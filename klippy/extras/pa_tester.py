# Pressure Advance calibration and tuning module
#
# Copyright (C) 2026  Dmitry Butyugin <dmbutyugin@google.com>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

import bisect, collections, importlib, logging, math, multiprocessing, traceback
import chelper, mathutil, queuelogger


MoveData = collections.namedtuple(
    'MoveData', ('print_time', 'move_t', 'start_v', 'accel'))

FitResult = collections.namedtuple(
    'FitResult', ('tau', 'status', 'mse', 'n_evals', 'n_fits'))

WindowData = collections.namedtuple(
    'WindowData', ('velocity', 'times', 'data'))



PA_WARNING_LOW = 0.005  # warn if PA is below this (possible detection issue)


class IntParam:
    def __init__(self, config, name, default=None, minval=None, maxval=None):
        self.name = name.upper()
        self._minval = minval
        self._maxval = maxval
        if default is None:
            self.config_value = config.getint(
                name, minval=self._minval, maxval=self._maxval)
        else:
            self.config_value = config.getint(
                name, default, minval=self._minval, maxval=self._maxval)
    def get(self, gcmd):
        return self.config_value if gcmd is None else \
                gcmd.get_int(self.name, self.config_value,
                             minval=self._minval, maxval=self._maxval)

class FloatParam:
    def __init__(self, config, name, default=None, minval=None, maxval=None,
                 above=None, below=None):
        self.name = name.upper()
        self._minval = minval
        self._maxval = maxval
        self._above = above
        self._below = below
        if default is None:
            self.config_value = config.getfloat(
                name, minval=self._minval, maxval=self._maxval,
                above=self._above, below=self._below)
        else:
            self.config_value = config.getfloat(
                name, default, minval=self._minval, maxval=self._maxval,
                above=self._above, below=self._below)
    def get(self, gcmd):
        return self.config_value if gcmd is None else \
                gcmd.get_float(self.name, self.config_value,
                               minval=self._minval, maxval=self._maxval,
                               above=self._above, below=self._below)

class PACalibrationExtruderTest:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.slow_flow = FloatParam(config, 'slow_flow', 2., above=0.)
        self.high_flow = FloatParam(config, 'high_flow', above=0.)
        self.purge_length = FloatParam(config, 'purge_length', above=0.)
        self.segment_time = FloatParam(config, 'segment_time', 0.5, above=0.)
        self.test_repetitions = IntParam(config, 'test_repetitions', 10,
                                         minval=1)

    def _flow_to_velocity(self, flow, filament_dia):
        area = math.pi * (filament_dia / 2.) ** 2
        return flow / area

    def gen_test(self, filament_dia, gcmd):
        slow_flow = self.slow_flow.get(gcmd)
        high_flow = self.high_flow.get(gcmd)
        if high_flow <= slow_flow:
            raise gcmd.error("high_flow must be greater than slow_flow=%.1f"
                             % slow_flow)
        if high_flow < 2. * slow_flow:
            gcmd.respond_info(
                "WARNING: high_flow is less than 2x slow_flow, PA "
                "estimation may be unreliable")
        slow_velocity = self._flow_to_velocity(slow_flow, filament_dia)
        fast_velocity = self._flow_to_velocity(high_flow, filament_dia)
        seg_time = self.segment_time.get(gcmd)
        # This is a bit approximate - we ignore accelerations and
        # decelerations of extruder here for distance calculation
        slow_dist = slow_velocity * seg_time
        fast_dist = fast_velocity * seg_time
        test_repetitions = self.test_repetitions.get(gcmd)
        segments = [(self.purge_length.get(gcmd), fast_velocity)]
        segments += [(slow_dist, slow_velocity),
                     (fast_dist, fast_velocity)] * test_repetitions
        segments.append((slow_dist, slow_velocity))
        return segments

    def run_test(self, segments):
        toolhead = self.printer.lookup_object('toolhead')
        tpos = toolhead.get_position()
        X, Y, Z = tpos[:3]
        E = tpos[3]
        for dist, velocity in segments:
            E += dist
            toolhead.move([X, Y, Z, E], velocity)
        toolhead.wait_moves()

class PACalibrationDataCollector:
    def __init__(self, printer, extruder, loadcell):
        self.printer = printer
        self.extruder = extruder
        self.loadcell = loadcell
        self.force_samples = []
        self.trapq_moves = []
        self._collecting = False

    def _on_loadcell_data(self, msg):
        if not self._collecting:
            return False
        data = msg.get('data')
        if data is not None:
            self.force_samples.extend(data)
        return True

    def _collect_trapq(self, start_time, end_time):
        ffi_main, ffi_lib = chelper.get_ffi()
        all_moves = []
        query_end = end_time
        while True:
            data = ffi_main.new('struct pull_move[128]')
            count = ffi_lib.trapq_extract_old(
                self.extruder.get_trapq(), data, len(data),
                start_time, query_end)
            if not count:
                break
            for i in range(count):
                m = data[i]
                all_moves.append(MoveData(m.print_time, m.move_t,
                                          m.start_v, m.accel))
            if count < len(data):
                break
            query_end = data[count - 1].print_time
        all_moves.reverse()
        return all_moves

    def collect(self, run_test_cb):
        toolhead = self.printer.lookup_object('toolhead')
        reactor = self.printer.get_reactor()
        start_time = toolhead.get_last_move_time()
        self._collecting = True
        self.loadcell.add_client(self._on_loadcell_data)
        run_test_cb()
        end_time = toolhead.get_last_move_time()
        timeout = reactor.monotonic() + 2.
        while not self.force_samples or self.force_samples[-1][0] < end_time:
            if reactor.monotonic() > timeout:
                break
            reactor.pause(reactor.monotonic() + 0.05)
        self._collecting = False
        self.trapq_moves = self._collect_trapq(start_time, end_time)
        return self.force_samples, self.trapq_moves


class TriangularWindowFilter:
    def __init__(self, numpy, filter_window):
        self.numpy = numpy
        self.filter_window = filter_window

    def smooth(self, times, values):
        np = self.numpy
        n = len(values)
        if n < 2:
            return values
        dt = np.mean(np.diff(times))
        half_w = 0.5 * self.filter_window
        h = max(math.ceil(half_w / dt), 1)
        if h < 2:
            return values
        frac = half_w / dt - h + 1.
        half_tri = np.arange(h, dtype=np.float64) + frac
        tri = np.concatenate([half_tri, half_tri[:-1][::-1]])
        m = len(tri)
        center = (m - 1) // 2
        conv = np.convolve(values, tri, mode='full')
        wsum = np.convolve(np.ones(n, dtype=np.float64), tri, mode='full')
        return conv[center:center + n] / wsum[center:center + n]


def velocity_key(v, places=4):
    return round(v, places)


def _match_samples_to_moves(moves, all_samples, move_window=0.):
    if not all_samples:
        return []
    all_samples.sort(key=lambda s: s[0])
    sample_times = [s[0] for s in all_samples]
    half_window = 0.5 * move_window
    results = []
    for move in moves:
        start_time = move.print_time
        end_time = move.print_time + move.move_t
        lo = bisect.bisect_left(sample_times, start_time - half_window)
        hi = bisect.bisect_right(sample_times, end_time + half_window)
        window = all_samples[lo:hi]
        if window:
            results.append(WindowData(
                velocity=move.start_v,
                times=[s[0] for s in window],
                data=[s[1] for s in window]))
    return results


def _background_process_exec(printer, method, args):
    if printer is None:
        return method(*args)
    parent_conn, child_conn = multiprocessing.Pipe()
    def wrapper():
        queuelogger.clear_bg_logging()
        try:
            res = method(*args)
        except:
            child_conn.send((True, traceback.format_exc()))
            child_conn.close()
            return
        child_conn.send((False, res))
        child_conn.close()
    calc_proc = multiprocessing.Process(target=wrapper)
    calc_proc.daemon = True
    calc_proc.start()
    reactor = printer.get_reactor()
    gcode = printer.lookup_object("gcode")
    eventtime = last_report_time = reactor.monotonic()
    while calc_proc.is_alive() and not parent_conn.poll():
        if eventtime > last_report_time + 5.:
            last_report_time = eventtime
            gcode.respond_info("Wait for calculations..", log=False)
        eventtime = reactor.pause(eventtime + .1)
    status, recv = parent_conn.recv()
    if status:
        raise printer.command_error(
            "Error in remote calculation: %s" % (recv,))
    calc_proc.join()
    parent_conn.close()
    return recv


class StepResponseTest:
    """Estimate PA by fitting A + B*exp(-C*t) + D*t to force data.

    Phase 1: fit each window of sensor data (with extruder constant
    velocity after flow change) independently to get tau=1/C estimate.

    Phase 2: jointly fit all windows per velocity, using the phase-1
    tau to determine the fit window, for a refined per-velocity tau.
    """

    PHASE1_FIT_FRACTION = 0.8  # fraction of window to use for phase 1 fit
    PHASE2_TAU_MULT = 4.       # phase 2 fit window in units of phase 1 tau
    MAX_FIT_MSE_PHASE1 = 0.001 # max normalized MSE to mark a fit as 'ok'
    WEIGHT_ALPHA = 10.         # relative weight of early vs. late samples
    WEIGHT_POWER = 2.          # exponent in the weighting function
    C_MIN = 0.1                # lower bound on C (gives PA_max = 1/C_MIN)
    C_MAX = 1000.              # upper bound on C (gives PA_min = 1/C_MAX)

    def __init__(self, numpy, config):
        self.numpy = numpy
        self.printer = config.get_printer()
        self.filter_window = FloatParam(config, 'filter_window', 0.01, above=0.)
        self.gcode = config.get_printer().lookup_object('gcode')

    def _filter_cv_moves(self, all_moves, skip_first=True):
        cv_moves = [m for m in all_moves if m.accel == 0.0]
        return cv_moves[(1 if skip_first else 0):]

    def _smooth_windows(self, windows, filt):
        np = self.numpy
        smoothed = []
        for window in windows:
            times = np.array(window.times, dtype=np.float64)
            forces = np.array(window.data, dtype=np.float64)
            s_data = filt.smooth(times, forces)
            smoothed.append(WindowData(
                velocity=window.velocity, times=times, data=s_data))
        return smoothed

    def _group_windows_by_velocity(self, windows):
        vel_groups = {}
        for window in windows:
            vk = velocity_key(window.velocity)
            vel_groups.setdefault(vk, []).append(window)
        return vel_groups

    def _prepare_window_data(self, windows, dt_limit):
        np = self.numpy
        window_data = []
        dt_max_global = 0.
        for window in windows:
            times = window.times
            data = window.data
            dt = times - times[0]
            skip_mask = dt < dt_limit
            dt_fit = dt[skip_mask]
            data_fit = data[skip_mask]
            if len(dt_fit) < 3:
                continue
            data_range = np.max(data_fit) - np.min(data_fit)
            if data_range < 1e-5:
                continue
            dt_max = dt_fit[-1]
            dt_max_global = max(dt_max_global, dt_max)
            p = self.WEIGHT_POWER
            gamma = self.WEIGHT_ALPHA ** (1. / p) - 1.
            sqrt_w = (gamma / (1. + gamma * dt_fit / dt_max)) ** (-0.5 * p)
            window_data.append((dt_fit, data_fit, sqrt_w, data_range))
        return window_data, dt_max_global

    def _make_cost(self, window_data):
        np = self.numpy
        def cost_func(C):
            if C <= self.C_MIN or C >= self.C_MAX:
                return 1e12
            try:
                total_cost = 0.
                for dt_fit, data_fit, sqrt_w, data_range in window_data:
                    exp_term = np.exp(-C * dt_fit)
                    sw = sqrt_w
                    X = np.column_stack((sw, sw * exp_term, sw * dt_fit))
                    yw = sw * data_fit
                    XtX = np.matmul(X.T, X)
                    Xty = np.matmul(X.T, yw)
                    sol = np.linalg.solve(XtX, Xty)
                    resid = yw - np.matmul(X, sol)
                    total_cost += np.sum(resid ** 2) / (data_range ** 2)
                return total_cost
            except Exception:
                return 1e12
        return cost_func

    def _build_c_candidates(self, C0, full=False):
        candidates = [C0, 2. * C0, 10. * C0]
        if full:
            candidates += [C0 / 10., C0 / 2.]
            candidates += [1., 10., 100.]
        distinct = set()
        for c in candidates:
            cr = round(c, 4)
            if cr > self.C_MIN and cr < self.C_MAX:
                distinct.add(cr)
        return sorted(distinct)

    def _run_fit(self, cost_func, n_total, c_candidates):
        adj_params = {'C'}
        best_params = None
        best_cost = float('inf')
        iter_count = [0]
        def counted_cost(params):
            iter_count[0] += 1
            return cost_func(params['C'])
        n_fits = len(c_candidates)
        for c_start in c_candidates:
            params = {'C': c_start}
            res = mathutil.coordinate_descent(adj_params, params, counted_cost)
            c_val = res['C']
            if c_val > self.C_MIN and c_val < self.C_MAX:
                cost = counted_cost(res)
                if cost < best_cost:
                    best_cost = cost
                    best_params = res
        if best_params is None:
            return FitResult(tau=None, status='c_out_of_bounds', mse=None,
                             n_evals=iter_count[0], n_fits=n_fits)
        C_opt = best_params['C']
        mse = best_cost / n_total
        tau = 1. / C_opt
        status = 'ok'
        if mse > self.MAX_FIT_MSE_PHASE1:
            status = 'high_mse'
        return FitResult(tau=tau, status=status, mse=mse,
                         n_evals=iter_count[0], n_fits=n_fits)

    def _fit_window_phase1(self, window):
        dt = window.times - window.times[0]
        dt_limit = dt[-1] * self.PHASE1_FIT_FRACTION
        window_data, dt_max = self._prepare_window_data([window], dt_limit)
        if not window_data:
            return FitResult(tau=None, status='too_few_after_skip',
                             mse=None, n_evals=0, n_fits=0)
        n_total = len(window_data[0][0])
        cost_func = self._make_cost(window_data)
        c_candidates = self._build_c_candidates(1. / dt_max, full=False)
        return self._run_fit(cost_func, n_total, c_candidates)

    def _build_phase1_report(self, vel_fits, velocity_taus, overall_tau,
                             n_evals, n_fits):
        ok_count = sum(1 for _, f in vel_fits if f.status == 'ok')
        fail_count = len(vel_fits) - ok_count
        lines = []
        lines.append("--- Phase 1 Fit Results (rough) ---")
        lines.append("  %-4s %-6s %-8s %-8s %-8s" % (
            'W#', 'v', 'tau', 'MSE', 'status'))
        lines.append("  " + "-" * 40)
        for idx, (vk, fit) in enumerate(vel_fits):
            tau_str = "%.4f" % fit.tau if fit.tau is not None else 'N/A'
            mse_str = "%.1e" % fit.mse if fit.mse is not None else 'N/A'
            lines.append("  %-4d %-6.2f %-8s %-8s %-8s" % (
                idx, vk, tau_str, mse_str, fit.status))
        lines.append("")
        lines.append("  Phase 1 fits succeeded: %d  failed: %d  total: %d" % (
            ok_count, fail_count, ok_count + fail_count))
        lines.append("")
        if velocity_taus:
            lines.append("  Phase 1 tau by velocity:")
            for v in sorted(velocity_taus.keys()):
                lines.append(
                    "    v=%.2f  mean_tau=%.6f" % (v, velocity_taus[v]))
            lines.append("")
        if overall_tau is not None:
            lines.append("  overall_tau=%.6f" % overall_tau)
            lines.append("")
        avg_evals = n_evals / n_fits if n_fits else 0
        lines.append("  Optimization func evals: %d total, %.1f avg/call" % (
            n_evals, avg_evals))
        lines.append("")
        return lines

    def _run_phase1(self, vel_groups):
        all_fits = []
        for vk, group in vel_groups.items():
            for window in group:
                fit = self._fit_window_phase1(window)
                all_fits.append((vk, fit))
        velocity_taus, overall_tau = self._aggregate_taus(all_fits)
        n_evals = sum(fit.n_evals for _, fit in all_fits)
        n_fits = sum(fit.n_fits for _, fit in all_fits)
        report_lines = self._build_phase1_report(
                all_fits, velocity_taus, overall_tau, n_evals, n_fits)
        return all_fits, report_lines

    def _fit_velocity_phase2(self, windows, approx_tau):
        dt_limit = self.PHASE2_TAU_MULT * approx_tau
        window_data, _ = self._prepare_window_data(windows, dt_limit)
        if not window_data:
            return FitResult(tau=None, status='no_valid_windows',
                             mse=None, n_evals=0, n_fits=0)
        n_total = sum(len(wd[0]) for wd in window_data)
        cost_func = self._make_cost(window_data)
        c_candidates = self._build_c_candidates(1. / approx_tau, full=True)
        return self._run_fit(cost_func, n_total, c_candidates)

    def _build_phase2_report(self, vel_fits, velocity_taus, overall_tau,
                             n_evals, n_fits):
        ok_count = sum(1 for _, f in vel_fits if f.status == 'ok')
        fail_count = len(vel_fits) - ok_count
        lines = []
        lines.append("--- Fit Results ---")
        lines.append("  %-4s %-6s %-8s %-8s %-8s" % (
            'V#', 'v', 'tau', 'MSE', 'status'))
        lines.append("  " + "-" * 40)
        for idx, (vk, fit) in enumerate(vel_fits):
            tau_str = "%.4f" % fit.tau if fit.tau is not None else 'N/A'
            mse_str = "%.1e" % fit.mse if fit.mse is not None else 'N/A'
            lines.append("  %-4d %-6.2f %-8s %-8s %-8s" % (
                idx, vk, tau_str, mse_str, fit.status))
        lines.append("")
        lines.append("  Fits succeeded: %d  failed: %d  total: %d" % (
            ok_count, fail_count, ok_count + fail_count))
        lines.append("")
        if velocity_taus:
            lines.append("  Phase 2 tau by velocity:")
            for v in sorted(velocity_taus.keys()):
                lines.append(
                        "    v=%.2f  tau=%.6f" % (v, velocity_taus[v]))
            lines.append("")
        if overall_tau is not None:
            lines.append("  overall_tau=%.6f" % overall_tau)
            lines.append("")
        avg_evals = n_evals / n_fits if n_fits else 0
        lines.append("  Optimization func evals: %d total, %.1f avg/call" % (
            n_evals, avg_evals))
        lines.append("")
        return lines

    def _run_phase2(self, vel_groups, velocity_taus):
        phase2_fits = []
        for vk, group in vel_groups.items():
            fit = self._fit_velocity_phase2(group, velocity_taus[vk])
            phase2_fits.append((vk, fit))
        p2_velocity_taus, p2_overall_tau = self._aggregate_taus(phase2_fits)
        n_evals = sum(fit.n_evals for _, fit in phase2_fits)
        n_fits = sum(fit.n_fits for _, fit in phase2_fits)
        report_lines = self._build_phase2_report(
                phase2_fits, p2_velocity_taus, p2_overall_tau, n_evals, n_fits)
        return phase2_fits, report_lines

    def _aggregate_taus(self, vel_fits):
        np = self.numpy
        vel_taus = {}
        for vk, fit in vel_fits:
            if fit.status == 'ok' and fit.tau is not None:
                vel_taus.setdefault(vk, []).append(fit.tau)
        agg_taus = {vk: np.mean(taus) for vk, taus in vel_taus.items()}
        overall_tau = np.mean(list(agg_taus.values())) if agg_taus else None
        return agg_taus, overall_tau

    def _respond(self, gcmd, velocity_taus, overall_tau, filament_diameter,
                 segment_time, filter_window):
        area = math.pi * (filament_diameter / 2.) ** 2
        lines = ["PA calibration results:"]
        for v in sorted(velocity_taus.keys()):
            flow = v * area
            lines.append("  flow=%5.2f  pressure_advance=%.3f" % (
                flow, velocity_taus[v]))
        if overall_tau is not None:
            lines.append("  overall pressure_advance=%.3f" % overall_tau)
            if overall_tau < PA_WARNING_LOW:
                lines.append(
                    "WARNING: detected PA is very low (%.3f). This may "
                    "indicate an issue with the sensor, hotend, or filament."
                    % overall_tau)
            if overall_tau < filter_window:
                lines.append(
                    "WARNING: detected PA (%.3f) is below filter_window "
                    "(%.4f). The signal filter may be distorting the "
                    "pressure response; try reducing filter_window."
                    % (overall_tau, filter_window))
            if overall_tau > 0.5 * segment_time:
                lines.append(
                    "WARNING: detected PA (%.3f) exceeds half of segment_time "
                    "(%.1f). The exponential decay may not have settled; "
                    "increase segment_time for a more reliable result."
                    % (overall_tau, segment_time))
        gcmd.respond_info("\n".join(lines))

    def run(self, gcmd, filament_diameter, extruder_test, data_collector):
        reactor = self.printer.get_reactor()
        self.gcode.run_script_from_command("SET_PRESSURE_ADVANCE ADVANCE=0")
        # Execute the extrusion test and collect force and move data
        segments = extruder_test.gen_test(filament_diameter, gcmd)
        all_samples, all_moves = data_collector.collect(
            lambda: extruder_test.run_test(segments))
        # Keep only constant-velocity moves and match force samples to them
        cv_moves = self._filter_cv_moves(all_moves)
        windows = _match_samples_to_moves(cv_moves, all_samples)
        # Smooth force data of each window to reduce sensor noise
        filt = TriangularWindowFilter(self.numpy,
                                      self.filter_window.get(gcmd))
        windows = self._smooth_windows(windows, filt)
        vel_groups = self._group_windows_by_velocity(windows)
        t0 = reactor.monotonic()
        # Phase 1: rough tau fit per window
        vel_fits, phase1_report_lines = _background_process_exec(
                self.printer, self._run_phase1, (vel_groups,))
        t1 = reactor.monotonic()
        for line in phase1_report_lines:
            logging.info(line)
        logging.info("Phase 1 wall time: %.3f s", t1 - t0)
        velocity_taus, _ = self._aggregate_taus(vel_fits)
        if len(velocity_taus) < len(vel_groups):
            raise self.printer.command_error(
                "PA calibration failed: force sensor data is noisy or "
                "unreliable. Check sensor connections, hotend mounting, "
                "recalibrate the sensor, and try again.")
        # Phase 2: refined joint fit per velocity group
        p2_fits, phase2_report_lines = _background_process_exec(
                self.printer, self._run_phase2,
                (vel_groups, velocity_taus))
        t2 = reactor.monotonic()
        for line in phase2_report_lines:
            logging.info(line)
        logging.info("Phase 2 wall time: %.3f s", t2 - t1)
        p2_velocity_taus, p2_overall_tau = self._aggregate_taus(p2_fits)
        segment_time = extruder_test.segment_time.get(gcmd)
        self._respond(gcmd, p2_velocity_taus, p2_overall_tau, filament_diameter,
                      segment_time, signal_filter.filter_window)


PA_CALIBRATION_METHODS_MAP = {
    'step_response': StepResponseTest,
}


class PATester:
    def __init__(self, config):
        self.printer = config.get_printer()
        try:
            self.numpy = importlib.import_module('numpy')
        except ImportError:
            raise self.printer.command_error(
                "Failed to import `numpy` module, make sure it was "
                "installed via `~/klippy-env/bin/pip install`")
        self.force_sensor = config.get('force_sensor', None)
        self.extruder = config.get('extruder', None)
        method = config.get('method', 'step_response').lower()
        if method not in PA_CALIBRATION_METHODS_MAP:
            raise config.error(
                    "Invalid method '%s' (available: %s)"
                    % (method, ', '.join(
                        sorted(PA_CALIBRATION_METHODS_MAP.keys()))))
        self.default_method = method
        self.extruder_test = PACalibrationExtruderTest(config)
        self.methods = {
            name: cls(self.numpy, config)
            for name, cls in PA_CALIBRATION_METHODS_MAP.items()
        }
        self.gcode = self.printer.lookup_object('gcode')
        self.gcode.register_command(
                "PA_CALIBRATE", self.cmd_PA_CALIBRATE,
                desc=self.cmd_PA_CALIBRATE_help)

    def _resolve_extruder(self, gcmd):
        toolhead = self.printer.lookup_object('toolhead')
        active_extruder = toolhead.get_extruder()
        if active_extruder is None:
            raise self.printer.command_error(
                "No active extruder, the tested extruder must be activated")
        active_name = active_extruder.get_name()
        name = gcmd.get('EXTRUDER', self.extruder)
        if name is None:
            name = active_name
        elif name != active_name:
            raise self.printer.command_error(
                "Another extruder '%s' is currently active" % (active_name,))
        extruder = self.printer.lookup_object(name, None)
        if extruder is None:
            raise self.printer.command_error(
                "Extruder '%s' not found" % (name,))
        return extruder

    def _resolve_loadcell(self, gcmd):
        name = gcmd.get('FORCE_SENSOR', self.force_sensor)
        if name is None:
            name = self._auto_detect_loadcell()
        if name is None:
            raise self.printer.command_error(
                "No loadcell: configure 'force_sensor' in [pa_tester] "
                "or pass FORCE_SENSOR= to PA_CALIBRATE")
        obj = self.printer.lookup_object(name, None)
        if obj is None:
            raise self.printer.command_error(
                "Force sensor '%s' not found" % (name,))
        if hasattr(obj, 'get_load_cell'):
            loadcell = obj.get_load_cell()
        else:
            loadcell = obj
        return loadcell, name

    def _auto_detect_loadcell(self):
        lc_objects = self.printer.lookup_objects("load_cell_probe")
        lc_objects.extend(self.printer.lookup_objects("load_cell"))
        if not lc_objects:
            return None
        if len(lc_objects) > 1:
            names = ', '.join(n for n, o in lc_objects)
            raise self.printer.command_error(
                "Multiple loadcells found (%s): configure 'force_sensor' "
                "in [pa_tester] or pass FORCE_SENSOR= to PA_CALIBRATE"
                % (names,))
        return lc_objects[0][0]

    cmd_PA_CALIBRATE_help = "Calibrate Pressure Advance using loadcell data"
    def cmd_PA_CALIBRATE(self, gcmd):
        method_name = gcmd.get('METHOD', self.default_method).lower()
        method = self.methods.get(method_name)
        if method is None:
            raise self.printer.command_error(
                "Unknown PA_CALIBRATE method '%s' (available: %s)"
                % (method_name, ', '.join(sorted(self.methods.keys()))))
        extruder = self._resolve_extruder(gcmd)
        loadcell, _ = self._resolve_loadcell(gcmd)
        systime = self.printer.get_reactor().monotonic()
        extruder_status = extruder.get_status(systime)
        if not extruder_status['can_extrude']:
            raise self.printer.command_error(
                "Extruder cannot extrude\n"
                "See the 'min_extrude_temp' config option for details")
        filament_diameter = extruder_status['filament_diameter']
        collector = PACalibrationDataCollector(self.printer, extruder, loadcell)
        method.run(gcmd, filament_diameter, self.extruder_test, collector)


def load_config(config):
    return PATester(config)
