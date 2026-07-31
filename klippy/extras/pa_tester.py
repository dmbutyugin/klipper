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
    'WindowData', ('velocity', 'times', 'data', 'accel'))

PATestVelocityResult = collections.namedtuple(
    'PATestVelocityResult', ('velocity', 'overshoot', 'ratios', 'n_overshoot'))

PATestResult = collections.namedtuple(
    'PATestResult', ('pressure_advance', 'overshoot', 'vel_results'))


PA_WARNING_LOW = 0.005  # warn if PA is below this (possible detection issue)
FORCE_EPS = 1e-5        # effective zero force signal


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

    def get_slow_velocity(self, gcmd, filament_dia):
        return self._flow_to_velocity(self.slow_flow.get(gcmd), filament_dia)

    def get_fast_velocity(self, gcmd, filament_dia):
        return self._flow_to_velocity(self.high_flow.get(gcmd), filament_dia)

    def gen_test(self, gcmd, filament_dia, full_purge=True,
                 test_repetitions=None, gen_force_sig_test=False):
        slow_velocity = self.get_slow_velocity(gcmd, filament_dia)
        fast_velocity = self.get_fast_velocity(gcmd, filament_dia)
        if fast_velocity <= slow_velocity:
            raise gcmd.error("high_flow must be greater than slow_flow")
        if fast_velocity < 2. * slow_velocity:
            gcmd.respond_info(
                "WARNING: high_flow is less than 2x slow_flow, PA "
                "estimation may be unreliable")
        seg_time = self.segment_time.get(gcmd)
        # This is a bit approximate - we ignore accelerations and
        # decelerations of extruder here for distance calculation
        slow_dist = slow_velocity * seg_time
        fast_dist = fast_velocity * seg_time
        reps = (test_repetitions if test_repetitions is not None
                else self.test_repetitions.get(gcmd))
        extruder_moves = []
        if full_purge:
            extruder_moves.append((self.purge_length.get(gcmd), fast_velocity))
        else:
            extruder_moves.append((fast_dist, fast_velocity))
        extruder_moves += [(slow_dist, slow_velocity),
                           (fast_dist, fast_velocity)] * reps
        if gen_force_sig_test:
            medium_velocity = 0.5 * (slow_velocity + fast_velocity)
            medium_dist = medium_velocity * seg_time
            extruder_moves.append((medium_dist, medium_velocity))
            extruder_moves.append((slow_dist, slow_velocity))
        return extruder_moves

    def run_test(self, extruder_moves):
        toolhead = self.printer.lookup_object('toolhead')
        tpos = toolhead.get_position()
        X, Y, Z = tpos[:3]
        E = tpos[3]
        for dist, velocity in extruder_moves:
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


def detect_force_sign(all_samples, all_moves):
    # Detect force sensor polarity from the two last CV moves
    # (the very last one is slower and should require less force)
    if not all_moves:
        return None
    cv_moves = [m for m in all_moves if m.accel == 0.0]
    last_moves = cv_moves[-2:]
    if not last_moves:
        return None
    times = (last_moves[0].print_time, last_moves[1].print_time,
             last_moves[1].print_time + last_moves[1].move_t)
    last_move_samples = [s[1] for s in all_samples
                         if times[1] <= s[0] < times[2]]
    prev_move_samples = [s[1] for s in all_samples
                         if times[0] <= s[0] < times[1]]
    if len(last_move_samples) < 2 or len(prev_move_samples) < 2:
        return None
    diff = (sum(prev_move_samples) / len(prev_move_samples)
            - sum(last_move_samples) / len(last_move_samples))
    if diff > FORCE_EPS:
        return 1
    if diff < -FORCE_EPS:
        return -1
    return None


def velocity_key(v, places=4):
    return round(v, places)


def _match_samples_to_moves(moves, all_samples, move_extra_window=0.):
    if not all_samples:
        return []
    all_samples.sort(key=lambda s: s[0])
    sample_times = [s[0] for s in all_samples]
    half_window = 0.5 * move_extra_window
    results = []
    for move, window_accel in moves:
        start_time = move.print_time
        end_time = move.print_time + move.move_t
        lo = bisect.bisect_left(sample_times, start_time - half_window)
        hi = bisect.bisect_right(sample_times, end_time + half_window)
        window = all_samples[lo:hi]
        if window:
            results.append(WindowData(
                velocity=move.start_v + move.accel * move.move_t,
                times=[s[0] for s in window],
                data=[s[1] for s in window],
                accel=window_accel))
    return results

def _smooth_windows(np, windows, filt):
    smoothed = []
    for window in windows:
        times = np.array(window.times, dtype=np.float64)
        forces = np.array(window.data, dtype=np.float64)
        s_data = filt.smooth(times, forces)
        smoothed.append(WindowData(
            velocity=window.velocity, times=times, data=s_data,
            accel=window.accel))
    return smoothed


def _group_windows_by_velocity(windows):
    vel_groups = {}
    for window in windows:
        vk = velocity_key(window.velocity)
        vel_groups.setdefault(vk, []).append(window)
    return vel_groups


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

    def __init__(self, printer, numpy, scipy=None):
        self.printer = printer
        self.numpy = numpy
        self.scipy = scipy
        self.gcode = printer.lookup_object('gcode')

    def _filter_cv_moves(self, all_moves):
        cv_moves = []
        prior_accel = 0.
        for move in all_moves:
            if move.accel == 0.0:
                cv_moves.append((move, prior_accel))
            prior_accel = move.accel
        return cv_moves

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
            if data_range < FORCE_EPS:
                continue
            dt_max = dt_fit[-1]
            dt_max_global = max(dt_max_global, dt_max)
            p = self.WEIGHT_POWER
            gamma = self.WEIGHT_ALPHA ** (1. / p) - 1.
            sqrt_w = (gamma / (1. + gamma * dt_fit / dt_max)) ** (-0.5 * p)
            window_data.append((dt_fit, data_fit, sqrt_w, data_range))
        return window_data, dt_max_global

    def _solve_wls(self, window_data, C):
        # Solve A, B, D for model A + B*exp(-C*t) + D*t via WLS
        np = self.numpy
        dt_fit, data_fit, sqrt_w, data_range = window_data
        exp_term = np.exp(-C * dt_fit)
        sw = sqrt_w
        X = np.column_stack((sw, sw * exp_term, sw * dt_fit))
        yw = sw * data_fit
        XtX = np.matmul(X.T, X)
        Xty = np.matmul(X.T, yw)
        sol = np.linalg.solve(XtX, Xty)
        resid = yw - np.matmul(X, sol)
        return sol, resid

    def _make_cost(self, window_data):
        np = self.numpy
        def cost_func(C):
            if C <= self.C_MIN or C >= self.C_MAX:
                return 1e12
            try:
                total_cost = 0.
                for wd in window_data:
                    _, resid = self._solve_wls(wd, C)
                    data_range = wd[3]
                    total_cost += np.sum(resid ** 2) / (data_range ** 2)
                return total_cost
            except Exception:
                return 1e12
        return cost_func

    def _build_c_candidates(self, C0, full=False):
        candidates = [10., C0, 2. * C0, 10. * C0]
        if full:
            candidates += [C0 / 10., C0 / 2.]
            candidates += [1., 100.]
        distinct = set()
        for c in candidates:
            cr = round(c, 4)
            if cr > self.C_MIN and cr < self.C_MAX:
                distinct.add(cr)
        return sorted(distinct)

    def _run_fit(self, cost_func, n_total, c_candidates):
        best_params = {}
        best_cost = float('inf')
        if self.scipy is not None:
            res = self.scipy.optimize.minimize_scalar(
                    cost_func, bounds=(self.C_MIN, self.C_MAX),
                    method='bounded', options={'xatol': 1e-8})
            n_evals, n_fits = res.nfev, 1
            best_params['C'] = res.x
            best_cost = res.fun
        else:
            adj_params = {'C'}
            iter_count = [0]
            def counted_cost(params):
                iter_count[0] += 1
                return cost_func(params['C'])
            n_fits = len(c_candidates)
            for c_start in c_candidates:
                res = mathutil.coordinate_descent(
                        adj_params, {'C': c_start}, counted_cost)
                cost = counted_cost(res)
                if cost < best_cost:
                    best_cost = cost
                    best_params = res
            n_evals = iter_count[0]
        C_opt = best_params['C']
        mse = best_cost / n_total
        tau = 1. / C_opt
        status = 'ok'
        if mse > self.MAX_FIT_MSE_PHASE1:
            status = 'high_mse'
        return FitResult(tau=tau, status=status, mse=mse,
                         n_evals=n_evals, n_fits=n_fits)

    def _fit_window_phase1(self, window, force_sign):
        dt = window.times - window.times[0]
        dt_limit = dt[-1] * self.PHASE1_FIT_FRACTION
        window_data, dt_max = self._prepare_window_data([window], dt_limit)
        if not window_data:
            return FitResult(tau=None, status='too_few_after_skip',
                             mse=None, n_evals=0, n_fits=0)
        n_total = len(window_data[0][0])
        cost_func = self._make_cost(window_data)
        c_candidates = self._build_c_candidates(1. / dt_max, full=False)
        fit = self._run_fit(cost_func, n_total, c_candidates)
        # Validate B decay factor sign: B should oppose force_sign and accel
        sol, _ = self._solve_wls(window_data[0], 1. / fit.tau)
        if force_sign * window.accel * sol[1] > 0:
            return FitResult(tau=fit.tau, status='wrong_b_sign', mse=fit.mse,
                             n_evals=fit.n_evals, n_fits=fit.n_fits)
        return fit

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

    def _run_phase1(self, vel_groups, force_sign):
        all_fits = []
        for vk, group in vel_groups.items():
            for window in group:
                fit = self._fit_window_phase1(window, force_sign)
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

    def run(self, gcmd, extruder_test, data_collector, extruder_status,
            signal_filter, force_sign):
        reactor = self.printer.get_reactor()
        self.gcode.run_script_from_command("SET_PRESSURE_ADVANCE ADVANCE=0")
        # Execute the extrusion test and collect force and move data
        filament_diameter = extruder_status['filament_diameter']
        segment_time = extruder_test.segment_time.get(gcmd)
        extruder_moves = extruder_test.gen_test(gcmd, filament_diameter,
                                                gen_force_sig_test=True)
        all_samples, all_moves = data_collector.collect(
                lambda: extruder_test.run_test(extruder_moves))
        # Detect force_sign from the the last two of the extruder_moves
        if not force_sign:
            force_sign = detect_force_sign(all_samples, all_moves)
            if not force_sign:
                raise gcmd.error("Could not detect the sensor force direction. "
                                 "Check sensor and filament and try again.")
            gcmd.respond_info("Detected extrude_force_sign=%d" % force_sign)
        # Keep only constant-velocity moves skipping the first one (purge) and
        # the two last ones (force sign test), and match force samples to them
        cv_moves = self._filter_cv_moves(all_moves)[1:-2]
        if not cv_moves:
            raise self.printer.command_error("No move data captured")
        windows = _match_samples_to_moves(cv_moves, all_samples)
        if not windows:
            raise self.printer.command_error(
                "No force samples matched to extruder moves")
        # Smooth force data of each window to reduce sensor noise
        smoothed_windows = _smooth_windows(self.numpy, windows, signal_filter)
        vel_groups = _group_windows_by_velocity(smoothed_windows)
        t0 = reactor.monotonic()
        # Phase 1: rough tau fit per window
        vel_fits, phase1_report_lines = _background_process_exec(
                self.printer, self._run_phase1, (vel_groups, force_sign))
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
        self._respond(gcmd, p2_velocity_taus, p2_overall_tau, filament_diameter,
                      segment_time, signal_filter.filter_window)


class SearchOvershootTest:
    """Find optimal PA via binary-like search for the overshoot threshold."""

    STEADY_STATE_TAIL = 0.2          # fraction of tail for steady-state mean
    OVERSHOOT_THRESHOLD = 0.99       # integral threshold value for an overshoot
    OVERSHOOT_VOTE_FRACTION = 0.4    # fraction of windows to declare overshoot
    FORCE_SMOOTH_TIME_SETTLE = 2.0   # extra window around flow change to test
    ZERO_VEL_THRESHOLD = 0.1         # min velocity to keep an accel move
    PA_SEARCH_START = 0.05           # starting PA for exponential bracketing
    PA_SEARCH_MAX = 2.0              # hard upper limit for PA search
    PA_SEARCH_MIN_INTERVAL = 0.002   # absolute convergence tolerance
    PA_SEARCH_REL_TOL = 0.01         # relative convergence tolerance
    PA_SEARCH_SHRINK_FACTOR = 0.33   # fraction of interval shifted each step
    PA_RANGE_WARNING_FRACTION = 0.10 # PA proximity to boundary for warning
    PA_WARNING_HIGH_MARGIN = 0.2     # warn if PA within this of PA_SEARCH_MAX
    ACCEL_DURATION_SMOOTH_TIME_MULT = 1. # accel duration, in smooth_time units
    PA_CORRECTED_REL_START = 0.5     # relative corrected PA search start
    PA_CORRECTED_REL_TOL = 0.001     # relative corrected PA search tolerance
    SIMULATION_SEG_TIME = 0.0001     # simulation time step (s)
    CORRECTION_MIN_FRACTION = 0.05   # min PA correction fraction to report

    def __init__(self, printer, numpy, scipy=None):
        self.printer = printer
        self.numpy = numpy
        self.scipy = scipy
        self.gcode = printer.lookup_object('gcode')

    def _filter_accel_moves(self, all_moves):
        accel_moves = []
        for m in all_moves:
            if m.accel == 0.:
                continue
            end_v = m.start_v + m.accel * m.move_t
            if min(abs(m.start_v), abs(end_v)) < self.ZERO_VEL_THRESHOLD:
                continue
            accel_moves.append((m, m.accel))
        return accel_moves

    def _get_force_settle_time(self, smooth_time, signal_filter):
        return self.FORCE_SMOOTH_TIME_SETTLE * max(smooth_time,
                                                   signal_filter.filter_window)

    def _estimate_steady_state(self, forces):
        np = self.numpy
        n = len(forces)
        tail_start = int(n * (1. - self.STEADY_STATE_TAIL))
        if tail_start < 1:
            tail_start = 1
        return np.mean(forces[tail_start:])

    def _check_window_overshoot(self, detrended, accel, force_sign):
        np = self.numpy
        n = len(detrended)
        if n < 2:
            return None
        # Ignore the last half of detrending window, helps detecting strong
        # undershoot when the pressure is still settling in the tail.
        data_end = int(n * (1. - self.STEADY_STATE_TAIL * 0.5))
        data = detrended[:data_end]
        mean_abs = np.mean(np.abs(data))
        if mean_abs < FORCE_EPS:
            return None
        sign = force_sign if accel < 0 else -force_sign
        signed_mean = sign * np.mean(data)
        ratio = signed_mean / mean_abs
        return ratio < self.OVERSHOOT_THRESHOLD, ratio

    def _check_overshoot(self, votes):
        valid = [v for v in votes if v is not None]
        if not valid:
            return False, [], 0
        overshoot_votes = [v[0] for v in valid]
        ratios = [v[1] for v in valid]
        n_overshoot = sum(overshoot_votes)
        has_overshoot = (n_overshoot / len(overshoot_votes)
                         >= self.OVERSHOOT_VOTE_FRACTION)
        return has_overshoot, ratios, n_overshoot

    def test_pa(self, gcmd, extruder_test, data_collector,
                extruder_status, signal_filter, pressure_advance, force_sign):
        gcmd.respond_info(
                "Test %d: pressure_advance=%.4f ..."
                % (self._test_count + 1, pressure_advance))
        filament_diameter = extruder_status['filament_diameter']
        smooth_time = extruder_status['smooth_time']
        extruder = data_collector.extruder
        self.gcode.run_script_from_command(
            "SET_PRESSURE_ADVANCE ADVANCE=%.6f" % pressure_advance)
        old_pa_flag = extruder.set_enable_pa_for_extrude_only_moves(True)
        extruder_moves = extruder_test.gen_test(
                gcmd, filament_diameter, full_purge=not self._extruder_purged)
        velocities = [v for _, v in extruder_moves]
        all_samples, all_moves = data_collector.collect(
            lambda: extruder_test.run_test(extruder_moves))
        self._extruder_purged = True
        extruder.set_enable_pa_for_extrude_only_moves(old_pa_flag)
        accel_moves = self._filter_accel_moves(all_moves)
        if not accel_moves:
            raise self.printer.command_error("No move data captured")
        match_window = 2. * self._get_force_settle_time(smooth_time,
                                                        signal_filter)
        windows = _match_samples_to_moves(
                accel_moves, all_samples, match_window)
        if not windows:
            raise self.printer.command_error(
                "No force samples matched to extruder moves")
        smoothed_windows = _smooth_windows(self.numpy, windows, signal_filter)
        vel_groups = _group_windows_by_velocity(smoothed_windows)
        vel_results = []
        detected_overshoot = False
        for vk, group in vel_groups.items():
            votes = []
            for window in group:
                steady = self._estimate_steady_state(window.data)
                detrended = window.data - steady
                votes.append(self._check_window_overshoot(
                    detrended, window.accel, force_sign))
            has_overshoot, ratios, n_overshoot = self._check_overshoot(votes)
            vel_results.append(PATestVelocityResult(
                velocity=vk, overshoot=has_overshoot,
                ratios=ratios, n_overshoot=n_overshoot))
            if has_overshoot:
                detected_overshoot = True
            logging.info(
                    "PA=%.4f vel=%.2f overshoot=%d n_overshoot=%d"
                    " ratios=[%.4f..%.4f]" % (pressure_advance, vk,
                                              has_overshoot, n_overshoot,
                                              min(ratios) if ratios else 0.,
                                              max(ratios) if ratios else 0.))
        self._test_count += 1
        gcmd.respond_info("Tested pressure_advance=%.4f -> %s" % (
            pressure_advance,
            'overshoot' if detected_overshoot else 'no overshoot'))
        return PATestResult(
            pressure_advance=pressure_advance,
            overshoot=detected_overshoot,
            vel_results=vel_results)

    def _bracket_pa(self, test_cb):
        pa = self.PA_SEARCH_START
        left = 0.
        res = test_cb(pa)
        overshoot = res.overshoot
        # Exponentially increase PA limit until overshoot is confirmed twice
        while pa <= self.PA_SEARCH_MAX:
            res = test_cb(2. * pa)
            if res.overshoot:
                if overshoot:
                    return left, 1.5 * pa
                overshoot = True
            else:
                overshoot = False
                left = pa
            pa *= 2.
        raise self.printer.command_error(
            "No overshoot detected up to pressure_advance=%.2f. "
            "Check sensor or try a higher flow rate."
            % self.PA_SEARCH_MAX)

    def _contracting_search(self, test_cb, left, right):
        while True:
            m = 0.5 * (left + right)
            width = right - left
            tol = max(self.PA_SEARCH_MIN_INTERVAL, self.PA_SEARCH_REL_TOL * m)
            if width <= tol:
                break
            adjust = width * self.PA_SEARCH_SHRINK_FACTOR
            res = test_cb(m)
            if res.overshoot:
                right -= adjust
            else:
                left += adjust
        return m

    def _search_optimal_pa(self, gcmd, extruder_test, data_collector,
                           extruder_status, signal_filter, force_sign,
                           pa_range=None):
        results = []
        def test_cb(pa):
            res = self.test_pa(
                gcmd, extruder_test, data_collector,
                extruder_status, signal_filter, pa, force_sign)
            results.append(res)
            return res
        # First detect the PA range (if not provided), then run a binary-like
        # search to determine an optimal PA value without overshoot
        if pa_range is not None:
            left, right = pa_range
        else:
            left, right = self._bracket_pa(test_cb)
        optimal_pa = self._contracting_search(test_cb, left, right)
        return optimal_pa, results

    def _log_search_stats(self, results, optimal_pa, corrected_pa, area):
        rows = []
        for idx, res in enumerate(results):
            for vr in res.vel_results:
                rows.append((
                    idx + 1, res.pressure_advance, vr.velocity,
                    min(vr.ratios) if vr.ratios else 0.,
                    max(vr.ratios) if vr.ratios else 0.,
                    vr.overshoot, vr.n_overshoot))
        rows.sort(key=lambda r: (r[1], r[2]))
        logging.info("--- PA Search Details ---")
        logging.info("  %-4s %-10s %-8s %-12s %-8s",
                      '#', 'PA', 'vel', 'ratios', 'overshoot')
        logging.info("  " + "-" * 44)
        prev_pa = None
        marked_optimal = False
        for tnum, pa, vk, min_r, max_r, is_overshoot, n_overshoot in rows:
            if not marked_optimal and pa >= optimal_pa:
                logging.info("  >> optimal PA ~ %.4f << " % optimal_pa)
                marked_optimal = True
            label = ('O (%d)' if is_overshoot else 'U (%d)') % n_overshoot
            if pa != prev_pa:
                logging.info("  %-4d %-10.4f %-8.2f %-12s %-8s" % (
                    tnum, pa, vk, "[%.4f..%.4f]" % (min_r, max_r), label))
                prev_pa = pa
            else:
                logging.info("  %-4s %-10s %-8.2f %-12s %-8s" % (
                    '', '', vk, "[%.4f..%.4f]" % (min_r, max_r), label))
        logging.info(
                "  Simulation-corrected PA: %.4f (raw: %.4f)" % (
                    corrected_pa, optimal_pa))

    # Simulation helpers for PA correction

    def _calc_pa(self, times, positions, extruder_pa, smooth_time):
        # Calculate PA-compensated extruder positions
        np = self.numpy
        pa_dt = extruder_pa / self.SIMULATION_SEG_TIME
        pos_diffs = np.diff(positions)
        pa_positions = positions + pa_dt * np.append(pos_diffs, pos_diffs[-1])
        pa_filter = TriangularWindowFilter(np, smooth_time)
        # Smooth velocity and not positions, since triangular filter will smooth
        # linear function incorrectly at boundaries (result is non-linear)
        return np.cumsum(np.append(
            pa_positions[0], pa_filter.smooth(times, np.diff(pa_positions))))

    def _calc_nozzle_flow(self, positions, system_pa):
        # Simulate nozzle first-order pressure response
        np = self.numpy
        dt_pa = self.SIMULATION_SEG_TIME / system_pa
        if self.scipy is not None:
            return self.scipy.signal.lfilter(
                    [dt_pa], [1., -(1. - dt_pa)], positions)
        n = len(positions)
        out = np.empty(n, dtype=np.float64)
        out[0] = last = positions[0]
        for i in range(1, n):
            last += (positions[i] - last) * dt_pa
            out[i] = last
        return out

    def _simulate_velocity_step(self, slow_vel, fast_vel, accel, smooth_time,
                                signal_filter):
        # Generate a synthetic velocity step from slow to fast velocity
        np = self.numpy
        seg_time = self.SIMULATION_SEG_TIME
        accel_time = (fast_vel - slow_vel) / accel
        settle_time = self._get_force_settle_time(smooth_time, signal_filter)
        total_duration = 2. * settle_time + accel_time
        n_samples = int(total_duration / seg_time + 0.5)
        times = np.arange(n_samples, dtype=np.float64) * seg_time
        accel_start = settle_time
        accel_end = accel_start + accel_time
        velocities = np.where(times > accel_end, fast_vel, slow_vel + np.where(
            times < accel_start, 0., accel * (times - accel_start)))
        return times, velocities

    def _simulate_force_signal(self, times, velocities, configured_pa,
                               system_pa, smooth_time, signal_filter):
        # Simulate the full PA chain and return detrended force signal
        np = self.numpy
        cmd_pos = np.cumsum(velocities) * self.SIMULATION_SEG_TIME
        smoothed_pa = self._calc_pa(times, cmd_pos, configured_pa, smooth_time)
        actual_pos = self._calc_nozzle_flow(smoothed_pa, system_pa)
        # Force in mm (without stiffness coefficient)
        force = smoothed_pa - actual_pos
        filtered_force = signal_filter.smooth(times, force)
        steady = self._estimate_steady_state(filtered_force)
        return filtered_force - steady

    def _simulate_check_overshoot(self, configured_pa, system_pa, smooth_time,
                                  signal_filter, times, velocities):
        # Simulate PA chain for a velocity step and check for overshoot
        detrended = self._simulate_force_signal(
            times, velocities, configured_pa, system_pa,
            smooth_time, signal_filter)
        # accel and force_sign values are just hardcoded conventions to
        # represent a generated acceleration direction and a maching force
        return self._check_window_overshoot(detrended, accel=1., force_sign=1.)

    def _find_corrected_pa(self, detected_pa, smooth_time, signal_filter,
                           slow_velocity, fast_velocity, accel):
        # Binary search over system_pa to find corrected PA value
        times, velocities = self._simulate_velocity_step(
            slow_velocity, fast_velocity, accel, smooth_time, signal_filter)
        def sim_test_cb(system_pa):
            res = self._simulate_check_overshoot(
                    detected_pa, system_pa, smooth_time, signal_filter,
                    times, velocities)
            if res is None:
                return False
            return res[0]
        left = detected_pa * self.PA_CORRECTED_REL_START
        right = detected_pa
        if not sim_test_cb(left) or sim_test_cb(right):
            return detected_pa
        while True:
            m = 0.5 * (left + right)
            tol = max(self.PA_SEARCH_MIN_INTERVAL,
                      self.PA_CORRECTED_REL_TOL * m)
            if right - left <= tol:
                return m
            if sim_test_cb(m):
                left = m
            else:
                right = m

    def _respond(self, gcmd, optimal_pa, corrected_pa, pa_range):
        lines = ["Search complete: %d tests, pressure_advance=%.4f" % (
            self._test_count, optimal_pa)]
        if optimal_pa < PA_WARNING_LOW:
            lines.append(
                "WARNING: tuned PA is very low (%.4f). This may indicate "
                "an issue with the sensor, hotend, or filament." % optimal_pa)
        if self.PA_WARNING_HIGH_MARGIN + optimal_pa > self.PA_SEARCH_MAX:
            lines.append(
                "WARNING: tuned PA is very high (%.4f, limit=%.2f). This "
                "may indicate an issue with the sensor, hotend, or filament."
                % (optimal_pa, self.PA_SEARCH_MAX))
        if (corrected_pa != optimal_pa
                and abs(corrected_pa - optimal_pa)
                    > self.CORRECTION_MIN_FRACTION * optimal_pa):
            lines.append("PA correction: pressure_advance=%.4f" % corrected_pa)
            lines.append("(accounts for force data filtering)")
        if pa_range is not None:
            lo, hi = pa_range
            warn_range = 0.5 * (lo + hi) * self.PA_RANGE_WARNING_FRACTION
            if hi < self.PA_SEARCH_MAX and hi - optimal_pa <= warn_range:
                lines.append(
                    "WARNING: tuned PA is too close to PA_RANGE upper "
                    "bound (%.4f). The optimal value may be higher." % hi)
            if lo > 0 and optimal_pa - lo <= warn_range:
                lines.append(
                    "WARNING: tuned PA is too close to PA_RANGE lower "
                    "bound (%.4f). The optimal value may be lower." % lo)
        gcmd.respond_info("\n".join(lines))

    def _parse_pa_range(self, gcmd):
        raw = gcmd.get('PA_RANGE', None)
        if raw is None:
            return None
        parts = [v.strip() for v in raw.split(',')]
        if len(parts) != 2:
            raise self.printer.command_error(
                "PA_RANGE requires two comma-separated values (min,max)")
        try:
            lo, hi = float(parts[0]), float(parts[1])
        except ValueError:
            raise self.printer.command_error(
                "PA_RANGE: invalid numeric value in '%s'" % raw)
        if lo < 0:
            raise self.printer.command_error(
                "PA_RANGE: min must be non-negative")
        if hi > self.PA_SEARCH_MAX:
            raise self.printer.command_error(
                "PA_RANGE: max exceeds search limit (%.2f)"
                % self.PA_SEARCH_MAX)
        if lo >= hi:
            raise self.printer.command_error(
                "PA_RANGE: min must be less than max")
        return (lo, hi)

    def _test_force_sign(self, gcmd, extruder_test, data_collector,
                         filament_diameter):
        seg_time = extruder_test.segment_time.get(gcmd)
        extruder_moves = extruder_test.gen_test(
                gcmd, filament_diameter, full_purge=True,
                test_repetitions=0, gen_force_sig_test=True)
        all_samples, all_moves = data_collector.collect(
                lambda: extruder_test.run_test(extruder_moves))
        force_sign = detect_force_sign(all_samples, all_moves)
        self._extruder_purged = True
        if not force_sign:
            raise gcmd.error("Could not detect the sensor force direction."
                            " Check sensor and filament and try again.")
        return force_sign

    def run(self, gcmd, extruder_test, data_collector, extruder_status,
            signal_filter, force_sign):
        self._extruder_purged = False
        self._test_count = 0
        filament_diameter = extruder_status['filament_diameter']
        area = math.pi * (filament_diameter / 2.) ** 2
        smooth_time = extruder_status['smooth_time']
        pa_range = self._parse_pa_range(gcmd)

        slow_velocity = extruder_test.get_slow_velocity(gcmd, filament_diameter)
        fast_velocity = extruder_test.get_fast_velocity(gcmd, filament_diameter)
        accel_limit = (fast_velocity - slow_velocity) / (
            self.ACCEL_DURATION_SMOOTH_TIME_MULT * smooth_time)
        extruder = data_collector.extruder
        old_e_accel = extruder.set_extrude_only_accel_limit(accel_limit)
        self.gcode.run_script_from_command("SET_PRESSURE_ADVANCE ADVANCE=0")

        try:
            if not force_sign:
                force_sign = self._test_force_sign(
                        gcmd, extruder_test, data_collector, filament_diameter)
                gcmd.respond_info("Detected extrude_force_sign=%d" % force_sign)
            gcmd.respond_info("Searching for optimal pressure_advance...")
            optimal_pa, results = self._search_optimal_pa(
                gcmd, extruder_test, data_collector,
                extruder_status, signal_filter, force_sign, pa_range)
        finally:
            extruder.set_extrude_only_accel_limit(old_e_accel)

        corrected_pa = self._find_corrected_pa(
                optimal_pa, smooth_time, signal_filter,
                slow_velocity, fast_velocity, accel_limit)
        self._log_search_stats(results, optimal_pa, corrected_pa, area)
        self._respond(gcmd, optimal_pa, corrected_pa, pa_range)


class PATester:
    def __init__(self, config):
        self.printer = config.get_printer()
        try:
            self.numpy = importlib.import_module('numpy')
        except ImportError:
            raise self.printer.command_error(
                "Failed to import `numpy` module, make sure it was "
                "installed via `~/klippy-env/bin/pip install`")
        use_scipy = config.getboolean('use_scipy', False)
        try:
            self.scipy = importlib.import_module('scipy') if use_scipy else None
        except ImportError:
            raise self.printer.command_error(
                "Failed to import `scipy` module, make sure it was "
                "installed via `~/klippy-env/bin/pip install`")
        self.force_sensor = config.get('force_sensor', None)
        self.extruder = config.get('extruder', None)
        self.filter_window = FloatParam(config, 'filter_window', 0.01,
                                        above=0.0, maxval=0.04)
        self.methods = {
                'step_response': StepResponseTest(self.printer, self.numpy,
                                                  self.scipy),
                'search_overshoot': SearchOvershootTest(self.printer,
                                                        self.numpy, self.scipy),
            }
        self.default_method = config.getchoice(
            'method', list(self.methods.keys()), 'step_response')
        self.extruder_test = PACalibrationExtruderTest(config)
        self.force_sign = IntParam(config, 'extrude_force_sign', 0,
                                   minval=-1, maxval=1)
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
        collector = PACalibrationDataCollector(self.printer, extruder, loadcell)
        signal_filter = TriangularWindowFilter(
                self.numpy, self.filter_window.get(gcmd))
        force_sign = self.force_sign.get(gcmd)
        method.run(gcmd, self.extruder_test, collector,
                   extruder_status, signal_filter, force_sign)


def load_config(config):
    return PATester(config)
