# Support for reading acceleration data from an adxl345 chip
#
# Copyright (C) 2020-2021  Kevin O'Connor <kevin@koconnor.net>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import logging, math, time, collections, threading, multiprocessing, os
from . import bus, motion_report, shaper_calibrate

# ADXL345 registers
REG_DEVID = 0x00
REG_BW_RATE = 0x2C
REG_POWER_CTL = 0x2D
REG_DATA_FORMAT = 0x31
REG_FIFO_CTL = 0x38
REG_MOD_READ = 0x80
REG_MOD_MULTI = 0x40

QUERY_RATES = {
    25: 0x8, 50: 0x9, 100: 0xa, 200: 0xb, 400: 0xc,
    800: 0xd, 1600: 0xe, 3200: 0xf,
}

ADXL345_DEV_ID = 0xe5
SET_FIFO_CTL = 0x90

FREEFALL_ACCEL = 9.80665 * 1000.
SCALE = 0.0039 * FREEFALL_ACCEL # 3.9mg/LSB * Earth gravity in mm/s**2

CALIBRATION_NOISE_THRESHOLD = 1e4
ACCELERATION_OUTLIER_THRESHOLD = FREEFALL_ACCEL * 5.

Accel_Measurement = collections.namedtuple(
    'Accel_Measurement', ('time', 'accel_x', 'accel_y', 'accel_z'))

# Helper class to obtain measurements
class ADXL345QueryHelper:
    def __init__(self, printer, cconn):
        self.printer = printer
        self.cconn = cconn
        print_time = printer.lookup_object('toolhead').get_last_move_time()
        self.request_start_time = self.request_end_time = print_time
        self.samples = self.raw_samples = []
    def finish_measurements(self):
        toolhead = self.printer.lookup_object('toolhead')
        self.request_end_time = toolhead.get_last_move_time()
        toolhead.wait_moves()
        self.cconn.finalize()
    def _get_raw_samples(self):
        raw_samples = self.cconn.get_messages()
        if raw_samples:
            self.raw_samples = raw_samples
        return self.raw_samples
    def has_valid_samples(self):
        raw_samples = self._get_raw_samples()
        for msg in raw_samples:
            data = msg['params']['data']
            first_sample_time = data[0][0]
            last_sample_time = data[-1][0]
            if (first_sample_time > self.request_end_time
                    or last_sample_time < self.request_start_time):
                continue
            # The time intervals [first_sample_time, last_sample_time]
            # and [request_start_time, request_end_time] have non-zero
            # intersection. It is still theoretically possible that none
            # of the samples from raw_samples fall into the time interval
            # [request_start_time, request_end_time] if it is too narrow
            # or on very heavy data losses. In practice, that interval
            # is at least 1 second, so this possibility is negligible.
            return True
        return False
    def get_samples(self):
        raw_samples = self._get_raw_samples()
        if not raw_samples:
            return self.samples
        total = sum([len(m['params']['data']) for m in raw_samples])
        count = 0
        self.samples = samples = [None] * total
        for msg in raw_samples:
            for samp_time, x, y, z in msg['params']['data']:
                if samp_time < self.request_start_time:
                    continue
                if samp_time > self.request_end_time:
                    break
                samples[count] = Accel_Measurement(samp_time, x, y, z)
                count += 1
        del samples[count:]
        return self.samples
    def get_measurements_range(self):
        return (self.request_start_time, self.request_end_time)
    def write_to_file(self, filename):
        def write_impl():
            try:
                # Try to re-nice writing process
                os.nice(20)
            except:
                pass
            f = open(filename, "w")
            f.write("#time,accel_x,accel_y,accel_z\n")
            samples = self.samples or self.get_samples()
            for t, accel_x, accel_y, accel_z in samples:
                f.write("%.6f,%.6f,%.6f,%.6f\n" % (
                    t, accel_x, accel_y, accel_z))
            f.close()
        write_proc = multiprocessing.Process(target=write_impl)
        write_proc.daemon = True
        write_proc.start()

class AccelerometerCalibrator:
    def __init__(self, printer, chip):
        self.printer = printer
        self.gcode = self.printer.lookup_object('gcode')
        self.chip = chip
        # Test parameters
        self.max_accel = 500.
        self.wait_test_sec = 5.
        self.move_test_runs = 15
        self.move_test_len = 4. # Should work well with GT-2 belts
    def _run_wait_test(self, toolhead):
        aclient = self.chip.start_internal_client()
        toolhead.dwell(self.wait_test_sec)
        aclient.finish_measurements()
        return aclient
    def _run_move_test(self, toolhead, axis_dir):
        X, Y, Z, E = toolhead.get_position()
        systime = self.printer.get_reactor().monotonic()
        toolhead_info = toolhead.get_status(systime)
        max_accel = toolhead_info['max_accel']
        max_accel_to_decel = toolhead_info['max_accel_to_decel']
        # The test runs as follows:
        # * accelerate for t_seg/2 time
        # * cruise for t_seg time
        # * decelerate for t_seg/2 time
        # * accelerate for t_seg/2 time in reverse direction
        # .....
        L = self.move_test_len
        accel = min(self.max_accel, 6. * max_accel_to_decel, max_accel)
        t_seg = math.sqrt(L / (.75 * accel))
        freq = .25 / t_seg
        max_v = .5 * t_seg * accel
        toolhead.cmd_M204(self.gcode.create_gcode_command(
            'M204', 'M204', {'S': accel}))
        nX = X + axis_dir[0] * L
        nY = Y + axis_dir[1] * L
        aclient = self.chip.start_internal_client()
        try:
            for i in range(self.move_test_runs):
                toolhead.move([nX, nY, Z, E], max_v)
                toolhead.move([X, Y, Z, E], max_v)
        finally:
            aclient.finish_measurements()
        toolhead.cmd_M204(self.gcode.create_gcode_command(
            'M204', 'M204', {'S': max_accel}))
        return (freq, accel, aclient)
    def _integrate(self, samples, t_range):
        avg = samples[(samples[:,0] >= t_range[0]) &
                      (samples[:,0] < t_range[1])].mean(axis=0)
        return avg[1:]
    def _compute_freefall_accel(self, data):
        samples = self.np.asarray(data.get_samples())
        g = self._integrate(samples, data.get_measurements_range())
        freefall_accel = self.np.linalg.norm(g)
        return freefall_accel, g
    def _compute_measured_accel(self, freq, data):
        start_t, end_t = data.get_measurements_range()
        samples = self.np.asarray(data.get_samples())
        quarter_period = .25 / freq
        a = self.np.zeros(shape=3)
        n = 0
        sign = 1
        while end_t - start_t > quarter_period:
            a += sign * self._integrate(
                    samples, (start_t, start_t + quarter_period))
            start_t += quarter_period
            a -= sign * self._integrate(
                    samples, (start_t, start_t + quarter_period))
            start_t += quarter_period
            sign = -sign
            n += 2
        # Acceleration is active only half of the time
        a *= 2. / n
        measured_acccel = self.np.linalg.norm(a)
        return measured_acccel, a
    def _calculate_gravity(self):
        self.gravity = [round(FREEFALL_ACCEL * coeff, 1)
                        for coeff in self.results['z']]
    def _calculate_axes_transform(self):
        linalg = self.np.linalg
        A = self.np.zeros(shape=(3, 3))
        if 'x' in self.results:
            A[:,0] = self.results['x']
        else:
            a_y = self.np.asarray(self.results['y'])
            a_z = self.np.asarray(self.results['z'])
            # Exact X axis direction does not matter, so
            # creating the standard right-handed coordinate system
            a_x = self.np.cross(a_y, a_z)
            # Normalize to the same scale as Y axis
            a_x *= linalg.norm(a_y) / linalg.norm(a_x)
            A[:,0] = a_x
        if 'y' in self.results:
            A[:,1] = self.results['y']
        else:
            a_x = self.np.asarray(self.results['x'])
            a_z = self.np.asarray(self.results['z'])
            # Exact Y axis direction does not matter, so
            # creating the standard right-handed coordinate system
            a_y = self.np.cross(a_z, a_x)
            # Normalize to the same scale as X axis
            a_y *= linalg.norm(a_x) / linalg.norm(a_y)
            A[:,1] = a_y
        A[:,2] = self.results['z']
        self.axes_transform = linalg.inv(A)
    def _get_chip_name(self):
        return self.chip.get_config().get_name()
    def _save_gravity(self, gcmd):
        chip_name = self._get_chip_name()
        str_val = ','.join(['%.1f' % (coeff,) for coeff in self.gravity])
        configfile = self.printer.lookup_object('configfile')
        configfile.set(chip_name, 'gravity', str_val)
        self.chip.set_transform(self.axes_transform, self.gravity)
        gcmd.respond_info(
                'SAVE_CONFIG command will update %s configuration with '
                'gravity = %s parameter' % (chip_name, str_val))
    def _save_axes_transform(self, gcmd):
        chip_name = self._get_chip_name()
        configfile = self.printer.lookup_object('configfile')
        if self.chip.get_config().get('axes_map', None, note_valid=False):
            configfile.set(chip_name, 'axes_map', '')
        str_val = '\n'.join([','.join(['%.9f' % (coeff,) for coeff in axis])
                             for axis in self.axes_transform])
        configfile.set(chip_name, 'axes_transform', '\n' + str_val)
        self.chip.set_transform(self.axes_transform, self.gravity)
        gcmd.respond_info(
                'SAVE_CONFIG command will also update %s configuration with '
                'axes_transform =\n%s' % (chip_name, str_val))
    def _calibrate_z_axis(self, toolhead, gcmd, output):
        gcmd.respond_info('Calibrating z axis')
        data = self._run_wait_test(toolhead)
        if not data.has_valid_samples():
            raise gcmd.error('No accelerometer measurements found')
        if output is not None:
            filename = "/tmp/%s-z.csv" % (output,)
            data.write_to_file(filename)
        for _, ax, ay, az in data.get_samples():
            if max(ax, ay, az) > ACCELERATION_OUTLIER_THRESHOLD:
                gcmd.respond_info(
                        'WARNING: large acceleration reading detected (%.1f '
                        'mm/sec^2). This generally indicates communication '
                        'errors with the accelerometer (e.g. because of '
                        'electromagnetic interference on the cables or '
                        'problems with the sensor). This may severerly impact '
                        'the calibration and resonance testing results.' % (
                            max(ax, ay, az),))
                break
        helper = shaper_calibrate.ShaperCalibrate(self.printer)
        processed_data = helper.process_accelerometer_data(data)
        psd = processed_data.get_psd()
        max_noise_ind = psd.argmax()
        if psd[max_noise_ind] > CALIBRATION_NOISE_THRESHOLD:
            gcmd.respond_info(
                    'WARNING: strong periodic noise detected at %.1f Hz. This '
                    'could be a loud unbalanced fan or some other devices '
                    'working nearby. Please turn off fans (e.g. hotend fan) '
                    'and other sources of vibrations for accelerometer '
                    'calibration and resonance testing for best results.' % (
                        processed_data.get_freq_bins()[max_noise_ind],))
        self.np = processed_data.get_numpy()
        freefall_accel, g = self._compute_freefall_accel(data)
        if abs(freefall_accel - FREEFALL_ACCEL) > .2 * FREEFALL_ACCEL:
            chip_name = self._get_chip_name()
            raise gcmd.error('%s is defunct: measured invalid freefall accel '
                             '%.3f (mm/sec^2) vs ~ %.3f (mm/sec^2)' % (
                                 chip_name, freefall_accel, FREEFALL_ACCEL))
        self.results['z'] = g / FREEFALL_ACCEL
        self._calculate_gravity()
        gcmd.respond_info(
                'Detected gravity direction: %s' % (', '.join(
                    ['%.6f' % (val,) for val in self.results['z']])))
        self._save_gravity(gcmd)
    def _calibrate_xy_axis(self, axis, axis_dir, toolhead, gcmd, output):
        gcmd.respond_info('Calibrating %s axis' % (axis,))
        chip_name = self._get_chip_name()
        freq, accel, data = self._run_move_test(toolhead, axis_dir)
        if not data.has_valid_samples():
            raise gcmd.error('No accelerometer measurements found')
        if output is not None:
            filename = "/tmp/%s-%s.csv" % (output, axis)
            data.write_to_file(filename)
        measured_acccel, a = self._compute_measured_accel(freq, data)
        if measured_acccel > .2 * accel:
            if abs(measured_acccel - accel) > .2 * accel:
                raise gcmd.error(
                        '%s measured spurious acceleration on %s axis: '
                        '%.3f vs %.3f (mm/sec^2)' % (chip_name, axis,
                                                     measured_acccel, accel))
            self.results[axis] = a / accel
            gcmd.respond_info(
                    'Detected %s direction: %s' % (axis, ', '.join(
                        ['%.6f' % (val,) for val in self.results[axis]])))
        else:
            gcmd.respond_info('%s is not kinematically connected to the '
                              'movement of %s axis' % (chip_name, axis))
    def calibrate(self, gcmd, output=None):
        toolhead = self.printer.lookup_object('toolhead')
        reactor = self.printer.get_reactor()
        # Reset adxl345 transformations
        self.axes_transform = [[1., 0., 0.],
                               [0., 1., 0.],
                               [0., 0., 1.]]
        self.gravity = [0., 0., 0.]
        self.chip.set_transform(self.axes_transform, self.gravity)
        self.results = {}
        self._calibrate_z_axis(toolhead, gcmd, output)
        reactor.pause(reactor.monotonic() + 0.1)
        self._calibrate_xy_axis('x', (1., 0.), toolhead, gcmd, output)
        reactor.pause(reactor.monotonic() + 0.1)
        self._calibrate_xy_axis('y', (0., 1.), toolhead, gcmd, output)
        reactor.pause(reactor.monotonic() + 0.1)
        if 'x' not in self.results and 'y' not in self.results:
            raise gcmd.error(
                    '%s is not kinematically connected to either of X or '
                    'Y printer axis, impossible to calibrate automatically. '
                    'Please manually set axes_map parameter.' % (
                        self._get_chip_name(),))
        gcmd.respond_info('Computing axes transform')
        self._calculate_axes_transform()
        self._save_axes_transform(gcmd)

# Helper class for G-Code commands
class ADXLCommandHelper:
    def __init__(self, config, chip):
        self.printer = config.get_printer()
        self.chip = chip
        self.bg_client = None
        self.name = config.get_name().split()[-1]
        self.register_commands(self.name)
        if self.name == "adxl345":
            self.register_commands(None)
    def register_commands(self, name):
        # Register commands
        gcode = self.printer.lookup_object('gcode')
        gcode.register_mux_command("ACCELEROMETER_MEASURE", "CHIP", name,
                                   self.cmd_ACCELEROMETER_MEASURE,
                                   desc=self.cmd_ACCELEROMETER_MEASURE_help)
        gcode.register_mux_command("ACCELEROMETER_QUERY", "CHIP", name,
                                   self.cmd_ACCELEROMETER_QUERY,
                                   desc=self.cmd_ACCELEROMETER_QUERY_help)
        gcode.register_mux_command("ACCELEROMETER_CALIBRATE", "CHIP", name,
                                   self.cmd_ACCELEROMETER_CALIBRATE,
                                   desc=self.cmd_ACCELEROMETER_CALIBRATE_help)
        gcode.register_mux_command("ACCELEROMETER_DEBUG_READ", "CHIP", name,
                                   self.cmd_ACCELEROMETER_DEBUG_READ,
                                   desc=self.cmd_ACCELEROMETER_DEBUG_READ_help)
        gcode.register_mux_command("ACCELEROMETER_DEBUG_WRITE", "CHIP", name,
                                   self.cmd_ACCELEROMETER_DEBUG_WRITE,
                                   desc=self.cmd_ACCELEROMETER_DEBUG_WRITE_help)
    cmd_ACCELEROMETER_MEASURE_help = "Start/stop accelerometer"
    def cmd_ACCELEROMETER_MEASURE(self, gcmd):
        if self.bg_client is None:
            # Start measurements
            self.bg_client = self.chip.start_internal_client()
            gcmd.respond_info("adxl345 measurements started")
            return
        # End measurements
        name = gcmd.get("NAME", time.strftime("%Y%m%d_%H%M%S"))
        if not name.replace('-', '').replace('_', '').isalnum():
            raise gcmd.error("Invalid adxl345 NAME parameter")
        bg_client = self.bg_client
        self.bg_client = None
        bg_client.finish_measurements()
        # Write data to file
        if self.name == "adxl345":
            filename = "/tmp/adxl345-%s.csv" % (name,)
        else:
            filename = "/tmp/adxl345-%s-%s.csv" % (self.name, name,)
        bg_client.write_to_file(filename)
        gcmd.respond_info("Writing raw accelerometer data to %s file"
                          % (filename,))
    cmd_ACCELEROMETER_QUERY_help = "Query accelerometer for the current values"
    def cmd_ACCELEROMETER_QUERY(self, gcmd):
        aclient = self.chip.start_internal_client()
        self.printer.lookup_object('toolhead').dwell(1.)
        aclient.finish_measurements()
        values = aclient.get_samples()
        if not values:
            raise gcmd.error("No adxl345 measurements found")
        _, accel_x, accel_y, accel_z = values[-1]
        gcmd.respond_info("adxl345 values (x, y, z): %.6f, %.6f, %.6f"
                          % (accel_x, accel_y, accel_z))
    cmd_ACCELEROMETER_CALIBRATE_help = "Automatically calibrate accelerometer"
    def cmd_ACCELEROMETER_CALIBRATE(self, gcmd):
        output = gcmd.get("OUTPUT", None)
        if not output.replace('-', '').replace('_', '').isalnum():
            raise gcmd.error("Invalid OUTPUT parameter")
        AccelerometerCalibrator(self.printer, self.chip).calibrate(gcmd, output)
    cmd_ACCELEROMETER_DEBUG_READ_help = "Query adxl345 register (for debugging)"
    def cmd_ACCELEROMETER_DEBUG_READ(self, gcmd):
        reg = gcmd.get("REG", minval=29, maxval=57, parser=lambda x: int(x, 0))
        val = self.chip.read_reg(reg)
        gcmd.respond_info("ADXL345 REG[0x%x] = 0x%x" % (reg, val))
    cmd_ACCELEROMETER_DEBUG_WRITE_help = "Set adxl345 register (for debugging)"
    def cmd_ACCELEROMETER_DEBUG_WRITE(self, gcmd):
        reg = gcmd.get("REG", minval=29, maxval=57, parser=lambda x: int(x, 0))
        val = gcmd.get("VAL", minval=0, maxval=255, parser=lambda x: int(x, 0))
        self.chip.set_reg(reg, val)

# Helper class for chip clock synchronization via linear regression
class ClockSyncRegression:
    def __init__(self, mcu, chip_clock_smooth, decay = 1. / 20.):
        self.mcu = mcu
        self.chip_clock_smooth = chip_clock_smooth
        self.decay = decay
        self.last_chip_clock = self.last_exp_mcu_clock = 0.
        self.mcu_clock_avg = self.mcu_clock_variance = 0.
        self.chip_clock_avg = self.chip_clock_covariance = 0.
    def reset(self, mcu_clock, chip_clock):
        self.mcu_clock_avg = self.last_mcu_clock = mcu_clock
        self.chip_clock_avg = chip_clock
        self.mcu_clock_variance = self.chip_clock_covariance = 0.
        self.last_chip_clock = self.last_exp_mcu_clock = 0.
    def update(self, mcu_clock, chip_clock):
        # Update linear regression
        decay = self.decay
        diff_mcu_clock = mcu_clock - self.mcu_clock_avg
        self.mcu_clock_avg += decay * diff_mcu_clock
        self.mcu_clock_variance = (1. - decay) * (
            self.mcu_clock_variance + diff_mcu_clock**2 * decay)
        diff_chip_clock = chip_clock - self.chip_clock_avg
        self.chip_clock_avg += decay * diff_chip_clock
        self.chip_clock_covariance = (1. - decay) * (
            self.chip_clock_covariance + diff_mcu_clock*diff_chip_clock*decay)
    def set_last_chip_clock(self, chip_clock):
        base_mcu, base_chip, inv_cfreq = self.get_clock_translation()
        self.last_chip_clock = chip_clock
        self.last_exp_mcu_clock = base_mcu + (chip_clock-base_chip) * inv_cfreq
    def get_clock_translation(self):
        inv_chip_freq = self.mcu_clock_variance / self.chip_clock_covariance
        if not self.last_chip_clock:
            return self.mcu_clock_avg, self.chip_clock_avg, inv_chip_freq
        # Find mcu clock associated with future chip_clock
        s_chip_clock = self.last_chip_clock + self.chip_clock_smooth
        scdiff = s_chip_clock - self.chip_clock_avg
        s_mcu_clock = self.mcu_clock_avg + scdiff * inv_chip_freq
        # Calculate frequency to converge at future point
        mdiff = s_mcu_clock - self.last_exp_mcu_clock
        s_inv_chip_freq = mdiff / self.chip_clock_smooth
        return self.last_exp_mcu_clock, self.last_chip_clock, s_inv_chip_freq
    def get_time_translation(self):
        base_mcu, base_chip, inv_cfreq = self.get_clock_translation()
        clock_to_print_time = self.mcu.clock_to_print_time
        base_time = clock_to_print_time(base_mcu)
        inv_freq = clock_to_print_time(base_mcu + inv_cfreq) - base_time
        return base_time, base_chip, inv_freq

MIN_MSG_TIME = 0.100

BYTES_PER_SAMPLE = 5
SAMPLES_PER_BLOCK = 10

# Printer class that controls ADXL345 chip
class ADXL345:
    def __init__(self, config):
        self.config = config
        self.printer = config.get_printer()
        ADXLCommandHelper(config, self)
        gravity = config.getfloatlist('gravity', (0., 0., 0.), count=3)
        axes_transform = config.getlists(
                'axes_transform', None, seps=(',', '\n'), parser=float, count=3)
        if axes_transform is None:
            axes_map = config.getlist('axes_map', ('x','y','z'), count=3)
            am = {'x': (0, 1.), 'y': (1, 1.), 'z': (2, 1.),
                  '-x': (0, -1.), '-y': (1, -1.), '-z': (2, -1.)}
            if any([a not in am for a in axes_map]):
                raise config.error("Invalid adxl345 axes_map parameter")
            axes_transform = [[0., 0., 0.],
                              [0., 0., 0.],
                              [0., 0., 0.]]
            for i, a in enumerate(axes_map):
                ind, val = am[a.strip()]
                axes_transform[i][ind] = val
        self.set_transform(axes_transform, gravity)
        self.data_rate = config.getint('rate', 3200)
        if self.data_rate not in QUERY_RATES:
            raise config.error("Invalid rate parameter: %d" % (self.data_rate,))
        # A lock to guard start/stop process and make it atomic
        self.start_stop_lock = threading.Lock()
        self.query_rate = 0
        # Measurement storage (accessed from background thread)
        self.lock = threading.Lock()
        self.raw_samples = []
        # Setup mcu sensor_adxl345 bulk query code
        self.spi = bus.MCU_SPI_from_config(config, 3, default_speed=5000000)
        self.mcu = mcu = self.spi.get_mcu()
        self.oid = oid = mcu.create_oid()
        self.query_adxl345_cmd = self.query_adxl345_end_cmd = None
        self.query_adxl345_status_cmd = None
        mcu.add_config_cmd("config_adxl345 oid=%d spi_oid=%d"
                           % (oid, self.spi.get_oid()))
        mcu.add_config_cmd("query_adxl345 oid=%d clock=0 rest_ticks=0"
                           % (oid,), on_restart=True)
        mcu.register_config_callback(self._build_config)
        mcu.register_response(self._handle_adxl345_data, "adxl345_data", oid)
        # Clock tracking
        self.last_sequence = self.max_query_duration = 0
        self.last_limit_count = self.last_error_count = 0
        self.clock_sync = ClockSyncRegression(self.mcu, 640)
        # API server endpoints
        self.api_dump = motion_report.APIDumpHelper(
            self.printer, self._api_update, self._api_startstop, 0.100)
        self.name = config.get_name().split()[-1]
        wh = self.printer.lookup_object('webhooks')
        wh.register_mux_endpoint("adxl345/dump_adxl345", "sensor", self.name,
                                 self._handle_dump_adxl345)
    def _build_config(self):
        cmdqueue = self.spi.get_command_queue()
        self.query_adxl345_cmd = self.mcu.lookup_command(
            "query_adxl345 oid=%c clock=%u rest_ticks=%u", cq=cmdqueue)
        self.query_adxl345_end_cmd = self.mcu.lookup_query_command(
            "query_adxl345 oid=%c clock=%u rest_ticks=%u",
            "adxl345_status oid=%c clock=%u query_ticks=%u next_sequence=%hu"
            " buffered=%c fifo=%c limit_count=%hu", oid=self.oid, cq=cmdqueue)
        self.query_adxl345_status_cmd = self.mcu.lookup_query_command(
            "query_adxl345_status oid=%c",
            "adxl345_status oid=%c clock=%u query_ticks=%u next_sequence=%hu"
            " buffered=%c fifo=%c limit_count=%hu", oid=self.oid, cq=cmdqueue)
    def read_reg(self, reg):
        params = self.spi.spi_transfer([reg | REG_MOD_READ, 0x00])
        response = bytearray(params['response'])
        return response[1]
    def set_reg(self, reg, val, minclock=0):
        self.spi.spi_send([reg, val & 0xFF], minclock=minclock)
        stored_val = self.read_reg(reg)
        if stored_val != val:
            raise self.printer.command_error(
                    "Failed to set ADXL345 register [0x%x] to 0x%x: got 0x%x. "
                    "This is generally indicative of connection problems "
                    "(e.g. faulty wiring) or a faulty adxl345 chip." % (
                        reg, val, stored_val))
    def get_config(self):
        return self.config
    def set_transform(self, axes_transform, gravity):
        self.gravity = [coeff / SCALE for coeff in gravity]
        self.axes_transform = [[coeff * SCALE for coeff in axis]
                               for axis in axes_transform]
    # Measurement collection
    def is_measuring(self):
        return self.query_rate > 0
    def _handle_adxl345_data(self, params):
        with self.lock:
            self.raw_samples.append(params)
    def _extract_samples(self, raw_samples):
        # Load variables to optimize inner loop below
        tr_x, tr_y, tr_z = self.axes_transform
        g_x, g_y, g_z = self.gravity
        last_sequence = self.last_sequence
        time_base, chip_base, inv_freq = self.clock_sync.get_time_translation()
        # Process every message in raw_samples
        count = seq = 0
        samples = [None] * (len(raw_samples) * SAMPLES_PER_BLOCK)
        for params in raw_samples:
            seq_diff = (last_sequence - params['sequence']) & 0xffff
            seq_diff -= (seq_diff & 0x8000) << 1
            seq = last_sequence - seq_diff
            d = bytearray(params['data'])
            msg_cdiff = seq * SAMPLES_PER_BLOCK - chip_base
            for i in range(len(d) // BYTES_PER_SAMPLE):
                d_xyz = d[i*BYTES_PER_SAMPLE:(i+1)*BYTES_PER_SAMPLE]
                xlow, ylow, zlow, xzhigh, yzhigh = d_xyz
                if yzhigh & 0x80:
                    self.last_error_count += 1
                    continue
                rx = ((xlow | ((xzhigh & 0x1f) << 8)) - ((xzhigh & 0x10) << 9)
                        - g_x)
                ry = ((ylow | ((yzhigh & 0x1f) << 8)) - ((yzhigh & 0x10) << 9)
                        - g_y)
                rz = ((zlow | ((xzhigh & 0xe0) << 3) | ((yzhigh & 0xe0) << 6))
                      - ((yzhigh & 0x40) << 7)) - g_z
                x = round(tr_x[0] * rx + tr_x[1] * ry + tr_x[2] * rz, 6)
                y = round(tr_y[0] * rx + tr_y[1] * ry + tr_y[2] * rz, 6)
                z = round(tr_z[0] * rx + tr_z[1] * ry + tr_z[2] * rz, 6)
                ptime = round(time_base + (msg_cdiff + i) * inv_freq, 6)
                samples[count] = (ptime, x, y, z)
                count += 1
        self.clock_sync.set_last_chip_clock(seq * SAMPLES_PER_BLOCK + i)
        del samples[count:]
        return samples
    def _update_clock(self, minclock=0):
        # Query current state
        for retry in range(5):
            params = self.query_adxl345_status_cmd.send([self.oid],
                                                        minclock=minclock)
            fifo = params['fifo'] & 0x7f
            if fifo <= 32:
                break
        else:
            raise self.printer.command_error("Unable to query adxl345 fifo")
        mcu_clock = self.mcu.clock32_to_clock64(params['clock'])
        sequence = (self.last_sequence & ~0xffff) | params['next_sequence']
        if sequence < self.last_sequence:
            sequence += 0x10000
        self.last_sequence = sequence
        buffered = params['buffered']
        limit_count = (self.last_limit_count & ~0xffff) | params['limit_count']
        if limit_count < self.last_limit_count:
            limit_count += 0x10000
        self.last_limit_count = limit_count
        duration = params['query_ticks']
        if duration > self.max_query_duration:
            # Skip measurement as a high query time could skew clock tracking
            self.max_query_duration = max(2 * self.max_query_duration,
                                          self.mcu.seconds_to_clock(.000005))
            return
        self.max_query_duration = 2 * duration
        msg_count = (sequence * SAMPLES_PER_BLOCK
                     + buffered // BYTES_PER_SAMPLE + fifo)
        # The "chip clock" is the message counter plus .5 for average
        # inaccuracy of query responses and plus .5 for assumed offset
        # of adxl345 hw processing time.
        chip_clock = msg_count + 1
        self.clock_sync.update(mcu_clock + duration // 2, chip_clock)
    def _start_measurements(self):
        if self.is_measuring():
            return
        # In case of miswiring, testing ADXL345 device ID prevents treating
        # noise or wrong signal as a correctly initialized device
        dev_id = self.read_reg(REG_DEVID)
        if dev_id != ADXL345_DEV_ID:
            raise self.printer.command_error(
                "Invalid adxl345 id (got %x vs %x).\n"
                "This is generally indicative of connection problems\n"
                "(e.g. faulty wiring) or a faulty adxl345 chip."
                % (dev_id, ADXL345_DEV_ID))
        # Setup chip in requested query rate
        self.set_reg(REG_POWER_CTL, 0x00)
        self.set_reg(REG_DATA_FORMAT, 0x0B)
        self.set_reg(REG_FIFO_CTL, 0x00)
        self.set_reg(REG_BW_RATE, QUERY_RATES[self.data_rate])
        self.set_reg(REG_FIFO_CTL, SET_FIFO_CTL)
        # Setup samples
        with self.lock:
            self.raw_samples = []
        # Start bulk reading
        systime = self.printer.get_reactor().monotonic()
        print_time = self.mcu.estimated_print_time(systime) + MIN_MSG_TIME
        reqclock = self.mcu.print_time_to_clock(print_time)
        rest_ticks = self.mcu.seconds_to_clock(4. / self.data_rate)
        self.query_rate = self.data_rate
        self.query_adxl345_cmd.send([self.oid, reqclock, rest_ticks],
                                    reqclock=reqclock)
        logging.info("ADXL345 starting '%s' measurements", self.name)
        # Initialize clock tracking
        self.last_sequence = 0
        self.last_limit_count = self.last_error_count = 0
        self.clock_sync.reset(reqclock, 0)
        self.max_query_duration = 1 << 31
        self._update_clock(minclock=reqclock)
        self.max_query_duration = 1 << 31
    def _finish_measurements(self):
        if not self.is_measuring():
            return
        # Halt bulk reading
        params = self.query_adxl345_end_cmd.send([self.oid, 0, 0])
        self.query_rate = 0
        with self.lock:
            self.raw_samples = []
        logging.info("ADXL345 finished '%s' measurements", self.name)
    # API interface
    def _api_update(self, eventtime):
        self._update_clock()
        with self.lock:
            raw_samples = self.raw_samples
            self.raw_samples = []
        if not raw_samples:
            return {}
        samples = self._extract_samples(raw_samples)
        if not samples:
            return {}
        return {'data': samples, 'errors': self.last_error_count,
                'overflows': self.last_limit_count}
    def _api_startstop(self, is_start):
        reactor = self.printer.get_reactor()
        # Do not block, especially the main Klipper thread, but instead give way
        # to other threads to prevent deadlocks
        while not self.start_stop_lock.acquire(False):
            reactor.pause(reactor.monotonic() + 0.1)
        try:
            if is_start:
                self._start_measurements()
            else:
                self._finish_measurements()
        finally:
            self.start_stop_lock.release()
    def _handle_dump_adxl345(self, web_request):
        self.api_dump.add_client(web_request)
        hdr = ('time', 'x_acceleration', 'y_acceleration', 'z_acceleration')
        web_request.send({'header': hdr})
    def start_internal_client(self):
        cconn = self.api_dump.add_internal_client()
        return ADXL345QueryHelper(self.printer, cconn)

def load_config(config):
    return ADXL345(config)

def load_config_prefix(config):
    return ADXL345(config)
