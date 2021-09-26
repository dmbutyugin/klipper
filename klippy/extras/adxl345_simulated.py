# Generate acceleration data by simulating adxl345 chip
#
# Copyright (C) 2020  Dmitry Butyugin <dmbutyugin@google.com>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import logging, time, multiprocessing, os
import chelper
from . import adxl345, motion_report, shaper_calibrate
from adxl345 import Accel_Measurement, ADXLCommandHelper
from adxl345 import write_accelerometer_data

NEVER_TIME = 9999999999999999.
UPDATE_INTERVAL = 0.1

# Results of the simulated measurements
class ADXL345SimulatedQueryHelper:
    def __init__(self, printer, cconn, data_rate):
        self.printer = printer
        self.cconn = cconn
        print_time = printer.lookup_object('toolhead').get_last_move_time()
        self.request_start_time = self.request_end_time = print_time
        self.samples = []
        self.data_rate = data_rate
        # Capture input shaping parameters
        self.shaper_x = self.shaper_y = self._unity_shaper()
        input_shaper = self.printer.lookup_object('input_shaper', None)
        if input_shaper is not None:
            shapers = input_shaper.get_shapers()
            for shaper in shapers:
                n, A, T = self._scale_and_shift(shaper.get_shaper())
                if not n:
                    continue
                if shaper.get_name() == 'shaper_x':
                    self.shaper_x = (n, A, T)
                elif shaper.get_name() == 'shaper_y':
                    self.shaper_y = (n, A, T)
    def finish_measurements(self):
        reactor = self.printer.get_reactor()
        toolhead = self.printer.lookup_object('toolhead')
        self.request_end_time = toolhead.get_last_move_time()
        toolhead.wait_moves()
        reactor.pause(reactor.monotonic() + UPDATE_INTERVAL)
        self.cconn.finalize()
        reactor.pause(reactor.monotonic() + UPDATE_INTERVAL)
    def _unity_shaper(self):
        return 1, [1.], [0.]
    def _scale_and_shift(self, shaper_vals):
        n, A, T = shaper_vals
        if not n:
            return shaper_vals
        inv_D = 1. / sum(A)
        # Calculate the input shaper shift
        ts = sum([A[i] * T[i] for i in range(n)]) * inv_D
        for i in range(n):
            A[i] *= inv_D
            T[i] -= ts
        return shaper_vals
    def _gen_timestamps(self):
        time = self.request_start_time
        time_per_sample = 1. / self.data_rate
        while time <= self.request_end_time:
            yield time
            time += time_per_sample
    def _gen_accel_bounds(self):
        next_time = NEVER_TIME
        accel_bounds = [Accel_Measurement(0., 0., 0., 0.)]
        for msg in self.raw_moves:
            for move in msg['params']['data']:
                move_time = move[0]
                if accel_bounds and accel_bounds[-1].time > move_time:
                    raise self.printer.command_error(
                            "Internal error: moves captured backwards")
                if next_time < move_time:
                    accel_bounds.append(
                            Accel_Measurement(next_time, 0., 0., 0.))
                ax = move[2] * move[3][0]
                ay = move[2] * move[3][1]
                az = move[2] * move[3][2]
                accel_bounds.append(Accel_Measurement(move_time, ax, ay, az))
                next_time = move_time + move[1]
        accel_bounds.append(Accel_Measurement(next_time, 0., 0., 0.))
        accel_bounds.append(Accel_Measurement(NEVER_TIME, 0., 0., 0.))
        return accel_bounds
    def _gen_samples(self, timestamps=None):
        if timestamps is None:
            timestamps = self._gen_timestamps()
        accel_bounds = self._gen_accel_bounds()
        nx, Ax, Tx = self.shaper_x
        ny, Ay, Ty = self.shaper_y
        nz, Az, Tz = self._unity_shaper()
        n = [nx, ny, nz]
        A = [Ax, Ay, Az]
        T = [Tx, Ty, Tz]
        idx = [[0] * n[0], [0] * n[1], [0] * n[2]]
        prev_time = None
        for time in timestamps:
            if prev_time and prev_time > time:
                raise self.printer.command_error(
                        "Internal error: timestamps must be increasing")
            accel = [None] * 3
            for j in range(3):
                a = A[j]
                t = T[j]
                ind = idx[j]
                res = 0.
                for k in range(n[j]):
                    ts = time - t[k]
                    i = ind[k]
                    while ts >= accel_bounds[i+1].time:
                        i += 1
                    # 0-th Accel_Measurement element is time
                    res += a[k] * accel_bounds[i][j + 1]
                    ind[k] = i
                accel[j] = res
            yield Accel_Measurement(time, accel[0], accel[1], accel[2])
            prev_time = time
    def get_samples(self):
        raw_moves = self.cconn.get_messages()
        if not raw_moves:
            return self.samples
        self.raw_moves = raw_moves
        samples = []
        for sample in self._gen_samples():
            samples.append(sample)
        self.samples = samples
        return samples
    def generate_samples(self, timestamps):
        raw_moves = self.cconn.get_messages()
        if raw_moves:
            self.raw_moves = raw_moves
        samples = []
        for sample in self._gen_samples(timestamps):
            samples.append(sample)
        return samples
    def write_to_file(self, filename):
        write_accelerometer_data(self.get_samples(), filename)

class AccelerometerDataDiff:
    def __init__(self, raw_data, simulated_data):
        self.raw_data = raw_data
        self.simulated_data = simulated_data
        self.samples = None
    def get_samples(self):
        if self.samples is not None:
            return self.samples
        raw_samples = self.raw_data.get_samples()
        timestamps = [rs.time for rs in raw_samples]
        expected_samples = self.simulated_data.generate_samples(timestamps)
        adjusted_samples = []
        for rs, es in zip(raw_samples, expected_samples):
            adjusted_samples.append(Accel_Measurement(rs.time,
                                                      rs.accel_x - es.accel_x,
                                                      rs.accel_y - es.accel_y,
                                                      rs.accel_z - es.accel_z))
        self.samples = adjusted_samples
        return self.samples
    def write_to_file(self, filename):
        write_accelerometer_data(self.get_samples(), filename)

# Printer class that controls measurments
class ADXL345Simulated:
    def __init__(self, config, printer=None):
        self.printer = config.get_printer() if printer is None else printer
        self.name = "simulated"
        if config:
            ADXLCommandHelper(config, self)
            self.data_rate = config.getint('rate', 3200, minval=1)
            if len(config.get_name().split()) > 1:
                self.name = config.get_name().split()[-1]
        else:
            self.data_rate = 3200
        # API server endpoints
        self.api_dump = motion_report.APIDumpHelper(
                self.printer, self._api_update,
                self._api_startstop, UPDATE_INTERVAL)
        self.last_api_msg = (0., 0.)
    def _start_measurements(self):
        toolhead = self.printer.lookup_object('toolhead')
        self.trapq = toolhead.get_trapq()
        logging.info("Simulated ADXL345 started measurements")
    def _finish_measurements(self):
        logging.info("Simulated ADXL345 finished measurements")
    def extract_trapq(self, start_time, end_time):
        ffi_main, ffi_lib = chelper.get_ffi()
        res = []
        while 1:
            data = ffi_main.new('struct pull_move[128]')
            count = ffi_lib.trapq_extract_old(self.trapq, data, len(data),
                                              start_time, end_time)
            if not count:
                break
            res.append((data, count))
            if count < len(data):
                break
            end_time = data[count-1].print_time
        res.reverse()
        return ([d[i] for d, cnt in res for i in range(cnt-1, -1, -1)], res)
    # API interface
    def _api_update(self, eventtime):
        qtime = self.last_api_msg[0] + min(self.last_api_msg[1], 0.100)
        data, cdata = self.extract_trapq(qtime, NEVER_TIME)
        d = [(m.print_time, m.move_t, m.accel, (m.x_r, m.y_r, m.z_r))
             for m in data]
        if d and d[0] == self.last_api_msg:
            d.pop(0)
        if not d:
            return {}
        self.last_api_msg = d[-1]
        return {"data": d}
    def _api_startstop(self, is_start):
        if is_start:
            self._start_measurements()
        else:
            self._finish_measurements()
    def start_internal_client(self):
        cconn = self.api_dump.add_internal_client()
        return ADXL345SimulatedQueryHelper(self.printer, cconn, self.data_rate)

def load_config(config):
    return ADXL345Simulated(config)
def load_config_prefix(config):
    return ADXL345Simulated(config)
