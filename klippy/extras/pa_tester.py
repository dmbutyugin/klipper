# Pressure Advance calibration and tuning module
#
# Copyright (C) 2026  Dmitry Butyugin <dmbutyugin@google.com>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

import math


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

    def run_test(self, segments, gcmd):
        toolhead = self.printer.lookup_object('toolhead')
        tpos = toolhead.get_position()
        X, Y, Z = tpos[:3]
        E = tpos[3]
        for dist, velocity in segments:
            E += dist
            toolhead.move([X, Y, Z, E], velocity)
        toolhead.wait_moves()


class PATester:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.force_sensor = config.get('force_sensor', None)
        self.extruder = config.get('extruder', None)
        self.extruder_test = PACalibrationExtruderTest(config)
        self.gcode = self.printer.lookup_object('gcode')
        self.gcode.register_command(
                "PA_CALIBRATE", self.cmd_PA_CALIBRATE,
                desc=self.cmd_PA_CALIBRATE_help)

    cmd_PA_CALIBRATE_help = "Calibrate Pressure Advance using loadcell data"
    def cmd_PA_CALIBRATE(self, gcmd):
        gcmd.respond_info("PA_CALIBRATE: force sensor integration pending")


def load_config(config):
    return PATester(config)
