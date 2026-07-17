# Pressure Advance calibration and tuning module
#
# Copyright (C) 2026  Dmitry Butyugin <dmbutyugin@google.com>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

import collections, math
import chelper

MoveData = collections.namedtuple(
    'MoveData', ['print_time', 'move_t', 'start_v', 'accel'])



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
        extruder = self._resolve_extruder(gcmd)
        loadcell, loadcell_name = self._resolve_loadcell(gcmd)
        systime = self.printer.get_reactor().monotonic()
        extruder_status = extruder.get_status(systime)
        if not extruder_status['can_extrude']:
            raise self.printer.command_error(
                "Extruder cannot extrude\n"
                "See the 'min_extrude_temp' config option for details")
        filament_diameter = extruder_status['filament_diameter']
        segments = self.extruder_test.gen_test(filament_diameter, gcmd)
        collector = PACalibrationDataCollector(self.printer, extruder, loadcell)
        force_samples, trapq_moves = collector.collect(
                run_test_cb=lambda: self.extruder_test.run_test(segments))
        gcmd.respond_info(
            "PA_CALIBRATE: extruder=%s loadcell=%s "
            "samples=%d moves=%d",
            extruder.get_name(), loadcell_name,
            len(force_samples), len(trapq_moves))


def load_config(config):
    return PATester(config)
