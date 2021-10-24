# Kinematic input shaper to minimize motion vibrations in XY plane
#
# Copyright (C) 2019-2020  Kevin O'Connor <kevin@koconnor.net>
# Copyright (C) 2020  Dmitry Butyugin <dmbutyugin@google.com>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import collections
import chelper
from . import shaper_defs

CUSTOM_SHAPER = 'custom'

def parse_custom_shaper(custom_a_str, custom_t_str, parse_error):
    def parse_str(s):
        res = []
        for line in s.split('\n'):
            for coeff in line.split(','):
                res.append(float(coeff.strip()))
        return res
    try:
        A = parse_str(custom_a_str)
    except:
        raise parse_error("Invalid shaper A string: '%s'" % (custom_a_str,))
    if min([abs(a) for a in A]) < 0.001:
        raise parse_error("All shaper A coefficients must be non-zero")
    if sum(A) < 0.001:
        raise parse_error("Shaper A parameter must sum up to a positive number")
    try:
        T = parse_str(custom_t_str)
    except:
        raise parse_error("Invalid shaper T string: '%s'" % (custom_t_str,))
    if T != sorted(T):
        raise parse_error("Shaper T parameter is not ordered: %s" % (T,))
    if len(A) != len(T):
        raise parse_error("Shaper A and T parameters must have the same length:"
                          " %d vs %d" % (len(A), len(T),))
    dur = T[-1] - T[0]
    if len(T) > 1 and dur < 0.001:
        raise parse_error("Shaper duration is too small (%.6f sec)" % (dur,))
    if dur > 0.2:
        raise parse_error("Shaper duration is too large (%.6f sec)" % (dur,))
    return len(A), A, T

class CustomInputShaperParams:
    def __init__(self, axis, n, A, T):
        self.axis = axis
        self.n, self.A, self.T = n, A, T
    @classmethod
    def init_from_cfg(cls, axis, config):
        shaper_a_str = config.get('shaper_a_' + axis)
        shaper_t_str = config.get('shaper_t_' + axis)
        n, A, T = parse_custom_shaper(
                shaper_a_str, shaper_t_str, config.error)
        return CustomInputShaperParams(axis, n, A, T)
    @classmethod
    def init_from_gcmd(cls, axis, gcmd):
        shaper_a_str = gcmd.get('SHAPER_A_' + axis.upper())
        shaper_t_str = gcmd.get('SHAPER_T_' + axis.upper())
        n, A, T = parse_custom_shaper(
                shaper_a_str, shaper_t_str, gcmd.error)
        return CustomInputShaperParams(axis, n, A, T)
    def is_compatible(self, shaper_type):
        return shaper_type == CUSTOM_SHAPER
    def update_from_gcmd(self, shaper_type, gcmd):
        axis = self.axis.upper()
        shaper_a_str = gcmd.get('SHAPER_A_' + axis, None)
        shaper_t_str = gcmd.get('SHAPER_T_' + axis, None)
        if (shaper_a_str is None) != (shaper_t_str is None):
            raise gcmd.error("Both SHAPER_A_%s and SHAPER_T_%s parameters"
                             " must be provided" % (axis, axis))
        if shaper_a_str is not None:
            self.n, self.A, self.T = parse_custom_shaper(
                    shaper_a_str, shaper_t_str, gcmd.error)
    def get_shaper(self):
        return self.n, self.A, self.T
    def get_shaper_type(self):
        return CUSTOM_SHAPER
    def get_status(self):
        if not self.n:
            return collections.OrderedDict([
                ('shaper_type', CUSTOM_SHAPER + '/disabled')])
        return collections.OrderedDict([
            ('shaper_type', CUSTOM_SHAPER),
            ('shaper_a', ','.join(['%.6f' % (a,) for a in self.A])),
            ('shaper_t', ','.join(['%.6f' % (t,) for t in self.T]))])

class InputShaperParams:
    shapers = {s.name : s.init_func for s in shaper_defs.INPUT_SHAPERS}
    def __init__(self, axis, shaper_type, shaper_freq, damping_ratio):
        self.axis = axis
        self.shaper_type = shaper_type
        self.shaper_freq = shaper_freq
        self.damping_ratio = damping_ratio
    @classmethod
    def init_from_cfg(cls, axis, shaper_type, config):
        if shaper_type == CUSTOM_SHAPER:
            return CustomInputShaperParams.init_from_cfg(axis, config)
        if shaper_type not in cls.shapers:
            raise config.error(
                    'Unsupported shaper type: %s' % (shaper_type,))
        damping_ratio = config.getfloat(
                'damping_ratio_' + axis, shaper_defs.DEFAULT_DAMPING_RATIO,
                minval=0., maxval=1.)
        shaper_freq = config.getfloat('shaper_freq_' + axis, 0., minval=0.)
        return InputShaperParams(axis, shaper_type, shaper_freq, damping_ratio)
    @classmethod
    def init_from_gcmd(cls, axis, shaper_type, gcmd):
        if shaper_type == CUSTOM_SHAPER:
            return CustomInputShaperParams.init_from_gcmd(axis, gcmd)
        if shaper_type not in cls.shapers:
            raise gcmd.error(
                    'Unsupported shaper type: %s' % (shaper_type,))
        damping_ratio = gcmd.get_float('DAMPING_RATIO_' + axis.upper(),
                                       shaper_defs.DEFAULT_DAMPING_RATIO,
                                       minval=0., maxval=1.)
        shaper_freq = gcmd.get_float('SHAPER_FREQ_' + axis.upper(),
                                     0., minval=0.)
        return InputShaperParams(axis, shaper_type, shaper_freq, damping_ratio)
    def update_from_gcmd(self, shaper_type, gcmd):
        if shaper_type not in self.shapers:
            raise gcmd.error(
                    'Unsupported shaper type: %s' % (shaper_type,))
        axis = self.axis.upper()
        self.damping_ratio = gcmd.get_float('DAMPING_RATIO_' + axis,
                                            self.damping_ratio,
                                            minval=0., maxval=1.)
        self.shaper_freq = gcmd.get_float('SHAPER_FREQ_' + axis,
                                          self.shaper_freq, minval=0.)
        self.shaper_type = shaper_type
    def is_compatible(self, shaper_type):
        return shaper_type != CUSTOM_SHAPER
    def get_shaper(self):
        if not self.shaper_freq:
            A, T = shaper_defs.get_none_shaper()
        else:
            A, T = self.shapers[self.shaper_type](self.shaper_freq,
                                                  self.damping_ratio)
        return len(A), A, T
    def get_shaper_type(self):
        return self.shaper_type
    def get_status(self):
        if not self.shaper_freq:
            return collections.OrderedDict([
                ('shaper_type', self.shaper_type + '/disabled')])
        return collections.OrderedDict([
            ('shaper_type', self.shaper_type),
            ('shaper_freq', '%.3f' % (self.shaper_freq,)),
            ('damping_ratio', '%.6f' % (self.damping_ratio,))])

class AxisInputShaper:
    def __init__(self, axis, config):
        self.axis = axis
        shaper_type = config.get('shaper_type', 'mzv')
        shaper_type = config.get('shaper_type_' + axis, shaper_type).lower()
        self.params = InputShaperParams.init_from_cfg(axis, shaper_type, config)
        self.n, self.A, self.T = self.params.get_shaper()
        self.saved = None
    def get_name(self):
        return 'shaper_' + self.axis
    def get_shaper(self):
        return self.n, self.A, self.T
    def update(self, gcmd):
        shaper_type = gcmd.get('SHAPER_TYPE', None)
        if shaper_type is None:
            shaper_type = gcmd.get('SHAPER_TYPE_' + self.axis.upper(),
                                   self.params.get_shaper_type())
        shaper_type = shaper_type.lower()
        if self.params.is_compatible(shaper_type):
            self.params.update_from_gcmd(shaper_type, gcmd)
        else:
            self.params = InputShaperParams.init_from_gcmd(self.axis,
                                                           shaper_type, gcmd)
        old_n, old_A, old_T = self.n, self.A, self.T
        self.n, self.A, self.T = self.params.get_shaper()
        return (old_n, old_A, old_T) != (self.n, self.A, self.T)
    def set_shaper_kinematics(self, sk):
        ffi_main, ffi_lib = chelper.get_ffi()
        success = ffi_lib.input_shaper_set_shaper_params(
                sk, self.axis.encode(), self.n, self.A, self.T) == 0
        if not success:
            self.disable_shaping()
            ffi_lib.input_shaper_set_shaper_params(
                    sk, self.axis.encode(), self.n, self.A, self.T)
        return success
    def get_step_generation_window(self):
        ffi_main, ffi_lib = chelper.get_ffi()
        return ffi_lib.input_shaper_get_step_generation_window(self.n,
                                                               self.A, self.T)
    def disable_shaping(self):
        if self.saved is None and self.n:
            self.saved = (self.n, self.A, self.T)
        A, T = shaper_defs.get_none_shaper()
        self.n, self.A, self.T = len(A), A, T
    def enable_shaping(self):
        if self.saved is None:
            # Input shaper was not disabled
            return
        self.n, self.A, self.T = self.saved
        self.saved = None
    def report(self, gcmd):
        info = ' '.join(["%s_%s:%s" % (key, self.axis, value)
                         for (key, value) in self.params.get_status().items()])
        gcmd.respond_info(info)

class InputShaper:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.printer.register_event_handler("klippy:connect", self.connect)
        self.toolhead = None
        self.shapers = [AxisInputShaper('x', config),
                        AxisInputShaper('y', config)]
        self.stepper_kinematics = []
        self.orig_stepper_kinematics = []
        # Register gcode commands
        gcode = self.printer.lookup_object('gcode')
        gcode.register_command("SET_INPUT_SHAPER",
                               self.cmd_SET_INPUT_SHAPER,
                               desc=self.cmd_SET_INPUT_SHAPER_help)
    def get_shapers(self):
        return self.shapers
    def connect(self):
        self.toolhead = self.printer.lookup_object("toolhead")
        kin = self.toolhead.get_kinematics()
        # Lookup stepper kinematics
        ffi_main, ffi_lib = chelper.get_ffi()
        steppers = kin.get_steppers()
        for s in steppers:
            sk = ffi_main.gc(ffi_lib.input_shaper_alloc(), ffi_lib.free)
            orig_sk = s.set_stepper_kinematics(sk)
            res = ffi_lib.input_shaper_set_sk(sk, orig_sk)
            if res < 0:
                s.set_stepper_kinematics(orig_sk)
                continue
            self.stepper_kinematics.append(sk)
            self.orig_stepper_kinematics.append(orig_sk)
        # Configure initial values
        self.old_delay = 0.
        self._update_input_shaping(error=self.printer.config_error)
    def _update_input_shaping(self, error=None):
        self.toolhead.flush_step_generation()
        new_delay = max([s.get_step_generation_window() for s in self.shapers])
        self.toolhead.note_step_generation_scan_time(new_delay,
                                                     old_delay=self.old_delay)
        failed = []
        for sk in self.stepper_kinematics:
            for shaper in self.shapers:
                if shaper in failed:
                    continue
                if not shaper.set_shaper_kinematics(sk):
                    failed.append(shaper)
        if failed:
            error = error or self.printer.command_error
            raise error("Failed to configure shaper(s) %s with given parameters"
                        % (', '.join([s.get_name() for s in failed])))
    def disable_shaping(self):
        for shaper in self.shapers:
            shaper.disable_shaping()
        self._update_input_shaping()
    def enable_shaping(self):
        for shaper in self.shapers:
            shaper.enable_shaping()
        self._update_input_shaping()
    cmd_SET_INPUT_SHAPER_help = "Set cartesian parameters for input shaper"
    def cmd_SET_INPUT_SHAPER(self, gcmd):
        updated = False
        for shaper in self.shapers:
            updated |= shaper.update(gcmd)
        if updated:
            self._update_input_shaping()
        for shaper in self.shapers:
            shaper.report(gcmd)

def load_config(config):
    return InputShaper(config)
