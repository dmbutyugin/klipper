# Code for handling printer nozzle extruders
#
# Copyright (C) 2016-2022  Kevin O'Connor <kevin@koconnor.net>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import math, logging
import stepper, chelper

class ExtruderStepper:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.name = config.get_name().split()[-1]
        self.linear_advance = self.linear_velocity = self.linear_offset = 0.
        self.pressure_advance_smooth_time = 0.
        self.pressure_advance_offset_time = 0.
        self.input_shaper_step_gen_window = 0.
        config_pa = config.getfloat('pressure_advance', None, minval=0.)
        if config_pa is not None:
            self.config_la = config_pa
            self.config_lo = self.config_lv = 0.
        else:
            self.config_la = config.getfloat('linear_advance', 0., minval=0.)
            self.config_lo = config.getfloat('linear_offset', 0., minval=0.)
            self.config_lv = config.getfloat('linear_velocity', 0., minval=0.,
                                             above=(0. if self.config_lo
                                                    else -1.))
        self.config_smooth_time = config.getfloat(
                'pressure_advance_smooth_time', 0.040, above=0., maxval=.200)
        self.config_offset_time = config.getfloat(
                'pressure_advance_offset_time', 0.0, minval=-0.2, maxval=0.2)
        # Setup stepper
        self.stepper = stepper.PrinterStepper(config)
        ffi_main, ffi_lib = chelper.get_ffi()
        self.sk_extruder = ffi_main.gc(ffi_lib.extruder_stepper_alloc(),
                                       ffi_lib.free)
        self.stepper.set_stepper_kinematics(self.sk_extruder)
        self.motion_queue = None
        # Register commands
        self.printer.register_event_handler("klippy:connect",
                                            self._handle_connect)
        gcode = self.printer.lookup_object('gcode')
        if self.name == 'extruder':
            gcode.register_mux_command("SET_PRESSURE_ADVANCE", "EXTRUDER", None,
                                       self.cmd_default_SET_PRESSURE_ADVANCE,
                                       desc=self.cmd_SET_PRESSURE_ADVANCE_help)
        gcode.register_mux_command("SET_PRESSURE_ADVANCE", "EXTRUDER",
                                   self.name, self.cmd_SET_PRESSURE_ADVANCE,
                                   desc=self.cmd_SET_PRESSURE_ADVANCE_help)
        gcode.register_mux_command("SET_EXTRUDER_ROTATION_DISTANCE", "EXTRUDER",
                                   self.name, self.cmd_SET_E_ROTATION_DISTANCE,
                                   desc=self.cmd_SET_E_ROTATION_DISTANCE_help)
        gcode.register_mux_command("SYNC_EXTRUDER_MOTION", "EXTRUDER",
                                   self.name, self.cmd_SYNC_EXTRUDER_MOTION,
                                   desc=self.cmd_SYNC_EXTRUDER_MOTION_help)
        gcode.register_mux_command("SET_EXTRUDER_STEP_DISTANCE", "EXTRUDER",
                                   self.name, self.cmd_SET_E_STEP_DISTANCE,
                                   desc=self.cmd_SET_E_STEP_DISTANCE_help)
        gcode.register_mux_command("SYNC_STEPPER_TO_EXTRUDER", "STEPPER",
                                   self.name, self.cmd_SYNC_STEPPER_TO_EXTRUDER,
                                   desc=self.cmd_SYNC_STEPPER_TO_EXTRUDER_help)
    def _handle_connect(self):
        toolhead = self.printer.lookup_object('toolhead')
        toolhead.register_step_generator(self.stepper.generate_steps)
        input_shaper = self.printer.lookup_object('input_shaper', None)
        if input_shaper:
            input_shaper.add_extruder(self)
        self._set_pressure_advance(self.config_lv, self.config_lo,
                                   self.config_la, self.config_smooth_time,
                                   self.config_offset_time)
    def get_status(self, eventtime):
        return {'pressure_advance': self.linear_advance,
                'linear_velocity': self.linear_velocity,
                'linear_offset': self.linear_offset,
                'smooth_time': self.pressure_advance_smooth_time,
                'offset_time': self.pressure_advance_offset_time,
                'motion_queue': self.motion_queue}
    def find_past_position(self, print_time):
        mcu_pos = self.stepper.get_past_mcu_position(print_time)
        return self.stepper.mcu_to_commanded_position(mcu_pos)
    def sync_to_extruder(self, extruder_name):
        toolhead = self.printer.lookup_object('toolhead')
        toolhead.flush_step_generation()
        if not extruder_name:
            self.stepper.set_trapq(None)
            self.motion_queue = None
            return
        extruder = self.printer.lookup_object(extruder_name, None)
        if extruder is None or not isinstance(extruder, PrinterExtruder):
            raise self.printer.command_error("'%s' is not a valid extruder."
                                             % (extruder_name,))
        self.stepper.set_position(extruder.last_position)
        self.stepper.set_trapq(extruder.get_trapq())
        self.motion_queue = extruder_name
    def _set_pressure_advance(self, linear_velocity, linear_offset,
                              linear_advance, smooth_time, offset_time):
        old_smooth_time = self.pressure_advance_smooth_time
        if not self.linear_advance and not self.linear_offset:
            old_smooth_time = 0.
        old_offset_time = self.pressure_advance_offset_time
        new_smooth_time = smooth_time
        if not linear_advance and not linear_offset:
            new_smooth_time = 0.
        is_step_gen_window = self.input_shaper_step_gen_window
        toolhead = self.printer.lookup_object("toolhead")
        toolhead.note_step_generation_scan_time(
                new_smooth_time*.5 + abs(offset_time) + is_step_gen_window,
                old_delay=(old_smooth_time*.5 + abs(old_offset_time) +
                           is_step_gen_window))
        ffi_main, ffi_lib = chelper.get_ffi()
        espa = ffi_lib.extruder_set_pressure_advance
        espa(self.sk_extruder, linear_velocity, linear_offset, linear_advance,
             new_smooth_time, offset_time)
        self.linear_velocity = linear_velocity
        self.linear_offset = linear_offset
        self.linear_advance = linear_advance
        self.pressure_advance_smooth_time = smooth_time
        self.pressure_advance_offset_time = offset_time
    def update_input_shaping(self, axis_shaper, input_shaper_step_gen_window):
        smooth_time = self.pressure_advance_smooth_time
        if not self.linear_advance and not self.linear_offset:
            smooth_time = 0.
        old_input_shaper_step_gen_window = self.input_shaper_step_gen_window
        toolhead = self.printer.lookup_object("toolhead")
        toolhead.note_step_generation_scan_time(
                smooth_time * .5 + input_shaper_step_gen_window,
                old_delay=smooth_time * .5 + old_input_shaper_step_gen_window)
        ffi_main, ffi_lib = chelper.get_ffi()
        axis, n, A, T = axis_shaper.get_shaper()
        success = ffi_lib.extruder_set_shaper_params(
                self.sk_extruder, axis.encode(), n, A, T) == 0
        if not success:
            axis_shaper.disable_shaping()
            axis, n, A, T = axis_shaper.get_shaper()
            ffi_lib.extruder_set_shaper_params(self.sk_extruder,
                                               axis.encode(), n, A, T)
        self.input_shaper_step_gen_window = input_shaper_step_gen_window
        return success
    cmd_SET_PRESSURE_ADVANCE_help = "Set pressure advance parameters"
    def cmd_default_SET_PRESSURE_ADVANCE(self, gcmd):
        extruder = self.printer.lookup_object('toolhead').get_extruder()
        if extruder.extruder_stepper is None:
            raise gcmd.error("Active extruder does not have a stepper")
        strapq = extruder.extruder_stepper.stepper.get_trapq()
        if strapq is not extruder.get_trapq():
            raise gcmd.error("Unable to infer active extruder stepper")
        extruder.extruder_stepper.cmd_SET_PRESSURE_ADVANCE(gcmd)
    def cmd_SET_PRESSURE_ADVANCE(self, gcmd):
        smooth_time = gcmd.get_float('SMOOTH_TIME',
                                     self.pressure_advance_smooth_time,
                                     minval=0., maxval=.200)
        offset_time = gcmd.get_float('OFFSET_TIME',
                                     self.pressure_advance_offset_time,
                                     minval=-0.2, maxval=0.2)
        pressure_advance = gcmd.get_float('ADVANCE', None, minval=0.)
        if pressure_advance is not None:
            self._set_pressure_advance(
                    0., 0., pressure_advance, smooth_time, offset_time)
            msg = ("pressure_advance: %.6f\n"
                   "pressure_advance_smooth_time: %.6f\n"
                   "pressure_advance_offset_time: %.6f"
                   % (pressure_advance, smooth_time, offset_time))
        else:
            linear_advance = gcmd.get_float('LINEAR_ADVANCE',
                                            self.linear_advance, minval=0.)
            linear_velocity = gcmd.get_float('LINEAR_VELOCITY',
                                             self.linear_velocity, minval=0.)
            linear_offset = gcmd.get_float('LINEAR_OFFSET',
                                           self.linear_offset, minval=0.)
            if linear_offset > 0. and linear_velocity <= 0.:
                raise gcmd.error('LINEAR_VELOCITY must be set to a positive '
                                 'value when LINEAR_OFFSET is configured')
            self._set_pressure_advance(linear_velocity, linear_offset,
                                       linear_advance, smooth_time, offset_time)
            msg = ("linear_advance: %.6f\n"
                   "linear_velocity: %.6f\n"
                   "linear_offset: %.6f\n"
                   "pressure_advance_smooth_time: %.6f\n"
                   "pressure_advance_offset_time: %.6f"
                   % (linear_advance, linear_velocity, linear_offset,
                      smooth_time, offset_time))
        self.printer.set_rollover_info(self.name, "%s: %s" % (self.name, msg))
        gcmd.respond_info(msg, log=False)
    cmd_SET_E_ROTATION_DISTANCE_help = "Set extruder rotation distance"
    def cmd_SET_E_ROTATION_DISTANCE(self, gcmd):
        rotation_dist = gcmd.get_float('DISTANCE', None)
        if rotation_dist is not None:
            if not rotation_dist:
                raise gcmd.error("Rotation distance can not be zero")
            invert_dir, orig_invert_dir = self.stepper.get_dir_inverted()
            next_invert_dir = orig_invert_dir
            if rotation_dist < 0.:
                next_invert_dir = not orig_invert_dir
                rotation_dist = -rotation_dist
            toolhead = self.printer.lookup_object('toolhead')
            toolhead.flush_step_generation()
            self.stepper.set_rotation_distance(rotation_dist)
            self.stepper.set_dir_inverted(next_invert_dir)
        else:
            rotation_dist, spr = self.stepper.get_rotation_distance()
        invert_dir, orig_invert_dir = self.stepper.get_dir_inverted()
        if invert_dir != orig_invert_dir:
            rotation_dist = -rotation_dist
        gcmd.respond_info("Extruder '%s' rotation distance set to %0.6f"
                          % (self.name, rotation_dist))
    cmd_SYNC_EXTRUDER_MOTION_help = "Set extruder stepper motion queue"
    def cmd_SYNC_EXTRUDER_MOTION(self, gcmd):
        ename = gcmd.get('MOTION_QUEUE')
        self.sync_to_extruder(ename)
        gcmd.respond_info("Extruder '%s' now syncing with '%s'"
                          % (self.name, ename))
    cmd_SET_E_STEP_DISTANCE_help = "Set extruder step distance"
    def cmd_SET_E_STEP_DISTANCE(self, gcmd):
        step_dist = gcmd.get_float('DISTANCE', None, above=0.)
        if step_dist is not None:
            toolhead = self.printer.lookup_object('toolhead')
            toolhead.flush_step_generation()
            rd, steps_per_rotation = self.stepper.get_rotation_distance()
            self.stepper.set_rotation_distance(step_dist * steps_per_rotation)
        else:
            step_dist = self.stepper.get_step_dist()
        gcmd.respond_info("Extruder '%s' step distance set to %0.6f"
                          % (self.name, step_dist))
    cmd_SYNC_STEPPER_TO_EXTRUDER_help = "Set extruder stepper"
    def cmd_SYNC_STEPPER_TO_EXTRUDER(self, gcmd):
        ename = gcmd.get('EXTRUDER')
        self.sync_to_extruder(ename)
        gcmd.respond_info("Extruder '%s' now syncing with '%s'"
                          % (self.name, ename))

# Tracking for hotend heater, extrusion motion queue, and extruder stepper
class PrinterExtruder:
    def __init__(self, config, extruder_num):
        self.printer = config.get_printer()
        self.name = config.get_name()
        self.last_position = [0., 0., 0.]
        # Setup hotend heater
        shared_heater = config.get('shared_heater', None)
        pheaters = self.printer.load_object(config, 'heaters')
        gcode_id = 'T%d' % (extruder_num,)
        if shared_heater is None:
            self.heater = pheaters.setup_heater(config, gcode_id)
        else:
            config.deprecate('shared_heater')
            self.heater = pheaters.lookup_heater(shared_heater)
        # Setup kinematic checks
        self.nozzle_diameter = config.getfloat('nozzle_diameter', above=0.)
        filament_diameter = config.getfloat(
            'filament_diameter', minval=self.nozzle_diameter)
        self.filament_area = math.pi * (filament_diameter * .5)**2
        def_max_cross_section = 4. * self.nozzle_diameter**2
        def_max_extrude_ratio = def_max_cross_section / self.filament_area
        max_cross_section = config.getfloat(
            'max_extrude_cross_section', def_max_cross_section, above=0.)
        self.max_extrude_ratio = max_cross_section / self.filament_area
        logging.info("Extruder max_extrude_ratio=%.6f", self.max_extrude_ratio)
        toolhead = self.printer.lookup_object('toolhead')
        max_velocity, max_accel = toolhead.get_max_velocity()
        self.max_e_velocity = config.getfloat(
            'max_extrude_only_velocity', max_velocity * def_max_extrude_ratio
            , above=0.)
        self.max_e_accel = config.getfloat(
            'max_extrude_only_accel', max_accel * def_max_extrude_ratio
            , above=0.)
        self.max_e_dist = config.getfloat(
            'max_extrude_only_distance', 50., minval=0.)
        self.instant_corner_v = config.getfloat(
            'instantaneous_corner_velocity', 1., minval=0.)
        # Setup extruder trapq (trapezoidal motion queue)
        ffi_main, ffi_lib = chelper.get_ffi()
        self.trapq = ffi_main.gc(ffi_lib.trapq_alloc(), ffi_lib.trapq_free)
        self.trapq_append = ffi_lib.trapq_append
        self.trapq_finalize_moves = ffi_lib.trapq_finalize_moves
        # Setup extruder stepper
        self.extruder_stepper = None
        if (config.get('step_pin', None) is not None
            or config.get('dir_pin', None) is not None
            or config.get('rotation_distance', None) is not None):
            self.extruder_stepper = ExtruderStepper(config)
            self.extruder_stepper.stepper.set_trapq(self.trapq)
        # Register commands
        gcode = self.printer.lookup_object('gcode')
        if self.name == 'extruder':
            toolhead.set_extruder(self, 0.)
            gcode.register_command("M104", self.cmd_M104)
            gcode.register_command("M109", self.cmd_M109)
        gcode.register_mux_command("ACTIVATE_EXTRUDER", "EXTRUDER",
                                   self.name, self.cmd_ACTIVATE_EXTRUDER,
                                   desc=self.cmd_ACTIVATE_EXTRUDER_help)
    def update_move_time(self, flush_time):
        self.trapq_finalize_moves(self.trapq, flush_time)
    def get_status(self, eventtime):
        sts = self.heater.get_status(eventtime)
        sts['can_extrude'] = self.heater.can_extrude
        if self.extruder_stepper is not None:
            sts.update(self.extruder_stepper.get_status(eventtime))
        return sts
    def get_name(self):
        return self.name
    def get_heater(self):
        return self.heater
    def get_trapq(self):
        return self.trapq
    def stats(self, eventtime):
        return self.heater.stats(eventtime)
    def check_move(self, move):
        axis_r = move.axes_r[3]
        if not self.heater.can_extrude:
            raise self.printer.command_error(
                "Extrude below minimum temp\n"
                "See the 'min_extrude_temp' config option for details")
        if (not move.axes_d[0] and not move.axes_d[1]) or axis_r < 0.:
            # Extrude only move (or retraction move) - limit accel and velocity
            if abs(move.axes_d[3]) > self.max_e_dist:
                raise self.printer.command_error(
                    "Extrude only move too long (%.3fmm vs %.3fmm)\n"
                    "See the 'max_extrude_only_distance' config"
                    " option for details" % (move.axes_d[3], self.max_e_dist))
            inv_extrude_r = 1. / abs(axis_r)
            move.limit_speed(self.max_e_velocity * inv_extrude_r,
                             self.max_e_accel * inv_extrude_r)
        elif axis_r > self.max_extrude_ratio:
            if move.axes_d[3] <= self.nozzle_diameter * self.max_extrude_ratio:
                # Permit extrusion if amount extruded is tiny
                return
            area = axis_r * self.filament_area
            logging.debug("Overextrude: %s vs %s (area=%.3f dist=%.3f)",
                          axis_r, self.max_extrude_ratio, area, move.move_d)
            raise self.printer.command_error(
                "Move exceeds maximum extrusion (%.3fmm^2 vs %.3fmm^2)\n"
                "See the 'max_extrude_cross_section' config option for details"
                % (area, self.max_extrude_ratio * self.filament_area))
    def calc_junction(self, prev_move, move):
        diff_r = move.axes_r[3] - prev_move.axes_r[3]
        if diff_r:
            return (self.instant_corner_v / abs(diff_r))**2
        return move.max_cruise_v2
    def move(self, print_time, move):
        axis_r = move.axes_r[3]
        accel = move.accel * abs(axis_r)
        start_v = move.start_v * abs(axis_r)
        cruise_v = move.cruise_v * abs(axis_r)
        extr_pos = self.last_position
        if move.is_kinematic_move:
            # Regular kinematic move with extrusion
            extr_r = [math.copysign(r * r, axis_r) for r in move.axes_r[:3]]
        else:
            # Extrude-only move, do not apply pressure advance
            extr_r = [0., 0., axis_r]
        self.trapq_append(self.trapq, print_time,
                          move.accel_t, move.cruise_t, move.decel_t,
                          extr_pos[0], extr_pos[1], extr_pos[2],
                          extr_r[0], extr_r[1], extr_r[2],
                          start_v, cruise_v, accel)
        for i in range(3):
            self.last_position[i] += move.axes_d[3] * abs(extr_r[i])
    def find_past_position(self, print_time):
        if self.extruder_stepper is None:
            return 0.
        return self.extruder_stepper.find_past_position(print_time)
    def cmd_M104(self, gcmd, wait=False):
        # Set Extruder Temperature
        temp = gcmd.get_float('S', 0.)
        index = gcmd.get_int('T', None, minval=0)
        if index is not None:
            section = 'extruder'
            if index:
                section = 'extruder%d' % (index,)
            extruder = self.printer.lookup_object(section, None)
            if extruder is None:
                if temp <= 0.:
                    return
                raise gcmd.error("Extruder not configured")
        else:
            extruder = self.printer.lookup_object('toolhead').get_extruder()
        pheaters = self.printer.lookup_object('heaters')
        pheaters.set_temperature(extruder.get_heater(), temp, wait)
    def cmd_M109(self, gcmd):
        # Set Extruder Temperature and Wait
        self.cmd_M104(gcmd, wait=True)
    cmd_ACTIVATE_EXTRUDER_help = "Change the active extruder"
    def cmd_ACTIVATE_EXTRUDER(self, gcmd):
        toolhead = self.printer.lookup_object('toolhead')
        if toolhead.get_extruder() is self:
            gcmd.respond_info("Extruder %s already active" % (self.name,))
            return
        gcmd.respond_info("Activating extruder %s" % (self.name,))
        toolhead.flush_step_generation()
        toolhead.set_extruder(self, sum(self.last_position))
        self.printer.send_event("extruder:activate_extruder")

# Dummy extruder class used when a printer has no extruder at all
class DummyExtruder:
    def __init__(self, printer):
        self.printer = printer
    def update_move_time(self, flush_time):
        pass
    def check_move(self, move):
        raise move.move_error("Extrude when no extruder present")
    def find_past_position(self, print_time):
        return 0.
    def calc_junction(self, prev_move, move):
        return move.max_cruise_v2
    def get_name(self):
        return ""
    def get_heater(self):
        raise self.printer.command_error("Extruder not configured")
    def get_trapq(self):
        raise self.printer.command_error("Extruder not configured")

def add_printer_objects(config):
    printer = config.get_printer()
    for i in range(99):
        section = 'extruder'
        if i:
            section = 'extruder%d' % (i,)
        if not config.has_section(section):
            break
        pe = PrinterExtruder(config.getsection(section), i)
        printer.add_object(section, pe)
