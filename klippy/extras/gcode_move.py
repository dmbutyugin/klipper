# G-Code G1 movement commands (and associated coordinate manipulation)
#
# Copyright (C) 2016-2025  Kevin O'Connor <kevin@koconnor.net>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import logging

class GCodeMove:
    def __init__(self, config):
        self.printer = printer = config.get_printer()
        self.is_printer_ready = False
        # Register g-code commands
        gcode = printer.lookup_object('gcode')
        handlers = [
            'G1', 'G20', 'G21',
            'M82', 'M83', 'G90', 'G91', 'G92', 'M220', 'M221',
            'SET_GCODE_OFFSET', 'SET_GCODE_POSITION',
            'SAVE_GCODE_STATE', 'RESTORE_GCODE_STATE',
        ]
        for cmd in handlers:
            func = getattr(self, 'cmd_' + cmd)
            desc = getattr(self, 'cmd_' + cmd + '_help', None)
            gcode.register_command(cmd, func, False, desc)
        gcode.register_command('G0', self.cmd_G1)
        gcode.register_command('M114', self.cmd_M114, True)
        gcode.register_command('GET_POSITION', self.cmd_GET_POSITION, True,
                               desc=self.cmd_GET_POSITION_help)
        self.Coord = gcode.Coord
        # G-Code coordinate manipulation
        self.absolute_coord = self.absolute_extrude = True
        self.base_position = [0.0, 0.0, 0.0, 0.0]
        self.last_position = [0.0, 0.0, 0.0, 0.0]
        self.homing_position = [0.0, 0.0, 0.0, 0.0]
        self.axis_map = {'X':0, 'Y': 1, 'Z': 2, 'E': 3}
        self.speed = 25.
        self.speed_factor = 1. / 60.
        self.extrude_factors = {'E': 1.}
        # G-Code state
        self.saved_states = {}
        self.move_transform = self.move_with_transform = None
        self.position_with_transform = (lambda: [0., 0., 0., 0.])
        # Register callbacks
        printer.register_event_handler("klippy:ready", self._handle_ready)
        printer.register_event_handler("klippy:shutdown", self._handle_shutdown)
        printer.register_event_handler("klippy:analyze_shutdown",
                                       self._handle_analyze_shutdown)
        printer.register_event_handler("toolhead:set_position",
                                       self.reset_last_position)
        printer.register_event_handler("toolhead:manual_move",
                                       self.reset_last_position)
        printer.register_event_handler("toolhead:update_extra_axes",
                                       self._update_extra_axes)
        printer.register_event_handler("gcode:command_error",
                                       self.reset_last_position)
        printer.register_event_handler("extruder:activate_extruder",
                                       self._handle_activate_extruder)
        printer.register_event_handler("homing:home_rails_end",
                                       self._handle_home_rails_end)
    def _handle_ready(self):
        self.is_printer_ready = True
        if self.move_transform is None:
            toolhead = self.printer.lookup_object('toolhead')
            self.move_with_transform = toolhead.move
            self.position_with_transform = toolhead.get_position
        self.reset_last_position()
    def _handle_shutdown(self):
        self.is_printer_ready = False
    def _handle_analyze_shutdown(self, msg, details):
        logging.info("gcode state: absolute_coord=%s absolute_extrude=%s"
                     " base_position=%s last_position=%s homing_position=%s"
                     " speed_factor=%s extrude_factors={%s} speed=%s",
                     self.absolute_coord, self.absolute_extrude,
                     self.base_position, self.last_position,
                     self.homing_position, self.speed_factor,
                     ','.join(["'%s':%f" % (a, f)
                               for a, f in self.extrude_factors.items()]),
                     self.speed)
    def _handle_activate_extruder(self, extruder):
        self.reset_last_position()
        gcode_id = extruder.get_axis_gcode_id()
        self.extrude_factors[gcode_id] = 1.
        pos = self.axis_map[gcode_id]
        self.base_position[pos] = self.last_position[pos]
    def _handle_home_rails_end(self, homing_state, rails):
        self.reset_last_position()
        for axis in homing_state.get_axes():
            self.base_position[axis] = self.homing_position[axis]
    def set_move_transform(self, transform, force=False):
        if self.move_transform is not None and not force:
            raise self.printer.config_error(
                "G-Code move transform already specified")
        old_transform = self.move_transform
        if old_transform is None:
            old_transform = self.printer.lookup_object('toolhead', None)
        self.move_transform = transform
        self.move_with_transform = transform.move
        self.position_with_transform = transform.get_position
        return old_transform
    def _get_gcode_position(self):
        p = [lp - bp for lp, bp in zip(self.last_position, self.base_position)]
        for axis, factor in self.extrude_factors.items():
            p[self.axis_map[axis]] /= factor
        return p
    def _get_gcode_speed(self):
        return self.speed / self.speed_factor
    def _get_gcode_speed_override(self):
        return self.speed_factor * 60.
    def get_status(self, eventtime=None):
        move_position = self._get_gcode_position()
        return {
            'speed_factor': self._get_gcode_speed_override(),
            'speed': self._get_gcode_speed(),
            'extrude_factor': self.extrude_factors['E'],
            'extrude_factors': { a.lower() : f
                                for a, f in self.extrude_factors.items()},
            'absolute_coordinates': self.absolute_coord,
            'absolute_extrude': self.absolute_extrude,
            'homing_origin': self.Coord(self.homing_position),
            'position': self.Coord(self.last_position),
            'gcode_position': self.Coord(move_position),
            'axis_map': self.axis_map,
        }
    def reset_last_position(self):
        if self.is_printer_ready:
            self.last_position = self.position_with_transform()
    def _update_extra_axes(self):
        toolhead = self.printer.lookup_object('toolhead')
        axis_map = {'X':0, 'Y': 1, 'Z': 2, 'E': 3}
        extra_axes = toolhead.get_extra_axes()
        for index, ea in enumerate(extra_axes):
            if ea is None:
                continue
            gcode_id = ea.get_axis_gcode_id()
            if (gcode_id is None or len(gcode_id) != 1 or not gcode_id.isupper()
                or gcode_id in axis_map or gcode_id in "FN"):
                continue
            axis_map[gcode_id] = index
        base_position = [0.] * len(extra_axes)
        for index, ea in enumerate(extra_axes):
            if ea is None:
                continue
            gcode_id = ea.get_axis_gcode_id()
            if gcode_id in self.axis_map:
                base_position[index] = self.base_position[
                        self.axis_map[gcode_id]]
        self.axis_map = axis_map
        self.base_position[3:] = base_position[3:]
        for axis in tuple(self.extrude_factors.keys()):
            if axis not in axis_map:
                del self.extrude_factors[axis]
        self.reset_last_position()
    # G-Code movement commands
    def cmd_G1(self, gcmd):
        # Move
        params = gcmd.get_command_parameters()
        try:
            for axis, pos in self.axis_map.items():
                if axis in params:
                    v = float(params[axis])
                    absolute_coord = self.absolute_coord
                    if axis in self.extrude_factors:
                        v *= self.extrude_factors[axis]
                        if not self.absolute_extrude:
                            absolute_coord = False
                    if not absolute_coord:
                        # value relative to position of last move
                        self.last_position[pos] += v
                    else:
                        # value relative to base coordinate position
                        self.last_position[pos] = v + self.base_position[pos]
            if 'F' in params:
                gcode_speed = float(params['F'])
                if gcode_speed <= 0.:
                    raise gcmd.error("Invalid speed in '%s'"
                                     % (gcmd.get_commandline(),))
                self.speed = gcode_speed * self.speed_factor
        except ValueError as e:
            raise gcmd.error("Unable to parse move '%s'"
                             % (gcmd.get_commandline(),))
        self.move_with_transform(self.last_position, self.speed)
    # G-Code coordinate manipulation
    def cmd_G20(self, gcmd):
        # Set units to inches
        raise gcmd.error('Machine does not support G20 (inches) command')
    def cmd_G21(self, gcmd):
        # Set units to millimeters
        pass
    def cmd_M82(self, gcmd):
        # Use absolute distances for extrusion
        self.absolute_extrude = True
    def cmd_M83(self, gcmd):
        # Use relative distances for extrusion
        self.absolute_extrude = False
    def cmd_G90(self, gcmd):
        # Use absolute coordinates
        self.absolute_coord = True
    def cmd_G91(self, gcmd):
        # Use relative coordinates
        self.absolute_coord = False
    def cmd_G92(self, gcmd):
        # Set position
        offsets = [ gcmd.get_float(a, None) for a in 'XYZE' ]
        for i, offset in enumerate(offsets):
            if offset is not None:
                if i == 3:
                    offset *= self.extrude_factors['E']
                self.base_position[i] = self.last_position[i] - offset
        if offsets == [None, None, None, None]:
            self.base_position[:4] = self.last_position[:4]
    def cmd_M114(self, gcmd):
        # Get Current Position
        p = self._get_gcode_position()
        gcmd.respond_raw("X:%.3f Y:%.3f Z:%.3f E:%.3f" % tuple(p[:4]))
    def cmd_M220(self, gcmd):
        # Set speed factor override percentage
        value = gcmd.get_float('S', 100., above=0.) / (60. * 100.)
        self.speed = self._get_gcode_speed() * value
        self.speed_factor = value
    def cmd_M221(self, gcmd):
        # Set extrude factor override percentage
        new_extrude_factor = gcmd.get_float('S', 100., above=0.) / 100.
        toolhead = self.printer.lookup_object('toolhead')
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
            extruder = toolhead.get_extruder()
        extra_axes = toolhead.get_extra_axes()
        if extruder not in extra_axes:
            return
        ea_index = extra_axes.index(extruder)
        gcode_id = extruder.get_axis_gcode_id()
        last_e_pos = self.last_position[ea_index]
        extrude_factor = self.extrude_factors[gcode_id]
        e_value = (last_e_pos - self.base_position[ea_index]) / extrude_factor
        self.base_position[ea_index] = last_e_pos - e_value * new_extrude_factor
        self.extrude_factors[gcode_id] = new_extrude_factor
    cmd_SET_GCODE_OFFSET_help = "Set a virtual offset to g-code positions"
    def cmd_SET_GCODE_OFFSET(self, gcmd):
        move_delta = [0.] * len(self.axis_map)
        for axis, pos in self.axis_map.items():
            homing_pos = self.homing_position[pos] if pos < 3 else 0.
            offset = gcmd.get_float(axis, None)
            if offset is None:
                offset = gcmd.get_float(axis + '_ADJUST', None)
                if offset is None:
                    continue
                offset += homing_pos
            delta = offset - homing_pos
            move_delta[pos] = delta
            self.base_position[pos] += delta
            if pos < 3:
                self.homing_position[pos] = offset
        # Move the toolhead the given offset if requested
        if gcmd.get_int('MOVE', 0):
            speed = gcmd.get_float('MOVE_SPEED', self.speed, above=0.)
            for axis, pos in self.axis_map.items():
                self.last_position[pos] += move_delta[pos]
            self.move_with_transform(self.last_position, speed)
    cmd_SET_GCODE_POSITION_help = "Set positions for g-code axes"
    def cmd_SET_GCODE_POSITION(self, gcmd):
        for axis, pos in self.axis_map.items():
            axis_position = gcmd.get_float(axis, None)
            if axis_position is None:
                continue
            if axis in self.extrude_factors:
                axis_position *= self.extrude_factors[axis]
            self.base_position[pos] = self.last_position[pos] - axis_position
    cmd_SAVE_GCODE_STATE_help = "Save G-Code coordinate state"
    def cmd_SAVE_GCODE_STATE(self, gcmd):
        state_name = gcmd.get('NAME', 'default')
        self.saved_states[state_name] = {
            'absolute_coord': self.absolute_coord,
            'absolute_extrude': self.absolute_extrude,
            'axis_map': dict(self.axis_map),
            'base_position': list(self.base_position),
            'last_position': list(self.last_position),
            'homing_position': list(self.homing_position),
            'speed': self.speed, 'speed_factor': self.speed_factor,
            'extrude_factors': dict(self.extrude_factors),
        }
    cmd_RESTORE_GCODE_STATE_help = "Restore a previously saved G-Code state"
    def cmd_RESTORE_GCODE_STATE(self, gcmd):
        state_name = gcmd.get('NAME', 'default')
        state = self.saved_states.get(state_name)
        if state is None:
            raise gcmd.error("Unknown g-code state: %s" % (state_name,))
        move_requested = gcmd.get_int('MOVE', 0)
        # Restore state
        self.absolute_coord = state['absolute_coord']
        self.absolute_extrude = state['absolute_extrude']
        self.base_position[:4] = state['base_position'][:4]
        self.homing_position = list(state['homing_position'])
        self.speed = state['speed']
        self.speed_factor = state['speed_factor']
        saved_axis_map = state['axis_map']
        extra_axes = self.printer.lookup_object('toolhead').get_extra_axes()
        for index, ea in enumerate(extra_axes):
            if ea is None:
                continue
            gcode_id = ea.get_axis_gcode_id()
            if gcode_id is None or gcode_id not in saved_axis_map:
                continue
            saved_index = saved_axis_map[gcode_id]
            if ea.is_extruder_axis():
                # Restore the relative E position
                e_diff = self.last_position[index] \
                        - state['last_position'][saved_index]
                self.base_position[index] += e_diff
                if gcode_id in state['extrude_factors']:
                    self.extrude_factors[gcode_id] = \
                            state['extrude_factors'][gcode_id]
                continue
            self.base_position[index] = state['base_position'][saved_index]
            if move_requested:
                self.last_position[index] = state['last_position'][saved_index]
        # Move the toolhead back if requested
        if move_requested:
            speed = gcmd.get_float('MOVE_SPEED', self.speed, above=0.)
            self.last_position[:3] = state['last_position'][:3]
            self.move_with_transform(self.last_position, speed)
    cmd_GET_POSITION_help = (
        "Return information on the current location of the toolhead")
    def cmd_GET_POSITION(self, gcmd):
        toolhead = self.printer.lookup_object('toolhead', None)
        if toolhead is None:
            raise gcmd.error("Printer not ready")
        kin = toolhead.get_kinematics()
        extra_axes = toolhead.get_extra_axes()
        steppers = kin.get_steppers()
        mcu_pos = " ".join(["%s:%d" % (s.get_name(), s.get_mcu_position())
                            for s in steppers])
        cinfo = [(s.get_name(), s.get_commanded_position()) for s in steppers]
        stepper_pos = " ".join(["%s:%.6f" % (a, v) for a, v in cinfo])
        kin_axes = ["X", "Y", "Z"] + [ea.get_axis_gcode_id()
                                      for ea in extra_axes
                                      if ea and ea.is_kinematic_axis()]
        kinfo = zip(kin_axes, kin.calc_position(dict(cinfo)))
        kin_pos = " ".join(["%s:%.6f" % (a, v) for a, v in kinfo])
        toolhead_axes = [("XYZ"[i] if i < 3 else ea.get_axis_gcode_id(), i)
                         for i, ea in enumerate(extra_axes)
                         if i < 3 or (ea and (ea.is_kinematic_axis()
                                              or ea.is_extruder_axis()))]
        tpos = toolhead.get_position()
        toolhead_pos = " ".join(["%s:%.6f" % (a, tpos[i])
                                 for a, i in toolhead_axes])
        gcode_axes = [" "] * len(self.axis_map)
        for a, i in self.axis_map.items():
            gcode_axes[i] = a
        gcode_pos = " ".join(["%s:%.6f"  % (a, v)
                              for a, v in zip(gcode_axes, self.last_position)])
        base_pos = " ".join(["%s:%.6f"  % (a, v)
                             for a, v in zip(gcode_axes, self.base_position)])
        homing_pos = " ".join(["%s:%.6f"  % (a, v)
                               for a, v in zip("XYZ", self.homing_position)])
        gcmd.respond_info("mcu: %s\n"
                          "stepper: %s\n"
                          "kinematic: %s\n"
                          "toolhead: %s\n"
                          "gcode: %s\n"
                          "gcode base: %s\n"
                          "gcode homing: %s"
                          % (mcu_pos, stepper_pos, kin_pos, toolhead_pos,
                             gcode_pos, base_pos, homing_pos))

def load_config(config):
    return GCodeMove(config)
