# Code for generic handling the kinematics of cartesian-like printers
#
# Copyright (C) 2024-2025  Dmitry Butyugin <dmbutyugin@google.com>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

import collections, copy, itertools, logging, math
import gcode, mathutil, stepper
from . import idex_modes
from . import kinematic_stepper as ks

VALID_AXES = ['x', 'y', 'z']

def mat_mul(a, b):
    if len(a[0]) != len(b):
        return None
    res = []
    for i in range(len(a)):
        res.append([])
        for j in range(len(b[0])):
            res[i].append(sum(a[i][k] * b[k][j] for k in range(len(b))))
    return res

def mat_transp(a):
    res = []
    for i in range(len(a[0])):
        res.append([a[j][i] for j in range(len(a))])
    return res

def mat_pseudo_inverse(m):
    mt = mat_transp(m)
    mtm = mat_mul(mt, m)
    pinv = mat_mul(mathutil.matrix_inv(mtm), mt)
    return pinv

class MainCarriage:
    def __init__(self, config):
        self.rail = stepper.GenericPrinterRail(config)
        self.homed = False
        carriage_name = self.rail.get_name(short=True)
        if carriage_name in VALID_AXES:
            self.axis_name = config.getchoice('axis', VALID_AXES, carriage_name)
        else:
            self.axis_name = config.getchoice('axis', VALID_AXES)
        self.axis = ord(self.axis_name) - ord('x')
        self.dual_carriage = None
    def get_name(self):
        return self.rail.get_name(short=True)
    def get_axis(self):
        return self.axis
    def get_rail(self):
        return self.rail
    def add_stepper(self, kin_stepper):
        self.rail.add_stepper(kin_stepper.get_stepper())
    def is_active(self):
        if self.dual_carriage is None:
            return True
        return self.dual_carriage.get_dc_module().is_active(self.rail)
    def is_homed(self):
        return self.homed
    def set_homed(self, homed=True):
        self.homed = homed
    def set_dual_carriage(self, carriage):
        self.dual_carriage = carriage
    def get_dual_carriage(self):
        return self.dual_carriage

class ExtraCarriage:
    def __init__(self, config, carriages):
        self.name = config.get_name().split()[-1]
        self.primary_carriage = config.getchoice('primary_carriage', carriages)
        self.endstop_pin = config.get('endstop_pin')
    def get_name(self):
        return self.name
    def get_axis(self):
        return self.primary_carriage.get_axis()
    def get_rail(self):
        return self.primary_carriage.get_rail()
    def add_stepper(self, kin_stepper):
        self.get_rail().add_stepper(kin_stepper.get_stepper(),
                                    self.endstop_pin, self.name)

class DualCarriage:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.rail = stepper.GenericPrinterRail(config)
        self.homed = False
        self.primary_carriage_name = config.get('primary_carriage', None)
        if self.primary_carriage_name is None:
            self.axis_name = config.getchoice('axis', VALID_AXES)
            self.axis = ord(self.axis_name) - ord('x')
            self.safe_dist = None
        else:
            self.axis_name = config.getchoice('axis', VALID_AXES + [None], None)
            self.safe_dist = config.getfloat('safe_distance', None, minval=0.)
        self.primary_carriage = self.dual_carriage = None
        self.config_error = config.error
    def resolve_primary_carriage(self, carriages):
        if self.primary_carriage_name is None:
            return
        if self.primary_carriage_name not in carriages:
            raise self.config_error(
                    "primary_carriage = '%s' for '%s' is not a valid choice"
                    % (self.primary_carriage_name, self.get_name()))
        self.primary_carriage = carriages[self.primary_carriage_name]
        axis_name = self.axis_name or self.primary_carriage.axis_name
        if axis_name != self.primary_carriage.axis_name:
            raise self.config_error("Mismatching axes between carriage '%s' "
                                    "(axis=%s) and dual_carriage '%s' (axis=%s)"
                                    % (self.primary_carriage.get_name(),
                                       self.primary_carriage.axis_name,
                                       self.get_name(), axis_name))
        self.axis = ord(axis_name) - ord('x')
        if self.primary_carriage.get_dual_carriage():
            raise self.config_error(
                    "Multiple dual carriages ('%s', '%s') for carriage '%s'" %
                    (self.primary_carriage.get_dual_carriage().get_name(),
                     self.get_name(), self.primary_carriage.get_name()))
        self.primary_carriage.set_dual_carriage(self)
        self.axis = self.primary_carriage.get_axis()
        if self.axis > 1:
            raise self.config_error("Invalid axis '%s' for dual_carriage '%s'" %
                                    ("xyz"[self.axis], self.get_name()))
        if self.safe_dist is None:
            dc_range = self.rail.get_range()
            pc_range = self.primary_carriage.get_rail().get_range()
            self.safe_dist = min([abs(l_pc - l_dc)
                                  for l_pc, l_dc in zip(pc_range, dc_range)])
    def get_name(self):
        return self.rail.get_name(short=True)
    def get_axis(self):
        return self.axis
    def get_rail(self):
        return self.rail
    def get_safe_dist(self):
        return self.safe_dist
    def get_dc_module(self):
        return self.printer.lookup_object('dual_carriage')
    def is_active(self):
        return self.get_dc_module().is_active(self.rail)
    def is_homed(self):
        return self.homed
    def set_homed(self, homed=True):
        self.homed = homed
    def set_dual_carriage(self, carriage):
        self.dual_carriage = carriage
    def get_dual_carriage(self):
        if self.dual_carriage is not None:
            return self.dual_carriage
        return self.primary_carriage
    def get_primary_carriage(self):
        return self.primary_carriage
    def add_stepper(self, kin_stepper):
        self.rail.add_stepper(kin_stepper.get_stepper())
    def get_axis_gcode_id(self):
        return self.get_name().upper()
    def is_kinematic_axis(self):
        return True
    def is_extruder_axis(self):
        return False
    def check_move(self, move, ea_index):
        axis = self.axis
        if not move.axes_d[ea_index] and not move.axes_d[axis]:
            return
        if self.safe_dist:
            start_pos = move.start_pos
            end_pos = move.end_pos
            ds = start_pos[axis] - start_pos[ea_index]
            de = end_pos[axis] - end_pos[ea_index]
            if min(abs(ds), abs(de)) < self.safe_dist:
                # Carriages get too close
                raise move.move_error()
            if (ds < 0) != (de < 0):
                # Carriages trajectories intersect
                raise move.move_error()
        if not move.axes_d[ea_index]:
            return
        if not self.is_homed():
            raise move.move_error(
                    "Must home carriage %s first" % self.get_axis_gcode_id())
        limits = self.get_rail().get_range()
        end_pos = move.end_pos[ea_index]
        if end_pos < limits[0] or end_pos > limits[1]:
            raise move.move_error()

class DCVirtualToolhead:
    def __init__(self, config, toolhead):
        self.printer = config.get_printer()
        self.toolhead = toolhead
        self.direct_dc_carriages = []
        self.motion_queuing = self.printer.load_object(config, 'motion_queuing')
        self.trapq = self.motion_queuing.allocate_trapq()
        self.trapq_append = self.motion_queuing.lookup_trapq_append()
    def get_active_carriages(self):
        return self.direct_dc_carriages
    def activate_direct_mode(self, dual_carriage, carriage_pos):
        self.toolhead.add_extra_axis(dual_carriage, carriage_pos)
        self.direct_dc_carriages.append(dual_carriage)
        dual_carriage.get_rail().set_trapq(self.trapq)
    def deactivate_direct_mode(self, dual_carriage):
        if dual_carriage not in self.direct_dc_carriages:
            return
        self.toolhead.remove_extra_axis(dual_carriage)
        self.direct_dc_carriages.remove(dual_carriage)
        # Restore the trapq of the main toolhead for the steppers no longer
        # associated with the carriages under direct control
        dual_carriage.get_rail().set_trapq(self.toolhead.get_trapq())
        for active_dc in self.direct_dc_carriages:
            active_dc.get_rail().set_trapq(self.trapq)
    def process_move(self, print_time, move):
        extra_axes = self.toolhead.get_extra_axes()
        axes_r = list(move.axes_r)
        start_pos = list(move.start_pos)
        for ea_index, ea in enumerate(extra_axes):
            if ea in self.direct_dc_carriages:
                dc_axis = ea.get_axis()
                start_pos[dc_axis] = start_pos[ea_index]
                axes_r[dc_axis] = axes_r[ea_index]
        if not any(axes_r):
            return
        self.trapq_append(
            self.trapq, print_time,
            move.accel_t, move.cruise_t, move.decel_t,
            start_pos[0], start_pos[1], start_pos[2],
            axes_r[0], axes_r[1], axes_r[2],
            move.start_v, move.cruise_v, move.accel)

class GenericCartesianKinematics:
    def __init__(self, toolhead, config):
        self.printer = config.get_printer()
        self.toolhead = toolhead
        self._load_kinematics(config)
        for s in self.get_steppers():
            s.set_trapq(toolhead.get_trapq())
        self.dc_module = self.dc_toolhead = None
        if self.dc_carriages:
            dc_axes = set(dc.get_axis() for dc in self.dc_carriages)
            pcs = ([pc for pc in self.primary_carriages
                    if pc.get_axis() in dc_axes] +
                   [dc for dc in self.dc_carriages
                    if dc.get_primary_carriage() is None])
            dcs = [pc.get_dual_carriage() for pc in pcs]
            primary_rails = [pc.get_rail() for pc in pcs]
            dual_rails = [dc.get_rail() if dc else None for dc in dcs]
            axes = [pc.get_axis() for pc in pcs]
            safe_dist = [dc.get_safe_dist() if dc else None for dc in dcs]
            self.dc_module = dc_module = idex_modes.DualCarriages(
                    self.printer, primary_rails, dual_rails, axes, safe_dist)
            self.dc_toolhead = DCVirtualToolhead(config, toolhead)
            zero_pos = (0., 0.)
            for acs in itertools.product(*zip(pcs, dcs)):
                for c in acs:
                    if c is None:
                        continue
                    dc_module.get_dc_rail_wrapper(c.get_rail()).activate(
                            idex_modes.PRIMARY, zero_pos)
                    dc = c.get_dual_carriage()
                    if dc is not None:
                        dc_module.get_dc_rail_wrapper(dc.get_rail()).inactivate(
                                zero_pos)
                self._check_kinematics(config.error)
            for dc in self.dc_carriages:
                dc_module.get_dc_rail_wrapper(dc.get_rail()).inactivate(
                        zero_pos)
            for pc in self.primary_carriages:
                if pc.get_axis() not in dc_axes:
                    continue
                dc_module.get_dc_rail_wrapper(pc.get_rail()).activate(
                        idex_modes.PRIMARY, zero_pos)
        else:
            self._check_kinematics(config.error)
        self.track_junction_carriages = []
        self._default_junction_carriages = [
                [c for c in self.carriages.values() if c.is_active()],
                # Dual carriages + remaining complements
                self.dc_carriages + [
                    c for c in self.carriages.values()
                    if c.is_active() and c.get_dual_carriage() is None]]
        self.test_junction_axes_set = [(0, 1, 2)]
        # Setup boundary checks
        max_velocity, max_accel = toolhead.get_max_velocity()
        self.max_z_velocity = config.getfloat('max_z_velocity', max_velocity,
                                              above=0., maxval=max_velocity)
        self.max_z_accel = config.getfloat('max_z_accel', max_accel,
                                           above=0., maxval=max_accel)
        self.limits = [None] * 3
        for carriage in self.carriages.values():
            if not carriage.is_active():
                continue
            self.limits[carriage.get_axis()] = carriage.get_rail().get_range()
        # Register gcode commands
        gcode = self.printer.lookup_object('gcode')
        gcode.register_command("SET_STEPPER_CARRIAGES",
                               self.cmd_SET_STEPPER_CARRIAGES,
                               desc=self.cmd_SET_STEPPER_CARRIAGES_help)
        gcode.register_command("SET_TRACK_CARRIAGES_JUNCTION",
                               self.cmd_SET_TRACK_CARRIAGES_JUNCTION,
                               desc=self.cmd_SET_TRACK_CARRIAGES_JUNCTION_help)
        self.printer.register_event_handler("toolhead:update_extra_axes",
                                            self._update_extra_axes)
    def _load_kinematics(self, config):
        primary_carriages = []
        for mcconfig in config.get_prefix_sections('carriage '):
            primary_carriages.append(MainCarriage(mcconfig))
        for axis, axis_name in enumerate(VALID_AXES):
            dups = [pc.get_name() for pc in primary_carriages
                    if pc.get_axis() == axis]
            if len(dups) > 1:
                raise config.error(
                        "Axis '%s' is set for multiple primary carriages (%s)"
                        % (axis_name, ', '.join(dups)))
            elif not dups:
                raise config.error(
                        "No carriage defined for axis '%s'" % axis_name)
        dc_carriages = []
        for dcconfig in config.get_prefix_sections('dual_carriage '):
            dc_carriages.append(DualCarriage(dcconfig))
        carriages = {}
        for carriage in primary_carriages + dc_carriages:
            name = carriage.get_name()
            if name in carriages:
                raise config.error("Redefinition of carriage %s" % name)
            carriages[name] = carriage
        for dc in dc_carriages:
            dc.resolve_primary_carriage(carriages)
        self.carriages = dict(carriages)
        self.primary_carriages = primary_carriages
        self.dc_carriages = dc_carriages
        ec_carriages = []
        for ecconfig in config.get_prefix_sections('extra_carriage '):
            ec_carriages.append(ExtraCarriage(ecconfig, carriages))
        for ec in ec_carriages:
            name = ec.get_name()
            if name in carriages:
                raise config.error("Redefinition of carriage %s" % name)
            carriages[name] = ec
        self.kin_steppers = self._load_steppers(config, carriages)
        self.all_carriages = carriages
        self._check_carriages_references(config.error)
        self._check_multi_mcu_homing(config.error)
    def _check_carriages_references(self, report_error):
        carriages = dict(self.all_carriages)
        for s in self.kin_steppers:
            for c in s.get_carriages():
                carriages.pop(c.get_name(), None)
        if carriages:
            raise report_error(
                    "Carriage(s) %s must be referenced by some "
                    "stepper(s)" % (", ".join(carriages),))
    def _check_multi_mcu_homing(self, report_error):
        for carriage in self.carriages.values():
            for es in carriage.get_rail().get_endstops():
                stepper_mcus = set([s.get_mcu() for s in es[0].get_steppers()
                                    if s in carriage.get_rail().get_steppers()])
                if len(stepper_mcus) > 1:
                    raise report_error("Multi-mcu homing not supported on"
                                       " multi-mcu shared carriage %s" % es[1])
    def _load_steppers(self, config, carriages):
        return [ks.KinematicStepper(c, carriages)
                for c in config.get_prefix_sections('stepper ')]
    def _update_extra_axes(self):
        test_junction_axes_set = []
        extra_axes = self.toolhead.get_extra_axes()
        for tc in self.track_junction_carriages \
                + self._default_junction_carriages:
            axes = [extra_axes.index(c) if c in extra_axes else c.get_axis()
                    for c in tc]
            axes.sort()
            if axes not in test_junction_axes_set:
                test_junction_axes_set.append(axes)
        self.test_junction_axes_set = test_junction_axes_set
    def get_steppers(self):
        return [s.get_stepper() for s in self.kin_steppers]
    def _get_kinematics_coeffs(self):
        matr = {s.get_name() : list(s.get_kin_coeffs())
                for s in self.kin_steppers}
        offs = {s.get_name() : [0.] * 3 for s in self.kin_steppers}
        if self.dc_module is None:
            return ([matr[s.get_name()] for s in self.kin_steppers],
                    [0. for s in self.kin_steppers])
        axes = [dc.get_axis() for dc in self.dc_carriages]
        orig_matr = copy.deepcopy(matr)
        for c in self.carriages.values():
            axis = c.get_axis()
            if axis in self.dc_module.get_axes():
                m, o = self.dc_module.get_transform(c.get_rail())
                for s in c.get_rail().get_steppers():
                    matr[s.get_name()][axis] *= m
                    offs[s.get_name()][axis] += o
        return ([matr[s.get_name()] for s in self.kin_steppers],
                [mathutil.matrix_dot(orig_matr[s.get_name()],
                                     offs[s.get_name()])
                 for s in self.kin_steppers])
    def _check_kinematics(self, report_error):
        matr, _ = self._get_kinematics_coeffs()
        det = mathutil.matrix_det(mat_mul(mat_transp(matr), matr))
        if abs(det) < 0.00001:
            raise report_error(
                    "Verify configured stepper(s) and their 'carriages' "
                    "specifications, the current configuration does not "
                    "allow independent movements of all printer axes.")
    def calc_position(self, stepper_positions):
        matr, offs = self._get_kinematics_coeffs()
        spos = [stepper_positions[s.get_name()] for s in self.kin_steppers]
        pinv = mat_pseudo_inverse(matr)
        pos = mat_mul([[sp-o for sp, o in zip(spos, offs)]], mat_transp(pinv))
        for i in range(3):
            if not any(pinv[i]):
                pos[0][i] = None
        return pos[0]
    def update_limits(self, i, range):
        self.limits[i] = range
    def set_position(self, newpos, homing_axes):
        for s in self.kin_steppers:
            s.set_position(newpos)
        for axis_name in homing_axes:
            axis = "xyz".index(axis_name)
            for c in self.carriages.values():
                if c.get_axis() == axis and c.is_active():
                    self.update_limits(axis, c.get_rail().get_range())
                    c.set_homed()
                    break
    def clear_homing_state(self, clear_axes):
        for c in self.carriages.values():
            axis = c.get_axis()
            axis_name = "xyz"[axis]
            if axis_name in clear_axes:
                c.set_homed(homed=False)
    def home_axis(self, homing_state, axis, rail):
        # Determine movement
        position_min, position_max = rail.get_range()
        hi = rail.get_homing_info()
        homepos = [None, None, None, None]
        homepos[axis] = hi.position_endstop
        forcepos = list(homepos)
        if hi.positive_dir:
            forcepos[axis] -= 1.5 * (hi.position_endstop - position_min)
        else:
            forcepos[axis] += 1.5 * (position_max - hi.position_endstop)
        # Perform homing
        homing_state.home_rails([rail], forcepos, homepos)
    def home(self, homing_state):
        self._check_kinematics(self.printer.command_error)
        primary_carriages = {c.get_axis(): c for c in self.primary_carriages}
        # Each axis is homed independently and in order
        for axis in homing_state.get_axes():
            if self.dc_module is not None and axis in self.dc_module.get_axes():
                self.dc_module.home(homing_state, axis)
            else:
                carriage = primary_carriages[axis]
                self.home_axis(homing_state, axis, carriage.get_rail())
    def check_move(self, move):
        limits = self.limits
        end_pos = move.end_pos
        for c in self.carriages.values():
            axis = c.get_axis()
            if not move.axes_d[axis]:
                continue
            if c.is_active() and not c.is_homed():
                raise move.move_error("Must home axis first")
            if (end_pos[axis] < self.limits[axis][0]
                or end_pos[axis] > self.limits[axis][1]):
                raise move.move_error()
            dc = c.get_dual_carriage()
            if dc is not None and dc in self.dc_toolhead.get_active_carriages():
                extra_axes = self.toolhead.get_extra_axes()
                # Always check with the dual carriage even if it is not moving
                dc.check_move(move, extra_axes.index(dc))
        if not move.axes_d[2]:
            # Normal XY move - use defaults
            return
        # Move with Z - update velocity and accel for slower Z axis
        z_ratio = move.move_d / abs(move.axes_d[2])
        move.limit_speed(
            self.max_z_velocity * z_ratio, self.max_z_accel * z_ratio)
    def calc_junction(self, prev_move, move):
        return min(move.calc_junction_v2(prev_move, axes, normalize=True)
                   for axes in self.test_junction_axes_set)
    def process_move(self, print_time, move):
        self.dc_toolhead.process_move(print_time, move)
    def activate_dc_direct_mode(self, carriage_name, carriage_pos):
        for dc in self.dc_carriages:
            if dc.get_name() == carriage_name:
                self.dc_toolhead.activate_direct_mode(dc, carriage_pos)
                return
        raise self.printer.command_error(
                "Cannot activate direct mode for carriage '%s'" % carriage_name)
    def deactivate_dc_direct_mode(self, carriage_name):
        for dc in self.dc_carriages:
            if dc.get_name() == carriage_name:
                self.dc_toolhead.deactivate_direct_mode(dc)
                return
        raise self.printer.command_error(
                "Cannot deactivate direct mode for carriage '%s'"
                % carriage_name)
    def get_status(self, eventtime):
        homed_axes = ["xyz"[c.get_axis()] for c in self.carriages.values()
                      if c.is_active() and c.is_homed()]
        ranges = [(min(c.get_rail().get_range()[0]
                       for c in self.carriages.values()
                       if c.get_axis() == axis),
                   max(c.get_rail().get_range()[1]
                       for c in self.carriages.values()
                       if c.get_axis() == axis))
                 for axis in range(3)]
        axes_min = gcode.Coord([r[0] for r in ranges])
        axes_max = gcode.Coord([r[1] for r in ranges])
        return {
            'homed_axes': "".join(homed_axes),
            'axis_minimum': axes_min,
            'axis_maximum': axes_max,
        }
    cmd_SET_TRACK_CARRIAGES_JUNCTION_help = "Set up cornering tracking"
    def cmd_SET_TRACK_CARRIAGES_JUNCTION(self, gcmd):
        self.toolhead.flush_step_generation()
        carriages_str = gcmd.get("CARRIAGES")
        carriages_list = carriages_str.split(',')
        if len(carriages_list) != 3:
            raise gcmd.error(
                    "CARRIAGES parameter must specify exactly 3 carriages")
        carriages = []
        for c_str in carriages_list:
            c_name = c_str.strip().lower()
            if c_name not in self.carriages:
                raise gcmd.error("Invalid carriage '%s' specified" % c_name)
            carriages.append(self.carriages[c_name])
        carriages.sort(key=lambda c: c.get_axis())
        for _, c_by_axis in itertools.groupby(carriages,
                                              key=lambda c: c.get_axis()):
            c_by_axis = list(c_by_axis)
            if len(c_by_axis) > 1:
                raise gcmd.error("Carriages %s share the same cartesian axis"
                                 % ', '.join(c.get_name() for c in c_by_axis))
        enable = gcmd.get_int("ENABLE", 1)
        if enable:
            if carriages not in self.track_junction_carriages and \
                    carriages not in self._default_junction_carriages:
                self.track_junction_carriages.append(carriages)
        else:
            if carriages in self.track_junction_carriages:
                self.track_junction_carriages.remove(carriages)
            else:
                raise gcmd.error(
                        "Cannot disable junction tracking for CARRIAGES=%s"
                        % carriages_str)
    cmd_SET_STEPPER_CARRIAGES_help = "Set stepper carriages"
    def cmd_SET_STEPPER_CARRIAGES(self, gcmd):
        self.toolhead.flush_step_generation()
        stepper_name = gcmd.get("STEPPER")
        steppers = [stepper for stepper in self.kin_steppers
                    if stepper.get_name() == stepper_name
                    or stepper.get_name(short=True) == stepper_name]
        if len(steppers) != 1:
            raise gcmd.error("Invalid STEPPER '%s' specified" % stepper_name)
        stepper = steppers[0]
        carriages_str = gcmd.get("CARRIAGES").lower()
        validate = not gcmd.get_int("DISABLE_CHECKS", 0)
        old_carriages = stepper.get_carriages()
        old_kin_coeffs = stepper.get_kin_coeffs()
        stepper.update_carriages(carriages_str, self.all_carriages, gcmd.error)
        new_carriages = stepper.get_carriages()
        if new_carriages != old_carriages:
            stepper.update_kin_coeffs(old_kin_coeffs)
            raise gcmd.error("SET_STEPPER_CARRIAGES cannot add or remove "
                             "carriages that the stepper controls")
        pos = self.toolhead.get_position()
        stepper.set_position(pos)
        if not validate:
            return
        if self.dc_module:
            dc_state = self.dc_module.save_dual_carriage_state()
            pcs = [dc.get_dual_carriage() for dc in self.dc_carriages]
            axes = [dc.get_axis() for dc in self.dc_carriages]
            for acs in itertools.product(*zip(pcs, self.dc_carriages)):
                for c in acs:
                    self.dc_module.get_dc_rail_wrapper(c.get_rail()).activate(
                            idex_modes.PRIMARY, pos)
                    self.dc_module.get_dc_rail_wrapper(
                            c.get_dual_carriage().get_rail()).inactivate(pos)
                self._check_kinematics(gcmd.error)
            self.dc_module.restore_dual_carriage_state(dc_state, move=0)
        else:
            self._check_kinematics(gcmd.error)

def load_kinematics(toolhead, config):
    return GenericCartesianKinematics(toolhead, config)
