#!/usr/bin/env python3
# A slicer post-processing script to enable fast IDEX tool swaps
#
# Copyright (C) 2025  Dmitry Butyugin <dmbutyugin@google.com>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

import datetime, math, optparse, re

PROCESSED_MARKER = '; Processed by fast tool swaps script'
EXTRUDER_SYNC_GCMD_TMPL = 'SYNC_EXTRUDER_MOTION EXTRUDER=%s MOTION_QUEUE=%s\n'
CARRIAGES_WIPE_PREPARE_SNIPPET = \
"""SET_DUAL_CARRIAGE CARRIAGE={next_carriage} MODE=DIRECT GCODE_AXIS={gcode_id}
G91
M83
M204 S{tool_wipe_accel:.0f}
G1 F{tool_wipe_feedrate:.0f}
"""
M104_GCMD_TMPL = 'M104 T%d S%.1f\n'

class GCodeState:
    def __init__(self, num_tools, other=None):
        self.num_tools = num_tools
        self.absolute_extrude = (True if other is None
                                 else other.absolute_extrude)
        self.print_time = 0. if other is None else other.print_time
        self.feed_rate = None if other is None else other.feed_rate
        self.acceleration = None if other is None else other.acceleration
        self.current_tool = 0 if other is None else other.current_tool
        self.extruder_temps = ([0] * num_tools if other is None
                               else other.extruder_temps)
        self.position = [0.] * 4 if other is None else list(other.position)
    def get_updated_from(self, gcmd):
        new = GCodeState(self.num_tools, other=self)
        gcmd_type = gcmd.get_type()
        match gcmd.get_type():
            case 'T':
                new.current_tool = int(gcmd.get_cmd()[1:])
                if new.current_tool < 0 or new.current_tool >= self.num_tools:
                    raise RuntimeError("Malformed '%s' encountered"
                                       % gcmd.get_line().strip())
            case 'G1':
                f = gcmd.get_param('F')
                if f is not None:
                    new.feed_rate = float(f)
                for i, a in enumerate("XYZ"):
                    new.position[i] = gcmd.get_param(a, float, self.position[i])
                if new.absolute_extrude:
                    new.position[3] = gcmd.get_param('E', float,
                                                     self.position[3])
                else:
                    new.position[3] += gcmd.get_param('E', float, 0.)
            case 'M104':
                tool = gcmd.get_param('T', int, self.current_tool)
                temp = gcmd.get_param('S', float, 0.)
                if tool < 0 or tool >= self.num_tools:
                    raise RuntimeError("Malformed '%s' encountered"
                                       % gcmd.get_line().strip())
                new.extruder_temps[tool] = temp
            case 'M204':
                new.acceleration = gcmd.get_param('S', float, self.acceleration)
            case 'SVL':
                new.acceleration = gcmd.get_param('ACCEL', float,
                                                  self.acceleration)
                new.feed_rate = gcmd.get_param('VELOCITY', float,
                                               self.feed_rate / 60.) * 60.
            case 'M82':
                new.absolute_extrude = gcmd.get_cmd() == 'M82'
            case 'G92':
                new.position[3] = gcmd.get_param('E', float, self.position[3])
        return new
    def update_print_time(self, gcmd):
        self.print_time += gcmd.estimate_print_time()

class GCodeCommand:
    args_r = re.compile('([A-Z_]+|[A-Z*])')
    def __init__(self, line, prev_gcode_state):
        self.line = line
        line = line.strip().upper()
        cpos = line.find(';')
        hpos = line.find('#')
        comment_pos = (max(cpos, hpos) if min(cpos, hpos) < 0
                       else min(cpos, hpos))
        if not line or comment_pos == 0:
            self.cmd = line
            self.params = {}
            self.gcode_state = prev_gcode_state
            return
        if comment_pos > 0:
            line = line[:comment_pos]
        parts = self.args_r.split(line)
        cmd = ''.join(parts[:3]).strip()
        # Build gcode "params"
        if len(cmd) > 1 and cmd[0].isalpha() and cmd[1].isdigit():
            params = { parts[i]: parts[i+1].strip()
                       for i in range(3, len(parts), 2) }
        else:
            parts = [s.split('=', 1)
                     for s in ''.join(parts[3:]).strip().split(' ') if s]
            params = { k: v for k, v in parts }
        self.cmd = cmd
        self.params = params
        self.prev_gcode_state = prev_gcode_state
        self.gcode_state = prev_gcode_state.get_updated_from(self)
        self.move_d = [
                self.gcode_state.position[i] - prev_gcode_state.position[i]
                for i in range(4)]
        if self.get_type() == 'G1' and not self.gcode_state.absolute_extrude:
            self.move_d[3] = self.get_param('E', float, 0.)
        self.gcode_state.update_print_time(self)
    def get_type(self):
        cmd = self.cmd
        if cmd in ('G0', 'G1'):
            return 'G1'
        if cmd in ('T0', 'T1'):
            return 'T'
        if cmd.startswith(';'):
            return 'C'
        if cmd in ('M104', 'M109'):
            return 'M104'
        if cmd == 'SET_VELOCITY_LIMIT':
            return 'SVL'
        if cmd == 'M204':
            return 'M204'
        if cmd in ('M82', 'M83'):
            return 'M82'
        if cmd == 'G92':
            return 'G92'
        return 'O'
    def get_cmd(self):
        return self.cmd
    def get_param(self, param_name, param_type=str, default=None):
        if param_name not in self.params:
            return default
        return param_type(self.params[param_name])
    def get_line(self):
        return self.line
    def get_gcode_state(self):
        return self.gcode_state
    def estimate_print_time(self, feed_rate=None):
        feed_rate = self.gcode_state.feed_rate
        if self.get_type() != 'G1' or not feed_rate:
            return 0.
        length = self.get_length()
        speed = feed_rate / 60.
        acceleration = self.gcode_state.acceleration
        t0 = length / speed
        if not acceleration or not t0:
            return t0
        t1 = t0 * (1. +  min(0.25 * speed / acceleration / t0, 1.0))
        return t1
    def get_length(self):
        if self.get_type() != 'G1':
            return 0.
        if any(self.move_d[:3]):
            return math.sqrt(sum(d**2 for d in self.move_d))
        return abs(self.get_extrusion_length())
    def get_extrusion_length(self):
        if self.get_type() != 'G1':
            return 0.
        return self.move_d[3]

class GCodeProcessor:
    def __init__(self, num_tools, optparser, options):
        self.num_tools = num_tools
        self.parse_per_tool_float_param(optparser, options,
                                        'inactive_tool_temp_reduction',
                                        default=1., min_val=0., max_val=1.)
        self.parse_per_tool_float_param(optparser, options, 'tool_shutoff_time',
                                        default=-1, min_val=0., max_val=3600.)
        self.enable_temp_control = [
                ittr < 1. or tst >= 0.
                for ittr, tst in zip(self.inactive_tool_temp_reduction,
                                     self.tool_shutoff_time)]
        self.parse_per_tool_float_param(
                optparser, options, 'min_preheat_time',
                default=(None if min(self.inactive_tool_temp_reduction) < 1.
                         else 0.),
                min_val=0., max_val=60.)
        self.parse_per_tool_float_param(
                optparser, options, 'min_purge_length',
                default=(None if min(self.inactive_tool_temp_reduction) < 1.
                         else -1),
                min_val=0., max_val=50.)
        self.parse_per_tool_float_param(
                optparser, options, 'long_preheat_time',
                default=(None if (any(tst >= 0. and ittr >= 1.
                                      for ittr, tst in zip(
                                          self.inactive_tool_temp_reduction,
                                          self.tool_shutoff_time)))
                         else -1),
                min_val=0., max_val=300.)
        self.parse_per_tool_float_param(
                optparser, options, 'long_preheat_purge_length',
                default=(None if (any(tst >= 0.
                                      for tst in self.tool_shutoff_time)
                                  and not any(mpl >= 0.
                                              for mpl in self.min_purge_length))
                         else -1),
                min_val=0., max_val=50.)
        for i in range(num_tools):
            if self.long_preheat_time[i] < 0 \
                    and self.inactive_tool_temp_reduction[i] < 1.:
                self.long_preheat_time[i] = (
                        self.min_preheat_time[i] /
                        (1. - self.inactive_tool_temp_reduction[i]))
                optparser.error(
            if self.long_preheat_time[i] >= 0 and \
                    self.long_preheat_time[i] < self.min_preheat_time[i]:
                        "Incorrect long_preheat_time=%.3f provided"
                        % self.long_preheat_time[i])
            if self.long_preheat_purge_length[i] < 0:
                self.long_preheat_purge_length[i] = self.min_purge_length[i]
            elif self.long_preheat_purge_length[i] < self.min_purge_length[i]:
                optparser.error(
                        "Incorrect long_preheat_purge_length=%.3f provided"
                        % self.long_preheat_purge_length[i])
        tools = []
        for i in range(num_tools):
            t_def = vars(options)['t%d' % i]
            t_parts = t_def.split(':')
            if len(t_parts) != 3 or t_parts[2] not in ('+', '-'):
                optparser.error("Incorrect t%d=%s provided" % (i, t_def))
            tools.append(t_parts)
        self.carriages = tuple(tools[i][0] for i in range(num_tools))
        self.extruders = tuple(tools[i][1] for i in range(num_tools))
        self.wipe_dirs = tuple(tools[i][2] for i in range(num_tools))
        self.extra_gcode_axis = options.extra_gcode_axis.upper()
        if len(self.extra_gcode_axis) != 1 or self.extra_gcode_axis in "XYZEFN":
            optparser.error("Incorrect extra_gcode_axis=%s provided"
                            % options.extra_gcode_axis)
        self.parse_per_tool_float_param(optparser, options, 'tool_wipe_dist',
                                        min_val=0., max_val=100.)
        self.wipe_enabled = [twd > 0. for twd in self.tool_wipe_dist]
        self.parse_per_tool_float_param(optparser, options, 'tool_wipe_speed',
                                        default=(None if any(self.wipe_enabled)
                                                 else 0.),
                                        min_val=1., max_val=10000.)
        self.parse_per_tool_float_param(optparser, options, 'tool_wipe_accel',
                                        default=(None if any(self.wipe_enabled)
                                                 else 0.),
                                        min_val=100., max_val=100000.)
        self.parse_per_tool_float_param(optparser, options, 'tool_wipe_z_hop',
                                        default=0., min_val=0., max_val=10.)
        self.parse_per_tool_float_param(optparser, options, 'tool_swap_speed',
                                        default=(0. if any(self.tool_wipe_speed)
                                                 else None),
                                        min_val=0., max_val=10000.)
        self.parse_per_tool_float_param(optparser, options, 'tool_swap_accel',
                                        default=(0. if any(self.tool_wipe_accel)
                                                 else None),
                                        min_val=0., max_val=100000.)
        for i, tws in enumerate(self.tool_wipe_speed):
            if not self.tool_swap_speed[i]:
                self.tool_swap_speed[i] = tws
        for i, twa in enumerate(self.tool_wipe_accel):
            if not self.tool_swap_accel[i]:
                self.tool_swap_accel[i] = twa

        for i in range(num_tools):
            if self.tool_shutoff_time[i] < 0:
                if self.enable_temp_control[i]:
                    self.tool_shutoff_time[i] = (
                            10. * self.min_preheat_time[i] /
                            (1. - self.inactive_tool_temp_reduction[i]))
                else:
                    self.tool_shutoff_time[i] = 0.
        self.tool_fast_change_gcode = options.tool_fast_change_gcode
        if not self.tool_fast_change_gcode:
            optparser.error("--tool_fast_change_gcode must be provided")
        if not "{next_tool}" in self.tool_fast_change_gcode:
            optparser.error("Malformed tool_fast_change_gcode='%s' provided"
                            % options.tool_fast_change_gcode)
        self.cur_tool = -1
    def parse_per_tool_float_param(self, optparser, options, name, default=None,
                                   min_val=None, max_val=None):
        value_str = vars(options)[name]
        if not value_str:
            if default is not None:
                vars(self)[name] = [default] * self.num_tools
                return
            else:
                optparser.error("Parameter %s must be provided" % name)
        elif ',' in value_str:
            try:
                values = [float(v.strip()) for v in value_str.split(',')]
            except ValueError:
                optparser.error("Invalid format for %s='%s'" % (name,
                                                                value_str))
        else:
            try:
                values = [float(value_str.strip())]
            except ValueError:
                optparser.error("Invalid format for %s='%s'" % (name,
                                                                value_str))
        for value in values:
            if min_val is not None and value < min_val:
                optparser.error("Too small %s=%f provided" % (name, value))
            if max_val is not None and value > max_val:
                optparser.error("Too large %s=%f provided" % (name, value))
        if len(values) == 1:
            values = values * self.num_tools
        elif len(values) != self.num_tools:
            optparser.error("Invalid number of values provided in %s='%s'"
                            % (name, value_str))
        vars(self)[name] = values
    def consume_next_line(self):
        while True:
            line = next(self.lines, None)
            if line is None:
                return None
            gcmd = GCodeCommand(line, self.gcode_state)
            self.gcode_state = gcmd.get_gcode_state()
            # Check if a GCode command should be skipped
            if gcmd.get_type() == 'M104' and \
                    self.enable_temp_control[self.gcode_state.current_tool]:
                T = gcmd.get_param('T', int, self.gcode_state.current_tool)
                if T != self.gcode_state.current_tool:
                    continue
            self.buffer.append(gcmd)
            return gcmd
    def process(self, lines):
        self.buffer = []
        self.cur_tool = -1
        self.tool_heated = [False, False]
        self.lines = lines
        self.gcode_state = GCodeState(num_tools=self.num_tools)
        yield PROCESSED_MARKER + datetime.datetime.now(datetime.UTC).strftime(
                ' on %Y-%m-%d at %H:%M:%S UTC\n')
        while True:
            gcmd = self.consume_next_line()
            if gcmd is None:
                break
            if gcmd.get_type() != 'T':
                continue
            for l in self.handle_tool_swap():
                yield l
        if self.buffer:
            cur_tool = self.buffer[0].get_gcode_state().current_tool
            for i, e in enumerate(
                    self.buffer[0].get_gcode_state().extruder_temps):
                if i != cur_tool and e and self.enable_temp_control[i]:
                    yield M104_GCMD_TMPL % (i, 0.)
                    self.tool_heated[i] = False
        for gcmd in self.buffer:
            yield gcmd.get_line()
    def handle_tool_swap(self):
        first_gcode_state = self.buffer[0].get_gcode_state()
        tool_swap_gcode_state = self.buffer[-1].get_gcode_state()
        next_tool = tool_swap_gcode_state.current_tool
        if self.cur_tool in (-1, next_tool):
            for gcmd in self.buffer[:-1]:
                yield gcmd.get_line()
            if self.cur_tool == -1:
                yield self.buffer[-1].get_line()
            del self.buffer[:]
            self.cur_tool = next_tool
            self.tool_heated[self.cur_tool] = True
            return
        min_purge_length = self.min_purge_length[next_tool]
        min_preheat_time = self.min_preheat_time[next_tool]
        tool_swap_print_time = (tool_swap_gcode_state.print_time
                                - first_gcode_state.print_time)
        inactive_tool_temp_reduction = \
                self.inactive_tool_temp_reduction[next_tool]
        long_preheat = self.enable_temp_control[next_tool] and (
                not self.tool_heated[next_tool]
                or (self.tool_shutoff_time[next_tool] and
                    tool_swap_print_time > max(
                        self.long_preheat_time[next_tool],
                        self.tool_shutoff_time[next_tool])))
        if long_preheat:
            min_preheat_time = self.long_preheat_time[next_tool]
            min_purge_length = self.long_preheat_purge_length[next_tool]
            if self.tool_heated[next_tool]:
                yield M104_GCMD_TMPL % (next_tool, 0.)
        min_purge_length = max(min_purge_length, 0)
        wipe_start_ind = wipe_end_ind = preheat_ind = purge_ind = None
        extrude_length = print_time = 0.
        tool_change_ind = len(self.buffer) - 1
        if not self.wipe_enabled[next_tool]:
            wipe_start_ind = wipe_end_ind = tool_change_ind - 1
        for i in range(tool_change_ind, -1, -1):
            gcmd = self.buffer[i]
            if wipe_start_ind is not None:
                extrude_length += gcmd.get_extrusion_length()
            if wipe_end_ind is None and gcmd.get_cmd().startswith(';WIPE_END'):
                wipe_end_ind = i
            elif (wipe_end_ind is not None and wipe_start_ind is None
                  and gcmd.get_cmd().startswith(';WIPE_START')):
                wipe_start_ind = i
            if (wipe_start_ind is not None and purge_ind is None
                  and extrude_length >= min_purge_length):
                if not min_purge_length:
                    purge_ind = i
                elif gcmd.get_length() and not gcmd.get_extrusion_length():
                    # Make sure to not insert the extruder sync command midway
                    # of the printing seqeunce of moves as this will cause a zit
                    purge_ind = i + 1
            if purge_ind is not None:
                print_time += gcmd.estimate_print_time()
            if (purge_ind is not None and preheat_ind is None
                  and print_time >= min_preheat_time):
                preheat_ind = i
                break
        if wipe_start_ind is None:
            raise RuntimeError(
                    "No annotated WIPE found before the tool change")
        preheat_ind = preheat_ind or 0
        purge_ind = purge_ind or 0
        need_z = False
        while True:
            gcmd = self.consume_next_line()
            if gcmd is None:
                break
            if gcmd.get_type() == 'G1' and any(a in gcmd.params for a in "XY"):
                if 'Z' in gcmd.params:
                    need_z = True
                break
        next_temp = self.buffer[-1].get_gcode_state().extruder_temps[next_tool]
        target_pos = self.buffer[-1].get_gcode_state().position[:3]
        for i in range(preheat_ind):
            yield self.buffer[i].get_line()
        if self.enable_temp_control[next_tool]:
            if not next_temp:
                raise RuntimeError(
                        "Must print with T%d, but its temperature is not set"
                        % next_tool)
            yield M104_GCMD_TMPL % (next_tool, next_temp)
        self.tool_heated[next_tool] = True
        for i in range(preheat_ind, purge_ind):
            yield self.buffer[i].get_line()
        actual_purge_length = 0.
        if min_purge_length:
            yield EXTRUDER_SYNC_GCMD_TMPL % (self.extruders[next_tool],
                                             self.extruders[self.cur_tool])
        for i in range(purge_ind, wipe_start_ind+1):
            actual_purge_length += self.buffer[i].get_extrusion_length()
            yield self.buffer[i].get_line()
        if inactive_tool_temp_reduction < 1.:
            yield M104_GCMD_TMPL % (
                    self.cur_tool,
                    tool_swap_gcode_state.extruder_temps[self.cur_tool]
                    * inactive_tool_temp_reduction)
            self.tool_heated[self.cur_tool] = True
        yield CARRIAGES_WIPE_PREPARE_SNIPPET.format(**{
            'next_carriage': self.carriages[next_tool],
            'prev_carriage': self.carriages[self.cur_tool],
            'gcode_id': self.extra_gcode_axis,
            'tool_wipe_feedrate': 60. * (self.tool_wipe_speed[next_tool] or
                                         self.tool_swap_speed[next_tool]),
            'tool_wipe_accel': (self.tool_wipe_accel[next_tool] or
                                self.tool_swap_accel[next_tool])
            })
        if any(gcmd.get_extrusion_length() > 1e-7
               for gcmd in self.buffer[wipe_start_ind+1:tool_change_ind]):
            raise RuntimeError("More extrusion after the last tool wipe"
                               + " before the tool change")
        half_retract = -0.5 * sum(
                gcmd.get_extrusion_length()
                for gcmd in self.buffer[wipe_start_ind+1:tool_change_ind])
        executed_retract = 0.
        wipe_dir = 1.0 if self.wipe_dirs[next_tool] == '+' else -1.0
        wipe_dist = self.tool_wipe_dist[next_tool]
        for i in range(wipe_start_ind+1, tool_change_ind):
            gcmd = self.buffer[i]
            if gcmd.get_type() != 'G1':
                yield self.buffer[i].get_line()
                continue
            if not any(gcmd.move_d):
                continue
            cur_retract = -gcmd.get_extrusion_length()
            while cur_retract > 0.:
                if executed_retract + cur_retract > half_retract:
                    retract = half_retract - executed_retract
                    move_d = [d * retract / -gcmd.move_d[3]
                              for d in gcmd.move_d]
                    executed_retract = 0.
                    cur_retract -= retract
                    if cur_retract < 0.000001:
                        cur_retract = 0.
                    next_wipe_dir = -wipe_dir
                else:
                    move_d = [d * cur_retract / -gcmd.move_d[3]
                              for d in gcmd.move_d]
                    executed_retract += cur_retract
                    cur_retract = 0.
                    next_wipe_dir = wipe_dir
                coord_str = ' '.join('%s%.6f' % (a, d)
                                     for a, d in zip("XYZE", move_d) if d)
                extra_axis_d = -wipe_dir * wipe_dist * move_d[3] / half_retract
                coord_str += ' %s%.6f' % (self.extra_gcode_axis, extra_axis_d)
                tool_wipe_z_hop = self.tool_wipe_z_hop[next_tool]
                if not gcmd.move_d[2] and tool_wipe_z_hop:
                    z_hop_d = -0.5 * tool_wipe_z_hop * move_d[3] / half_retract
                    coord_str += ' Z%.6f' % z_hop_d
                yield 'G1 %s\n' % coord_str
                wipe_dir = next_wipe_dir
        if min_purge_length:
            yield EXTRUDER_SYNC_GCMD_TMPL % (self.extruders[next_tool],
                                             self.extruders[next_tool])
        yield 'G90\n'
        if self.wipe_enabled[next_tool]:
            yield 'M204 S%.0f\n' % self.tool_swap_accel[next_tool]
            yield 'G1 F%.0f\n' % (60. * self.tool_swap_speed[next_tool])
        yield self.tool_fast_change_gcode.format(**{
            'next_tool': next_tool,
            'xpos': '%.6f' % target_pos[0],
            'ypos': '%.6f' % target_pos[1],
            }) + '\n'
        if need_z:
            yield 'G1 Z%.6f\n' % target_pos[2]
        yield 'G1 F%.0f\n' % tool_swap_gcode_state.feed_rate
        yield 'M204 S%.f\n' % tool_swap_gcode_state.acceleration
        if tool_swap_gcode_state.absolute_extrude:
            yield 'M82\n'
        self.cur_tool = next_tool
        del self.buffer[:tool_change_ind+1]

def main():
    # Parse command-line arguments
    usage = "%prog [options] <input>"
    opts = optparse.OptionParser(usage)
    opts.add_option("-o", "--output", type="string", dest="output",
                    default=None, help="Filename of the output. If not set,"
                    + " the script edits the input file in-place.")
    opts.add_option("--inactive_tool_temp_reduction", type="string",
                    default=None, help="If set, defines a ratio [0, 1] by how"
                    + " much to reduce the temperature during tool inactivity, "
                    + " either a single number or a comma-separated list"
                    + " of values, one per tool. If unset or set to 1, disables"
                    + " temperature control of the tools by the script.")
    opts.add_option("-t", "--min_preheat_time", type="string",
                    dest="min_preheat_time", default=None,
                    help="Minimum extruder preheat time before activation (s),"
                    + " either a single number or a comma-separated list"
                    + " of values, one per tool. Must be provided if the"
                    + " temperature control is enabled via"
                    + " --inactive_tool_temp_reduction option.")
    opts.add_option("-l", "--min_purge_length", type="string",
                    dest="min_purge_length", default=None,
                    help="Minimum extruder purge length before activation (mm),"
                    + " either a single number or a comma-separated list"
                    + " of values, one per tool. Must be provided if"
                    + " --inactive_tool_temp_reduction is configured.")
    opts.add_option("--tool_shutoff_time", type="string", default=None,
                    help="A time of tool inactivity to turn it off (s). If not"
                    + " specified, will be inferred from --min_preheat_time."
                    + " Set this option to 0 to disable this feature. The tool"
                    + " shutoff can be enalbed even if temperature reduction"
                    + " is not configured, which still requires either of"
                    + " --min_purge_length or --long_preheat_purge_length"
                    + " options to be specified to set a purge length. Can be"
                    + " either a single number or a comma-separated list"
                    + " of values, one per tool.")
    opts.add_option("--long_preheat_time", type="string", default=None,
                    help="Min extruder preheat tiem after long inactivity (s),"
                    + " either a single number or a comma-separated list"
                    + " of values, one per tool. If not set, will be inferred"
                    + " from --min_preheat_time option.")
    opts.add_option("--long_preheat_purge_length", type="string", default=None,
                    help="Min extruder purge length after long inactivity (mm),"
                    + " either a single number or a comma-separated list"
                    + " of values, one per tool. If not set, defaults to the"
                    + " values of --min_purge_length option.")
    opts.add_option("--t0", type="string", dest="t0", default="xc:extruder:+",
                    help="T0 definition in the format 'carriage:extruder:[+|-]'"
                    + " with 'carriage' being a name of generic cartesian"
                    + " carriage name of T0, 'extruder' being a name of the"
                    + " corresponding extruder, and a sign indicating the"
                    + " direction of the wipe move along the carriage axis."
                    + " This option must be provided.")
    opts.add_option("--t1", type="string", dest="t1", default="uc:extruder1:-",
                    help="T1 definition in the format 'carriage:extruder:[+|-]'"
                    + " with 'carriage' being a name of generic cartesian"
                    + " carriage name of T1, 'extruder' being a name of the"
                    + " corresponding extruder, and a sign indicating the"
                    + " direction of the wipe move along the carriage axis."
                    + " This option must be provided.")
    opts.add_option("--extra_gcode_axis", type="string", default="U",
                    help="Extra GCode axis to use for tool swaps."
                    + " If not set, defaults to 'U' GCode axis.")
    opts.add_option("--tool_wipe_dist", type="string", default=None,
                    help="Tool wipe distance during tool swaps (mm)."
                    + " This option defines the distance the tool will be"
                    + " moved from its parked position forth and back during"
                    + " tool activation. Should either be a single number or"
                    + " a comma-separated list of values, one per tool."
                    + " This option must be provided.")
    opts.add_option("--tool_wipe_speed", type="string", default=None,
                    help="Tool wipe speed during tool swaps (mm/s),"
                    + " either a single number or a comma-separated list of"
                    + " values, one per tool. This option must be provided.")
    opts.add_option("--tool_wipe_accel", type="string", default=None,
                    help="Tool wipe accel during tool swaps (mm/s^2),"
                    + " either a single number or a comma-separated list of"
                    + " values, one per tool. This option must be provided.")
    opts.add_option("--tool_fast_change_gcode", type="string", default=None,
                    help="A GCode template to call for the tool swap. Available"
                    + " template substitutions are {next_tool}, and {xpos} and"
                    + " {ypos} for the target tool position. The {next_tool}"
                    + " substitution must be used by the template.")
    opts.add_option("--tool_swap_speed", type="string", default=None,
                    help="Tool swap speed during tool swaps (mm/s) when"
                    + " the tool moves to its designated position (only if"
                    + " TOOL_FAST_CHANGE_GCODE uses {xpos} and {ypos}"
                    + " substitutions, otherwise has no effect and regular"
                    + " travel move speed from the slicer is used),"
                    + " either a single number or a comma-separated list of"
                    + " values, one per tool. If unset, defaults to"
                    + " --tool_wipe_speed option.")
    opts.add_option("--tool_swap_accel", type="string", default=None,
                    help="Tool swap accel during tool swaps (mm/s^2) when"
                    + " the tool moves to its designated position (only if"
                    + " TOOL_FAST_CHANGE_GCODE uses {xpos} and {ypos}"
                    + " substitutions, otherwise has no effect and regular"
                    + " travel move acceleration from the slicer is used),"
                    + " either a single number or a comma-separated list of"
                    + " values, one per tool. If unset, defaults to"
                    + " --tool_wipe_accel option.")
    opts.add_option("--tool_wipe_z_hop", type="string", default=None,
                    help="Tool wipe z-hop (mm), either a single number or"
                    + " a comma-separated list of values, one per tool."
                    + " When not set, explicit z-hops are disabled.")
    options, args = opts.parse_args()
    if len(args) != 1:
        opts.error("Incorrect number of arguments")

    gcode_processor = GCodeProcessor(num_tools=2,
                                     optparser=opts, options=options)
    with open(args[0], encoding='utf-8') as fi:
        input_lines = [l for l in fi]
    if any(True for line in input_lines
           if line.startswith(PROCESSED_MARKER)):
        opts.error('The input file was already processed by the fast tool'
                   ' swaps script once. You must use the original file'
                   ' for processing, or slice GCode anew.')
    with open(options.output or args[0], mode='wt', encoding='utf-8') as fo:
        try:
            for out_line in gcode_processor.process(iter(input_lines)):
                fo.write(out_line)
        except Exception as e:
            opts.error(str(e))

if __name__ == '__main__':
    main()
