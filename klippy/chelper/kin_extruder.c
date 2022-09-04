// Extruder stepper pulse time generation
//
// Copyright (C) 2018-2019  Kevin O'Connor <kevin@koconnor.net>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <math.h> // sqrt
#include <stddef.h> // offsetof
#include <stdlib.h> // malloc
#include <string.h> // memset
#include "compiler.h" // __visible
#include "itersolve.h" // struct stepper_kinematics
#include "kin_shaper.h" // struct shaper_pulses
#include "pyhelper.h" // errorf
#include "trapq.h" // move_get_distance

struct extruder_stepper {
    struct stepper_kinematics sk;
    struct shaper_pulses sp[3];
    double pressure_advance, half_smooth_time, inv_half_smooth_time2;
};

// Without pressure advance, the extruder stepper position is:
//     extruder_position(t) = nominal_position(t)
// When pressure advance is enabled, additional filament is pushed
// into the extruder during acceleration (and retracted during
// deceleration). The formula is:
//     pa_position(t) = (nominal_position(t)
//                       + pressure_advance * nominal_velocity(t))
// Which is then "smoothed" using a weighted average:
//     smooth_position(t) = (
//         definitive_integral(pa_position(x) * (smooth_time/2 - abs(t-x)) * dx,
//                             from=t-smooth_time/2, to=t+smooth_time/2)
//         / ((smooth_time/2)**2))

// Calculate the definitive integral of the motion formula:
//   position(t) = base + t * (start_v + t * half_accel)
static double
extruder_integrate(double base, double start_v, double half_accel
                   , double start, double end)
{
    double half_v = .5 * start_v, sixth_a = (1. / 3.) * half_accel;
    double si = start * (base + start * (half_v + start * sixth_a));
    double ei = end * (base + end * (half_v + end * sixth_a));
    return ei - si;
}

// Calculate the definitive integral of time weighted position:
//   weighted_position(t) = t * (base + t * (start_v + t * half_accel))
static double
extruder_integrate_time(double base, double start_v, double half_accel
                        , double start, double end)
{
    double half_b = .5 * base, third_v = (1. / 3.) * start_v;
    double eighth_a = .25 * half_accel;
    double si = start * start * (half_b + start * (third_v + start * eighth_a));
    double ei = end * end * (half_b + end * (third_v + end * eighth_a));
    return ei - si;
}

// Calculate the definitive integral of extruder for a given move
static double
pa_move_integrate(struct move *m, int axis, double move_dir_r
                  , double pressure_advance, double base
                  , double start, double end, double time_offset)
{
    if (start < 0.)
        start = 0.;
    if (end > m->move_t)
        end = m->move_t;
    // Calculate base position and velocity with pressure advance
    int can_pressure_advance = m->axes_r.x > 0. || m->axes_r.y > 0.;
    if (!can_pressure_advance)
        pressure_advance = 0.;
    double axis_r = move_dir_r * m->axes_r.axis[axis - 'x'];
    double start_v = m->start_v * axis_r;
    double ha = m->half_accel * axis_r;
    base += pressure_advance * start_v;
    start_v += pressure_advance * 2. * ha;
    // Calculate definitive integral
    double iext = extruder_integrate(base, start_v, ha, start, end);
    double wgt_ext = extruder_integrate_time(base, start_v, ha, start, end);
    return wgt_ext - time_offset * iext;
}

// Calculate the definitive integral of the extruder over a range of moves
static double
pa_range_integrate(struct move *m, int axis, double move_dir_r, double move_time
                   , double pressure_advance, double hst)
{
    while (unlikely(move_time < 0.)) {
        m = list_prev_entry(m, node);
        move_time += m->move_t;
    }
    while (unlikely(move_time > m->move_t)) {
        move_time -= m->move_t;
        m = list_next_entry(m, node);
    }
    // Calculate integral for the current move
    double res = 0., start = move_time - hst, end = move_time + hst;
    double start_base = m->start_pos.axis[axis - 'x'];
    res += pa_move_integrate(m, axis, move_dir_r, pressure_advance, 0.
                             , start, move_time, start);
    res -= pa_move_integrate(m, axis, move_dir_r, pressure_advance, 0.
                             , move_time, end, end);
    // Integrate over previous moves
    struct move *prev = m;
    while (unlikely(start < 0.)) {
        prev = list_prev_entry(prev, node);
        start += prev->move_t;
        double base = prev->start_pos.axis[axis - 'x'] - start_base;
        res += pa_move_integrate(prev, axis, move_dir_r, pressure_advance, base
                                 , start, prev->move_t, start);
    }
    // Integrate over future moves
    while (unlikely(end > m->move_t)) {
        end -= m->move_t;
        m = list_next_entry(m, node);
        double base = m->start_pos.axis[axis - 'x'] - start_base;
        res -= pa_move_integrate(m, axis, move_dir_r, pressure_advance, base
                                 , 0., end, end);
    }
    return res + start_base * hst * hst;
}

/****************************************************************
 * Extruder per-axis position calculation via shaper convolution
 ****************************************************************/

static inline double
get_axis_position(struct move *m, int axis, double move_dir_r, double move_time)
{
    double axis_r = move_dir_r * m->axes_r.axis[axis - 'x'];
    double start_pos = m->start_pos.axis[axis - 'x'];
    double move_dist = move_get_distance(m, move_time);
    return start_pos + axis_r * move_dist;
}

static inline double
get_axis_position_across_moves(struct move *m, int axis, double move_dir_r
                               , double time)
{
    while (likely(time < 0.)) {
        m = list_prev_entry(m, node);
        time += m->move_t;
    }
    while (likely(time > m->move_t)) {
        time -= m->move_t;
        m = list_next_entry(m, node);
    }
    return get_axis_position(m, axis, move_dir_r, time);
}

// Calculate the position from the convolution of the shaper with input signal
static inline double
shaper_calc_position(struct move *m, int axis, double move_dir_r
                              , double move_time, struct shaper_pulses *sp)
{
    int num_pulses = sp->num_pulses, i;
    if (!num_pulses) {
        return get_axis_position(m, axis, move_dir_r, move_time);
    }
    double res = 0.;
    for (i = 0; i < num_pulses; ++i) {
        double t = sp->pulses[i].t, a = sp->pulses[i].a;
        res += a * get_axis_position_across_moves(m, axis, move_dir_r
                                                  , move_time + t);
    }
    return res;
}

/****************************************************************
 * Toolhead velocity direction calculation via shaper convolution
 ****************************************************************/

static inline double
get_axis_velocity(struct move *m, int axis, double move_time)
{
    double axis_r = m->axes_r.axis[axis - 'x'];
    return axis_r * (m->start_v + 2. * m->half_accel * move_time);
}

static inline double
get_axis_velocity_across_moves(struct move *m, int axis, double time)
{
    while (likely(time < 0.)) {
        m = list_prev_entry(m, node);
        time += m->move_t;
    }
    while (likely(time > m->move_t)) {
        time -= m->move_t;
        m = list_next_entry(m, node);
    }
    return get_axis_velocity(m, axis, time);
}

// Calculate the velocity from the convolution of the shaper with input signal
static inline double
calc_velocity(struct move *m, int axis, double move_time
              , struct shaper_pulses *sp)
{
    int num_pulses = sp->num_pulses, i;
    if (!num_pulses) {
        return get_axis_velocity(m, axis, move_time);
    }
    double res = 0.;
    for (i = 0; i < num_pulses; ++i) {
        double t = sp->pulses[i].t, a = sp->pulses[i].a;
        res += a * get_axis_velocity_across_moves(m, axis, move_time + t);
    }
    return res;
}

static struct coord
get_move_dir(struct move *m, double move_time, struct extruder_stepper *es)
{
    struct coord move_dir;
    int i;
    double norm2 = 0.;
    for (i = 0; i < 3; ++i) {
        int axis = 'x' + i;
        struct shaper_pulses* sp = &es->sp[i];
        move_dir.axis[i] = calc_velocity(m, axis, move_time, sp);
        norm2 += move_dir.axis[i] * move_dir.axis[i];
    }
    double inv_norm = 1. / sqrt(norm2);
    for (i = 0; i < 3; ++i) move_dir.axis[i] *= inv_norm;
    return move_dir;
}

/****************************************************************
 * Extruder PA calculation via shaper convolution
 ****************************************************************/

static double
shaper_pa_range_integrate(struct move *m, int axis, double move_dir_r
                          , double move_time, double pressure_advance
                          , double hst, struct shaper_pulses *sp)
{
    double res = 0.;
    int num_pulses = sp->num_pulses, i;
    if (!num_pulses) {
        return pa_range_integrate(m, axis, move_dir_r, move_time
                                  , pressure_advance, hst);
    }
    for (i = 0; i < num_pulses; ++i) {
        double t = sp->pulses[i].t, a = sp->pulses[i].a;
        res += a * pa_range_integrate(m, axis, move_dir_r, move_time + t
                                      , pressure_advance, hst);
    }
    return res;
}

static double
extruder_calc_position(struct stepper_kinematics *sk, struct move *m
                       , double move_time)
{
    struct extruder_stepper *es = container_of(sk, struct extruder_stepper, sk);
    double hst = es->half_smooth_time;
    int i;
    struct coord e_pos;
    struct coord move_dir = get_move_dir(m, move_time, es);
    for (i = 0; i < 3; ++i) {
        int axis = 'x' + i;
        struct shaper_pulses* sp = &es->sp[i];
        if (!hst) {
            e_pos.axis[i] = shaper_calc_position(
                    m, axis, move_dir.axis[i], move_time, sp);
        } else {
            double area = shaper_pa_range_integrate(
                    m, axis, move_dir.axis[i], move_time,
                    es->pressure_advance, hst, sp);
            e_pos.axis[i] = area * es->inv_half_smooth_time2;
        }
    }
    return e_pos.x + e_pos.y + e_pos.z;
}

static void
extruder_note_generation_time(struct extruder_stepper *es)
{
    double pre_active = 0., post_active = 0.;
    int i;
    for (i = 0; i < 2; ++i) {
        struct shaper_pulses* sp = &es->sp[i];
        if (!es->sp[i].num_pulses) continue;
        pre_active = sp->pulses[sp->num_pulses-1].t > pre_active
            ? sp->pulses[sp->num_pulses-1].t : pre_active;
        post_active = -sp->pulses[0].t > post_active
            ? -sp->pulses[0].t : post_active;
    }
    if (es->half_smooth_time) {
        pre_active += es->half_smooth_time;
        post_active += es->half_smooth_time;
    }
    es->sk.gen_steps_pre_active = pre_active;
    es->sk.gen_steps_post_active = post_active;
}

void __visible
extruder_set_pressure_advance(struct stepper_kinematics *sk
                              , double pressure_advance, double smooth_time)
{
    struct extruder_stepper *es = container_of(sk, struct extruder_stepper, sk);
    double hst = smooth_time * .5;
    es->half_smooth_time = hst;
    extruder_note_generation_time(es);
    if (! hst)
        return;
    es->inv_half_smooth_time2 = 1. / (hst * hst);
    es->pressure_advance = pressure_advance;
}

int __visible
extruder_set_shaper_params(struct stepper_kinematics *sk, char axis
                           , int n, double a[], double t[])
{
    if (axis != 'x' && axis != 'y')
        return -1;
    struct extruder_stepper *es = container_of(sk, struct extruder_stepper, sk);
    struct shaper_pulses *sp = &es->sp[axis-'x'];
    int status = init_shaper(n, a, t, sp);
    extruder_note_generation_time(es);
    return status;
}

struct stepper_kinematics * __visible
extruder_stepper_alloc(void)
{
    struct extruder_stepper *es = malloc(sizeof(*es));
    memset(es, 0, sizeof(*es));
    es->sk.calc_position_cb = extruder_calc_position;
    es->sk.active_flags = AF_X | AF_Y | AF_Z;
    return &es->sk;
}
