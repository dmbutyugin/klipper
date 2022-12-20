// Extruder stepper pulse time generation
//
// Copyright (C) 2018-2019  Kevin O'Connor <kevin@koconnor.net>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <stddef.h> // offsetof
#include <stdlib.h> // malloc
#include <string.h> // memset
#include "compiler.h" // __visible
#include "itersolve.h" // struct stepper_kinematics
#include "kin_shaper.h" // struct shaper_pulses
#include "pyhelper.h" // errorf
#include "trapq.h" // move_get_distance

// Without pressure advance, the extruder stepper position is:
//     extruder_position(t) = nominal_position(t)
// When pressure advance is enabled, additional filament is pushed
// into the extruder during acceleration (and retracted during
// deceleration). The formula is:
//     pa_position(t) = (nominal_position(t)
//                       + advance * nominal_velocity(t)
//                       [ + advance2 * nominal_velocity(t)^2])
// Where the nominal position and velocitiy are actually "smoothed" using
// a weighted average prior to using PA formula:
//     smooth_position(t) = (
//         definitive_integral(nominal_position(x) *
//                                 (smooth_time/2 - abs(t-x)) * dx,
//                             from=t-smooth_time/2, to=t+smooth_time/2)
//         / ((smooth_time/2)**2))
// and similarly for nominal_velocity(x).

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

// Calculate the definitive integrals of extruder for a given move
static void
pa_move_integrate(struct move *m, int axis, double linear_velocity
                  , double linear_offset, double linear_advance
                  , double base, double start, double end, double time_offset
                  , double *pos_integral, double *pa_vel_integral)
{
    if (start < 0.)
        start = 0.;
    if (end > m->move_t)
        end = m->move_t;
    double axis_r = m->axes_r.axis[axis - 'x'];
    double start_v = m->start_v * axis_r;
    double ha = m->half_accel * axis_r;
    // Calculate base position and velocity for pressure advance
    int can_pressure_advance = m->axes_r.x > 0. || m->axes_r.y > 0.;
    // Calculate definitive integrals
    double iext = extruder_integrate(base, start_v, ha, start, end);
    double wgt_ext = extruder_integrate_time(base, start_v, ha, start, end);
    *pos_integral = wgt_ext - time_offset * iext;
    if (!can_pressure_advance) {
        *pa_vel_integral = 0.;
    } else {
        double ivel = extruder_integrate(start_v, 2. * ha, 0., start, end);
        double wgt_vel = extruder_integrate_time(start_v, 2. * ha, 0., start, end);
        *pa_vel_integral = wgt_vel - time_offset * ivel;
    }
}

// Calculate the definitive integrals of the extruder over a range of moves
static void
pa_range_integrate(struct move *m, int axis, double move_time
                   , double linear_velocity, double linear_offset
                   , double linear_advance, double hst
                   , double *pos_integral, double *pa_vel_integral)
{
    while (unlikely(move_time < 0.)) {
        m = list_prev_entry(m, node);
        move_time += m->move_t;
    }
    while (unlikely(move_time > m->move_t)) {
        move_time -= m->move_t;
        m = list_next_entry(m, node);
    }
    // Calculate integrals for the current move
    double start = move_time - hst, end = move_time + hst;
    double start_base = m->start_pos.axis[axis - 'x'];
    *pos_integral = *pa_vel_integral = 0.;
    double move_pos_int, move_pa_vel_int;
    pa_move_integrate(
            m, axis, linear_velocity, linear_offset, linear_advance, 0.,
            start, move_time, start, &move_pos_int, &move_pa_vel_int);
    *pos_integral += move_pos_int;
    *pa_vel_integral += move_pa_vel_int;
    pa_move_integrate(
            m, axis, linear_velocity, linear_offset, linear_advance, 0.,
            move_time, end, end, &move_pos_int, &move_pa_vel_int);
    *pos_integral -= move_pos_int;
    *pa_vel_integral -= move_pa_vel_int;
    // Integrate over previous moves
    struct move *prev = m;
    while (unlikely(start < 0.)) {
        prev = list_prev_entry(prev, node);
        start += prev->move_t;
        double base = prev->start_pos.axis[axis - 'x'] - start_base;
        pa_move_integrate(prev, axis, linear_velocity, linear_offset,
                          linear_advance, base, start, prev->move_t,
                          start, &move_pos_int, &move_pa_vel_int);
        *pos_integral += move_pos_int;
        *pa_vel_integral += move_pa_vel_int;
    }
    // Integrate over future moves
    while (unlikely(end > m->move_t)) {
        end -= m->move_t;
        m = list_next_entry(m, node);
        double base = m->start_pos.axis[axis - 'x'] - start_base;
        pa_move_integrate(m, axis, linear_velocity, linear_offset,
                          linear_advance, base, 0., end, end,
                          &move_pos_int, &move_pa_vel_int);
        *pos_integral -= move_pos_int;
        *pa_vel_integral -= move_pa_vel_int;
    }
    *pos_integral += start_base * hst * hst;
}

static void
shaper_pa_range_integrate(struct move *m, int axis, double move_time
                          , double linear_velocity, double linear_offset
                          , double linear_advance, double hst
                          , struct shaper_pulses *sp
                          , double *pos_integral, double *pa_vel_integral)
{
    int num_pulses = sp->num_pulses, i;
    *pos_integral = *pa_vel_integral = 0.;
    for (i = 0; i < num_pulses; ++i) {
        double t = sp->pulses[i].t, a = sp->pulses[i].a;
        double pulse_pos_int, pulse_pa_vel_int;
        pa_range_integrate(
                m, axis, move_time + t, linear_velocity, linear_offset,
                linear_advance, hst, &pulse_pos_int, &pulse_pa_vel_int);
        *pos_integral += a * pulse_pos_int;
        *pa_vel_integral += a * pulse_pa_vel_int;
    }
}

struct extruder_stepper {
    struct stepper_kinematics sk;
    struct shaper_pulses sp[3];
    double linear_velocity, linear_offset, linear_advance;
    double half_smooth_time, inv_half_smooth_time2;
};

static void
extruder_calc_position_and_pa_velocity(struct extruder_stepper *es
                                       , struct move *m, double move_time
                                       , double *pos, double *pa_vel)
{
    double hst = es->half_smooth_time;
    int i;
    struct coord e_pos, e_vel;
    double move_dist = move_get_distance(m, move_time);
    for (i = 0; i < 3; ++i) {
        int axis = 'x' + i;
        struct shaper_pulses* sp = &es->sp[i];
        int num_pulses = sp->num_pulses;
        if (!hst) {
            e_pos.axis[i] = num_pulses
                ? shaper_calc_position(m, axis, move_time, sp)
                : m->start_pos.axis[i] + m->axes_r.axis[i] * move_dist;
            e_vel.axis[i] = 0.;
        } else {
            if (!num_pulses) {
                pa_range_integrate(
                        m, axis, move_time, es->linear_velocity,
                        es->linear_offset, es->linear_advance, hst,
                        e_pos.axis + i, e_vel.axis + i);
            } else {
                shaper_pa_range_integrate(
                        m, axis, move_time, es->linear_velocity,
                        es->linear_offset, es->linear_advance, hst, sp,
                        e_pos.axis + i, e_vel.axis + i);
            }
            e_pos.axis[i] *= es->inv_half_smooth_time2;
            e_vel.axis[i] *= es->inv_half_smooth_time2;
        }
    }
    *pos = e_pos.x + e_pos.y + e_pos.z;
    *pa_vel = e_vel.x + e_vel.y + e_vel.z;
}

static double
extruder_calc_position(struct stepper_kinematics *sk, struct move *m
                       , double move_time)
{
    struct extruder_stepper *es = container_of(sk, struct extruder_stepper, sk);
    double pos, pa_vel;
    extruder_calc_position_and_pa_velocity(es, m, move_time, &pos, &pa_vel);
    if (!es->half_smooth_time)
        return pos;
    pos += es->linear_advance * pa_vel;
    if (pa_vel >= es->linear_velocity) {
        pos += es->linear_offset;
    } else {
        double rel_vel = pa_vel / es->linear_velocity;
        pos += es->linear_offset * rel_vel * (2. - rel_vel);
    }
    return pos;
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
                              , double linear_velocity, double linear_offset
                              , double linear_advance, double smooth_time)
{
    struct extruder_stepper *es = container_of(sk, struct extruder_stepper, sk);
    double hst = smooth_time * .5;
    es->half_smooth_time = hst;
    extruder_note_generation_time(es);
    if (! hst)
        return;
    es->inv_half_smooth_time2 = 1. / (hst * hst);
    es->linear_velocity = linear_velocity;
    es->linear_offset = linear_offset;
    es->linear_advance = linear_advance;
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
