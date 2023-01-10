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
#include "pyhelper.h" // errorf
#include "trapq.h" // move_get_distance

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
static void
pa_move_integrate(struct move *m, double base, double start, double end
                  , double time_offset, double *pos_integral
                  , double *pa_velocity_integral)
{
    if (start < 0.)
        start = 0.;
    if (end > m->move_t)
        end = m->move_t;
    // Calculate base position and velocity with pressure advance
    int can_pressure_advance = m->axes_r.y != 0.;
    // Calculate definitive integral
    double start_v = m->start_v;
    double ha = m->half_accel;
    double iext = extruder_integrate(base, start_v, ha, start, end);
    double wgt_ext = extruder_integrate_time(base, start_v, ha, start, end);
    *pos_integral = wgt_ext - time_offset * iext;
    if (!can_pressure_advance) {
        *pa_velocity_integral = 0.;
    } else {
        double ivel = extruder_integrate(start_v, 2. * ha, 0., start, end);
        double wgt_vel = extruder_integrate(0., start_v, 2. * ha, start, end);
        *pa_velocity_integral = wgt_vel - time_offset * ivel;
    }
}

// Calculate the definitive integral of the extruder over a range of moves
static void
pa_range_integrate(struct move *m, double move_time, double hst
                   , double *pos_integral, double *pa_velocity_integral)
{
    // Calculate integral for the current move
    *pos_integral = *pa_velocity_integral = 0.;
    double m_pos_int, m_pa_vel_int;
    double start = move_time - hst, end = move_time + hst;
    double start_base = m->start_pos.x;
    pa_move_integrate(m, 0., start, move_time, start,
                      &m_pos_int, &m_pa_vel_int);
    *pos_integral += m_pos_int;
    *pa_velocity_integral += m_pa_vel_int;
    pa_move_integrate(m, 0., move_time, end, end, &m_pos_int, &m_pa_vel_int);
    *pos_integral -= m_pos_int;
    *pa_velocity_integral -= m_pa_vel_int;
    // Integrate over previous moves
    struct move *prev = m;
    while (unlikely(start < 0.)) {
        prev = list_prev_entry(prev, node);
        start += prev->move_t;
        double base = prev->start_pos.x - start_base;
        pa_move_integrate(prev, base, start, prev->move_t, start,
                          &m_pos_int, &m_pa_vel_int);
        *pos_integral += m_pos_int;
        *pa_velocity_integral += m_pa_vel_int;
    }
    // Integrate over future moves
    while (unlikely(end > m->move_t)) {
        end -= m->move_t;
        m = list_next_entry(m, node);
        double base = m->start_pos.x - start_base;
        pa_move_integrate(m, base, 0., end, end, &m_pos_int, &m_pa_vel_int);
        *pos_integral -= m_pos_int;
        *pa_velocity_integral -= m_pa_vel_int;
    }
}

struct extruder_stepper {
    struct stepper_kinematics sk;
    double linear_velocity, linear_offset, linear_advance;
    double half_smooth_time, inv_half_smooth_time2;
};

static double
extruder_calc_position(struct stepper_kinematics *sk, struct move *m
                       , double move_time)
{
    struct extruder_stepper *es = container_of(sk, struct extruder_stepper, sk);
    double hst = es->half_smooth_time;
    if (!hst)
        // Pressure advance not enabled
        return m->start_pos.x + move_get_distance(m, move_time);
    // Apply pressure advance and average over smooth_time
    double pa_pos, pa_velocity;
    pa_range_integrate(m, move_time, hst, &pa_pos, &pa_velocity);
    pa_pos *= es->inv_half_smooth_time2;
    pa_velocity *= es->inv_half_smooth_time2;
    pa_pos += es->linear_advance * pa_velocity;
    if (pa_velocity < 0.) pa_velocity = 0.;
    if (pa_velocity < es->linear_velocity) {
        double rel_vel = pa_velocity / es->linear_velocity;
        pa_pos += es->linear_offset * (2. * sqrt(rel_vel) - rel_vel);
    } else {
        pa_pos += es->linear_offset;
    }
    return m->start_pos.x + pa_pos;
}

void __visible
extruder_set_pressure_advance(struct stepper_kinematics *sk
                              , double linear_velocity, double linear_offset
                              , double linear_advance, double smooth_time)
{
    struct extruder_stepper *es = container_of(sk, struct extruder_stepper, sk);
    double hst = smooth_time * .5;
    es->half_smooth_time = hst;
    es->sk.gen_steps_pre_active = es->sk.gen_steps_post_active = hst;
    if (! hst)
        return;
    es->inv_half_smooth_time2 = 1. / (hst * hst);
    es->linear_advance = linear_advance;
    if (linear_velocity > 0.) {
        es->linear_velocity = linear_velocity;
        es->linear_offset = linear_offset;
    } else {
        es->linear_velocity = 0.;
        es->linear_offset = 0.;
    }
}

struct stepper_kinematics * __visible
extruder_stepper_alloc(void)
{
    struct extruder_stepper *es = malloc(sizeof(*es));
    memset(es, 0, sizeof(*es));
    es->sk.calc_position_cb = extruder_calc_position;
    es->sk.active_flags = AF_X;
    return &es->sk;
}
