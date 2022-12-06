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
pa_move_integrate(struct move *m, double linear_velocity, double linear_offset
                  , double linear_advance, double base, double start, double end
                  , double time_offset)
{
    if (start < 0.)
        start = 0.;
    if (end > m->move_t)
        end = m->move_t;
    double start_v = m->start_v;
    double ha = m->half_accel;
    // Calculate base position and velocity with pressure advance
    int can_pressure_advance = m->axes_r.y != 0.;
    if (can_pressure_advance) {
        double switch_t =
            !ha || !linear_velocity ?  0.
                                    : .5 * (linear_velocity - start_v) / ha;
        if (switch_t > start && switch_t < end) {
            double res = pa_move_integrate(
                    m, linear_velocity, linear_offset, linear_advance,
                    base, start, switch_t, time_offset);
            res += pa_move_integrate(
                    m, linear_velocity, linear_offset, linear_advance,
                    base, switch_t, end, time_offset);
            return res;
        } else if (start_v + ha * (start + end) < linear_velocity) {
            double recip_lv = 1. / linear_velocity;
            base += (linear_advance
                    + linear_offset * recip_lv * (2. - recip_lv * start_v)) * start_v;
            start_v += (linear_advance
                        + 2. * linear_offset * recip_lv * (1. - recip_lv * start_v)) * 2. * ha;
            ha -= 4. * linear_offset * recip_lv * recip_lv * ha * ha;
        } else {
            base += linear_offset + linear_advance * start_v;
            start_v += 2. * linear_advance * ha;
        }
    }
    // Calculate definitive integral
    double iext = extruder_integrate(base, start_v, ha, start, end);
    double wgt_ext = extruder_integrate_time(base, start_v, ha, start, end);
    return wgt_ext - time_offset * iext;
}

// Calculate the definitive integral of the extruder over a range of moves
static double
pa_range_integrate(struct move *m, double move_time, double linear_velocity
                   , double linear_offset, double linear_advance, double hst)
{
    // Calculate integral for the current move
    double res = 0., start = move_time - hst, end = move_time + hst;
    double start_base = m->start_pos.x;
    res += pa_move_integrate(m, linear_velocity, linear_offset, linear_advance,
                             0., start, move_time, start);
    res -= pa_move_integrate(m, linear_velocity, linear_offset, linear_advance,
                             0., move_time, end, end);
    // Integrate over previous moves
    struct move *prev = m;
    while (unlikely(start < 0.)) {
        prev = list_prev_entry(prev, node);
        start += prev->move_t;
        double base = prev->start_pos.x - start_base;
        res += pa_move_integrate(
                prev, linear_velocity, linear_offset, linear_advance,
                base, start, prev->move_t, start);
    }
    // Integrate over future moves
    while (unlikely(end > m->move_t)) {
        end -= m->move_t;
        m = list_next_entry(m, node);
        double base = m->start_pos.x - start_base;
        res -= pa_move_integrate(
                m, linear_velocity, linear_offset, linear_advance,
                base, 0., end, end);
    }
    return res;
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
    double area = pa_range_integrate(
            m, move_time, es->linear_velocity, es->linear_offset,
            es->linear_advance, hst);
    return m->start_pos.x + area * es->inv_half_smooth_time2;
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
    es->linear_velocity = linear_velocity;
    es->linear_offset = linear_offset;
    es->linear_advance = linear_advance;
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
