// Bezier curve acceleration
//
// Copyright (C) 2018-2020  Kevin O'Connor <kevin@koconnor.net>
// Copyright (C) 2019-2022  Dmitry Butyugin <dmbutyugin@google.com>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <math.h> // fabs
#include <string.h> // memset

#include "compiler.h" // likely
#include "scurve.h" // scurve

// Find the distance travel on an S-Curve
double
scurve_eval(const struct scurve *s, double time)
{
    double v = s->c6;
    v = s->c5 + v * time;
    v = s->c4 + v * time;
    v = s->c3 + v * time;
    v = s->c2 + v * time;
    v = s->c1 + v * time;
    return v * time;
}

double
scurve_velocity(const struct scurve *s, double time)
{
    double v = 6. * s->c6;
    v = 5. * s->c5 + v * time;
    v = 4. * s->c4 + v * time;
    v = 3. * s->c3 + v * time;
    v = 2. * s->c2 + v * time;
    v = s->c1 + v * time;
    return v;
}

static void
scurve_fill_bezier2(struct scurve *s, double start_accel_v
                  , double effective_accel, double accel_offset_t)
{
    s->c2 = .5 * effective_accel;
    s->c1 = start_accel_v + effective_accel * accel_offset_t;
}

// Determine the coefficients for a 4th order bezier position function
static void
scurve_fill_bezier4(struct scurve *s, double start_accel_v
        , double effective_accel, double total_accel_t, double accel_offset_t)
{
    if (!total_accel_t)
        return;
    double inv_accel_t = 1. / total_accel_t;
    double accel_div_accel_t = effective_accel * inv_accel_t;
    double accel_div_accel_t2 = accel_div_accel_t * inv_accel_t;
    s->c4 = -.5 * accel_div_accel_t2;
    s->c3 = accel_div_accel_t;
    s->c1 = start_accel_v;

    scurve_offset(s, accel_offset_t);
}

// Determine the coefficients for a 6th order bezier position function
static void
scurve_fill_bezier6(struct scurve *s, double start_accel_v
        , double effective_accel, double total_accel_t, double accel_offset_t)
{
    if (!total_accel_t)
        return;
    double inv_accel_t = 1. / total_accel_t;
    double accel_div_accel_t2 = effective_accel * inv_accel_t * inv_accel_t;
    double accel_div_accel_t3 = accel_div_accel_t2 * inv_accel_t;
    double accel_div_accel_t4 = accel_div_accel_t3 * inv_accel_t;
    s->c6 = accel_div_accel_t4;
    s->c5 = -3. * accel_div_accel_t3;
    s->c4 = 2.5 * accel_div_accel_t2;
    s->c1 = start_accel_v;

    scurve_offset(s, accel_offset_t);
}

double scurve_get_time(const struct scurve *s, double distance)
{
    double low = 0;
    double high = s->total_accel_t;
    if (scurve_eval(s, high) <= distance) return high;
    if (scurve_eval(s, low) > distance) return low;
    while (likely(high - low > .000000001)) {
        double guess_time = (high + low) * .5;
        if (scurve_eval(s, guess_time) > distance) {
            high = guess_time;
        } else {
            low = guess_time;
        }
    }
    return (high + low) * .5;
}

void
scurve_fill(struct scurve *s, int accel_order
            , double accel_t, double accel_offset_t, double total_accel_t
            , double start_accel_v, double effective_accel)
{
    memset(s, 0, sizeof(*s));
    s->total_accel_t = total_accel_t;
    if (accel_order == 4) {
        scurve_fill_bezier4(s, start_accel_v, effective_accel
                , total_accel_t, accel_offset_t);
    } else if (accel_order == 6) {
        scurve_fill_bezier6(s, start_accel_v, effective_accel
                , total_accel_t, accel_offset_t);
    } else {
        scurve_fill_bezier2(s, start_accel_v, effective_accel, accel_offset_t);
    }
}

void
scurve_offset(struct scurve *s, double offset_t)
{
    s->c1 += ((((6. * s->c6 * offset_t + 5. * s->c5) * offset_t
                    + 4. * s->c4) * offset_t
                + 3. * s->c3) * offset_t + 2. * s->c2) * offset_t;
    s->c2 += (((15. * s->c6 * offset_t + 10. * s->c5) * offset_t
                + 6. * s->c4) * offset_t + 3. * s->c3) * offset_t;
    s->c3 += ((20. * s->c6 * offset_t
                + 10. * s->c5) * offset_t + 4. * s->c4) * offset_t;
    s->c4 += (15. * s->c6 * offset_t + 5. * s->c5) * offset_t;
    s->c5 += 6. * s->c6 * offset_t;
}

void
scurve_copy_scaled(const struct scurve *src, double ratio, struct scurve *dst)
{
    dst->total_accel_t = src->total_accel_t;

    dst->c6 = src->c6 * ratio;
    dst->c5 = src->c5 * ratio;
    dst->c4 = src->c4 * ratio;
    dst->c3 = src->c3 * ratio;
    dst->c2 = src->c2 * ratio;
    dst->c1 = src->c1 * ratio;
}

double
scurve_add_deriv(const struct scurve *src, double ratio, struct scurve *dst)
{
    dst->c5 += 6. * src->c6 * ratio;
    dst->c4 += 5. * src->c5 * ratio;
    dst->c3 += 4. * src->c4 * ratio;
    dst->c2 += 3. * src->c3 * ratio;
    dst->c1 += 2. * src->c2 * ratio;
    return src->c1 * ratio;
}

double
scurve_add_2nd_deriv(const struct scurve *src, double ratio, struct scurve *dst)
{
    dst->c4 += 30. * src->c6 * ratio;
    dst->c3 += 20. * src->c5 * ratio;
    dst->c2 += 12. * src->c4 * ratio;
    dst->c1 += 6. * src->c3 * ratio;
    return 2. * src->c2 * ratio;
}

static const double scurve_antideriv_coeffs[] = {
    1./1., 1./2., 1./3., 1./4., 1./5., 1./6., 1./7., 1./8., 1./9., 1./10.,
    1./11., 1./12., 1./13., 1./14., 1./15.,
};

double
scurve_tn_antiderivative(const struct scurve *s, int n, double time)
{
    const double *coeffs = scurve_antideriv_coeffs + n;
    double v = s->c6 * coeffs[6];
    v = s->c5 * coeffs[5] + v * time;
    v = s->c4 * coeffs[4] + v * time;
    v = s->c3 * coeffs[3] + v * time;
    v = s->c2 * coeffs[2] + v * time;
    v = s->c1 * coeffs[1] + v * time;
    for (; likely(n >= 0); --n)
        v *= time;
    return v * time;
}

double
scurve_tn_integrate(const struct scurve *s, int n, double start, double end)
{
    const double *coeffs = scurve_antideriv_coeffs + n;
    double delta_t = end - start;
    double sum = end + start, e_n = end;
    for (; unlikely(n > 0); --n) {
        e_n *= end;
        sum = sum * start + e_n;
    }
    double avg_val = coeffs[1] * s->c1 * sum;
    e_n *= end;
    sum = sum * start + e_n;
    avg_val += coeffs[2] * s->c2 * sum;
    e_n *= end;
    sum = sum * start + e_n;
    avg_val += coeffs[3] * s->c3 * sum;
    e_n *= end;
    sum = sum * start + e_n;
    avg_val += coeffs[4] * s->c4 * sum;
    e_n *= end;
    sum = sum * start + e_n;
    avg_val += coeffs[5] * s->c5 * sum;
    e_n *= end;
    sum = sum * start + e_n;
    avg_val += coeffs[6] * s->c6 * sum;
    return delta_t * avg_val;
}
