#ifndef SCURVE_H
#define SCURVE_H

struct scurve {
    double c1, c2, c3, c4, c5, c6;
    double total_accel_t;
};

double scurve_eval(const struct scurve *s, double time);
double scurve_velocity(const struct scurve *s, double time);
double scurve_tn_antiderivative(const struct scurve *s, int n, double time);
double scurve_tn_integrate(const struct scurve *s, int n
                           , double start, double end);

void scurve_offset(struct scurve *s, double offset_t);
void scurve_copy_scaled(const struct scurve *src, double ratio
                        , struct scurve *dst);
double scurve_add_deriv(const struct scurve *src, double ratio
                        , struct scurve *dst);
double scurve_add_2nd_deriv(const struct scurve *src, double ratio
                            , struct scurve *dst);
void scurve_fill(struct scurve *s, int accel_order
                 , double accel_t, double accel_offset_t, double total_accel_t
                 , double start_accel_v, double effective_accel);
double scurve_get_time(const struct scurve *s, double distance);

#endif // scurve.h
