#ifndef __KIN_SHAPER_H
#define __KIN_SHAPER_H

#include "trapq.h" // struct move

struct shaper_pulses {
    int num_pulses;
    struct {
        double t, a;
    } pulses[5];
};

int init_shaper(int n, double a[], double t[], struct shaper_pulses *sp);

#endif  // kin_shaper.h
