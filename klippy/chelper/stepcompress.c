// Stepper pulse schedule compression
//
// Copyright (C) 2016-2021  Kevin O'Connor <kevin@koconnor.net>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

// The goal of this code is to take a series of scheduled stepper
// pulse times and compress them into a handful of commands that can
// be efficiently transmitted and executed on a microcontroller (mcu).
// The mcu accepts step pulse commands that take interval, count, and
// add parameters such that 'count' pulses occur, with each step event
// calculating the next step event time using:
//  next_wake_time = last_wake_time + interval; interval += add
// This code is written in C (instead of python) for processing
// efficiency - the repetitive integer math is vastly faster in C.

#include <math.h> // sqrt
#include <stddef.h> // offsetof
#include <stdint.h> // uint32_t
#include <stdio.h> // fprintf
#include <stdlib.h> // malloc
#include <string.h> // memset
#include "compiler.h" // DIV_ROUND_UP
#include "pyhelper.h" // errorf
#include "serialqueue.h" // struct queue_message
#include "stepcompress.h" // stepcompress_alloc

#define CHECK_LINES 1
#define QUEUE_START_SIZE 1024

// Maximum error between compressed and actual step,
// in the form of (step_i - step_{i-1}) >> MAX_ERR_2P
#define MAX_ERR_2P 6
// Limits on step_move values
#define MAX_COUNT 256
// Limits below are optimized for "message blocks" encoding
// and underlying storage types and MCU-side implementation
#define MAX_INTRVL 0x3FFFFFF
#define MAX_ADD 0x7FFF
#define MAX_ADD2 0xFFF
#define MAX_SHIFT 16
#define MIN_SHIFT -8

struct stepcompress {
    // Buffer management
    uint32_t *queue, *queue_end, *queue_pos, *queue_next;
    // Internal tracking
    uint32_t max_error;
    double mcu_time_offset, mcu_freq, last_step_print_time;
    // Message generation
    uint64_t last_step_clock;
    struct list_head msg_queue;
    uint32_t oid;
    int32_t queue_step_msgtag, set_next_step_dir_msgtag;
    int sdir, invert_sdir;
    // Step+dir+step filter
    uint64_t next_step_clock;
    int next_step_dir;
    // History tracking
    int64_t last_position;
    struct list_head history_list;
};

struct step_move {
    uint32_t interval;
    uint16_t count;
    int16_t add;
    int16_t add2;
    int8_t shift;
    uint64_t first_step, last_step;
};

#define HISTORY_EXPIRE (30.0)

struct history_steps {
    struct list_node node;
    uint64_t first_clock, last_clock;
    int64_t start_position;
    int step_count, interval, add, add2, shift;
};


/****************************************************************
 * Step compression
 ****************************************************************/

// A bias of the first step in the generated step sequence step_move towards
// the first requires step in the least squares method. The zero value means
// no bias: the first step error is minimized as all other step errors, and
// values larger than zero make the method reduce the first step error more
// than other step errors. This helps pushing back the maximum errors towards
// the end of the generated step sequence.
#define FIRST_STEP_BIAS 1.0

struct matrix_3x3 {
    double a00, a10, a11, a20, a21, a22;
};

struct rhs_3 {
    double b0, b1, b2;
};

struct matrix_3x3 least_squares_ldl[MAX_COUNT] = {0};

static void
fill_least_squares_matrix_3x3(uint16_t count, struct matrix_3x3 *m)
{
    int32_t i;
    memset(m, 0, sizeof(*m));
    for (i = 0; i < count; ++i) {
        int32_t c0 = i+1;
        int32_t c1 = c0 * i / 2;
        int32_t c2 = c1 * (i-1) / 3;

        m->a00 += (double)c0 * c0;
        m->a10 += (double)c1 * c0;
        m->a11 += (double)c1 * c1;
        m->a20 += (double)c2 * c0;
        m->a21 += (double)c2 * c1;
        m->a22 += (double)c2 * c2;
    }
    m->a00 += FIRST_STEP_BIAS;
    if (i < 2) m->a11 = 1.;
    if (i < 3) m->a22 = 1.;
}

static void
compute_ldl_3x3(struct matrix_3x3 *m)
{
    double d0 = m->a00;
    m->a00 = 1. / d0;
    m->a10 *= m->a00;
    m->a20 *= m->a00;

    double d1 = m->a11 - d0 * m->a10 * m->a10;
    m->a11 = 1. / d1;
    m->a21 -= m->a20 * m->a10 * d0;
    m->a21 *= m->a11;

    double d2 = m->a22 - d0 * m->a20 * m->a20 - d1 * m->a21 * m->a21;
    m->a22 = 1. / d2;
}

static void
compute_rhs_3(struct stepcompress *sc, uint16_t count, struct rhs_3 *f)
{
    memset(f, 0, sizeof(*f));
    uint32_t lsc = sc->last_step_clock;
    f->b0 += FIRST_STEP_BIAS * (*sc->queue_pos - lsc);
    for (uint16_t i = 0; i < count; ++i) {
        double d = sc->queue_pos[i] - lsc;
        int32_t c = i+1;
        f->b0 += d * c;
        c = c * i / 2;
        f->b1 += d * c;
        c = c * (i-1) / 3;
        f->b2 += d * c;
    }
}

static void
solve_3x3(struct matrix_3x3 *m, struct rhs_3 *f)
{
    f->b1 -= f->b0 * m->a10;
    f->b2 -= f->b0 * m->a20 + f->b1 * m->a21;

    f->b0 *= m->a00;
    f->b1 *= m->a11;
    f->b2 *= m->a22;

    f->b1 -= f->b2 * m->a21;
    f->b0 -= f->b1 * m->a10 + f->b2 * m->a20;
}

static struct matrix_3x3*
get_least_squares_ldl_3x3(uint16_t count)
{
    if (count > MAX_COUNT) return NULL;
    struct matrix_3x3 *m = &least_squares_ldl[count-1];
    if (!m->a00) {
        fill_least_squares_matrix_3x3(count, m);
        compute_ldl_3x3(m);
    }
    return m;
}

static struct step_move
step_move_encode(uint16_t count, struct rhs_3* f)
{
    struct step_move res;
    memset(&res, 0, sizeof(res));
    double interval = f->b0, add = f->b1, add2 = f->b2;
    if (interval < 0.)
        return res;
    if (count <= 1) {
        res.count = count;
        res.interval = round(interval);
        return res;
    }
    double end_add = add + add2 * count;
    double max_int_inc = count * (fabs(add) > fabs(end_add) ?
                                  fabs(add) : fabs(end_add));
    double max_end_int = interval + max_int_inc;
    if (fabs(add) > MAX_ADD || fabs(end_add) > MAX_ADD ||
            fabs(add2) > MAX_ADD2 || max_end_int > MAX_INTRVL) {
        while (res.shift >= MIN_SHIFT && (
                    fabs(add) > MAX_ADD || fabs(end_add) > MAX_ADD ||
                    fabs(add2) > MAX_ADD2 || max_end_int > MAX_INTRVL)) {
            interval *= 0.5;
            add *= 0.5;
            add2 *= 0.5;
            end_add *= 0.5;
            max_end_int *= 0.5;
            --res.shift;
        }
        if (res.shift < MIN_SHIFT)
            // Cannot encode the current rhs_3, the values are too large
            return res;
    } else if (max_int_inc >= 0.5 ||
               count * fabs(interval-round(interval)) >= 0.5) {
        while (res.shift < MAX_SHIFT &&
                fabs(add * 2.) <= MAX_ADD && fabs(end_add * 2.) <= MAX_ADD &&
                fabs(add2 * 2.) <= MAX_ADD2 && max_end_int * 2. <= MAX_INTRVL) {
            interval *= 2.;
            add *= 2.;
            add2 *= 2.;
            end_add *= 2.;
            max_end_int *= 2.;
            ++res.shift;
        }
    }
    res.count = count;
    res.interval = round(interval);
    res.add = round(add);
    res.add2 = round(add2);
    return res;
}

struct points {
    int32_t minp, maxp;
};

// Given a requested step time, return the minimum and maximum
// acceptable times
static inline struct points
minmax_point(struct stepcompress *sc, uint32_t *pos)
{
    uint32_t lsc = sc->last_step_clock, point = *pos - lsc;
    uint32_t prevpoint = pos > sc->queue_pos ? *(pos-1) - lsc : 0;
    uint32_t nextpoint = pos + 1 < sc->queue_next ? *(pos+1) - lsc : point;
    uint32_t max_bck_error = (point - prevpoint) >> MAX_ERR_2P;
    uint32_t max_frw_error = (nextpoint - point) >> MAX_ERR_2P;
    if (max_bck_error > sc->max_error)
        max_bck_error = sc->max_error;
    if (max_frw_error > sc->max_error)
        max_frw_error = sc->max_error;
    return (struct points){ point - max_bck_error, point + max_frw_error };
}

struct stepper_moves {
    uint32_t interval;
    uint16_t int_low;
    int32_t int_low_acc;
    int32_t add;
    int32_t add2;
    uint32_t count;
};

static inline void
add_interval(uint32_t* time, struct stepper_moves *s)
{
    uint32_t next_time = *time + s->interval;
    if (likely(s->int_low)) {
        int32_t int_low_acc = s->int_low_acc + s->int_low;
        if (unlikely(int_low_acc >= 0 && ((int_low_acc & 0x00008000) ||
                                          ((int_low_acc >> 16) & 0xff)))) {
            ++next_time;
            int_low_acc -= 0x10000;
        }
        s->int_low_acc = int_low_acc;
    }
    *time = next_time;
}

static inline void
inc_interval(struct stepper_moves *s)
{
    uint32_t int_add = (uint32_t)(s->int_low + s->add);
    s->int_low = int_add & 0xffff;
    uint16_t int_add_high = int_add >> 16;
    if (!(int_add_high & 0x8000))
        s->interval += int_add_high;
    else
        s->interval -= (uint16_t)~int_add_high + 1;
    s->add += s->add2;
}

static void
fill_stepper_moves(struct step_move *m, struct stepper_moves *s)
{
    s->count = m->count;
    uint32_t interval = m->interval;
    int16_t add = m->add;
    int16_t add2 = m->add2;
    int8_t shift = m->shift;

    if (shift > 0) {
        s->interval = interval >> shift;
        s->int_low = (interval << (16 - shift)) & 0xFFFF;
    } else {
        s->interval = interval << -shift;
        s->int_low = 0;
    }
    // Left shift of the signed int is undefined behavior,
    // use addition instead.
    int32_t add_shifted = add, add2_shifted = add2;
    for (uint_fast8_t i = 16 - shift; i > 0; --i) {
        add_shifted += add_shifted;
        add2_shifted += add2_shifted;
    }
    s->add = add_shifted;
    s->add2 = add2_shifted;
    s->int_low_acc = 0;
}

static int
test_step_move(struct stepcompress *sc, struct step_move *m, int report_errors)
{
    if (!m->count || (!m->interval && !m->add && !m->add2 && m->count > 1)
        || m->interval >= 0x80000000
        || m->add < -MAX_ADD || m->add > MAX_ADD
        || m->add2 < -MAX_ADD2 || m->add2 > MAX_ADD2
        || m->shift > MAX_SHIFT || m->shift < MIN_SHIFT) {
        if (report_errors)
            errorf("stepcompress o=%d i=%d c=%d a=%d, a2=%d, s=%d:"
                   " Invalid sequence"
                   , sc->oid, m->interval, m->count, m->add, m->add2, m->shift);
        m->count = 0;
        return ERROR_RET;
    }
    struct stepper_moves s;
    fill_stepper_moves(m, &s);
    uint16_t i;
    uint32_t cur_step = 0;
    for (i = 0; i < m->count; ++i) {
        add_interval(&cur_step, &s);
        struct points point = minmax_point(sc, sc->queue_pos + i);
        if (cur_step < point.minp || cur_step > point.maxp) {
            if (report_errors)
                errorf("stepcompress o=%d i=%d c=%d a=%d, a2=%d, s=%d:"
                       " Point %u: %d not in %d:%d"
                       , sc->oid, m->interval, m->count, m->add, m->add2
                       , m->shift, i+1, cur_step, point.minp, point.maxp);
            // The least squares method does not minimize the maximum error
            // in the generated step sequence, but rather the total error.
            // However, we can still use the longest good generated prefix.
            m->count = i;
            return ERROR_RET;
        }
        inc_interval(&s);
        if (s.interval >= 0x80000000) {
            if (report_errors)
                errorf("stepcompress o=%d i=%d c=%d a=%d, a2=%d, s=%d:"
                       " Point %d: interval overflow %d"
                       , sc->oid, m->interval, m->count, m->add, m->add2
                       , m->shift, i+1, s.interval);
            m->count = i;
            return ERROR_RET;
        }
        m->last_step = cur_step;
        if (!m->first_step)
            m->first_step = cur_step;
    }
    return 0;
}

static struct step_move
test_step_count(struct stepcompress *sc, uint16_t count)
{
    struct matrix_3x3 *m = get_least_squares_ldl_3x3(count);
    struct rhs_3 f;
    compute_rhs_3(sc, count, &f);
    solve_3x3(m, &f);
    struct step_move res = step_move_encode(count, &f);
    test_step_move(sc, &res, /*report_errors=*/0);
    return res;
}

// Find a 'step_move' that covers a series of step times
static struct step_move
compress_bisect_count(struct stepcompress *sc)
{
    uint16_t left = 0, right = 8;
    struct step_move cur, best;
    best.count = 0;
    for (; right <= MAX_COUNT; right <<= 1) {
        cur = test_step_count(sc, right);
        if (cur.count > left) {
            best = cur;
            left = cur.count;
        }
        else break;
    }
    if (right > MAX_COUNT) right = MAX_COUNT + 1;
    while (right - left > 1) {
        uint16_t count = (left + right) / 2;
        cur = test_step_count(sc, count);
        if (cur.count > best.count) best = cur;
        if (cur.count < count) right = count;
        if (cur.count > left) left = count;
    }
    if (best.count <= 1) {
        uint32_t interval = *sc->queue_pos - (uint32_t)sc->last_step_clock;
        return (struct step_move){ interval, 1, 0, 0, 0, interval, interval };
    }
    return best;
}


/****************************************************************
 * Step compress checking
 ****************************************************************/

// Verify that a given 'step_move' matches the actual step times
static inline int
check_line(struct stepcompress *sc, struct step_move move)
{
    if (!CHECK_LINES)
        return 0;
    return test_step_move(sc, &move, /*report_errors=*/1);
}


/****************************************************************
 * Step compress interface
 ****************************************************************/

// Allocate a new 'stepcompress' object
struct stepcompress * __visible
stepcompress_alloc(uint32_t oid)
{
    struct stepcompress *sc = malloc(sizeof(*sc));
    memset(sc, 0, sizeof(*sc));
    list_init(&sc->msg_queue);
    list_init(&sc->history_list);
    sc->oid = oid;
    sc->sdir = -1;
    return sc;
}

// Fill message id information
void __visible
stepcompress_fill(struct stepcompress *sc, uint32_t max_error
                  , int32_t queue_step_msgtag, int32_t set_next_step_dir_msgtag)
{
    sc->max_error = max_error;
    sc->queue_step_msgtag = queue_step_msgtag;
    sc->set_next_step_dir_msgtag = set_next_step_dir_msgtag;
}

// Set the inverted stepper direction flag
void __visible
stepcompress_set_invert_sdir(struct stepcompress *sc, uint32_t invert_sdir)
{
    invert_sdir = !!invert_sdir;
    if (invert_sdir != sc->invert_sdir) {
        sc->invert_sdir = invert_sdir;
        if (sc->sdir >= 0)
            sc->sdir ^= 1;
    }
}

// Helper to free items from the history_list
static void
free_history(struct stepcompress *sc, uint64_t end_clock)
{
    while (!list_empty(&sc->history_list)) {
        struct history_steps *hs = list_last_entry(
            &sc->history_list, struct history_steps, node);
        if (hs->last_clock > end_clock)
            break;
        list_del(&hs->node);
        free(hs);
    }
}

// Free memory associated with a 'stepcompress' object
void __visible
stepcompress_free(struct stepcompress *sc)
{
    if (!sc)
        return;
    free(sc->queue);
    message_queue_free(&sc->msg_queue);
    free_history(sc, UINT64_MAX);
    free(sc);
}

uint32_t
stepcompress_get_oid(struct stepcompress *sc)
{
    return sc->oid;
}

int
stepcompress_get_step_dir(struct stepcompress *sc)
{
    return sc->next_step_dir;
}

// Determine the "print time" of the last_step_clock
static void
calc_last_step_print_time(struct stepcompress *sc)
{
    double lsc = sc->last_step_clock;
    sc->last_step_print_time = sc->mcu_time_offset + (lsc - .5) / sc->mcu_freq;

    if (lsc > sc->mcu_freq * HISTORY_EXPIRE)
        free_history(sc, lsc - sc->mcu_freq * HISTORY_EXPIRE);
}

// Set the conversion rate of 'print_time' to mcu clock
static void
stepcompress_set_time(struct stepcompress *sc
                      , double time_offset, double mcu_freq)
{
    sc->mcu_time_offset = time_offset;
    sc->mcu_freq = mcu_freq;
    calc_last_step_print_time(sc);
}

// Maximium clock delta between messages in the queue
#define CLOCK_DIFF_MAX (3<<28)

// Helper to create a queue_step command from a 'struct step_move'
static void
add_move(struct stepcompress *sc, uint64_t first_clock, struct step_move *move)
{
    uint64_t last_clock = sc->last_step_clock + move->last_step;

    // Create and queue a queue_step command
    uint32_t msg[7] = {
        sc->queue_step_msgtag, sc->oid, move->interval, move->count,
        move->add, move->add2, move->shift
    };
    struct queue_message *qm = message_alloc_and_encode(msg, 7);
    qm->min_clock = qm->req_clock = sc->last_step_clock;
    if (move->count == 1 && first_clock >= sc->last_step_clock + CLOCK_DIFF_MAX)
        qm->req_clock = first_clock;
    list_add_tail(&qm->node, &sc->msg_queue);
    sc->last_step_clock = last_clock;

    // Create and store move in history tracking
    struct history_steps *hs = malloc(sizeof(*hs));
    hs->first_clock = first_clock;
    hs->last_clock = last_clock;
    hs->start_position = sc->last_position;
    hs->interval = move->interval;
    hs->add = move->add;
    hs->add2 = move->add2;
    hs->shift = move->shift;
    hs->step_count = sc->sdir ? move->count : -move->count;
    sc->last_position += hs->step_count;
    list_add_head(&hs->node, &sc->history_list);
}

// Convert previously scheduled steps into commands for the mcu
static int
queue_flush(struct stepcompress *sc, uint64_t move_clock)
{
    if (sc->queue_pos >= sc->queue_next)
        return 0;
    while (sc->last_step_clock < move_clock) {
        struct step_move move = compress_bisect_count(sc);
        int ret = check_line(sc, move);
        if (ret)
            return ret;

        add_move(sc, sc->last_step_clock + move.first_step, &move);

        if (sc->queue_pos + move.count >= sc->queue_next) {
            sc->queue_pos = sc->queue_next = sc->queue;
            break;
        }
        sc->queue_pos += move.count;
    }
    calc_last_step_print_time(sc);
    return 0;
}

// Generate a queue_step for a step far in the future from the last step
static int
stepcompress_flush_far(struct stepcompress *sc, uint64_t abs_step_clock)
{
    uint64_t interval = abs_step_clock - sc->last_step_clock;
    struct step_move move = { interval, 1, 0, 0, 0, interval, interval };
    add_move(sc, abs_step_clock, &move);
    calc_last_step_print_time(sc);
    return 0;
}

// Send the set_next_step_dir command
static int
set_next_step_dir(struct stepcompress *sc, int sdir)
{
    if (sc->sdir == sdir)
        return 0;
    int ret = queue_flush(sc, UINT64_MAX);
    if (ret)
        return ret;
    sc->sdir = sdir;
    uint32_t msg[3] = {
        sc->set_next_step_dir_msgtag, sc->oid, sdir ^ sc->invert_sdir
    };
    struct queue_message *qm = message_alloc_and_encode(msg, 3);
    qm->req_clock = sc->last_step_clock;
    list_add_tail(&qm->node, &sc->msg_queue);
    return 0;
}

// Slow path for queue_append() - handle next step far in future
static int
queue_append_far(struct stepcompress *sc)
{
    uint64_t step_clock = sc->next_step_clock;
    sc->next_step_clock = 0;
    int ret = queue_flush(sc, step_clock - CLOCK_DIFF_MAX + 1);
    if (ret)
        return ret;
    if (step_clock >= sc->last_step_clock + CLOCK_DIFF_MAX)
        return stepcompress_flush_far(sc, step_clock);
    *sc->queue_next++ = step_clock;
    return 0;
}

// Slow path for queue_append() - expand the internal queue storage
static int
queue_append_extend(struct stepcompress *sc)
{
    if (sc->queue_next - sc->queue_pos > 65535 + 2000) {
        // No point in keeping more than 64K steps in memory
        uint32_t flush = (*(sc->queue_next-65535)
                          - (uint32_t)sc->last_step_clock);
        int ret = queue_flush(sc, sc->last_step_clock + flush);
        if (ret)
            return ret;
    }

    if (sc->queue_next >= sc->queue_end) {
        // Make room in the queue
        int in_use = sc->queue_next - sc->queue_pos;
        if (sc->queue_pos > sc->queue) {
            // Shuffle the internal queue to avoid having to allocate more ram
            memmove(sc->queue, sc->queue_pos, in_use * sizeof(*sc->queue));
        } else {
            // Expand the internal queue of step times
            int alloc = sc->queue_end - sc->queue;
            if (!alloc)
                alloc = QUEUE_START_SIZE;
            while (in_use >= alloc)
                alloc *= 2;
            sc->queue = realloc(sc->queue, alloc * sizeof(*sc->queue));
            sc->queue_end = sc->queue + alloc;
        }
        sc->queue_pos = sc->queue;
        sc->queue_next = sc->queue + in_use;
    }

    *sc->queue_next++ = sc->next_step_clock;
    sc->next_step_clock = 0;
    return 0;
}

// Add a step time to the queue (flushing the queue if needed)
static int
queue_append(struct stepcompress *sc)
{
    if (unlikely(sc->next_step_dir != sc->sdir)) {
        int ret = set_next_step_dir(sc, sc->next_step_dir);
        if (ret)
            return ret;
    }
    if (unlikely(sc->next_step_clock >= sc->last_step_clock + CLOCK_DIFF_MAX))
        return queue_append_far(sc);
    if (unlikely(sc->queue_next >= sc->queue_end))
        return queue_append_extend(sc);
    *sc->queue_next++ = sc->next_step_clock;
    sc->next_step_clock = 0;
    return 0;
}

#define SDS_FILTER_TIME .000750

// Add next step time
int
stepcompress_append(struct stepcompress *sc, int sdir
                    , double print_time, double step_time)
{
    // Calculate step clock
    double offset = print_time - sc->last_step_print_time;
    double rel_sc = (step_time + offset) * sc->mcu_freq;
    uint64_t step_clock = sc->last_step_clock + (uint64_t)rel_sc;
    // Flush previous pending step (if any)
    if (sc->next_step_clock) {
        if (unlikely(sdir != sc->next_step_dir)) {
            double diff = (int64_t)(step_clock - sc->next_step_clock);
            if (diff < SDS_FILTER_TIME * sc->mcu_freq) {
                // Rollback last step to avoid rapid step+dir+step
                sc->next_step_clock = 0;
                sc->next_step_dir = sdir;
                return 0;
            }
        }
        int ret = queue_append(sc);
        if (ret)
            return ret;
    }
    // Store this step as the next pending step
    sc->next_step_clock = step_clock;
    sc->next_step_dir = sdir;
    return 0;
}

// Commit next pending step (ie, do not allow a rollback)
int
stepcompress_commit(struct stepcompress *sc)
{
    if (sc->next_step_clock)
        return queue_append(sc);
    return 0;
}

// Flush pending steps
static int
stepcompress_flush(struct stepcompress *sc, uint64_t move_clock)
{
    if (sc->next_step_clock && move_clock >= sc->next_step_clock) {
        int ret = queue_append(sc);
        if (ret)
            return ret;
    }
    return queue_flush(sc, move_clock);
}

// Reset the internal state of the stepcompress object
int __visible
stepcompress_reset(struct stepcompress *sc, uint64_t last_step_clock)
{
    int ret = stepcompress_flush(sc, UINT64_MAX);
    if (ret)
        return ret;
    sc->last_step_clock = last_step_clock;
    sc->sdir = -1;
    calc_last_step_print_time(sc);
    return 0;
}

// Set last_position in the stepcompress object
int __visible
stepcompress_set_last_position(struct stepcompress *sc, uint64_t clock
                               , int64_t last_position)
{
    int ret = stepcompress_flush(sc, UINT64_MAX);
    if (ret)
        return ret;
    sc->last_position = last_position;

    // Add a marker to the history list
    struct history_steps *hs = malloc(sizeof(*hs));
    memset(hs, 0, sizeof(*hs));
    hs->first_clock = hs->last_clock = clock;
    hs->start_position = last_position;
    list_add_head(&hs->node, &sc->history_list);
    return 0;
}

// Search history of moves to find a past position at a given clock
int64_t __visible
stepcompress_find_past_position(struct stepcompress *sc, uint64_t clock)
{
    int64_t last_position = sc->last_position;
    struct history_steps *hs;
    list_for_each_entry(hs, &sc->history_list, node) {
        if (clock < hs->first_clock) {
            last_position = hs->start_position;
            continue;
        }
        if (clock >= hs->last_clock)
            return hs->start_position + hs->step_count;
        int64_t ticks = clock - hs->first_clock;
        int64_t interval = hs->interval, add = hs->add, add2 = hs->add2;
        int count = hs->step_count, shift = hs->shift;
        if (count < 0) count = -count;
        if (shift <= 0) {
            int mul = 1 << -shift;
            interval <<= -shift;
            add *= mul;
            add2 *= mul;
        } else {
            ticks *= 1 << shift;
        }
        // When clock == hs->first_clock, offset == 1
        ticks += interval;
        int left = 0, right = count;
        while (right - left > 1) {
            int cnt = (left + right) / 2;
            int64_t step_clock = cnt * interval + cnt * (cnt-1) / 2 * add +
                cnt * (cnt-1) * (cnt-2) / 6 * add2;
            if (step_clock <= ticks) left = cnt;
            else right = cnt;
        }
        int64_t clock_left = left * interval + left * (left-1) / 2 * add +
            left * (left-1) * (left-2) / 6 * add2;
        int64_t clock_right = right * interval + right * (right-1) / 2 * add +
            right * (right-1) * (right-2) / 6 * add2;
        int offset = ticks - clock_left <= clock_right - ticks ? left : right;
        if (hs->step_count < 0)
            return hs->start_position - offset;
        return hs->start_position + offset;
    }
    return last_position;
}

// Queue an mcu command to go out in order with stepper commands
int __visible
stepcompress_queue_msg(struct stepcompress *sc, uint32_t *data, int len)
{
    int ret = stepcompress_flush(sc, UINT64_MAX);
    if (ret)
        return ret;

    struct queue_message *qm = message_alloc_and_encode(data, len);
    qm->req_clock = sc->last_step_clock;
    list_add_tail(&qm->node, &sc->msg_queue);
    return 0;
}

// Return history of queue_step commands
int __visible
stepcompress_extract_old(struct stepcompress *sc, struct pull_history_steps *p
                         , int max, uint64_t start_clock, uint64_t end_clock)
{
    int res = 0;
    struct history_steps *hs;
    list_for_each_entry(hs, &sc->history_list, node) {
        if (start_clock >= hs->last_clock || res >= max)
            break;
        if (end_clock <= hs->first_clock)
            continue;
        p->first_clock = hs->first_clock;
        p->last_clock = hs->last_clock;
        p->start_position = hs->start_position;
        p->step_count = hs->step_count;
        p->interval = hs->interval;
        p->add = hs->add;
        p->add2 = hs->add2;
        p->shift = hs->shift;
        p++;
        res++;
    }
    return res;
}


/****************************************************************
 * Step compress synchronization
 ****************************************************************/

// The steppersync object is used to synchronize the output of mcu
// step commands.  The mcu can only queue a limited number of step
// commands - this code tracks when items on the mcu step queue become
// free so that new commands can be transmitted.  It also ensures the
// mcu step queue is ordered between steppers so that no stepper
// starves the other steppers of space in the mcu step queue.

struct steppersync {
    // Serial port
    struct serialqueue *sq;
    struct command_queue *cq;
    // Storage for associated stepcompress objects
    struct stepcompress **sc_list;
    int sc_num;
    // Storage for list of pending move clocks
    uint64_t *move_clocks;
    int num_move_clocks;
};

// Allocate a new 'steppersync' object
struct steppersync * __visible
steppersync_alloc(struct serialqueue *sq, struct stepcompress **sc_list
                  , int sc_num, int move_num)
{
    struct steppersync *ss = malloc(sizeof(*ss));
    memset(ss, 0, sizeof(*ss));
    ss->sq = sq;
    ss->cq = serialqueue_alloc_commandqueue();

    ss->sc_list = malloc(sizeof(*sc_list)*sc_num);
    memcpy(ss->sc_list, sc_list, sizeof(*sc_list)*sc_num);
    ss->sc_num = sc_num;

    ss->move_clocks = malloc(sizeof(*ss->move_clocks)*move_num);
    memset(ss->move_clocks, 0, sizeof(*ss->move_clocks)*move_num);
    ss->num_move_clocks = move_num;

    return ss;
}

// Free memory associated with a 'steppersync' object
void __visible
steppersync_free(struct steppersync *ss)
{
    if (!ss)
        return;
    free(ss->sc_list);
    free(ss->move_clocks);
    serialqueue_free_commandqueue(ss->cq);
    free(ss);
}

// Set the conversion rate of 'print_time' to mcu clock
void __visible
steppersync_set_time(struct steppersync *ss, double time_offset
                     , double mcu_freq)
{
    int i;
    for (i=0; i<ss->sc_num; i++) {
        struct stepcompress *sc = ss->sc_list[i];
        stepcompress_set_time(sc, time_offset, mcu_freq);
    }
}

// Implement a binary heap algorithm to track when the next available
// 'struct move' in the mcu will be available
static void
heap_replace(struct steppersync *ss, uint64_t req_clock)
{
    uint64_t *mc = ss->move_clocks;
    int nmc = ss->num_move_clocks, pos = 0;
    for (;;) {
        int child1_pos = 2*pos+1, child2_pos = 2*pos+2;
        uint64_t child2_clock = child2_pos < nmc ? mc[child2_pos] : UINT64_MAX;
        uint64_t child1_clock = child1_pos < nmc ? mc[child1_pos] : UINT64_MAX;
        if (req_clock <= child1_clock && req_clock <= child2_clock) {
            mc[pos] = req_clock;
            break;
        }
        if (child1_clock < child2_clock) {
            mc[pos] = child1_clock;
            pos = child1_pos;
        } else {
            mc[pos] = child2_clock;
            pos = child2_pos;
        }
    }
}

// Find and transmit any scheduled steps prior to the given 'move_clock'
int __visible
steppersync_flush(struct steppersync *ss, uint64_t move_clock)
{
    // Flush each stepcompress to the specified move_clock
    int i;
    for (i=0; i<ss->sc_num; i++) {
        int ret = stepcompress_flush(ss->sc_list[i], move_clock);
        if (ret)
            return ret;
    }

    // Order commands by the reqclock of each pending command
    struct list_head msgs;
    list_init(&msgs);
    for (;;) {
        // Find message with lowest reqclock
        uint64_t req_clock = MAX_CLOCK;
        struct queue_message *qm = NULL;
        for (i=0; i<ss->sc_num; i++) {
            struct stepcompress *sc = ss->sc_list[i];
            if (!list_empty(&sc->msg_queue)) {
                struct queue_message *m = list_first_entry(
                    &sc->msg_queue, struct queue_message, node);
                if (m->req_clock < req_clock) {
                    qm = m;
                    req_clock = m->req_clock;
                }
            }
        }
        if (!qm || (qm->min_clock && req_clock > move_clock))
            break;

        uint64_t next_avail = ss->move_clocks[0];
        if (qm->min_clock)
            // The qm->min_clock field is overloaded to indicate that
            // the command uses the 'move queue' and to store the time
            // that move queue item becomes available.
            heap_replace(ss, qm->min_clock);
        // Reset the min_clock to its normal meaning (minimum transmit time)
        qm->min_clock = next_avail;

        // Batch this command
        list_del(&qm->node);
        list_add_tail(&qm->node, &msgs);
    }

    // Transmit commands
    if (!list_empty(&msgs))
        serialqueue_send_batch(ss->sq, ss->cq, &msgs);
    return 0;
}
