/*
 * impulsive_types.c -- Core type lifecycle for impulsive systems.
 * Implements creation, destruction, validation, and manipulation of
 * ImpTimeSeq, ImpSystem, ImpSolution.  Key formulas:
 *   Periodic: tau_k = t0 + k*period
 *   Dwell:    dwell_k = tau_{k+1} - tau_k
 * Complexity: O(1) create/free, O(log n) find_index.
 */
#include "impulsive_types.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

static int cmp_dbl(const void *a, const void *b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da < db) ? -1 : (da > db) ? 1 : 0;
}

/* ---- ImpTimeSeq ---- */

ImpTimeSeq* imp_time_seq_create(int capacity, double t0, double T)
{
    if (capacity < 1 || capacity > IMP_MAX_IMPULSES) return NULL;
    if (!isfinite(t0) || !isfinite(T) || T <= t0) return NULL;
    ImpTimeSeq *seq = (ImpTimeSeq*)calloc(1, sizeof(ImpTimeSeq));
    if (!seq) return NULL;
    seq->times = (double*)malloc((size_t)capacity * sizeof(double));
    if (!seq->times) { free(seq); return NULL; }
    seq->capacity = capacity; seq->count = 0;
    seq->t0 = t0; seq->T = T;
    seq->is_periodic = false; seq->period = 0.0;
    return seq;
}

void imp_time_seq_free(ImpTimeSeq *seq)
{
    if (!seq) return;
    free(seq->times); free(seq);
}

int imp_time_seq_add(ImpTimeSeq *seq, double tau)
{
    if (!seq || !seq->times || !isfinite(tau)) return -1;
    if (tau < seq->t0 - IMP_EPS || tau > seq->T + IMP_EPS) return -1;
    if (seq->count > 0 && tau <= seq->times[seq->count - 1] + IMP_EPS) return -1;
    if (seq->count >= seq->capacity) {
        int new_cap = seq->capacity * 2;
        if (new_cap > IMP_MAX_IMPULSES) return -2;
        double *nt = (double*)realloc(seq->times, (size_t)new_cap * sizeof(double));
        if (!nt) return -2;
        seq->times = nt; seq->capacity = new_cap;
    }
    seq->times[seq->count++] = tau;
    return 0;
}

int imp_time_seq_complete(ImpTimeSeq *seq)
{
    if (!seq || !seq->times || seq->count < 1) return -1;
    qsort(seq->times, (size_t)seq->count, sizeof(double), cmp_dbl);
    for (int i = 1; i < seq->count; i++)
        if (seq->times[i] <= seq->times[i-1] + IMP_EPS) return -1;
    if (seq->times[0] < seq->t0 - IMP_EPS) return -1;
    if (seq->times[seq->count-1] > seq->T + IMP_EPS) return -1;
    return 0;
}

ImpTimeSeq* imp_time_seq_create_periodic(double t0, double T, double period)
{
    if (!isfinite(t0) || !isfinite(T) || !isfinite(period)) return NULL;
    if (period <= IMP_EPS || T <= t0) return NULL;
    int count = (int)floor((T - t0) / period) + 1;
    if (count < 1 || count > IMP_MAX_IMPULSES) return NULL;
    ImpTimeSeq *seq = imp_time_seq_create(count, t0, T);
    if (!seq) return NULL;
    for (int k = 0; k < count; k++)
        seq->times[k] = t0 + (double)k * period;
    seq->count = count; seq->is_periodic = true; seq->period = period;
    return seq;
}

ImpTimeSeq* imp_time_seq_create_aperiodic(const double *times, int count,
                                            double t0, double T)
{
    if (!times || count < 1 || count > IMP_MAX_IMPULSES) return NULL;
    if (!isfinite(t0) || !isfinite(T) || T <= t0) return NULL;
    for (int i = 0; i < count; i++) {
        if (!isfinite(times[i])) return NULL;
        if (i > 0 && times[i] <= times[i-1] + IMP_EPS) return NULL;
    }
    if (times[0] < t0 - IMP_EPS || times[count-1] > T + IMP_EPS) return NULL;
    ImpTimeSeq *seq = imp_time_seq_create(count, t0, T);
    if (!seq) return NULL;
    memcpy(seq->times, times, (size_t)count * sizeof(double));
    seq->count = count;
    return seq;
}

int imp_time_seq_get_dwell_stats(const ImpTimeSeq *seq, ImpDwellStats *stats)
{
    if (!seq || !stats || seq->count < 2) return -1;
    stats->num_dwells = seq->count - 1;
    stats->min_dwell = IMP_HUGE; stats->max_dwell = 0.0; stats->total_time = 0.0;
    for (int i = 0; i < stats->num_dwells; i++) {
        double dwell = seq->times[i+1] - seq->times[i];
        if (dwell < stats->min_dwell) stats->min_dwell = dwell;
        if (dwell > stats->max_dwell) stats->max_dwell = dwell;
        stats->total_time += dwell;
    }
    stats->avg_dwell = stats->total_time / (double)stats->num_dwells;
    return 0;
}

int imp_time_seq_find_index(const ImpTimeSeq *seq, double t)
{
    if (!seq || seq->count < 1) return -1;
    if (t < seq->times[0] - IMP_EPS) return -1;
    if (t >= seq->times[seq->count - 1] - IMP_EPS) return seq->count - 1;
    int lo = 0, hi = seq->count - 1;
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (seq->times[mid] <= t + IMP_EPS) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

ImpTimeSeq* imp_time_seq_clone(const ImpTimeSeq *src)
{
    if (!src) return NULL;
    ImpTimeSeq *cl = imp_time_seq_create(src->capacity, src->t0, src->T);
    if (!cl) return NULL;
    memcpy(cl->times, src->times, (size_t)src->count * sizeof(double));
    cl->count = src->count; cl->is_periodic = src->is_periodic;
    cl->period = src->period;
    return cl;
}

/* ---- ImpSystem ---- */

ImpSystem* imp_system_create(int n, ImpVectorField f, ImpJumpMap I,
                              ImpTimeSeq *seq, void *ctx)
{
    if (n < 1 || n > IMP_MAX_DIM || !f || !I || !seq) return NULL;
    ImpSystem *sys = (ImpSystem*)calloc(1, sizeof(ImpSystem));
    if (!sys) return NULL;
    sys->n = n; sys->f = f; sys->I = I; sys->seq = seq;
    sys->trigger = IMP_TRIGGER_TIME;
    sys->guard = NULL; sys->event = NULL; sys->ctx = ctx;
    snprintf(sys->name, sizeof(sys->name), "imp_sys_n%d", n);
    return sys;
}

void imp_system_free(ImpSystem *sys) { free(sys); }

bool imp_system_validate(const ImpSystem *sys)
{
    if (!sys || sys->n < 1 || sys->n > IMP_MAX_DIM) return false;
    if (!sys->f || !sys->I || !sys->seq || sys->seq->count < 1) return false;
    if (sys->trigger == IMP_TRIGGER_STATE && !sys->guard) return false;
    if (sys->trigger == IMP_TRIGGER_EVENT && !sys->event) return false;
    return true;
}

/* ---- ImpSolution ---- */

ImpSolution* imp_solution_create(int npts_cap, int nimp_cap, int n)
{
    if (npts_cap < 1 || nimp_cap < 0 || n < 1 || n > IMP_MAX_DIM) return NULL;
    ImpSolution *sol = (ImpSolution*)calloc(1, sizeof(ImpSolution));
    if (!sol) return NULL;
    sol->t = (double*)malloc((size_t)npts_cap * sizeof(double));
    sol->x = (double*)malloc((size_t)npts_cap * n * sizeof(double));
    if (!sol->t || !sol->x) { imp_solution_free(sol); return NULL; }
    if (nimp_cap > 0) {
        size_t jb = (size_t)nimp_cap * n * sizeof(double);
        sol->x_jump_before = (double*)malloc(jb);
        sol->x_jump_after  = (double*)malloc(jb);
        if (!sol->x_jump_before || !sol->x_jump_after) {
            imp_solution_free(sol); return NULL; }
    }
    sol->npts_cap = npts_cap; sol->nimp_cap = nimp_cap;
    sol->n = n; sol->npts = 0; sol->nimp = 0;
    return sol;
}

void imp_solution_free(ImpSolution *sol)
{
    if (!sol) return;
    free(sol->t); free(sol->x);
    free(sol->x_jump_before); free(sol->x_jump_after);
    free(sol);
}

int imp_solution_resize(ImpSolution *sol)
{
    if (!sol) return -1;
    int nnc = sol->npts_cap * 2;
    double *nt = (double*)realloc(sol->t, (size_t)nnc * sizeof(double));
    double *nx = (double*)realloc(sol->x, (size_t)nnc * sol->n * sizeof(double));
    if (!nt || !nx) return -2;
    sol->t = nt; sol->x = nx; sol->npts_cap = nnc;
    if (sol->nimp_cap > 0) {
        int nic = sol->nimp_cap * 2;
        size_t jb = (size_t)nic * sol->n * sizeof(double);
        double *njb = (double*)realloc(sol->x_jump_before, jb);
        double *nja = (double*)realloc(sol->x_jump_after, jb);
        if (!njb || !nja) return -2;
        sol->x_jump_before = njb; sol->x_jump_after = nja;
        sol->nimp_cap = nic;
    }
    return 0;
}

int imp_solution_add_point(ImpSolution *sol, double t, const double *x)
{
    if (!sol || !x) return -1;
    if (sol->npts >= sol->npts_cap && imp_solution_resize(sol) != 0) return -2;
    sol->t[sol->npts] = t;
    memcpy(&sol->x[(size_t)sol->npts * sol->n], x, (size_t)sol->n * sizeof(double));
    sol->npts++; sol->T = t;
    return 0;
}

int imp_solution_add_jump(ImpSolution *sol, const double *x_before,
                           const double *x_after)
{
    if (!sol || !x_before || !x_after || sol->nimp_cap == 0) return -1;
    if (sol->nimp >= sol->nimp_cap) {
        int nc = sol->nimp_cap * 2;
        size_t jb = (size_t)nc * sol->n * sizeof(double);
        double *njb = (double*)realloc(sol->x_jump_before, jb);
        double *nja = (double*)realloc(sol->x_jump_after, jb);
        if (!njb || !nja) return -2;
        sol->x_jump_before = njb; sol->x_jump_after = nja;
        sol->nimp_cap = nc;
    }
    size_t off = (size_t)sol->nimp * sol->n;
    memcpy(&sol->x_jump_before[off], x_before, (size_t)sol->n * sizeof(double));
    memcpy(&sol->x_jump_after[off], x_after, (size_t)sol->n * sizeof(double));
    sol->nimp++;
    return 0;
}
