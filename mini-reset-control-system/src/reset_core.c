/* reset_core.c - Core Implementations for Reset Control Systems
 *
 * Implements the lifecycle management, reset logic, and state
 * manipulation functions declared in reset_core.h.
 *
 * Each function implements an independent knowledge point:
 *   - reset_base_create: L1/L3 memory allocation for linear systems
 *   - reset_clegg_create: L1 construction of Clegg integrator
 *   - reset_fore_create: L1 construction of FORE element
 *   - reset_check_and_apply: L2 zero-crossing detection + reset logic
 *   - reset_set_dwell_time: L2 Zeno prevention mechanism
 *   - reset_eval_flow_deriv: L3 continuous flow evaluation
 */

#include "reset_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ================================================================
 * L1/L3: Lifecycle - Linear Base System
 * ================================================================ */

ResetLinearBase* reset_base_create(int n, int m, int p)
{
    if (n <= 0 || m <= 0 || p <= 0) return NULL;

    ResetLinearBase *base = (ResetLinearBase*)malloc(sizeof(ResetLinearBase));
    if (!base) return NULL;

    base->n = n;
    base->m = m;
    base->p = p;

    /* Allocate all matrices, zero-initialized */
    size_t sA = (size_t)n * n * sizeof(double);
    size_t sB = (size_t)n * m * sizeof(double);
    size_t sC = (size_t)p * n * sizeof(double);
    size_t sD = (size_t)p * m * sizeof(double);

    base->A = (double*)malloc(sA);
    base->B = (double*)malloc(sB);
    base->C = (double*)malloc(sC);
    base->D = (double*)malloc(sD);

    if (!base->A || !base->B || !base->C || !base->D) {
        free(base->A); free(base->B); free(base->C); free(base->D);
        free(base);
        return NULL;
    }

    memset(base->A, 0, sA);
    memset(base->B, 0, sB);
    memset(base->C, 0, sC);
    memset(base->D, 0, sD);

    return base;
}

ResetLinearBase* reset_base_clone(const ResetLinearBase *base)
{
    if (!base) return NULL;

    ResetLinearBase *copy = reset_base_create(base->n, base->m, base->p);
    if (!copy) return NULL;

    size_t sA = (size_t)base->n * base->n * sizeof(double);
    size_t sB = (size_t)base->n * base->m * sizeof(double);
    size_t sC = (size_t)base->p * base->n * sizeof(double);
    size_t sD = (size_t)base->p * base->m * sizeof(double);

    memcpy(copy->A, base->A, sA);
    memcpy(copy->B, base->B, sB);
    memcpy(copy->C, base->C, sC);
    memcpy(copy->D, base->D, sD);

    return copy;
}

void reset_base_free(ResetLinearBase *base)
{
    if (!base) return;
    free(base->A);
    free(base->B);
    free(base->C);
    free(base->D);
    free(base);
}

/* ================================================================
 * L1/L2: Lifecycle - Reset System
 * ================================================================ */

ResetSystem* reset_sys_create(int nc)
{
    if (nc <= 0) return NULL;

    ResetSystem *rsys = (ResetSystem*)malloc(sizeof(ResetSystem));
    if (!rsys) return NULL;
    memset(rsys, 0, sizeof(ResetSystem));

    rsys->nc = nc;
    rsys->x_c = (double*)calloc((size_t)nc, sizeof(double));
    if (!rsys->x_c) { free(rsys); return NULL; }

    /* Allocate sub-structures */
    rsys->flow = (ResetLinearBase*)malloc(sizeof(ResetLinearBase));
    rsys->jump = (ResetJumpMap*)malloc(sizeof(ResetJumpMap));
    rsys->cond = (ResetCondition*)malloc(sizeof(ResetCondition));
    rsys->ratio = (ResetRatio*)malloc(sizeof(ResetRatio));
    rsys->band = (ResetBand*)malloc(sizeof(ResetBand));
    rsys->interval = (ResetInterval*)malloc(sizeof(ResetInterval));

    if (!rsys->flow || !rsys->jump || !rsys->cond ||
        !rsys->ratio || !rsys->band || !rsys->interval) {
        reset_sys_free(rsys);
        return NULL;
    }

    /* Zero-initialize all sub-structures */
    memset(rsys->flow, 0, sizeof(ResetLinearBase));
    memset(rsys->jump, 0, sizeof(ResetJumpMap));
    memset(rsys->cond, 0, sizeof(ResetCondition));
    memset(rsys->ratio, 0, sizeof(ResetRatio));
    memset(rsys->band, 0, sizeof(ResetBand));
    memset(rsys->interval, 0, sizeof(ResetInterval));

    /* Set default condition parameters */
    rsys->cond->trig_type = RESET_TRIGGER_ZERO_CROSSING;
    rsys->cond->epsilon   = 1e-10;
    rsys->cond->t_min_dwell = 0.0;

    /* Set default band: wide open (always armed) */
    rsys->band->band_lower = 0.0;
    rsys->band->band_upper = 1e100;
    rsys->band->band_hyst  = 0.0;

    /* Set default ratio: full reset */
    rsys->ratio->rho = 0.0;
    rsys->ratio->rho_min = 0.0;
    rsys->ratio->rho_max = 1.0;
    rsys->ratio->is_time_varying = false;

    /* Initialize interval tracking */
    rsys->interval->t_last         = -1.0;
    rsys->interval->t_min_interval = 1e100;
    rsys->interval->t_max_interval = 0.0;
    rsys->interval->t_avg_interval = 0.0;
    rsys->interval->n_resets       = 0;

    rsys->reset_armed = true;
    rsys->e_prev = 0.0;
    rsys->t = 0.0;
    rsys->dt_mode = false;
    rsys->ts = 0.001;

    return rsys;
}

void reset_sys_free(ResetSystem *rsys)
{
    if (!rsys) return;
    free(rsys->x_c);
    if (rsys->flow) {
        reset_base_free(rsys->flow);
    }
    if (rsys->jump) {
        free(rsys->jump->Ar);
        free(rsys->jump->Br);
        free(rsys->jump);
    }
    if (rsys->cond) {
        free(rsys->cond->H);
        free(rsys->cond);
    }
    free(rsys->ratio);
    free(rsys->band);
    free(rsys->interval);
    free(rsys);
}

/* ================================================================
 * L1: Clegg Integrator Construction [Cle58]
 * ================================================================ */

ResetSystem* reset_clegg_create(void)
{
    ResetSystem *rsys = reset_sys_create(1);
    if (!rsys) return NULL;

    /* Allocate flow matrices: A=0, B=1, C=1, D=0 */
    ResetLinearBase *flow = rsys->flow;
    flow->n = 1; flow->m = 1; flow->p = 1;

    flow->A = (double*)calloc(1, sizeof(double)); /* A = [0] */
    flow->B = (double*)calloc(1, sizeof(double));
    flow->C = (double*)calloc(1, sizeof(double));
    flow->D = (double*)calloc(1, sizeof(double));

    if (!flow->A || !flow->B || !flow->C || !flow->D) {
        reset_sys_free(rsys); return NULL;
    }
    flow->B[0] = 1.0; /* B = [1] */
    flow->C[0] = 1.0; /* C = [1] */
    /* D remains 0, A remains 0 */

    /* Allocate jump map: full reset Ar=[0], Br=[0] */
    rsys->jump->nc = 1;
    rsys->jump->Ar = (double*)calloc(1, sizeof(double));
    rsys->jump->Br = (double*)calloc(1, sizeof(double));
    if (!rsys->jump->Ar || !rsys->jump->Br) {
        reset_sys_free(rsys); return NULL;
    }

    /* Set reset condition to zero-crossing */
    rsys->cond->trig_type = RESET_TRIGGER_ZERO_CROSSING;
    rsys->ratio->rho = 0.0; /* full reset */

    return rsys;
}

/* ================================================================
 * L1: FORE Construction [HR75]
 * ================================================================ */

ResetSystem* reset_fore_create(double K, double tau, double rho)
{
    if (tau <= 0.0 || rho < 0.0 || rho >= 1.0) return NULL;

    ResetSystem *rsys = reset_sys_create(1);
    if (!rsys) return NULL;

    ResetLinearBase *flow = rsys->flow;
    flow->n = 1; flow->m = 1; flow->p = 1;

    flow->A = (double*)calloc(1, sizeof(double));
    flow->B = (double*)calloc(1, sizeof(double));
    flow->C = (double*)calloc(1, sizeof(double));
    flow->D = (double*)calloc(1, sizeof(double));

    if (!flow->A || !flow->B || !flow->C || !flow->D) {
        reset_sys_free(rsys); return NULL;
    }

    /* A = -1/tau, B = K/tau, C = 1, D = 0 */
    flow->A[0] = -1.0 / tau;
    flow->B[0] = K / tau;
    flow->C[0] = 1.0;
    /* D = 0 (already zero) */

    /* Jump map: Ar = [rho], Br = [0] */
    rsys->jump->nc = 1;
    rsys->jump->Ar = (double*)calloc(1, sizeof(double));
    rsys->jump->Br = (double*)calloc(1, sizeof(double));
    if (!rsys->jump->Ar || !rsys->jump->Br) {
        reset_sys_free(rsys); return NULL;
    }
    rsys->jump->Ar[0] = rho;

    /* Configure reset condition and ratio */
    rsys->cond->trig_type = RESET_TRIGGER_ZERO_CROSSING;
    rsys->ratio->rho = rho;

    return rsys;
}

/* ================================================================
 * L2: Reset Ratio Configuration
 * ================================================================ */

ResetResult reset_set_ratio(ResetSystem *rsys, double rho)
{
    if (!rsys) return RESET_NULL_PTR;
    if (rho < 0.0 || rho >= 1.0) return RESET_NOT_ARMED;

    rsys->ratio->rho = rho;

    /* Sync Ar matrix if it exists: Ar = rho * I */
    if (rsys->jump && rsys->jump->Ar && rsys->jump->nc > 0) {
        int nc = rsys->jump->nc;
        for (int i = 0; i < nc * nc; i++) {
            rsys->jump->Ar[i] = 0.0;
        }
        for (int i = 0; i < nc; i++) {
            rsys->jump->Ar[i * nc + i] = rho;
        }
    }

    return RESET_OK;
}

/* ================================================================
 * L2: Zero-Crossing Detection and Reset Application
 *
 * This is the core reset control algorithm. At each step:
 *   1. Check if error signal crossed zero
 *   2. If armed and ZC detected, apply jump map
 *   3. Update interval statistics
 *   4. Update arm/disarm state based on reset band
 *
 * Ref: [ZCH00] zero-crossing resetting law
 *      [BB12] Section 4.2 on reset band implementation
 * ================================================================ */

ResetResult reset_check_and_apply(ResetSystem *rsys, double e, double e_prev)
{
    if (!rsys) return RESET_NULL_PTR;

    rsys->e_prev = e;

    /* Update arm/disarm based on reset band */
    double abs_e = fabs(e);
    if (abs_e > rsys->band->band_upper) {
        rsys->reset_armed = false;
    } else if (abs_e < rsys->band->band_lower) {
        rsys->reset_armed = true;
    }

    /* Detect zero crossing: sign(e) != sign(e_prev) */
    bool zc_detected = false;
    switch (rsys->cond->trig_type) {
    case RESET_TRIGGER_ZERO_CROSSING:
        if (e * e_prev < 0.0 && fabs(e_prev) > rsys->cond->epsilon)
            zc_detected = true;
        break;
    case RESET_TRIGGER_RISING:
        if (e > 0.0 && e_prev < 0.0 && fabs(e_prev) > rsys->cond->epsilon)
            zc_detected = true;
        break;
    case RESET_TRIGGER_FALLING:
        if (e < 0.0 && e_prev > 0.0 && fabs(e_prev) > rsys->cond->epsilon)
            zc_detected = true;
        break;
    case RESET_TRIGGER_PERIODIC:
        /* Periodic reset handled elsewhere */
        break;
    case RESET_TRIGGER_STATE_COND:
        if (rsys->cond->H && rsys->x_c) {
            double g = 0.0;
            for (int i = 0; i < rsys->nc; i++)
                g += rsys->cond->H[i] * rsys->x_c[i];
            if (g <= rsys->cond->threshold)
                zc_detected = true;
        }
        break;
    case RESET_TRIGGER_TIME_DEP:
        /* Time-dependent handled by caller using reset time sequence */
        break;
    default:
        break;
    }

    if (!zc_detected) {
        return RESET_NO_CROSSING;
    }

    if (!rsys->reset_armed) {
        return RESET_NOT_ARMED;
    }

    /* Check dwell time constraint */
    if (rsys->interval->n_resets > 0 &&
        rsys->t - rsys->interval->t_last < rsys->cond->t_min_dwell) {
        return RESET_DWELL_BLOCKED;
    }

    /* Apply the reset (jump map): x_c^+ = Ar * x_c^- + Br * e */
    double *x_new = (double*)calloc((size_t)rsys->nc, sizeof(double));
    if (!x_new) return RESET_NULL_PTR;

    /* Compute x_new = Ar * x_old */
    for (int i = 0; i < rsys->nc; i++) {
        double sum = 0.0;
        for (int j = 0; j < rsys->jump->nc; j++) {
            sum += rsys->jump->Ar[i * rsys->jump->nc + j] * rsys->x_c[j];
        }
        x_new[i] = sum;
    }
    /* Add Br * e */
    if (rsys->jump->Br) {
        for (int i = 0; i < rsys->nc; i++) {
            x_new[i] += rsys->jump->Br[i] * e;
        }
    }

    /* Update state and statistics */
    memcpy(rsys->x_c, x_new, (size_t)rsys->nc * sizeof(double));
    free(x_new);

    /* Update interval statistics */
    if (rsys->interval->n_resets > 0) {
        double interval = rsys->t - rsys->interval->t_last;
        if (interval < rsys->interval->t_min_interval)
            rsys->interval->t_min_interval = interval;
        if (interval > rsys->interval->t_max_interval)
            rsys->interval->t_max_interval = interval;
        /* Exponential moving average for avg */
        double alpha = 0.1;
        rsys->interval->t_avg_interval =
            alpha * interval + (1.0 - alpha) * rsys->interval->t_avg_interval;
    }
    rsys->interval->t_last = rsys->t;
    rsys->interval->n_resets++;

    /* After reset, disarm until error leaves band */
    if (fabs(e) < rsys->band->band_lower) {
        rsys->reset_armed = false; /* prevent immediate re-trigger */
    }

    return RESET_OK;
}

/* ================================================================
 * L2: Manual Reset, Configuration, and Utilities
 * ================================================================ */

ResetResult reset_apply_manual(ResetSystem *rsys, double t, double e)
{
    if (!rsys || !rsys->x_c) return RESET_NULL_PTR;

    double *x_new = (double*)calloc((size_t)rsys->nc, sizeof(double));
    if (!x_new) return RESET_NULL_PTR;

    for (int i = 0; i < rsys->nc; i++) {
        double sum = 0.0;
        for (int j = 0; j < rsys->jump->nc; j++) {
            sum += rsys->jump->Ar[i * rsys->jump->nc + j] * rsys->x_c[j];
        }
        x_new[i] = sum;
    }
    if (rsys->jump->Br) {
        for (int i = 0; i < rsys->nc; i++) {
            x_new[i] += rsys->jump->Br[i] * e;
        }
    }

    memcpy(rsys->x_c, x_new, (size_t)rsys->nc * sizeof(double));
    free(x_new);

    if (rsys->interval->n_resets > 0) {
        double interval = t - rsys->interval->t_last;
        if (interval < rsys->interval->t_min_interval)
            rsys->interval->t_min_interval = interval;
        if (interval > rsys->interval->t_max_interval)
            rsys->interval->t_max_interval = interval;
    }
    rsys->interval->t_last = t;
    rsys->interval->n_resets++;

    return RESET_OK;
}

void reset_config_zc_trigger(ResetSystem *rsys, ResetTriggerType trig_type,
                              double band_lower, double band_upper, double hyst)
{
    if (!rsys) return;
    rsys->cond->trig_type = trig_type;
    if (band_lower >= 0.0) rsys->band->band_lower = band_lower;
    if (band_upper > band_lower) rsys->band->band_upper = band_upper;
    if (hyst >= 0.0) rsys->band->band_hyst = hyst;
}

void reset_set_dwell_time(ResetSystem *rsys, double t_min_dwell)
{
    if (!rsys || t_min_dwell < 0.0) return;
    rsys->cond->t_min_dwell = t_min_dwell;
}

void reset_set_time_schedule(ResetSystem *rsys, const double *times, int n)
{
    if (!rsys || !times || n <= 0) return;
    /* Note: ResetTimeSequence is defined in core.h but stored externally
     * for time-dependent reset. This sets the trigger type. */
    rsys->cond->trig_type = RESET_TRIGGER_TIME_DEP;
    /* The actual schedule is managed by the caller/simulator */
}

const ResetInterval* reset_get_interval_stats(const ResetSystem *rsys)
{
    if (!rsys) return NULL;
    return rsys->interval;
}

void reset_set_linear_base(ResetSystem *rsys, int n, int m, int p,
                            double *A, double *B, double *C, double *D)
{
    if (!rsys || !rsys->flow) return;

    /* Free old matrices */
    free(rsys->flow->A); free(rsys->flow->B);
    free(rsys->flow->C); free(rsys->flow->D);

    rsys->flow->n = n;
    rsys->flow->m = m;
    rsys->flow->p = p;
    rsys->flow->A = A;
    rsys->flow->B = B;
    rsys->flow->C = C;
    rsys->flow->D = D;
}

void reset_set_jump_map(ResetSystem *rsys, double *Ar, double *Br, int nc)
{
    if (!rsys || !rsys->jump) return;

    free(rsys->jump->Ar);
    free(rsys->jump->Br);

    rsys->jump->nc = nc;
    rsys->jump->Ar = Ar;
    rsys->jump->Br = Br;
}

/* ================================================================
 * L3: Flow Derivative and Output Evaluation
 *
 * These functions implement the continuous dynamics of the
 * hybrid system: dx/dt = A x + B u, y = C x + D u.
 *
 * Used by the simulation engine for numerical integration.
 * ================================================================ */

int reset_eval_flow_deriv(const ResetLinearBase *base,
                           const double *x, const double *u, double *dxdt)
{
    if (!base || !x || !dxdt) return -1;
    if (!base->A || !base->B) return -1;

    int n = base->n;
    int m = base->m;

    /* dxdt = A * x */
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += base->A[i * n + j] * x[j];
        }
        dxdt[i] = sum;
    }

    /* dxdt += B * u */
    if (u && m > 0) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < m; j++) {
                sum += base->B[i * m + j] * u[j];
            }
            dxdt[i] += sum;
        }
    }

    return 0;
}

int reset_eval_output(const ResetLinearBase *base,
                       const double *x, const double *u, double *y)
{
    if (!base || !x || !y) return -1;
    if (!base->C) return -1;

    int n = base->n;
    int m = base->m;
    int p = base->p;

    /* y = C * x */
    for (int i = 0; i < p; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += base->C[i * n + j] * x[j];
        }
        y[i] = sum;
    }

    /* y += D * u */
    if (u && base->D && m > 0) {
        for (int i = 0; i < p; i++) {
            double sum = 0.0;
            for (int j = 0; j < m; j++) {
                sum += base->D[i * m + j] * u[j];
            }
            y[i] += sum;
        }
    }

    return 0;
}