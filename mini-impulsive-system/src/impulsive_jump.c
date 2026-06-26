/*
 * impulsive_jump.c -- Impulsive reset / jump map implementations
 *
 * Implements: linear, affine, projection, hard-reset, impulsive
 * control, nonlinear, time-varying, and composed jump maps.
 *
 * Each jump map computes:  x(tau_k^+) = x(tau_k^-) + I_k(x(tau_k^-))
 *
 * Key formulas:
 *   Linear:   x^+ = (I+J) * x^-
 *   Affine:   x^+ = (I+J) * x^- + d
 *   Project:  x_i^+ = clamp(x_i^-, lb_i, ub_i)
 *   Control:  x^+ = x^- + B * u
 *
 * References:
 *   Yang (2001) "Impulsive Control Theory"
 *   Bainov & Simeonov (1989) "Systems with Impulse Effect"
 */
#include "impulsive_jump.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- Linear Jump ---- */

ImpJumpLinear* imp_jump_linear_create(const double *J, int n)
{
    if (!J || n < 1 || n > IMP_MAX_DIM) return NULL;
    ImpJumpLinear *jl = (ImpJumpLinear*)calloc(1, sizeof(ImpJumpLinear));
    if (!jl) return NULL;
    jl->J = (double*)malloc((size_t)n * n * sizeof(double));
    if (!jl->J) { free(jl); return NULL; }
    for (int i = 0; i < n * n; i++) jl->J[i] = J[i];
    jl->n = n;
    return jl;
}

void imp_jump_linear_free(ImpJumpLinear *jl)
{
    if (!jl) return;
    free(jl->J); free(jl);
}

int imp_jump_linear_eval(double t_k, const double *x_before, int n,
                          double *x_after, void *ctx)
{
    (void)t_k;
    ImpJumpLinear *jl = (ImpJumpLinear*)ctx;
    if (!jl || !x_before || !x_after || n != jl->n) return -1;
    /* x^+ = (I+J) * x^- + d = x^- + J*x^- */
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
            sum += jl->J[i * n + j] * x_before[j];
        x_after[i] = x_before[i] + sum;
    }
    return 0;
}

/* ---- Affine Jump ---- */

ImpJumpAffine* imp_jump_affine_create(const double *J, const double *d, int n)
{
    if (!J || !d || n < 1 || n > IMP_MAX_DIM) return NULL;
    ImpJumpAffine *ja = (ImpJumpAffine*)calloc(1, sizeof(ImpJumpAffine));
    if (!ja) return NULL;
    ja->J = (double*)malloc((size_t)n * n * sizeof(double));
    ja->d = (double*)malloc((size_t)n * sizeof(double));
    if (!ja->J || !ja->d) { imp_jump_affine_free(ja); return NULL; }
    memcpy(ja->J, J, (size_t)n * n * sizeof(double));
    memcpy(ja->d, d, (size_t)n * sizeof(double));
    ja->n = n;
    return ja;
}

void imp_jump_affine_free(ImpJumpAffine *ja)
{
    if (!ja) return;
    free(ja->J); free(ja->d); free(ja);
}

int imp_jump_affine_eval(double t_k, const double *x_before, int n,
                          double *x_after, void *ctx)
{
    (void)t_k;
    ImpJumpAffine *ja = (ImpJumpAffine*)ctx;
    if (!ja || !x_before || !x_after || n != ja->n) return -1;
    for (int i = 0; i < n; i++) {
        double sum = ja->d[i];
        for (int j = 0; j < n; j++)
            sum += ja->J[i * n + j] * x_before[j];
        x_after[i] = x_before[i] + sum;
    }
    return 0;
}

/* ---- Projection Jump ---- */

ImpJumpProject* imp_jump_project_create(const double *lb, const double *ub,
                                          int n)
{
    if (!lb || !ub || n < 1 || n > IMP_MAX_DIM) return NULL;
    /* Validate lb[i] <= ub[i] */
    for (int i = 0; i < n; i++)
        if (lb[i] > ub[i]) return NULL;
    ImpJumpProject *jp = (ImpJumpProject*)calloc(1, sizeof(ImpJumpProject));
    if (!jp) return NULL;
    jp->lb = (double*)malloc((size_t)n * sizeof(double));
    jp->ub = (double*)malloc((size_t)n * sizeof(double));
    if (!jp->lb || !jp->ub) { imp_jump_project_free(jp); return NULL; }
    memcpy(jp->lb, lb, (size_t)n * sizeof(double));
    memcpy(jp->ub, ub, (size_t)n * sizeof(double));
    jp->n = n;
    return jp;
}

void imp_jump_project_free(ImpJumpProject *jp)
{
    if (!jp) return;
    free(jp->lb); free(jp->ub); free(jp);
}

int imp_jump_project_eval(double t_k, const double *x_before, int n,
                           double *x_after, void *ctx)
{
    (void)t_k;
    ImpJumpProject *jp = (ImpJumpProject*)ctx;
    if (!jp || !x_before || !x_after || n != jp->n) return -1;
    for (int i = 0; i < n; i++) {
        if (x_before[i] < jp->lb[i])
            x_after[i] = jp->lb[i];
        else if (x_before[i] > jp->ub[i])
            x_after[i] = jp->ub[i];
        else
            x_after[i] = x_before[i];
    }
    return 0;
}

/* ---- Hard Reset Jump ---- */

ImpJumpHardReset* imp_jump_hard_reset_create(const double *x_target, int n)
{
    if (!x_target || n < 1 || n > IMP_MAX_DIM) return NULL;
    ImpJumpHardReset *hr = (ImpJumpHardReset*)calloc(1, sizeof(ImpJumpHardReset));
    if (!hr) return NULL;
    hr->x_target = (double*)malloc((size_t)n * sizeof(double));
    if (!hr->x_target) { free(hr); return NULL; }
    memcpy(hr->x_target, x_target, (size_t)n * sizeof(double));
    hr->n = n;
    return hr;
}

void imp_jump_hard_reset_free(ImpJumpHardReset *hr)
{
    if (!hr) return;
    free(hr->x_target); free(hr);
}

int imp_jump_hard_reset_eval(double t_k, const double *x_before, int n,
                              double *x_after, void *ctx)
{
    (void)t_k; (void)x_before;
    ImpJumpHardReset *hr = (ImpJumpHardReset*)ctx;
    if (!hr || !x_after || n != hr->n) return -1;
    memcpy(x_after, hr->x_target, (size_t)n * sizeof(double));
    return 0;
}

/* ---- Impulsive Control Jump ---- */

ImpJumpControl* imp_jump_control_create(const double *B, int n, int m)
{
    if (!B || n < 1 || n > IMP_MAX_DIM || m < 1 || m > IMP_MAX_DIM) return NULL;
    ImpJumpControl *jc = (ImpJumpControl*)calloc(1, sizeof(ImpJumpControl));
    if (!jc) return NULL;
    jc->B = (double*)malloc((size_t)n * m * sizeof(double));
    jc->u = (double*)calloc((size_t)m, sizeof(double));
    if (!jc->B || !jc->u) { imp_jump_control_free(jc); return NULL; }
    memcpy(jc->B, B, (size_t)n * m * sizeof(double));
    jc->n = n; jc->m = m;
    return jc;
}

void imp_jump_control_free(ImpJumpControl *jc)
{
    if (!jc) return;
    free(jc->B); free(jc->u); free(jc);
}

void imp_jump_control_set_u(ImpJumpControl *jc, const double *u_new)
{
    if (!jc || !u_new) return;
    memcpy(jc->u, u_new, (size_t)jc->m * sizeof(double));
}

int imp_jump_control_eval(double t_k, const double *x_before, int n,
                           double *x_after, void *ctx)
{
    (void)t_k;
    ImpJumpControl *jc = (ImpJumpControl*)ctx;
    if (!jc || !x_before || !x_after || n != jc->n) return -1;
    /* x^+ = x^- + B * u */
    for (int i = 0; i < n; i++) {
        double Bu = 0.0;
        for (int j = 0; j < jc->m; j++)
            Bu += jc->B[i * jc->m + j] * jc->u[j];
        x_after[i] = x_before[i] + Bu;
    }
    return 0;
}

/* ---- Nonlinear Jump ---- */

ImpJumpNonlinear* imp_jump_nonlinear_create(ImpJumpMap g, void *params, int n)
{
    if (!g || n < 1 || n > IMP_MAX_DIM) return NULL;
    ImpJumpNonlinear *jn = (ImpJumpNonlinear*)calloc(1, sizeof(ImpJumpNonlinear));
    if (!jn) return NULL;
    jn->g = g; jn->params = params; jn->n = n;
    return jn;
}

void imp_jump_nonlinear_free(ImpJumpNonlinear *jn) { free(jn); }

int imp_jump_nonlinear_eval(double t_k, const double *x_before, int n,
                             double *x_after, void *ctx)
{
    (void)ctx;
    ImpJumpNonlinear *jn = (ImpJumpNonlinear*)ctx;
    if (!jn || !x_before || !x_after || n != jn->n) return -1;
    return jn->g(t_k, x_before, n, x_after, jn->params);
}

/* ---- Time-Varying Jump ---- */

ImpJumpTimeVarying* imp_jump_time_varying_create(int n)
{
    if (n < 1 || n > IMP_MAX_DIM) return NULL;
    ImpJumpTimeVarying *jtv = (ImpJumpTimeVarying*)calloc(1, sizeof(ImpJumpTimeVarying));
    if (!jtv) return NULL;
    jtv->J_k = (double*)calloc((size_t)n * n, sizeof(double));
    if (!jtv->J_k) { free(jtv); return NULL; }
    jtv->n = n;
    return jtv;
}

void imp_jump_time_varying_free(ImpJumpTimeVarying *jtv)
{
    if (!jtv) return;
    free(jtv->J_k); free(jtv);
}

void imp_jump_time_varying_set_gain(ImpJumpTimeVarying *jtv,
                                     const double *J_k)
{
    if (!jtv || !J_k) return;
    memcpy(jtv->J_k, J_k, (size_t)jtv->n * jtv->n * sizeof(double));
}

int imp_jump_time_varying_eval(double t_k, const double *x_before, int n,
                                double *x_after, void *ctx)
{
    (void)t_k;
    ImpJumpTimeVarying *jtv = (ImpJumpTimeVarying*)ctx;
    if (!jtv || !x_before || !x_after || n != jtv->n) return -1;
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
            sum += jtv->J_k[i * n + j] * x_before[j];
        x_after[i] = x_before[i] + sum;
    }
    return 0;
}

/* ---- Jump Map Composition ---- */

int imp_jump_compose(ImpJumpMap J1, void *ctx1,
                      ImpJumpMap J2, void *ctx2,
                      double t_k, const double *x_before, int n,
                      double *x_after)
{
    if (!J1 || !J2 || !x_before || !x_after || n < 1) return -1;
    double *x_mid = (double*)malloc((size_t)n * sizeof(double));
    if (!x_mid) return -2;
    int ret = J1(t_k, x_before, n, x_mid, ctx1);
    if (ret != 0) { free(x_mid); return ret; }
    ret = J2(t_k, x_mid, n, x_after, ctx2);
    free(x_mid);
    return ret;
}

int imp_jump_sequence(ImpJumpMap *J, void **ctx_arr, int num,
                       double t_k, const double *x_before, int n,
                       double *x_after)
{
    if (!J || !ctx_arr || num < 1 || !x_before || !x_after || n < 1)
        return -1;
    double *buf1 = (double*)malloc((size_t)n * sizeof(double));
    double *buf2 = (double*)malloc((size_t)n * sizeof(double));
    if (!buf1 || !buf2) { free(buf1); free(buf2); return -2; }

    double *src = (double*)x_before;
    double *dst = buf1;
    for (int k = 0; k < num; k++) {
        int ret = J[k](t_k, src, n, dst, ctx_arr[k]);
        if (ret != 0) { free(buf1); free(buf2); return ret; }
        /* Swap buffers for next iteration */
        double *tmp = src; src = dst; dst = tmp;
    }
    memcpy(x_after, src, (size_t)n * sizeof(double));
    free(buf1); free(buf2);
    return 0;
}
