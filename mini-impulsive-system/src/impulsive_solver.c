/*
 * impulsive_solver.c -- Numerical solvers for impulsive ODEs
 *
 * Implements: Euler, Heun, Midpoint, RK4, and adaptive RK45 methods
 * for solving impulsive differential equations.
 *
 * Algorithm: Between impulse times, integrate dx/dt = f(t,x) using
 * chosen method. At impulse times, apply jump map x^+ = I(x^-).
 * For state-driven impulses, use zero-crossing event detection.
 *
 * Error control via Richardson extrapolation (RK4 + RK5 embedded pair).
 *
 * References:
 *   Hairer, Norsett, Wanner (1993) "Solving ODE I"
 *   Shampine & Reichelt (1997) "The MATLAB ODE Suite"
 *   Dormand & Prince (1980) "A family of embedded Runge-Kutta formulae"
 */
#include "impulsive_solver.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- Default Config ---- */

ImpSolverConfig imp_solver_config_default(void)
{
    ImpSolverConfig cfg;
    cfg.method   = IMP_SOLVER_RK4;
    cfg.h_init   = 0.01;
    cfg.h_min    = 1e-8;
    cfg.h_max    = 0.1;
    cfg.atol     = 1e-8;
    cfg.rtol     = 1e-6;
    cfg.max_steps = IMP_MAX_ITER;
    cfg.detect_events = false;
    cfg.event_tol = 1e-4;
    return cfg;
}

/* ---- Euler Step ---- */

int imp_solver_euler_step(ImpVectorField f, void *ctx,
                           double t, const double *x, int n,
                           double h, double *x_next)
{
    if (!f || !x || !x_next || n < 1 || n > IMP_MAX_DIM) return -1;
    double *dxdt = (double*)malloc((size_t)n * sizeof(double));
    if (!dxdt) return -2;
    int ret = f(t, x, n, dxdt, ctx);
    if (ret != 0) { free(dxdt); return ret; }
    for (int i = 0; i < n; i++)
        x_next[i] = x[i] + h * dxdt[i];
    free(dxdt);
    return 0;
}

/* ---- Heun (Improved Euler) Step ---- */

int imp_solver_heun_step(ImpVectorField f, void *ctx,
                          double t, const double *x, int n,
                          double h, double *x_next)
{
    if (!f || !x || !x_next || n < 1) return -1;
    double *k1 = (double*)malloc((size_t)n * sizeof(double));
    double *k2 = (double*)malloc((size_t)n * sizeof(double));
    double *xt = (double*)malloc((size_t)n * sizeof(double));
    if (!k1 || !k2 || !xt) { free(k1); free(k2); free(xt); return -2; }

    f(t, x, n, k1, ctx);
    for (int i = 0; i < n; i++) xt[i] = x[i] + h * k1[i];
    f(t + h, xt, n, k2, ctx);
    for (int i = 0; i < n; i++)
        x_next[i] = x[i] + 0.5 * h * (k1[i] + k2[i]);

    free(k1); free(k2); free(xt);
    return 0;
}

/* ---- Midpoint Step ---- */

int imp_solver_midpoint_step(ImpVectorField f, void *ctx,
                              double t, const double *x, int n,
                              double h, double *x_next)
{
    if (!f || !x || !x_next || n < 1) return -1;
    double *k1 = (double*)malloc((size_t)n * sizeof(double));
    double *k2 = (double*)malloc((size_t)n * sizeof(double));
    double *xt = (double*)malloc((size_t)n * sizeof(double));
    if (!k1 || !k2 || !xt) { free(k1); free(k2); free(xt); return -2; }

    f(t, x, n, k1, ctx);
    for (int i = 0; i < n; i++) xt[i] = x[i] + 0.5 * h * k1[i];
    f(t + 0.5 * h, xt, n, k2, ctx);
    for (int i = 0; i < n; i++)
        x_next[i] = x[i] + h * k2[i];

    free(k1); free(k2); free(xt);
    return 0;
}

/* ---- RK4 Step ---- */

int imp_solver_rk4_step(ImpVectorField f, void *ctx,
                          double t, const double *x, int n,
                          double h, double *x_next)
{
    if (!f || !x || !x_next || n < 1 || n > IMP_MAX_DIM) return -1;
    double *k1 = (double*)malloc((size_t)4 * n * sizeof(double));
    double *k2 = k1 + n, *k3 = k1 + 2*n, *k4 = k1 + 3*n;
    double *xt = (double*)malloc((size_t)n * sizeof(double));
    if (!k1 || !xt) { free(k1); free(xt); return -2; }

    double h2 = 0.5 * h;

    f(t, x, n, k1, ctx);
    for (int i = 0; i < n; i++) xt[i] = x[i] + h2 * k1[i];
    f(t + h2, xt, n, k2, ctx);
    for (int i = 0; i < n; i++) xt[i] = x[i] + h2 * k2[i];
    f(t + h2, xt, n, k3, ctx);
    for (int i = 0; i < n; i++) xt[i] = x[i] + h * k3[i];
    f(t + h, xt, n, k4, ctx);

    for (int i = 0; i < n; i++)
        x_next[i] = x[i] + (h / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);

    free(k1); free(xt);
    return 0;
}

/* ---- Error Estimation ---- */

double imp_solver_estimate_error(const double *x1, const double *x2,
                                  int n, double atol, double rtol)
{
    if (!x1 || !x2 || n < 1) return IMP_HUGE;
    double err = 0.0;
    for (int i = 0; i < n; i++) {
        double scale = atol + rtol * fmax(fabs(x1[i]), fabs(x2[i]));
        double diff = fabs(x1[i] - x2[i]) / scale;
        err += diff * diff;
    }
    return sqrt(err / (double)n);
}

/* ---- Adaptive Step (RK45 embedded) ---- */

int imp_solver_adaptive_step(ImpVectorField f, void *ctx,
                              double t, const double *x, int n,
                              double *h, double *x_next,
                              double *err_out,
                              const ImpSolverConfig *cfg)
{
    if (!f || !x || !h || !x_next || !cfg || n < 1) return -1;

    /* Dormand-Prince 5(4) coefficients */
    static const double a21 = 1.0/5.0;
    static const double a31 = 3.0/40.0,  a32 = 9.0/40.0;
    static const double a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0;
    static const double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0;
    static const double a53 = 64448.0/6561.0, a54 = -212.0/729.0;
    static const double a61 = 9017.0/3168.0, a62 = -355.0/33.0;
    static const double a63 = 46732.0/5247.0, a64 = 49.0/176.0;
    static const double a65 = -5103.0/18656.0;
    static const double b1  = 35.0/384.0; /* b2=0 unused */
    static const double b3  = 500.0/1113.0, b4 = 125.0/192.0;
    static const double b5  = -2187.0/6784.0, b6 = 11.0/84.0;
    static const double d1  = 5179.0/57600.0; /* d2=0 unused */
    static const double d3  = 7571.0/16695.0, d4 = 393.0/640.0;
    static const double d5  = -92097.0/339200.0, d6 = 187.0/2100.0;
    static const double d7  = 1.0/40.0;

    int ndim = n;
    double *k1 = (double*)malloc((size_t)7 * ndim * sizeof(double));
    double *k2 = k1 + ndim, *k3 = k2 + ndim, *k4 = k3 + ndim;
    double *k5 = k4 + ndim, *k6 = k5 + ndim, *k7 = k6 + ndim;
    double *yt = (double*)malloc((size_t)ndim * sizeof(double));
    if (!k1 || !yt) { free(k1); free(yt); return -2; }

    double hh = *h;

    f(t, x, ndim, k1, ctx);
    for (int i = 0; i < ndim; i++) yt[i] = x[i] + hh * a21 * k1[i];
    f(t + 0.2*hh, yt, ndim, k2, ctx);
    for (int i = 0; i < ndim; i++) yt[i] = x[i] + hh*(a31*k1[i] + a32*k2[i]);
    f(t + 0.3*hh, yt, ndim, k3, ctx);
    for (int i = 0; i < ndim; i++) yt[i] = x[i] + hh*(a41*k1[i] + a42*k2[i] + a43*k3[i]);
    f(t + 0.8*hh, yt, ndim, k4, ctx);
    for (int i = 0; i < ndim; i++) yt[i] = x[i] + hh*(a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
    f(t + 8.0/9.0*hh, yt, ndim, k5, ctx);
    for (int i = 0; i < ndim; i++) yt[i] = x[i] + hh*(a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
    f(t + hh, yt, ndim, k6, ctx);

    /* 5th order solution */
    for (int i = 0; i < ndim; i++)
        x_next[i] = x[i] + hh*(b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i]);

    /* 4th order solution for error estimation */
    f(t + hh, x_next, ndim, k7, ctx);
    double err = 0.0;
    for (int i = 0; i < ndim; i++) {
        double e = hh*(d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*k7[i]);
        double scale = cfg->atol + cfg->rtol * fmax(fabs(x[i]), fabs(x_next[i]));
        err += (e / scale) * (e / scale);
    }
    err = sqrt(err / (double)ndim);
    if (err_out) *err_out = err;

    /* Step size control */
    double safety = 0.9;
    if (err > 1.0) {
        *h = fmax(cfg->h_min, safety * hh * pow(err, -0.2));
        free(k1); free(yt); return 1;  /* reject step */
    }
    *h = fmin(cfg->h_max, safety * hh * pow(err, -0.25));
    free(k1); free(yt);
    return 0;
}

/* ---- Apply Impulse ---- */

int imp_solver_apply_impulse(ImpJumpMap I, void *ctx,
                              double t_k, double *x, int n)
{
    if (!I || !x || n < 1) return -1;
    double *x_new = (double*)malloc((size_t)n * sizeof(double));
    if (!x_new) return -2;
    int ret = I(t_k, x, n, x_new, ctx);
    if (ret == 0) memcpy(x, x_new, (size_t)n * sizeof(double));
    free(x_new);
    return ret;
}

/* ---- Dense Output (Linear Interpolation) ---- */

int imp_solver_interpolate(const ImpSolution *sol, double t,
                            double *x_out, int n)
{
    if (!sol || !x_out || n != sol->n || sol->npts < 2) return -1;
    if (t < sol->t[0] || t > sol->T) return -1;

    /* Find interval */
    int idx = 0;
    for (int i = 0; i < sol->npts - 1; i++) {
        if (sol->t[i] <= t && t <= sol->t[i+1]) { idx = i; break; }
    }

    double dt = sol->t[idx+1] - sol->t[idx];
    if (dt < IMP_EPS) {
        memcpy(x_out, &sol->x[(size_t)idx * n], (size_t)n * sizeof(double));
        return 0;
    }
    double theta = (t - sol->t[idx]) / dt;
    const double *x0 = &sol->x[(size_t)idx * n];
    const double *x1 = &sol->x[(size_t)(idx+1) * n];
    for (int i = 0; i < n; i++)
        x_out[i] = (1.0 - theta) * x0[i] + theta * x1[i];
    return 0;
}

/* ---- Fixed-Step Simulation with Known Impulse Times ---- */

int imp_solver_simulate_fixed(const ImpSystem *sys,
                               const double *x0, int n,
                               const ImpSolverConfig *cfg,
                               ImpSolution *sol)
{
    if (!sys || !x0 || !cfg || !sol || n != sys->n) return -1;

    double t = sys->seq->t0;
    double h = cfg->h_init;
    double *x = (double*)malloc((size_t)n * sizeof(double));
    double *x_tmp = (double*)malloc((size_t)n * sizeof(double));
    if (!x || !x_tmp) { free(x); free(x_tmp); return -2; }
    memcpy(x, x0, (size_t)n * sizeof(double));

    sol->t0 = t;
    imp_solution_add_point(sol, t, x);

    int next_imp_idx = 0;
    int (*step_fn)(ImpVectorField, void*, double, const double*, int, double, double*);
    switch (cfg->method) {
        case IMP_SOLVER_EULER:    step_fn = imp_solver_euler_step; break;
        case IMP_SOLVER_HEUN:     step_fn = imp_solver_heun_step; break;
        case IMP_SOLVER_MIDPOINT: step_fn = imp_solver_midpoint_step; break;
        default:                  step_fn = imp_solver_rk4_step; break;
    }

    int nsteps = 0;
    while (t < sys->seq->T && nsteps < cfg->max_steps) {
        /* Check if next step crosses an impulse time */
        double t_next = t + h;
        bool impulse_due = false;
        double tau_k = t_next;

        while (next_imp_idx < sys->seq->count &&
               sys->seq->times[next_imp_idx] <= t + IMP_EPS)
            next_imp_idx++;

        if (next_imp_idx < sys->seq->count) {
            tau_k = sys->seq->times[next_imp_idx];
            if (t_next >= tau_k - IMP_EPS) {
                impulse_due = true;
                h = tau_k - t;
            }
        }

        /* Integrate one step */
        if (h < cfg->h_min) h = cfg->h_min;
        int ret = step_fn(sys->f, sys->ctx, t, x, n, h, x_tmp);
        if (ret != 0) { free(x); free(x_tmp); return ret; }

        t += h;
        memcpy(x, x_tmp, (size_t)n * sizeof(double));

        /* Apply impulse if needed */
        if (impulse_due) {
            double *x_before = (double*)malloc((size_t)n * sizeof(double));
            memcpy(x_before, x, (size_t)n * sizeof(double));
            imp_solver_apply_impulse(sys->I, sys->ctx, tau_k, x, n);
            imp_solution_add_jump(sol, x_before, x);
            free(x_before);
            next_imp_idx++;
            h = cfg->h_init;  /* reset step size after impulse */
        }

        imp_solution_add_point(sol, t, x);
        nsteps++;
    }

    free(x); free(x_tmp);
    return 0;
}

/* ---- Adaptive + Event Detection Simulation ---- */

int imp_solver_simulate_adaptive(const ImpSystem *sys,
                                  const double *x0, int n,
                                  const ImpSolverConfig *cfg,
                                  ImpSolution *sol)
{
    if (!sys || !x0 || !cfg || !sol || n != sys->n) return -1;

    double t = sys->seq->t0;
    double h = cfg->h_init;
    double *x = (double*)malloc((size_t)n * sizeof(double));
    double *x_tmp = (double*)malloc((size_t)n * sizeof(double));
    if (!x || !x_tmp) { free(x); free(x_tmp); return -2; }
    memcpy(x, x0, (size_t)n * sizeof(double));

    sol->t0 = t;
    imp_solution_add_point(sol, t, x);

    int next_imp_idx = 0;
    int nsteps = 0;
    while (t < sys->seq->T && nsteps < cfg->max_steps) {
        /* Clamp step to next impulse */
        while (next_imp_idx < sys->seq->count &&
               sys->seq->times[next_imp_idx] <= t + IMP_EPS)
            next_imp_idx++;

        if (next_imp_idx < sys->seq->count) {
            double tau = sys->seq->times[next_imp_idx];
            if (t + h > tau) h = tau - t;
            if (h < cfg->h_min) h = cfg->h_min;
        }

        double err;
        int ret = imp_solver_adaptive_step(sys->f, sys->ctx, t, x, n,
                                            &h, x_tmp, &err, cfg);
        if (ret == 1) continue;  /* step rejected, retry */
        if (ret != 0) break;

        t += h;
        memcpy(x, x_tmp, (size_t)n * sizeof(double));

        /* Check for impulse */
        if (next_imp_idx < sys->seq->count &&
            fabs(t - sys->seq->times[next_imp_idx]) < cfg->event_tol) {
            double *x_before = (double*)malloc((size_t)n * sizeof(double));
            memcpy(x_before, x, (size_t)n * sizeof(double));
            imp_solver_apply_impulse(sys->I, sys->ctx,
                                      sys->seq->times[next_imp_idx], x, n);
            imp_solution_add_jump(sol, x_before, x);
            free(x_before);
            next_imp_idx++;
            h = cfg->h_init;
        }

        imp_solution_add_point(sol, t, x);
        nsteps++;
    }

    free(x); free(x_tmp);
    return 0;
}

/* ---- Main Simulate (dispatches to fixed or adaptive) ---- */

int imp_solver_simulate(const ImpSystem *sys,
                         const double *x0, int n,
                         double t0, double T,
                         const ImpSolverConfig *cfg,
                         ImpSolution *sol)
{
    (void)t0; (void)T;
    if (!sys || !x0 || !cfg || !sol || n != sys->n) return -1;
    if (!imp_system_validate(sys)) return -1;

    if (cfg->method == IMP_SOLVER_RK45)
        return imp_solver_simulate_adaptive(sys, x0, n, cfg, sol);
    else
        return imp_solver_simulate_fixed(sys, x0, n, cfg, sol);
}
