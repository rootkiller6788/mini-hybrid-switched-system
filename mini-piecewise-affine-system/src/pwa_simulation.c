/**
 * @file pwa_simulation.c
 * @brief PWA System Simulation Engine — L5 Algorithms
 *
 * Simulates piecewise affine systems with event detection,
 * numerical integration, Zeno detection, and trajectory
 * recording and export.
 *
 * Knowledge coverage:
 *   L5: RK4 integration, adaptive RKF45, event detection via
 *       bisection, Zeno detection and regularization,
 *       trajectory simulation, CSV export
 *
 * References:
 *   Hairer, Nørsett, Wanner (1993). "Solving ODEs I." Springer.
 *   Cellier & Kofman (2006). "Continuous System Simulation." Springer.
 *   Zhang et al. (2001). "Zeno hybrid systems." IJRNLC, 11(5):435-451.
 */

#include "pwa_simulation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

/*===========================================================================
 * L5: Default Configuration
 *===========================================================================*/

PWASimConfig pwa_sim_config_default(void)
{
    PWASimConfig cfg;
    cfg.t_start = 0.0;
    cfg.t_end = 10.0;
    cfg.dt_max = 0.01;
    cfg.dt_min = 1e-8;
    cfg.max_steps = 100000;
    cfg.event_tol = 1e-8;
    cfg.max_events = 1000;
    cfg.store_all = 1;
    cfg.verbose = 0;
    return cfg;
}

/*===========================================================================
 * L5: Vector Field Evaluation
 *===========================================================================*/

int pwa_vector_field(const PWASystem *sys, int region,
                      const double *x, const double *u, double *dxdt)
{
    if (!sys || !x || !dxdt) return -1;
    if (region < 0 || region >= sys->n_regions) return -1;

    const PWAAffineDynamics *dyn = &sys->dynamics[region];
    int n = sys->n_state;
    int m = sys->n_input;

    /* dx/dt = A x + B u + f */
    for (int i = 0; i < n; i++) {
        dxdt[i] = dyn->f ? dyn->f[i] : 0.0;

        if (dyn->A) {
            for (int j = 0; j < n; j++) {
                dxdt[i] += dyn->A[i * n + j] * x[j];
            }
        }

        if (dyn->B && u) {
            for (int j = 0; j < m; j++) {
                dxdt[i] += dyn->B[i * m + j] * u[j];
            }
        }
    }

    return 0;
}

/*===========================================================================
 * L5: RK4 Integration Step
 *===========================================================================*/

int pwa_rk4_step(const PWASystem *sys, int region,
                  double t, const double *x, const double *u,
                  double h, double *x_next)
{
    (void)t;
    if (!sys || !x || !x_next) return -1;

    int n = sys->n_state;
    double *k1 = (double*)malloc((size_t)(4 * n) * sizeof(double));
    if (!k1) return -1;
    double *k2 = k1 + n;
    double *k3 = k1 + 2 * n;
    double *k4 = k1 + 3 * n;
    double *xtmp = (double*)malloc((size_t)n * sizeof(double));
    if (!xtmp) { free(k1); return -1; }

    /* k1 = f(t, x) */
    pwa_vector_field(sys, region, x, u, k1);

    /* k2 = f(t + h/2, x + h*k1/2) */
    for (int i = 0; i < n; i++) xtmp[i] = x[i] + 0.5 * h * k1[i];
    pwa_vector_field(sys, region, xtmp, u, k2);

    /* k3 = f(t + h/2, x + h*k2/2) */
    for (int i = 0; i < n; i++) xtmp[i] = x[i] + 0.5 * h * k2[i];
    pwa_vector_field(sys, region, xtmp, u, k3);

    /* k4 = f(t + h, x + h*k3) */
    for (int i = 0; i < n; i++) xtmp[i] = x[i] + h * k3[i];
    pwa_vector_field(sys, region, xtmp, u, k4);

    /* x_next = x + h*(k1 + 2*k2 + 2*k3 + k4)/6 */
    for (int i = 0; i < n; i++) {
        x_next[i] = x[i] + (h / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

    free(k1);
    free(xtmp);
    return 0;
}

/*===========================================================================
 * L5: RKF45 Adaptive Step
 *===========================================================================*/

int pwa_rkf45_step(const PWASystem *sys, int region,
                    double t, const double *x, const double *u,
                    double *h, double *x_next, double *err, double tol)
{
    (void)t;
    if (!sys || !x || !u || !h || !x_next || !err) return -1;

    int n = sys->n_state;
    double *k1 = (double*)malloc((size_t)(6 * n) * sizeof(double));
    if (!k1) return -1;
    double *k2 = k1 + n, *k3 = k1 + 2*n, *k4 = k1 + 3*n;
    double *k5 = k1 + 4*n, *k6 = k1 + 5*n;
    double *xtmp = (double*)malloc((size_t)n * sizeof(double));
    if (!xtmp) { free(k1); return -1; }

    double hi = *h;

    /* Butcher tableau for RKF45 */
    /* k1 = f(t, x) */
    pwa_vector_field(sys, region, x, u, k1);

    /* k2 = f(t + 1/4 h, x + 1/4 h k1) */
    for (int i = 0; i < n; i++) xtmp[i] = x[i] + 0.25 * hi * k1[i];
    pwa_vector_field(sys, region, xtmp, u, k2);

    /* k3 = f(t + 3/8 h, x + 3/32 h k1 + 9/32 h k2) */
    for (int i = 0; i < n; i++)
        xtmp[i] = x[i] + hi * (3.0/32.0 * k1[i] + 9.0/32.0 * k2[i]);
    pwa_vector_field(sys, region, xtmp, u, k3);

    /* k4 = f(t + 12/13 h, x + 1932/2197 h k1 - 7200/2197 h k2
     *                       + 7296/2197 h k3) */
    for (int i = 0; i < n; i++)
        xtmp[i] = x[i] + hi * (1932.0/2197.0 * k1[i]
                              - 7200.0/2197.0 * k2[i]
                              + 7296.0/2197.0 * k3[i]);
    pwa_vector_field(sys, region, xtmp, u, k4);

    /* k5 = f(t + h, x + 439/216 h k1 - 8 h k2 + 3680/513 h k3
     *                    - 845/4104 h k4) */
    for (int i = 0; i < n; i++)
        xtmp[i] = x[i] + hi * (439.0/216.0 * k1[i] - 8.0 * k2[i]
                              + 3680.0/513.0 * k3[i]
                              - 845.0/4104.0 * k4[i]);
    pwa_vector_field(sys, region, xtmp, u, k5);

    /* k6 = f(t + 1/2 h, x - 8/27 h k1 + 2 h k2 - 3544/2565 h k3
     *                         + 1859/4104 h k4 - 11/40 h k5) */
    for (int i = 0; i < n; i++)
        xtmp[i] = x[i] + hi * (-8.0/27.0 * k1[i] + 2.0 * k2[i]
                              - 3544.0/2565.0 * k3[i]
                              + 1859.0/4104.0 * k4[i]
                              - 0.275 * k5[i]);
    pwa_vector_field(sys, region, xtmp, u, k6);

    /* 4th order estimate */
    for (int i = 0; i < n; i++)
        x_next[i] = x[i] + hi * (25.0/216.0 * k1[i]
                               + 1408.0/2565.0 * k3[i]
                               + 2197.0/4104.0 * k4[i]
                               - 0.2 * k5[i]);

    /* 5th order estimate for error */
    double max_err = 0.0;
    for (int i = 0; i < n; i++) {
        double e5 = hi * (16.0/135.0 * k1[i]
                        + 6656.0/12825.0 * k3[i]
                        + 28561.0/56430.0 * k4[i]
                        - 0.18 * k5[i]
                        + 2.0/55.0 * k6[i]);
        double e4 = x_next[i] - x[i];
        double e = fabs(e5 - e4);
        err[i] = e;
        if (e > max_err) max_err = e;
    }

    /* Adjust step size */
    if (max_err > tol) {
        double factor = 0.9 * pow(tol / max_err, 0.2);
        if (factor < 0.1) factor = 0.1;
        *h *= factor;
    } else if (max_err < tol * 0.1) {
        double factor = 0.9 * pow(tol / max_err, 0.25);
        if (factor > 5.0) factor = 5.0;
        *h *= factor;
    }

    free(k1);
    free(xtmp);
    return 0;
}

/*===========================================================================
 * L5: Boundary Crossing Detection
 *===========================================================================*/

int pwa_detect_boundary_crossing(const PWASystem *sys, int region,
                                  const double *x0, const double *u0,
                                  const double *x1, const double *u1,
                                  int *crossed_halfspace)
{
    if (!sys || !x0 || !x1) return 0;
    if (region < 0 || region >= sys->n_regions) return 0;

    const PWARegion *reg = &sys->regions[region];
    int nz = sys->n_state + sys->n_input;

    /* Build extended vectors */
    double *z0 = (double*)malloc((size_t)nz * sizeof(double));
    double *z1 = (double*)malloc((size_t)nz * sizeof(double));
    if (!z0 || !z1) { free(z0); free(z1); return 0; }

    memcpy(z0, x0, (size_t)sys->n_state * sizeof(double));
    memcpy(z1, x1, (size_t)sys->n_state * sizeof(double));
    if (u0 && sys->n_input > 0)
        memcpy(z0 + sys->n_state, u0, (size_t)sys->n_input * sizeof(double));
    if (u1 && sys->n_input > 0)
        memcpy(z1 + sys->n_state, u1, (size_t)sys->n_input * sizeof(double));

    *crossed_halfspace = -1;
    int crossing_detected = 0;

    for (int c = 0; c < reg->n_constraints; c++) {
        double g0 = -reg->K[c];
        double g1 = -reg->K[c];
        const double *Hc = &reg->H[c * nz];

        for (int j = 0; j < nz; j++) {
            g0 += Hc[j] * z0[j];
            g1 += Hc[j] * z1[j];
        }

        /* Sign change: g0 ≤ 0 (inside), g1 > 0 (outside) */
        if (g0 <= 1e-10 && g1 > 1e-10) {
            *crossed_halfspace = c;
            crossing_detected = 1;
            break;
        }
    }

    free(z0);
    free(z1);
    return crossing_detected;
}

/*===========================================================================
 * L5: Bisection for Exact Crossing Time
 *===========================================================================*/

int pwa_find_crossing_time(const PWASystem *sys, int region,
                            double ta, const double *xa,
                            const double *u_val,
                            double tb, const double *xb,
                            int halfspace,
                            double *t_cross, double *x_cross, double tol)
{
    if (!sys || !xa || !xb || !t_cross || !x_cross) return -1;

    int n = sys->n_state;
    int nz = n + sys->n_input;
    const PWARegion *reg = &sys->regions[region];
    const double *Hc = &reg->H[halfspace * nz];
    double Kc = reg->K[halfspace];
    (void)ta; (void)tb;

    /* Constraint: g(z) = Hc·z - Kc = 0 at crossing */
    double *x_L = (double*)malloc((size_t)n * sizeof(double));
    double *x_R = (double*)malloc((size_t)n * sizeof(double));
    double *x_M = (double*)malloc((size_t)n * sizeof(double));
    double *z_M = (double*)malloc((size_t)nz * sizeof(double));

    if (!x_L || !x_R || !x_M || !z_M) {
        free(x_L); free(x_R); free(x_M); free(z_M);
        return -1;
    }

    memcpy(x_L, xa, (size_t)n * sizeof(double));
    memcpy(x_R, xb, (size_t)n * sizeof(double));

    /* Constraint function g(z) = Hc·[x;u_val] - Kc */
    /* Inline: evaluate constraint at a state */
    double g_L = -Kc;
    {
        memcpy(z_M, xa, (size_t)n * sizeof(double));
        if (u_val) memcpy(z_M + n, u_val, (size_t)sys->n_input * sizeof(double));
        else memset(z_M + n, 0, (size_t)sys->n_input * sizeof(double));
        for (int j = 0; j < nz; j++) g_L += Hc[j] * z_M[j];
    }

    double g_R = -Kc;
    {
        memcpy(z_M, xb, (size_t)n * sizeof(double));
        if (u_val) memcpy(z_M + n, u_val, (size_t)sys->n_input * sizeof(double));
        else memset(z_M + n, 0, (size_t)sys->n_input * sizeof(double));
        for (int j = 0; j < nz; j++) g_R += Hc[j] * z_M[j];
    }
    (void)ta; (void)tb;

    /* Ensure we bracket the root */
    if (g_L * g_R > 0) {
        /* No sign change within tolerance; just use midpoint */
        *t_cross = (ta + tb) * 0.5;
        for (int i = 0; i < n; i++)
            x_cross[i] = (xa[i] + xb[i]) * 0.5;
        free(x_L); free(x_R); free(x_M); free(z_M);
        return 0;
    }

    /* Bisection */
    double t_L = ta, t_R = tb;
    for (int iter = 0; iter < 60; iter++) {
        double t_M = (t_L + t_R) * 0.5;
        for (int i = 0; i < n; i++)
            x_M[i] = (x_L[i] + x_R[i]) * 0.5;

        double g_M = -Kc;
        {
            memcpy(z_M, x_M, (size_t)n * sizeof(double));
            if (u_val) memcpy(z_M + n, u_val, (size_t)sys->n_input * sizeof(double));
            else memset(z_M + n, 0, (size_t)sys->n_input * sizeof(double));
            for (int j = 0; j < nz; j++) g_M += Hc[j] * z_M[j];
        }

        if (fabs(g_M) < tol) {
            *t_cross = t_M;
            memcpy(x_cross, x_M, (size_t)n * sizeof(double));
            free(x_L); free(x_R); free(x_M); free(z_M);
            return 0;
        }

        if (g_L * g_M < 0) {
            t_R = t_M;
            memcpy(x_R, x_M, (size_t)n * sizeof(double));
            g_R = g_M;
        } else {
            t_L = t_M;
            memcpy(x_L, x_M, (size_t)n * sizeof(double));
            g_L = g_M;
        }
    }

    *t_cross = (t_L + t_R) * 0.5;
    for (int i = 0; i < n; i++)
        x_cross[i] = (x_L[i] + x_R[i]) * 0.5;

    free(x_L); free(x_R); free(x_M); free(z_M);
    return 0;
}

/*===========================================================================
 * L5: Discrete-Time PWA Simulation
 *===========================================================================*/

int pwa_simulate_dt(PWASystem *sys,
                     const double *x0, const double *u_seq,
                     const PWASimConfig *config, PWATrajectory *traj)
{
    if (!sys || !x0 || !config || !traj) return -1;
    if (!sys->is_continuous && sys->dt <= 0) return -1;

    int n = sys->n_state;
    int m = sys->n_input;
    int p = sys->n_output;

    double dt = sys->dt;
    int max_steps = config->max_steps;
    double t_end = config->t_end;

    /* Initialize trajectory with initial state */
    traj->n_steps = 0;
    traj->n_state = n;
    traj->n_input = m;
    traj->n_output = p;
    traj->n_events = 0;

    double *x_cur = (double*)malloc((size_t)n * sizeof(double));
    double *u_cur = (double*)calloc((size_t)m, sizeof(double));
    if (!x_cur || !u_cur) { free(x_cur); free(u_cur); return -1; }

    memcpy(x_cur, x0, (size_t)n * sizeof(double));

    double t = config->t_start;
    int step = 0;

    /* Store initial state */
    traj->t_hist[step] = t;
    memcpy(traj->x_hist + step * n, x_cur, (size_t)n * sizeof(double));
    if (u_seq) memcpy(traj->u_hist + step * m, u_seq + step * m,
                       (size_t)m * sizeof(double));
    traj->region_hist[step] = -1;

    step++;
    traj->n_steps = step;

    for (; step < max_steps && t < t_end; step++) {
        t += dt;

        /* Get input for this step */
        if (u_seq) {
            memcpy(u_cur, u_seq + step * m, (size_t)m * sizeof(double));
        }

        /* Find active region */
        int region = pwa_point_location(sys, x_cur, u_cur);
        if (region < 0) {
            /* Try with pure state */
            region = pwa_point_location(sys, x_cur, NULL);
        }
        if (region < 0) {
            /* State outside all regions - use last known region */
            if (traj->region_hist[step-1] >= 0) {
                region = traj->region_hist[step-1];
            } else {
                break;  /* Invalid */
            }
        }

        /* Apply affine dynamics */
        const PWAAffineDynamics *dyn = &sys->dynamics[region];
        double *x_next = (double*)calloc((size_t)n, sizeof(double));
        if (!x_next) break;

        /* x_next = A x + B u + f */
        for (int i = 0; i < n; i++) {
            x_next[i] = dyn->f ? dyn->f[i] : 0.0;
            if (dyn->A) {
                for (int j = 0; j < n; j++)
                    x_next[i] += dyn->A[i * n + j] * x_cur[j];
            }
            if (dyn->B) {
                for (int j = 0; j < m; j++)
                    x_next[i] += dyn->B[i * m + j] * u_cur[j];
            }
        }

        /* Detect mode switch */
        if (region != traj->region_hist[step-1] && traj->region_hist[step-1] >= 0) {
            if (traj->n_events < config->max_events) {
                PWAEvent *ev = &traj->events[traj->n_events];
                ev->t = t;
                ev->from_region = traj->region_hist[step-1];
                ev->to_region = region;
                ev->type = PWA_SWITCH_AUTONOMOUS;
                ev->x_at_event = (double*)malloc((size_t)n * sizeof(double));
                if (ev->x_at_event) memcpy(ev->x_at_event, x_cur,
                                             (size_t)n * sizeof(double));
                traj->n_events++;
            }
        }

        /* Store trajectory point */
        traj->t_hist[step] = t;
        memcpy(traj->x_hist + step * n, x_next, (size_t)n * sizeof(double));
        if (u_seq) memcpy(traj->u_hist + step * m, u_seq + step * m,
                           (size_t)m * sizeof(double));
        traj->region_hist[step] = region;

        /* Update state */
        memcpy(x_cur, x_next, (size_t)n * sizeof(double));
        free(x_next);

        traj->n_steps = step + 1;
    }

    free(x_cur);
    free(u_cur);
    return 0;
}

/*===========================================================================
 * L5: Continuous-Time PWA Simulation
 *===========================================================================*/

int pwa_simulate_ct(PWASystem *sys,
                     const double *x0,
                     double (*u_func)(double t, int dim, void *ctx),
                     void *u_func_ctx,
                     const PWASimConfig *config, PWATrajectory *traj)
{
    if (!sys || !x0 || !config || !traj) return -1;

    int n = sys->n_state;
    int m = sys->n_input;
    int p = sys->n_output;

    double t = config->t_start;
    double t_end = config->t_end;
    double h = config->dt_max;
    double dt_min = config->dt_min;

    /* Initialize trajectory */
    traj->n_steps = 0;
    traj->n_state = n;
    traj->n_input = m;
    traj->n_output = p;
    traj->n_events = 0;

    double *x_cur = (double*)malloc((size_t)n * sizeof(double));
    double *x_try = (double*)malloc((size_t)n * sizeof(double));
    double *u_cur = (double*)calloc((size_t)m, sizeof(double));
    double *dxdt = (double*)calloc((size_t)n, sizeof(double));
    double *y_cur = (double*)calloc((size_t)p, sizeof(double));

    if (!x_cur || !x_try || !u_cur) {
        free(x_cur); free(x_try); free(u_cur); free(dxdt); free(y_cur);
        return -1;
    }

    memcpy(x_cur, x0, (size_t)n * sizeof(double));

    /* Evaluate initial input */
    if (u_func && m > 0) {
        for (int i = 0; i < m; i++)
            u_cur[i] = u_func(t, i, u_func_ctx);
    }

    /* Find initial region */
    int cur_region = pwa_point_location(sys, x_cur, u_cur);
    if (cur_region < 0) {
        /* Force: search all regions */
        cur_region = 0;
    }

    /* Store initial state */
    traj->t_hist[0] = t;
    memcpy(traj->x_hist, x_cur, (size_t)n * sizeof(double));
    memcpy(traj->u_hist, u_cur, (size_t)m * sizeof(double));
    traj->region_hist[0] = cur_region;
    traj->n_steps = 1;

    int step = 1;
    int event_count = 0;

    while (t < t_end - 1e-10 && step < config->max_steps) {
        if (h < dt_min) {
            h = dt_min;  /* Prevent infinite looping */
        }

        double h_use = h;
        if (t + h_use > t_end) h_use = t_end - t;

        /* Take RK4 step */
        int ret = pwa_rk4_step(sys, cur_region, t, x_cur, u_cur, h_use, x_try);
        if (ret < 0) break;

        /* Check for boundary crossing */
        int new_region;
        int crossed_hs;
        int has_crossing = pwa_detect_boundary_crossing(sys, cur_region,
                                 x_cur, u_cur, x_try, u_cur, &crossed_hs);

        if (has_crossing && crossed_hs >= 0) {
            /* Find exact crossing time */
            double t_cross;
            double *x_cross = (double*)malloc((size_t)n * sizeof(double));
            if (x_cross) {
                pwa_find_crossing_time(sys, cur_region, t, x_cur, u_cur,
                                        t + h_use, x_try, crossed_hs,
                                        &t_cross, x_cross, config->event_tol);

                /* Store state at crossing */
                traj->t_hist[step] = t_cross;
                memcpy(traj->x_hist + step * n, x_cross, (size_t)n * sizeof(double));
                memcpy(traj->u_hist + step * m, u_cur, (size_t)m * sizeof(double));
                traj->region_hist[step] = cur_region;
                step++;
                traj->n_steps = step;
                t = t_cross;
                memcpy(x_cur, x_cross, (size_t)n * sizeof(double));
                free(x_cross);

                /* Find new region */
                new_region = pwa_point_location(sys, x_cur, u_cur);
                if (new_region < 0) new_region = 0;

                /* Record event */
                if (new_region != cur_region && event_count < config->max_events) {
                    PWAEvent *ev = &traj->events[event_count];
                    ev->t = t_cross;
                    ev->from_region = cur_region;
                    ev->to_region = new_region;
                    ev->type = PWA_SWITCH_AUTONOMOUS;
                    ev->x_at_event = (double*)malloc((size_t)n * sizeof(double));
                    if (ev->x_at_event)
                        memcpy(ev->x_at_event, x_cur, (size_t)n * sizeof(double));
                    event_count++;
                    traj->n_events = event_count;
                }

                cur_region = new_region;
                h = config->dt_max;  /* Reset step size */
            }
        } else {
            /* No crossing: accept step */
            t += h_use;
            memcpy(x_cur, x_try, (size_t)n * sizeof(double));

            /* Update input */
            if (u_func && m > 0) {
                for (int i = 0; i < m; i++)
                    u_cur[i] = u_func(t, i, u_func_ctx);
            }

            /* Find region (may have changed without crossing detection) */
            new_region = pwa_point_location(sys, x_cur, u_cur);
            if (new_region >= 0 && new_region != cur_region) {
                if (event_count < config->max_events) {
                    PWAEvent *ev = &traj->events[event_count];
                    ev->t = t;
                    ev->from_region = cur_region;
                    ev->to_region = new_region;
                    ev->type = PWA_SWITCH_AUTONOMOUS;
                    ev->x_at_event = (double*)malloc((size_t)n * sizeof(double));
                    if (ev->x_at_event)
                        memcpy(ev->x_at_event, x_cur, (size_t)n * sizeof(double));
                    event_count++;
                    traj->n_events = event_count;
                }
                cur_region = new_region;
            }
            if (new_region < 0) new_region = cur_region;

            /* Store */
            traj->t_hist[step] = t;
            memcpy(traj->x_hist + step * n, x_cur, (size_t)n * sizeof(double));
            memcpy(traj->u_hist + step * m, u_cur, (size_t)m * sizeof(double));
            traj->region_hist[step] = new_region >= 0 ? new_region : cur_region;
            step++;
            traj->n_steps = step;
        }
    }

    free(x_cur); free(x_try); free(u_cur); free(dxdt); free(y_cur);
    return 0;
}

/*===========================================================================
 * L5: Zeno Detection
 *===========================================================================*/

int pwa_detect_zeno(const PWATrajectory *traj, double *tb)
{
    if (!traj || !tb) return 0;

    if (traj->n_events < 10) {
        *tb = 0.0;
        return 0;
    }

    /* Check if inter-event times are decreasing and accumulating */
    double *intervals = (double*)malloc((size_t)(traj->n_events - 1) * sizeof(double));
    if (!intervals) return 0;

    for (int i = 0; i < traj->n_events - 1; i++) {
        intervals[i] = traj->events[i+1].t - traj->events[i].t;
    }

    /* Exponential fit: interval ≈ c * λ^i */
    /* Sum of intervals: convergence check */
    double total = 0.0;
    for (int i = 0; i < traj->n_events - 1; i++) {
        total += intervals[i];
    }

    /* Check if interval ratio is decreasing */
    double prev = intervals[0];
    int decreasing_count = 0;
    for (int i = 1; i < traj->n_events - 1; i++) {
        if (intervals[i] < prev * 0.95) {
            decreasing_count++;
        }
        prev = intervals[i];
    }

    free(intervals);

    if (decreasing_count > traj->n_events / 4) {
        *tb = traj->events[0].t + total;
        return 1;
    }

    *tb = 0.0;
    return 0;
}

int pwa_regularize_zeno(PWASystem *sys, double min_dwell, double hysteresis)
{
    if (!sys || min_dwell < 0.0) return -1;

    /* Add hysteresis band to region boundaries:
     * For each constraint H_i z ≤ K_i, replace with:
     *   H_i z ≤ K_i + hysteresis  (outer boundary)
     *
     * During simulation, use modified constraints that require
     * crossing an extended boundary before a switch is triggered.
     * After switching, the new region's constraints are tightened
     * to prevent immediate switching back.
     *
     * For simplicity, store a parallel set of "relaxed" constraints. */

    for (int i = 0; i < sys->n_regions; i++) {
        if (!sys->regions[i].is_active) continue;
        PWARegion *reg = &sys->regions[i];

        /* Add hysteresis: tighten all constraints temporarily
         * For now, just modify K values by hysteresis */
        double *K_hyst = (double*)malloc((size_t)reg->n_constraints * sizeof(double));
        if (!K_hyst) continue;

        memcpy(K_hyst, reg->K, (size_t)reg->n_constraints * sizeof(double));
        for (int c = 0; c < reg->n_constraints; c++) {
            K_hyst[c] -= hysteresis;  /* Tighten to prevent chatter */
        }

        /* Save original K and replace with tightened version */
        double *K_orig = reg->K;
        reg->K = K_hyst;
        free(K_orig);  /* Old constraints freed; new tighter ones active */
    }

    return 0;
}

/*===========================================================================
 * L5: Trajectory Management
 *===========================================================================*/

PWATrajectory* pwa_trajectory_create(int n_state, int n_input,
                                      int n_output, int n_max)
{
    if (n_max <= 0) return NULL;

    PWATrajectory *traj = (PWATrajectory*)calloc(1, sizeof(PWATrajectory));
    if (!traj) return NULL;

    traj->n_state = n_state;
    traj->n_input = n_input;
    traj->n_output = n_output;
    traj->n_max = n_max;
    traj->n_steps = 0;
    traj->n_events = 0;

    traj->t_hist = (double*)calloc((size_t)n_max, sizeof(double));
    traj->x_hist = (double*)calloc((size_t)(n_max * n_state), sizeof(double));
    traj->u_hist = (double*)calloc((size_t)(n_max * n_input), sizeof(double));
    traj->y_hist = (double*)calloc((size_t)(n_max * n_output), sizeof(double));
    traj->region_hist = (int*)calloc((size_t)n_max, sizeof(int));

    traj->events = (PWAEvent*)calloc((size_t)n_max, sizeof(PWAEvent));

    if (!traj->t_hist || !traj->x_hist || !traj->u_hist ||
        !traj->y_hist || !traj->region_hist || !traj->events) {
        pwa_trajectory_destroy(traj);
        return NULL;
    }

    return traj;
}

void pwa_trajectory_destroy(PWATrajectory *traj)
{
    if (!traj) return;

    for (int i = 0; i < traj->n_events; i++) {
        free(traj->events[i].x_at_event);
        free(traj->events[i].boundary_dist);
    }

    free(traj->t_hist);
    free(traj->x_hist);
    free(traj->u_hist);
    free(traj->y_hist);
    free(traj->region_hist);
    free(traj->events);
    free(traj);
}

void pwa_trajectory_print(const PWATrajectory *traj)
{
    if (!traj) {
        printf("Trajectory: NULL\n");
        return;
    }

    printf("=== PWA Trajectory ===\n");
    printf("Steps: %d (max %d)\n", traj->n_steps, traj->n_max);
    printf("Dimensions: n_state=%d, n_input=%d, n_output=%d\n",
           traj->n_state, traj->n_input, traj->n_output);
    printf("Events: %d\n", traj->n_events);

    if (traj->n_steps > 0) {
        printf("Time range: [%.4f, %.4f]\n",
               traj->t_hist[0], traj->t_hist[traj->n_steps - 1]);

        int regions_visited[100] = {0};
        int n_regions_visited = 0;
        for (int i = 0; i < traj->n_steps; i++) {
            int r = traj->region_hist[i];
            int found = 0;
            for (int j = 0; j < n_regions_visited; j++) {
                if (regions_visited[j] == r) { found = 1; break; }
            }
            if (!found && n_regions_visited < 100) {
                regions_visited[n_regions_visited++] = r;
            }
        }
        printf("Regions visited: %d (", n_regions_visited);
        for (int j = 0; j < n_regions_visited; j++) {
            printf("%d%s", regions_visited[j],
                   j < n_regions_visited - 1 ? ", " : "");
        }
        printf(")\n");
    }
}

int pwa_trajectory_export_csv(const PWATrajectory *traj, const char *filename)
{
    if (!traj || !filename) return -1;

    FILE *fp = fopen(filename, "w");
    if (!fp) return -1;

    /* Header */
    fprintf(fp, "time");
    for (int i = 0; i < traj->n_state; i++) fprintf(fp, ",x%d", i);
    for (int i = 0; i < traj->n_input; i++) fprintf(fp, ",u%d", i);
    for (int i = 0; i < traj->n_output; i++) fprintf(fp, ",y%d", i);
    fprintf(fp, ",region\n");

    /* Data rows */
    for (int step = 0; step < traj->n_steps; step++) {
        fprintf(fp, "%.10f", traj->t_hist[step]);
        for (int i = 0; i < traj->n_state; i++)
            fprintf(fp, ",%.10f", traj->x_hist[step * traj->n_state + i]);
        for (int i = 0; i < traj->n_input; i++)
            fprintf(fp, ",%.10f", traj->u_hist[step * traj->n_input + i]);
        for (int i = 0; i < traj->n_output; i++)
            fprintf(fp, ",%.10f", traj->y_hist[step * traj->n_output + i]);
        fprintf(fp, ",%d\n", traj->region_hist[step]);
    }

    fclose(fp);
    return 0;
}
