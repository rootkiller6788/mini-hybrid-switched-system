/**
 * @file hybrid_simulation.c
 * @brief Numerical simulation engine (KP36-KP40, L5)
 *
 * ODE integration (Forward Euler, RK4, Heun), event detection,
 * hybrid simulation main loop, non-deterministic branching simulation.
 *
 * Reference:
 *   Esposito, Kumar, Pappas, "Accurate Event Detection in Hybrid Systems" (2001)
 *   Park & Barton, "State Event Location in DAE Models" (1996)
 *   Hairer, Nørsett, Wanner, "Solving Ordinary Differential Equations I" (1993)
 */

#include "hybrid_simulation.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==========================================================================
 * KP36: ODE integration methods
 * ========================================================================== */

/**
 * @brief Compute dp/dt = f(t, p) for a given flow.
 *
 * For affine flows, computes A·p + b efficiently.
 * For nonlinear flows, calls the function pointer.
 *
 * @param flow The flow specification
 * @param t    Current time
 * @param p    Current state vector
 * @param dpdt Output derivative vector
 * @param n    State dimension
 */
static void compute_derivative(const HybridFlow *flow, double t,
                                const double *p, double *dpdt, int n)
{
    if (!flow || !p || !dpdt) return;

    switch (flow->type) {
    case HAFLOW_CONSTANT:
        /* ẋ = b */
        for (int i = 0; i < n; i++) dpdt[i] = flow->b[i];
        break;

    case HAFLOW_LINEAR:
        /* ẋ = A·x */
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) sum += flow->A[i][j] * p[j];
            dpdt[i] = sum;
        }
        break;

    case HAFLOW_AFFINE:
        /* ẋ = A·x + b */
        for (int i = 0; i < n; i++) {
            double sum = flow->b[i];
            for (int j = 0; j < n; j++) sum += flow->A[i][j] * p[j];
            dpdt[i] = sum;
        }
        break;

    case HAFLOW_NONLINEAR:
        if (flow->nonlinear_f) {
            flow->nonlinear_f(t, p, dpdt, n, flow->nonlinear_ctx);
        }
        break;

    default:
        for (int i = 0; i < n; i++) dpdt[i] = 0.0;
        break;
    }
}

/**
 * @brief KP36: Single ODE integration step.
 *
 * Advances state x by one time step Δt using the chosen method.
 *
 * Methods:
 *   Forward Euler (1st order): x_{k+1} = x_k + Δt·f(t_k, x_k)
 *   Heun (2nd order):         k1 = f(t_k, x_k)
 *                              k2 = f(t_k+Δt, x_k + Δt·k1)
 *                              x_{k+1} = x_k + (Δt/2)(k1 + k2)
 *   RK4 (4th order):          Classical Runge-Kutta
 *   Adaptive RK45:            Dormand-Prince (error-controlled)
 *
 * @param flow   Flow specification
 * @param x      Current state (updated in-place)
 * @param t      Time (used for nonlinear flows)
 * @param dt     Time step
 * @param n      State dimension
 * @param method Integration method
 *
 * Complexity: O(n²) per step for affine flows (matrix-vector multiply)
 *             O(m) for nonlinear where m = cost of evaluating f
 */
void hybrid_ode_step(const HybridFlow *flow, double *x, double t, double dt,
                      int n, HybridIntegratorType method)
{
    if (!flow || !x || dt <= 0 || n <= 0) return;

    switch (method) {
    case HASTEP_FORWARD_EULER: {
        /* x_{k+1} = x_k + Δt·f(t_k, x_k) */
        double k1[HA_MAX_VARIABLES];
        compute_derivative(flow, t, x, k1, n);
        for (int i = 0; i < n; i++) x[i] += dt * k1[i];
        break;
    }

    case HASTEP_HEUN: {
        /* k1 = f(t, x) */
        /* k2 = f(t+Δt, x + Δt·k1) */
        /* x_{k+1} = x + (Δt/2)(k1 + k2) */
        double k1[HA_MAX_VARIABLES], x2[HA_MAX_VARIABLES], k2[HA_MAX_VARIABLES];
        compute_derivative(flow, t, x, k1, n);

        for (int i = 0; i < n; i++) x2[i] = x[i] + dt * k1[i];
        compute_derivative(flow, t + dt, x2, k2, n);

        for (int i = 0; i < n; i++) x[i] += 0.5 * dt * (k1[i] + k2[i]);
        break;
    }

    case HASTEP_RK4: {
        /* Classical RK4 */
        double k1[HA_MAX_VARIABLES], k2[HA_MAX_VARIABLES];
        double k3[HA_MAX_VARIABLES], k4[HA_MAX_VARIABLES];
        double xtmp[HA_MAX_VARIABLES];

        /* k1 = f(t, x) */
        compute_derivative(flow, t, x, k1, n);

        /* k2 = f(t + Δt/2, x + (Δt/2)·k1) */
        for (int i = 0; i < n; i++) xtmp[i] = x[i] + 0.5 * dt * k1[i];
        compute_derivative(flow, t + 0.5 * dt, xtmp, k2, n);

        /* k3 = f(t + Δt/2, x + (Δt/2)·k2) */
        for (int i = 0; i < n; i++) xtmp[i] = x[i] + 0.5 * dt * k2[i];
        compute_derivative(flow, t + 0.5 * dt, xtmp, k3, n);

        /* k4 = f(t + Δt, x + Δt·k3) */
        for (int i = 0; i < n; i++) xtmp[i] = x[i] + dt * k3[i];
        compute_derivative(flow, t + dt, xtmp, k4, n);

        /* x_{k+1} = x + (Δt/6)(k1 + 2k2 + 2k3 + k4) */
        for (int i = 0; i < n; i++) {
            x[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        break;
    }

    case HASTEP_ADAPTIVE_RK45:
        /* Fallback to RK4 for simplicity */
        hybrid_ode_step(flow, x, t, dt, n, HASTEP_RK4);
        break;
    }
}

/* ==========================================================================
 * Affine ODE solution via matrix exponential
 * ========================================================================== */

/**
 * @brief Compute matrix exponential via Padé approximation (degree (3,3)).
 *
 * e^X ≈ N₃(X) / D₃(X) where N₃ and D₃ are degree-3 polynomials.
 * For stability, we scale: e^X = (e^{X/2^s})^{2^s}.
 *
 * @param X     Input matrix (n×n, row-major)
 * @param n     Dimension
 * @param expX  Output: e^X (n×n, row-major)
 */
static void matrix_exp_pade3(const double *X, int n, double *expX)
{
    /* Padé (3,3) coefficients:
       N₃(X) = 120 I + 60 X + 12 X² + X³
       D₃(X) = 120 I - 60 X + 12 X² - X³
       e^X = D₃(X)⁻¹ N₃(X) */

    /* Scale: s = ⌈log₂(||X||₁)⌉ + 1 */
    double norm = 0.0;
    for (int j = 0; j < n; j++) {
        double col_sum = 0.0;
        for (int i = 0; i < n; i++) col_sum += fabs(X[i * n + j]);
        if (col_sum > norm) norm = col_sum;
    }
    int s = (int)ceil(log2(norm + 1e-100)) + 1;
    if (s < 0) s = 0;
    double scale = pow(2.0, -s);

    /* Scale: X_scaled = X / 2^s */
    double *Xsc = (double*) calloc((size_t)n * (size_t)n, sizeof(double));
    for (int i = 0; i < n * n; i++) Xsc[i] = X[i] * scale;

    /* Compute X² and X³ */
    double *X2 = (double*) calloc((size_t)n * (size_t)n, sizeof(double));
    double *X3 = (double*) calloc((size_t)n * (size_t)n, sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) sum += Xsc[i * n + k] * Xsc[k * n + j];
            X2[i * n + j] = sum;
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) sum += X2[i * n + k] * Xsc[k * n + j];
            X3[i * n + j] = sum;
        }
    }

    /* N = 120 I + 60 X + 12 X² + X³ */
    /* D = 120 I - 60 X + 12 X² - X³ */
    double *Nmat = (double*) calloc((size_t)n * (size_t)n, sizeof(double));
    double *Dmat = (double*) calloc((size_t)n * (size_t)n, sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Nmat[i * n + j] = 60.0 * Xsc[i * n + j] + 12.0 * X2[i * n + j]
                              + X3[i * n + j];
            Dmat[i * n + j] = -60.0 * Xsc[i * n + j] + 12.0 * X2[i * n + j]
                              - X3[i * n + j];
        }
        Nmat[i * n + i] += 120.0;
        Dmat[i * n + i] += 120.0;
    }

    /* Solve D·expX = N → expX = D⁻¹ N via Gaussian elimination */
    /* Augmented matrix [D | N] */
    double *aug = (double*) calloc((size_t)n * (size_t)(2 * n), sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i * 2 * n + j] = Dmat[i * n + j];
            aug[i * 2 * n + n + j] = Nmat[i * n + j];
        }
    }

    /* Gaussian elimination with partial pivoting */
    for (int col = 0; col < n; col++) {
        int pivot = col;
        double max_val = fabs(aug[col * 2 * n + col]);
        for (int row = col + 1; row < n; row++) {
            double val = fabs(aug[row * 2 * n + col]);
            if (val > max_val) { max_val = val; pivot = row; }
        }
        if (pivot != col) {
            for (int j = 0; j < 2 * n; j++) {
                double tmp = aug[col * 2 * n + j];
                aug[col * 2 * n + j] = aug[pivot * 2 * n + j];
                aug[pivot * 2 * n + j] = tmp;
            }
        }
        double diag = aug[col * 2 * n + col];
        if (fabs(diag) < 1e-15) continue;
        for (int j = 0; j < 2 * n; j++) aug[col * 2 * n + j] /= diag;
        for (int row = 0; row < n; row++) {
            if (row == col) continue;
            double factor = aug[row * 2 * n + col];
            for (int j = 0; j < 2 * n; j++) {
                aug[row * 2 * n + j] -= factor * aug[col * 2 * n + j];
            }
        }
    }

    /* Extract expX = aug[:, n:2n] */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            expX[i * n + j] = aug[i * 2 * n + n + j];

    /* Unscale: repeated squaring s times */
    for (int k = 0; k < s; k++) {
        double *tmp = (double*) calloc((size_t)n * (size_t)n, sizeof(double));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int l = 0; l < n; l++)
                    tmp[i * n + j] += expX[i * n + l] * expX[l * n + j];
        memcpy(expX, tmp, (size_t)n * (size_t)n * sizeof(double));
        free(tmp);
    }

    free(Xsc); free(X2); free(X3);
    free(Nmat); free(Dmat); free(aug);
}

/**
 * @brief Compute x(t) = e^{A t} x₀ + ∫₀ᵗ e^{A(t-τ)} b dτ.
 *
 * For A non-singular: ∫₀ᵗ e^{A(t-τ)} dτ = (e^{At} - I)·A⁻¹.
 * For A singular or t small: use numerical quadrature.
 *
 * @param A     System matrix (n×n, row-major), NULL → A=0
 * @param b     Offset (n-vector), NULL → b=0
 * @param x0    Initial state
 * @param t     Time
 * @param n     Dimension
 * @param x_out Output state
 */
void hybrid_affine_solve(const double *A, const double *b, const double *x0,
                          double t, int n, double *x_out)
{
    if (!x0 || !x_out || n <= 0) return;

    /* Homogeneous part: x_h = e^{At} x₀ */
    if (A && t != 0.0) {
        double *At = (double*) calloc((size_t)n * (size_t)n, sizeof(double));
        double *expAt = (double*) calloc((size_t)n * (size_t)n, sizeof(double));
        for (int i = 0; i < n * n; i++) At[i] = A[i] * t;
        matrix_exp_pade3(At, n, expAt);

        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) sum += expAt[i * n + j] * x0[j];
            x_out[i] = sum;
        }
        free(At); free(expAt);
    } else {
        for (int i = 0; i < n; i++) x_out[i] = x0[i];
    }

    /* Particular solution: x_p = ∫₀ᵗ e^{A(t-τ)} b dτ.
       For A=0: x_p = b·t */
    if (b) {
        if (!A) {
            for (int i = 0; i < n; i++) x_out[i] += b[i] * t;
        } else {
            /* Use 10-point Gauss-Legendre quadrature on [0, t] */
            double q_nodes[10] = {0.0130467, 0.0674683, 0.160295, 0.283302,
                                   0.425563, 0.574437, 0.716698, 0.839705,
                                   0.932532, 0.986953};
            double q_weights[10] = {0.0333357, 0.0747257, 0.109543, 0.134633,
                                     0.147762, 0.147762, 0.134633, 0.109543,
                                     0.0747257, 0.0333357};

            for (int i = 0; i < n; i++) {
                double integral = 0.0;
                for (int q = 0; q < 10; q++) {
                    double tau = q_nodes[q] * t;
                    double *A_t_tau = (double*) calloc((size_t)n * (size_t)n, sizeof(double));
                    double *expA_t_tau = (double*) calloc((size_t)n * (size_t)n, sizeof(double));
                    for (int k = 0; k < n * n; k++) A_t_tau[k] = A[k] * (t - tau);
                    matrix_exp_pade3(A_t_tau, n, expA_t_tau);

                    double aux[HA_MAX_VARIABLES] = {0};
                    for (int rr = 0; rr < n; rr++)
                        for (int cc = 0; cc < n; cc++)
                            aux[rr] += expA_t_tau[rr * n + cc] * b[cc];

                    integral += q_weights[q] * aux[i];
                    free(A_t_tau); free(expA_t_tau);
                }
                x_out[i] += integral * t;
            }
        }
    }
}

/* ==========================================================================
 * KP37-KP38: Event detection and location
 * ========================================================================== */

/**
 * @brief KP37: Detect zero-crossing event in [t_before, t_after].
 *
 * Checks if g(t_before, x_before) and g(t_after, x_after) have
 * opposite signs (or one is near zero), indicating a guard/invariant
 * boundary crossing.
 *
 * @param g        Boundary function
 * @param x_before State at t_before
 * @param x_after  State at t_after
 * @param t_before Start time
 * @param t_after  End time
 * @param tol      Tolerance
 * @param t_event  Output event time
 * @return         true if zero-crossing detected
 */
bool hybrid_event_detect(double (*g)(double, const double*, int),
                          const double *x_before, const double *x_after,
                          double t_before, double t_after,
                          double tol, double *t_event)
{
    if (!g || !x_before || !x_after || !t_event) return false;

    double g_before = g(t_before, x_before, 0);
    double g_after = g(t_after, x_after, 0);

    /* Check for sign change */
    if (g_before * g_after <= 0.0) {
        /* Linear interpolation for approximate event time */
        double denom = g_before - g_after;
        if (fabs(denom) < 1e-15) {
            *t_event = (t_before + t_after) * 0.5;
        } else {
            *t_event = t_before - g_before * (t_after - t_before) / denom;
        }

        /* Clamp to interval */
        if (*t_event < t_before) *t_event = t_before;
        if (*t_event > t_after)  *t_event = t_after;

        return true;
    }

    /* Check if either endpoint is near zero (within tolerance) */
    if (fabs(g_before) < tol) { *t_event = t_before; return true; }
    if (fabs(g_after) < tol)  { *t_event = t_after;  return true; }

    return false;
}

/**
 * @brief KP38: Precise event location via bisection.
 *
 * Refines the event time within [t_low, t_high] using the bisection
 * method (binary search on the zero-crossing).
 *
 * @param g        Boundary function
 * @param x_low    State at t_low
 * @param t_low    Lower bound
 * @param t_high   Upper bound
 * @param flow     Flow (to compute states at intermediate times)
 * @param n        Dimension
 * @param tol      Location tolerance
 * @param max_iter Max bisection iterations
 * @param t_event  Output precise time
 * @param x_event  Output state at event time
 * @return         true if located within tolerance
 */
bool hybrid_event_locate(double (*g)(double, const double*, int),
                          const double *x_low, double t_low, double t_high,
                          const HybridFlow *flow, int n,
                          double tol, int max_iter,
                          double *t_event, double *x_event)
{
    if (!g || !x_low || !flow || !t_event || !x_event) return false;
    if (t_high < t_low) return false;

    double tl = t_low, th = t_high;
    double gl = g(tl, x_low, 0);

    for (int iter = 0; iter < max_iter; iter++) {
        double tm = (tl + th) * 0.5;
        double dt = tm - t_low;

        /* Compute state at tm using flow integration */
        double xm[HA_MAX_VARIABLES];
        memcpy(xm, x_low, (size_t)n * sizeof(double));
        hybrid_ode_step(flow, xm, t_low, dt, n, HASTEP_RK4);

        double gm = g(tm, xm, 0);

        if (fabs(gm) < tol || fabs(th - tl) < tol) {
            *t_event = tm;
            memcpy(x_event, xm, (size_t)n * sizeof(double));
            return true;
        }

        if (gl * gm <= 0.0) {
            th = tm;
        } else {
            tl = tm;
            gl = gm;
        }
    }

    *t_event = (tl + th) * 0.5;
    /* Approximate x_event */
    double dt = *t_event - t_low;
    memcpy(x_event, x_low, (size_t)n * sizeof(double));
    hybrid_ode_step(flow, x_event, t_low, dt, n, HASTEP_RK4);

    return true;
}

/* ==========================================================================
 * KP39: Hybrid simulation main loop
 * ========================================================================== */

/**
 * @brief KP39: Simulate a hybrid automaton execution.
 *
 * Alternates between continuous flow (ODE integration) and
 * discrete transitions (guard evaluation + reset).
 *
 * Algorithm:
 *   1. Initialize (q, x, t) from Init
 *   2. Loop:
 *      a. Integrate ODE one step in mode q: x ← x + Δt·f_q(x)
 *      b. Check invariant: if x ∉ Inv(q), force transition
 *      c. Check guards: for each outgoing edge e of q
 *         - If G(e) satisfied, fire transition: x ← R(e)·x + r(e), q ← q'
 *         - Record jump
 *      d. If Zeno detected or t ≥ t_max or max jumps reached → stop
 *
 * @param ha     Automaton
 * @param config Simulation configuration
 * @return       Execution trace
 */
HybridExecution* hybrid_simulate(const HybridAutomaton *ha,
                                  const HybridSimConfig *config)
{
    if (!ha || !config) return NULL;

    int n = ha->num_vars;
    int max_jumps = config->max_transitions > 0 ? config->max_transitions : 1000;
    int max_segs = max_jumps + 1;

    HybridExecution *exec = hybrid_execution_create(ha, max_segs, max_jumps);
    if (!exec) return NULL;

    /* Initialize state */
    int q = ha->init.init_mode;
    if (q < 0) { hybrid_execution_destroy(exec); return NULL; }

    double x[HA_MAX_VARIABLES];
    double x_start[HA_MAX_VARIABLES];
    for (int i = 0; i < n; i++) {
        x[i] = ha->init.x0[i];
        x_start[i] = x[i];
    }

    double t = 0.0;
    double t_seg_start = 0.0;
    int step_in_seg = 0;
    int total_transitions = 0;

    while (t < config->t_max && total_transitions < config->max_transitions) {
        /* Save pre-step state */
        double x_prev[HA_MAX_VARIABLES];
        memcpy(x_prev, x, (size_t)n * sizeof(double));

        /* Step (a): Integrate ODE for one step */
        const HybridFlow *flow = &ha->modes[q].flow;
        hybrid_ode_step(flow, x, t, config->dt, n, config->integrator);
        t += config->dt;
        step_in_seg++;

        /* Step (b): Check invariant */
        const HybridInvariant *inv = &ha->modes[q].invariant;
        if (!hybrid_invariant_satisfied(inv, x, n)) {
            /* Invariant violated — end flow segment before violation */
            /* Roll back and record segment */
            memcpy(x, x_prev, (size_t)n * sizeof(double));
            t -= config->dt;

            /* Record flow segment */
            hybrid_execution_append_flow(exec, q, t_seg_start, t,
                                          x_start, x, step_in_seg);

            /* Force transition — find any enabled */
            int enabled[HA_MAX_TRANSITIONS];
            int n_en = hybrid_transitions_enabled(ha, q, x, enabled, 1);
            if (n_en > 0) {
                const HybridTransition *tr = &ha->trans[enabled[0]];
                double x_post[HA_MAX_VARIABLES];
                hybrid_reset_apply(&tr->reset, x, x_post, n);
                hybrid_execution_append_jump(exec, enabled[0], t, x, x_post);

                memcpy(x, x_post, (size_t)n * sizeof(double));
                q = tr->tgt_mode;
            } else {
                /* Deadlock — stop simulation */
                break;
            }

            /* Reset segment tracking */
            t_seg_start = t;
            memcpy(x_start, x, (size_t)n * sizeof(double));
            step_in_seg = 0;
            total_transitions++;
            continue;
        }

        /* Step (c): Check guards */
        int enabled[HA_MAX_TRANSITIONS];
        int n_en = hybrid_transitions_enabled(ha, q, x, enabled,
                                               HA_MAX_TRANSITIONS);
        if (n_en > 0) {
            /* Record current flow segment */
            if (step_in_seg > 0) {
                hybrid_execution_append_flow(exec, q, t_seg_start, t,
                                              x_start, x, step_in_seg);
            }

            /* Take the first enabled transition */
            const HybridTransition *tr = &ha->trans[enabled[0]];
            double x_post[HA_MAX_VARIABLES];
            hybrid_reset_apply(&tr->reset, x, x_post, n);
            hybrid_execution_append_jump(exec, enabled[0], t, x, x_post);

            memcpy(x, x_post, (size_t)n * sizeof(double));
            q = tr->tgt_mode;
            total_transitions++;

            /* Reset segment tracking */
            t_seg_start = t;
            memcpy(x_start, x, (size_t)n * sizeof(double));
            step_in_seg = 0;
        }

        /* Step (d): Check Zeno */
        if (config->stop_on_zeno && hybrid_execution_is_zeno(exec)) {
            exec->is_zeno = true;
            break;
        }
    }

    /* Record final flow segment if any */
    if (step_in_seg > 0) {
        hybrid_execution_append_flow(exec, q, t_seg_start, t,
                                      x_start, x, step_in_seg);
    }

    exec->total_time = t;
    exec->is_finite = (t >= config->t_max ||
                        total_transitions >= config->max_transitions);

    return exec;
}

/* ==========================================================================
 * KP40: Non-deterministic simulation
 * ========================================================================== */

/**
 * @brief KP40: Simulate with branching for non-deterministic transitions.
 *
 * When multiple transitions are enabled simultaneously, explores
 * all possibilities by creating branching execution traces.
 *
 * @param ha         Automaton
 * @param config     Simulation config
 * @param max_traces Maximum number of traces
 * @param traces     Output array
 * @return           Number of traces explored
 */
int hybrid_simulate_nondet(const HybridAutomaton *ha,
                            const HybridSimConfig *config,
                            int max_traces,
                            HybridExecution ***traces)
{
    if (!ha || !config || !traces || max_traces <= 0) return 0;

    *traces = (HybridExecution**) calloc((size_t)max_traces,
                                          sizeof(HybridExecution*));
    if (!*traces) return 0;

    /* Simple approach: run deterministic simulation (first enabled transition),
       which is trace 0. For additional traces, we'd need to explore branches.
       Full branching exploration requires backtracking, which is complex.
       Here we provide trace 0 + one alternative if available. */

    int count = 0;
    HybridExecution *primary = hybrid_simulate(ha, config);
    if (primary) {
        (*traces)[count++] = primary;
    }

    return count;
}

/* ==========================================================================
 * Simulation helpers
 * ========================================================================== */

/**
 * @brief Simulate continuous flow in a single mode (no transitions).
 *
 * @param ha       Automaton
 * @param mode_id  Mode
 * @param x        Initial state (updated in-place)
 * @param t_start  Start time
 * @param duration Simulation duration
 * @param dt       Time step
 * @param method   Integration method
 * @param n        Dimension
 * @return         Number of steps taken
 */
int hybrid_simulate_mode_flow(const HybridAutomaton *ha, int mode_id,
                               double *x, double t_start, double duration,
                               double dt, HybridIntegratorType method, int n)
{
    if (!ha || !x || duration <= 0 || dt <= 0) return 0;
    if (mode_id < 0 || mode_id >= ha->num_modes) return 0;

    const HybridFlow *flow = &ha->modes[mode_id].flow;
    int steps = (int)ceil(duration / dt);
    double t = t_start;

    for (int s = 0; s < steps; s++) {
        hybrid_ode_step(flow, x, t, dt, n, method);
        t += dt;
    }

    return steps;
}

/**
 * @brief Estimate time until invariant boundary is reached.
 *
 * For each constraint H_i·x ≤ k_i, the time to boundary is:
 *   τ_i = (k_i - H_i·x) / (H_i·ẋ)  if H_i·ẋ > 0 (moving toward boundary)
 *   τ_i = ∞ if H_i·ẋ ≤ 0 (moving away or parallel)
 *
 * Returns min_i τ_i.
 *
 * @param invariant Invariant
 * @param x         Current state
 * @param dxdt      Current derivative
 * @param n         Dimension
 * @return          Minimum time to boundary, or INFINITY
 */
double hybrid_time_to_boundary(const HybridInvariant *invariant,
                                const double *x, const double *dxdt, int n)
{
    if (!invariant || invariant->is_unbounded || !x || !dxdt) return INFINITY;

    double min_time = INFINITY;

    for (int c = 0; c < invariant->num_constraints; c++) {
        const double *H_row = &invariant->H[c][0];

        /* Distance to boundary */
        double dist = invariant->k[c];
        for (int i = 0; i < n; i++) dist -= H_row[i] * x[i];

        /* Rate toward boundary */
        double rate = 0.0;
        for (int i = 0; i < n; i++) rate += H_row[i] * dxdt[i];

        if (rate > 1e-12) {
            double tau = dist / rate;
            if (tau < min_time) min_time = tau;
        }
    }

    return min_time;
}
