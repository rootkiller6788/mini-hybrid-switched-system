#include "etc_self.h"
#include "etc_trigger.h"
#include "etc_dynamics.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Self-Triggered Control (STC)
 *
 * STC pre-computes the next event time at the current event, eliminating
 * the need for continuous state monitoring. The predicted inter-event time
 * τ_k = f(x(t_k)) must guarantee that |e(t)| ≤ σ|x(t)| for all t ∈ [0, τ_k].
 *
 * The computation of τ_k involves:
 *   1. Predicting x(t) = e^{Acl t} x_k under zero-order hold
 *   2. Computing e(t) = x_k − x(t) (error evolution)
 *   3. Finding the first t where |e(t)| = σ|x(t)|
 *
 * For linear systems, this reduces to a 1D root-finding problem.
 * ============================================================================ */

/* Internal: evaluate the trigger at a predicted future time under STC */
static double etc_stc_predict_trigger(const ETCVector* xk,
                                       const ETCMatrix* Acl,
                                       const ETCMatrix* BK,
                                       double t, double sigma) {
    int n = Acl->rows;
    /* Predict x(t) = e^{Acl t} x_k + ∫₀ᵗ e^{Acl(t−s)} ds · BK x_k */
    /* For small time increments, we use the flow map approximation */

    /* Compute the flow map for both the autonomous and forced parts.
     * Since the control is u = K x_k (constant), the forced response is:
     *   x(t) = e^{Acl t} x_k  [since e^{At} B K x_k ds integral
     *                          with Acl = A + BK gives exactly this] */
    ETCVector x_t = etc_vector_create(n);
    ETCVector u_const = etc_vector_create(n);
    etc_matrix_vec_mul(BK, xk, &u_const);

    etc_flow_map(Acl, BK, xk, xk, t, n, n, &x_t);
    /* Note: for the flow map with Acl, the input is K x_k acting through BK.
     * The actual dynamics are ẋ = Acl x − BK e, which under zero error
     * becomes ẋ = Acl x. With e ≠ 0, we need a more careful analysis.
     *
     * Simplified approach: use the closed-form from the linear dynamics:
     *   ẋ = A x + B K x_k = Acl x − BK (x − x_k) = Acl x − BK e
     *
     * We compute e(t) ≈ −∫₀ᵗ Acl e^{Acl s} x_k ds (first-order approx) */

    ETCVector e_t = etc_vector_create(n);
    /* e(t) = x_k − x(t) */
    for (int i = 0; i < n; i++)
        e_t.data[i] = xk->data[i] - x_t.data[i];

    double trigger_val = etc_trigger_quadratic(&x_t, &e_t, sigma, 0.0);

    etc_vector_free(&x_t);
    etc_vector_free(&u_const);
    etc_vector_free(&e_t);

    return trigger_val;
}

/* ============================================================================
 * STC Public API
 * ============================================================================ */

STCConfig etc_stc_config_create(void) {
    STCConfig cfg;
    cfg.tau_min = 0.001;
    cfg.tau_max = 1.0;
    cfg.sigma = 0.1;
    cfg.prediction_horizon = 5.0;
    cfg.n_prediction_steps = 100;
    cfg.is_conservative = true;
    cfg.Acl = etc_matrix_create(1, 1); /* initial; resized in etc_stc_init */
    cfg.BK = etc_matrix_create(1, 1);
    return cfg;
}

void etc_stc_config_free(STCConfig* cfg) {
    if (cfg) {
        etc_matrix_free(&cfg->Acl);
        etc_matrix_free(&cfg->BK);
    }
}

void etc_stc_init(ETCSystem* sys, STCConfig* cfg) {
    if (!sys || !cfg) return;
    int n = sys->n_states;

    /* Cache Acl and BK */
    etc_matrix_free(&cfg->Acl);
    etc_matrix_free(&cfg->BK);
    cfg->Acl = etc_matrix_create(n, n);
    cfg->BK = etc_matrix_create(n, n);

    for (int i = 0; i < n * n; i++) cfg->Acl.data[i] = sys->Acl.data[i];
    etc_matrix_mul(&sys->B, &sys->K, &cfg->BK);

    cfg->sigma = sys->sigma;

    /* Compute worst-case tau_min from ETC lower bound */
    double tau_wc = etc_inter_event_time_bound(sys);
    if (tau_wc < cfg->tau_min) cfg->tau_min = tau_wc * 0.5;
    if (cfg->tau_min < 1e-6) cfg->tau_min = 1e-6;
}

double etc_stc_next_interval(const ETCVector* xk, const STCConfig* cfg) {
    if (!xk || !cfg) return cfg->tau_max;

    double sigma = cfg->sigma;
    double tau_min = cfg->tau_min;
    double tau_max = cfg->tau_max;

    /* Search for the first time t ∈ [tau_min, tau_max] where trigger fires.
     * Use binary search / grid search for robustness. */
    int n_steps = cfg->n_prediction_steps;
    double tau = tau_max;
    bool found = false;

    /* Grid search from tau_min to tau_max */
    for (int k = 0; k <= n_steps; k++) {
        double t = tau_min + (tau_max - tau_min) * (double)k / (double)n_steps;
        double gamma = etc_stc_predict_trigger(xk, &cfg->Acl, &cfg->BK, t, sigma);
        if (gamma >= 0.0) {
            tau = t;
            found = true;
            break;
        }
    }

    /* If no trigger found within tau_max, use tau_max */
    if (!found) return tau_max;

    /* Refine with bisection around the found tau */
    double t_lo = (tau > tau_min + 1e-8) ? tau_min : tau_min;
    double t_hi = tau;
    for (int iter = 0; iter < 20; iter++) {
        double t_mid = 0.5 * (t_lo + t_hi);
        double gamma = etc_stc_predict_trigger(xk, &cfg->Acl, &cfg->BK, t_mid, sigma);
        if (gamma >= 0.0)
            t_hi = t_mid;
        else
            t_lo = t_mid;
    }
    return 0.5 * (t_lo + t_hi);
}

double etc_stc_worst_case_interval(const STCConfig* cfg) {
    if (!cfg) return 0.001;
    /* Worst-case bound from Tabuada (2007):
     * τ_wc = σ / (||Acl|| + σ||BK||) */
    double norm_Acl = etc_matrix_norm_2(&cfg->Acl);
    double norm_BK = etc_matrix_norm_2(&cfg->BK);
    double denom = norm_Acl + cfg->sigma * norm_BK;
    if (denom < 1e-15) return cfg->tau_max;
    return cfg->sigma / denom;
}

void etc_stc_simulate(ETCSystem* sys, const STCConfig* cfg,
                       double t_max, double dt) {
    if (!sys || !cfg || t_max <= 0.0 || dt <= 0.0) return;

    int n = sys->n_states;
    sys->dt = dt;

    /* Initial event (k = 0) */
    etc_vector_copy(&sys->x, &sys->x_hat);
    etc_system_compute_control(sys);
    etc_system_trigger_event(sys);

    /* Main STC loop */
    double tau = etc_stc_next_interval(&sys->x_hat, cfg);
    double t_next_event = sys->t + tau;

    while (sys->t < t_max) {
        /* Integrate with fixed control u = K x_hat */
        if (sys->t + dt >= t_next_event) {
            /* Step exactly to event time */
            double dt_to_event = t_next_event - sys->t;
            if (dt_to_event > 0.0) {
                sys->dt = dt_to_event;
                etc_system_step_rk4(sys);
                sys->dt = dt; /* restore */
            }

            /* Fire event */
            etc_vector_copy(&sys->x, &sys->x_hat);
            etc_system_compute_control(sys);
            etc_system_trigger_event(sys);

            /* Compute next interval */
            tau = etc_stc_next_interval(&sys->x_hat, cfg);
            t_next_event = sys->t + tau;

            if (t_next_event > t_max) break;
        } else {
            /* Normal integration step */
            if (!etc_system_step_rk4(sys)) break;
            /* Update error for tracking */
            for (int i = 0; i < n; i++)
                sys->e.data[i] = sys->x_hat.data[i] - sys->x.data[i];
        }
    }
}

void etc_stc_compare(ETCSystem* sys, STCConfig* cfg, double t_max,
                      double dt, double* overhead,
                      int* stc_events, int* etc_events) {
    if (!sys || !cfg || !overhead || !stc_events || !etc_events) return;

    /* Run STC simulation */
    ETCSystem* sys_stc = etc_system_create(
        sys->A.data, sys->B.data, sys->K.data,
        sys->n_states, sys->n_inputs);
    if (!sys_stc) return;

    ETCVector x0 = etc_vector_create(sys->n_states);
    for (int i = 0; i < sys->n_states; i++) x0.data[i] = sys->x.data[i];
    etc_system_set_initial_state(sys_stc, x0.data);
    etc_system_set_trigger(sys_stc, sys->trigger_type,
                            sys->sigma, sys->epsilon, sys->trigger_fn);
    etc_vector_free(&x0);

    STCConfig stc_cfg = etc_stc_config_create();
    etc_stc_init(sys_stc, &stc_cfg);
    stc_cfg.sigma = sys->sigma;

    etc_stc_simulate(sys_stc, &stc_cfg, t_max, dt);
    *stc_events = sys_stc->event_count;

    /* Estimate continuous ETC events */
    double tau_min = etc_inter_event_time_bound(sys);
    *etc_events = (int)(t_max / tau_min) + 1;

    *overhead = (*etc_events > 0)
                ? (double)(*stc_events) / (double)(*etc_events)
                : 1.0;

    etc_stc_config_free(&stc_cfg);
    etc_system_free(sys_stc);
}

void etc_stc_adapt(ETCSystem* sys, STCConfig* cfg,
                    double current_time, double current_x_norm,
                    double initial_x_norm) {
    if (!sys || !cfg) return;
    (void)current_time;

    /* Adaptation logic:
     * As the state converges (|x| → 0), events become less frequent.
     * We can increase tau_min to reduce the number of events.
     *
     * Convergence ratio ∈ [0, 1]: 0 = initial, 1 = converged */
    double convergence = (initial_x_norm > 1e-10)
                         ? 1.0 - current_x_norm / initial_x_norm
                         : 1.0;

    /* Scale tau_min: increase as system converges */
    double scale = 1.0 + 9.0 * convergence; /* 1x → 10x */
    cfg->tau_min = cfg->tau_min * scale;
    if (cfg->tau_min > cfg->tau_max * 0.5)
        cfg->tau_min = cfg->tau_max * 0.5;
}

double etc_stc_robust_interval(const ETCVector* xk, const STCConfig* cfg,
                                double W, double Bw_norm) {
    if (!xk || !cfg) return cfg->tau_min;

    /* With bounded disturbance |w(t)| ≤ W, the error dynamics become:
     *
     *   ė = −Acl x − BK e − Bw w
     *
     * The disturbance adds a term Bw w that accelerates error growth.
     * A conservative bound on the inter-event time:
     *
     *   τ_robust ≤ σ · (|xk|) / (||Acl|| |xk| + σ||BK|| |xk| + ||Bw|| W)
     *
     * which is smaller than the nominal τ.
     */
    double xk_norm = etc_vector_norm(xk);
    if (xk_norm < 1e-12) return cfg->tau_max;

    double norm_Acl = etc_matrix_norm_2(&cfg->Acl);
    double norm_BK = etc_matrix_norm_2(&cfg->BK);
    double sigma = cfg->sigma;

    double denom = norm_Acl * xk_norm + sigma * norm_BK * xk_norm + Bw_norm * W;
    if (denom < 1e-15) return cfg->tau_max;

    double tau_robust = sigma * xk_norm / denom;

    if (tau_robust < cfg->tau_min) tau_robust = cfg->tau_min;
    if (tau_robust > cfg->tau_max) tau_robust = cfg->tau_max;

    return tau_robust;
}
