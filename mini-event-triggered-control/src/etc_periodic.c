#include "etc_periodic.h"
#include "etc_trigger.h"
#include "etc_dynamics.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * Periodic Event-Triggered Control (PETC)
 *
 * PETC evaluates the triggering condition only at periodic instants
 * t = k·h (k ∈ ℕ). This guarantees by construction that the minimum
 * inter-event time is at least h, eliminating Zeno behavior.
 *
 * Key difference from ETC:
 *   - ETC: continuous monitoring, event when Γ ≥ 0
 *   - PETC: periodic check at multiples of h, event when Γ ≥ 0
 *
 * Trade-off:
 *   + Zeno-free by construction (MIET ≥ h)
 *   + Implementable on digital hardware (fixed sampling clock)
 *   − Slightly more events than continuous ETC (delayed detection)
 *   − Requires h small enough to preserve stability
 * ============================================================================ */

PETCConfig etc_petc_config_create(double h) {
    PETCConfig cfg;
    cfg.h = h > 0.0 ? h : 0.01;
    cfg.t_last_check = 0.0;
    cfg.check_count = 0;
    cfg.event_count = 0;
    cfg.event_pending = false;
    cfg.last_trigger_value = -1.0;
    return cfg;
}

void etc_system_set_petc(ETCSystem* sys, double h) {
    /* PETC does not require structural changes to the ETC system.
     * The sampling period is enforced externally by the simulation loop.
     * However, we ensure the trigger type is compatible. */
    if (!sys) return;
    /* PETC works with any trigger type — the difference is in evaluation timing */
    (void)h;
}

bool etc_petc_should_check(const ETCSystem* sys, double h) {
    if (!sys || h <= 0.0) return false;
    return (sys->t - sys->t_last_event) >= h;
}

bool etc_petc_evaluate(ETCSystem* sys, PETCConfig* cfg) {
    if (!sys || !cfg) return false;

    cfg->check_count++;
    cfg->last_trigger_value = etc_system_eval_trigger(sys);

    if (cfg->last_trigger_value >= 0.0) {
        /* Trigger condition met: fire event */
        etc_system_trigger_event(sys);
        cfg->event_count++;
        cfg->event_pending = false;
        cfg->t_last_check = sys->t;
        return true;
    }

    cfg->event_pending = false;
    cfg->t_last_check = sys->t;
    return false;
}

void etc_petc_simulate(ETCSystem* sys, PETCConfig* cfg,
                        double t_max, double dt) {
    if (!sys || !cfg || t_max <= 0.0 || dt <= 0.0) return;
    if (dt > cfg->h) dt = cfg->h; /* dt must be ≤ h */

    sys->dt = dt;
    /* Record initial event */
    etc_system_trigger_event(sys);
    cfg->t_last_check = 0.0;

    while (sys->t < t_max) {
        /* Integrate one step */
        if (!etc_system_step_rk4(sys)) break;

        /* Check if it's time for a PETC check */
        if (etc_petc_should_check(sys, cfg->h)) {
            etc_petc_evaluate(sys, cfg);
        }
    }
}

double etc_petc_max_sampling_period(const ETCSystem* sys, double sigma) {
    if (!sys || sigma <= 0.0) return 0.0;

    int n = sys->n_states;
    ETCMatrix BK = etc_matrix_create(n, n);
    etc_matrix_mul(&sys->B, &sys->K, &BK);
    double norm_Acl = etc_matrix_norm_2(&sys->Acl);
    double norm_BK = etc_matrix_norm_2(&BK);
    etc_matrix_free(&BK);

    /* Heemels et al. (2013, Theorem III.2):
     * h_max = (1/||Acl||) · ln(1 + σ||Acl|| / (||Acl|| + σ||BK||))
     */
    if (norm_Acl < 1e-15) return INFINITY;

    double arg1 = sigma * norm_Acl / (norm_Acl + sigma * norm_BK);
    double h_max = (1.0 / norm_Acl) * log(1.0 + arg1);

    return h_max;
}

bool etc_petc_design(const ETCSystem* sys, double h, double* sigma) {
    if (!sys || !sigma || h <= 0.0) return false;

    int n = sys->n_states;
    ETCMatrix BK = etc_matrix_create(n, n);
    etc_matrix_mul(&sys->B, &sys->K, &BK);
    double norm_Acl = etc_matrix_norm_2(&sys->Acl);
    double norm_BK = etc_matrix_norm_2(&BK);
    etc_matrix_free(&BK);

    /* Inverse of h_max formula:
     * Given h, find σ such that h = (1/||Acl||) ln(1 + σ||Acl|| / (||Acl|| + σ||BK||))
     *
     * Let u = exp(h ||Acl||) − 1, then:
     *   σ = u ||Acl|| / (||Acl|| − u ||BK|| + u ||Acl||)
     *   σ = u ||Acl|| / (||Acl||(1+u) − u ||BK||)\n     *
     * This requires u ||BK|| < ||Acl||(1+u) for σ > 0.
     */
    double u = exp(h * norm_Acl) - 1.0;
    double denom = norm_Acl * (1.0 + u) - u * norm_BK;

    if (denom <= 0.0) {
        /* h is too large — no stabilizing σ exists */
        *sigma = 1.0;
        return false;
    }

    *sigma = u * norm_Acl / denom;

    /* Clamp to valid range */
    if (*sigma <= 0.0 || *sigma >= 1.0) {
        *sigma = 1.0;
        return false;
    }

    return true;
}

void etc_petc_compare(ETCSystem* sys, double h, double t_max, double dt,
                       double* ratio, int* petc_events, int* cont_events) {
    if (!sys || !ratio || !petc_events || !cont_events) return;

    /* Run PETC simulation */
    ETCSystem* sys_copy = etc_system_create(
        sys->A.data, sys->B.data, sys->K.data,
        sys->n_states, sys->n_inputs);
    if (!sys_copy) return;

    /* Copy state */
    ETCVector x0 = etc_vector_create(sys->n_states);
    for (int i = 0; i < sys->n_states; i++) x0.data[i] = sys->x.data[i];
    etc_system_set_initial_state(sys_copy, x0.data);
    etc_system_set_trigger(sys_copy, sys->trigger_type, sys->sigma,
                            sys->epsilon, sys->trigger_fn);
    etc_vector_free(&x0);

    PETCConfig cfg = etc_petc_config_create(h);
    etc_petc_simulate(sys_copy, &cfg, t_max, dt);
    *petc_events = cfg.event_count;

    /* Estimate continuous ETC events:
     * The expected number of events for continuous ETC is
     * approximately t_max / τ_min, where τ_min = σ / (||Acl|| + σ||BK||). */
    double tau_min = etc_inter_event_time_bound(sys);
    *cont_events = (int)(t_max / tau_min) + 1;

    *ratio = (*cont_events > 0) ? (double)(*petc_events) / (double)(*cont_events) : 1.0;

    etc_system_free(sys_copy);
}

bool etc_petc_is_zeno_free(double h) {
    /* PETC is Zeno-free by construction: events are separated by at least h.
     * The trigger is only checked at multiples of h, so even if the continuous-time
     * trigger would fire arbitrarily fast, only one event per period is possible. */
    return h > 0.0;
}
