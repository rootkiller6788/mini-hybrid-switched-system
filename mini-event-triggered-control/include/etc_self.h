#ifndef ETC_SELF_H
#define ETC_SELF_H

#include "etc_core.h"

/* ============================================================================
 * Self-Triggered Control (STC)
 *
 * Self-triggered control computes the next event time t_{k+1} = t_k + τ_k
 * at the current event k, using only the sampled state x(t_k) and the
 * system model. No continuous monitoring of the trigger condition is needed.
 *
 * Key advantages over ETC:
 *   — No continuous state measurement required
 *   — Fully predictable event schedule (offline verifiable)
 *   — Naturally robust to temporary sensor failures
 *
 * Disadvantage:
 *   — τ_k must be conservative (worst-case) ⇒ more events than ETC
 *
 * Core mechanism:
 *   At t = t_k:
 *     1. Measure state x(t_k)
 *     2. Compute τ_k = f(x(t_k)) using model-based prediction
 *     3. Schedule next event at t_{k+1} = t_k + τ_k
 *     4. Apply control u(t) = K x(t_k) for t ∈ [t_k, t_{k+1})
 *
 * References:
 *   Velasco, M. et al. (2003). Self-triggered control of networked
 *     control systems. IEEE CDC.
 *   Mazo, M. Jr. et al. (2010). On self-triggered control for linear
 *     systems: guarantees and complexity. ECC.
 *   Anta, A. & Tabuada, P. (2010). To sample or not to sample:
 *     self-triggered control for nonlinear systems. IEEE TAC.
 *   Almeida, J. et al. (2014). Self-triggered state feedback control
 *     of linear plants under bounded disturbances. IEEE TAC.
 * ============================================================================ */

/** Self-triggered control configuration. */
typedef struct {
    double tau_min;          /* Minimum inter-event time (MIET guarantee) */
    double tau_max;          /* Maximum inter-event time (safety bound) */
    double sigma;            /* Virtual threshold for τ computation */
    double prediction_horizon;
    int n_prediction_steps;  /* Discrete steps for τ search */
    ETCMatrix Acl;           /* Closed-loop matrix (cached) */
    ETCMatrix BK;            /* BK product (cached) */
    bool is_conservative;    /* True if using worst-case bound */
} STCConfig;

/** Create a self-triggered control configuration. */
STCConfig etc_stc_config_create(void);

/** Free resources associated with STC configuration. */
void etc_stc_config_free(STCConfig* cfg);

/** Initialize STC with system matrices.
 *  Pre-computes Acl = A + BK and BK product. */
void etc_stc_init(ETCSystem* sys, STCConfig* cfg);

/** Compute the next inter-event time τ_k = f(x_k) for STC.
 *
 *  Method: Predict state evolution under zero-order hold,
 *  find the smallest t > 0 such that |e(t)| = σ|x(t)|.
 *
 *  For linear systems, the error dynamics are:
 *     ė(t) = −ẋ(t) = −(Acl x(t) − BK e(t))
 *     e(0) = 0
 *  with x(t) = e^{Acl t} x_k + ∫₀ᵗ e^{Acl(t−s)} BK ds · x_k.
 *
 *  This function searches over [tau_min, tau_max] for the
 *  first time the trigger condition is met.
 *
 *  @param xk       Sampled state at current event
 *  @param cfg      STC configuration (with Acl, BK cached)
 *  @return         τ_k ∈ [tau_min, tau_max] */
double etc_stc_next_interval(const ETCVector* xk, const STCConfig* cfg);

/** Compute a worst-case (conservative) STC inter-event time.
 *  Uses global Lipschitz bounds instead of state-dependent prediction.
 *
 *  τ_wc = σ / (||Acl|| + σ||BK||)  (same as ETC lower bound)
 *
 *  This is always safe but generates more events than state-dependent STC.
 *
 *  @param cfg  STC configuration
 *  @return     τ_wc */
double etc_stc_worst_case_interval(const STCConfig* cfg);

/** Simulate a self-triggered control system.
 *
 *  At each event k:
 *    1. Sample x(t_k)
 *    2. Compute τ_k = etc_stc_next_interval(x(t_k), cfg)
 *    3. Integrate open-loop with u = K x(t_k) over [t_k, t_k + τ_k]
 *    4. t_{k+1} = t_k + τ_k, k ← k+1
 *
 *  @param sys     ETC system
 *  @param cfg     STC configuration
 *  @param t_max   Simulation horizon
 *  @param dt      Integration step size (≪ τ_min) */
void etc_stc_simulate(ETCSystem* sys, const STCConfig* cfg,
                       double t_max, double dt);

/** Compare STC against continuous ETC.
 *  STC is more conservative (more events) but requires no monitoring.
 *  This function quantifies the overhead.
 *
 *  @param sys           ETC system
 *  @param cfg           STC configuration
 *  @param t_max         Simulation horizon
 *  @param dt            Integration step
 *  @param overhead      Output: STC events / ETC events ratio (≥ 1)
 *  @param stc_events    Output: number of STC events
 *  @param etc_events    Output: estimated ETC events */
void etc_stc_compare(ETCSystem* sys, STCConfig* cfg, double t_max,
                      double dt, double* overhead,
                      int* stc_events, int* etc_events);

/** Adaptive self-triggered control: adjust τ_min online based on
 *  observed performance (state convergence rate).
 *
 *  @param sys           ETC system (must be in simulation)
 *  @param cfg           STC configuration (tau_min updated)
 *  @param current_time  Current simulation time
 *  @param current_x_norm Current state norm
 *  @param initial_x_norm Initial state norm */
void etc_stc_adapt(ETCSystem* sys, STCConfig* cfg,
                    double current_time, double current_x_norm,
                    double initial_x_norm);

/** STC with disturbance estimation.
 *
 *  When bounded disturbances w(t) are present (|w(t)| ≤ W),
 *  the self-triggered interval must be shortened to account for
 *  worst-case disturbance effects.
 *
 *  τ_k = max{ τ > 0 : |e(t)| ≤ σ|x(t)| for all t ∈ [0, τ]
 *                      under worst-case w(t) }
 *
 *  @param xk         Sampled state
 *  @param cfg        STC configuration
 *  @param W          Disturbance bound |w(t)| ≤ W
 *  @param Bw_norm    ||Bw|| (disturbance input matrix norm)
 *  @return           Robust τ_k */
double etc_stc_robust_interval(const ETCVector* xk, const STCConfig* cfg,
                                double W, double Bw_norm);

#endif /* ETC_SELF_H */
