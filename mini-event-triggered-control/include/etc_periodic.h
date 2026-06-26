#ifndef ETC_PERIODIC_H
#define ETC_PERIODIC_H

#include "etc_core.h"

/* ============================================================================
 * Periodic Event-Triggered Control (PETC)
 *
 * PETC evaluates the triggering condition only at periodic sampling
 * instants t = k·h (k ∈ N). This guarantees:
 *   — No Zeno behavior (sampling period h > 0 is a hard lower bound)
 *   — Implementable on digital platforms with fixed clock
 *   — Verifiable using discrete-time Lyapunov methods
 *
 * The PETC scheme:
 *   Given sampling period h > 0:
 *   t_{k+1} = min{ t = j·h > t_k : Γ(x(t), e(t)) ≥ 0 }
 *
 * Key references:
 *   Heemels, W.P.M.H. et al. (2013). Periodic event-triggered control
 *     based on state feedback. IEEE TAC, 58(4), 847–861.
 *   Postoyan, R. et al. (2015). Periodic event-triggered control for
 *     nonlinear systems. Automatica, 51, 200–206.
 *   Borgers, D.P. & Heemels, W.P.M.H. (2014). Event-separation
 *     properties of event-triggered control systems. IEEE TAC.
 * ============================================================================ */

/** PETC configuration.
 *  Extends basic ETC with a mandatory sampling period. */
typedef struct {
    double h;                  /* Sampling period (h > 0) */
    double t_last_check;       /* Last time trigger was evaluated */
    int check_count;           /* Number of trigger evaluations performed */
    int event_count;           /* Number of events triggered */
    bool event_pending;        /* True if trigger condition met at last check */
    double last_trigger_value; /* Γ value at last check */
} PETCConfig;

/** Create a PETC configuration with given sampling period. */
PETCConfig etc_petc_config_create(double h);

/** Set the PETC configuration on an ETC system.
 *  Converts the continuous trigger evaluation to periodic. */
void etc_system_set_petc(ETCSystem* sys, double h);

/** Check if it is time to evaluate the trigger (t ≥ t_last_check + h). */
bool etc_petc_should_check(const ETCSystem* sys, double h);

/** Evaluate the trigger at a PETC sampling instant.
 *  If Γ(x,e) ≥ 0, fires event and resets t_last_check.
 *  Otherwise, does nothing and schedules next check at t + h.
 *
 *  @param sys   ETC system
 *  @param cfg   PETC configuration (updated with check count, events)
 *  @return      True if event was triggered this check */
bool etc_petc_evaluate(ETCSystem* sys, PETCConfig* cfg);

/** PETC simulation loop: integrate with step dt, evaluate trigger every h.
 *
 *  @param sys     ETC system
 *  @param cfg     PETC configuration
 *  @param t_max   Simulation horizon
 *  @param dt      Integration step (dt ≤ h) */
void etc_petc_simulate(ETCSystem* sys, PETCConfig* cfg,
                        double t_max, double dt);

/** Compute the maximum allowable sampling period h for PETC.
 *  The PETC with sampling period h guarantees the same stability
 *  as continuous ETC if h ≤ h_max, where h_max depends on the
 *  system dynamics and Lyapunov function.
 *
 *  From Heemels et al. (2013, Theorem III.2):
 *     h_max = (1/||Acl||) · ln(1 + σ||Acl|| / (||Acl|| + σ||BK||))
 *
 *  @param sys   ETC system with A, B, K set
 *  @param sigma Trigger threshold
 *  @return      h_max (maximum sampling period preserving stability) */
double etc_petc_max_sampling_period(const ETCSystem* sys, double sigma);

/** Design a PETC with guaranteed stability.
 *  Given desired sampling period h, computes the required σ
 *  to ensure asymptotic stability.
 *
 *  @param sys   ETC system
 *  @param h     Desired sampling period
 *  @param sigma Output: required threshold σ
 *  @return      True if a stabilizing σ ∈ (0, 1) exists for given h */
bool etc_petc_design(const ETCSystem* sys, double h, double* sigma);

/** Compare PETC behavior against continuous ETC.
 *  Computes the ratio of events triggered under PETC vs continuous ETC.
 *
 *  @param sys     ETC system
 *  @param h       PETC sampling period
 *  @param t_max   Simulation horizon
 *  @param dt      Integration step
 *  @param ratio   Output: PETC events / continuous ETC events
 *  @param petc_events    Output: number of PETC events
 *  @param cont_events    Output: estimated continuous ETC events */
void etc_petc_compare(ETCSystem* sys, double h, double t_max, double dt,
                       double* ratio, int* petc_events, int* cont_events);

/** Check if PETC induces Zeno behavior.
 *  By construction, PETC cannot be Zeno since events are separated
 *  by at least h. Returns false always for valid h > 0. */
bool etc_petc_is_zeno_free(double h);

#endif /* ETC_PERIODIC_H */
