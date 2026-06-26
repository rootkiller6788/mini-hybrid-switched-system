#ifndef ETC_DYNAMICS_H
#define ETC_DYNAMICS_H

#include "etc_core.h"

/* ============================================================================
 * Event-Triggered Control — Dynamics & Simulation
 *
 * This module implements the numerical integration of event-triggered
 * control systems. Key aspects:
 *
 * 1. Continuous-time plant dynamics between events (RK4 integration)
 * 2. Trigger condition monitoring at each integration step
 * 3. Event detection and state reset (zero-order hold update)
 * 4. Inter-event time statistics collection
 * 5. Flow/impulse decomposition of hybrid trajectory
 *
 * References:
 *   Heemels, W.P.M.H. et al. (2012). Periodic event-triggered control.
 *   Dolk, V.S. et al. (2017). Event-triggered control systems under
 *     Denial-of-Service attacks. IEEE TAC.
 *   Postoyan, R. et al. (2015). Event-triggered tracking control.
 *     Automatica, 53, 275–281.
 * ============================================================================ */

/* --- Simulation Configuration --- */

/** Simulation configuration parameters. */
typedef struct {
    double t_start;          /* Start time */
    double t_end;            /* End time (horizon) */
    double dt;               /* Integration step size */
    double dt_max;           /* Maximum allowed step (adaptive) */
    double dt_min;           /* Minimum allowed step */
    int max_events;          /* Maximum number of events before abort */
    int max_steps;           /* Maximum number of integration steps */
    bool detect_zeno;        /* Enable Zeno detection */
    double zeno_threshold;   /* Inter-event time below which Zeno suspected */
    double event_tolerance;  /* Tolerance for trigger zero-crossing */
} ETCSimConfig;

/** Simulation results summary. */
typedef struct {
    double t_final;          /* Final simulation time */
    int n_steps;             /* Total integration steps taken */
    int n_events;            /* Total events triggered */
    int n_zeno_detections;   /* Number of Zeno-like clusters detected */
    double min_iet;          /* Minimum observed inter-event time */
    double max_iet;          /* Maximum observed inter-event time */
    double avg_iet;          /* Average inter-event time */
    double std_iet;          /* Standard deviation of inter-event times */
    double total_control_updates;
    double avg_error_norm;   /* Average measurement error norm */
    double avg_state_norm;   /* Average state norm */
    double settling_time;    /* Time to reach 2% of initial state */
    bool converged;          /* True if state converged to origin */
    bool zeno_free;          /* True if no Zeno detected */
} ETCSimResult;

/* --- Dynamics Kernel --- */

/** Compute closed-loop dynamics ẋ = Acl x − BK e.
 *  Separates nominal (Acl x) and error-driven (−BK e) components.
 *
 *  @param sys  ETC system with Acl = A+BK set
 *  @param dx   Output: state derivative (n_states elements) */
void etc_dynamics_closed_loop(const ETCSystem* sys, double* dx);

/** Compute open-loop plant dynamics ẋ = A x + B u.
 *
 *  @param sys  ETC system
 *  @param dx   Output: state derivative */
void etc_dynamics_open_loop(const ETCSystem* sys, double* dx);

/** Compute the flow map Φ_t(x₀, u) for the linear system under
 *  constant input u (zero-order hold).
 *
 *  Φ_t(x₀, u) = e^{A t} x₀ + ∫₀ᵗ e^{A(t−s)} B ds · u
 *
 *  Uses matrix exponential via Padé approximation (degree 6).
 *
 *  @param A       State matrix (n×n)
 *  @param B       Input matrix (n×m)
 *  @param x0      Initial state (n×1)
 *  @param u_const Constant control input (m×1)
 *  @param t       Time horizon
 *  @param n       State dimension
 *  @param m       Input dimension
 *  @param x_t     Output: Φ_t(x₀, u) */
void etc_flow_map(const ETCMatrix* A, const ETCMatrix* B,
                   const ETCVector* x0, const ETCVector* u_const,
                   double t, int n, int m, ETCVector* x_t);

/** Compute matrix exponential e^{A t} using truncated Taylor series.
 *
 *  e^{A t} = Σ_{k=0}^{N} (A t)^k / k!
 *
 *  Uses scaling-and-squaring for numerical stability.
 *
 *  @param A       Square matrix (n×n)
 *  @param t       Time scalar
 *  @param n       Dimension
 *  @param expAt   Output: e^{A t} (n×n), must be pre-allocated */
void etc_matrix_exponential(const ETCMatrix* A, double t, int n,
                             ETCMatrix* expAt);

/* --- Simulation Engine --- */

/** Initialize simulation configuration with sensible defaults. */
void etc_sim_config_init(ETCSimConfig* cfg);

/** Run event-triggered simulation with detailed logging.
 *
 *  At each integration step:
 *    1. Integrate plant dynamics (RK4)
 *    2. Compute measurement error e = x̂ − x
 *    3. Evaluate trigger Γ(x, e)
 *    4. If Γ ≥ 0: fire event (x̂ ← x, e ← 0, record event)
 *    5. Collect statistics
 *
 *  @param sys    ETC system (must be fully configured)
 *  @param cfg    Simulation configuration
 *  @param result Output: simulation results summary */
void etc_simulate_full(const ETCSystem* sys, const ETCSimConfig* cfg,
                        ETCSimResult* result);

/** Run a single integration time step for the event-triggered system.
 *  Evaluates trigger after integration and fires event if needed.
 *
 *  @param sys    ETC system
 *  @param dt     Step size
 *  @param event_fired  Output: true if event was triggered this step
 *  @return            True if step successful */
bool etc_step_integrate(ETCSystem* sys, double dt, bool* event_fired);

/** Detect zero-crossing of trigger function within interval [t, t+dt].
 *  Uses linear interpolation for efficient root-finding.
 *
 *  @param sys         ETC system at current state
 *  @param dt          Step size
 *  @param event_time  Output: precise event time (interpolated)
 *  @return            True if zero-crossing detected */
bool etc_detect_event(ETCSystem* sys, double dt, double* event_time);

/* --- Inter-Event Time Analysis --- */

/** Compute the theoretical lower bound on inter-event time.
 *
 *  For static trigger Γ = |e| − σ|x|:
 *     τ_min = σ / (||Acl|| + σ||BK||)
 *
 *  This is a global lower bound valid for all events.
 *
 *  @param sys  ETC system with Acl=BK set
 *  @return     Positive lower bound τ_min */
double etc_inter_event_time_bound(const ETCSystem* sys);

/** Compute inter-event time as a function of the sampled state.
 *
 *  For linear systems, the next event time can be estimated from
 *  the evolution of the error dynamics: ė = −Acl x − BK e.
 *
 *  @param sys   ETC system
 *  @param xk    Sampled state at event k
 *  @return      Estimated next inter-event time */
double etc_iet_estimate(const ETCSystem* sys, const ETCVector* xk);

/** Compute the average inter-event time from event history.
 *  Useful for communication rate estimation.
 *
 *  @param hist  Event history (populated via simulation)
 *  @return      Average Δ_k = t_{k} − t_{k-1} */
double etc_average_iet(const ETCEventHistory* hist);

/** Compute the event-triggering communication ratio.
 *  Ratio = (number of events) / (number of periodic samples that would
 *  be needed at the same integration step size).
 *
 *  @param n_events   Number of events
 *  @param t_total    Total simulation time
 *  @param dt         Integration step size
 *  @return           Communication ratio in [0, 1] */
double etc_communication_ratio(int n_events, double t_total, double dt);

/** Print simulation results summary to stdout. */
void etc_sim_result_print(const ETCSimResult* result);

/* --- Adaptive Step-Size Control --- */

/** Adjust integration step size based on trigger margin.
 *  When Γ ≪ 0 (far from triggering), step size can increase.
 *  When Γ ≈ 0 (near triggering), step size must decrease.
 *
 *  @param sys        ETC system
 *  @param dt_current Current step size
 *  @param dt_next    Output: suggested next step size
 *  @return           True if step size adjustment applied */
bool etc_adaptive_step(const ETCSystem* sys, double dt_current,
                        double* dt_next);

#endif /* ETC_DYNAMICS_H */
