/**
 * @file hybrid_simulation.h
 * @brief Simulation engine for hybrid automata (L5-L6)
 *
 * Numerical simulation of hybrid automaton executions using ODE
 * integration within modes and discrete transition detection.
 *
 * Reference:
 *   Esposito, Kumar, Pappas, "Accurate Event Detection for Simulating
 *     Hybrid Systems" (2001), HSCC
 *   Park, Barton, "State Event Location in Differential-Algebraic Models"
 *     (1996), ACM TOMS
 *
 * Course Mapping:
 *   MIT 6.841    - Numerical simulation
 *   Berkeley EECS 291E - Simulation-based analysis
 *   CMU 15-455   - Simulation algorithms
 *
 * Knowledge points:
 *   KP36: Time-triggered ODE integration (fixed step)
 *   KP37: Event detection (zero-crossing of guards/invariants)
 *   KP38: Event location (precise transition time)
 *   KP39: Hybrid simulation main loop
 *   KP40: Fixed-point iteration for transition enablement
 */

#ifndef HYBRID_SIMULATION_H
#define HYBRID_SIMULATION_H

#include "hybrid_automaton.h"
#include "hybrid_execution.h"
#include <stdbool.h>

/* ==========================================================================
 * ODE integrator types
 * ========================================================================== */

typedef enum {
    HASTEP_FORWARD_EULER,    /**< Forward Euler: x_{k+1} = x_k + Δt·f(t_k, x_k) */
    HASTEP_RK4,              /**< Classical Runge-Kutta 4th order */
    HASTEP_HEUN,             /**< Heun's method (2nd order) */
    HASTEP_ADAPTIVE_RK45     /**< Dormand-Prince adaptive step (RK45) */
} HybridIntegratorType;

/* ==========================================================================
 * Simulation configuration
 * ========================================================================== */

/**
 * @brief Configuration for hybrid simulation.
 */
typedef struct {
    HybridIntegratorType integrator;  /**< ODE integration method */
    double  dt;                       /**< Base time step for integration */
    double  t_max;                    /**< Maximum simulation time */
    int     max_transitions;          /**< Maximum discrete transitions (safety limit) */
    double  guard_tolerance;          /**< Tolerance for guard zero-crossing detection */
    double  invariant_tolerance;      /**< Tolerance for invariant boundary detection */
    bool    stop_on_zeno;             /**< Stop simulation if Zeno detected */
    int     verbosity;               /**< Verbosity level */
} HybridSimConfig;

/** Default simulation configuration */
#define HYBRID_SIM_CONFIG_DEFAULT { \
    .integrator = HASTEP_RK4, \
    .dt = 0.001, \
    .t_max = 100.0, \
    .max_transitions = 10000, \
    .guard_tolerance = 1e-6, \
    .invariant_tolerance = 1e-6, \
    .stop_on_zeno = true, \
    .verbosity = 0 \
}

/* ==========================================================================
 * ODE integration (KP36)
 * ========================================================================== */

/**
 * @brief KP36: Single ODE integration step.
 *
 * Advances the continuous state by one time step Δt using the
 * specified integration method.
 *
 * For affine ODE ẋ = A·x + b, this computes:
 *   Forward Euler: x(t+Δt) = x(t) + Δt·(A·x(t) + b)
 *   RK4: k1 = Δt·f(t, x), k2 = Δt·f(t+Δt/2, x+k1/2), ...
 *
 * @param flow    The flow specification (affine or nonlinear)
 * @param x       Current state (n-dimensional, updated in-place)
 * @param t       Current time (updated to t+dt on nonlinear flows)
 * @param dt      Time step
 * @param n       State dimension
 * @param method  Integration method
 */
void hybrid_ode_step(const HybridFlow *flow, double *x, double t, double dt,
                      int n, HybridIntegratorType method);

/**
 * @brief Compute affine ODE solution x(t) = e^{A t} x₀ + ∫₀ᵗ e^{A(t-τ)} b dτ.
 *
 * Uses matrix exponential (Padé approximation) for the homogeneous part
 * and numerical quadrature for the particular solution.
 *
 * @param A     System matrix (n×n, row-major), may be NULL (zero matrix)
 * @param b     Affine offset (n-vector), may be NULL (zero vector)
 * @param x0    Initial state (n-vector)
 * @param t     Time to advance
 * @param n     State dimension
 * @param x_out Output state at time t (n-vector)
 *
 * Complexity: O(n³) for matrix exponential (Padé)
 */
void hybrid_affine_solve(const double *A, const double *b, const double *x0,
                          double t, int n, double *x_out);

/* ==========================================================================
 * Event detection (KP37-KP38)
 * ========================================================================== */

/**
 * @brief KP37: Detect a zero-crossing event within a time interval.
 *
 * A zero-crossing occurs when a guard or invariant boundary function
 * g(t) changes sign, indicating that a discrete transition may be enabled
 * or that the invariant may be violated.
 *
 * Uses linear interpolation between integration steps.
 *
 * @param g          Boundary function: g(t, x) ≤ 0 means constraint holds
 * @param x_before   State before the step
 * @param x_after    State after the step
 * @param t_before   Time before the step
 * @param t_after    Time after the step
 * @param tol        Detection tolerance
 * @param t_event    Output: estimated event time (if event detected)
 * @return           true if a zero-crossing was detected
 */
bool hybrid_event_detect(double (*g)(double, const double*, int),
                          const double *x_before, const double *x_after,
                          double t_before, double t_after,
                          double tol, double *t_event);

/**
 * @brief KP38: Precise event location via bisection or secant method.
 *
 * Given a zero-crossing detected in [t_low, t_high], refine the
 * exact event time using iterative root-finding.
 *
 * @param g        Boundary function
 * @param x_low    State at t_low
 * @param t_low    Lower bound of interval
 * @param t_high   Upper bound of interval
 * @param flow     Flow specification (to compute intermediate states)
 * @param n        State dimension
 * @param tol      Location tolerance
 * @param max_iter Maximum iterations
 * @param t_event  Output: precise event time
 * @param x_event  Output: state at event time (n-vector)
 * @return         true if event located within tolerance
 */
bool hybrid_event_locate(double (*g)(double, const double*, int),
                          const double *x_low, double t_low, double t_high,
                          const HybridFlow *flow, int n,
                          double tol, int max_iter,
                          double *t_event, double *x_event);

/* ==========================================================================
 * Simulation engine (KP39-KP40)
 * ========================================================================== */

/**
 * @brief KP39: Main hybrid simulation loop.
 *
 * Executes the hybrid automaton from its initial state, alternating
 * between continuous flow (ODE integration) and discrete transitions
 * (guard evaluation and reset application).
 *
 * Algorithm:
 *   1. (q, x, t) ← initial state from Init
 *   2. While t < t_max and transitions < max_transitions:
 *      a. Integrate ODE ẋ = f_q(x) for one time step
 *      b. Check invariants: if x ∉ Inv(q) → error (must force transition earlier)
 *      c. Check guards: if any e = (q → q') enabled:
 *         - Apply reset: x ← R(e)·x + r(e)
 *         - Update mode: q ← q'
 *         - Increment transition counter
 *      d. If Zeno detected → stop
 *   3. Return execution trace
 *
 * @param ha     The hybrid automaton to simulate
 * @param config Simulation configuration
 * @return       Execution trace, or NULL on error
 */
HybridExecution* hybrid_simulate(const HybridAutomaton *ha,
                                  const HybridSimConfig *config);

/**
 * @brief KP40: Simulate with fixed-point iteration for transition enablement.
 *
 * For non-deterministic automata, when multiple transitions are enabled,
 * this variant explores all enabled transitions by creating branching
 * execution trees.
 *
 * @param ha       The automaton
 * @param config   Simulation config
 * @param max_traces Maximum number of execution traces to explore
 * @param traces   Output: array of execution traces
 * @return         Number of traces actually explored
 */
int  hybrid_simulate_nondet(const HybridAutomaton *ha,
                             const HybridSimConfig *config,
                             int max_traces,
                             HybridExecution ***traces);

/**
 * @brief Simulate a single mode's continuous flow for a given duration.
 *
 * Useful for analyzing a specific mode without transition logic.
 *
 * @param ha       Automaton (for flow definition)
 * @param mode_id  Mode to simulate
 * @param x        Initial state (updated in-place)
 * @param t_start  Start time
 * @param duration How long to simulate
 * @param dt       Time step
 * @param method   Integration method
 * @param n        State dimension
 * @return         Number of integration steps taken
 */
int  hybrid_simulate_mode_flow(const HybridAutomaton *ha, int mode_id,
                                double *x, double t_start, double duration,
                                double dt, HybridIntegratorType method, int n);

/**
 * @brief Check if an invariant boundary is approached (predictive check).
 *
 * Returns the estimated time until the invariant boundary is reached,
 * based on the current velocity and distance to the boundary.
 *
 * @param invariant The invariant region
 * @param x         Current state
 * @param dxdt      Current derivative
 * @param n         State dimension
 * @return          Estimated time to boundary, or INFINITY if receding
 */
double hybrid_time_to_boundary(const HybridInvariant *invariant,
                                const double *x, const double *dxdt, int n);

#endif /* HYBRID_SIMULATION_H */
