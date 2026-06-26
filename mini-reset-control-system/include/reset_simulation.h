/* reset_simulation.h - Hybrid Simulation Engine with Reset Detection
 *
 * Implements the numerical simulation of reset control systems as
 * hybrid dynamical systems with continuous flow and discrete jumps.
 *
 * Key challenge: detecting zero-crossings accurately within a
 * numerical integration framework, since the reset instant t_k
 * is generally not aligned with integration mesh points.
 *
 * Knowledge coverage:
 *   L2: Hybrid time domain (t,j), Zeno detection, dwell time
 *   L3: Hybrid automaton simulation, event detection
 *   L5: Zero-crossing bracketing, bisection, dense output
 *   L8: Variable-step integration with event detection
 *
 * Ref: [GST12] Ch.4 on hybrid simulation;
 *      Shampine & Thompson (2000), "Event location for ODEs",
 *      Computers & Mathematics with Applications
 */

#ifndef RESET_SIMULATION_H
#define RESET_SIMULATION_H

#include "reset_core.h"
#include "reset_system.h"
#include <stdbool.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L3: Simulation Configuration and State
 * ================================================================ */

/** SimMethod - numerical integration method for continuous flow. */
typedef enum {
    SIM_METHOD_EULER        = 0,  /* Forward Euler (1st order)        */
    SIM_METHOD_HEUN         = 1,  /* Heun's method (RK2, 2nd order)  */
    SIM_METHOD_RK4          = 2,  /* Classical Runge-Kutta (4th order)*/
    SIM_METHOD_RKF45        = 3   /* Runge-Kutta-Fehlberg (adaptive)  */
} SimMethod;

/** SimConfig - configuration for hybrid simulation. */
typedef struct {
    SimMethod    method;          /**< integration method             */
    double       dt_min;          /**< minimum time step (adaptive)   */
    double       dt_max;          /**< maximum time step (adaptive)   */
    double       dt_fixed;        /**< fixed step size                */
    double       t_start;         /**< simulation start time          */
    double       t_end;           /**< simulation end time            */
    double       zc_tol;          /**< zero-crossing detection tolerance */
    int          max_resets;      /**< maximum resets before abort    */
    int          max_steps;       /**< maximum integration steps      */
    double       zc_refine_tol;   /**< tolerance for ZC refinement    */
    int          zc_max_iter;     /**< max iterations for ZC bisection */
    bool         adaptive;        /**< use adaptive step size?        */
    bool         detect_zeno;     /**< check for Zeno behavior?       */
    double       zeno_threshold;  /**< min interval for Zeno warning  */
} SimConfig;

/** SimState - runtime state of the hybrid simulation. */
typedef struct {
    ResetSystem     *rsys;         /**< the reset control system      */
    double           t;            /**< current simulation time        */
    double          *x;            /**< current continuous state       */
    double           u;            /**< current input (scalar for SISO)*/
    double          *u_vec;        /**< input vector (MIMO)            */
    int              step_count;   /**< number of integration steps    */
    int              reset_count;  /**< number of reset events         */
    double           t_last_reset; /**< time of last reset             */
    double           dt_current;   /**< current step size (adaptive)   */
    bool             in_flow;      /**< true if in continuous flow     */
    bool             zeno_warning; /**< Zeno behavior suspected        */
    ResetMemory     *reset_log;    /**< log of reset events            */
    int              log_capacity; /**< capacity of reset_log          */
    int              log_count;    /**< number of entries in log       */
} SimState;

/** SimResult - outcome of a complete simulation run. */
typedef struct {
    int        n_steps;         /**< total integration steps taken     */
    int        n_resets;        /**< total reset events occurred       */
    double     cpu_time_ms;     /**< approximate CPU time in ms       */
    double     t_final;         /**< final simulation time             */
    bool       completed;       /**< true if reached t_end             */
    bool       zeno_detected;   /**< true if Zeno behavior detected   */
    double    *t_history;       /**< time history [n_steps+1]          */
    double    *x_history;       /**< state history [n_steps+1 * nx]   */
    double    *y_history;       /**< output history [n_steps+1 * ny]  */
    int        history_len;     /**< length of history arrays         */
} SimResult;

/* ================================================================
 * L5: Simulation Engine Functions
 * ================================================================ */

/** Create a default simulation configuration.
 *  Sets: RK4 method, dt_fixed=0.001, t_end=10.0, max_resets=10000. */
SimConfig sim_config_default(void);

/** Create a simulation state for a given reset system. */
SimState* sim_state_create(ResetSystem *rsys);

/** Free a simulation state and all allocated memory. */
void sim_state_free(SimState *state);

/** Set the initial continuous state of the simulation. */
void sim_set_initial_state(SimState *state, const double *x0);

/** Perform one integration step of the continuous flow dynamics.
 *  Uses the configured integration method (Euler, Heun, RK4, RKF45).
 *  Does NOT check for reset conditions - call sim_check_reset() after.
 *  Complexity: O(n^2) to O(n^3) depending on method.
 *  Returns the step size actually used. */
double sim_integrate_flow(SimState *state, const SimConfig *cfg,
                           double u);

/** Check for zero-crossing between the last two states.
 *  If a ZC is detected within tolerance zc_tol:
 *    - Uses bisection to refine the crossing time
 *    - Applies the reset jump map
 *    - Resumes flow integration from the after-reset state
 *  Returns true if a reset was triggered and applied. */
bool sim_check_reset(SimState *state, const SimConfig *cfg,
                      double e_now, double e_prev);

/** Run a complete hybrid simulation from t_start to t_end.
 *  Handles the loop: integrate flow -> check reset -> integrate flow.
 *  Stops early if max_resets or max_steps is exceeded.
 *  Returns a SimResult with simulation history.
 *  Caller must call sim_result_free() to release history arrays. */
SimResult* sim_run(SimState *state, const SimConfig *cfg,
                    const double *u_vals, int n_u);

/** Run simulation with error-based feedback: u = -(Cx + Du) equivalent.
 *  r_vals: reference input at each time point.
 *  n_r: length of reference array (also determines number of steps). */
SimResult* sim_run_feedback(SimState *state, const SimConfig *cfg,
                             const double *r_vals, int n_r);

/** Free a simulation result and all history arrays. */
void sim_result_free(SimResult *result);

/* ---- Zero-Crossing Detection (L5 Algorithm) ---- */

/** Detect a zero crossing between (t0, x0) and (t1, x1).
 *  Uses linear interpolation to estimate ZC time.
 *  zc_time: output, estimated ZC time.
 *  Returns true if ZC detected. */
bool sim_detect_zc(double e0, double e1, double t0, double t1,
                    double *zc_time);

/** Refine a zero-crossing time using the secant method.
 *  f_eval: callback that evaluates the error signal at time t.
 *  t_a, t_b: bracketing interval (f(t_a)*f(t_b) < 0).
 *  tol: convergence tolerance. max_iter: maximum iterations.
 *  Returns the refined ZC time. */
typedef double (*sim_error_func)(double t, void *ctx);

double sim_refine_zc(sim_error_func f_eval, void *ctx,
                      double t_a, double t_b, double tol, int max_iter);

/** Count resets in a time window to detect Zeno behavior.
 *  If (n_resets_in_window) / (window_size) > zeno_threshold,
 *  the system likely exhibits Zeno behavior.
 *  Returns estimated Zeno time if detected, INFINITY otherwise. */
double sim_detect_zeno(const SimState *state, double window_size,
                        double zeno_threshold);

/* ---- Output and Logging ---- */

/** Print the simulation history to stdout in CSV format.
 *  Useful for piping to plotting tools. */
void sim_print_history(const SimResult *result, FILE *fp);

/** Get the reset event log. Returns read-only pointer. */
const ResetMemory* sim_get_reset_log(const SimState *state, int *count);

/** Compute RMS error between output history and reference. */
double sim_compute_rms_error(const SimResult *result,
                              const double *ref, int n_ref);

/** Compute settling time (time to stay within tol% of final value). */
double sim_settling_time(const SimResult *result, double tol_pct);

#ifdef __cplusplus
}
#endif

#endif /* RESET_SIMULATION_H */