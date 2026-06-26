/**
 * @file pwa_simulation.h
 * @brief PWA System Simulation Engine — L5 Algorithms
 *
 * Simulates piecewise affine systems in both continuous-time (CT)
 * and discrete-time (DT). The main challenge in CT simulation is
 * event detection: precisely locating when the state trajectory
 * crosses a region boundary.
 *
 * For DT-PWA systems, simulation is straightforward: at each step
 * apply the active region's affine map, then re-evaluate point
 * location to determine the region for the next step.
 *
 * For CT-PWA systems, we use numerical ODE integration (RK4) with
 * event detection via root-finding on the boundary constraint
 * functions. When a boundary crossing is detected, we use bisection
 * or the secant method to locate the exact crossing time.
 *
 * References:
 *   Hairer, E., Nørsett, S. P., & Wanner, G. (1993).
 *     "Solving Ordinary Differential Equations I: Nonstiff Problems."
 *     Springer-Verlag. Chapter II.6 (Event Location).
 *   Cellier, F. E. & Kofman, E. (2006). "Continuous System Simulation."
 *     Springer. Chapter 8 (Hybrid Systems).
 *   Mosterman, P. J. (1999). "An overview of hybrid simulation
 *     phenomena and their support by simulation packages."
 *     HSCC '99, LNCS 1569:165-177.
 *
 * Knowledge coverage:
 *   L5: Event detection, RK4 integration, Zeno detection,
 *       bisection root-finding, PWA trajectory simulation
 *
 * Nine-school course alignment:
 *   MIT 6.241J — Lec 21 Numerical simulation
 *   Stanford AA203 — Lec 5 ODE solvers
 *   Berkeley EE221A — Lec 16 Simulation
 *   CMU 24-677 — Lec 20 Hybrid simulation
 *   ETH 227-0216 — Lec 18 Numerical methods
 */

#ifndef PWA_SIMULATION_H
#define PWA_SIMULATION_H

#include "pwa_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
 * L5: Simulation Parameters and Configuration
 *===========================================================================*/

/**
 * @brief Simulation configuration.
 */
typedef struct {
    double  t_start;        /**< Start time */
    double  t_end;          /**< End time */
    double  dt_max;         /**< Maximum integration step (CT) or fixed step (DT) */
    double  dt_min;         /**< Minimum step size (prevents Zeno lock) */
    int     max_steps;      /**< Maximum number of steps */
    double  event_tol;      /**< Tolerance for event detection (boundary crossing) */
    int     max_events;     /**< Maximum number of switch events */
    int     store_all;      /**< 1 = store all steps, 0 = store final only */
    int     verbose;        /**< Verbosity level (0=quiet, 1=summary, 2=detailed) */
} PWASimConfig;

/**
 * @brief ODE integration method for CT-PWA simulation.
 */
typedef enum {
    PWA_ODE_EULER       = 0,  /**< Forward Euler (1st order) */
    PWA_ODE_HEUN        = 1,  /**< Heun's method (2nd order, explicit trapezoidal) */
    PWA_ODE_RK4         = 2,  /**< Classical Runge-Kutta 4th order */
    PWA_ODE_RKF45       = 3   /**< Runge-Kutta-Fehlberg 4(5) adaptive */
} PWAODEMethod;

/*===========================================================================
 * L5: ODE Integration Functions
 *===========================================================================*/

/**
 * @brief Evaluate the PWA vector field at a point (x,u) in a region.
 *
 * f(x,u) = A_i * x + B_i * u + f_i   for [x;u] ∈ region i
 *
 * @param sys    PWA system
 * @param region Region index
 * @param x      State vector (n_state)
 * @param u      Input vector (n_input)
 * @param dxdt   Output: time derivative (n_state)
 * @return 0 on success, -1 if region invalid
 *
 * Complexity: O(n_state^2 + n_state * n_input)
 */
int pwa_vector_field(const PWASystem *sys, int region,
                      const double *x, const double *u, double *dxdt);

/**
 * @brief Single RK4 integration step for PWA system.
 *
 * k1 = h · f(t_n, x_n)
 * k2 = h · f(t_n + h/2, x_n + k1/2)
 * k3 = h · f(t_n + h/2, x_n + k2/2)
 * k4 = h · f(t_n + h, x_n + k3)
 * x_{n+1} = x_n + (k1 + 2k2 + 2k3 + k4) / 6
 *
 * @param sys    PWA system
 * @param region Current region
 * @param t      Current time
 * @param x      Current state (n_state)
 * @param u      Current input (n_input)
 * @param h      Step size
 * @param x_next Output: next state (n_state)
 * @return 0 on success, -1 on error
 *
 * Complexity: O(n_state^2 + n_state * n_input)
 * Local truncation error: O(h^5)
 */
int pwa_rk4_step(const PWASystem *sys, int region,
                  double t, const double *x, const double *u,
                  double h, double *x_next);

/**
 * @brief Adaptive step Runge-Kutta-Fehlberg 4(5) step.
 *
 * Estimates error and adjusts step size to meet tolerance.
 *
 * @param sys       PWA system
 * @param region    Current region
 * @param t         Current time
 * @param x         Current state (n_state)
 * @param u         Current input (n_input)
 * @param h         Input/output: step size (adjusted)
 * @param x_next    Output: next state (n_state)
 * @param err       Output: estimated local error
 * @param tol       Desired tolerance
 * @return 0 on success, -1 on error
 */
int pwa_rkf45_step(const PWASystem *sys, int region,
                    double t, const double *x, const double *u,
                    double *h, double *x_next, double *err, double tol);

/*===========================================================================
 * L5: Event Detection (Boundary Crossing)
 *===========================================================================*/

/**
 * @brief Check if a state has crossed any region boundary.
 *
 * Evaluates all constraint functions g_j(x) = H_j·z - K_j
 * for the current region. A crossing is detected when any
 * g_j changes sign between two consecutive states.
 *
 * @param sys    PWA system
 * @param region Current region
 * @param x0     Previous state
 * @param u0     Previous input
 * @param x1     Current state
 * @param u1     Current input
 * @param crossed_halfspace Output: index of crossed half-space, -1 if none
 * @return 1 if boundary crossing detected, 0 otherwise
 */
int pwa_detect_boundary_crossing(const PWASystem *sys, int region,
                                  const double *x0, const double *u0,
                                  const double *x1, const double *u1,
                                  int *crossed_halfspace);

/**
 * @brief Find the exact time of boundary crossing using bisection.
 *
 * Given that a crossing occurred between states x(t_a) and x(t_b),
 * find t* such that the constraint g_j(x(t*)) = 0.
 *
 * @param sys      PWA system
 * @param region   Current region
 * @param ta       Start time
 * @param xa       State at start time
 * @param u        Input (constant over interval)
 * @param tb       End time
 * @param xb       State at end time
 * @param halfspace Half-space index that was crossed
 * @param t_cross  Output: exact crossing time
 * @param x_cross  Output: state at crossing time
 * @param tol      Tolerance for root finding
 * @return 0 on success, -1 on failure
 *
 * Complexity: O(log((tb-ta)/tol) · n_state^2)
 */
int pwa_find_crossing_time(const PWASystem *sys, int region,
                            double ta, const double *xa,
                            const double *u,
                            double tb, const double *xb,
                            int halfspace,
                            double *t_cross, double *x_cross, double tol);

/*===========================================================================
 * L5: PWA Simulation Engine
 *===========================================================================*/

/**
 * @brief Simulate a discrete-time PWA system.
 *
 * Iterates: find region, apply affine map, repeat.
 * Detects region switches at each step.
 *
 * Stores the complete trajectory including state, input, output,
 * region, and event history.
 *
 * @param sys           PWA system
 * @param x0            Initial state (n_state)
 * @param u_seq         Input sequence (n_steps × n_input) or NULL for u=0
 * @param config        Simulation configuration
 * @param traj          Output: trajectory (caller allocates via pwa_trajectory_create)
 * @return 0 on success, -1 on error
 *
 * Complexity: O(n_steps · n_regions · n_cons · (n_state+n_input))
 *
 * Reference: Bemporad (2004), "Hybrid Toolbox - User's Guide."
 */
int pwa_simulate_dt(PWASystem *sys,
                     const double *x0, const double *u_seq,
                     const PWASimConfig *config, PWATrajectory *traj);

/**
 * @brief Simulate a continuous-time PWA system with event detection.
 *
 * Uses adaptive-step RK4(5) integration with boundary crossing
 * detection and bisection-based event location.
 *
 * Handles mode switches by detecting boundary crossings,
 * locating exact switching times, and resuming integration
 * in the new region.
 *
 * @param sys           PWA system
 * @param x0            Initial state (n_state)
 * @param u_func        Input as function of time (NULL for zero input)
 * @param u_func_ctx    Context for input function
 * @param config        Simulation configuration
 * @param traj          Output: trajectory
 * @return 0 on success, -1 on error
 */
int pwa_simulate_ct(PWASystem *sys,
                     const double *x0,
                     double (*u_func)(double t, int dim, void *ctx),
                     void *u_func_ctx,
                     const PWASimConfig *config, PWATrajectory *traj);

/**
 * @brief Default simulation configuration.
 *
 * @return PWASimConfig with safe defaults
 */
PWASimConfig pwa_sim_config_default(void);

/*===========================================================================
 * L5: Zeno Detection
 *===========================================================================*/

/**
 * @brief Check for Zeno behavior in a simulation trajectory.
 *
 * Zeno behavior: infinite number of discrete transitions in finite time.
 * Detected when switching times accumulate with sum converging to
 * a finite limit.
 *
 * @param traj Simulation trajectory
 * @param tb   Output: Zeno accumulation time (if detected)
 * @return 1 if Zeno detected, 0 otherwise
 *
 * Reference: Zhang, J., Johansson, K. H., Lygeros, J., & Sastry, S. (2001).
 *   "Zeno hybrid systems." Int. J. Robust Nonlinear Control, 11(5):435-451.
 */
int pwa_detect_zeno(const PWATrajectory *traj, double *tb);

/**
 * @brief Regularize a PWA system to prevent Zeno behavior.
 *
 * Adds hysteresis or minimum dwell time constraints to prevent
 * chattering at region boundaries.
 *
 * @param sys          PWA system
 * @param min_dwell    Minimum dwell time in each region
 * @param hysteresis   Hysteresis band width
 * @return 0 on success, -1 on error
 */
int pwa_regularize_zeno(PWASystem *sys, double min_dwell, double hysteresis);

/*===========================================================================
 * L5: Trajectory Operations
 *===========================================================================*/

/**
 * @brief Create a trajectory structure for PWA simulation output.
 *
 * @param n_state  State dimension
 * @param n_input  Input dimension
 * @param n_output Output dimension
 * @param n_max    Maximum number of time steps
 * @return New trajectory, or NULL on error
 */
PWATrajectory* pwa_trajectory_create(int n_state, int n_input,
                                      int n_output, int n_max);

/**
 * @brief Destroy a trajectory and free all memory.
 *
 * @param traj Trajectory to destroy
 */
void pwa_trajectory_destroy(PWATrajectory *traj);

/**
 * @brief Print trajectory summary to stdout.
 *
 * @param traj Trajectory to print
 */
void pwa_trajectory_print(const PWATrajectory *traj);

/**
 * @brief Export trajectory to CSV format.
 *
 * @param traj      Trajectory to export
 * @param filename  Output filename
 * @return 0 on success, -1 on error
 */
int pwa_trajectory_export_csv(const PWATrajectory *traj, const char *filename);

#ifdef __cplusplus
}
#endif

#endif /* PWA_SIMULATION_H */
