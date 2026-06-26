/**
 * @file pwa_defs.h
 * @brief Piecewise Affine (PWA) System — L1 Core Definitions
 *
 * A Piecewise Affine (PWA) system is a hybrid dynamical system where
 * the state-input space is partitioned into a finite set of polyhedral
 * regions, and within each region the dynamics are affine (linear + offset).
 *
 * Discrete-time PWA (DT-PWA):
 *   x(k+1) = A_i x(k) + B_i u(k) + f_i   for  [x(k); u(k)] ∈ R_i
 *   y(k)   = C_i x(k) + D_i u(k) + g_i   for  [x(k); u(k)] ∈ R_i
 *
 * Continuous-time PWA (CT-PWA):
 *   dx/dt  = A_i x(t) + B_i u(t) + f_i   for  [x(t); u(t)] ∈ R_i
 *   y(t)   = C_i x(t) + D_i u(t) + g_i   for  [x(t); u(t)] ∈ R_i
 *
 * Each region R_i is a convex polyhedron:
 *   R_i = { z ∈ R^{n+m} | H_i z ≤ K_i }
 *
 * References:
 *   Sontag, E. D. (1981). "Nonlinear Regulation: The Piecewise Linear Approach."
 *     IEEE TAC, 26(2):346-358.
 *   Johansson, M. (2003). "Piecewise Linear Control Systems."
 *     Springer-Verlag. Lecture Notes in Control and Information Sciences.
 *   Bemporad, A. & Morari, M. (1999). "Control of systems integrating logic,
 *     dynamics, and constraints." Automatica, 35(3):407-427.
 *
 * Knowledge coverage:
 *   L1: PWA system definition, region types, affine dynamics, mode types
 *   L2: state-space partitioning, mode switching, well-posedness
 *
 * Nine-school course alignment:
 *   MIT 6.241J Dynamic Systems — Lec 20 Hybrid/PWA systems
 *   Stanford AA203 Optimal Control — Lec 15 Hybrid MPC
 *   Berkeley EE222 Nonlinear Systems — Ch 7 PWA systems
 *   CMU 24-677 Nonlinear Ctrl — Lec 18 PWA/MLD
 *   Princeton MAE 546 Optimal Ctrl — Lec 19 Hybrid systems
 *   Caltech CDS140 Nonlinear Dynamics — Lec 16 PWA models
 *   Cambridge 4F3 Nonlinear Ctrl — Lec 8 Switched/PWA
 *   Oxford B4 Predictive Control — Lec 9 Hybrid MPC
 *   ETH 227-0220 Model Reduction — Lec 12 PWA approximations
 */

#ifndef PWA_DEFS_H
#define PWA_DEFS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
 * L1: Core PWA Type Definitions
 *===========================================================================*/

/**
 * @brief Affine dynamics within a single region.
 *
 * Represents: dx/dt = A*x + B*u + f  (CT) or x^+ = A*x + B*u + f (DT)
 *           y = C*x + D*u + g
 *
 * A ∈ R^{n×n}, B ∈ R^{n×m}, f ∈ R^n
 * C ∈ R^{p×n}, D ∈ R^{p×m}, g ∈ R^p
 */
typedef struct {
    int     n_state;      /**< State dimension n */
    int     n_input;      /**< Input dimension m */
    int     n_output;     /**< Output dimension p */
    double *A;            /**< State matrix, n_state × n_state */
    double *B;            /**< Input matrix, n_state × n_input */
    double *f;            /**< Affine offset, n_state */
    double *C;            /**< Output matrix, n_output × n_state */
    double *D;            /**< Feedthrough matrix, n_output × n_input */
    double *g;            /**< Output offset, n_output */
} PWAAffineDynamics;

/**
 * @brief A single polyhedral region of a PWA partition.
 *
 * R_i = { x ∈ R^n, u ∈ R^m | H_i * [x;u] ≤ K_i }
 *
 * H_i ∈ R^{n_constraints × (n+m)}, K_i ∈ R^{n_constraints}
 */
typedef struct {
    int     id;           /**< Region identifier */
    int     n_state;      /**< State dimension n */
    int     n_input;      /**< Input dimension m */
    int     n_constraints; /**< Number of half-space constraints */
    double *H;            /**< Constraint matrix, n_constraints × (n+m) */
    double *K;            /**< Constraint right-hand side, n_constraints */
    int     dynamics_id;  /**< Index into PWA system dynamics array */
    int     is_active;    /**< Whether this region is active */
} PWARegion;

/**
 * @brief Complete Piecewise Affine (PWA) System.
 *
 * A PWA system is defined by a set of regions that partition the
 * state-input space, each with its own affine dynamics.
 *
 * Mode switches occur when the state-input vector crosses region
 * boundaries. The system is well-posed if the regions form a partition
 * (disjoint interiors, union covers the domain of interest).
 */
typedef struct {
    int              n_state;        /**< State dimension */
    int              n_input;        /**< Input dimension */
    int              n_output;       /**< Output dimension */
    int              n_regions;      /**< Number of regions */
    int              n_allocated;    /**< Allocated capacity for regions */
    PWARegion       *regions;        /**< Array of n_regions regions */
    PWAAffineDynamics *dynamics;     /**< Array of n_regions dynamics */
    int              is_continuous;  /**< 1 = continuous-time, 0 = discrete-time */
    double           dt;             /**< Sampling time for discrete-time (0 for CT) */
    int              current_region; /**< Current active region (-1 if unknown) */
    double          *current_state;  /**< Current state vector (n_state) */
    double          *current_input;  /**< Current input vector (n_input) */
    double          *x_min;          /**< State lower bounds (n_state) */
    double          *x_max;          /**< State upper bounds (n_state) */
    double          *u_min;          /**< Input lower bounds (n_input) */
    double          *u_max;          /**< Input upper bounds (n_input) */
} PWASystem;

/**
 * @brief Switching type classification.
 *
 * PWA systems can exhibit different types of switches:
 * - EXOGENOUS: External signal triggers switch
 * - AUTONOMOUS: State crossing boundary triggers switch
 * - CONTROLLED: Control input choice triggers switch
 */
typedef enum {
    PWA_SWITCH_EXOGENOUS  = 0,  /**< Externally triggered */
    PWA_SWITCH_AUTONOMOUS = 1,  /**< State-dependent */
    PWA_SWITCH_CONTROLLED = 2   /**< Control-dependent */
} PWASwitchType;

/**
 * @brief Event record for PWA simulation.
 *
 * Captures mode switches and boundary crossing events.
 */
typedef struct {
    double          t;              /**< Time of event */
    int             from_region;    /**< Previous region index */
    int             to_region;      /**< New region index */
    PWASwitchType   type;           /**< Type of switch */
    double         *x_at_event;     /**< State at event (n_state) */
    double         *boundary_dist;  /**< Signed distance to each boundary */
} PWAEvent;

/**
 * @brief PWA simulation trajectory.
 *
 * Stores the complete state, input, output, and event history
 * of a PWA system simulation.
 */
typedef struct {
    int             n_steps;        /**< Number of time steps stored */
    int             n_state;        /**< State dimension */
    int             n_input;        /**< Input dimension */
    int             n_output;       /**< Output dimension */
    int             n_max;          /**< Maximum allocated steps */
    double         *t_hist;         /**< Time history (n_max) */
    double         *x_hist;         /**< State history, row-major (n_max × n_state) */
    double         *u_hist;         /**< Input history (n_max × n_input) */
    double         *y_hist;         /**< Output history (n_max × n_output) */
    int            *region_hist;    /**< Active region history (n_max) */
    int             n_events;       /**< Number of events */
    PWAEvent       *events;         /**< Event list */
} PWATrajectory;

/**
 * @brief PWA Mode (discrete state) definition.
 *
 * A mode is a discrete state with associated continuous dynamics.
 * In the hybrid automaton view, modes are locations and region
 * boundaries are guards.
 */
typedef struct {
    int             mode_id;        /**< Mode identifier */
    int             dynamics_id;    /**< Index into dynamics array */
    int             n_invariants;   /**< Number of invariant constraints */
    double         *H_inv;          /**< Invariant: H_inv * [x;u] ≤ K_inv */
    double         *K_inv;          /**< Invariant RHS */
    int             n_transitions;  /**< Number of outgoing transitions */
    int            *target_modes;   /**< Target mode for each transition */
    int            *guard_regions;  /**< Guard region index for each transition */
} PWAMode;

/**
 * @brief PWA Hybrid Automaton representation.
 *
 * A hybrid automaton with modes (locations) and transitions (edges).
 * Equivalent representation to the region-based PWA system.
 */
typedef struct {
    int             n_state;        /**< State dimension */
    int             n_input;        /**< Input dimension */
    int             n_output;       /**< Output dimension */
    int             n_modes;        /**< Number of discrete modes */
    PWAMode        *modes;          /**< Array of modes */
    PWAAffineDynamics *dynamics;    /**< Shared dynamics array */
    int             initial_mode;   /**< Initial mode */
} PWAHybridAutomaton;

/**
 * @brief Well-posedness status for a PWA system.
 *
 * A PWA system is well-posed if:
 * 1. Regions have disjoint interiors
 * 2. Union of regions covers the domain of interest
 * 3. No Zeno behavior (in CT) - infinite switches in finite time
 */
typedef enum {
    PWA_WELLPOSED_OK       = 0,  /**< System is well-posed */
    PWA_WELLPOSED_OVERLAP  = 1,  /**< Regions overlap */
    PWA_WELLPOSED_GAPS     = 2,  /**< Regions have gaps */
    PWA_WELLPOSED_ZENO     = 3,  /**< Zeno behavior detected */
    PWA_WELLPOSED_INVALID  = 4   /**< Invalid definition */
} PWAWellPosedStatus;

/*===========================================================================
 * L1+L2: Core PWA Operations
 *===========================================================================*/

/**
 * @brief Create a new PWA system.
 *
 * Allocates and initializes a PWASystem with given dimensions.
 * Regions and dynamics arrays are pre-allocated to n_regions_max.
 *
 * @param n_state       State dimension
 * @param n_input       Input dimension
 * @param n_output      Output dimension
 * @param n_regions_max Maximum number of regions
 * @param is_continuous 1 for continuous-time, 0 for discrete-time
 * @param dt            Sampling time (0 for continuous-time)
 * @return Pointer to new PWASystem, or NULL on allocation failure
 *
 * Complexity: O(n_state^2 + n_regions_max * (n_state * (n_state + n_input)))
 */
PWASystem* pwa_system_create(int n_state, int n_input, int n_output,
                              int n_regions_max, int is_continuous, double dt);

/**
 * @brief Destroy a PWA system and free all allocated memory.
 *
 * @param sys PWA system to destroy
 */
void pwa_system_destroy(PWASystem *sys);

/**
 * @brief Add an affine dynamics to a PWA system.
 *
 * @param sys   PWA system
 * @param A     State matrix (n×n, row-major)
 * @param B     Input matrix (n×m, row-major)
 * @param f     Affine offset vector (n)
 * @param C     Output matrix (p×n, row-major)
 * @param D     Feedthrough matrix (p×m, row-major)
 * @param g     Output offset vector (p)
 * @return Dynamics index, or -1 on error
 */
int pwa_add_dynamics(PWASystem *sys,
                      const double *A, const double *B, const double *f,
                      const double *C, const double *D, const double *g);

/**
 * @brief Add a polyhedral region to a PWA system.
 *
 * The region is defined by H*z ≤ K where z = [x; u].
 *
 * @param sys           PWA system
 * @param H             Constraint matrix (n_cons × (n+m))
 * @param K             Constraint RHS (n_cons)
 * @param n_cons        Number of constraints
 * @param dynamics_id   Index of dynamics for this region
 * @return Region index, or -1 on error
 */
int pwa_add_region(PWASystem *sys,
                    const double *H, const double *K, int n_cons,
                    int dynamics_id);

/**
 * @brief Validate PWA system for consistency.
 *
 * Checks:
 * - All regions referenced dynamics exist
 * - Region dimensions match system dimensions
 * - Partition covers domain (no gaps in bounding box)
 *
 * @param sys PWA system to validate
 * @return 0 on success, negative on error
 */
int pwa_system_validate(const PWASystem *sys);

/**
 * @brief Check if a point [x;u] lies within a given region.
 *
 * Evaluates all half-space constraints H*z ≤ K.
 *
 * @param region Region to test
 * @param z      Point vector [x; u] of length n_state + n_input
 * @return 1 if z ∈ region, 0 otherwise
 *
 * Theorem (Polyhedron Membership): z ∈ R_i ⇔ H_i z ≤ K_i (component-wise)
 * Complexity: O(n_constraints * (n_state + n_input))
 */
int pwa_point_in_region(const PWARegion *region, const double *z);

/**
 * @brief Point location: find which region contains [x;u].
 *
 * Searches through all regions to find the active one.
 * If multiple regions contain the point (overlap), returns the first.
 *
 * @param sys PWA system
 * @param x   State vector (n_state)
 * @param u   Input vector (n_input)
 * @return Region index, or -1 if no region contains the point
 *
 * Complexity: O(n_regions * n_constraints_max * (n_state + n_input))
 */
int pwa_point_location(const PWASystem *sys, const double *x, const double *u);

/**
 * @brief Check well-posedness of a PWA system.
 *
 * Tests for overlapping regions and coverage gaps in the
 * specified domain [x_min, x_max] × [u_min, u_max].
 *
 * @param sys PWA system to check
 * @return PWAWellPosedStatus
 *
 * Reference: Johansson (2003), Definition 2.1
 */
PWAWellPosedStatus pwa_check_wellposed(const PWASystem *sys);

/**
 * @brief Compute region adjacency: find which regions share a boundary.
 *
 * Two regions are adjacent if their boundaries intersect in a
 * (n+m-1)-dimensional facet.
 *
 * @param sys          PWA system
 * @param adjacency    Output: adjacency matrix (n_regions × n_regions),
 *                     adjacency[i * n_regions + j] = 1 if regions i,j adjacent
 * @return 0 on success, -1 on error
 *
 * Complexity: O(n_regions^2 * n_cons^2)
 */
int pwa_compute_adjacency(const PWASystem *sys, int *adjacency);

/**
 * @brief Print a PWA system description to stdout (for debugging).
 *
 * @param sys PWA system to print
 * @param verbose 0=summary, 1=full details
 */
void pwa_system_print(const PWASystem *sys, int verbose);

#ifdef __cplusplus
}
#endif

#endif /* PWA_DEFS_H */
