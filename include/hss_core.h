/**
 * @file hss_core.h
 * @brief Hybrid Switched System (HSS) — Unified Core Definitions (L1-L2)
 *
 * This header provides the unified theoretical framework for hybrid switched
 * systems, integrating concepts from:
 *
 *   - Hybrid Automata (Alur/Henzinger 1993–1996)
 *   - Switched Systems Stability (Liberzon 2003; Morse 1996)
 *   - Dwell-Time Analysis (Hespanha & Morse 1999)
 *   - Event-Triggered Control (Tabuada 2007)
 *   - Impulsive Systems (Lakshmikantham/Bainov/Simeonov 1989)
 *   - Piecewise Affine Systems (Sontag 1981; Johansson 2003)
 *   - Reset Control (Clegg 1958; Horowitz/Rosenbaum 1975; Beker/Hollot/Chait 2004)
 *   - Supervisory Control (Ramadge & Wonham 1987)
 *
 * A hybrid switched system is formally defined as:
 *   H = (Q, X, Σ, F, Init, Inv, E, G, R)
 * where:
 *   Q   = finite set of discrete modes (locations)
 *   X   = continuous state space (R^n)
 *   Σ   = set of switching signals σ: R≥0 → Q
 *   F   = family of vector fields {f_q : X → TX}
 *   Init = initial condition set ⊆ Q × X
 *   Inv  = invariant region for each mode
 *   E    = set of discrete transitions (edges)
 *   G    = guard conditions per edge
 *   R    = reset maps per edge
 *
 * ============================================================================
 * Knowledge Coverage (L1-L2)
 * ============================================================================
 *
 * L1 — Core Definitions:
 *   KP1: HSS_System — unified hybrid switched system struct
 *   KP2: HSS_Mode — discrete operating mode with continuous dynamics
 *   KP3: HSS_SwitchingSignal — piecewise-constant signal σ(t)
 *   KP4: HSS_Transition — guarded discrete transition between modes
 *   KP5: HSS_ContinuousState — continuous state vector x ∈ R^n
 *   KP6: HSS_SystemClass — system classification (linear/affine/nonlinear)
 *   KP7: HSS_StabilityConcept — stability definitions for hybrid systems
 *   KP8: HSS_SolverType — numerical integration method taxonomy
 *   KP9: HSS_Event — discrete event triggering transition
 *   KP10: HSS_ResetMap — state discontinuity at switching instants
 *
 * L2 — Core Concepts:
 *   KC1: Hybrid time domain — (t, j) ∈ R≥0 × N
 *   KC2: Flow set and jump set partitioning
 *   KC3: Zeno behavior detection and exclusion
 *   KC4: Dwell-time constrained switching
 *   KC5: Average dwell-time switching
 *   KC6: State-dependent switching logic
 *   KC7: Common vs. Multiple Lyapunov Functions
 *   KC8: Input-to-State Stability (ISS) for hybrid systems
 *
 * Course Mapping:
 *   MIT 6.241/6.824 — Hybrid & Distributed Systems
 *   Stanford CS359 — Hybrid Systems
 *   Berkeley EECS 291E — Hybrid Systems
 *   CMU 15-424 — Foundations of Cyber-Physical Systems
 *   ETH 227-0690 — Hybrid System Theory
 *   Cambridge Part II — Control Theory
 *   Georgia Tech CS 6290 — HPCA + advanced control
 *
 * References:
 *   [1] Liberzon, D. "Switching in Systems and Control." Birkhäuser, 2003.
 *   [2] Goebel, R., Sanfelice, R.G., Teel, A.R. "Hybrid Dynamical Systems."
 *       Princeton, 2012.
 *   [3] Tabuada, P. "Verification and Control of Hybrid Systems." Springer, 2009.
 *   [4] van der Schaft, A., Schumacher, J.M. "An Introduction to Hybrid
 *       Dynamical Systems." Springer, 2000.
 *   [5] Alur, R. "Principles of Cyber-Physical Systems." MIT Press, 2015.
 */

#ifndef HSS_CORE_H
#define HSS_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Maximum number of discrete modes in a hybrid switched system */
#define HSS_MAX_MODES            128

/** Maximum number of continuous state variables */
#define HSS_MAX_STATE_DIM         64

/** Maximum number of control inputs */
#define HSS_MAX_INPUT_DIM         32

/** Maximum number of output variables */
#define HSS_MAX_OUTPUT_DIM        32

/** Maximum number of discrete transitions between modes */
#define HSS_MAX_TRANSITIONS       512

/** Maximum number of guards per transition */
#define HSS_MAX_GUARDS            16

/** Maximum switching signal history length */
#define HSS_MAX_SWITCH_HISTORY   2048

/** Maximum length of mode/state/variable names */
#define HSS_NAME_LEN              64

/** Convergence tolerance for numerical methods */
#define HSS_EPSILON              1e-9

/** Maximum iterations for iterative algorithms */
#define HSS_MAX_ITER             5000

/* ============================================================================
 * L1 KP1-KP3: System Classification Enums
 * ============================================================================ */

/**
 * @brief System dynamics class (L1 KP6)
 *
 * Determines the form of the vector field f_q(x,u):
 *   LINEAR:     ẋ = A_q x + B_q u
 *   AFFINE:     ẋ = A_q x + B_q u + c_q
 *   NONLINEAR:  ẋ = f_q(x, u)  [general nonlinear]
 *   BILINEAR:   ẋ = A_q x + Σ u_i N_i x + B_q u
 *   LPV:        ẋ = A_q(ρ) x + B_q(ρ) u  [linear parameter varying]
 *   POLYNOMIAL: ẋ = Σ a_q,α x^α  [polynomial vector field]
 */
typedef enum {
    HSS_CLASS_LINEAR      = 0,
    HSS_CLASS_AFFINE      = 1,
    HSS_CLASS_NONLINEAR   = 2,
    HSS_CLASS_BILINEAR    = 3,
    HSS_CLASS_LPV         = 4,
    HSS_CLASS_POLYNOMIAL  = 5
} HSS_SystemClass;

/**
 * @brief Stability concept taxonomy (L1 KP7)
 *
 * Different notions of stability for hybrid/switched systems.
 * Based on Goebel/Sanfelice/Teel (2012) and Liberzon (2003).
 */
typedef enum {
    HSS_STAB_LYAPUNOV             = 0,   /* ||x(0)|| < δ ⇒ ||x(t)|| < ε    */
    HSS_STAB_ASYMPTOTIC           = 1,   /* x(t) → 0 as t → ∞              */
    HSS_STAB_EXPONENTIAL          = 2,   /* ||x(t)|| ≤ c·e^{-λt}·||x(0)|| */
    HSS_STAB_GLOBAL_ASYMPTOTIC    = 3,   /* GAS: global asymptotic stable  */
    HSS_STAB_GLOBAL_EXPONENTIAL   = 4,   /* GES: global exponential stable */
    HSS_STAB_UNIFORM_ASYMPTOTIC   = 5,   /* UAS: uniform w.r.t. switching   */
    HSS_STAB_UNIFORM_EXPONENTIAL  = 6,   /* UES: unif. exponential stable   */
    HSS_STAB_INPUT_TO_STATE       = 7,   /* ISS: stable under bounded input  */
    HSS_STAB_PRACTICAL            = 8,   /* Practical: stable to bounded set */
    HSS_STAB_INCREMENTAL          = 9,   /* Incremental: contraction property */
    HSS_STAB_STRING_PASSIVITY     = 10,  /* Passivity-based stability        */
    HSS_STAB_UNSTABLE             = 11   /* Not stable                       */
} HSS_StabilityConcept;

/**
 * @brief Numerical integration method (L1 KP8)
 *
 * Taxonomy of solvers for simulating continuous dynamics within each mode.
 */
typedef enum {
    HSS_SOLVER_EULER_FORWARD      = 0,    /* Forward Euler: O(h)            */
    HSS_SOLVER_EULER_BACKWARD     = 1,    /* Backward Euler: O(h), A-stable */
    HSS_SOLVER_HEUN               = 2,    /* Heun's method: O(h²)           */
    HSS_SOLVER_RK4                = 3,    /* Classical RK4: O(h⁴)           */
    HSS_SOLVER_RKF45              = 4,    /* Runge-Kutta-Fehlberg: adaptive */
    HSS_SOLVER_DORMAND_PRINCE     = 5,    /* DOPRI5(4): adaptive O(h⁵)      */
    HSS_SOLVER_IMPLICIT_MIDPOINT  = 6,    /* Implicit midpoint: symplectic  */
    HSS_SOLVER_TRAPEZOIDAL        = 7,    /* Trapezoidal: O(h²), A-stable   */
    HSS_SOLVER_RADAU_IIA          = 8     /* Radau IIA: fully implicit O(h⁵)*/
} HSS_SolverType;

/**
 * @brief Event type for hybrid system transitions (L1 KP9)
 *
 * Events trigger discrete transitions. Classification follows
 * the CPS event taxonomy from Alur (2015).
 */
typedef enum {
    HSS_EVENT_GUARD_CROSSING      = 0,    /* State enters guard region       */
    HSS_EVENT_INVARIANT_VIOLATION  = 1,   /* State leaves invariant region   */
    HSS_EVENT_TIME_TRIGGERED      = 2,    /* Periodic/time-based event       */
    HSS_EVENT_EXTERNAL_INPUT      = 3,    /* External command/input event    */
    HSS_EVENT_STATE_THRESHOLD     = 4,    /* State variable crosses threshold*/
    HSS_EVENT_RESET_TRIGGERED     = 5,    /* Reset condition activated       */
    HSS_EVENT_SAFETY_VIOLATION    = 6,    /* Unsafe region entered           */
    HSS_EVENT_ZENO_DETECTED       = 7,    /* Zeno behavior detected           */
    HSS_EVENT_SIMULATION_END      = 8     /* Maximum time reached            */
} HSS_EventType;

/**
 * @brief Reset map type (L1 KP10)
 *
 * Describes how the continuous state changes at switching instants.
 * From the hybrid automata formalism.
 */
typedef enum {
    HSS_RESET_IDENTITY            = 0,    /* x⁺ = x⁻ (continuous)           */
    HSS_RESET_ZERO                = 1,    /* x⁺ = 0 (full reset)             */
    HSS_RESET_LINEAR              = 2,    /* x⁺ = M x⁻ (linear reset)        */
    HSS_RESET_AFFINE              = 3,    /* x⁺ = M x⁻ + b (affine reset)   */
    HSS_RESET_PROJECTION          = 4,    /* x⁺ = Π(x⁻) (projection)         */
    HSS_RESET_STATE_DEPENDENT     = 5,    /* x⁺ = ρ(x⁻) (general nonlinear) */
    HSS_RESET_PARTIAL             = 6,    /* Only some states reset           */
    HSS_RESET_FILTER              = 7     /* Filter-style reset (Clegg)       */
} HSS_ResetType;

/* ============================================================================
 * L1 KP1-KP2: Core Data Structures
 * ============================================================================ */

/**
 * @brief Continuous state vector (L1 KP5)
 *
 * x ∈ R^n: the continuous component of the hybrid state.
 * Hybrid state = (q, x) ∈ Q × R^n.
 */
typedef struct {
    double  *data;      /**< State values x[0]...x[n-1]             */
    int      dim;        /**< Dimension n of the state space         */
    char    *labels;     /**< Optional labels (n × HSS_NAME_LEN)     */
    double   time;       /**< Current continuous time t ∈ R≥0        */
    int64_t  jumps;      /**< Number of discrete jumps (j index)      */
} HSS_ContinuousState;

/**
 * @brief System matrix for linear modes (L1)
 *
 * For a linear subsystem: ẋ = A x + B u
 * Stores A ∈ R^{n×n} and B ∈ R^{n×m}.
 */
typedef struct {
    double  *A;          /**< State matrix (n×n, row-major)          */
    double  *B;          /**< Input matrix (n×m, row-major)          */
    double  *c;          /**< Affine term c ∈ R^n (NULL if linear)   */
    int      n;          /**< State dimension                        */
    int      m;          /**< Input dimension                        */
    bool     has_affine; /**< True if affine term c is present       */
} HSS_SystemMatrix;

/**
 * @brief Guard condition for a discrete transition (L1)
 *
 * A guard is a predicate over (q, x) that enables a transition.
 * Guard: state x satisfies linear constraint aᵀx ≤ b.
 * Multiple guards are AND-combined to form the transition guard set.
 */
typedef struct {
    double  *a;          /**< Normal vector a ∈ R^n                  */
    double   b;          /**< Threshold scalar b                      */
    int      dim;        /**< Dimension n of a                        */
    bool     is_active;  /**< Whether this guard condition is active  */
} HSS_Guard;

/**
 * @brief Discrete transition between modes (L1 KP4)
 *
 * Edge e = (q_src, guard, reset, q_dst) ∈ E
 * When the continuous state satisfies the guard condition,
 * the system may transition from mode q_src to q_dst,
 * applying the reset map to the continuous state.
 */
typedef struct {
    int             src_mode;    /**< Source mode index q_src          */
    int             dst_mode;    /**< Destination mode index q_dst     */
    HSS_Guard      *guards;      /**< Guard conditions (AND-combined) */
    int             num_guards;  /**< Number of guard conditions       */
    HSS_ResetType   reset_type;  /**< Type of reset map               */
    double         *reset_M;     /**< Reset matrix M (if linear reset)*/
    double         *reset_b;     /**< Reset offset b (if affine)      */
    int             reset_dim;   /**< Dimension of reset map           */
    double          priority;    /**< Transition priority (if multiple)*/
    bool            is_urgent;   /**< Urgent: must fire when enabled   */
    char            label[HSS_NAME_LEN]; /**< Human-readable label    */
} HSS_Transition;

/**
 * @brief Discrete operating mode (L1 KP2)
 *
 * A mode q ∈ Q represents a discrete operating state.
 * Within each mode, continuous dynamics evolve according to f_q.
 * The system must satisfy the mode invariant while in mode q.
 */
typedef struct {
    int              id;              /**< Mode identifier              */
    char             name[HSS_NAME_LEN]; /**< Mode name                */
    HSS_SystemClass  dynamics_class;  /**< Type of vector field        */
    HSS_SystemMatrix matrix;          /**< Linear/affine dynamics data  */

    /** Nonlinear dynamics: ẋ = f(x, u, t); NULL for linear modes */
    void (*nonlinear_flow)(const double *x, const double *u,
                           double t, double *dxdt, int n, void *params);
    void            *flow_params;      /**< Parameters for nonlinear flow */

    /** Mode invariant: set of linear constraints aᵀx ≤ b */
    HSS_Guard       *invariants;      /**< Invariant constraints        */
    int              num_invariants;  /**< Number of invariant constraints */

    /** Lyapunov function for this mode (for stability analysis) */
    double          *lyapunov_P;       /**< Quadratic LF matrix P       */
    double           lyapunov_decay;   /**< Decay rate: V̇ ≤ -α V       */
    bool             is_stable;        /**< Mode stability status        */

    /** Output equation: y = C x + D u */
    double          *C;                /**< Output matrix (p×n)         */
    double          *D;                /**< Feedthrough matrix (p×m)    */
    int              output_dim;       /**< Output dimension p          */
} HSS_Mode;

/**
 * @brief Switching signal (L1 KP3)
 *
 * σ(t): R≥0 → Q is a piecewise-constant function.
 * It determines which mode q ∈ Q is active at time t.
 * The switching signal may be:
 *   - Arbitrary (no constraints)
 *   - State-dependent: σ(t) = φ(x(t))
 *   - Time-dependent: σ(t) = ψ(t)
 *   - Constrained by dwell time: t_{k+1} - t_k ≥ τ_d
 *   - Constrained by average dwell time
 */
typedef struct {
    int             *mode_sequence;   /**< q_0, q_1, ... sequence       */
    double          *switch_times;    /**< t_0, t_1, ... switch instants*/
    int              num_switches;    /**< Number of switches κ          */
    int              max_sequence;    /**< Allocated sequence length     */
    double           dwell_time_min;  /**< Minimum dwell time τ_d       */
    double           dwell_time_max;  /**< Maximum dwell time (or INF)  */
    double           average_dwell;   /**< Average dwell time τ_a       */
    int              chatter_bound;   /**< Max switches per time N₀     */
    char             signal_name[HSS_NAME_LEN]; /**< Identifier        */
} HSS_SwitchingSignal;

/**
 * @brief Unified Hybrid Switched System (L1 KP1)
 *
 * H = (Q, X, Σ, F, Init, Inv, E, G, R)
 *
 * This is the central data structure representing a complete
 * hybrid switched system. It combines:
 *   - Discrete automaton topology (modes + transitions)
 *   - Continuous dynamics per mode
 *   - Switching signal constraints
 *   - Initial conditions
 */
typedef struct {
    /* == Discrete structure == */
    HSS_Mode          *modes;          /**< Array of discrete modes Q    */
    int                num_modes;       /**< Number of modes |Q|         */
    HSS_Transition    *transitions;    /**< Array of edges E            */
    int                num_transitions; /**< Number of transitions |E|   */

    /* == Continuous structure == */
    int                state_dim;       /**< Dimension n of state space  */
    int                input_dim;       /**< Dimension m of input space  */
    int                output_dim;      /**< Dimension p of output space */
    HSS_ContinuousState state;          /**< Current hybrid state (q, x)*/

    /* == Switching signal == */
    HSS_SwitchingSignal signal;         /**< Switching signal σ(t)      */
    int                active_mode;     /**< Currently active mode q     */

    /* == System properties == */
    HSS_SystemClass    system_class;    /**< Overall system class        */
    char               name[HSS_NAME_LEN]; /**< System identifier       */
    double             time_horizon;    /**< Simulation/synthesis horizon*/

    /* == Solver configuration == */
    HSS_SolverType     solver;          /**< Numerical integration method*/
    double             step_size;       /**< Integration step size h     */
    double             tolerance;       /**< Numerical tolerance ε       */

    /* == Zeno detection == */
    bool               zeno_check;      /**< Enable Zeno detection       */
    double             zeno_min_time;   /**< Minimum inter-event time    */
    int                zeno_max_jumps;  /**< Max jumps before Zeno alarm */
} HSS_System;

/**
 * @brief Hybrid execution trace (L2 KC1)
 *
 * Records the evolution of the hybrid system over hybrid time.
 * Hybrid time: τ = ∪_{j=0}^{J} [t_j, t_{j+1}] × {j}
 * An execution is a sequence (τ, q, x) over the hybrid time domain.
 */
typedef struct {
    double  *times;          /**< Continuous time stamps t_k             */
    int     *modes;          /**< Active mode at each step               */
    double  *states;         /**< State vectors (num_steps × state_dim)  */
    int      num_steps;      /**< Number of time points recorded         */
    int      max_steps;      /**< Allocated capacity                     */
    int      num_jumps;      /**< Number of discrete jumps J             */
    double   total_time;     /**< Total continuous time T                */
} HSS_ExecutionTrace;

/**
 * @brief Hybrid system simulation configuration
 */
typedef struct {
    HSS_SolverType  solver;          /**< Integration method             */
    double          t_start;         /**< Start time t₀                  */
    double          t_end;           /**< End time T                     */
    double          dt;              /**< Fixed step size (ignored if    */
                                      /*   adaptive solver is chosen)    */
    double          dt_min;          /**< Minimum step size (adaptive)   */
    double          dt_max;          /**< Maximum step size (adaptive)   */
    double          abstol;          /**< Absolute error tolerance       */
    double          reltol;          /**< Relative error tolerance       */
    int             max_steps;       /**< Maximum simulation steps       */
    bool            detect_events;   /**< Enable event detection         */
    bool            record_trace;    /**< Record execution trace         */
    double          event_tol;       /**< Event detection tolerance      */
} HSS_SimConfig;

/* ============================================================================
 * L1 API: Core Functions
 * ============================================================================ */

/**
 * @brief Allocate and initialize a hybrid switched system
 * @param name System identifier
 * @param num_modes Number of discrete modes
 * @param state_dim Dimension of continuous state space
 * @param input_dim Dimension of input space
 * @return Newly allocated HSS_System, or NULL on failure
 */
HSS_System *hss_system_create(const char *name, int num_modes,
                               int state_dim, int input_dim);

/**
 * @brief Deallocate a hybrid switched system and all its resources
 * @param sys System to destroy
 */
void hss_system_destroy(HSS_System *sys);

/**
 * @brief Set the dynamics for a specific mode
 * @param sys HSS system
 * @param mode_id Mode index
 * @param cls System class (linear/affine/nonlinear)
 * @param A State matrix (n×n, row-major), or NULL
 * @param B Input matrix (n×m, row-major), or NULL
 * @param c Affine term (n-vector), or NULL
 * @return 0 on success, -1 on error
 *
 * For nonlinear modes, use hss_mode_set_nonlinear_flow() instead.
 */
int hss_mode_set_dynamics(HSS_System *sys, int mode_id,
                           HSS_SystemClass cls,
                           const double *A, const double *B,
                           const double *c);

/**
 * @brief Set nonlinear flow function for a mode
 * @param sys HSS system
 * @param mode_id Mode index
 * @param flow_func Function computing ẋ = f(x, u, t)
 * @param params User-defined parameters
 * @return 0 on success, -1 on error
 */
int hss_mode_set_nonlinear_flow(HSS_System *sys, int mode_id,
                                 void (*flow_func)(const double *x,
                                     const double *u, double t,
                                     double *dxdt, int n, void *params),
                                 void *params);

/**
 * @brief Add a discrete transition between modes
 * @param sys HSS system
 * @param src Source mode index
 * @param dst Destination mode index
 * @param guards Array of guard conditions
 * @param num_guards Number of guards
 * @param reset_type Type of reset map
 * @param reset_M Reset matrix (NULL for identity/no reset)
 * @param reset_b Reset offset (NULL if not affine)
 * @param label Human-readable transition label
 * @return Transition index on success, -1 on error
 */
int hss_add_transition(HSS_System *sys, int src, int dst,
                        const HSS_Guard *guards, int num_guards,
                        HSS_ResetType reset_type,
                        const double *reset_M, const double *reset_b,
                        const char *label);

/**
 * @brief Add an invariant constraint to a mode
 * @param sys HSS system
 * @param mode_id Mode index
 * @param a Normal vector a ∈ R^n
 * @param b Threshold scalar
 * @return 0 on success, -1 on error
 */
int hss_mode_add_invariant(HSS_System *sys, int mode_id,
                            const double *a, double b);

/**
 * @brief Set initial hybrid state (q₀, x₀)
 * @param sys HSS system
 * @param initial_mode Starting mode index
 * @param x0 Initial continuous state vector
 * @return 0 on success, -1 on error
 */
int hss_set_initial_state(HSS_System *sys, int initial_mode,
                           const double *x0);

/**
 * @brief Set Lyapunov matrix for a mode
 * @param sys HSS system
 * @param mode_id Mode index
 * @param P Lyapunov matrix P ∈ R^{n×n} (symmetric positive definite)
 * @param decay_rate Decay rate α > 0 such that V̇ ≤ -α V
 * @return 0 on success, -1 on error
 */
int hss_mode_set_lyapunov(HSS_System *sys, int mode_id,
                           const double *P, double decay_rate);

/**
 * @brief Validate system consistency
 *
 * Checks: mode count, transition connectivity, invariant feasibility,
 * guard well-formedness, dimension consistency.
 *
 * @param sys HSS system
 * @return 0 if valid, -1 with diagnostic if invalid
 */
int hss_system_validate(const HSS_System *sys);

/**
 * @brief Print a human-readable summary of the system
 * @param sys HSS system
 * @param fp Output stream
 */
void hss_system_print_summary(const HSS_System *sys, FILE *fp);

/* ---- L2: Hybrid time and invariant checking ---- */

/**
 * @brief Advance hybrid time by dt
 * @param state Continuous state to advance
 * @param dt Time increment
 */
void hss_hybrid_time_advance(HSS_ContinuousState *state, double dt);

/**
 * @brief Increment jump counter in hybrid time
 * @param state Continuous state
 */
void hss_hybrid_time_jump(HSS_ContinuousState *state);

/**
 * @brief Check if state satisfies all invariants of current mode
 * @param sys HSS system
 * @return true if state is inside invariant set
 */
bool hss_check_invariant(const HSS_System *sys);

/**
 * @brief Check if a specific guard condition is satisfied
 * @param guard Guard condition
 * @param x State vector
 * @return true if guard is enabled
 */
bool hss_check_guard_condition(const HSS_Guard *guard, const double *x);

/**
 * @brief Check if system is Zeno-prone
 * @param sys HSS system
 * @return true if Zeno behavior suspected
 */
bool hss_is_zeno_prone(const HSS_System *sys);

/**
 * @brief Estimate Zeno time for geometric sequence of inter-event times
 * @param first_interval Time of first interval
 * @param decay_rate Rate of decrease of intervals
 * @return Estimated Zeno accumulation time (INFINITY if no Zeno)
 */
double hss_estimate_zeno_time(double first_interval, double decay_rate);

#endif /* HSS_CORE_H */
