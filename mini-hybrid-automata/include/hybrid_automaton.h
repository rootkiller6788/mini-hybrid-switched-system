/**
 * @file hybrid_automaton.h
 * @brief Core hybrid automaton data structures and API (L1 Definitions)
 *
 * A hybrid automaton H = (Q, X, Init, Flow, Inv, E, Guard, Reset) models
 * systems with both discrete state transitions and continuous dynamics.
 *
 * Reference: Thomas A. Henzinger, "The Theory of Hybrid Automata" (1996)
 *            Alur, Courcoubetis, Henzinger, Ho, "Hybrid Automata: An
 *            Algorithmic Approach to the Specification and Verification
 *            of Hybrid Systems" (1993)
 *
 * Course Mapping (L1 Definitions):
 *   MIT 6.841   - Hybrid systems modeling
 *   Stanford CS359 - Hybrid automata fundamentals
 *   Berkeley EECS 291E - Hybrid systems
 *   ETH 227-0690 - Hybrid system theory
 *
 * Knowledge points implemented:
 *   KP1: Hybrid automaton tuple definition (struct HybridAutomaton)
 *   KP2: Discrete mode/location definition (struct HybridMode)
 *   KP3: Continuous variable definition (struct HybridVariable)
 *   KP4: Transition/edge definition (struct HybridTransition)
 *   KP5: Guard condition definition (struct HybridGuard)
 *   KP6: Reset map definition (struct HybridReset)
 *   KP7: Invariant condition definition (struct HybridInvariant)
 *   KP8: Continuous flow/dynamics definition (struct HybridFlow)
 *   KP9: Initial condition definition (struct HybridInit)
 */

#ifndef HYBRID_AUTOMATON_H
#define HYBRID_AUTOMATON_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/* ==========================================================================
 * Configuration constants
 * ========================================================================== */

/** Maximum number of discrete modes (locations) in a hybrid automaton */
#define HA_MAX_MODES        64

/** Maximum number of continuous state variables */
#define HA_MAX_VARIABLES    32

/** Maximum number of discrete transitions (edges) between modes */
#define HA_MAX_TRANSITIONS  256

/** Maximum number of linear constraints per invariant or guard */
#define HA_MAX_CONSTRAINTS  64

/** Maximum length of a mode name string */
#define HA_NAME_LEN         64

/** Maximum length of a variable name string */
#define HA_VARNAME_LEN      32

/* ==========================================================================
 * Core type definitions
 * ========================================================================== */

/**
 * @brief KP1: Hybrid Automaton tuple H = (Q, X, Init, Flow, Inv, E, G, R)
 *
 * This struct captures the complete formal definition of a hybrid automaton
 * as given in Henzinger (1996). The automaton transitions between discrete
 * modes (locations) while continuous variables evolve according to mode-
 * specific differential equations.
 */
typedef struct HybridAutomaton HybridAutomaton;

/**
 * @brief KP2: Discrete mode (location) q ∈ Q
 *
 * Each mode represents a discrete operating state. Within a mode,
 * continuous variables evolve according to the mode's flow function.
 * The system must remain within the mode's invariant region.
 *
 * Mathematical structure:
 *   Mode q has: name, invariant Inv(q) ⊆ ℝⁿ, flow f_q(x) = ẋ
 */
typedef struct HybridMode HybridMode;

/**
 * @brief KP3: Continuous state variable x_i ∈ X
 *
 * Variables represent physical quantities (position, velocity, temperature,
 * etc.) that evolve continuously over time.
 *
 * Domains supported:
 *   - Real-valued (ℝ)
 *   - Bounded interval (ℝ on [lo, hi])
 *   - Discrete-valued (ℤ, used in timed automata clocks)
 */
typedef struct HybridVariable HybridVariable;

/**
 * @brief KP4: Discrete transition (edge) e ∈ E ⊆ Q × Q
 *
 * A transition from source mode to target mode, enabled when its guard
 * condition is satisfied. Upon firing, the reset map updates continuous
 * variables to new values.
 *
 * Mathematical structure:
 *   Edge e = (q_src → q_tgt, Guard(e), Reset(e))
 */
typedef struct HybridTransition HybridTransition;

/**
 * @brief KP5: Guard condition G(e) ⊆ ℝⁿ for a transition
 *
 * The guard is a predicate on continuous state. When the continuous state
 * satisfies the guard, the transition is _enabled_. Guards are typically
 * conjunctions of linear inequalities: A·x ≤ b.
 */
typedef struct HybridGuard HybridGuard;

/**
 * @brief KP6: Reset map R(e): ℝⁿ → ℝⁿ for a transition
 *
 * The reset map specifies how continuous variables are updated when a
 * discrete transition fires. Supports affine resets: x' = R·x + r,
 * where R ∈ ℝ^{n×n} and r ∈ ℝⁿ.
 */
typedef struct HybridReset HybridReset;

/**
 * @brief KP7: Invariant condition Inv(q) ⊆ ℝⁿ for a mode
 *
 * The invariant specifies the region of the continuous state space within
 * which the system must remain while in mode q. If the state reaches the
 * boundary of the invariant, a discrete transition must occur (forced
 * transition).
 *
 * Represents linear constraints: H·x ≤ k.
 */
typedef struct HybridInvariant HybridInvariant;

/**
 * @brief KP8: Continuous flow/dynamics f_q: ℝⁿ → ℝⁿ for a mode
 *
 * The flow defines how continuous variables evolve over time in mode q.
 * Supports:
 *   - Affine ODE:  ẋ = A_q·x + b_q
 *   - Linear ODE:  ẋ = A_q·x
 *   - Constant:    ẋ = b_q
 *   - General ODE: ẋ = f_q(x) via function pointer
 *
 * The solution x(t) satisfies x(0) = x₀ and ẋ(t) = f_q(x(t)).
 */
typedef struct HybridFlow HybridFlow;

/**
 * @brief KP9: Initial condition Init ⊆ Q₀ × X₀
 *
 * Specifies the set of initial (mode, continuous_state) pairs from which
 * executions may begin. For deterministic initialization, a single
 * (q₀, x₀) pair is given.
 */
typedef struct HybridInit HybridInit;

/* ==========================================================================
 * Variable domain type
 * ========================================================================== */

typedef enum {
    HAVAR_REAL,          /**< Unbounded real variable */
    HAVAR_REAL_BOUNDED,  /**< Real variable with [lo, hi] domain */
    HAVAR_INTEGER,       /**< Discrete integer variable */
    HAVAR_CLOCK           /**< Timed automaton clock (rate 1, reset to 0) */
} HybridVariableType;

/* ==========================================================================
 * Mode type classification
 * ========================================================================== */

typedef enum {
    HAMODE_NORMAL,       /**< Standard operating mode */
    HAMODE_INITIAL,      /**< Initial mode (system starts here) */
    HAMODE_TERMINAL,     /**< Terminal/accepting mode */
    HAMODE_URGENT        /**< Urgent mode (no time can elapse) */
} HybridModeType;

/* ==========================================================================
 * Transition trigger classification
 * ========================================================================== */

typedef enum {
    HATRIG_AUTONOMOUS,   /**< Transition taken automatically when guard true */
    HATRIG_CONTROLLED,   /**< Transition taken by controller decision */
    HATRIG_SYNCHRONIZED, /**< Synchronized transition (input event) */
    HATRIG_FORCED        /**< Forced transition (invariant about to be violated) */
} HybridTransitionType;

/* ==========================================================================
 * Flow type classification
 * ========================================================================== */

typedef enum {
    HAFLOW_CONSTANT,     /**< ẋ = b (constant derivative) */
    HAFLOW_LINEAR,       /**< ẋ = A·x (linear ODE, origin-centered) */
    HAFLOW_AFFINE,       /**< ẋ = A·x + b (affine ODE) */
    HAFLOW_NONLINEAR,    /**< ẋ = f(x) (general nonlinear ODE) */
    HAFLOW_DIFF_INCL     /**< ẋ ∈ F(x) (differential inclusion) */
} HybridFlowType;

/* ==========================================================================
 * Reset type classification
 * ========================================================================== */

typedef enum {
    HARESET_IDENTITY,    /**< x' = x (continuous state unchanged) */
    HARESET_CONSTANT,    /**< x' = c (jump to constant value) */
    HARESET_AFFINE,      /**< x' = R·x + r (affine reset) */
    HARESET_CLOCK_RESET  /**< x' = 0 (clock reset in timed automata) */
} HybridResetType;

/* ==========================================================================
 * Core struct definitions
 * ========================================================================== */

/**
 * @brief HybridVariable: continuous state variable
 *
 * KP3 Implementation: Each variable has a name, type, and optional bounds.
 * Clock variables (used in timed automata) have rate = 1 and natural bound = [0,∞).
 */
struct HybridVariable {
    char   name[HA_VARNAME_LEN];   /**< Variable name (e.g., "x", "temp", "clock1") */
    HybridVariableType type;       /**< Domain type */
    double lo;                     /**< Lower bound (if bounded/clock) */
    double hi;                     /**< Upper bound (if bounded) */
    double rate;                   /**< Derivative (1.0 for clock variables) */
};

/**
 * @brief HybridGuard: guard condition for a transition
 *
 * KP5 Implementation: A guard is a conjunction of linear constraints G(e) = {x | A·x ≤ b}.
 * Guards may also represent boolean combinations via multiple constraint groups.
 */
struct HybridGuard {
    int      num_constraints;                        /**< Number of linear constraints */
    double   A[HA_MAX_CONSTRAINTS][HA_MAX_VARIABLES]; /**< Constraint coefficient matrix */
    double   b[HA_MAX_CONSTRAINTS];                  /**< Constraint upper bounds */
    bool     is_trivially_true;                      /**< Guard ≡ true (always enabled) */
};

/**
 * @brief HybridReset: reset map for a transition
 *
 * KP6 Implementation: x' = R·x + r where R ∈ ℝ^{n×n}, r ∈ ℝⁿ.
 * Supports identity reset (x' = x), clock reset (x' = 0),
 * constant reset (x' = c), and general affine reset.
 */
struct HybridReset {
    HybridResetType type;                          /**< Reset classification */
    double   R[HA_MAX_VARIABLES][HA_MAX_VARIABLES]; /**< Reset matrix (n×n) */
    double   r[HA_MAX_VARIABLES];                   /**< Reset offset vector */
};

/**
 * @brief HybridInvariant: invariant region for a mode
 *
 * KP7 Implementation: Inv(q) = {x | H·x ≤ k}. The system must remain
 * within this region while in mode q. Convex polyhedron in H-representation.
 */
struct HybridInvariant {
    int      num_constraints;                        /**< Number of half-space constraints */
    double   H[HA_MAX_CONSTRAINTS][HA_MAX_VARIABLES]; /**< Constraint normal matrix */
    double   k[HA_MAX_CONSTRAINTS];                  /**< Constraint bound vector */
    bool     is_unbounded;                           /**< Inv ≡ ℝⁿ (no constraint) */
};

/**
 * @brief HybridFlow: continuous dynamics for a mode
 *
 * KP8 Implementation: ẋ = f_q(x). Supports affine ODE ẋ = A·x + b,
 * linear ODE ẋ = A·x, constant derivative ẋ = b, and nonlinear function pointer.
 * Matrix A may be NULL for nonlinear flows (function pointer used instead).
 */
struct HybridFlow {
    HybridFlowType type;                          /**< Flow classification */
    bool     has_A;                               /**< Whether A matrix is active */
    double   A[HA_MAX_VARIABLES][HA_MAX_VARIABLES]; /**< System matrix (n×n, affine ODE) */
    double   b[HA_MAX_VARIABLES];                  /**< Affine offset / constant derivative */
    /** Nonlinear flow function pointer: computes ẋ = f(t, x) */
    void   (*nonlinear_f)(double t, const double *x, double *dxdt, int n, void *ctx);
    void   *nonlinear_ctx;                         /**< Context for nonlinear flow */
};

/**
 * @brief HybridInit: initial condition
 *
 * KP9 Implementation: Init ⊆ Q₀ × X₀. Specifies which mode(s) and
 * continuous state(s) are valid initial points for executions.
 */
struct HybridInit {
    int      init_mode;     /**< Initial discrete mode index, -1 if set of modes */
    bool     has_point;     /**< Whether x0 is a single point */
    double   x0[HA_MAX_VARIABLES]; /**< Initial continuous state (if has_point) */
    /** Set-based initialization (if !has_point): rectangular region [lo_i, hi_i] */
    double   init_lo[HA_MAX_VARIABLES];
    double   init_hi[HA_MAX_VARIABLES];
    bool     init_has_bounds; /**< Whether init_lo/hi are set */
};

/**
 * @brief HybridMode: discrete location/mode
 *
 * KP2 Implementation: Each mode packages its invariant, flow, and
 * optional name/type metadata.
 */
struct HybridMode {
    int      id;              /**< Unique mode identifier (index) */
    char     name[HA_NAME_LEN]; /**< Human-readable name */
    HybridModeType type;      /**< Mode classification */
    HybridInvariant invariant; /**< Inv(q) - region where mode can stay */
    HybridFlow     flow;       /**< f_q - continuous dynamics */
};

/**
 * @brief HybridTransition: discrete edge between modes
 *
 * KP4 Implementation: Edge e = (src → tgt, Guard(e), Reset(e)).
 * A transition fires when the guard is satisfied (and for forced transitions,
 * when the source mode's invariant boundary is reached).
 */
struct HybridTransition {
    int      id;              /**< Unique transition identifier (index) */
    char     name[HA_NAME_LEN]; /**< Transition event name (e.g., "switch_on") */
    int      src_mode;        /**< Source mode index (q_src) */
    int      tgt_mode;        /**< Target mode index (q_tgt) */
    HybridTransitionType type; /**< Trigger classification */
    HybridGuard    guard;      /**< Guard condition G(e) */
    HybridReset    reset;      /**< Reset map R(e) */
    bool     is_active;       /**< Whether this transition is enabled */
};

/**
 * @brief HybridAutomaton: full hybrid automaton
 *
 * KP1 Implementation: H = (Q, X, Init, Flow, Inv, E, Guard, Reset).
 * This is the top-level struct aggregating all components.
 */
struct HybridAutomaton {
    char     name[HA_NAME_LEN];  /**< Automaton name/identifier */
    int      num_modes;          /**< |Q| - number of discrete modes */
    int      num_vars;           /**< |X| = n - number of continuous variables */
    int      num_transitions;    /**< |E| - number of edges */
    HybridMode       modes[HA_MAX_MODES];        /**< Set of discrete modes Q */
    HybridVariable   vars[HA_MAX_VARIABLES];      /**< Set of continuous variables X */
    HybridTransition trans[HA_MAX_TRANSITIONS];   /**< Set of edges E */
    HybridInit       init;                        /**< Initial condition */
    /** Mode-to-transition adjacency: trans_from[q] = list of transition IDs */
    int      trans_from[HA_MAX_MODES][HA_MAX_TRANSITIONS];
    int      trans_from_count[HA_MAX_MODES];
    /** Mode-to-transition adjacency: trans_to[q] = list of incoming transition IDs */
    int      trans_to[HA_MAX_MODES][HA_MAX_TRANSITIONS];
    int      trans_to_count[HA_MAX_MODES];
};

/* ==========================================================================
 * API: Construction and configuration
 * ========================================================================== */

/**
 * @brief Create a new hybrid automaton with given name and variable count.
 *
 * KP1: Instantiates the formal tuple H = (Q, X, Init, Flow, Inv, E, G, R).
 *
 * @param name         Identifier for the automaton
 * @param num_vars     Number of continuous state variables (dimension of X)
 * @return             Pointer to newly allocated HybridAutomaton, or NULL on error
 *
 * Complexity: O(1) allocation
 */
HybridAutomaton* hybrid_automaton_create(const char *name, int num_vars);

/**
 * @brief Destroy a hybrid automaton and free all associated resources.
 *
 * @param ha  Pointer to the automaton to destroy
 *
 * Complexity: O(|Q| + |E|)
 */
void hybrid_automaton_destroy(HybridAutomaton *ha);

/**
 * @brief KP2: Add a new discrete mode (location) to the automaton.
 *
 * Each mode q receives a unique ID, a name, and is initialized with
 * an unbounded invariant and zero flow (steady state).
 *
 * @param ha    The automaton to modify
 * @param name  Human-readable mode name
 * @param type  Mode classification (NORMAL, INITIAL, TERMINAL, URGENT)
 * @return      Mode index (≥0) on success, -1 if max modes reached
 *
 * Complexity: O(1) insertion
 * Theorem reference: mode set Q is finite by definition (Henzinger 1996, Def 2.1)
 */
int  hybrid_automaton_add_mode(HybridAutomaton *ha, const char *name, HybridModeType type);

/**
 * @brief KP3: Set a continuous variable's name and domain type.
 *
 * Each variable x_i ∈ X represents a physical quantity tracked across modes.
 *
 * @param ha         The automaton
 * @param var_index  Variable index (0 to num_vars-1)
 * @param name       Variable name
 * @param type       Domain type
 * @param lo         Lower bound (if bounded)
 * @param hi         Upper bound (if bounded)
 * @return           true on success, false on invalid index
 */
bool hybrid_automaton_set_variable(HybridAutomaton *ha, int var_index,
                                    const char *name, HybridVariableType type,
                                    double lo, double hi);

/**
 * @brief KP4: Add a discrete transition (edge) between two modes.
 *
 * Transition e = (src → tgt, guard, reset). The transition is enabled
 * when the guard predicate is satisfied by the current continuous state.
 *
 * @param ha     The automaton
 * @param src    Source mode index
 * @param tgt    Target mode index
 * @param type   Transition trigger classification
 * @param name   Transition event name
 * @return       Transition index on success, -1 if max exceeded or invalid modes
 */
int  hybrid_automaton_add_transition(HybridAutomaton *ha, int src, int tgt,
                                      HybridTransitionType type, const char *name);

/**
 * @brief KP5: Configure the guard condition for a transition.
 *
 * Guard G(e) = {x ∈ ℝⁿ | A·x ≤ b}. A conjunction of linear constraints.
 * For trivially true guards, pass num_constraints = 0.
 *
 * @param ha             The automaton
 * @param trans_id       Transition index
 * @param num_constraints Number of linear inequalities
 * @param A              Constraint matrix (num_constraints × n), may be NULL if 0
 * @param b              Constraint bounds (length num_constraints), may be NULL if 0
 * @return               true on success
 *
 * Theorem: Linear guard guards preserve the semi-algebraic property of
 *          reachable sets in initialized rectangular automata.
 */
bool hybrid_guard_set(HybridAutomaton *ha, int trans_id,
                      int num_constraints, const double *A, const double *b);

/**
 * @brief KP6: Configure the reset map for a transition.
 *
 * Reset R(e): ℝⁿ → ℝⁿ, x' = R·x + r. Supports identity, constant,
 * clock reset, and affine reset types.
 *
 * @param ha       The automaton
 * @param trans_id Transition index
 * @param type     Reset classification
 * @param R        Reset matrix (n×n in row-major), may be NULL for identity/constant
 * @param r        Reset offset vector (length n), may be NULL
 * @return         true on success
 */
bool hybrid_reset_set(HybridAutomaton *ha, int trans_id,
                      HybridResetType type, const double *R, const double *r);

/**
 * @brief KP7: Configure the invariant region for a mode.
 *
 * Inv(q) = {x ∈ ℝⁿ | H·x ≤ k}. The system must stay within this
 * polyhedron while in mode q. For unbounded modes, pass num_constraints = 0.
 *
 * @param ha             The automaton
 * @param mode_id        Mode index
 * @param num_constraints Number of half-space constraints
 * @param H              Normal matrix (num_constraints × n)
 * @param k              Bound vector (length num_constraints)
 * @return               true on success
 *
 * Theorem: Invariant violation forces a discrete transition (safety requirement).
 *          Alur et al. (1995), "The Algorithmic Analysis of Hybrid Systems"
 */
bool hybrid_invariant_set(HybridAutomaton *ha, int mode_id,
                          int num_constraints, const double *H, const double *k);

/**
 * @brief KP8: Configure the continuous flow for a mode.
 *
 * Flow f_q: ẋ = A·x + b (affine ODE) or ẋ = f(t,x) (nonlinear, via function pointer).
 *
 * @param ha      The automaton
 * @param mode_id Mode index
 * @param type    Flow classification
 * @param A       System matrix (n×n, row-major), may be NULL
 * @param b       Affine offset / constant vector (length n), may be NULL
 * @return        true on success
 */
bool hybrid_flow_set(HybridAutomaton *ha, int mode_id,
                     HybridFlowType type, const double *A, const double *b);

/**
 * @brief Set a nonlinear flow function for a mode.
 *
 * KP8 Extension: For general ODE ẋ = f(t, x), provides a callback that
 * computes the derivative at any given (t, x) point.
 *
 * @param ha      The automaton
 * @param mode_id Mode index
 * @param f       Function pointer computing dx/dt = f(t, x, dxdt, n, ctx)
 * @param ctx     Opaque context pointer passed to f
 * @return        true on success
 */
bool hybrid_flow_set_nonlinear(HybridAutomaton *ha, int mode_id,
                                void (*f)(double, const double*, double*, int, void*),
                                void *ctx);

/**
 * @brief KP9: Set the initial condition of the automaton.
 *
 * Init ⊆ Q₀ × X₀ specifies valid starting states. For point initialization,
 * provides a single (q₀, x₀). For set initialization, provides bounds.
 *
 * @param ha          The automaton
 * @param init_mode   Initial mode index
 * @param x0          Initial continuous state (n-dimensional), may be NULL
 * @return            true on success
 */
bool hybrid_init_set(HybridAutomaton *ha, int init_mode, const double *x0);

/**
 * @brief Set the initial condition as a rectangular region.
 *
 * @param ha        The automaton
 * @param init_mode Initial mode index
 * @param lo        Lower bounds (n-dimensional)
 * @param hi        Upper bounds (n-dimensional)
 * @return          true on success
 */
bool hybrid_init_set_rect(HybridAutomaton *ha, int init_mode,
                           const double *lo, const double *hi);

/* ==========================================================================
 * API: Query and introspection
 * ========================================================================== */

/** Get the number of modes |Q| */
int  hybrid_automaton_mode_count(const HybridAutomaton *ha);

/** Get the number of continuous variables |X| = n */
int  hybrid_automaton_var_count(const HybridAutomaton *ha);

/** Get the number of transitions |E| */
int  hybrid_automaton_trans_count(const HybridAutomaton *ha);

/** Get a mode by index */
const HybridMode* hybrid_automaton_get_mode(const HybridAutomaton *ha, int mode_id);

/** Get a transition by index */
const HybridTransition* hybrid_automaton_get_transition(const HybridAutomaton *ha, int trans_id);

/** Get a variable by index */
const HybridVariable* hybrid_automaton_get_variable(const HybridAutomaton *ha, int var_idx);

/** Get number of outgoing transitions from a mode (out-degree) */
int  hybrid_automaton_out_degree(const HybridAutomaton *ha, int mode_id);

/** Get number of incoming transitions to a mode (in-degree) */
int  hybrid_automaton_in_degree(const HybridAutomaton *ha, int mode_id);

/**
 * @brief Check whether the automaton is well-formed.
 *
 * Validates: (1) all transitions reference valid modes, (2) at least
 * one initial mode, (3) no mode has self-contradictory invariants.
 *
 * @param ha  The automaton
 * @return    true if well-formed
 */
bool hybrid_automaton_is_well_formed(const HybridAutomaton *ha);

/**
 * @brief Print a human-readable description of the automaton.
 *
 * Outputs the automaton's graph structure, modes, transitions, invariants,
 * guards, resets, and initial conditions to stdout for debugging.
 *
 * @param ha  The automaton to print
 */
void hybrid_automaton_print(const HybridAutomaton *ha);

#endif /* HYBRID_AUTOMATON_H */
