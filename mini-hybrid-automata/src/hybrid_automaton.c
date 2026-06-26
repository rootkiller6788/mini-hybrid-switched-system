/**
 * @file hybrid_automaton.c
 * @brief Core hybrid automaton implementation (KP1-KP9, L1-L2)
 *
 * Implements the fundamental data structures and operations for
 * hybrid automata: construction, mode/variable/transition management,
 * guard/reset/invariant/flow/init configuration, and introspection.
 *
 * Reference:
 *   Henzinger, "The Theory of Hybrid Automata" (1996), LICS
 *   Alur et al., "Hybrid Automata: An Algorithmic Approach..." (1993)
 *   Lygeros et al., "Dynamical Properties of Hybrid Automata" (2003), IEEE TAC
 */

#include "hybrid_automaton.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==========================================================================
 * Memory allocation helpers
 * ========================================================================== */

/**
 * @brief Safe string copy with truncation guarantee.
 *
 * Copies at most n-1 characters and always null-terminates.
 *
 * @param dst  Destination buffer
 * @param src  Source string
 * @param n    Buffer size
 */
static void safe_strncpy(char *dst, const char *src, size_t n)
{
    if (!dst || !src || n == 0) return;
    size_t i;
    for (i = 0; i < n - 1 && src[i]; i++) {
        dst[i] = src[i];
    }
    dst[i] = '\0';
}

/**
 * @brief Allocate a zero-initialized memory block.
 *
 * Wraps calloc with NULL check.
 */
static void* safe_calloc(size_t nmemb, size_t size)
{
    void *ptr = calloc(nmemb, size);
    if (!ptr) {
        fprintf(stderr, "[hybrid_automaton] ERROR: calloc(%zu, %zu) failed\n",
                nmemb, size);
    }
    return ptr;
}

/* ==========================================================================
 * Identity matrix initialization
 * ========================================================================== */

/**
 * @brief Initialize the n×n identity matrix in a flat array (row-major).
 *
 * Used for identity reset maps and default system matrices.
 *
 * @param mat  Flat array of size n*n
 * @param n    Dimension
 */
static void matrix_set_identity(double *mat, int n)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mat[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

/**
 * @brief Set all entries of an n×n matrix to zero.
 */
static void matrix_set_zero(double *mat, int n)
{
    for (int i = 0; i < n * n; i++) {
        mat[i] = 0.0;
    }
}

/**
 * @brief Copy an n×n matrix from src to dst (both row-major).
 */
static void matrix_copy(double *dst, const double *src, int n)
{
    for (int i = 0; i < n * n; i++) {
        dst[i] = src[i];
    }
}

/**
 * @brief Copy an n-vector from src to dst.
 */
static void vector_copy(double *dst, const double *src, int n)
{
    for (int i = 0; i < n; i++) {
        dst[i] = src[i];
    }
}

/* ==========================================================================
 * Core creation and destruction (KP1)
 * ========================================================================== */

/**
 * @brief KP1: Create a new hybrid automaton.
 *
 * Instantiates H = (Q, X, Init, Flow, Inv, E, Guard, Reset) with the
 * given name and continuous variable count. All modes, transitions,
 * guards, resets, invariants, flows, and initial conditions are
 * initialized to default (safe) values.
 *
 * Mathematical definition (Henzinger 1996, Definition 2.1):
 *   A hybrid automaton H consists of:
 *   - A finite set Q = {q_1, ..., q_m} of discrete states (modes)
 *   - A finite set X = {x_1, ..., x_n} of real-valued variables
 *   - Initial condition Init ⊆ Q × ℝⁿ
 *   - For each q ∈ Q: a flow condition f_q defining ẋ
 *   - For each q ∈ Q: an invariant condition Inv(q) ⊆ ℝⁿ
 *   - A finite set E ⊆ Q × Q of edges (transitions)
 *   - For each e ∈ E: a guard condition G(e) ⊆ ℝⁿ
 *   - For each e ∈ E: a reset relation R(e) ⊆ ℝⁿ × ℝⁿ
 *
 * @param name      Identifier for this automaton
 * @param num_vars  Number of continuous variables (dimension n)
 * @return          Newly allocated and zero-initialized automaton, or NULL
 *
 * Complexity: O(1) — single allocation with fixed-size internal arrays
 */
HybridAutomaton* hybrid_automaton_create(const char *name, int num_vars)
{
    if (!name || num_vars < 0 || num_vars > HA_MAX_VARIABLES) {
        fprintf(stderr, "[hybrid_automaton_create] Invalid parameters: "
                "name=%p, num_vars=%d (max %d)\n",
                (void*)name, num_vars, HA_MAX_VARIABLES);
        return NULL;
    }

    HybridAutomaton *ha = (HybridAutomaton*) safe_calloc(1, sizeof(HybridAutomaton));
    if (!ha) return NULL;

    /* KP1: Initialize the automaton tuple */
    safe_strncpy(ha->name, name, HA_NAME_LEN);
    ha->num_modes = 0;
    ha->num_vars = num_vars;
    ha->num_transitions = 0;

    /* Initialize mode array: all modes start as unused */
    for (int i = 0; i < HA_MAX_MODES; i++) {
        ha->modes[i].id = -1;
        ha->modes[i].type = HAMODE_NORMAL;
        ha->modes[i].invariant.is_unbounded = true;
        ha->modes[i].invariant.num_constraints = 0;
        ha->modes[i].flow.type = HAFLOW_CONSTANT;
        ha->modes[i].flow.has_A = false;
        for (int j = 0; j < HA_MAX_VARIABLES; j++) {
            ha->modes[i].flow.b[j] = 0.0;
        }
    }

    /* Initialize variable array (KP3) */
    for (int i = 0; i < HA_MAX_VARIABLES; i++) {
        ha->vars[i].type = HAVAR_REAL;
        ha->vars[i].lo = -INFINITY;
        ha->vars[i].hi = INFINITY;
        ha->vars[i].rate = 0.0;
        ha->vars[i].name[0] = '\0';
    }

    /* Initialize transition adjacency lists (KP4) */
    for (int i = 0; i < HA_MAX_MODES; i++) {
        ha->trans_from_count[i] = 0;
        ha->trans_to_count[i] = 0;
        for (int j = 0; j < HA_MAX_TRANSITIONS; j++) {
            ha->trans_from[i][j] = -1;
            ha->trans_to[i][j] = -1;
        }
    }

    /* Initialize transitions */
    for (int i = 0; i < HA_MAX_TRANSITIONS; i++) {
        ha->trans[i].id = -1;
        ha->trans[i].src_mode = -1;
        ha->trans[i].tgt_mode = -1;
        ha->trans[i].is_active = true;
        ha->trans[i].guard.is_trivially_true = true;
        ha->trans[i].guard.num_constraints = 0;
        ha->trans[i].reset.type = HARESET_IDENTITY;
    }

    /* Initialize initial condition (KP9) */
    ha->init.init_mode = -1;
    ha->init.has_point = false;
    ha->init.init_has_bounds = false;
    for (int i = 0; i < HA_MAX_VARIABLES; i++) {
        ha->init.x0[i] = 0.0;
        ha->init.init_lo[i] = -INFINITY;
        ha->init.init_hi[i] = INFINITY;
    }

    return ha;
}

/**
 * @brief Destroy a hybrid automaton and release all memory.
 *
 * Frees the automaton struct. Since all components are inline arrays
 * (not dynamically allocated pointers), a single free suffices.
 *
 * @param ha Automaton to destroy (may be NULL)
 *
 * Complexity: O(1)
 */
void hybrid_automaton_destroy(HybridAutomaton *ha)
{
    free(ha);
}

/* ==========================================================================
 * Mode management (KP2)
 * ========================================================================== */

/**
 * @brief KP2: Add a discrete mode (location) to the automaton.
 *
 * Each mode represents a discrete state with its own continuous dynamics,
 * invariant region, and outgoing/incoming transitions.
 *
 * Mathematical definition:
 *   q ∈ Q (finite set), with associated:
 *   - Inv(q) ⊆ ℝⁿ (invariant: where the system can stay)
 *   - f_q: ℝⁿ → ℝⁿ (flow: how continuous state evolves)
 *
 * @param ha    The automaton
 * @param name  Human-readable mode name
 * @param type  Mode classification (NORMAL, INITIAL, TERMINAL, URGENT)
 * @return      Mode index on success, -1 if full
 *
 * Complexity: O(1) amortized insertion into pre-allocated array
 * Theorem reference: Q must be finite (Henzinger 1996, Def 2.1, condition 1)
 */
int hybrid_automaton_add_mode(HybridAutomaton *ha, const char *name,
                               HybridModeType type)
{
    if (!ha) return -1;
    if (ha->num_modes >= HA_MAX_MODES) {
        fprintf(stderr, "[hybrid_automaton_add_mode] Max modes (%d) reached\n",
                HA_MAX_MODES);
        return -1;
    }

    int idx = ha->num_modes;
    HybridMode *m = &ha->modes[idx];

    m->id = idx;
    safe_strncpy(m->name, name ? name : "", HA_NAME_LEN);
    m->type = type;

    /* Default to unbounded invariant (no constraints) */
    m->invariant.is_unbounded = true;
    m->invariant.num_constraints = 0;

    /* Default to zero flow: ẋ = 0 (steady state) */
    m->flow.type = HAFLOW_CONSTANT;
    m->flow.has_A = false;
    for (int i = 0; i < ha->num_vars; i++) {
        m->flow.b[i] = 0.0;
    }

    ha->num_modes++;

    /* If this is the first mode marked INITIAL, set it */
    if (type == HAMODE_INITIAL && ha->init.init_mode < 0) {
        ha->init.init_mode = idx;
    }

    return idx;
}

/* ==========================================================================
 * Variable management (KP3)
 * ========================================================================== */

/**
 * @brief KP3: Configure a continuous variable.
 *
 * Each x_i ∈ X represents a physical quantity tracked continuously.
 * Variables can be real-valued, bounded, integer-valued, or clocks.
 *
 * Mathematical definition:
 *   X = {x₁, ..., xₙ} is a finite set of variables taking values in ℝ.
 *   A valuation is a function ν: X → ℝ.
 *
 * @param ha        The automaton
 * @param var_index Variable index (0 ≤ var_index < num_vars)
 * @param name      Descriptive name
 * @param type      Domain type
 * @param lo        Lower bound (for bounded/clock types)
 * @param hi        Upper bound (for bounded/clock types)
 * @return          true on success
 *
 * Complexity: O(1)
 */
bool hybrid_automaton_set_variable(HybridAutomaton *ha, int var_index,
                                    const char *name, HybridVariableType type,
                                    double lo, double hi)
{
    if (!ha || var_index < 0 || var_index >= ha->num_vars) {
        return false;
    }

    HybridVariable *v = &ha->vars[var_index];
    safe_strncpy(v->name, name ? name : "", HA_VARNAME_LEN);
    v->type = type;
    v->lo = lo;
    v->hi = hi;

    /* Clock variables always have rate = 1 */
    if (type == HAVAR_CLOCK) {
        v->rate = 1.0;
        v->lo = 0.0;
    }

    return true;
}

/* ==========================================================================
 * Transition management (KP4)
 * ========================================================================== */

/**
 * @brief KP4: Add a discrete transition (edge) between two modes.
 *
 * Transition e = (q_src → q_tgt, G(e), R(e)).
 * A transition fires instantaneously when its guard is satisfied.
 *
 * The adjacency lists trans_from[q_src] and trans_to[q_tgt] are updated.
 *
 * @param ha   The automaton
 * @param src  Source mode index
 * @param tgt  Target mode index
 * @param type Transition trigger classification
 * @param name Transition event name
 * @return     Transition index on success, -1 on error
 *
 * Complexity: O(1) insertion + O(1) adjacency update
 */
int hybrid_automaton_add_transition(HybridAutomaton *ha, int src, int tgt,
                                     HybridTransitionType type, const char *name)
{
    if (!ha) return -1;
    if (src < 0 || src >= ha->num_modes) {
        fprintf(stderr, "[hybrid_automaton_add_transition] Invalid src mode %d\n", src);
        return -1;
    }
    if (tgt < 0 || tgt >= ha->num_modes) {
        fprintf(stderr, "[hybrid_automaton_add_transition] Invalid tgt mode %d\n", tgt);
        return -1;
    }
    if (ha->num_transitions >= HA_MAX_TRANSITIONS) {
        fprintf(stderr, "[hybrid_automaton_add_transition] Max transitions (%d) reached\n",
                HA_MAX_TRANSITIONS);
        return -1;
    }

    int idx = ha->num_transitions;
    HybridTransition *t = &ha->trans[idx];

    t->id = idx;
    safe_strncpy(t->name, name ? name : "", HA_NAME_LEN);
    t->src_mode = src;
    t->tgt_mode = tgt;
    t->type = type;
    t->is_active = true;

    /* Default: trivially true guard (always enabled) */
    t->guard.is_trivially_true = true;
    t->guard.num_constraints = 0;

    /* Default: identity reset (x' = x) */
    t->reset.type = HARESET_IDENTITY;

    /* Update adjacency: outgoing from src */
    int fc = ha->trans_from_count[src];
    if (fc < HA_MAX_TRANSITIONS) {
        ha->trans_from[src][fc] = idx;
        ha->trans_from_count[src] = fc + 1;
    }

    /* Update adjacency: incoming to tgt */
    int tc = ha->trans_to_count[tgt];
    if (tc < HA_MAX_TRANSITIONS) {
        ha->trans_to[tgt][tc] = idx;
        ha->trans_to_count[tgt] = tc + 1;
    }

    ha->num_transitions++;

    return idx;
}

/* ==========================================================================
 * Guard configuration (KP5)
 * ========================================================================== */

/**
 * @brief KP5: Configure a linear guard condition for a transition.
 *
 * Guard G(e) = {x ∈ ℝⁿ | A·x ≤ b}, where A is a num_constraints × n matrix
 * and b is a vector of length num_constraints.  A conjunction of
 * linear inequality constraints.
 *
 * For trivially true guards (always enabled), set num_constraints = 0.
 *
 * Mathematical definition (Henzinger 1996, Def 2.1, condition 7):
 *   For each edge e ∈ E, the guard G(e) ⊆ ℝⁿ defines when the
 *   edge may be taken. G(e) is a predicate on the continuous state.
 *
 * @param ha              The automaton
 * @param trans_id        Transition index
 * @param num_constraints Number of constraints (0 = trivially true)
 * @param A               Constraint matrix (num_constraints × n, row-major)
 * @param b               Constraint bounds (length num_constraints)
 * @return                true on success
 *
 * Complexity: O(num_constraints · n) for copying
 */
bool hybrid_guard_set(HybridAutomaton *ha, int trans_id,
                      int num_constraints, const double *A, const double *b)
{
    if (!ha || trans_id < 0 || trans_id >= ha->num_transitions) return false;
    if (num_constraints < 0 || num_constraints > HA_MAX_CONSTRAINTS) return false;

    HybridTransition *t = &ha->trans[trans_id];
    HybridGuard *g = &t->guard;

    if (num_constraints == 0) {
        /* Trivially true guard: transition always enabled */
        g->is_trivially_true = true;
        g->num_constraints = 0;
        return true;
    }

    if (!A || !b) return false;

    g->is_trivially_true = false;
    g->num_constraints = num_constraints;

    /* Copy constraint matrix A (num_constraints × n) */
    for (int i = 0; i < num_constraints; i++) {
        for (int j = 0; j < ha->num_vars; j++) {
            g->A[i][j] = A[i * ha->num_vars + j];
        }
        g->b[i] = b[i];
    }

    return true;
}

/* ==========================================================================
 * Reset configuration (KP6)
 * ========================================================================== */

/**
 * @brief KP6: Configure the reset map for a transition.
 *
 * Reset R(e): ℝⁿ → ℝⁿ maps the pre-transition continuous state to
 * the post-transition continuous state: x' = R·x + r.
 *
 * Supported reset types:
 *   HARESET_IDENTITY  — x' = x (continuous state unchanged)
 *   HARESET_CONSTANT  — x' = c (jump to fixed value, via r)
 *   HARESET_AFFINE    — x' = R·x + r (general affine reset)
 *   HARESET_CLOCK_RESET — x_i' = 0 for specific clock variables
 *
 * Mathematical definition:
 *   Reset(e) ⊆ ℝⁿ × ℝⁿ relates pre-jump and post-jump states.
 *   For deterministic resets, it's a function ℝⁿ → ℝⁿ.
 *
 * @param ha        The automaton
 * @param trans_id  Transition index
 * @param type      Reset classification
 * @param R         Reset matrix (n×n, row-major), NULL for non-affine
 * @param r         Reset offset (n-vector), NULL for identity
 * @return          true on success
 *
 * Complexity: O(n²) for affine copy, O(1) otherwise
 */
bool hybrid_reset_set(HybridAutomaton *ha, int trans_id,
                      HybridResetType type, const double *R, const double *r)
{
    if (!ha || trans_id < 0 || trans_id >= ha->num_transitions) return false;

    HybridTransition *t = &ha->trans[trans_id];
    HybridReset *reset = &t->reset;

    reset->type = type;

    switch (type) {
    case HARESET_IDENTITY:
        /* x' = x: R = I, r = 0 */
        matrix_set_identity(reset->R[0], ha->num_vars);
        for (int i = 0; i < ha->num_vars; i++) reset->r[i] = 0.0;
        break;

    case HARESET_CONSTANT:
        /* x' = c: R = 0, r = c */
        if (!r) return false;
        matrix_set_zero(reset->R[0], ha->num_vars);
        vector_copy(reset->r, r, ha->num_vars);
        break;

    case HARESET_AFFINE:
        /* x' = R·x + r */
        if (R) {
            matrix_copy(reset->R[0], R, ha->num_vars);
        } else {
            matrix_set_identity(reset->R[0], ha->num_vars);
        }
        if (r) {
            vector_copy(reset->r, r, ha->num_vars);
        } else {
            for (int i = 0; i < ha->num_vars; i++) reset->r[i] = 0.0;
        }
        break;

    case HARESET_CLOCK_RESET:
        /* x_i' = 0 for designated clocks */
        matrix_set_identity(reset->R[0], ha->num_vars);
        if (r) {
            for (int i = 0; i < ha->num_vars; i++) {
                reset->r[i] = r[i];
            }
        } else {
            for (int i = 0; i < ha->num_vars; i++) reset->r[i] = 0.0;
        }
        break;

    default:
        return false;
    }

    return true;
}

/* ==========================================================================
 * Invariant configuration (KP7)
 * ========================================================================== */

/**
 * @brief KP7: Configure the invariant region for a mode.
 *
 * Inv(q) = {x ∈ ℝⁿ | H·x ≤ k}, a convex polyhedron in H-representation.
 * While in mode q, the continuous state must remain within Inv(q).
 *
 * Mathematical definition:
 *   Inv: Q → 2^{ℝⁿ} assigns to each mode an invariant region.
 *   During a continuous flow in mode q, x(t) ∈ Inv(q) ∀t.
 *   If x(t) reaches ∂Inv(q), a discrete transition MUST occur.
 *
 * @param ha              The automaton
 * @param mode_id         Mode index
 * @param num_constraints Number of half-space constraints
 * @param H               Normal matrix (num_constraints × n, row-major)
 * @param k               Bound vector (length num_constraints)
 * @return                true on success
 *
 * Complexity: O(num_constraints · n)
 * Theorem: Invariant boundary crossings force discrete transitions
 *          (Alur et al., "Algorithmic Analysis of Hybrid Systems", 1995)
 */
bool hybrid_invariant_set(HybridAutomaton *ha, int mode_id,
                          int num_constraints, const double *H, const double *k)
{
    if (!ha || mode_id < 0 || mode_id >= ha->num_modes) return false;
    if (num_constraints < 0 || num_constraints > HA_MAX_CONSTRAINTS) return false;

    HybridMode *m = &ha->modes[mode_id];
    HybridInvariant *inv = &m->invariant;

    if (num_constraints == 0) {
        inv->is_unbounded = true;
        inv->num_constraints = 0;
        return true;
    }

    if (!H || !k) return false;

    inv->is_unbounded = false;
    inv->num_constraints = num_constraints;

    for (int i = 0; i < num_constraints; i++) {
        for (int j = 0; j < ha->num_vars; j++) {
            inv->H[i][j] = H[i * ha->num_vars + j];
        }
        inv->k[i] = k[i];
    }

    return true;
}

/* ==========================================================================
 * Flow configuration (KP8)
 * ========================================================================== */

/**
 * @brief KP8: Configure the continuous flow for a mode.
 *
 * Flow f_q: ẋ = A·x + b (affine ODE).
 *
 * Mathematical definition:
 *   Flow: Q → (ℝⁿ → ℝⁿ) assigns to each mode a vector field
 *   defining the continuous evolution: ẋ(t) = f_q(x(t)).
 *
 * For affine flows: ẋ(t) = A_q·x(t) + b_q.
 * The solution (for constant A_q, b_q):
 *   x(t) = e^{A_q t} x₀ + (∫₀ᵗ e^{A_q(t-τ)} dτ) b_q
 *
 * Special case A_q = 0: x(t) = x₀ + b_q·t (constant derivative).
 * Special case b_q = 0: x(t) = e^{A_q t} x₀ (linear ODE).
 *
 * @param ha      The automaton
 * @param mode_id Mode index
 * @param type    Flow classification
 * @param A       System matrix (n×n, row-major), NULL for non-matrical flows
 * @param b       Affine offset (n-vector), NULL for zero offset
 * @return        true on success
 *
 * Complexity: O(n²) for matrix copy
 */
bool hybrid_flow_set(HybridAutomaton *ha, int mode_id,
                     HybridFlowType type, const double *A, const double *b)
{
    if (!ha || mode_id < 0 || mode_id >= ha->num_modes) return false;

    HybridMode *m = &ha->modes[mode_id];
    HybridFlow *flow = &m->flow;

    flow->type = type;
    flow->has_A = (A != NULL);

    /* Copy system matrix A if provided */
    if (A) {
        for (int i = 0; i < ha->num_vars; i++) {
            for (int j = 0; j < ha->num_vars; j++) {
                flow->A[i][j] = A[i * ha->num_vars + j];
            }
        }
    } else {
        for (int i = 0; i < ha->num_vars; i++) {
            for (int j = 0; j < ha->num_vars; j++) {
                flow->A[i][j] = 0.0;
            }
        }
    }

    /* Copy offset vector b if provided */
    if (b) {
        for (int i = 0; i < ha->num_vars; i++) {
            flow->b[i] = b[i];
        }
    } else {
        for (int i = 0; i < ha->num_vars; i++) {
            flow->b[i] = 0.0;
        }
    }

    /* Clear nonlinear flow function for non-NONLINEAR types */
    if (type != HAFLOW_NONLINEAR) {
        flow->nonlinear_f = NULL;
        flow->nonlinear_ctx = NULL;
    }

    return true;
}

/**
 * @brief Set a nonlinear flow function for a mode.
 *
 * For general ODE ẋ = f(t, x), provides a callback that computes the
 * derivative. The system matrix A and offset b are ignored when using
 * a nonlinear flow.
 *
 * @param ha      The automaton
 * @param mode_id Mode index
 * @param f       Function computing dx/dt: f(t, x, dxdt, n, ctx)
 * @param ctx     Opaque context pointer
 * @return        true on success
 *
 * Complexity: O(1) — sets function pointers only
 */
bool hybrid_flow_set_nonlinear(HybridAutomaton *ha, int mode_id,
                                void (*f)(double, const double*, double*, int, void*),
                                void *ctx)
{
    if (!ha || mode_id < 0 || mode_id >= ha->num_modes) return false;
    if (!f) return false;

    HybridMode *m = &ha->modes[mode_id];
    m->flow.type = HAFLOW_NONLINEAR;
    m->flow.nonlinear_f = f;
    m->flow.nonlinear_ctx = ctx;
    m->flow.has_A = false;

    return true;
}

/* ==========================================================================
 * Initial condition (KP9)
 * ========================================================================== */

/**
 * @brief KP9: Set the initial condition as a single point.
 *
 * Init = {(q₀, x₀)}: the system starts in mode q₀ with continuous
 * state x₀. This is the most common initialization for deterministic
 * simulation.
 *
 * @param ha        The automaton
 * @param init_mode Initial mode index
 * @param x0        Initial continuous state (n-vector), NULL = zeros
 * @return          true on success
 *
 * Complexity: O(n) for copying
 */
bool hybrid_init_set(HybridAutomaton *ha, int init_mode, const double *x0)
{
    if (!ha || init_mode < 0 || init_mode >= ha->num_modes) return false;

    ha->init.init_mode = init_mode;
    ha->init.has_point = true;
    ha->init.init_has_bounds = false;

    if (x0) {
        vector_copy(ha->init.x0, x0, ha->num_vars);
    } else {
        for (int i = 0; i < ha->num_vars; i++) {
            ha->init.x0[i] = 0.0;
        }
    }

    /* Mark the mode as INITIAL type */
    ha->modes[init_mode].type = HAMODE_INITIAL;

    return true;
}

/**
 * @brief Set the initial condition as a rectangular region.
 *
 * Init ⊆ {q₀} × [lo, hi]: the system starts in mode q₀ with
 * continuous state anywhere in the axis-aligned box [lo_i, hi_i].
 *
 * @param ha        The automaton
 * @param init_mode Initial mode index
 * @param lo        Lower bounds (n-vector)
 * @param hi        Upper bounds (n-vector)
 * @return          true on success
 */
bool hybrid_init_set_rect(HybridAutomaton *ha, int init_mode,
                           const double *lo, const double *hi)
{
    if (!ha || init_mode < 0 || init_mode >= ha->num_modes) return false;
    if (!lo || !hi) return false;

    ha->init.init_mode = init_mode;
    ha->init.has_point = false;
    ha->init.init_has_bounds = true;

    vector_copy(ha->init.init_lo, lo, ha->num_vars);
    vector_copy(ha->init.init_hi, hi, ha->num_vars);

    ha->modes[init_mode].type = HAMODE_INITIAL;

    return true;
}

/* ==========================================================================
 * Query functions
 * ========================================================================== */

int hybrid_automaton_mode_count(const HybridAutomaton *ha)
{
    return ha ? ha->num_modes : 0;
}

int hybrid_automaton_var_count(const HybridAutomaton *ha)
{
    return ha ? ha->num_vars : 0;
}

int hybrid_automaton_trans_count(const HybridAutomaton *ha)
{
    return ha ? ha->num_transitions : 0;
}

const HybridMode* hybrid_automaton_get_mode(const HybridAutomaton *ha, int mode_id)
{
    if (!ha || mode_id < 0 || mode_id >= ha->num_modes) return NULL;
    return &ha->modes[mode_id];
}

const HybridTransition* hybrid_automaton_get_transition(const HybridAutomaton *ha, int trans_id)
{
    if (!ha || trans_id < 0 || trans_id >= ha->num_transitions) return NULL;
    return &ha->trans[trans_id];
}

const HybridVariable* hybrid_automaton_get_variable(const HybridAutomaton *ha, int var_idx)
{
    if (!ha || var_idx < 0 || var_idx >= ha->num_vars) return NULL;
    return &ha->vars[var_idx];
}

int hybrid_automaton_out_degree(const HybridAutomaton *ha, int mode_id)
{
    if (!ha || mode_id < 0 || mode_id >= ha->num_modes) return 0;
    return ha->trans_from_count[mode_id];
}

int hybrid_automaton_in_degree(const HybridAutomaton *ha, int mode_id)
{
    if (!ha || mode_id < 0 || mode_id >= ha->num_modes) return 0;
    return ha->trans_to_count[mode_id];
}

/* ==========================================================================
 * Well-formedness check
 * ========================================================================== */

/**
 * @brief Check whether the hybrid automaton is well-formed.
 *
 * Well-formedness conditions:
 *   1. At least one mode exists
 *   2. All transitions reference valid source and target modes
 *   3. The initial mode is set and valid
 *   4. No invariants are self-contradictory (trivially empty)
 *
 * This is a syntactic check; semantic properties like non-blocking
 * and determinism are checked elsewhere.
 *
 * @param ha The automaton
 * @return   true if well-formed
 *
 * Complexity: O(|Q| + |E|)
 */
bool hybrid_automaton_is_well_formed(const HybridAutomaton *ha)
{
    if (!ha) return false;

    /* Condition 1: At least one mode */
    if (ha->num_modes < 1) {
        fprintf(stderr, "[well-formed] FAIL: no modes defined\n");
        return false;
    }

    /* Condition 2: Valid transition endpoints */
    for (int i = 0; i < ha->num_transitions; i++) {
        const HybridTransition *t = &ha->trans[i];
        if (t->src_mode < 0 || t->src_mode >= ha->num_modes) {
            fprintf(stderr, "[well-formed] FAIL: transition %d has invalid src %d\n",
                    i, t->src_mode);
            return false;
        }
        if (t->tgt_mode < 0 || t->tgt_mode >= ha->num_modes) {
            fprintf(stderr, "[well-formed] FAIL: transition %d has invalid tgt %d\n",
                    i, t->tgt_mode);
            return false;
        }
    }

    /* Condition 3: Initial mode set and valid */
    if (ha->init.init_mode < 0 || ha->init.init_mode >= ha->num_modes) {
        fprintf(stderr, "[well-formed] FAIL: invalid or missing initial mode %d\n",
                ha->init.init_mode);
        return false;
    }

    /* Condition 4: No trivially contradictory invariants */
    /* An invariant with constraints H·x ≤ k is trivially empty if for some j,
       the constraint 0 ≤ -1 appears (H_row all zeros, k negative).
       We flag this but don't fail — it's theoretically valid (empty mode). */
    for (int m = 0; m < ha->num_modes; m++) {
        const HybridInvariant *inv = &ha->modes[m].invariant;
        if (inv->is_unbounded) continue;
        for (int c = 0; c < inv->num_constraints; c++) {
            bool all_zero = true;
            for (int j = 0; j < ha->num_vars; j++) {
                if (fabs(inv->H[c][j]) > 1e-12) { all_zero = false; break; }
            }
            if (all_zero && inv->k[c] < -1e-12) {
                fprintf(stderr, "[well-formed] WARNING: mode %d has "
                        "contradictory constraint (0 ≤ %g)\n", m, inv->k[c]);
            }
        }
    }

    return true;
}

/* ==========================================================================
 * Print / introspection
 * ========================================================================== */

/**
 * @brief Print a human-readable description of the automaton.
 *
 * Outputs to stdout the full specification of the hybrid automaton:
 * name, modes with their invariants and flows, transitions with
 * their guards and resets, and initial condition.
 *
 * @param ha The automaton to print
 */
void hybrid_automaton_print(const HybridAutomaton *ha)
{
    if (!ha) {
        printf("NULL automaton\n");
        return;
    }

    printf("============================================\n");
    printf("Hybrid Automaton: %s\n", ha->name);
    printf("============================================\n");
    printf("Continuous variables: %d\n", ha->num_vars);
    printf("Discrete modes:       %d\n", ha->num_modes);
    printf("Transitions:          %d\n", ha->num_transitions);
    printf("\n");

    /* Print variables */
    printf("--- Variables ---\n");
    for (int i = 0; i < ha->num_vars; i++) {
        const HybridVariable *v = &ha->vars[i];
        const char *type_str = "REAL";
        switch (v->type) {
        case HAVAR_REAL:         type_str = "REAL"; break;
        case HAVAR_REAL_BOUNDED: type_str = "REAL_BOUNDED"; break;
        case HAVAR_INTEGER:      type_str = "INTEGER"; break;
        case HAVAR_CLOCK:        type_str = "CLOCK"; break;
        }
        printf("  x[%d]: %s (%s) [%g, %g] rate=%g\n",
               i, v->name[0] ? v->name : "(unnamed)", type_str,
               v->lo, v->hi, v->rate);
    }
    printf("\n");

    /* Print modes */
    printf("--- Modes ---\n");
    for (int i = 0; i < ha->num_modes; i++) {
        const HybridMode *m = &ha->modes[i];
        const char *type_str = "NORMAL";
        switch (m->type) {
        case HAMODE_INITIAL:  type_str = "INITIAL"; break;
        case HAMODE_TERMINAL: type_str = "TERMINAL"; break;
        case HAMODE_URGENT:   type_str = "URGENT"; break;
        default: break;
        }
        printf("  Mode %d: %s (%s)\n", i, m->name[0] ? m->name : "(unnamed)", type_str);

        /* Invariant */
        if (m->invariant.is_unbounded) {
            printf("    Inv: unbounded (R^n)\n");
        } else {
            printf("    Inv: %d constraints\n", m->invariant.num_constraints);
            for (int c = 0; c < m->invariant.num_constraints; c++) {
                printf("      ");
                for (int j = 0; j < ha->num_vars; j++) {
                    if (fabs(m->invariant.H[c][j]) > 1e-12) {
                        printf("%+.2g*x%d ", m->invariant.H[c][j], j);
                    }
                }
                printf("<= %g\n", m->invariant.k[c]);
            }
        }

        /* Flow */
        const char *flow_str = "CONSTANT";
        switch (m->flow.type) {
        case HAFLOW_CONSTANT:  flow_str = "CONSTANT"; break;
        case HAFLOW_LINEAR:    flow_str = "LINEAR"; break;
        case HAFLOW_AFFINE:    flow_str = "AFFINE"; break;
        case HAFLOW_NONLINEAR: flow_str = "NONLINEAR"; break;
        case HAFLOW_DIFF_INCL: flow_str = "DIFF_INCLUSION"; break;
        }
        printf("    Flow: %s", flow_str);
        if (m->flow.has_A) {
            printf(" A=[%dx%d]", ha->num_vars, ha->num_vars);
        }
        printf("\n");

        /* Outgoing transitions */
        int out = ha->trans_from_count[i];
        if (out > 0) {
            printf("    Outgoing: %d transitions", out);
            for (int t = 0; t < out; t++) {
                printf(" ->%d", ha->modes[ha->trans[ha->trans_from[i][t]].tgt_mode].id);
            }
            printf("\n");
        }
    }
    printf("\n");

    /* Print transitions */
    printf("--- Transitions ---\n");
    for (int i = 0; i < ha->num_transitions; i++) {
        const HybridTransition *t = &ha->trans[i];
        printf("  T%d: %s (%d -> %d)", i,
               t->name[0] ? t->name : "(unnamed)",
               t->src_mode, t->tgt_mode);

        if (t->guard.is_trivially_true) {
            printf(" guard:TRUE");
        } else {
            printf(" guard:%d constr", t->guard.num_constraints);
        }

        const char *rst_str = "ID";
        switch (t->reset.type) {
        case HARESET_IDENTITY:   rst_str = "ID"; break;
        case HARESET_CONSTANT:   rst_str = "CONST"; break;
        case HARESET_AFFINE:     rst_str = "AFFINE"; break;
        case HARESET_CLOCK_RESET: rst_str = "CLK_RST"; break;
        }
        printf(" reset:%s", rst_str);

        if (!t->is_active) printf(" [INACTIVE]");
        printf("\n");
    }
    printf("\n");

    /* Print initial condition */
    printf("--- Initial Condition ---\n");
    if (ha->init.init_mode >= 0) {
        printf("  Mode: %d (%s)\n", ha->init.init_mode,
               ha->modes[ha->init.init_mode].name);
    }
    if (ha->init.has_point) {
        printf("  State: [");
        for (int i = 0; i < ha->num_vars; i++) {
            printf("%g%s", ha->init.x0[i], i < ha->num_vars - 1 ? ", " : "");
        }
        printf("]\n");
    } else if (ha->init.init_has_bounds) {
        printf("  Region: [");
        for (int i = 0; i < ha->num_vars; i++) {
            printf("%g", ha->init.init_lo[i]);
            if (i < ha->num_vars - 1) printf(", ");
        }
        printf("] -> [");
        for (int i = 0; i < ha->num_vars; i++) {
            printf("%g", ha->init.init_hi[i]);
            if (i < ha->num_vars - 1) printf(", ");
        }
        printf("]\n");
    }
    printf("============================================\n");
}
