/**
 * @file hybrid_execution.c
 * @brief Execution semantics implementation (KP10-KP19, L2)
 *
 * Implements execution construction, semantic predicates (determinism,
 * non-blocking, Zeno detection), guard/invariant evaluation, reset
 * application, and parallel composition of hybrid automata.
 *
 * Reference:
 *   Lygeros et al., "Dynamical Properties of Hybrid Automata" (2003)
 *   Goebel, Sanfelice, Teel, "Hybrid Dynamical Systems" (2009)
 *   Lynch, Segala, Vaandrager, "Hybrid I/O Automata" (2003)
 */

#include "hybrid_execution.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==========================================================================
 * Execution construction (KP11-KP13)
 * ========================================================================== */

/**
 * @brief KP11: Create a new execution for a given automaton.
 *
 * An execution χ is a sequence of alternating continuous flow segments
 * and discrete jumps, starting from the automaton's initial state.
 *
 * Execution χ = τ₀ → (jump₀) → τ₁ → (jump₁) → τ₂ → ...
 *
 * where τ_i = (q_i, [t_i, t_i'], ξ_i) is a flow segment in mode q_i
 * with trajectory ξ_i, and jump_i is an instantaneous transition.
 *
 * @param ha                The automaton to execute
 * @param max_flow_segments Initial capacity for flow segments
 * @param max_jumps         Initial capacity for jumps
 * @return                  New execution, or NULL on error
 *
 * Complexity: O(max_segments + max_jumps) allocation
 */
HybridExecution* hybrid_execution_create(const HybridAutomaton *ha,
                                          int max_flow_segments, int max_jumps)
{
    if (!ha || max_flow_segments < 1 || max_jumps < 0) return NULL;

    HybridExecution *exec = (HybridExecution*) calloc(1, sizeof(HybridExecution));
    if (!exec) return NULL;

    /* Cast away const — the execution doesn't modify the automaton,
       but we need a non-const reference for compatibility */
    exec->automaton = (HybridAutomaton*) ha;
    exec->num_segments = 0;
    exec->num_jumps = 0;
    exec->max_segments = max_flow_segments;
    exec->max_jumps = max_jumps;
    exec->is_finite = false;
    exec->is_zeno = false;
    exec->total_time = 0.0;

    exec->flow_segments = (HybridFlowSegment*) calloc(
        (size_t)max_flow_segments, sizeof(HybridFlowSegment));
    exec->jump_points = (HybridJumpPoint*) calloc(
        (size_t)max_jumps, sizeof(HybridJumpPoint));

    if (!exec->flow_segments || !exec->jump_points) {
        free(exec->flow_segments);
        free(exec->jump_points);
        free(exec);
        return NULL;
    }

    return exec;
}

/**
 * @brief Destroy an execution and free resources.
 */
void hybrid_execution_destroy(HybridExecution *exec)
{
    if (!exec) return;
    free(exec->flow_segments);
    free(exec->jump_points);
    free(exec);
}

/**
 * @brief KP12: Append a continuous flow segment to the execution.
 *
 * Records that the system evolved continuously in a specific mode
 * over time interval [t_start, t_end], reaching x_end from x_start.
 *
 * The ODE solution connecting x_start to x_end must satisfy:
 *   ẋ(t) = f_q(x(t)) for all t ∈ [t_start, t_end]
 *   x(t) ∈ Inv(q) for all t ∈ [t_start, t_end]
 *
 * @param exec       The execution
 * @param mode_index Mode during flow
 * @param t_start    Start time
 * @param t_end      End time (must be ≥ t_start)
 * @param x_start    State at t_start
 * @param x_end      State at t_end
 * @param num_steps  Number of integration steps
 * @return           true on success
 *
 * Complexity: O(1) append, O(n) copy
 */
bool hybrid_execution_append_flow(HybridExecution *exec, int mode_index,
                                   double t_start, double t_end,
                                   const double *x_start, const double *x_end,
                                   int num_steps)
{
    if (!exec || !x_start || !x_end) return false;
    if (t_end < t_start) return false;
    if (mode_index < 0 || mode_index >= exec->automaton->num_modes) return false;

    /* Expand capacity if needed */
    if (exec->num_segments >= exec->max_segments) {
        int new_max = exec->max_segments * 2;
        HybridFlowSegment *new_segs = (HybridFlowSegment*) realloc(
            exec->flow_segments, (size_t)new_max * sizeof(HybridFlowSegment));
        if (!new_segs) return false;
        exec->flow_segments = new_segs;
        exec->max_segments = new_max;
    }

    HybridFlowSegment *seg = &exec->flow_segments[exec->num_segments];
    seg->mode_index = mode_index;
    seg->t_start = t_start;
    seg->t_end = t_end;
    seg->num_steps = num_steps;

    int n = exec->automaton->num_vars;
    memcpy(seg->x_start, x_start, (size_t)n * sizeof(double));
    memcpy(seg->x_end, x_end, (size_t)n * sizeof(double));

    exec->num_segments++;
    exec->total_time += (t_end - t_start);

    return true;
}

/**
 * @brief KP13: Append a discrete jump to the execution.
 *
 * Records an instantaneous mode switch at time t_jump.
 * The transition fires when its guard is satisfied at x_pre,
 * and the continuous state is reset to x_post = R(e)·x_pre + r(e).
 *
 * @param exec           The execution
 * @param transition_id  Index of the transition
 * @param t_jump         Time of the jump
 * @param x_pre          State before jump
 * @param x_post         State after reset
 * @return               true on success
 *
 * Complexity: O(1) append, O(n) copy
 */
bool hybrid_execution_append_jump(HybridExecution *exec, int transition_id,
                                   double t_jump,
                                   const double *x_pre, const double *x_post)
{
    if (!exec || !x_pre || !x_post) return false;
    if (transition_id < 0 || transition_id >= exec->automaton->num_transitions)
        return false;

    /* Expand capacity if needed */
    if (exec->num_jumps >= exec->max_jumps) {
        int new_max = exec->max_jumps * 2;
        HybridJumpPoint *new_jumps = (HybridJumpPoint*) realloc(
            exec->jump_points, (size_t)new_max * sizeof(HybridJumpPoint));
        if (!new_jumps) return false;
        exec->jump_points = new_jumps;
        exec->max_jumps = new_max;
    }

    HybridJumpPoint *jp = &exec->jump_points[exec->num_jumps];
    const HybridTransition *t = &exec->automaton->trans[transition_id];

    jp->transition_id = transition_id;
    jp->src_mode = t->src_mode;
    jp->tgt_mode = t->tgt_mode;
    jp->t_jump = t_jump;

    int n = exec->automaton->num_vars;
    memcpy(jp->x_pre, x_pre, (size_t)n * sizeof(double));
    memcpy(jp->x_post, x_post, (size_t)n * sizeof(double));

    exec->num_jumps++;

    return true;
}

/* ==========================================================================
 * KP14: Determinism check
 * ========================================================================== */

/**
 * @brief KP14: Check if the hybrid automaton is deterministic.
 *
 * Determinism conditions for hybrid automata:
 *   1. Initial condition is a single point (not a set)
 *   2. For each mode q, Flow(q) is a function (not differential inclusion)
 *   3. For any state (q, x), at most one outgoing transition guard is true
 *      (guard regions from same mode must be mutually exclusive)
 *
 * Non-determinism sources:
 *   - Guards overlap: ∃e₁≠e₂ from q s.t. G(e₁) ∩ G(e₂) ≠ ∅
 *   - Differential inclusions (set-valued flows)
 *   - Set-based initial conditions
 *
 * @param ha The automaton
 * @return   true if deterministic
 *
 * Complexity: O(|E|² · n) for guard overlap checking
 * Theorem: Determinism is semi-decidable for linear HA (Henzinger et al., 1995)
 */
bool hybrid_is_deterministic(const HybridAutomaton *ha)
{
    if (!ha) return false;

    /* Condition 1: Point initial condition (or not set) */
    /* We consider point initial as deterministic; set-based as non-deterministic */
    if (!ha->init.has_point) return false;

    /* Condition 2: No differential inclusions (each mode has a single flow) */
    for (int m = 0; m < ha->num_modes; m++) {
        if (ha->modes[m].flow.type == HAFLOW_DIFF_INCL) {
            return false;
        }
    }

    /* Condition 3: Mutual exclusion of guard regions from same mode */
    /* For each pair of distinct transitions from the same source mode,
       check if their guards could simultaneously be true. This is
       computationally hard in general; we use a simple check: if both
       guards are trivially true, they overlap. */
    for (int q = 0; q < ha->num_modes; q++) {
        int out_count = ha->trans_from_count[q];
        if (out_count <= 1) continue;

        for (int i = 0; i < out_count; i++) {
            for (int j = i + 1; j < out_count; j++) {
                int ti = ha->trans_from[q][i];
                int tj = ha->trans_from[q][j];
                const HybridGuard *gi = &ha->trans[ti].guard;
                const HybridGuard *gj = &ha->trans[tj].guard;

                /* If both are trivially true, they overlap (non-det) */
                if (gi->is_trivially_true && gj->is_trivially_true) {
                    return false;
                }
            }
        }
    }

    return true;
}

/* ==========================================================================
 * KP15: Non-blocking check
 * ========================================================================== */

/**
 * @brief KP15: Check if the automaton is non-blocking (deadlock-free).
 *
 * Deadlock occurs when the system reaches a state (q, x) where:
 *   - x is on the boundary of Inv(q) (cannot flow forward)
 *   - No guard is satisfied (cannot jump)
 *
 * This function checks the structural property: for every mode,
 * either the invariant is unbounded (can always flow) or there
 * exists at least one outgoing transition with a guard that
 * intersects the invariant boundary.
 *
 * @param ha The automaton
 * @return   true if structurally non-blocking
 */
bool hybrid_is_nonblocking(const HybridAutomaton *ha)
{
    if (!ha) return false;

    for (int q = 0; q < ha->num_modes; q++) {
        const HybridMode *mode = &ha->modes[q];

        /* Mode with unbounded invariant can flow forever (non-blocking) */
        if (mode->invariant.is_unbounded) continue;

        /* Mode with no outgoing transitions — check if flow can stop
           before hitting invariant boundary */
        int out_degree = ha->trans_from_count[q];
        if (out_degree == 0) {
            /* Terminal mode — acceptable as non-blocking if it's marked TERMINAL */
            if (mode->type != HAMODE_TERMINAL) {
                /* Check if all flows point away from boundary — simplified */
                /* If flow is zero (steady state), system might be stuck at boundary */
                bool all_zero_flow = true;
                for (int v = 0; v < ha->num_vars; v++) {
                    if (fabs(mode->flow.b[v]) > 1e-12) { all_zero_flow = false; break; }
                }
                if (all_zero_flow && !mode->flow.has_A) {
                    return false; /* Deadlock: stuck in mode with invariant, no exit */
                }
            }
            continue;
        }

        /* Mode with outgoing transitions — check at least one is trivially
           true or guards collectively cover the invariant boundary.
           Simplified: if any outgoing transition has trivially true guard,
           the mode is non-blocking. */
        bool has_escape = false;
        for (int t = 0; t < out_degree; t++) {
            int tid = ha->trans_from[q][t];
            if (ha->trans[tid].guard.is_trivially_true) {
                has_escape = true;
                break;
            }
        }

        if (!has_escape) {
            /* Conservative: we cannot guarantee non-blocking */
            /* In a full implementation, we'd check coverage of Inv boundary
               by the union of all guard regions. Here we flag as potentially
               blocking. */
            return false;
        }
    }

    return true;
}

/* ==========================================================================
 * KP16: Zeno detection
 * ========================================================================== */

/**
 * @brief KP16: Check for Zeno behavior in an execution.
 *
 * A Zeno execution has infinitely many discrete transitions
 * accumulating in finite total time. Formally:
 *   ∃T_∞ < ∞ such that lim_{j→∞} t_j = T_∞
 *
 * Detection heuristics:
 *   - If last N transitions have very short inter-transition times
 *   - If the time interval between consecutive jumps is monotonically
 *     decreasing and approaches zero
 *
 * @param exec The execution
 * @return     true if Zeno behavior detected
 *
 * Theorem: Zeno detection is undecidable for general HA
 *          but semi-decidable for timed/rectangular subclasses
 */
bool hybrid_execution_is_zeno(const HybridExecution *exec)
{
    if (!exec || exec->num_jumps < 10) return false;

    /* Check trend of inter-jump times over the last 10 jumps */
    int check_window = 10;
    if (exec->num_jumps < check_window) check_window = exec->num_jumps;
    if (check_window < 3) return false;

    int start = exec->num_jumps - check_window;
    double prev_dt = INFINITY;
    int decreasing_count = 0;
    int very_small_count = 0;

    for (int i = start + 1; i < exec->num_jumps; i++) {
        double dt = exec->jump_points[i].t_jump -
                     exec->jump_points[i - 1].t_jump;
        if (dt < prev_dt * 0.95) decreasing_count++;
        if (dt < 1e-9) very_small_count++;
        prev_dt = dt;
    }

    /* Zeno indicators */
    if (decreasing_count >= check_window - 2) {
        /* Monotonically decreasing jump intervals → Zeno */
        return true;
    }

    if (very_small_count >= check_window / 2) {
        /* Many near-zero jump intervals → Zeno */
        return true;
    }

    return false;
}

/**
 * @brief Estimate the Zeno accumulation time.
 *
 * Uses a geometric series approximation: if inter-jump times
 * decrease geometrically with ratio r < 1, then
 * T_∞ = t_N + d_N / (1 - r) where d_N is the last inter-jump time.
 *
 * @param exec The execution
 * @return     Estimated Zeno time, or INFINITY if not Zeno
 */
double hybrid_execution_zeno_time(const HybridExecution *exec)
{
    if (!exec || exec->num_jumps < 3) return INFINITY;

    /* Estimate the geometric ratio from last 3 inter-jump times */
    int nj = exec->num_jumps;
    double d1 = exec->jump_points[nj - 2].t_jump -
                 exec->jump_points[nj - 3].t_jump;
    double d2 = exec->jump_points[nj - 1].t_jump -
                 exec->jump_points[nj - 2].t_jump;

    if (d1 <= 1e-15 || fabs(d1 - d2) < 1e-15) return INFINITY;

    double ratio = d2 / d1;
    if (ratio >= 1.0 || ratio <= 0.0) return INFINITY;

    double t_last = exec->jump_points[nj - 1].t_jump;
    if (fabs(1.0 - ratio) < 1e-15) return INFINITY;

    return t_last + d2 / (1.0 - ratio);
}

/* ==========================================================================
 * KP17: Transition enabled check
 * ========================================================================== */

/**
 * @brief KP17: Check if a transition is enabled at (mode, x).
 *
 * A transition e = (q → q', G(e), R(e)) is enabled if:
 *   1. Current mode matches q (unless mode = -1 for "any mode")
 *   2. Guard G(e) is satisfied by x: A·x ≤ b (∀ constraints)
 *   3. After reset, the resulting state is in Inv(q')
 *
 * @param ha        The automaton
 * @param trans_id  Transition index
 * @param x         Current state
 * @param mode      Current mode (-1 to skip mode check)
 * @return          true if transition can fire
 *
 * Complexity: O(n · max(num_constraints, n)) for guard + invariant check
 */
bool hybrid_transition_enabled(const HybridAutomaton *ha, int trans_id,
                                const double *x, int mode)
{
    if (!ha || !x || trans_id < 0 || trans_id >= ha->num_transitions)
        return false;

    const HybridTransition *t = &ha->trans[trans_id];
    int n = ha->num_vars;

    /* Check 1: mode match */
    if (mode >= 0 && t->src_mode != mode) return false;

    /* Check 2: guard satisfaction */
    if (!t->guard.is_trivially_true) {
        for (int c = 0; c < t->guard.num_constraints; c++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += t->guard.A[c][j] * x[j];
            }
            if (sum > t->guard.b[c] + 1e-12) {
                return false; /* Constraint violated */
            }
        }
    }

    /* Check 3: post-reset state satisfies target invariant */
    const HybridInvariant *inv_tgt = &ha->modes[t->tgt_mode].invariant;
    if (!inv_tgt->is_unbounded) {
        /* Compute x' = R·x + r */
        double x_post[HA_MAX_VARIABLES];
        hybrid_reset_apply(&t->reset, x, x_post, n);

        /* Check invariant */
        if (!hybrid_invariant_satisfied(inv_tgt, x_post, n)) {
            return false;
        }
    }

    return true;
}

/* ==========================================================================
 * Guard evaluation
 * ========================================================================== */

/**
 * @brief Evaluate a guard condition at state x.
 *
 * Returns true if and only if A·x ≤ b for all constraints.
 * Trivially true guards always return true.
 *
 * @param guard The guard to evaluate
 * @param x     Continuous state
 * @param n     State dimension
 * @return      true if x satisfies the guard
 */
bool hybrid_guard_evaluate(const HybridGuard *guard, const double *x, int n)
{
    if (!guard) return true;
    if (guard->is_trivially_true) return true;
    if (!x || n <= 0) return false;

    for (int c = 0; c < guard->num_constraints; c++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += guard->A[c][j] * x[j];
        }
        if (sum > guard->b[c] + 1e-12) {
            return false;
        }
    }
    return true;
}

/* ==========================================================================
 * Reset application
 * ========================================================================== */

/**
 * @brief Apply a reset map to a state.
 *
 * Computes x_out = R(e)·x + r(e).
 * x and x_out may be the same array (in-place update).
 *
 * @param reset The reset map
 * @param x     Input state
 * @param x_out Output state
 * @param n     Dimension
 *
 * Complexity: O(n²) for affine, O(n) for identity/constant
 */
void hybrid_reset_apply(const HybridReset *reset, const double *x,
                         double *x_out, int n)
{
    if (!reset || !x || !x_out || n <= 0) return;

    /* Use a temporary to handle in-place aliasing */
    double tmp[HA_MAX_VARIABLES];

    switch (reset->type) {
    case HARESET_IDENTITY:
        /* x_out = x */
        if (x_out != x) {
            for (int i = 0; i < n; i++) x_out[i] = x[i];
        }
        break;

    case HARESET_CONSTANT:
        /* x_out = r */
        for (int i = 0; i < n; i++) x_out[i] = reset->r[i];
        break;

    case HARESET_AFFINE:
        /* x_out = R·x + r */
        for (int i = 0; i < n; i++) {
            double sum = reset->r[i];
            for (int j = 0; j < n; j++) {
                sum += reset->R[i][j] * x[j];
            }
            tmp[i] = sum;
        }
        for (int i = 0; i < n; i++) x_out[i] = tmp[i];
        break;

    case HARESET_CLOCK_RESET:
        /* x_out = R·x + r, where R is modified identity
           (some diagonal entries zeroed for clock reset) */
        for (int i = 0; i < n; i++) {
            double sum = reset->r[i];
            for (int j = 0; j < n; j++) {
                sum += reset->R[i][j] * x[j];
            }
            tmp[i] = sum;
        }
        for (int i = 0; i < n; i++) x_out[i] = tmp[i];
        break;

    default:
        /* Fallback: identity */
        if (x_out != x) {
            for (int i = 0; i < n; i++) x_out[i] = x[i];
        }
        break;
    }
}

/* ==========================================================================
 * KP18: Invariant satisfaction
 * ========================================================================== */

/**
 * @brief KP18: Check if state x satisfies an invariant condition.
 *
 * Returns true if H·x ≤ k for all constraints (or if invariant is unbounded).
 *
 * @param invariant The invariant to check
 * @param x         State vector
 * @param n         Dimension
 * @return          true if x ∈ Inv
 *
 * Complexity: O(num_constraints · n)
 */
bool hybrid_invariant_satisfied(const HybridInvariant *invariant,
                                 const double *x, int n)
{
    if (!invariant) return false;
    if (invariant->is_unbounded) return true;
    if (!x || n <= 0) return false;

    for (int c = 0; c < invariant->num_constraints; c++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += invariant->H[c][j] * x[j];
        }
        if (sum > invariant->k[c] + 1e-12) {
            return false;
        }
    }
    return true;
}

/* ==========================================================================
 * Find all enabled transitions
 * ========================================================================== */

/**
 * @brief Find all transitions from a mode that are enabled at state x.
 *
 * @param ha          The automaton
 * @param mode        Current mode
 * @param x           Current state
 * @param enabled_ids Output array
 * @param max_count   Max entries
 * @return            Number of enabled transitions
 */
int hybrid_transitions_enabled(const HybridAutomaton *ha, int mode,
                                const double *x, int *enabled_ids, int max_count)
{
    if (!ha || !x || !enabled_ids || max_count <= 0) return 0;
    if (mode < 0 || mode >= ha->num_modes) return 0;

    int count = 0;
    int out_degree = ha->trans_from_count[mode];

    for (int i = 0; i < out_degree && count < max_count; i++) {
        int tid = ha->trans_from[mode][i];
        if (!ha->trans[tid].is_active) continue;

        if (hybrid_transition_enabled(ha, tid, x, mode)) {
            enabled_ids[count++] = tid;
        }
    }

    return count;
}

/* ==========================================================================
 * KP19: Parallel composition
 * ========================================================================== */

/**
 * @brief KP19: Parallel composition H₁ || H₂.
 *
 * Two hybrid automata are composed via interleaving on continuous flow
 * and synchronization on shared event labels.
 *
 * Composition rules (from Lynch, Segala, Vaandrager 2003):
 *   - Q = Q₁ × Q₂ (Cartesian product of mode sets)
 *   - X = X₁ ∪ X₂ (union of variables)
 *   - Init = Init₁ × Init₂
 *   - Flow: for mode (q₁, q₂), the flow is the product flow:
 *           ẋ₁ = f_{q₁}(x₁), ẋ₂ = f_{q₂}(x₂)
 *   - Invariant: Inv(q₁, q₂) = Inv₁(q₁) × Inv₂(q₂)
 *   - Transitions:
 *     * Async move of H₁: ((q₁,q₂) → (q₁',q₂)) if q₁ → q₁' in E₁
 *       and guard and reset only affect X₁ variables
 *     * Async move of H₂: ((q₁,q₂) → (q₁,q₂')) if q₂ → q₂' in E₂
 *     * Synchronization on shared events (not implemented — simplified)
 *
 * @param h1 First automaton
 * @param h2 Second automaton
 * @return   Composed automaton, or NULL on error
 *
 * Complexity: O(|Q₁|·|Q₂| + |E₁|·|Q₂| + |Q₁|·|E₂|)
 */
HybridAutomaton* hybrid_compose_parallel(const HybridAutomaton *h1,
                                          const HybridAutomaton *h2)
{
    if (!h1 || !h2) return NULL;

    /* Check dimensions: variables should be independent (union) */
    int n1 = h1->num_vars;
    int n2 = h2->num_vars;
    int n_total = n1 + n2;

    if (n_total > HA_MAX_VARIABLES) {
        fprintf(stderr, "[compose] Total variables %d exceed max %d\n",
                n_total, HA_MAX_VARIABLES);
        return NULL;
    }

    /* Product modes: |Q| = |Q₁| × |Q₂| */
    int num_modes = h1->num_modes * h2->num_modes;
    if (num_modes > HA_MAX_MODES) {
        fprintf(stderr, "[compose] Product modes %d exceed max %d\n",
                num_modes, HA_MAX_MODES);
        return NULL;
    }

    /* Create the composed automaton */
    char comp_name[HA_NAME_LEN];
    snprintf(comp_name, HA_NAME_LEN, "%s||%s", h1->name, h2->name);

    HybridAutomaton *h = hybrid_automaton_create(comp_name, n_total);
    if (!h) return NULL;

    /* Create product modes */
    for (int q1 = 0; q1 < h1->num_modes; q1++) {
        for (int q2 = 0; q2 < h2->num_modes; q2++) {
            char mname[HA_NAME_LEN];
            snprintf(mname, HA_NAME_LEN, "(%s,%s)",
                     h1->modes[q1].name[0] ? h1->modes[q1].name : "?",
                     h2->modes[q2].name[0] ? h2->modes[q2].name : "?");
            HybridModeType mtype = HAMODE_NORMAL;
            if (h1->modes[q1].type == HAMODE_INITIAL &&
                h2->modes[q2].type == HAMODE_INITIAL) {
                mtype = HAMODE_INITIAL;
            }
            int mid = hybrid_automaton_add_mode(h, mname, mtype);
            if (mid < 0) { hybrid_automaton_destroy(h); return NULL; }

            /* Combine invariants: we extend H1's constraints to n_total dim
               and append H2's constraints (offset to X₂). Since we're using
               fixed-size arrays, we place the product invariant manually. */

            /* Copy H1 invariant (extended with zeros for X₂ dimensions) */
            int total_c = 0;
            const HybridInvariant *inv1 = &h1->modes[q1].invariant;
            if (!inv1->is_unbounded) {
                for (int c = 0; c < inv1->num_constraints && total_c < HA_MAX_CONSTRAINTS; c++) {
                    for (int j = 0; j < n1; j++) {
                        h->modes[mid].invariant.H[total_c][j] = inv1->H[c][j];
                    }
                    for (int j = n1; j < n_total; j++) {
                        h->modes[mid].invariant.H[total_c][j] = 0.0;
                    }
                    h->modes[mid].invariant.k[total_c] = inv1->k[c];
                    total_c++;
                }
            }
            /* Copy H2 invariant (offset to X₂ position) */
            const HybridInvariant *inv2 = &h2->modes[q2].invariant;
            if (!inv2->is_unbounded) {
                for (int c = 0; c < inv2->num_constraints && total_c < HA_MAX_CONSTRAINTS; c++) {
                    for (int j = 0; j < n1; j++) {
                        h->modes[mid].invariant.H[total_c][j] = 0.0;
                    }
                    for (int j = 0; j < n2; j++) {
                        h->modes[mid].invariant.H[total_c][n1 + j] = inv2->H[c][j];
                    }
                    h->modes[mid].invariant.k[total_c] = inv2->k[c];
                    total_c++;
                }
            }
            if (total_c > 0) {
                h->modes[mid].invariant.num_constraints = total_c;
                h->modes[mid].invariant.is_unbounded = false;
            }

            /* Copy combined flows: place A1 in upper-left, A2 in lower-right */
            if (h1->modes[q1].flow.has_A || h2->modes[q2].flow.has_A) {
                double A_comp[HA_MAX_VARIABLES * HA_MAX_VARIABLES] = {0};
                /* Upper-left: A1 */
                for (int i = 0; i < n1; i++)
                    for (int j = 0; j < n1; j++)
                        A_comp[i * n_total + j] = h1->modes[q1].flow.A[i][j];
                /* Lower-right: A2 */
                for (int i = 0; i < n2; i++)
                    for (int j = 0; j < n2; j++)
                        A_comp[(n1 + i) * n_total + (n1 + j)] = h2->modes[q2].flow.A[i][j];

                HybridFlowType ftype = (h1->modes[q1].flow.type == HAFLOW_AFFINE ||
                                         h2->modes[q2].flow.type == HAFLOW_AFFINE)
                                        ? HAFLOW_AFFINE : HAFLOW_LINEAR;
                double b_comp[HA_MAX_VARIABLES] = {0};
                for (int i = 0; i < n1; i++) b_comp[i] = h1->modes[q1].flow.b[i];
                for (int i = 0; i < n2; i++) b_comp[n1 + i] = h2->modes[q2].flow.b[i];

                hybrid_flow_set(h, mid, ftype, A_comp, b_comp);
            }
        }
    }

    /* Create async transitions for H1's edges */
    for (int t1 = 0; t1 < h1->num_transitions; t1++) {
        const HybridTransition *tr1 = &h1->trans[t1];
        for (int q2 = 0; q2 < h2->num_modes; q2++) {
            int src_prod = tr1->src_mode * h2->num_modes + q2;
            int tgt_prod = tr1->tgt_mode * h2->num_modes + q2;
            int tid = hybrid_automaton_add_transition(h, src_prod, tgt_prod,
                                                       tr1->type, tr1->name);
            if (tid < 0) continue;

            /* Extend guard to n_total dimensions (X₂ part trivially true) */
            if (!tr1->guard.is_trivially_true) {
                double A_ext[HA_MAX_CONSTRAINTS * HA_MAX_VARIABLES] = {0};
                double b_ext[HA_MAX_CONSTRAINTS];
                for (int c = 0; c < tr1->guard.num_constraints; c++) {
                    for (int j = 0; j < n1; j++)
                        A_ext[c * n_total + j] = tr1->guard.A[c][j];
                    b_ext[c] = tr1->guard.b[c];
                }
                hybrid_guard_set(h, tid, tr1->guard.num_constraints, A_ext, b_ext);
            }

            /* Extend reset: identity on X₂ */
            if (tr1->reset.type != HARESET_IDENTITY) {
                double R_ext[HA_MAX_VARIABLES * HA_MAX_VARIABLES];
                double r_ext[HA_MAX_VARIABLES] = {0};
                for (int i = 0; i < n_total; i++)
                    for (int j = 0; j < n_total; j++)
                        R_ext[i * n_total + j] = (i == j) ? 1.0 : 0.0;
                for (int i = 0; i < n1; i++) {
                    for (int j = 0; j < n1; j++)
                        R_ext[i * n_total + j] = tr1->reset.R[i][j];
                    r_ext[i] = tr1->reset.r[i];
                }
                hybrid_reset_set(h, tid, tr1->reset.type, R_ext, r_ext);
            }
        }
    }

    /* Create async transitions for H2's edges */
    for (int t2 = 0; t2 < h2->num_transitions; t2++) {
        const HybridTransition *tr2 = &h2->trans[t2];
        for (int q1 = 0; q1 < h1->num_modes; q1++) {
            int src_prod = q1 * h2->num_modes + tr2->src_mode;
            int tgt_prod = q1 * h2->num_modes + tr2->tgt_mode;
            int tid = hybrid_automaton_add_transition(h, src_prod, tgt_prod,
                                                       tr2->type, tr2->name);
            if (tid < 0) continue;

            if (!tr2->guard.is_trivially_true) {
                double A_ext[HA_MAX_CONSTRAINTS * HA_MAX_VARIABLES] = {0};
                double b_ext[HA_MAX_CONSTRAINTS];
                for (int c = 0; c < tr2->guard.num_constraints; c++) {
                    for (int j = 0; j < n2; j++)
                        A_ext[c * n_total + (n1 + j)] = tr2->guard.A[c][j];
                    b_ext[c] = tr2->guard.b[c];
                }
                hybrid_guard_set(h, tid, tr2->guard.num_constraints, A_ext, b_ext);
            }

            if (tr2->reset.type != HARESET_IDENTITY) {
                double R_ext[HA_MAX_VARIABLES * HA_MAX_VARIABLES];
                double r_ext[HA_MAX_VARIABLES] = {0};
                for (int i = 0; i < n_total; i++)
                    for (int j = 0; j < n_total; j++)
                        R_ext[i * n_total + j] = (i == j) ? 1.0 : 0.0;
                for (int i = 0; i < n2; i++) {
                    for (int j = 0; j < n2; j++)
                        R_ext[(n1 + i) * n_total + (n1 + j)] = tr2->reset.R[i][j];
                    r_ext[n1 + i] = tr2->reset.r[i];
                }
                hybrid_reset_set(h, tid, tr2->reset.type, R_ext, r_ext);
            }
        }
    }

    /* Set initial condition: product of individual initial states */
    int q1_init = h1->init.init_mode;
    int q2_init = h2->init.init_mode;
    if (q1_init >= 0 && q2_init >= 0) {
        int init_prod = q1_init * h2->num_modes + q2_init;
        double x0_comp[HA_MAX_VARIABLES] = {0};
        for (int i = 0; i < n1; i++) x0_comp[i] = h1->init.x0[i];
        for (int i = 0; i < n2; i++) x0_comp[n1 + i] = h2->init.x0[i];
        hybrid_init_set(h, init_prod, x0_comp);
    }

    return h;
}
