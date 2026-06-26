/**
 * @file hybrid_execution.h
 * @brief Execution semantics for hybrid automata (L2 Core Concepts)
 *
 * Defines trajectories, executions, hybrid time domains, and the
 * formal semantics of hybrid automaton evolution. Covers the continuous
 * flow within modes and discrete jumps between modes.
 *
 * Reference:
 *   Lygeros, Johansson, Simic, Zhang, Sastry, "Dynamical Properties of
 *     Hybrid Automata" (2003), IEEE TAC
 *   Goebel, Sanfelice, Teel, "Hybrid Dynamical Systems" (2009)
 *
 * Course Mapping (L2 Core Concepts):
 *   MIT 6.841    - Execution semantics, hybrid time
 *   Stanford CS359 - Trajectories and reachability
 *   Berkeley EECS 291E - Hybrid trajectories
 *   Cambridge Part III - Hybrid dynamical systems
 *
 * Knowledge points implemented:
 *   KP10: Hybrid time domain definition (HybridTimeSet)
 *   KP11: Execution trajectory definition (HybridExecution)
 *   KP12: Continuous flow segment (HybridFlowSegment)
 *   KP13: Discrete jump point (HybridJumpPoint)
 *   KP14: Determinism check for hybrid automata
 *   KP15: Non-blocking check (no deadlock)
 *   KP16: Zeno execution detection
 *   KP17: Transition enabled check (guard evaluation)
 *   KP18: Invariant satisfaction check
 *   KP19: Parallel composition of hybrid automata
 */

#ifndef HYBRID_EXECUTION_H
#define HYBRID_EXECUTION_H

#include "hybrid_automaton.h"
#include <stddef.h>
#include <stdbool.h>

/* ==========================================================================
 * Hybrid time domain
 * ========================================================================== */

/**
 * @brief A single time interval [t_low, t_high] within one mode.
 *
 * During a flow segment, time advances continuously from t_low to t_high
 * while the system stays in a fixed discrete mode.
 *
 * KP10: Hybrid time is a subset of ℕ × ℝ, structured as a sequence of
 * intervals. This struct is one "step" in that sequence.
 */
typedef struct {
    double t_start;    /**< Start time of this interval */
    double t_end;      /**< End time of this interval */
    int    mode_index; /**< The discrete mode during this interval */
} HybridTimeInterval;

/**
 * @brief A hybrid time set is a finite or infinite sequence of time intervals.
 *
 * KP10 Implementation: H = ∪_{j=0}^{J} [t_j, t_{j+1}] × {j} where J ∈ ℕ ∪ {∞}.
 * Each step j has a continuous time interval and a fixed discrete mode.
 *
 * Properties:
 *   - Intervals are non-overlapping and ordered by time
 *   - t_{j+1} ≥ t_j (monotonic)
 *   - If J is finite, the last interval may be infinite (infinite horizon)
 */
typedef struct {
    HybridTimeInterval *intervals; /**< Array of time intervals */
    int    num_intervals;          /**< Number of intervals in the sequence */
    int    max_intervals;          /**< Allocated capacity */
    bool   is_infinite;            /**< Whether the execution is infinite */
    bool   is_zeno;                /**< Whether Zeno phenomenon detected */
} HybridTimeSet;

/* ==========================================================================
 * Execution state
 * ========================================================================== */

/**
 * @brief State at a point in the execution: (mode, continuous_state, time).
 *
 * KP11: A hybrid state s = (q, x) ∈ Q × ℝⁿ. This is the fundamental
 * state of the hybrid system at any instant.
 */
typedef struct {
    int    mode;                         /**< Current discrete mode q */
    double x[HA_MAX_VARIABLES];          /**< Current continuous state x ∈ ℝⁿ */
    double t;                            /**< Current time */
} HybridState;

/**
 * @brief KP12: A continuous flow segment in one mode.
 *
 * During [t_start, t_end], the system stays in a single mode q,
 * and the continuous state evolves according to ẋ = f_q(x).
 *
 * The segment includes:
 *   - Start state s_start = (q, x_start) at time t_start
 *   - End state s_end = (q, x_end) at time t_end
 *   - ODE solution connecting x_start to x_end
 *
 * For affine flows ẋ = A·x + b, the solution is:
 *   x(t) = e^{A(t-t₀)} x₀ + ∫_{t₀}^{t} e^{A(t-τ)} b dτ
 */
typedef struct {
    int    mode_index;                   /**< Discrete mode during this segment */
    double t_start;                      /**< Start time */
    double t_end;                        /**< End time */
    double x_start[HA_MAX_VARIABLES];    /**< State at t_start */
    double x_end[HA_MAX_VARIABLES];      /**< State at t_end */
    int    num_steps;                    /**< Number of integration steps taken */
} HybridFlowSegment;

/**
 * @brief KP13: A discrete jump point between modes.
 *
 * At time t_jump, the system instantaneously transitions from
 * mode q_src to mode q_tgt. The continuous state is updated
 * according to the reset map: x' = R(e)·x + r(e).
 *
 * The jump is instant (zero duration in continuous time).
 */
typedef struct {
    int    transition_id;                /**< Index of the transition fired */
    int    src_mode;                     /**< Source mode before jump */
    int    tgt_mode;                     /**< Target mode after jump */
    double t_jump;                       /**< Time of the jump */
    double x_pre[HA_MAX_VARIABLES];      /**< State just before jump */
    double x_post[HA_MAX_VARIABLES];     /**< State just after jump (reset applied) */
} HybridJumpPoint;

/**
 * @brief Full execution of a hybrid automaton.
 *
 * KP11 Implementation: An execution χ is a (possibly infinite) sequence
 *   χ = τ₀ → (jump₀) → τ₁ → (jump₁) → τ₂ → ...
 *
 * where each τ_i is a continuous flow segment (HybridFlowSegment) and
 * each jump_i is a discrete transition (HybridJumpPoint).
 *
 * The execution must satisfy:
 *   1. τ₀ starts from an initial state (q₀, x₀) ∈ Init
 *   2. For each flow segment τ_i in mode q_i, ẋ(t) = f_{q_i}(x(t))
 *      and x(t) ∈ Inv(q_i) for all t in [t_i, t_i']
 *   3. For each jump from q_i to q_{i+1} with transition e:
 *      x_pre = x(t_i') satisfies Guard(e), and
 *      x_post = R(e)·x_pre + r(e) ∈ Inv(q_{i+1})
 */
typedef struct {
    HybridAutomaton *automaton;          /**< The automaton being executed */
    HybridFlowSegment *flow_segments;    /**< Array of continuous flow segments */
    HybridJumpPoint  *jump_points;       /**< Array of discrete jumps */
    int    num_segments;                 /**< Number of flow segments */
    int    num_jumps;                    /**< Number of jumps */
    int    max_segments;                 /**< Allocated capacity for segments */
    int    max_jumps;                    /**< Allocated capacity for jumps */
    bool   is_finite;                    /**< Whether execution has terminated */
    bool   is_zeno;                      /**< Whether Zeno detected (infinite jumps in finite time) */
    double total_time;                   /**< Total elapsed hybrid time */
} HybridExecution;

/* ==========================================================================
 * API: Execution construction
 * ========================================================================== */

/**
 * @brief Create a new execution for a given automaton.
 *
 * Allocates an execution with initial capacity for flow segments and jumps.
 * The execution starts at the automaton's initial state.
 *
 * @param ha               The hybrid automaton
 * @param max_flow_segments Initial capacity for flow segments
 * @param max_jumps        Initial capacity for discrete jumps
 * @return                 Newly allocated execution, or NULL on error
 */
HybridExecution* hybrid_execution_create(const HybridAutomaton *ha,
                                          int max_flow_segments, int max_jumps);

/**
 * @brief Destroy an execution and free resources.
 */
void hybrid_execution_destroy(HybridExecution *exec);

/**
 * @brief Append a continuous flow segment to the execution.
 *
 * KP12: Records that the system flowed in a specific mode for
 * a given time interval, reaching x_end from x_start.
 *
 * @param exec       The execution
 * @param mode_index Mode during the flow
 * @param t_start    Start time
 * @param t_end      End time
 * @param x_start    State at t_start (n-dimensional)
 * @param x_end      State at t_end (n-dimensional)
 * @param num_steps  Number of integration steps
 * @return           true on success
 */
bool hybrid_execution_append_flow(HybridExecution *exec, int mode_index,
                                   double t_start, double t_end,
                                   const double *x_start, const double *x_end,
                                   int num_steps);

/**
 * @brief Append a discrete jump to the execution.
 *
 * KP13: Records an instantaneous mode switch at time t_jump.
 *
 * @param exec           The execution
 * @param transition_id  Index of the transition taken
 * @param t_jump         Time of the jump
 * @param x_pre          State before jump (n-dimensional)
 * @param x_post         State after reset (n-dimensional)
 * @return               true on success
 */
bool hybrid_execution_append_jump(HybridExecution *exec, int transition_id,
                                   double t_jump,
                                   const double *x_pre, const double *x_post);

/* ==========================================================================
 * API: Semantic predicates (L2 Core Concepts)
 * ========================================================================== */

/**
 * @brief KP14: Check whether the automaton is deterministic.
 *
 * A hybrid automaton is deterministic at state (q, x) if:
 *   - For each mode q, Flow(q) is a single-valued function (not diff inclusion)
 *   - At any state (q, x), at most one transition guard is satisfied
 *   - Guard regions are mutually exclusive for transitions from the same mode
 *
 * Non-determinism arises from:
 *   - Overlapping guard conditions (multiple transitions enabled simultaneously)
 *   - Differential inclusions (set-valued flows)
 *   - Initial condition sets (rather than points)
 *
 * @param ha The automaton to check
 * @return   true if deterministic
 *
 * Theorem: Determinism is semi-decidable for linear hybrid automata
 *          (Henzinger et al., 1995)
 */
bool hybrid_is_deterministic(const HybridAutomaton *ha);

/**
 * @brief KP15: Check whether the automaton is non-blocking (deadlock-free).
 *
 * A hybrid automaton is non-blocking if from every reachable state (q, x),
 * either:
 *   (a) There exists δ > 0 such that the continuous flow can proceed for
 *       time δ while staying within Inv(q), OR
 *   (b) There exists at least one transition whose guard is satisfied
 *
 * Deadlock occurs when the invariant boundary is reached but no transition
 * guard is enabled — the system cannot flow further and cannot jump.
 *
 * @param ha The automaton to check
 * @return   true if non-blocking for all reachable states
 */
bool hybrid_is_nonblocking(const HybridAutomaton *ha);

/**
 * @brief KP16: Check for Zeno behavior potential.
 *
 * A Zeno execution is one with infinitely many discrete transitions
 * occurring in finite total time. This is a well-known pathological
 * behavior in hybrid systems.
 *
 * Detection conditions (sufficient):
 *   - A cycle of modes exists where time cannot advance (guards trigger
 *     immediately, resets re-enable the same guards)
 *   - Total accumulated flow time → 0 as transitions accumulate
 *
 * @param exec The execution to check
 * @return     true if Zeno behavior is detected
 *
 * Theorem: Zeno detection is undecidable for general hybrid automata
 *          but semi-decidable for specific subclasses (Zhang et al., 2001)
 */
bool hybrid_execution_is_zeno(const HybridExecution *exec);

/**
 * @brief Estimate Zeno accumulation time (pre-Zeno time).
 *
 * For a Zeno execution, returns the finite accumulation time T_∞
 * beyond which no further discrete transitions are possible.
 *
 * @param exec The execution
 * @return     Estimated Zeno time, or INFINITY if not Zeno
 */
double hybrid_execution_zeno_time(const HybridExecution *exec);

/**
 * @brief KP17: Check whether a specific transition is enabled at state (q, x).
 *
 * A transition e = (q_src → q_tgt, G(e), R(e)) is enabled if:
 *   1. The current mode matches q_src
 *   2. The guard G(e) is satisfied by x: A·x ≤ b
 *   3. The target invariant Inv(q_tgt) is satisfied after reset: H·(R·x+r) ≤ k
 *
 * @param ha        The automaton
 * @param trans_id  Transition index
 * @param x         Current continuous state (n-dimensional)
 * @param mode      Current mode (if -1, check regardless of mode)
 * @return          true if the transition can fire
 */
bool hybrid_transition_enabled(const HybridAutomaton *ha, int trans_id,
                                const double *x, int mode);

/**
 * @brief Evaluate a guard condition at a given continuous state.
 *
 * Returns true if A·x ≤ b for all constraints in the guard.
 *
 * @param guard The guard to evaluate
 * @param x     Continuous state (n-dimensional)
 * @param n     Dimension of x
 * @return      true if guard is satisfied
 */
bool hybrid_guard_evaluate(const HybridGuard *guard, const double *x, int n);

/**
 * @brief Apply a reset map to a continuous state.
 *
 * Computes x' = R(e)·x + r(e) and stores result in x_out.
 *
 * @param reset  The reset map
 * @param x      Input state (n-dimensional)
 * @param x_out  Output state (n-dimensional, may alias x)
 * @param n      Dimension
 */
void hybrid_reset_apply(const HybridReset *reset, const double *x,
                         double *x_out, int n);

/**
 * @brief KP18: Check whether a state satisfies a mode's invariant.
 *
 * Returns true if H·x ≤ k for all constraints in the invariant
 * (or if the invariant is unbounded).
 *
 * @param invariant The invariant to check
 * @param x         Continuous state (n-dimensional)
 * @param n         Dimension
 * @return          true if x ∈ Inv(q)
 */
bool hybrid_invariant_satisfied(const HybridInvariant *invariant,
                                 const double *x, int n);

/**
 * @brief Find all transitions enabled from a given state.
 *
 * Fills the provided array with transition IDs that are currently enabled
 * from (mode, x). Returns the number of enabled transitions.
 *
 * @param ha            The automaton
 * @param mode          Current mode
 * @param x             Current state
 * @param enabled_ids   Output array of enabled transition IDs
 * @param max_count     Max entries in enabled_ids
 * @return              Number of enabled transitions (may exceed max_count)
 */
int  hybrid_transitions_enabled(const HybridAutomaton *ha, int mode,
                                 const double *x, int *enabled_ids, int max_count);

/**
 * @brief KP19: Parallel composition of two hybrid automata.
 *
 * H = H₁ || H₂. The parallel composition synchronizes on shared
 * transition labels (event-based synchronization) while allowing
 * independent continuous evolution.
 *
 * The composed automaton has:
 *   Q = Q₁ × Q₂ (product modes)
 *   X = X₁ ∪ X₂ (union of variables)
 *   Sync transitions: (q₁, q₂) → (q₁', q₂') on shared event a
 *   Async transitions: (q₁, q₂) → (q₁', q₂) on H₁ event (H₂ stays)
 *
 * @param h1  First automaton
 * @param h2  Second automaton
 * @return    Composed automaton, or NULL on error
 *
 * Complexity: O(|Q₁|·|Q₂|·|E₁|·|E₂|)
 * Theorem: Bisimulation is a congruence for parallel composition
 *          in hybrid I/O automata (Lynch, Segala, Vaandrager, 2003)
 */
HybridAutomaton* hybrid_compose_parallel(const HybridAutomaton *h1,
                                          const HybridAutomaton *h2);

#endif /* HYBRID_EXECUTION_H */
