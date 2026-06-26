/**
 * @file hybrid_subclass.h
 * @brief Subclass specializations: timed, rectangular, initialized automata (L3-L4)
 *
 * Certain subclasses of hybrid automata have important decidability
 * properties and efficient algorithms that the general class lacks.
 *
 * Reference:
 *   Alur & Dill, "A Theory of Timed Automata" (1994), TCS
 *   Henzinger, Kopke, Puri, Varaiya, "What's Decidable About Hybrid
 *     Automata?" (1995), STOC
 *   Lafferriere, Pappas, Sastry, "O-Minimal Hybrid Systems" (2000)
 *
 * Course Mapping:
 *   MIT 6.841    - Decidability frontiers
 *   Stanford CS359 - Subclass analysis
 *   Berkeley EECS 291E - Timed and rectangular HA
 *   ETH 263-4650 - Decidability in hybrid systems
 *
 * Knowledge points:
 *   KP41: Timed automaton representation and regions
 *   KP42: Region equivalence for timed automata
 *   KP43: Rectangular hybrid automaton classification
 *   KP44: Initialized rectangular automata decidability
 *   KP45: Linear hybrid automaton (LHA) representation
 *   KP46: Piecewise affine system (PWA) representation
 *   KP47: Mixed Logical Dynamical (MLD) system representation
 *   KP48: O-minimal hybrid automata
 */

#ifndef HYBRID_SUBCLASS_H
#define HYBRID_SUBCLASS_H

#include "hybrid_automaton.h"
#include "hybrid_reachability.h"
#include <stdbool.h>

/* ==========================================================================
 * KP41: Timed Automaton
 * ========================================================================== */

/**
 * @brief KP41: A timed automaton is a special hybrid automaton where:
 *
 *   - All continuous variables are CLOCKS (rate = 1, start at 0)
 *   - Flows: ẋ_i = 1 for all clocks x_i
 *   - Invariants: conjunctions of x_i ≤ c or x_i < c (clock upper bounds)
 *   - Guards: conjunctions of x_i ≥ c or x_i > c (clock lower bounds)
 *   - Resets: x_i' = 0 (clock reset) or x_i' = x_i (identity)
 *
 * Timed automata have DECIDABLE reachability (Alur & Dill, 1994).
 * The decision procedure uses region equivalence.
 */

/** Maximum value a clock is compared to (defines region graph size) */
#define TA_MAX_CONSTANT 100

/**
 * @brief Clock constraint type.
 */
typedef enum {
    TACC_LT,    /**< x < c */
    TACC_LE,    /**< x ≤ c */
    TACC_GT,    /**< x > c */
    TACC_GE,    /**< x ≥ c */
    TACC_EQ     /**< x = c */
} TimedConstraintType;

/**
 * @brief A single clock constraint.
 */
typedef struct {
    int    clock_index;           /**< Which clock this constrains */
    TimedConstraintType type;    /**< Comparison operator */
    int    value;                 /**< Integer constant c ∈ ℕ₀ */
} TimedConstraint;

/**
 * @brief Region equivalence class for timed automata.
 *
 * Two clock valuations ν, ν' are region-equivalent (ν ≈ ν') iff:
 *   1. For each clock x, either ⌊ν(x)⌋ = ⌊ν(x')⌋ or both > c_max
 *   2. For clocks with ν(x) ≤ c_max, frac(ν(x)) = 0 ⇔ frac(ν(x')) = 0
 *   3. For clocks with ν(x) ≤ c_max and ν(y) ≤ c_max:
 *      frac(ν(x)) ≤ frac(ν(y)) ⇔ frac(ν(x')) ≤ frac(ν(y'))
 *
 * The number of regions is finite: O(|C|! · 2^{|C|} · c_max^{|C|})
 */
typedef struct {
    int    num_clocks;
    int    max_constant;
    int    *int_parts;            /**< Integer parts ⌊x_i⌋ for each clock */
    int    *frac_order;           /**< Fractional ordering permutation */
    int    *zero_frac;            /**< Which clocks have zero fractional part */
} TimedRegion;

/**
 * @brief Convert a hybrid automaton to a timed automaton if possible.
 *
 * Checks whether all variables are clocks, all flows are ẋ = 1,
 * all invariants/guards/resets satisfy the timed automaton constraints.
 *
 * @param ha            The hybrid automaton to check/convert
 * @return              true if the automaton is a valid timed automaton
 */
bool hybrid_is_timed_automaton(const HybridAutomaton *ha);

/**
 * @brief Get the maximum constant appearing in constraints for a clock.
 *
 * @param ha          The automaton
 * @param clock_index Clock variable index
 * @return            Maximum constant c such that x_i ∼ c appears
 */
int  hybrid_clock_max_constant(const HybridAutomaton *ha, int clock_index);

/**
 * @brief KP42: Compute the region of a clock valuation.
 *
 * Given a clock valuation ν, returns its region equivalence class.
 *
 * @param valuation    Clock values (num_clocks-dimensional)
 * @param num_clocks   Number of clocks
 * @param max_constant Maximum constant
 * @param region       Output: filled region struct
 *
 * Theorem: The region graph of a timed automaton is finite and
 *          preserves reachability properties (Alur & Dill, 1994).
 */
void hybrid_region_compute(const double *valuation, int num_clocks,
                            int max_constant, TimedRegion *region);

/**
 * @brief Check if two clock valuations are region-equivalent.
 *
 * @param v1, v2       Clock valuations
 * @param num_clocks   Number of clocks
 * @param max_constant Maximum constant
 * @return             true if v1 ≈ v2 (same region)
 */
bool hybrid_region_equivalent(const double *v1, const double *v2,
                               int num_clocks, int max_constant);

/**
 * @brief Compute the successor region of a given region.
 *
 * Advances time until a clock hits its next integer boundary
 * or until all clocks exceed max_constant.
 *
 * @param region     Current region (updated in-place to successor)
 * @param num_clocks Number of clocks
 * @param max_constant Maximum constant
 * @return           true if successor exists (finite)
 */
bool hybrid_region_successor(TimedRegion *region, int num_clocks, int max_constant);

/**
 * @brief Build the region graph for a timed automaton.
 *
 * The region graph nodes are (mode, region) pairs. An edge exists
 * if a successor can be reached by time elapse or a discrete transition.
 *
 * @param ha              The automaton (must be timed)
 * @param num_regions_out Output: number of region graph nodes
 * @param max_constant    Maximum clock constant
 * @return                Array of region pointers (length *num_regions_out)
 */
TimedRegion** hybrid_region_graph_build(const HybridAutomaton *ha,
                                          int max_constant,
                                          int *num_regions_out);

/* ==========================================================================
 * KP43-KP44: Rectangular Hybrid Automata
 * ========================================================================== */

/**
 * @brief KP43: A rectangular hybrid automaton (RHA) has:
 *
 *   - Flow inclusions of the form: ẋ_i ∈ [a_{q,i}, b_{q,i}]
 *     (constant rate bounds that depend only on mode q)
 *   - Invariants: x_i ∈ [lo_{q,i}, hi_{q,i}] (axis-aligned boxes)
 *   - Guards: x_i ∈ [g_lo_{e,i}, g_hi_{e,i}] (axis-aligned intervals)
 *   - Resets: x_i' ∈ [r_lo_{e,i}, r_hi_{e,i}] or x_i' = x_i
 *
 * This struct captures the interval bounds for rectangular dynamics.
 */
typedef struct {
    double *rate_lo;     /**< Lower rate bounds: num_modes × dim */
    double *rate_hi;     /**< Upper rate bounds: num_modes × dim */
    double *inv_lo;      /**< Invariant lower bounds: num_modes × dim */
    double *inv_hi;      /**< Invariant upper bounds: num_modes × dim */
    double *guard_lo;    /**< Guard lower bounds: num_trans × dim */
    double *guard_hi;    /**< Guard upper bounds: num_trans × dim */
    double *reset_lo;    /**< Reset lower bounds: num_trans × dim */
    double *reset_hi;    /**< Reset upper bounds: num_trans × dim */
    int    num_modes;
    int    num_trans;
    int    dim;
} RectangularHAData;

/**
 * @brief Check if a hybrid automaton is rectangular.
 *
 * @param ha The automaton
 * @return   true if all dynamics, invariants, guards, resets are rectangular
 */
bool hybrid_is_rectangular(const HybridAutomaton *ha);

/**
 * @brief Extract rectangular HA data from a hybrid automaton.
 *
 * @param ha   The automaton (must be rectangular)
 * @param data Output: filled rectangular HA data (caller allocates)
 * @return     true if extraction successful
 */
bool hybrid_rectangular_extract(const HybridAutomaton *ha, RectangularHAData *data);

/**
 * @brief KP44: Check if a rectangular HA is initialized.
 *
 * An initialized rectangular HA requires that after any discrete transition
 * that resets a variable x_i, the continuous dynamics of x_i in the target
 * mode are uniquely determined (i.e., rate_lo = rate_hi for that variable).
 *
 * Initialized RHA have DECIDABLE reachability (Henzinger et al., 1995).
 *
 * @param ha The automaton
 * @return   true if initialized
 */
bool hybrid_is_initialized_rectangular(const HybridAutomaton *ha);

/* ==========================================================================
 * KP45-KP47: LHA, PWA, MLD
 * ========================================================================== */

/**
 * @brief KP45: Linear Hybrid Automaton (LHA).
 *
 * LHA constraints:
 *   - Flows: ẋ = c_q (constant derivative, possibly different per mode)
 *   - Invariants, Guards: linear inequalities H·x ≤ k
 *   - Resets: affine x' = R·x + r
 *
 * LHA reachability is undecidable in general, but decidable for
 * various subclasses (e.g., timed, initialized rectangular).
 */
bool hybrid_is_linear_hybrid(const HybridAutomaton *ha);

/**
 * @brief KP46: Check if automaton represents a Piecewise Affine (PWA) system.
 *
 * A PWA system has:
 *   - State space partitioned into polyhedral regions
 *   - In each region, dynamics are affine: ẋ = A_i x + b_i
 *   - No reset maps (continuous state is preserved across region boundaries)
 *
 * PWA systems are a common modeling framework for hybrid MPC.
 *
 * @param ha The automaton
 * @return   true if the automaton is equivalent to a PWA system
 *
 * Reference: Sontag (1981), "Nonlinear Regulation: The Piecewise Linear Approach"
 */
bool hybrid_is_pwa(const HybridAutomaton *ha);

/**
 * @brief KP47: Convert a hybrid automaton to Mixed Logical Dynamical (MLD) form.
 *
 * MLD systems (Bemporad & Morari, 1999) model hybrid dynamics as:
 *   x(t+1) = A x(t) + B₁ u(t) + B₂ δ(t) + B₃ z(t)
 *   y(t)   = C x(t) + D₁ u(t) + D₂ δ(t) + D₃ z(t)
 *   E₂ δ(t) + E₃ z(t) ≤ E₁ u(t) + E₄ x(t) + E₅
 *
 * where δ ∈ {0,1} are binary auxiliary variables and z are continuous
 * auxiliary variables.
 *
 * MLD is used for hybrid MPC with MIQP/MILP solvers.
 *
 * @param ha             The automaton
 * @param Ts             Sampling time for discretization
 * @param A_out, B_out   Output matrices (caller allocates)
 * @param n              State dimension
 * @return               true if conversion successful
 */
bool hybrid_to_mld(const HybridAutomaton *ha, double Ts,
                    double *A_out, double *B_out, int n);

/* ==========================================================================
 * KP48: O-minimal hybrid automata
 * ========================================================================== */

/**
 * @brief KP48: Check if the automaton's definable sets are o-minimal.
 *
 * O-minimality is a property of the first-order theory of the reals
 * that guarantees finiteness of certain topological operations.
 *
 * O-minimal hybrid automata have decidable reachability for a broad
 * class of dynamics (Lafferriere, Pappas, Sastry, 2000).
 *
 * The key requirement is that all guards, invariants, resets, and
 * flows are definable in an o-minimal expansion of the real field
 * (e.g., ℝ_{exp}, ℝ_{an}).
 *
 * @param ha The automaton
 * @return   true if the automaton is in an o-minimal fragment
 */
bool hybrid_is_ominimal(const HybridAutomaton *ha);

/**
 * @brief Estimate the number of connected components of the reachable set.
 *
 * For o-minimal systems, this number is finite. This function
 * provides a coarse upper bound based on the formula complexity.
 *
 * @param ha The automaton
 * @return   Upper bound on number of connected components, or -1 if unknown
 */
int  hybrid_reachable_components_bound(const HybridAutomaton *ha);

/* ==========================================================================
 * Utility: subclass classification
 * ========================================================================== */

/**
 * @brief Classify a hybrid automaton into its most specific subclass.
 */
typedef enum {
    HACLASS_GENERAL,            /**< General hybrid automaton */
    HACLASS_LINEAR_HYBRID,      /**< Linear Hybrid Automaton (LHA) */
    HACLASS_RECTANGULAR,        /**< Rectangular Hybrid Automaton (RHA) */
    HACLASS_INIT_RECT,          /**< Initialized Rectangular HA */
    HACLASS_TIMED,              /**< Timed Automaton (TA) */
    HACLASS_PWA,                /**< Piecewise Affine System */
    HACLASS_FINITE_STATE        /**< Finite State Machine (no continuous dynamics) */
} HybridAutomatonClass;

/**
 * @brief Classify the given automaton.
 *
 * @param ha The automaton
 * @return   The most specific applicable subclass
 */
HybridAutomatonClass hybrid_classify(const HybridAutomaton *ha);

#endif /* HYBRID_SUBCLASS_H */
