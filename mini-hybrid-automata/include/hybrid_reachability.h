/**
 * @file hybrid_reachability.h
 * @brief Reachability analysis for hybrid automata (L3-L5)
 *
 * Reachability is the fundamental verification problem for hybrid
 * systems: given an initial set X₀ and a target/unsafe set X_unsafe,
 * determine whether any execution starting from X₀ can reach X_unsafe.
 *
 * Reference:
 *   Alur, Dang, Ivančić, "Progress on Reachability Analysis of Hybrid
 *     Systems Using Predicate Abstraction" (2003), HSCC
 *   Frehse et al., "SpaceEx: Scalable Verification of Hybrid Systems"
 *     (2011), CAV
 *   Girard, "Reachability of Uncertain Linear Systems Using Zonotopes"
 *     (2005), HSCC
 *
 * Course Mapping:
 *   MIT 6.841    - Reachability algorithms
 *   Berkeley EECS 291E - Reachability for hybrid systems
 *   CMU 15-855   - Verification algorithms
 *   Oxford        - Hybrid system verification
 *
 * Knowledge points:
 *   KP20: Forward reachability computation
 *   KP21: Backward reachability computation
 *   KP22: Flowpipe computation via time discretization
 *   KP23: Zonotope representation for reachable sets
 *   KP24: Polyhedral reachability (H-representation)
 *   KP25: Support function based reachability
 *   KP26: Bisimulation quotient computation
 *   KP27: Abstraction refinement (CEGAR) loop
 *   KP28: Onion-peeling iterative reachability
 */

#ifndef HYBRID_REACHABILITY_H
#define HYBRID_REACHABILITY_H

#include "hybrid_automaton.h"
#include "hybrid_execution.h"
#include <stdbool.h>

/* ==========================================================================
 * Reachability results and options
 * ========================================================================== */

/**
 * @brief Result of a reachability query.
 */
typedef enum {
    HAREACH_SAFE,         /**< No execution reaches the unsafe set */
    HAREACH_UNSAFE,       /**< Some execution reaches the unsafe set */
    HAREACH_UNKNOWN,      /**< Could not determine (abstraction too coarse) */
    HAREACH_TIMEOUT       /**< Computation exceeded resource bounds */
} HybridReachabilityResult;

/**
 * @brief Reachability computation options.
 */
typedef struct {
    double time_horizon;      /**< Maximum time to explore (0 = unlimited) */
    int    max_iterations;    /**< Maximum discrete transition steps */
    double time_step;         /**< Time discretization step for flowpipe */
    bool   use_zonotope;      /**< Use zonotope representation */
    bool   use_support_func;  /**< Use support function method */
    int    verbosity;         /**< Verbosity level (0=silent, 1=progress, 2=detailed) */
} HybridReachOptions;

/**
 * @brief Default reachability options.
 */
#define HYBRID_REACH_OPTIONS_DEFAULT { \
    .time_horizon = 0.0, \
    .max_iterations = 1000, \
    .time_step = 0.01, \
    .use_zonotope = true, \
    .use_support_func = false, \
    .verbosity = 0 \
}

/* ==========================================================================
 * Set representations for reachability
 * ========================================================================== */

/**
 * @brief KP23: Zonotope — center + linear combination of generators.
 *
 * Z = { c + Σ_{i=1}^{p} β_i g_i | β_i ∈ [-1, 1] }
 *
 * Zonotopes are closed under Minkowski sum and linear transformation,
 * making them efficient for reachability of affine dynamics.
 * They are centrally symmetric convex polytopes.
 *
 * Reference: Girard (2005), Kuhn (1998)
 */
typedef struct {
    int    dim;              /**< Dimension of the state space */
    int    num_generators;   /**< Number of generators p */
    int    max_generators;   /**< Allocated capacity */
    double *center;          /**< Center point c ∈ ℝⁿ (length dim) */
    double *generators;      /**< Generator matrix G ∈ ℝ^{dim × p}, column-major */
} HybridZonotope;

/**
 * @brief KP24: Convex polyhedron in H-representation.
 *
 * P = { x ∈ ℝⁿ | H·x ≤ k }
 *
 * The H-representation uses half-spaces defined by normal vectors
 * and bound values. This is the most general convex polyhedral
 * representation for reachability.
 */
typedef struct {
    int    dim;              /**< Dimension of state space */
    int    num_halfspaces;   /**< Number of half-space constraints */
    int    max_halfspaces;   /**< Allocated capacity */
    double *H;               /**< Constraint normals (max_hs × dim, row-major) */
    double *k;               /**< Constraint bounds (length max_hs) */
} HybridPolyhedron;

/**
 * @brief Hyper-rectangle (axis-aligned box).
 *
 * B = { x ∈ ℝⁿ | lo_i ≤ x_i ≤ hi_i, i = 1..n }
 *
 * Simplest set representation — efficient for rectangular hybrid automata.
 */
typedef struct {
    int    dim;
    double *lo;
    double *hi;
} HybridHyperRect;

/* ==========================================================================
 * Set operations for reachability
 * ========================================================================== */

/** Create a zonotope with given dimension and generator capacity */
HybridZonotope* hybrid_zono_create(int dim, int max_generators);
void hybrid_zono_destroy(HybridZonotope *z);

/** Set the center of a zonotope */
void hybrid_zono_set_center(HybridZonotope *z, const double *c);
/** Add a generator vector to a zonotope */
bool hybrid_zono_add_generator(HybridZonotope *z, const double *g);
/** Apply affine transformation: Z' = A·Z + b (Minkowski operations) */
bool hybrid_zono_affine_transform(HybridZonotope *z, const double *A,
                                   const double *b, int A_rows, int A_cols);
/** Over-approximate Minkowski sum of two zonotopes */
HybridZonotope* hybrid_zono_minkowski_sum(const HybridZonotope *z1,
                                           const HybridZonotope *z2);
/** Check if zonotope intersects a polyhedron (support function based) */
bool hybrid_zono_intersects_poly(const HybridZonotope *z,
                                  const HybridPolyhedron *p);
/** Enclose zonotope in a hyper-rectangle (interval hull) */
void hybrid_zono_to_hyperrect(const HybridZonotope *z, HybridHyperRect *r);

/** Create polyhedron with given dimension */
HybridPolyhedron* hybrid_poly_create(int dim, int max_halfspaces);
void hybrid_poly_destroy(HybridPolyhedron *p);
/** Add a half-space constraint H_row·x ≤ k_val */
bool hybrid_poly_add_halfspace(HybridPolyhedron *p, const double *H_row, double k_val);
/** Check if two polyhedra intersect */
bool hybrid_poly_intersects(const HybridPolyhedron *p1, const HybridPolyhedron *p2);
/** Check if point x satisfies all constraints (x ∈ P) */
bool hybrid_poly_contains(const HybridPolyhedron *p, const double *x);

/** Create a hyper-rectangle */
HybridHyperRect* hybrid_hyperrect_create(int dim);
void hybrid_hyperrect_destroy(HybridHyperRect *r);
/** Check if a point is inside the hyper-rectangle */
bool hybrid_hyperrect_contains(const HybridHyperRect *r, const double *x);

/* ==========================================================================
 * Forward reachability (KP20)
 * ========================================================================== */

/**
 * @brief KP20: Forward reachability analysis for a hybrid automaton.
 *
 * Starting from the initial set Init, iteratively computes the set of
 * states reachable through continuous flow and discrete transitions.
 *
 * Algorithm (overview):
 *   1. R₀ ← Init  (initial reachable set)
 *   2. For each iteration k:
 *      a. Continuous post:  Post_C(R_k) = states reachable by flowing
 *         from R_k within each mode's invariant
 *      b. Discrete post:    Post_D(R_k) = states reachable by taking
 *         one transition from R_k
 *      c. R_{k+1} ← R_k ∪ Post_C(R_k) ∪ Post_D(R_k)
 *   3. Fixed point:  R_∞ = lim_{k→∞} R_k
 *
 * This implements the classical reachability algorithm from
 * Alur, Courcoubetis, Henzinger, Ho (1995).
 *
 * @param ha       The hybrid automaton
 * @param unsafe   The unsafe set (polyhedron), NULL if just computing reach set
 * @param options  Computation options
 * @param result   Output: reachability result
 * @return         The final reachable set (as a polyhedron), or NULL on error
 *
 * Complexity: PSPACE-hard in general (Henzinger et al., 1995)
 *             EXPTIME for timed automata (Alur & Dill, 1994)
 */
HybridPolyhedron* hybrid_reachable_forward(const HybridAutomaton *ha,
                                            const HybridPolyhedron *unsafe,
                                            const HybridReachOptions *options,
                                            HybridReachabilityResult *result);

/**
 * @brief KP21: Backward reachability analysis.
 *
 * Starting from the unsafe/target set, computes the set of states
 * that can reach it. More efficient than forward when the target
 * set is small relative to the state space.
 *
 * Algorithm: Pre-image computation through flows and transitions.
 *
 * @param ha       The hybrid automaton
 * @param target   The target set to reach
 * @param options  Computation options
 * @param result   Output: reachability result
 * @return         The backward reachable set (polyhedron), or NULL
 */
HybridPolyhedron* hybrid_reachable_backward(const HybridAutomaton *ha,
                                             const HybridPolyhedron *target,
                                             const HybridReachOptions *options,
                                             HybridReachabilityResult *result);

/**
 * @brief KP22: Flowpipe computation for a single mode.
 *
 * Given an initial set X₀ ⊆ ℝⁿ and an affine flow ẋ = A·x + b,
 * compute the set of states reachable within time [0, T] while
 * staying within invariant Inv.
 *
 * The flowpipe is the union over t ∈ [0,T] of the reachable states
 * at time t:  F_{[0,T]} = ⋃_{t∈[0,T]} e^{At} X₀ + ∫₀ᵗ e^{A(t-τ)} b dτ
 *
 * @param ha        The automaton (for variable dimension)
 * @param mode_id   The mode whose flow is used
 * @param X0        Initial set (zonotope)
 * @param T         Time horizon
 * @param time_step Discretization step
 * @return          The flowpipe as an array of zonotopes (one per time step)
 *                  Caller receives count via *num_steps out parameter.
 *
 * Complexity: O(T/Δt · (dim³ + dim²·p)) for affine flows with p generators
 */
HybridZonotope** hybrid_flowpipe_compute(const HybridAutomaton *ha, int mode_id,
                                          const HybridZonotope *X0, double T,
                                          double time_step, int *num_steps);

/**
 * @brief KP25: Support function based overapproximation.
 *
 * For a convex set S, the support function ρ_S(d) = max_{x∈S} d^T x.
 * This provides a tight overapproximation using template directions.
 *
 * Used for efficient reachability with large linear systems.
 *
 * @param z       The zonotope to bound
 * @param dirs    Array of direction vectors
 * @param num_dirs Number of directions
 * @param bounds  Output: ρ_S(d_i) for each direction
 */
void hybrid_support_function(const HybridZonotope *z, const double **dirs,
                              int num_dirs, double *bounds);

/**
 * @brief KP28: Onion-peeling iterative reachability.
 *
 * Alternates between continuous post and discrete post computation,
 * "peeling layers" of the reachable set outward from the initial region.
 * Each layer adds the states reachable in one more discrete step.
 *
 * @param ha           The automaton
 * @param max_discrete_steps Maximum number of discrete transition layers
 * @param options      Reachability options
 * @param layers       Output: number of layers computed
 * @return             The layered reachable sets (array of polyhedra, size *layers)
 */
HybridPolyhedron** hybrid_onion_peeling(const HybridAutomaton *ha,
                                         int max_discrete_steps,
                                         const HybridReachOptions *options,
                                         int *layers);

/* ==========================================================================
 * Bisimulation and abstraction (KP26-KP27)
 * ========================================================================== */

/**
 * @brief KP26: Compute the bisimulation quotient of a hybrid automaton.
 *
 * A bisimulation relation ≡ is an equivalence relation on states such that
 * if s₁ ≡ s₂, then for every transition/flow from s₁, there exists a
 * corresponding transition/flow from s₂ leading to equivalent states.
 *
 * The quotient automaton H/≡ has equivalence classes as its modes.
 * Bisimulation preserves all CTL* properties (and thus safety/liveness).
 *
 * @param ha     The automaton
 * @param max_eq_classes Maximum number of equivalence classes (for bounded comp.)
 * @return       The quotient automaton, or NULL if not computable
 *
 * Theorem: For timed automata, region equivalence is a finite bisimulation
 *          (Alur & Dill, 1994). For general HA, finite bisimulation may not exist.
 */
HybridAutomaton* hybrid_bisimulation_quotient(const HybridAutomaton *ha,
                                               int max_eq_classes);

/**
 * @brief KP27: Counter-Example Guided Abstraction Refinement (CEGAR).
 *
 * CEGAR iteratively refines an abstract model based on spurious
 * counterexamples found during model checking:
 *
 *   1. Build initial abstraction (coarse partition of state space)
 *   2. Model-check the abstraction
 *   3. If SAFE → original system is SAFE (abstraction is conservative)
 *   4. If counterexample found → check if concretizable:
 *      a. If real → return UNSAFE with witness execution
 *      b. If spurious → refine abstraction to eliminate it, goto 2
 *
 * @param ha            The automaton
 * @param unsafe        The unsafe set
 * @param options       Reachability options
 * @param max_refinements Maximum CEGAR iterations
 * @param result        Output: SAFE, UNSAFE, or UNKNOWN
 * @param witness       Output: witness execution (if UNSAFE), NULL otherwise
 * @return              true if the algorithm converged (SAFE or UNSAFE)
 *
 * Reference: Clarke, Grumberg, Jha, Lu, Veith, "Counterexample-Guided
 *            Abstraction Refinement for Symbolic Model Checking" (2003), JACM
 */
bool hybrid_cegar(const HybridAutomaton *ha, const HybridPolyhedron *unsafe,
                  const HybridReachOptions *options, int max_refinements,
                  HybridReachabilityResult *result, HybridExecution **witness);

#endif /* HYBRID_REACHABILITY_H */
