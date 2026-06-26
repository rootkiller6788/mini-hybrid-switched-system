/**
 * @file hss_analysis.h
 * @brief Hybrid Switched System Analysis Framework (L4-L5)
 *
 * This module provides analytical methods for hybrid switched systems:
 *
 * L4 — Standards / Theorems (with code verification):
 *   KP1: Common Lyapunov Function Theorem — CLF ⇒ GUAS for arbitrary switching
 *   KP2: Multiple Lyapunov Function Theorem — MLF + decrease at switches ⇒ GAS
 *   KP3: Dwell-Time Stability Theorem — τ_d > ln(μ)/λ ⇒ GAS
 *   KP4: Average Dwell-Time Theorem — τ_a > ln(μ)/λ₀ ⇒ GAS
 *   KP5: LaSalle Invariance Principle for hybrid systems
 *   KP6: Matrosov's Theorem for time-varying/hybrid systems
 *   KP7: Small-Gain Theorem for interconnected hybrid systems
 *
 * L5 — Algorithms / Methods:
 *   KP1: CLF search via convex optimization (LMI)
 *   KP2: MLF computation with edge constraints
 *   KP3: Minimum dwell-time computation
 *   KP4: Switching signal synthesis for stabilization
 *   KP5: Barrier certificate generation for safety
 *   KP6: Reachable set approximation via flowpipe
 *   KP7: Zeno detection via inter-event time analysis
 *   KP8: Bisimulation quotient computation
 *
 * Course Mapping:
 *   MIT 6.241 — Lyapunov stability theory
 *   Stanford CS359 — Hybrid system verification
 *   Berkeley EECS 291E — Safety & reachability
 *   CMU 15-424 — CPS verification
 *   ETH 227-0690 — Stability of hybrid systems
 *   Cambridge Part II — Nonlinear & hybrid control
 *   Georgia Tech CS 7641 — ML + control
 */

#ifndef HSS_ANALYSIS_H
#define HSS_ANALYSIS_H

#include "hss_core.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * L4 KP1-KP2: Lyapunov Function Structures
 * ============================================================================ */

/**
 * @brief Common Lyapunov Function (CLF)
 *
 * V(x) = xᵀ P x where P = Pᵀ > 0 and P is common across ALL modes.
 * Theorem (Liberzon 2003): If ∃ P > 0 such that A_qᵀ P + P A_q < 0
 * for all q ∈ Q, then the switched system is GUAS for arbitrary switching.
 */
typedef struct {
    double *P;              /**< Common Lyapunov matrix (n×n)            */
    int     n;              /**< State dimension                         */
    double  max_eigenvalue; /**< λ_max(P)                                */
    double  min_eigenvalue; /**< λ_min(P)                                */
    double  condition_num;  /**< κ(P) = λ_max / λ_min                    */
    bool    is_valid;       /**< True if P satisfies CLF conditions      */
    double  margin;         /**< Stability margin: -max eig(A_qᵀP+PA_q) */
} HSS_CommonLyapunov;

/**
 * @brief Multiple Lyapunov Functions (MLF)
 *
 * For each mode q, V_q(x) = xᵀ P_q x with P_q > 0.
 * Theorem (Branicky 1998): If V_q decreases while mode q is active
 * AND V_{q'}(x(t_k)) ≤ V_q(x(t_k)) at switching instants t_k,
 * then the system is GAS.
 */
typedef struct {
    double **P;             /**< Mode-dependent Lyapunov matrices P_q    */
    double  *decay_rates;   /**< Decay rates α_q: V̇_q ≤ -α_q V_q       */
    int      num_modes;     /**< Number of modes                         */
    int      n;             /**< State dimension                         */
    bool    *mode_valid;     /**< LMI feasibility per mode               */
    bool     overall_valid;  /**< True if all modes valid + switch cond. */
    double  *mu;             /**< Switching ratio: V_i ≤ μ V_j (if >1)   */
    double   max_mu;         /**< Maximum switching ratio                 */
} HSS_MultipleLyapunov;

/* ============================================================================
 * L4 KP3-KP4: Dwell-Time Stability Structures
 * ============================================================================ */

/**
 * @brief Dwell-time stability result (L4 KP3)
 *
 * For a switched system with stable subsystems (all A_q Hurwitz),
 * let λ_min = min_q (-α(A_q)) be the minimum decay rate.
 * If the dwell time τ_d > ln(μ) / λ_min, then GAS holds,
 * where μ = max_{i,j} λ_max(P_i)/λ_min(P_j) is the "overshoot" ratio.
 */
typedef struct {
    double   tau_d_min;          /**< Minimum dwell time for stability   */
    double   tau_d_computed;     /**< Computed dwell time bound          */
    double   lambda_min;         /**< Min decay rate across all modes    */
    double   mu_overshoot;       /**< Overshoot ratio μ                  */
    double  *mode_decay_rates;   /**< Decay rate for each mode           */
    int      num_modes;
    bool     is_stable;          /**< True if τ_d > τ_d_min              */
    bool     all_modes_stable;   /**< True if all subsystems are stable  */
    HSS_StabilityConcept verdict;
} HSS_DwellTimeResult;

/**
 * @brief Average dwell-time stability result (L4 KP4)
 *
 * If τ_a > ln(μ) / λ₀ for some λ₀ > 0, and the chatter bound N₀ ≥ 1,
 * then the switched system is GUES.
 * This allows occasional fast switching as long as the average
 * interval between switches is long enough.
 */
typedef struct {
    double   tau_a_min;          /**< Minimum average dwell time         */
    double   tau_a_computed;     /**< Actual average dwell time          */
    double   lambda_0;           /**< Decay parameter λ₀                 */
    double   mu_overshoot;       /**< Overshoot ratio μ                  */
    double   N0;                 /**< Chatter bound N₀                   */
    int      num_switches;       /**< Number of switches analyzed        */
    double   time_interval;      /**< Analysis interval length           */
    bool     is_stable;          /**< True if τ_a satisfies condition    */
    HSS_StabilityConcept verdict;
} HSS_AverageDwellResult;

/* ============================================================================
 * L4 KP5-KP6: Invariance Principles
 * ============================================================================ */

/**
 * @brief LaSalle Invariance for hybrid systems (L4 KP5)
 *
 * Extended LaSalle's invariance principle to hybrid systems.
 * Theorem (Sanfelice/Goebel/Teel 2007): If V is non-increasing along
 * flows and non-increasing at jumps, then trajectories converge to
 * the largest weakly invariant set in {x | V̇(x) = 0}.
 */
typedef struct {
    double  *omega_limit_set;    /**< Approximated ω-limit set           */
    int      set_size;           /**< Size of the limit set              */
    double   convergence_rate;   /**< Rate of convergence to limit set   */
    int      iterations;         /**< Iterations to convergence           */
    bool     is_convergent;      /**< Whether convergence was detected   */
} HSS_LaSalleResult;

/**
 * @brief Matrosov's theorem result (L4 KP6)
 *
 * For non-autonomous/hybrid systems where Lyapunov function is not
 * strictly decreasing (V̇ ≤ 0), Matrosov's theorem provides an
 * auxiliary function W to establish asymptotic stability.
 */
typedef struct {
    double  *W;              /**< Auxiliary function evaluated at points */
    bool     is_W_bounded;   /**< W is bounded along trajectories       */
    bool     is_asymptotic;  /**< Asymptotic stability established      */
    double   V_scalar;       /**< Lyapunov function value               */
    double   W_upper_bound;  /**< Upper bound on |W|                   */
    int      check_points;   /**< Number of verification points         */
} HSS_MatrosovResult;

/* ============================================================================
 * L4 KP7: Small-Gain Theorem
 * ============================================================================ */

/**
 * @brief Small-gain interconnection analysis (L4 KP7)
 *
 * For two interconnected hybrid subsystems S₁ and S₂ with ISS gains
 * γ₁, γ₂: if γ₁ · γ₂ < 1 (in the nonlinear case, if the gain
 * composition is a contraction), then the interconnection is ISS.
 */
typedef struct {
    double   gain_12;           /**< Gain from S₂ output to S₁ state     */
    double   gain_21;           /**< Gain from S₁ output to S₂ state     */
    double   composite_gain;    /**< Composite gain γ₁₂ · γ₂₁           */
    bool     small_gain_holds;  /**< True if composite gain < 1          */
    bool     is_interconnection_ISS; /**< Overall ISS stability          */
    double   iss_margin;        /**< 1 - composite_gain (positive⇒ISS)  */
} HSS_SmallGainResult;

/* ============================================================================
 * L5 KP1-KP4: Algorithm Results
 * ============================================================================ */

/**
 * @brief CLF search result (L5 KP1)
 *
 * Solves the LMI feasibility problem:
 *   Find P = Pᵀ > 0 s.t. A_qᵀ P + P A_q < 0 for all q ∈ Q.
 * Uses iterative eigenvalue-based method (simpler than full SDP solver).
 */
typedef struct {
    HSS_CommonLyapunov clf;     /**< Found CLF (if any)                 */
    int    iterations;          /**< Iterations used in search          */
    double residual;            /**< Final residual norm                */
    bool   found;               /**< True if CLF was found              */
    char   method[64];          /**< Method used for the search          */
} HSS_CLFSearchResult;

/**
 * @brief MLF computation result (L5 KP2)
 *
 * Computes mode-dependent Lyapunov functions P_q satisfying:
 *   A_qᵀ P_q + P_q A_q < -α_q P_q  (decay in mode q)
 *   xᵀ P_{q'} x ≤ μ xᵀ P_q x     (at switch q → q')
 */
typedef struct {
    HSS_MultipleLyapunov mlf;   /**< Computed MLF                       */
    int    iterations;           /**< Total iterations                  */
    int    modes_converged;     /**< Modes that converged               */
    double max_residual;        /**< Maximum residual across modes      */
    bool   solved;              /**< True if MLF was found              */
} HSS_MLFComputationResult;

/**
 * @brief Minimum dwell-time computation (L5 KP3)
 *
 * Algorithm: For each mode pair (q₁, q₂), compute the minimum time
 * that must elapse before switching from q₁ to q₂ maintains stability.
 * Uses eigenvalue analysis of A_q and Lyapunov overshoot ratios.
 */
typedef struct {
    HSS_DwellTimeResult dwell_result;     /**< Dwell-time analysis      */
    double  *pairwise_dwell;              /**< τ_d(qi,qj) matrix        */
    int      num_pairs;                   /**< Number of mode pairs     */
    double   global_min_dwell;            /**< min_{i,j} τ_d(i,j)      */
    bool     mode_pair_stable[128][128];  /**< Per-pair stability flag  */
} HSS_DwellComputation;

/**
 * @brief Switching signal synthesis (L5 KP4)
 *
 * Given a switched system with unstable modes, synthesizes a
 * switching signal σ(t) that renders the system stable.
 * Uses: state-dependent switching (min rule), hysteresis switching,
 * or dwell-time switching.
 *
 * Theory (Liberzon 2003, Chapter 3): If there exists a convex
 * combination of unstable subsystems that is stable, a stabilizing
 * switching signal can be constructed.
 */
typedef struct {
    HSS_SwitchingSignal signal;       /**< Synthesized switching signal */
    bool    is_stabilizing;           /**< True if signal stabilizes    */
    double  decay_rate_achieved;      /**< Achieved decay rate          */
    char    method[64];               /**< Synthesis method used         */
    int     method_enum;              /**< 0=min-rule, 1=hysteresis,
                                            2=dwell-time, 3=MPC-like    */
} HSS_SwitchingSynthesis;

/* ============================================================================
 * L5 KP5-KP8: Safety & Reachability
 * ============================================================================ */

/**
 * @brief Barrier certificate (L5 KP5)
 *
 * A barrier certificate B(x) satisfies:
 *   B(x) ≤ 0 for all x ∈ Init (initial states)
 *   B(x) > 0 for all x ∈ Unsafe (unsafe states)
 *   Ḃ(x) ≤ 0 along flows (non-increasing)
 *   B(x⁺) ≤ B(x⁻) at jumps
 * If such B exists, the unsafe region is unreachable from Init.
 */
typedef struct {
    double *coeffs;           /**< Polynomial coefficients of B(x)       */
    int     degree;           /**< Degree of barrier polynomial          */
    int     n;                /**< State dimension                       */
    bool    is_valid;         /**< True if barrier conditions hold       */
    double  init_max;         /**< max_{Init} B(x)                       */
    double  unsafe_min;       /**< min_{Unsafe} B(x)                     */
    double  margin;           /**< Separation margin: unsafe_min-init_max*/
} HSS_BarrierCertificate;

/**
 * @brief Reachable set (L5 KP6)
 *
 * Approximates the set of states reachable from Init within time T.
 * Uses flowpipe construction: sequence of polyhedra over-approximating
 * the reachable set over time intervals.
 */
typedef struct {
    double  *vertices;           /**< Polyhedron vertices for flowpipe   */
    int      num_vertices;       /**< Total vertices across flowpipe     */
    int      num_segments;       /**< Number of flowpipe segments       */
    double  *segment_times;      /**< Time interval per segment          */
    double   time_horizon;       /**< Total reachability horizon T       */
    bool     safe;               /**< True if no unsafe state reachable  */
    bool     over_approximate;   /**< True if result is over-approximation*/
} HSS_ReachableSet;

/**
 * @brief Zeno detection analysis (L5 KP7)
 *
 * Detects potential Zeno behavior: infinite discrete transitions
 * in finite time. Uses inter-event time analysis and Lyapunov-based
 * sufficient conditions to rule out Zeno.
 */
typedef struct {
    bool     zeno_detected;         /**< Zeno behavior found             */
    double   zeno_time;             /**< Accumulation time (if zeno)     */
    int      jump_count;            /**< Jumps analyzed                  */
    double   min_inter_event;       /**< Minimum inter-event time        */
    double   avg_inter_event;       /**< Average inter-event time        */
    bool     exclusion_possible;    /**< Zeno exclusion can be proved    */
    char     exclusion_method[128]; /**< Method used for exclusion      */
} HSS_ZenoAnalysis;

/**
 * @brief Bisimulation quotient (L5 KP8)
 *
 * Computes the coarsest bisimulation equivalence on the hybrid
 * state space, yielding a finite quotient automaton that is
 * bisimilar to the original infinite-state system.
 *
 * Theory (Alur/Henzinger/Lafferriere/Pappas 2000):
 * For o-minimal hybrid systems, the bisimulation quotient
 * exists and is computable.
 */
typedef struct {
    HSS_System *quotient;        /**< Finite quotient automaton          */
    int        *state_partition;  /**< Partition map: state → block      */
    int         num_blocks;      /**< Number of equivalence classes      */
    bool        is_finite;       /**< True if quotient is finite         */
    bool        bisimilar;       /**< True if bisimulation verified      */
} HSS_BisimulationQuotient;

/* ============================================================================
 * L4-L5 API: Analysis Functions
 * ============================================================================ */

/* ---- L4 Theorem Verifications ---- */

/**
 * @brief Verify Common Lyapunov Function theorem (L4 KP1)
 *
 * Checks: Does there exist P > 0 such that A_qᵀP + PA_q < 0 for all q?
 *
 * @param sys HSS system (must be linear)
 * @param tolerance Numerical tolerance
 * @return CLF result with validity flag
 */
HSS_CommonLyapunov hss_verify_clf_theorem(const HSS_System *sys,
                                           double tolerance);

/**
 * @brief Verify Multiple Lyapunov Function theorem (L4 KP2)
 *
 * Checks Branicky's MLF conditions for stability under switching.
 *
 * @param sys HSS system
 * @param P_array Array of candidate P_q matrices (or NULL to search)
 * @param mu_expected Expected overshoot ratio
 * @return MLF result with validity flags
 */
HSS_MultipleLyapunov hss_verify_mlf_theorem(const HSS_System *sys,
                                              double **P_array,
                                              double mu_expected);

/**
 * @brief Compute dwell-time stability (L4 KP3)
 *
 * Computes the minimum dwell time τ_d needed for GAS.
 *
 * @param sys HSS system
 * @return Dwell-time stability result
 */
HSS_DwellTimeResult hss_compute_dwell_time(const HSS_System *sys);

/**
 * @brief Compute average dwell-time stability (L4 KP4)
 *
 * Verifies if the system is GAS under average dwell-time constraints.
 *
 * @param sys HSS system
 * @param signal Switching signal to analyze
 * @return Average dwell-time result
 */
HSS_AverageDwellResult hss_compute_average_dwell(
    const HSS_System *sys, const HSS_SwitchingSignal *signal);

/**
 * @brief Apply LaSalle invariance principle (L4 KP5)
 *
 * Computes the largest invariant set for hybrid systems.
 *
 * @param sys HSS system
 * @param trace Execution trace to analyze
 * @return LaSalle invariance result
 */
HSS_LaSalleResult hss_lasalle_invariance(const HSS_System *sys,
                                          const HSS_ExecutionTrace *trace);

/**
 * @brief Verify Matrosov's theorem (L4 KP6)
 *
 * Checks asymptotic stability using auxiliary function W.
 *
 * @param sys HSS system
 * @param V Lyapunov function values along trajectory
 * @param num_points Number of trajectory points
 * @return Matrosov result
 */
HSS_MatrosovResult hss_verify_matrosov(const HSS_System *sys,
                                         const double *V, int num_points);

/**
 * @brief Small-gain theorem for interconnected hybrid systems (L4 KP7)
 *
 * @param sys1 First subsystem
 * @param sys2 Second subsystem
 * @param gain_12 Gain bound from S₂→S₁
 * @param gain_21 Gain bound from S₁→S₂
 * @return Small-gain analysis result
 */
HSS_SmallGainResult hss_small_gain_analysis(
    const HSS_System *sys1, const HSS_System *sys2,
    double gain_12, double gain_21);

/* ---- L5 Algorithms ---- */

/**
 * @brief Search for a Common Lyapunov Function (L5 KP1)
 *
 * Uses gradient-based iterative method to find P > 0 satisfying
 * LMI constraints for all modes.
 *
 * Complexity: O(num_modes · n³ · max_iter)
 *
 * @param sys HSS system
 * @param max_iterations Maximum iterations
 * @param step_size Gradient step size
 * @return CLF search result
 */
HSS_CLFSearchResult hss_search_clf(const HSS_System *sys,
                                     int max_iterations, double step_size);

/**
 * @brief Compute Multiple Lyapunov Functions (L5 KP2)
 *
 * Sequential computation of P_q for each mode with switching constraints.
 *
 * @param sys HSS system
 * @param mu_desired Desired overshoot ratio
 * @return MLF computation result
 */
HSS_MLFComputationResult hss_compute_mlf(const HSS_System *sys,
                                           double mu_desired);

/**
 * @brief Compute minimum dwell time for stability (L5 KP3)
 *
 * Algorithm based on eigenvalue analysis and overshoot computation.
 *
 * @param sys HSS system
 * @return Dwell-time computation result
 */
HSS_DwellComputation hss_compute_min_dwell(const HSS_System *sys);

/**
 * @brief Synthesize stabilizing switching signal (L5 KP4)
 *
 * Constructs σ(t) that stabilizes the switched system.
 *
 * @param sys HSS system
 * @param method 0=min-rule, 1=hysteresis, 2=dwell-time
 * @param time_horizon Time horizon for synthesis
 * @return Switching synthesis result
 */
HSS_SwitchingSynthesis hss_synthesize_switching(
    const HSS_System *sys, int method, double time_horizon);

/**
 * @brief Generate barrier certificate for safety (L5 KP5)
 *
 * Attempts to find a polynomial barrier certificate B(x) that
 * separates initial states from unsafe states.
 *
 * @param sys HSS system
 * @param unsafe_set Vertices of unsafe polytope
 * @param num_vertices Number of vertices in unsafe polytope
 * @param degree Desired degree of barrier polynomial
 * @return Barrier certificate result
 */
HSS_BarrierCertificate hss_generate_barrier_certificate(
    const HSS_System *sys,
    const double *unsafe_set, int num_vertices, int degree);

/**
 * @brief Compute over-approximate reachable set (L5 KP6)
 *
 * Flowpipe construction for hybrid reachability analysis.
 *
 * @param sys HSS system
 * @param config Simulation configuration
 * @return Reachable set (over-approximation)
 */
HSS_ReachableSet hss_compute_reachable_set(
    const HSS_System *sys, const HSS_SimConfig *config);

/**
 * @brief Detect and characterize Zeno behavior (L5 KP7)
 *
 * Analyzes an execution trace for Zeno accumulation points.
 *
 * @param trace Execution trace
 * @param zeno_threshold Inter-event time threshold for Zeno suspicion
 * @return Zeno analysis result
 */
HSS_ZenoAnalysis hss_detect_zeno(const HSS_ExecutionTrace *trace,
                                   double zeno_threshold);

/**
 * @brief Compute bisimulation quotient (L5 KP8)
 *
 * Attempts to construct a finite bisimulation quotient.
 *
 * @param sys HSS system
 * @param grid_size Discretization grid size for state space
 * @return Bisimulation quotient result
 */
HSS_BisimulationQuotient hss_compute_bisimulation(
    const HSS_System *sys, double grid_size);

/* ---- Auxiliary Functions ---- */

/**
 * @brief Check if a matrix is Hurwitz (all eigenvalues have negative real part)
 * @param A Square matrix (row-major)
 * @param n Matrix dimension
 * @return true if Hurwitz
 */
bool is_hurwitz(const double *A, int n);

/**
 * @brief Approximate spectral abscissa (maximum real part of eigenvalues)
 * @param A Square matrix (row-major)
 * @param n Matrix dimension
 * @return Approximate spectral abscissa
 */
double spectral_abscissa(const double *A, int n);

#endif /* HSS_ANALYSIS_H */
