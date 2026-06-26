/**
 * @file hybrid_safety.h
 * @brief Safety verification for hybrid automata (L4-L5)
 *
 * Safety verification determines whether a hybrid system can reach
 * an unsafe state. Methods include barrier certificates, Lyapunov
 * functions for hybrid systems, and inductive invariants.
 *
 * Reference:
 *   Prajna, Jadbabaie, "Safety Verification of Hybrid Systems Using
 *     Barrier Certificates" (2004), HSCC
 *   Branicky, "Multiple Lyapunov Functions and Other Analysis Tools
 *     for Switched and Hybrid Systems" (1998), IEEE TAC
 *   Platzer, "Logical Analysis of Hybrid Systems" (2010), Springer
 *
 * Course Mapping:
 *   MIT 6.841    - Safety verification, barrier certificates
 *   Stanford CS359 - Hybrid system verification
 *   CMU 15-855   - Inductive invariants
 *   Caltech CS 154 - Safety in hybrid systems
 *
 * Knowledge points:
 *   KP29: Safety property specification
 *   KP30: Barrier certificate definition and checking
 *   KP31: Common Lyapunov function for hybrid systems
 *   KP32: Multiple Lyapunov functions (MLF)
 *   KP33: Inductive invariant generation
 *   KP34: Dwell-time based safety
 *   KP35: Event-triggered safety verification
 */

#ifndef HYBRID_SAFETY_H
#define HYBRID_SAFETY_H

#include "hybrid_automaton.h"
#include "hybrid_execution.h"
#include "hybrid_reachability.h"
#include <stdbool.h>

/* ==========================================================================
 * Safety property specification
 * ========================================================================== */

/**
 * @brief KP29: Safety property for a hybrid automaton.
 *
 * A safety property φ_safe specifies a set S ⊆ Q × ℝⁿ of "safe" states.
 * The system is safe if all reachable states remain within S.
 *
 * Typically, S is defined by its complement — the "unsafe" set U.
 * Safety ≡ □(¬U) in LTL: "always not unsafe."
 */
typedef struct {
    char     name[HA_NAME_LEN];          /**< Property name */
    HybridPolyhedron *unsafe_set;        /**< Unsafe region in continuous space */
    int      unsafe_mode;                /**< Mode where unsafe set is defined, -1 = all modes */
    char     description[256];           /**< Human-readable description */
    bool     is_global;                  /**< true = unsafe applies in all modes */
} HybridSafetyProperty;

/* ==========================================================================
 * Barrier certificates (KP30)
 * ========================================================================== */

/**
 * @brief KP30: Barrier certificate for hybrid systems.
 *
 * A barrier certificate B: Q × ℝⁿ → ℝ is a function satisfying:
 *
 *   (I)   B(q, x) ≤ 0  for all (q, x) ∈ Init     (initial states safe)
 *   (II)  B(q, x) > 0   for all (q, x) ∈ U       (unsafe states excluded)
 *   (III) Ḃ(q, x) ≤ 0   for all (q, x) with B(q, x) = 0 and x ∈ Inv(q)
 *                      (flow cannot cross from safe to unsafe)
 *   (IV)  B(q', R(e)·x + r(e)) ≤ 0  for all (q, x) with B(q, x) ≤ 0
 *                      (transitions preserve safety)
 *
 * If such a B exists, the system is safe — no execution from Init can reach U.
 *
 * This struct represents a specific barrier certificate candidate.
 *
 * Reference: Prajna & Jadbabaie (2004), HSCC
 */
typedef struct {
    int    num_modes;              /**< Number of modes in the automaton */
    int    dim;                    /**< Continuous state dimension */

    /** Quadratic barriers: B_q(x) = x^T P_q x + 2 q_q^T x + r_q */
    /** For each mode q, P_q is a symmetric n×n matrix, q_q is n-vector, r_q scalar */
    double *P;                    /**< P matrices: num_modes × dim × dim (row-major per mode) */
    double *q;                    /**< q vectors: num_modes × dim */
    double *r;                    /**< r scalars: num_modes */

    /** Whether each condition (I-IV) is validated */
    bool   condition_I_holds;
    bool   condition_II_holds;
    bool   condition_III_holds;
    bool   condition_IV_holds;

    /** Verification tolerance */
    double tolerance;
} HybridBarrierCertificate;

/**
 * @brief Allocate a barrier certificate for a hybrid automaton.
 *
 * @param num_modes Number of modes
 * @param dim       State dimension
 * @return          New barrier certificate (zero-initialized)
 */
HybridBarrierCertificate* hybrid_barrier_create(int num_modes, int dim);

/** Destroy a barrier certificate */
void hybrid_barrier_destroy(HybridBarrierCertificate *bc);

/**
 * @brief Check condition I: B(q, x) ≤ 0 for all initial states.
 *
 * @param bc  Barrier certificate
 * @param ha  Automaton (for initial condition)
 * @return    true if condition I holds (numerically, within tolerance)
 */
bool hybrid_barrier_check_init(const HybridBarrierCertificate *bc,
                                const HybridAutomaton *ha);

/**
 * @brief Check condition II: B(q, x) > 0 for all unsafe states.
 *
 * @param bc   Barrier certificate
 * @param sp   Safety property (defines unsafe set)
 * @return     true if condition II holds
 */
bool hybrid_barrier_check_unsafe(const HybridBarrierCertificate *bc,
                                  const HybridSafetyProperty *sp);

/**
 * @brief Check condition III: Ḃ(q, x) ≤ 0 on the zero sublevel set boundary.
 *
 * For affine flows ẋ = A_q x + b_q:
 *   Ḃ_q(x) = 2 x^T P_q (A_q x + b_q) + 2 q_q^T (A_q x + b_q)
 *
 * This must be ≤ 0 whenever B_q(x) = 0 and x ∈ Inv(q).
 *
 * @param bc  Barrier certificate
 * @param ha  Automaton (for flows and invariants)
 * @return    true if condition III holds for all modes
 */
bool hybrid_barrier_check_flow(const HybridBarrierCertificate *bc,
                                const HybridAutomaton *ha);

/**
 * @brief Check condition IV: transitions preserve safety.
 *
 * For each transition e = (q → q', G(e), R(e)):
 *   If B_q(x) ≤ 0 and x satisfies G(e), then B_{q'}(R(e)·x + r(e)) ≤ 0.
 *
 * @param bc  Barrier certificate
 * @param ha  Automaton (for transitions)
 * @return    true if condition IV holds for all transitions
 */
bool hybrid_barrier_check_transitions(const HybridBarrierCertificate *bc,
                                       const HybridAutomaton *ha);

/**
 * @brief Verify the complete barrier certificate (all four conditions).
 *
 * @param bc  Barrier certificate
 * @param ha  Automaton
 * @param sp  Safety property
 * @return    true if all conditions I-IV hold → system is SAFE
 */
bool hybrid_barrier_verify(const HybridBarrierCertificate *bc,
                            const HybridAutomaton *ha,
                            const HybridSafetyProperty *sp);

/* ==========================================================================
 * Lyapunov functions for hybrid safety (KP31-KP32)
 * ========================================================================== */

/**
 * @brief KP31: Common Lyapunov function for all modes.
 *
 * A function V: ℝⁿ → ℝ is a Common Lyapunov Function (CLF) if:
 *   (1) V(x) > 0 for all x ≠ 0, V(0) = 0 (positive definite)
 *   (2) V̇_q(x) = ∇V(x)·f_q(x) < 0 for all modes q (strict decrease)
 *
 * If a CLF exists, the origin is globally asymptotically stable
 * under arbitrary switching. This implies safety for a sufficiently
 * large sublevel set.
 *
 * Quadratic CLF: V(x) = x^T P x with P ≻ 0 and P A_q + A_q^T P ≺ 0 ∀q.
 *
 * @param ha    The automaton (for affine flow matrices)
 * @param n     State dimension
 * @param P_out Output: Lyapunov matrix P (n×n, row-major), caller allocates
 * @return      true if a common Lyapunov function was found
 *
 * Theorem: Existence of CLF ⇔ GAS under arbitrary switching for
 *          switched linear systems (Liberzon & Morse, 1999)
 */
bool hybrid_common_lyapunov(const HybridAutomaton *ha, int n, double *P_out);

/**
 * @brief KP32: Multiple Lyapunov Functions (MLF) for hybrid systems.
 *
 * When a CLF does not exist, MLF methods assign a separate Lyapunov
 * function V_q to each mode q. Stability requires:
 *   (a) V_q decreases during continuous flow in mode q
 *   (b) V_{q'}(x') ≤ V_q(x) for each transition q → q'
 *       (Lyapunov function does not increase on jumps)
 *
 * This is the classic MLF theorem from Branicky (1998).
 *
 * @param ha         The automaton
 * @param n          State dimension
 * @param P_out      Output: P matrices (num_modes × n × n, row-major per mode)
 * @param num_modes  Number of modes
 * @return           true if MLF found for all modes
 */
bool hybrid_multiple_lyapunov(const HybridAutomaton *ha, int n,
                               double *P_out, int num_modes);

/**
 * @brief Check MLF condition (b): V decreases or stays same on transitions.
 *
 * For transition e = (q → q'), checks:
 *   x'^T P_{q'} x' ≤ x^T P_q x   for all x satisfying G(e) and in Inv(q)
 *
 * @param ha         The automaton
 * @param P          Array of P matrices (num_modes × n^2)
 * @param n          State dimension
 * @return           true if all transitions satisfy the MLF condition
 */
bool hybrid_mlf_check_transitions(const HybridAutomaton *ha,
                                   const double *P, int n);

/* ==========================================================================
 * Inductive invariants (KP33)
 * ========================================================================== */

/**
 * @brief KP33: Inductive invariant for hybrid safety.
 *
 * An inductive invariant I ⊆ Q × ℝⁿ satisfies:
 *   (1) Init ⊆ I                     (initial states in I)
 *   (2) I ∩ U = ∅                    (I is disjoint from unsafe)
 *   (3) For any (q, x) ∈ I and any δ > 0:
 *       if flow stays in Inv(q), then (q, ξ_q(x, δ)) ∈ I
 *       (I is preserved by continuous flow)
 *   (4) For any (q, x) ∈ I and transition e = (q → q'):
 *       if x ∈ G(e), then (q', R(e)·x + r(e)) ∈ I
 *       (I is preserved by discrete transitions)
 *
 * This struct represents a candidate inductive invariant as a polyhedral set.
 */
typedef struct {
    int    num_modes;
    HybridPolyhedron **mode_invariants;  /**< I_q for each mode q (NULL = all of ℝⁿ) */
} HybridInductiveInvariant;

/** Allocate an inductive invariant structure */
HybridInductiveInvariant* hybrid_inductive_invariant_create(int num_modes);

/** Destroy an inductive invariant */
void hybrid_inductive_invariant_destroy(HybridInductiveInvariant *inv);

/**
 * @brief Check whether a given set is an inductive invariant.
 *
 * Verifies conditions (1)-(4) for the candidate invariant.
 *
 * @param inv Candidate invariant
 * @param ha  Automaton
 * @param sp  Safety property (for unsafe check)
 * @return    true if the candidate is a valid inductive invariant
 */
bool hybrid_inductive_invariant_check(const HybridInductiveInvariant *inv,
                                       const HybridAutomaton *ha,
                                       const HybridSafetyProperty *sp);

/**
 * @brief KP34: Dwell-time based safety verification.
 *
 * If the automaton has a minimum dwell time τ_d > 0 in each mode
 * (i.e., transitions cannot fire within τ_d of entering a mode),
 * then Lyapunov decrease during flow can compensate for possible
 * increases at transitions.
 *
 * @param ha         The automaton
 * @param dwell_time Minimum dwell time τ_d
 * @param n          State dimension
 * @return           true if dwell time guarantees safety
 *
 * Theorem: Switched linear systems with dwell time τ_d are GAS if
 *          each mode is stable and τ_d is sufficiently large
 *          (Morse, 1996)
 */
bool hybrid_dwell_time_safety(const HybridAutomaton *ha, double dwell_time, int n);

/**
 * @brief KP35: Event-triggered safety check.
 *
 * For automata with controlled transitions, verify that a given
 * event-triggering policy maintains safety. The policy specifies
 * when transition e should be taken based on state x.
 *
 * @param ha            The automaton
 * @param sp            Safety property
 * @param trigger_policy Array of mode→transition mappings (policy[q] = trans_id or -1)
 * @param num_modes     Number of modes
 * @return              true if the policy is safe
 */
bool hybrid_event_triggered_safety(const HybridAutomaton *ha,
                                    const HybridSafetyProperty *sp,
                                    const int *trigger_policy, int num_modes);

/**
 * @brief Create a default safety property that marks a given polyhedron as unsafe.
 *
 * @param name        Property name
 * @param unsafe      The unsafe polyhedron
 * @param unsafe_mode Mode where unsafe is defined (-1 = global)
 * @return            New safety property
 */
HybridSafetyProperty* hybrid_safety_property_create(const char *name,
                                                     const HybridPolyhedron *unsafe,
                                                     int unsafe_mode);

/** Destroy a safety property */
void hybrid_safety_property_destroy(HybridSafetyProperty *sp);

#endif /* HYBRID_SAFETY_H */
