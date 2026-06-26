/**
 * @file hss_integrator.h
 * @brief Hybrid Switched System Integration Layer (L3-L4)
 *
 * This module provides the bridge between the unified HSS framework
 * and the individual sub-module representations:
 *
 *   hss_hybrid_automata_bridge  ↔ mini-hybrid-automata
 *   hss_switched_stability_bridge ↔ mini-switched-stability
 *   hss_dwell_time_bridge       ↔ mini-dwell-time-analysis
 *   hss_event_triggered_bridge  ↔ mini-event-triggered-control
 *   hss_impulsive_bridge        ↔ mini-impulsive-system
 *   hss_pwa_bridge              ↔ mini-piecewise-affine-system
 *   hss_reset_bridge            ↔ mini-reset-control-system
 *   hss_supervisory_bridge      ↔ mini-supervisory-control
 *
 * Knowledge Points:
 *   L3 KP1: System morphism — structure-preserving maps between formalisms
 *   L3 KP2: Composition — series/parallel/feedback composition of HSS
 *   L3 KP3: Abstraction hierarchy — lift/lower between levels
 *   L3 KP4: Translation verification — validate mapping correctness
 *   L4 KP1: Equivalence notions — simulation, bisimulation, trace equivalence
 *   L4 KP2: Stability preservation under composition
 *   L4 KP3: Detectability of switching events
 *
 * Course Mapping:
 *   CMU 15-424 — CPS composition & refinement
 *   Stanford CS359 — Hybrid automata composition
 *   Berkeley EECS 291E — System equivalence
 *   MIT 6.241 — System interconnections
 */

#ifndef HSS_INTEGRATOR_H
#define HSS_INTEGRATOR_H

#include "hss_core.h"

/* ============================================================================
 * L3 KP1: System Morphism Types
 * ============================================================================ */

/**
 * @brief Morphism direction between formalisms
 */
typedef enum {
    HSS_MORPH_TO_HYBRID_AUTOMATA    = 0,
    HSS_MORPH_TO_SWITCHED_SYSTEM    = 1,
    HSS_MORPH_TO_DWELL_TIME         = 2,
    HSS_MORPH_TO_EVENT_TRIGGERED    = 3,
    HSS_MORPH_TO_IMPULSIVE          = 4,
    HSS_MORPH_TO_PWA               = 5,
    HSS_MORPH_TO_RESET_CONTROL      = 6,
    HSS_MORPH_TO_SUPERVISORY        = 7,
    HSS_MORPH_FROM_HYBRID_AUTOMATA = 8,
    HSS_MORPH_FROM_SWITCHED_SYSTEM = 9,
    HSS_MORPH_FROM_DWELL_TIME      = 10,
    HSS_MORPH_FROM_EVENT_TRIGGERED = 11,
    HSS_MORPH_FROM_IMPULSIVE       = 12,
    HSS_MORPH_FROM_PWA             = 13,
    HSS_MORPH_FROM_RESET_CONTROL   = 14,
    HSS_MORPH_FROM_SUPERVISORY     = 15
} HSS_MorphismDirection;

/**
 * @brief Equivalence level between system representations (L4 KP1)
 *
 * From the theory of system equivalence (van der Schaft 2004):
 *   - TRACE: same sets of execution traces
 *   - BISIMULATION: mutual step-wise simulation
 *   - SIMULATION: one system simulates another
 *   - APPROXIMATE: systems are ε-close
 *   - BEHAVIORAL: same input-output behavior
 */
typedef enum {
    HSS_EQUIV_TRACE         = 0,
    HSS_EQUIV_BISIMULATION  = 1,
    HSS_EQUIV_SIMULATION    = 2,
    HSS_EQUIV_APPROXIMATE   = 3,
    HSS_EQUIV_BEHAVIORAL    = 4,
    HSS_EQUIV_LYAPUNOV      = 5,
    HSS_EQUIV_NONE          = 6
} HSS_EquivalenceLevel;

/* ============================================================================
 * L3 KP2: Composition Types
 * ============================================================================ */

/**
 * @brief System composition topology (L3 KP2)
 *
 * Standard interconnection topologies from systems theory:
 *   SERIES:    u₁ → S₁ → y₁ → S₂ → y₂
 *   PARALLEL:  u → S₁ || S₂ → y
 *   FEEDBACK:  u → [+]→ S₁ → y → S₂ → [-] (back to +)
 *   CASCADE:   Multi-stage series
 *   MESH:      Arbitrary interconnection via adjacency matrix
 *   HIERARCHICAL: Supervisor-subordinate structure
 */
typedef enum {
    HSS_COMP_SERIES      = 0,
    HSS_COMP_PARALLEL    = 1,
    HSS_COMP_FEEDBACK    = 2,
    HSS_COMP_CASCADE     = 3,
    HSS_COMP_MESH        = 4,
    HSS_COMP_HIERARCHICAL = 5
} HSS_CompositionType;

/**
 * @brief Composed hybrid system result
 *
 * Stores the product automaton of two (or more) hybrid systems.
 * For systems with |Q₁| modes and |Q₂| modes, the product has
 * up to |Q₁| × |Q₂| modes (reduced by synchronization constraints).
 */
typedef struct {
    HSS_System  *product;         /**< Composed HSS                     */
    int         *mode_map;        /**< (q₁,q₂) → product_mode mapping  */
    int          num_components;  /**< Number of composed subsystems    */
    HSS_System **components;      /**< Original subsystems (borrowed)   */
} HSS_Composition;

/* ============================================================================
 * L3 KP3: Abstraction Level Types
 * ============================================================================ */

/**
 * @brief Abstraction level for hybrid systems
 *
 * Different levels of modeling resolution:
 *   L0_SIGNAL:    Input-output signal level
 *   L1_CONTINUOUS: Pure continuous dynamics (ODE)
 *   L2_DISCRETE:   Pure discrete automaton
 *   L3_HYBRID:     Combined continuous + discrete
 *   L4_TIMED:      Timed automaton (continuous = clock)
 *   L5_STOCHASTIC: Stochastic hybrid system
 */
typedef enum {
    HSS_LEVEL_SIGNAL      = 0,
    HSS_LEVEL_CONTINUOUS  = 1,
    HSS_LEVEL_DISCRETE    = 2,
    HSS_LEVEL_HYBRID      = 3,
    HSS_LEVEL_TIMED       = 4,
    HSS_LEVEL_STOCHASTIC  = 5
} HSS_AbstractionLevel;

/* ============================================================================
 * L3-L4 API: Integration Functions
 * ============================================================================ */

/**
 * @brief Compose two hybrid systems (L3 KP2)
 *
 * Computes the synchronous product or interconnection of two HSS.
 * For series composition: output of S₁ feeds input of S₂.
 * For parallel composition: states are stacked, independent evolution.
 * For feedback composition: S₁ output feeds back through S₂.
 *
 * Theory: The synchronous product of hybrid automata is a hybrid
 * automaton. Theorem (Alur/Henzinger 1996): reachability of the
 * product automaton implies reachability of both components.
 *
 * @param comp_type Composition topology
 * @param sys1 First HSS
 * @param sys2 Second HSS
 * @return Composed system (caller must free via hss_composition_free)
 */
HSS_Composition *hss_compose(HSS_CompositionType comp_type,
                              const HSS_System *sys1,
                              const HSS_System *sys2);

/**
 * @brief Compose multiple systems in cascade (L3 KP2)
 * @param systems Array of HSS pointers
 * @param count Number of systems
 * @return Cascade composition
 */
HSS_Composition *hss_compose_cascade(HSS_System **systems, int count);

/**
 * @brief Free a composition result
 * @param comp Composition to free
 */
void hss_composition_free(HSS_Composition *comp);

/**
 * @brief Check equivalence between two system representations (L4 KP1)
 *
 * Verifies whether two hybrid system representations are equivalent
 * at the specified level. Uses the theoretical framework from:
 *   - Milner (1989): Communication and Concurrency (bisimulation)
 *   - van der Schaft (2004): Equivalence of dynamical systems
 *
 * @param sys1 First HSS
 * @param sys2 Second HSS
 * @param level Equivalence level to check
 * @return true if equivalent at the given level
 */
bool hss_are_equivalent(const HSS_System *sys1, const HSS_System *sys2,
                         HSS_EquivalenceLevel level);

/**
 * @brief Verify that a translation preserves stability (L4 KP2)
 *
 * Theorem: If the original system is GAS and the translation
 * preserves Lyapunov function values at switching instants,
 * then the translated system is also GAS.
 *
 * @param original Source HSS
 * @param translated Target HSS
 * @param tolerance Numerical tolerance for LF comparison
 * @return true if stability is preserved
 */
bool hss_verify_stability_preservation(const HSS_System *original,
                                        const HSS_System *translated,
                                        double tolerance);

/**
 * @brief Detect if switching events are observable (L4 KP3)
 *
 * Observability of discrete mode changes from continuous output.
 * Based on the concept of mode distinguishability.
 *
 * @param sys HSS system
 * @param mode_a First mode
 * @param mode_b Second mode
 * @param threshold Distinguishability threshold
 * @return true if modes are distinguishable from output
 */
bool hss_modes_distinguishable(const HSS_System *sys,
                                int mode_a, int mode_b,
                                double threshold);

/**
 * @brief Lift a system to a higher abstraction level (L3 KP3)
 *
 * Abstraction by adding detail: e.g., discrete → hybrid
 * by refining discrete modes with continuous dynamics.
 *
 * @param sys Original system at lower level
 * @param target_level Desired abstraction level
 * @return Lifted system (caller owns), or NULL if impossible
 */
HSS_System *hss_lift_abstraction(const HSS_System *sys,
                                  HSS_AbstractionLevel target_level);

/**
 * @brief Lower (abstract) a system to a coarser level (L3 KP3)
 *
 * Abstraction by removing detail: e.g., hybrid → discrete
 * by quotienting out continuous dynamics (bisimulation quotient).
 *
 * @param sys Original system at higher level
 * @param target_level Desired abstraction level
 * @return Abstracted system (caller owns), or NULL if impossible
 */
HSS_System *hss_lower_abstraction(const HSS_System *sys,
                                   HSS_AbstractionLevel target_level);

/**
 * @brief Merge two hybrid systems sharing state (L3)
 *
 * Creates a unified system where certain modes of sys1 and sys2
 * are identified (merged). Useful for integrating independently
 * modeled subsystems.
 *
 * @param sys1 First system
 * @param sys2 Second system
 * @param mode_identifications Array of (mode_in_sys1, mode_in_sys2) pairs
 * @param num_pairs Number of identification pairs
 * @return Merged system
 */
HSS_System *hss_merge_systems(const HSS_System *sys1,
                               const HSS_System *sys2,
                               const int *mode_identifications,
                               int num_pairs);

#endif /* HSS_INTEGRATOR_H */
