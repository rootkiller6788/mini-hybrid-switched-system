/**
 * @file synthesis.h
 * @brief Supervisor Synthesis Algorithms
 *
 * Computes the supremal controllable sublanguage supC(K, G).
 * Algorithm:
 *   1. Compute the product G || K (synchronous composition)
 *   2. Iteratively remove "bad" states (states where uncontrollable
 *      events lead outside the specification)
 *   3. The remaining states form supC(K, G)
 *   4. If the initial state survives, a supervisor exists.
 *
 * This is essentially the Kumar-Garg algorithm for computing the
 * supremal controllable sublanguage of regular languages.
 *
 * References:
 * - Wonham & Ramadge (1987) "On the supremal controllable sublanguage"
 * - Kumar & Garg (1995) "Modeling and Control of Logical DES"
 */

#ifndef SYNTHESIS_H
#define SYNTHESIS_H

#include "des_automaton.h"
#include "controllability.h"
#include <stddef.h>
#include <stdint.h>

#define SYNTH_MAX_ITER      100
#define SYNTH_MAX_QUEUE    2048

typedef struct {
    int     feasible;
    uint16_t removed_states[DES_MAX_STATES];
    uint16_t n_removed;
    uint16_t iterations;
    uint16_t final_nstates;
} synthesis_result_t;

/** Standard fixpoint algorithm for supC(K, G). */
void synthesis_supremal_controllable(const des_automaton_t *G,
    const des_automaton_t *K, des_automaton_t *supC,
    synthesis_result_t *result);

/** Compute the supremal controllable and nonblocking sublanguage. */
void synthesis_supremal_nonblocking(const des_automaton_t *G,
    const des_automaton_t *K, des_automaton_t *supCNB,
    synthesis_result_t *result);

/** Compute the supremal controllable and normal sublanguage
 *  (for partial observation case). */
void synthesis_supremal_normal(const des_automaton_t *G,
    const des_automaton_t *K, des_automaton_t *supN,
    synthesis_result_t *result);

/** State-based synthesis: return 1 if state q should be in supC. */
int synthesis_is_legal_state(const des_automaton_t *G,
    const des_automaton_t *K, uint16_t qG, uint16_t qK);

/** Remove a state and all transitions to/from it in-place. */
void synthesis_remove_state(des_automaton_t *A, uint16_t q);

/** Compute the set of controllable events that must be disabled
 *  at each state to enforce the specification. */
uint16_t synthesis_disable_map(const des_automaton_t *G,
    const des_automaton_t *K, uint32_t *disable_map, uint16_t cap);

/** Modular synthesis: synthesize a supervisor for each specification
 *  component and then combine. */
int synthesis_modular(const des_automaton_t *G,
    const des_automaton_t **specs, uint16_t nspecs,
    des_automaton_t **supervisors, synthesis_result_t *results);

/** Print synthesis results in human-readable form. */
void synthesis_print_result(const synthesis_result_t *result, FILE *fp);

/** Compute supervisor automaton from supC: the automaton that
 *  implements the supervisory control law. */
void synthesis_build_supervisor_automaton(const des_automaton_t *supC,
    const des_automaton_t *G, des_automaton_t *supervisor);

/** Reduced supervisor: minimize the supervisor automaton. */
void synthesis_reduce_supervisor(des_automaton_t *supervisor);

#endif /* SYNTHESIS_H */
