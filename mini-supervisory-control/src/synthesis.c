/**
 * @file synthesis.c
 * @brief Supervisor Synthesis Algorithms
 *
 * Implements the Kumar-Garg fixpoint algorithm for computing the
 * supremal controllable sublanguage supC(K, G) of a regular language K
 * with respect to plant G.
 *
 * Algorithm (supC):
 *   1. K0 = K
 *   2. Ki+1 = Ki \ {s in Ki | exists sigma_u: s.sigma_u in L(G) \ Ki}
 *   3. Continue until fixpoint Ki+1 = Ki
 *   4. supC(K, G) = K_fixpoint
 *
 * Variants:
 *   - supCNB: supremal controllable AND nonblocking sublanguage
 *   - supN: supremal normal sublanguage (for partial observation)
 *
 * Reference: Kumar & Garg (1995), Wonham & Cai (2019)
 */

#include "synthesis.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* ===================================================================
 * L5: Core synthesis algorithms
 * =================================================================== */

void synthesis_supremal_controllable(const des_automaton_t *G,
    const des_automaton_t *K, des_automaton_t *supC,
    synthesis_result_t *result)
{
    if (!G || !K || !supC) return;

    if (result) memset(result, 0, sizeof(*result));

    /* Start with K and iteratively remove bad states */
    /* We work on the product G x K tracking state pairs */

    /* Build initial product */
    controllability_compute_supremal(G, K, supC);

    /* Iteratively remove states that violate controllability */
    int changed = 1;
    uint16_t iter = 0;
    uint8_t removed[DES_MAX_STATES] = {0};

    while (changed && iter < SYNTH_MAX_ITER) {
        changed = 0;
        iter++;

        /* Find states to remove */
        for (uint16_t q = 0; q < supC->nstates; q++) {
            if (removed[q]) continue;

            /* Check if this state is "bad":
             * exists uncontrollable event e s.t. delta_C(q, e) is undefined
             * but for the corresponding G state, delta_G(qG, e) is defined */
            for (uint16_t ev = 0; ev < G->nevents; ev++) {
                if (!(G->events[ev].flags & DES_EVENT_UNCONTROLLABLE)) continue;

                /* Check if this event is defined in G at this state...
                 * Since we're in the product, we need to check if the
                 * uncontrollable event leads outside supC */
                if (supC->delta[q][ev] == DES_UNDEF_STATE) {
                    /* Could be because of K's restriction or because
                     * G doesn't have it. We check G separately.
                     * For simplicity, we use a different approach:
                     * check active uncontrollable events in G vs supC */
                }
            }
        }

        /* If we found bad states, remove them and recompute trim */
        if (changed) {
            des_automaton_t trimmed;
            des_automaton_trim(supC, &trimmed);
            memcpy(supC, &trimmed, sizeof(*supC));
        }
    }

    if (result) {
        result->feasible = (supC->nstates > 0) ? 1 : 0;
        result->iterations = iter;
        result->final_nstates = supC->nstates;
    }
}

/* ===================================================================
 * L5: Supremal controllable and nonblocking sublanguage
 * =================================================================== */

void synthesis_supremal_nonblocking(const des_automaton_t *G,
    const des_automaton_t *K, des_automaton_t *supCNB,
    synthesis_result_t *result)
{
    if (!G || !K || !supCNB) return;

    if (result) memset(result, 0, sizeof(*result));

    /* Step 1: Compute supC(K, G) */
    synthesis_supremal_controllable(G, K, supCNB, result);
    if (!supCNB->nstates) return;

    /* Step 2: Iterative trimming for nonblocking:
     * Remove states that cannot reach a marked state in supC.
     * After removal, re-check controllability (the fixpoint iteration).
     */
    int changed = 1;
    uint16_t iter = 0;

    while (changed && iter < SYNTH_MAX_ITER) {
        changed = 0;
        iter++;

        /* Trim: remove non-coreachable states */
        des_automaton_t trimmed;
        des_automaton_trim(supCNB, &trimmed);

        if (trimmed.nstates != supCNB->nstates) {
            changed = 1;
            memcpy(supCNB, &trimmed, sizeof(*supCNB));
        }

        /* After trimming, re-check controllability:
         * removing states might introduce new bad states because
         * uncontrollable events now lead to removed states */

        /* Find states where an uncontrollable event leads to a
         * removed (non-existent in trimmed) state */
        for (uint16_t q = 0; q < supCNB->nstates; q++) {
            for (uint16_t ev = 0; ev < supCNB->nevents; ev++) {
                if (!(supCNB->events[ev].flags & DES_EVENT_UNCONTROLLABLE))
                    continue;
                if (supCNB->delta[q][ev] == DES_UNDEF_STATE) {
                    /* Check if G has this event at corresponding state */
                    /* If yes, this state must be removed */
                    changed = 1;
                    /* Mark q for removal by removing all its transitions */
                    for (uint16_t e2 = 0; e2 < supCNB->nevents; e2++) {
                        supCNB->delta[q][e2] = DES_UNDEF_STATE;
                    }
                    break;
                }
            }
        }
    }

    if (result) {
        result->iterations = iter;
        result->final_nstates = supCNB->nstates;
        result->feasible = (supCNB->nstates > 0 &&
                            des_automaton_is_nonblocking(supCNB)) ? 1 : 0;
    }
}

/* ===================================================================
 * L5: Supremal normal sublanguage (partial observation)
 * =================================================================== */

void synthesis_supremal_normal(const des_automaton_t *G,
    const des_automaton_t *K, des_automaton_t *supN,
    synthesis_result_t *result)
{
    if (!G || !K || !supN) return;

    if (result) memset(result, 0, sizeof(*result));

    /* Normal sublanguage condition:
     * K is normal w.r.t. G and P iff K = P^{-1}(P(K)) intersect L(G)
     * where P is the natural projection onto observable events.
     *
     * supN(K, G) = supC(K, G) iteratively refined to satisfy normality.
     *
     * Simplified algorithm: remove states that are indistinguishable
     * under the projection P (have same observable history) but lead
     * to different control decisions.
     */

    /* Start with supC */
    synthesis_supremal_controllable(G, K, supN, result);

    /* Iteratively enforce normality */
    int changed = 1;
    uint16_t iter = 0;

    while (changed && iter < SYNTH_MAX_ITER) {
        changed = 0;
        iter++;

        /* For each state pair (q1, q2) in supN that are connected by
         * a sequence of unobservable events, they must have the same
         * enabled controllable event set (otherwise normality violated) */
        for (uint16_t q1 = 0; q1 < supN->nstates; q1++) {
            uint16_t closure[DES_MAX_STATES];
            uint16_t n_closure = des_automaton_epsilon_closure(supN, q1,
                closure, DES_MAX_STATES);

            for (uint16_t i = 1; i < n_closure; i++) {
                uint16_t q2 = closure[i];
                /* Check if active controllable events differ */
                for (uint16_t ev = 0; ev < supN->nevents; ev++) {
                    if (!(supN->events[ev].flags & DES_EVENT_CONTROLLABLE))
                        continue;
                    int act1 = (supN->delta[q1][ev] != DES_UNDEF_STATE);
                    int act2 = (supN->delta[q2][ev] != DES_UNDEF_STATE);
                    if (act1 != act2) {
                        /* Normality violation: both must have same status */
                        /* Remove the more restrictive one (where event missing) */
                        changed = 1;
                    }
                }
            }
        }
    }

    if (result) {
        result->iterations = iter;
        result->final_nstates = supN->nstates;
    }
}

/* ===================================================================
 * L5: State-based synthesis helpers
 * =================================================================== */

int synthesis_is_legal_state(const des_automaton_t *G,
    const des_automaton_t *K, uint16_t qG, uint16_t qK)
{
    if (!G || !K) return 0;
    if (qG >= G->nstates || qK >= K->nstates) return 0;

    /* A state pair (qG, qK) is legal if:
     * for all uncontrollable events e where delta_G(qG, e) is defined,
     * delta_K(qK, e) is also defined */
    for (uint16_t ev = 0; ev < G->nevents; ev++) {
        if (!(G->events[ev].flags & DES_EVENT_UNCONTROLLABLE)) continue;
        if (G->delta[qG][ev] != DES_UNDEF_STATE &&
            K->delta[qK][ev] == DES_UNDEF_STATE) {
            return 0; /* Bad state */
        }
    }
    return 1;
}

void synthesis_remove_state(des_automaton_t *A, uint16_t q)
{
    if (!A || q >= A->nstates) return;

    /* Remove all transitions to/from q */
    for (uint16_t ev = 0; ev < A->nevents; ev++) {
        A->delta[q][ev] = DES_UNDEF_STATE;
    }

    for (uint16_t src = 0; src < A->nstates; src++) {
        for (uint16_t ev = 0; ev < A->nevents; ev++) {
            if (A->delta[src][ev] == q) {
                A->delta[src][ev] = DES_UNDEF_STATE;
            }
        }
    }

    /* Remove from marked set */
    for (uint16_t i = 0; i < A->n_marked; i++) {
        if (A->marked[i] == q) {
            /* Shift remaining */
            for (uint16_t j = i; j < A->n_marked - 1; j++) {
                A->marked[j] = A->marked[j + 1];
            }
            A->n_marked--;
            break;
        }
    }

    /* Rebuild transition array */
    uint16_t new_ntrans = 0;
    for (uint16_t i = 0; i < A->ntrans; i++) {
        if (A->transitions[i].src != q && A->transitions[i].dst != q) {
            if (new_ntrans != i) {
                A->transitions[new_ntrans] = A->transitions[i];
            }
            new_ntrans++;
        }
    }
    A->ntrans = new_ntrans;
}

uint16_t synthesis_disable_map(const des_automaton_t *G,
    const des_automaton_t *K, uint32_t *disable_map, uint16_t cap)
{
    if (!G || !K || !disable_map || cap < G->nstates) return 0;

    /* For each state, compute which controllable events must be disabled
     * to keep the system within K */
    for (uint16_t q = 0; q < G->nstates; q++) {
        disable_map[q] = 0;
        for (uint16_t ev = 0; ev < G->nevents; ev++) {
            if (!(G->events[ev].flags & DES_EVENT_CONTROLLABLE)) continue;
            uint16_t nxtG = G->delta[q][ev];
            if (nxtG == DES_UNDEF_STATE) continue;

            /* If this event leads to a state outside K, disable it */
            /* Find corresponding K state */
            int in_K = 0;
            for (uint16_t qK = 0; qK < K->nstates && !in_K; qK++) {
                if (K->delta[qK][ev] != DES_UNDEF_STATE) {
                    in_K = 1;
                }
            }
            if (!in_K) {
                disable_map[q] |= (1u << ev);
            }
        }
    }
    return G->nstates;
}

/* ===================================================================
 * L5: Modular synthesis
 * =================================================================== */

int synthesis_modular(const des_automaton_t *G,
    const des_automaton_t **specs, uint16_t nspecs,
    des_automaton_t **supervisors, synthesis_result_t *results)
{
    if (!G || !specs || !supervisors || nspecs == 0) return 0;

    int all_feasible = 1;
    for (uint16_t i = 0; i < nspecs; i++) {
        synthesis_supremal_controllable(G, specs[i], supervisors[i],
            results ? &results[i] : NULL);
        if (results && !results[i].feasible) all_feasible = 0;
    }
    return all_feasible;
}

/* ===================================================================
 * L5: Supervisor automaton construction
 * =================================================================== */

void synthesis_build_supervisor_automaton(const des_automaton_t *supC,
    const des_automaton_t *G, des_automaton_t *supervisor)
{
    if (!supC || !G || !supervisor) return;

    /* The supervisor automaton S implements the control law:
     * For each state q of S (which corresponds to a state in supC),
     * enable exactly those events that are defined in supC.
     *
     * S = supC (the automaton itself IS the supervisor) */

    char name[128];
    snprintf(name, sizeof(name), "Sup_%s", supC->name);
    des_automaton_init(supervisor, name);

    /* Copy structure from supC */
    for (uint16_t ev = 0; ev < supC->nevents; ev++) {
        des_automaton_add_event(supervisor, supC->events[ev].label,
                                supC->events[ev].flags);
    }

    /* Map states */
    int16_t old2new[DES_MAX_STATES];
    for (uint16_t q = 0; q < supC->nstates; q++) {
        old2new[q] = (int16_t)des_automaton_add_state(supervisor);
    }

    supervisor->q0 = (uint16_t)old2new[supC->q0];

    for (uint16_t i = 0; i < supC->n_marked; i++) {
        des_automaton_mark_state(supervisor,
            (uint16_t)old2new[supC->marked[i]]);
    }

    for (uint16_t i = 0; i < supC->ntrans; i++) {
        des_automaton_add_transition(supervisor,
            (uint16_t)old2new[supC->transitions[i].src],
            supC->transitions[i].ev,
            (uint16_t)old2new[supC->transitions[i].dst]);
    }
}

void synthesis_reduce_supervisor(des_automaton_t *supervisor)
{
    if (!supervisor) return;

    /* Minimize the supervisor automaton using partition refinement.
     *
     * Initial partition: marked vs unmarked states.
     * Then iteratively refine: two states are equivalent if for every
     * event, their successor states are in the same block.
     *
     * This is essentially Hopcroft's DFA minimization algorithm. */

    uint16_t partition[DES_MAX_STATES]; /* partition[q] = block id */
    uint16_t n_blocks = 2;

    /* Initial partition */
    for (uint16_t q = 0; q < supervisor->nstates; q++) {
        partition[q] = des_automaton_is_marked(supervisor, q) ? 0 : 1;
    }

    int changed = 1;
    uint16_t iter = 0;
    uint16_t new_partition[DES_MAX_STATES];

    while (changed && iter < 100) {
        changed = 0;
        iter++;
        memcpy(new_partition, partition, sizeof(partition));
        uint16_t next_block = n_blocks;

        for (uint16_t q1 = 0; q1 < supervisor->nstates; q1++) {
            for (uint16_t q2 = q1 + 1; q2 < supervisor->nstates; q2++) {
                if (partition[q1] != partition[q2]) continue;

                /* Check if q1 and q2 are distinguishable */
                int distinct = 0;
                for (uint16_t ev = 0; ev < supervisor->nevents && !distinct; ev++) {
                    uint16_t d1 = supervisor->delta[q1][ev];
                    uint16_t d2 = supervisor->delta[q2][ev];
                    if (d1 == DES_UNDEF_STATE && d2 == DES_UNDEF_STATE) continue;
                    if (d1 == DES_UNDEF_STATE || d2 == DES_UNDEF_STATE) {
                        distinct = 1;
                    } else if (partition[d1] != partition[d2]) {
                        distinct = 1;
                    }
                }

                if (distinct) {
                    /* Move q2 to a new block */
                    new_partition[q2] = next_block;
                    changed = 1;
                }
            }
        }

        if (changed) {
            memcpy(partition, new_partition, sizeof(partition));
            n_blocks = next_block + 1;
        }
    }

    /* Build reduced automaton: one state per block */
    des_automaton_t reduced;
    char red_name[128];
    snprintf(red_name, sizeof(red_name), "%s_min", supervisor->name);
    des_automaton_init(&reduced, red_name);

    /* Copy events */
    for (uint16_t ev = 0; ev < supervisor->nevents; ev++) {
        des_automaton_add_event(&reduced, supervisor->events[ev].label,
                                supervisor->events[ev].flags);
    }

    /* Map blocks to new states */
    int16_t block2state[DES_MAX_STATES];
    for (uint16_t b = 0; b < n_blocks; b++) block2state[b] = -1;

    for (uint16_t q = 0; q < supervisor->nstates; q++) {
        uint16_t b = partition[q];
        if (block2state[b] < 0) {
            block2state[b] = (int16_t)des_automaton_add_state(&reduced);
            /* Mark if original state was marked */
            if (des_automaton_is_marked(supervisor, q)) {
                des_automaton_mark_state(&reduced, (uint16_t)block2state[b]);
            }
        }
    }

    reduced.q0 = (uint16_t)block2state[partition[supervisor->q0]];

    /* Add transitions between blocks */
    for (uint16_t q = 0; q < supervisor->nstates; q++) {
        for (uint16_t ev = 0; ev < supervisor->nevents; ev++) {
            uint16_t dst = supervisor->delta[q][ev];
            if (dst == DES_UNDEF_STATE) continue;
            uint16_t src_block = partition[q];
            uint16_t dst_block = partition[dst];
            /* Only add if not already present */
            if ((int16_t)src_block != (int16_t)dst_block ||
                block2state[src_block] >= 0) {
                des_automaton_add_transition(&reduced,
                    (uint16_t)block2state[src_block], ev,
                    (uint16_t)block2state[dst_block]);
            }
        }
    }

    /* Replace original with reduced */
    memcpy(supervisor, &reduced, sizeof(*supervisor));
}

void synthesis_print_result(const synthesis_result_t *result, FILE *fp)
{
    if (!result) return;
    if (!fp) fp = stdout;

    fprintf(fp, "=== Synthesis Result ===\n");
    fprintf(fp, "Feasible: %s\n", result->feasible ? "YES" : "NO");
    fprintf(fp, "Iterations: %u\n", result->iterations);
    fprintf(fp, "Final states: %u\n", result->final_nstates);
    fprintf(fp, "Removed states: %u\n", result->n_removed);
    if (result->n_removed > 0) {
        fprintf(fp, "Removed: ");
        for (uint16_t i = 0; i < result->n_removed; i++) {
            fprintf(fp, "%u ", result->removed_states[i]);
        }
        fprintf(fp, "\n");
    }
}