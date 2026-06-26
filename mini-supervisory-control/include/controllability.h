/**
 * @file controllability.h
 * @brief Controllability Analysis
 *
 * Central theorem: K is controllable wrt G iff
 *   forall s in prefix(K), forall sigma in Sigma_u:
 *   s.sigma in L(G) => s.sigma in prefix(K)
 *
 * References:
 * - Ramadge & Wonham (1987) The control of discrete event systems
 * - Cassandras & Lafortune (2008) Intro to Discrete Event Systems
 */

#ifndef CONTROLLABILITY_H
#define CONTROLLABILITY_H

#include "des_automaton.h"
#include <stddef.h>
#include <stdint.h>

#define CTRL_MAX_CEX_LEN   256

typedef struct {
    uint16_t counterexample[CTRL_MAX_CEX_LEN];
    size_t   cex_len;
    uint16_t violation_event;
    uint16_t violation_state;
    int      is_controllable;
} controllability_result_t;

int controllability_check(const des_automaton_t *G, const des_automaton_t *K,
                           controllability_result_t *result);
int controllability_check_language(const des_automaton_t *G,
    const uint16_t *K_strings, const uint8_t *K_lens, uint16_t K_count,
    controllability_result_t *result);
void controllability_compute_supremal(const des_automaton_t *G,
    const des_automaton_t *K, des_automaton_t *supC);
int controllability_is_sublanguage(const des_automaton_t *H,
    const des_automaton_t *G, int check_marked);
int controllability_observer_property(const des_automaton_t *G);
uint16_t controllability_bad_states(const des_automaton_t *G,
    const des_automaton_t *K, uint16_t *bad_states, uint16_t cap);
int controllability_synthesize(const des_automaton_t *G,
    const des_automaton_t *K, des_automaton_t *supervisor_automaton);
void controllability_print_result(const controllability_result_t *result,
                                   FILE *fp);

#endif /* CONTROLLABILITY_H */
