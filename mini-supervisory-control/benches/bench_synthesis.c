/**
 * @file bench_synthesis.c
 * @brief Benchmark: supervisor synthesis performance
 */
#include "des_automaton.h"
#include "supervisor.h"
#include "controllability.h"
#include "synthesis.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

int main(void)
{
    printf("=== Supervisor Synthesis Benchmark ===\n");
    clock_t start = clock();

    /* Build a moderate-sized plant */
    des_automaton_t G;
    des_automaton_init(&G, "bench_plant");
    uint16_t states[20];
    states[0] = G.q0;
    for (int i = 1; i < 20; i++) states[i] = des_automaton_add_state(&G);
    uint16_t ev[5];
    for (int i = 0; i < 5; i++) {
        char lbl[8]; snprintf(lbl, sizeof(lbl), "e%d", i);
        uint8_t fl = (i < 3) ? (DES_EVENT_CONTROLLABLE|DES_EVENT_OBSERVABLE)
                             : (DES_EVENT_UNCONTROLLABLE|DES_EVENT_OBSERVABLE);
        ev[i] = des_automaton_add_event(&G, lbl, fl);
    }
    for (int i = 0; i < 19; i++)
        for (int j = 0; j < 5; j++)
            des_automaton_add_transition(&G, states[i], ev[j], states[i+1]);
    des_automaton_mark_state(&G, states[19]);

    des_automaton_t K;
    des_automaton_init(&K, "bench_spec");
    uint16_t ks[20];
    ks[0] = K.q0;
    for (int i = 1; i < 20; i++) ks[i] = des_automaton_add_state(&K);
    for (int i = 0; i < 5; i++) {
        char lbl[8]; snprintf(lbl, sizeof(lbl), "e%d", i);
        uint8_t fl = DES_EVENT_CONTROLLABLE|DES_EVENT_OBSERVABLE;
        des_automaton_add_event(&K, lbl, fl);
    }
    for (int i = 0; i < 19; i++)
        for (int j = 0; j < 3; j++)
            des_automaton_add_transition(&K, ks[i], ev[j], ks[i+1]);
    des_automaton_mark_state(&K, ks[19]);

    double t_build = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("Build time: %.4f s\n", t_build);

    start = clock();
    controllability_result_t cr;
    int c = controllability_check(&G, &K, &cr);
    double t_check = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("Controllability: %s (%.4f s)\n", c?"YES":"NO", t_check);

    start = clock();
    des_automaton_t supC;
    synthesis_result_t sr;
    synthesis_supremal_controllable(&G, &K, &supC, &sr);
    double t_synth = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("Synthesis: %u states (%.4f s)\n", supC.nstates, t_synth);

    printf("Total: %.4f s\n", t_build + t_check + t_synth);
    return 0;
}