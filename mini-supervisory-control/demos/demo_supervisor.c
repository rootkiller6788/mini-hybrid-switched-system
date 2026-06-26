/**
 * @file demo_supervisor.c
 * @brief Interactive demo: Supervisory control visualization
 *
 * Demonstrates:
 * - Building a plant automaton
 * - Defining safety specification
 * - Checking controllability
 * - Synthesizing supervisor
 * - Simulating closed-loop behavior
 */
#include "des_automaton.h"
#include "supervisor.h"
#include "controllability.h"
#include "synthesis.h"
#include <stdio.h>
#include <string.h>

int main(void)
{
    printf("\n");
    printf("  ╔══════════════════════════════════════════╗\n");
    printf("  ║   Supervisory Control Demo               ║\n");
    printf("  ║   Ramadge-Wonham Framework               ║\n");
    printf("  ╚══════════════════════════════════════════╝\n\n");

    /* Build a simple machine with start/stop/fault */
    des_automaton_t plant;
    des_automaton_init(&plant, "DemoMachine");

    uint16_t idle = plant.q0;
    uint16_t running = des_automaton_add_state(&plant);
    uint16_t stopped = des_automaton_add_state(&plant);

    uint16_t ev_start = des_automaton_add_event(&plant, "START",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t ev_stop = des_automaton_add_event(&plant, "STOP",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t ev_fault = des_automaton_add_event(&plant, "FAULT",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);

    des_automaton_add_transition(&plant, idle, ev_start, running);
    des_automaton_add_transition(&plant, running, ev_stop, idle);
    des_automaton_add_transition(&plant, running, ev_fault, stopped);
    des_automaton_mark_state(&plant, idle);
    des_automaton_mark_state(&plant, running);

    printf("Plant: %s\n", plant.name);
    des_automaton_print_table(&plant, stdout);

    /* Specification: never enter stopped state */
    des_automaton_t spec;
    des_automaton_init(&spec, "NoStop");
    uint16_t s0 = spec.q0;
    uint16_t s1 = des_automaton_add_state(&spec);

    uint16_t se_start = des_automaton_add_event(&spec, "START",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t se_stop = des_automaton_add_event(&spec, "STOP",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);

    des_automaton_add_transition(&spec, s0, se_start, s1);
    des_automaton_add_transition(&spec, s1, se_stop, s0);
    des_automaton_mark_state(&spec, s0);
    des_automaton_mark_state(&spec, s1);

    printf("\nSpecification: %s\n", spec.name);

    /* Check controllability */
    controllability_result_t result;
    int ctrl = controllability_check(&plant, &spec, &result);
    printf("\nControllability: %s\n", ctrl ? "YES" : "NO");

    /* Synthesize */
    des_automaton_t supC;
    synthesis_result_t syn_res;
    synthesis_supremal_controllable(&plant, &spec, &supC, &syn_res);
    printf("Supremal controllable: %u states (feasible=%s)\n",
           supC.nstates, syn_res.feasible ? "yes" : "no");

    /* Supervisor */
    supervisor_t S;
    supervisor_init(&S, &plant, &spec, SUP_TYPE_MONOLITHIC, "DemoSupervisor");
    for (uint16_t q = 0; q < supC.nstates; q++) {
        uint32_t mask = 0xFFFFFFFF;
        supervisor_set_pattern(&S, q, mask);
    }

    printf("\n=== Supervisor Status ===\n");
    supervisor_print_status(&S, stdout);

    /* Simulate a run */
    printf("\n=== Simulation ===\n");
    supervisor_reset(&S);
    printf("Observing START... %s\n",
           supervisor_observe(&S, ev_start) ? "OK" : "VIOLATION");
    printf("Observing STOP... %s\n",
           supervisor_observe(&S, ev_stop) ? "OK" : "VIOLATION");
    printf("Final state: %u\n", S.current_state);

    printf("\nDemo complete.\n");
    return 0;
}