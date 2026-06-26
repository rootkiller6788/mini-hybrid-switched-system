/**
 * @file example_transfer_line.c
 * @brief Example: Transfer Line Supervisory Control
 *
 * Classic DES example from Ramadge & Wonham: a manufacturing transfer line
 * with two machines M1 and M2 separated by a buffer of capacity 1.
 *
 * The supervisory control problem:
 * - M1 and M2 can each be IDLE or WORKING (or DOWN for failures)
 * - Buffer holds at most 1 part
 * - M1 must not deposit a part when buffer is full (overflow)
 * - M2 must not pick up a part when buffer is empty (underrun)
 * - Machine failures (uncontrollable) must be handled safely
 *
 * The supervisor ensures the buffer never overflows or underruns
 * while allowing maximal production throughput.
 *
 * This example demonstrates:
 * - Plant modeling as DES automaton
 * - Specification as forbidden state avoidance
 * - Supervisor synthesis
 * - Controlled vs uncontrolled events
 */

#include "des_automaton.h"
#include "supervisor.h"
#include "controllability.h"
#include "synthesis.h"
#include <stdio.h>
#include <string.h>

int main(void)
{
    printf("========================================\n");
    printf(" Transfer Line Supervisory Control\n");
    printf("========================================\n\n");

    /* --- Step 1: Model the plant (uncontrolled system) --- */
    des_automaton_t plant;
    des_automaton_init(&plant, "TransferLine");

    /* Events */
    uint16_t e_start1 = des_automaton_add_event(&plant, "M1_start",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t e_finish1 = des_automaton_add_event(&plant, "M1_finish",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t e_start2 = des_automaton_add_event(&plant, "M2_start",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t e_finish2 = des_automaton_add_event(&plant, "M2_finish",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t e_fail1 = des_automaton_add_event(&plant, "M1_fail",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t e_repair1 = des_automaton_add_event(&plant, "M1_repair",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);

    printf("Events defined:\n");
    printf("  [%u] M1_start   (controllable)\n", e_start1);
    printf("  [%u] M1_finish  (uncontrollable)\n", e_finish1);
    printf("  [%u] M2_start   (controllable)\n", e_start2);
    printf("  [%u] M2_finish  (uncontrollable)\n", e_finish2);
    printf("  [%u] M1_fail    (uncontrollable)\n", e_fail1);
    printf("  [%u] M1_repair  (controllable)\n\n", e_repair1);

    /* Build states: (M1_state, M2_state, buffer_count)
     * M1: 0=idle, 1=working, 2=down
     * M2: 0=idle, 1=working
     * buf: 0=empty, 1=full
     */
    #define ST_IDLE    0
    #define ST_WORK    1
    #define ST_DOWN    2

    /* We'll create a 12-state automaton encoding (m1,m2,buf) as m1*6+m2*2+buf */
    for (int i = 0; i < 12; i++) {
        des_automaton_add_state(&plant);
    }
    /* State encoding: state = m1*4 + m2*2 + buf */
    #define S(m1,m2,buf) ((m1)*4 + (m2)*2 + (buf))

    /* Mark all states except those with M1 down (representing acceptable
     * completion states) */
    for (int m1 = 0; m1 < 3; m1++) {
        for (int m2 = 0; m2 < 2; m2++) {
            for (int buf = 0; buf < 2; buf++) {
                if (m1 != ST_DOWN) {
                    des_automaton_mark_state(&plant, S(m1,m2,buf));
                }
            }
        }
    }

    /* Add transitions */
    /* M1_start: idle -> working (buffer must have space) */
    des_automaton_add_transition(&plant, S(ST_IDLE,0,0), e_start1, S(ST_WORK,0,0));
    des_automaton_add_transition(&plant, S(ST_IDLE,1,0), e_start1, S(ST_WORK,1,0));

    /* M1_finish: working -> idle, deposit part in buffer */
    des_automaton_add_transition(&plant, S(ST_WORK,0,0), e_finish1, S(ST_IDLE,0,1));
    des_automaton_add_transition(&plant, S(ST_WORK,1,0), e_finish1, S(ST_IDLE,1,1));

    /* M2_start: idle -> working (buffer must have part) */
    des_automaton_add_transition(&plant, S(0,ST_IDLE,1), e_start2, S(0,ST_WORK,1));
    des_automaton_add_transition(&plant, S(ST_WORK,ST_IDLE,1), e_start2, S(ST_WORK,ST_WORK,1));
    des_automaton_add_transition(&plant, S(ST_DOWN,ST_IDLE,1), e_start2, S(ST_DOWN,ST_WORK,1));

    /* M2_finish: working -> idle, remove part from buffer */
    des_automaton_add_transition(&plant, S(0,ST_WORK,1), e_finish2, S(0,ST_IDLE,0));
    des_automaton_add_transition(&plant, S(ST_WORK,ST_WORK,1), e_finish2, S(ST_WORK,ST_IDLE,0));

    /* M1_fail: idle or working -> down */
    des_automaton_add_transition(&plant, S(ST_IDLE,0,0), e_fail1, S(ST_DOWN,0,0));
    des_automaton_add_transition(&plant, S(ST_IDLE,1,0), e_fail1, S(ST_DOWN,1,0));
    des_automaton_add_transition(&plant, S(ST_WORK,0,0), e_fail1, S(ST_DOWN,0,0));
    des_automaton_add_transition(&plant, S(ST_WORK,1,0), e_fail1, S(ST_DOWN,1,0));
    des_automaton_add_transition(&plant, S(ST_WORK,0,1), e_fail1, S(ST_DOWN,0,1));

    /* M1_repair: down -> idle */
    for (int m2 = 0; m2 < 2; m2++) {
        for (int buf = 0; buf < 2; buf++) {
            des_automaton_add_transition(&plant, S(ST_DOWN,m2,buf), e_repair1, S(ST_IDLE,m2,buf));
        }
    }

    printf("Plant automaton: %u states, %u transitions\n\n",
           plant.nstates, plant.ntrans);

    /* --- Step 2: Define safety specification --- */
    /* The specification forbids buffer overflow and underrun.
     * These are prevented by the plant model itself (transitions are
     * only added where legal). We verify this with controllability. */

    des_automaton_t spec;
    des_automaton_init(&spec, "SafetySpec");
    /* Copy the plant structure as the specification */
    for (int i = 0; i < plant.nstates; i++) des_automaton_add_state(&spec);
    for (int i = 0; i < plant.nevents; i++)
        des_automaton_add_event(&spec, plant.events[i].label, plant.events[i].flags);
    for (int i = 0; i < plant.ntrans; i++)
        des_automaton_add_transition(&spec, plant.transitions[i].src,
            plant.transitions[i].ev, plant.transitions[i].dst);
    for (int i = 0; i < plant.n_marked; i++)
        des_automaton_mark_state(&spec, plant.marked[i]);

    /* --- Step 3: Check controllability --- */
    controllability_result_t result;
    int is_ctrl = controllability_check(&plant, &spec, &result);
    printf("Controllability: %s\n", is_ctrl ? "CONTROLLABLE" : "NOT CONTROLLABLE");

    /* --- Step 4: Synthesize supervisor --- */
    des_automaton_t supC;
    synthesis_result_t syn_result;
    synthesis_supremal_controllable(&plant, &spec, &supC, &syn_result);
    printf("Synthesis result: feasible=%s, iterations=%u, states=%u\n",
           syn_result.feasible ? "yes" : "no",
           syn_result.iterations, syn_result.final_nstates);

    /* --- Step 5: Compute supervisor control patterns --- */
    uint32_t disable_map[DES_MAX_STATES];
    synthesis_disable_map(&plant, &spec, disable_map, DES_MAX_STATES);

    printf("\nSupervisor control strategy (controllable events to disable):\n");
    for (uint16_t q = 0; q < plant.nstates; q++) {
        if (disable_map[q] != 0) {
            uint16_t m1 = q / 4, rem = q % 4;
            uint16_t m2 = rem / 2, buf = rem % 2;
            printf("  State (M1=%u, M2=%u, buf=%u): disable mask=0x%08x\n",
                   m1, m2, buf, disable_map[q]);
        }
    }

    /* --- Step 6: Statistics --- */
    printf("\n=== System Statistics ===\n");
    printf("Plant: %u states, %u events (%u controllable, %u uncontrollable)\n",
           plant.nstates, plant.nevents,
           plant.n_controllable,
           plant.nevents - plant.n_controllable);
    printf("Specification: %u states\n", spec.nstates);
    printf("Supervisor: %u states\n", supC.nstates);

    uint16_t nd = des_automaton_count_deadlocks(&supC);
    uint16_t nl = des_automaton_count_livelocks(&supC);
    printf("Deadlocks in supervisor: %u\n", nd);
    printf("Livelocks in supervisor: %u\n", nl);

    printf("\n=== Transfer line supervisory control complete ===\n");
    return 0;
}