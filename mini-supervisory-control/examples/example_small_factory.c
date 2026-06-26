/**
 * @file example_small_factory.c
 * @brief Example: Small Factory (Wonham's classic example)
 *
 * Classic DES example from Wonham's supervisory control textbook.
 * Two machines M1, M2 with a shared resource (robot).
 *
 * The Small Factory problem:
 * - Two machines each process parts independently
 * - A single robot serves both machines (shared resource)
 * - Machines can break down (uncontrollable)
 * - The supervisor must prevent both machines from using the
 *   robot simultaneously (mutual exclusion)
 *
 * L1: event types, automaton
 * L2: controllability concept
 * L4: controllability theorem (K controllable <=> supervisor exists)
 * L5: supremal controllable sublanguage
 * L6: small factory canonical problem
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
    printf(" Small Factory Supervisory Control\n");
    printf(" (Ramadge-Wonham / Wonham textbook)\n");
    printf("========================================\n\n");

    /* Machine 1 automaton: IDLE -> WORKING -> IDLE, can fail */
    des_automaton_t M1;
    des_automaton_init(&M1, "Machine1");
    uint16_t m1_idle = M1.q0;
    uint16_t m1_work = des_automaton_add_state(&M1);
    uint16_t m1_down = des_automaton_add_state(&M1);

    uint16_t m1_take = des_automaton_add_event(&M1, "M1_take_robot",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t m1_release = des_automaton_add_event(&M1, "M1_release",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t m1_fail = des_automaton_add_event(&M1, "M1_fail",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t m1_repair = des_automaton_add_event(&M1, "M1_repair",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);

    des_automaton_add_transition(&M1, m1_idle, m1_take, m1_work);
    des_automaton_add_transition(&M1, m1_work, m1_release, m1_idle);
    des_automaton_add_transition(&M1, m1_idle, m1_fail, m1_down);
    des_automaton_add_transition(&M1, m1_work, m1_fail, m1_down);
    des_automaton_add_transition(&M1, m1_down, m1_repair, m1_idle);
    des_automaton_mark_state(&M1, m1_idle);
    des_automaton_mark_state(&M1, m1_work);

    /* Machine 2 automaton: analogous */
    des_automaton_t M2;
    des_automaton_init(&M2, "Machine2");
    uint16_t m2_idle = M2.q0;
    uint16_t m2_work = des_automaton_add_state(&M2);
    uint16_t m2_down = des_automaton_add_state(&M2);

    uint16_t m2_take = des_automaton_add_event(&M2, "M2_take_robot",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t m2_release = des_automaton_add_event(&M2, "M2_release",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t m2_fail = des_automaton_add_event(&M2, "M2_fail",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t m2_repair = des_automaton_add_event(&M2, "M2_repair",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);

    des_automaton_add_transition(&M2, m2_idle, m2_take, m2_work);
    des_automaton_add_transition(&M2, m2_work, m2_release, m2_idle);
    des_automaton_add_transition(&M2, m2_idle, m2_fail, m2_down);
    des_automaton_add_transition(&M2, m2_work, m2_fail, m2_down);
    des_automaton_add_transition(&M2, m2_down, m2_repair, m2_idle);
    des_automaton_mark_state(&M2, m2_idle);
    des_automaton_mark_state(&M2, m2_work);

    printf("M1: %u states, M2: %u states\n", M1.nstates, M2.nstates);

    /* Plant = M1 || M2 (synchronous product) */
    des_automaton_t plant;
    des_automaton_product(&M1, &M2, &plant);
    printf("Plant (M1||M2): %u states, %u events, %u transitions\n",
           plant.nstates, plant.nevents, plant.ntrans);

    /* Mutual exclusion specification: M1_work and M2_work must not
     * both be true simultaneously. This means the "take" events
     * are mutually exclusive.
     *
     * We build a specification automaton that forbids the state
     * where both machines are working. */

    des_automaton_t spec;
    des_automaton_init(&spec, "MutualExclusion");

    /* States of the specification: track which machine has the robot.
     * 0: robot free, 1: M1 has robot, 2: M2 has robot */
    uint16_t s_free = spec.q0;
    uint16_t s_m1 = des_automaton_add_state(&spec);
    uint16_t s_m2 = des_automaton_add_state(&spec);

    /* Copy relevant events from plant */
    uint16_t se_m1t = des_automaton_add_event(&spec, "M1_take_robot",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t se_m1r = des_automaton_add_event(&spec, "M1_release",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t se_m2t = des_automaton_add_event(&spec, "M2_take_robot",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t se_m2r = des_automaton_add_event(&spec, "M2_release",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);

    /* Transitions */
    des_automaton_add_transition(&spec, s_free, se_m1t, s_m1);
    des_automaton_add_transition(&spec, s_free, se_m2t, s_m2);
    des_automaton_add_transition(&spec, s_m1, se_m1r, s_free);
    des_automaton_add_transition(&spec, s_m2, se_m2r, s_free);

    /* Self-loops for non-conflicting events */
    des_automaton_add_transition(&spec, s_m1, se_m1t, s_m1); /* M1 can take again? No, self-loop prevents undefined */
    des_automaton_add_transition(&spec, s_m2, se_m2t, s_m2);

    /* Mark all states */
    des_automaton_mark_state(&spec, s_free);
    des_automaton_mark_state(&spec, s_m1);
    des_automaton_mark_state(&spec, s_m2);

    printf("Spec (mutex): %u states, %u events\n", spec.nstates, spec.nevents);

    /* Check controllability */
    controllability_result_t ctrl_result;
    int is_controllable = controllability_check(&plant, &spec, &ctrl_result);
    printf("\nControllability of K wrt G: %s\n",
           is_controllable ? "YES (supervisor exists)" : "NO");

    if (is_controllable) {
        /* Synthesize supervisor */
        des_automaton_t supervisor_automaton;
        int synth_ok = controllability_synthesize(&plant, &spec,
                                                   &supervisor_automaton);
        printf("Synthesis: %s\n", synth_ok ? "SUCCESS" : "FAILED");

        if (synth_ok) {
            printf("Supervisor automaton: %u states\n",
                   supervisor_automaton.nstates);

            /* Build supervisor from automaton */
            supervisor_t S;
            supervisor_init(&S, &plant, &spec, SUP_TYPE_MONOLITHIC,
                            "SmallFactorySup");

            /* Set control patterns: for each state in supervisor,
             * enable events that keep us within the legal states */
            for (uint16_t q = 0; q < supervisor_automaton.nstates; q++) {
                uint32_t mask = 0xFFFFFFFF;
                supervisor_set_pattern(&S, q, mask);
            }

            /* Verify the supervisor is proper */
            int proper = supervisor_is_proper(&S);
            printf("Supervisor proper: %s\n", proper ? "YES" : "NO");

            /* Enumerate enabled events at initial state */
            uint16_t enabled[16];
            uint16_t n_enabled = supervisor_enabled_set(&S, enabled, 16);
            printf("Enabled events at initial state: %u\n", n_enabled);
            for (uint16_t i = 0; i < n_enabled; i++) {
                printf("  [%u] %s\n", enabled[i],
                       plant.events[enabled[i]].label);
            }
        }
    } else {
        printf("No supervisor exists for this specification.\n");
        controllability_print_result(&ctrl_result, stdout);

        /* Compute supremal controllable sublanguage */
        des_automaton_t supC;
        synthesis_result_t syn_result;
        synthesis_supremal_controllable(&plant, &spec, &supC, &syn_result);

        printf("\nsupC(K,G): %u states (feasible=%s)\n",
               supC.nstates, syn_result.feasible ? "yes" : "no");
    }

    printf("\n=== Small factory supervisory control complete ===\n");
    return 0;
}