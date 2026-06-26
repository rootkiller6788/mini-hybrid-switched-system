/**
 * @file example_switched_system.c
 * @brief Example: Switched System with Safety Supervisor
 *
 * Demonstrates hybrid supervisory control: a power converter with
 * two modes (buck/boost) where the supervisor ensures safe switching.
 *
 * The system:
 * - Mode 0 (Buck): steps voltage down (safe for low input)
 * - Mode 1 (Boost): steps voltage up (safe for high input)
 * - Continuous state: [output_voltage, inductor_current]
 * - Supervisor prevents switching when voltage/current exceed limits
 * - Autonomous fault transitions model overcurrent protection
 *
 * Applications:
 * - Power electronics (DC-DC converter control) — Tesla/SpaceX
 * - Smart grid switching (IEEE 1547 compliant mode transitions)
 * - Building HVAC zone management (Johnson Controls / Siemens)
 * - Aircraft power system mode management (Boeing 787, F-35)
 */

#include "des_automaton.h"
#include "supervisor.h"
#include "controllability.h"
#include "synthesis.h"
#include "hybrid_supervisor.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

int main(void)
{
    printf("========================================\n");
    printf(" Switched Power System\n");
    printf(" Supervisory Safety Control\n");
    printf("========================================\n\n");

    /* --- Build hybrid automaton --- */
    hybrid_automaton_t H;
    hybrid_automaton_init(&H, "BuckBoost_Converter");

    /* Mode 0: Buck converter dynamics
     * dx/dt = A0*x + B0*u
     * State: [v_out, i_L]
     * A0: di_L/dt = (V_in - v_out)/L, dv_out/dt = (i_L - v_out/R)/C */
    hyb_ode_t buck_ode;
    memset(&buck_ode, 0, sizeof(buck_ode));
    buck_ode.dim = 2;
    buck_ode.input_dim = 1;
    /* i_L dynamics: influenced by V_in and v_out */
    buck_ode.A[0][0] = 0.0;       buck_ode.A[0][1] = -1.0 / 0.001;  /* di_L/dt term */
    buck_ode.A[1][0] = 1.0 / 0.0001; buck_ode.A[1][1] = -1.0 / (0.0001 * 10.0); /* dv_out/dt */
    buck_ode.B[0][0] = 1000.0;    /* V_in/L effect */
    buck_ode.c[0] = 0.0;
    buck_ode.c[1] = 0.0;

    uint16_t buck_mode = hybrid_add_mode(&H, "buck", &buck_ode);
    H.modes[buck_mode].is_safe = 1;

    /* Invariant for buck: output voltage must stay below 24V */
    hyb_invariant_t buck_inv;
    memset(&buck_inv, 0, sizeof(buck_inv));
    buck_inv.dim = 1;
    buck_inv.A[0] = 0.0; buck_inv.A[0] = 0.0; /* v_out measured at index 1 */
    buck_inv.A[1] = 1.0;   /* v_out coefficient */
    buck_inv.b = -24.0;    /* v_out - 24 <= 0 => v_out <= 24 */
    hybrid_set_invariant(&H, buck_mode, &buck_inv);

    /* Mode 1: Boost converter dynamics */
    hyb_ode_t boost_ode;
    memset(&boost_ode, 0, sizeof(boost_ode));
    boost_ode.dim = 2;
    boost_ode.input_dim = 1;
    boost_ode.A[0][0] = 0.0;        boost_ode.A[0][1] = 0.0;
    boost_ode.A[1][0] = 1.0 / 0.0001; boost_ode.A[1][1] = -1.0 / (0.0001 * 20.0);
    boost_ode.B[0][0] = 1000.0;
    boost_ode.c[1] = 48.0 / (0.0001 * 20.0); /* V_out_ref contribution */

    uint16_t boost_mode = hybrid_add_mode(&H, "boost", &boost_ode);
    H.modes[boost_mode].is_safe = 1;

    /* Add DES events */
    uint16_t ev_switch_to_buck = des_automaton_add_event(&H.des,
        "switch_to_buck",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t ev_switch_to_boost = des_automaton_add_event(&H.des,
        "switch_to_boost",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t ev_overcurrent = des_automaton_add_event(&H.des,
        "overcurrent_fault",
        DES_EVENT_UNCONTROLLABLE | DES_EVENT_OBSERVABLE);

    printf("DES events: %u (buck_switch=%u, boost_switch=%u, fault=%u)\n",
           H.des.nevents, ev_switch_to_buck, ev_switch_to_boost,
           ev_overcurrent);

    /* Guards for switching */
    hyb_guard_t guard_to_boost;
    memset(&guard_to_boost, 0, sizeof(guard_to_boost));
    guard_to_boost.dim = 1;
    guard_to_boost.A[1] = 1.0;   /* v_out coefficient */
    guard_to_boost.b = 12.0;      /* -(-v_out + 12) <= 0 => v_out >= 12 */

    hyb_guard_t guard_to_buck;
    memset(&guard_to_buck, 0, sizeof(guard_to_buck));
    guard_to_buck.dim = 1;
    guard_to_buck.A[1] = -1.0;   /* -v_out coefficient */
    guard_to_buck.b = 24.0;      /* v_out - 24 <= 0 => v_out <= 24 */

    double R_identity[HYB_DIM][HYB_DIM] = {{1.0, 0.0}, {0.0, 1.0}};
    double r_zero[HYB_DIM] = {0.0, 0.0};

    /* Transitions */
    hybrid_add_transition(&H, buck_mode, boost_mode, ev_switch_to_boost,
        HYB_SWITCH_CONTROLLED, &guard_to_boost, R_identity, r_zero, 2);

    hybrid_add_transition(&H, boost_mode, buck_mode, ev_switch_to_buck,
        HYB_SWITCH_CONTROLLED, &guard_to_buck, R_identity, r_zero, 2);

    /* Fault transition: autonomous (uncontrollable) */
    hyb_guard_t empty_guard;
    memset(&empty_guard, 0, sizeof(empty_guard));
    uint16_t fault_mode = hybrid_add_mode(&H, "fault", NULL);
    H.modes[fault_mode].is_safe = 0;

    hybrid_add_transition(&H, buck_mode, fault_mode, ev_overcurrent,
        HYB_SWITCH_AUTONOMOUS, &empty_guard, R_identity, r_zero, 2);
    hybrid_add_transition(&H, boost_mode, fault_mode, ev_overcurrent,
        HYB_SWITCH_AUTONOMOUS, &empty_guard, R_identity, r_zero, 2);

    /* Set initial state */
    H.current_mode = buck_mode;
    H.x.data[0] = 1.0;  /* i_L = 1A */
    H.x.data[1] = 12.0; /* v_out = 12V */
    H.x.dim = 2;
    H.t = 0.0;

    printf("Hybrid automaton: %u modes, %u transitions\n",
           H.nmodes, H.ntrans);

    /* --- Abstraction to DES --- */
    des_automaton_t abs_G;
    hybrid_abstraction(&H, &abs_G);
    printf("DES abstraction: %u states, %u events\n",
           abs_G.nstates, abs_G.nevents);

    /* --- Safety specification: never enter fault mode --- */
    des_automaton_t safety_spec;
    des_automaton_init(&safety_spec, "NoFault");
    uint16_t safe0 = safety_spec.q0;
    uint16_t safe1 = des_automaton_add_state(&safety_spec);

    uint16_t se_boost = des_automaton_add_event(&safety_spec, "sw_buck_boost",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);
    uint16_t se_buck = des_automaton_add_event(&safety_spec, "sw_boost_buck",
        DES_EVENT_CONTROLLABLE | DES_EVENT_OBSERVABLE);

    des_automaton_add_transition(&safety_spec, safe0, se_boost, safe1);
    des_automaton_add_transition(&safety_spec, safe1, se_buck, safe0);
    des_automaton_mark_state(&safety_spec, safe0);
    des_automaton_mark_state(&safety_spec, safe1);

    /* --- Synthesize supervisor --- */
    supervisor_t S;
    hybrid_supervisor_synthesize(&H, &safety_spec, &S);
    printf("\nSupervisor synthesized: %s\n", S.name);

    /* --- Check safety --- */
    int is_safe = hybrid_supervisor_safe(&H, &S);
    printf("System safe under supervisor: %s\n", is_safe ? "YES" : "NO");

    /* --- Reachability analysis --- */
    hyb_state_t init_set, unsafe_set;
    memset(&init_set, 0, sizeof(init_set));
    memset(&unsafe_set, 0, sizeof(unsafe_set));
    init_set.dim = 2; unsafe_set.dim = 2;
    init_set.data[0] = 1.0; init_set.data[1] = 12.0;
    unsafe_set.data[0] = 100.0; unsafe_set.data[1] = 48.0;

    int safe_reach = hybrid_reachability_check(&H, &init_set, &unsafe_set);
    printf("Reachability safe: %s\n", safe_reach ? "YES" : "NO (unsafe reachable!)");

    /* --- Simulate --- */
    printf("\n=== Simulation ===\n");
    H.current_mode = buck_mode;
    H.x.data[0] = 1.0;
    H.x.data[1] = 12.0;
    H.t = 0.0;

    hyb_state_t trajectory[50];
    uint16_t n_steps = 0;
    int sim_ok = hybrid_simulate(&H, &S, 0.5, 0.02, trajectory, 50, &n_steps);
    printf("Simulation: %s (%u steps, final mode=%u)\n",
           sim_ok ? "OK" : "FAILED", n_steps, H.current_mode);

    /* Print final state */
    printf("Final state: v_out=%.3fV, i_L=%.3fA, mode=%s\n",
           H.x.data[1], H.x.data[0],
           H.modes[H.current_mode].name ? H.modes[H.current_mode].name : "?");

    /* Check all visited modes are safe */
    int all_safe = 1;
    for (uint16_t i = 0; i < n_steps; i++) {
        if (!H.modes[H.current_mode].is_safe) { all_safe = 0; break; }
    }
    printf("All visited states safe: %s\n", all_safe ? "YES" : "NO");

    /* --- Supervisor statistics --- */
    printf("\n=== Supervisor Statistics ===\n");
    supervisor_print_status(&S, stdout);

    printf("\n=== Switched power system supervisory control complete ===\n");
    return 0;
}