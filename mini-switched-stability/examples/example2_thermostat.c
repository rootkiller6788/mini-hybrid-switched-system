#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "switched_types.h"
#include "switched_applications.h"
#include "switched_stability.h"
#include "switched_dwell_time.h"

/**
 * Example 2: Thermostat as a Switched Control System (L6)
 *
 * Demonstrates a bang-bang thermostat with hysteresis as a
 * switched system with three modes (OFF, HEATING, COOLING).
 *
 * The hysteresis creates a minimum dwell time naturally:
 * the temperature must traverse the deadband before switching.
 *
 * This example:
 *   1. Creates a thermostat with 20 +/- 1 degree deadband
 *   2. Simulates temperature evolution
 *   3. Analyzes switching statistics
 *   4. Discusses the equivalent dwell-time analysis
 */
int main(void) {
    printf("=== Example 2: Thermostat Bang-Bang Control ===\n\n");

    /* Thermostat: setpoint 20C, deadband +/-1C */
    ThermostatSystem *thermo = thermo_create(
        20.0,   /* setpoint */
        1.0,    /* deadband */
        0.05,   /* heating rate (C/s) */
        0.03,   /* cooling rate (C/s) */
        0.001   /* ambient loss rate */
    );

    printf("Setpoint: %.1f C, Deadband: +/- %.1f C\n", thermo->setpoint, thermo->deadband);
    printf("Heating rate: %.3f C/s, Cooling rate: %.3f C/s\n",
           thermo->heating_rate, thermo->cooling_rate);

    /* Set initial temp below deadband to trigger heating */
    thermo->temp = 18.0;
    printf("Initial temperature: %.1f C\n\n", thermo->temp);

    /* Simulate for 10 minutes (600 seconds) */
    printf("Simulating 10 minutes...\n");
    thermo_simulate(thermo, 600.0, 0.1);

    /* Switched system analysis */
    printf("\n--- Switched System Analysis ---\n");

    /* Create equivalent switched linear system */
    SwitchedSystem *sys = sss_create("Thermostat", 1, 3);

    /* Mode 0: OFF, dynamics: dT/dt = -k*(T - T_set) */
    SwitchedMatrix A0 = sm_create(1, 1);
    sm_set(&A0, 0, 0, -thermo->ambient_loss);
    sss_add_subsystem(sys, 0, &A0);

    /* Mode 1: HEATING, dynamics: dT/dt = -k*(T-T_set) + r_h */
    /* (Affine, linearized around setpoint: same A, different B) */
    SwitchedMatrix A1 = sm_create(1, 1);
    sm_set(&A1, 0, 0, -thermo->ambient_loss);
    sss_add_subsystem(sys, 1, &A1);

    /* Mode 2: COOLING */
    SwitchedMatrix A2 = sm_create(1, 1);
    sm_set(&A2, 0, 0, -thermo->ambient_loss);
    sss_add_subsystem(sys, 2, &A2);

    /* All modes share the same stable dynamics (A = -k),
     * so a common Lyapunov function V = (T-20)^2 exists.
     * System is trivially GUES under arbitrary switching. */
    printf("All A matrices = [%.4f] (stable)\n", -thermo->ambient_loss);
    printf("Common Lyapunov function: V(T) = (T - %.1f)^2\n", thermo->setpoint);

    /* Dwell-time analysis */
    DwellTimeAnalysis dta;
    sss_dwell_time_theorem(sys, &dta);

    printf("\nStability margin lambda_0: %.6f\n", dta.stability_margin);
    printf("Minimum dwell time tau_d: %.4f s\n", dta.tau_d);
    printf("Actual avg dwell time tau_a: %.4f s\n",
           sdt_actual_avg_dwell(sys->signal));

    /* Since all modes share the same stable dynamics (decoupled from switching),
     * the system is exponentially stable regardless of switching pattern. */
    printf("\nVerdict: The thermostat is G.U.E.S. under arbitrary switching\n");
    printf("since all subsystems are individually exponentially stable\n");
    printf("and share a common quadratic Lyapunov function.\n");

    /* Switching frequency analysis */
    thermo_switched_analysis(thermo);

    sm_free(&A0); sm_free(&A1); sm_free(&A2);
    sss_free(sys);
    thermo_free(thermo);

    printf("\n=== Example 2 Complete ===\n");
    return 0;
}