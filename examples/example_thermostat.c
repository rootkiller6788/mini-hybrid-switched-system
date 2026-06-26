/**
 * @file example_thermostat.c
 * @brief L6 KP2: Thermostat with Hysteresis — Hybrid Control Example
 *
 * Demonstrates a two-mode ON/OFF control system with hysteresis.
 * The thermostat maintains room temperature within a deadband
 * by switching the heater on and off.
 *
 * This is the canonical example of a bang-bang controller with
 * state-dependent switching. It exhibits a stable limit cycle:
 * temperature oscillates between T_low and T_high.
 */
#include "hss_simulation.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Thermostat with Hysteresis (L6) ===\n\n");

    HSS_Thermostat t = hss_thermostat_init(
        20.0,    /* T_start: 20°C */
        18.0,    /* T_low:   heat ON below 18°C */
        22.0,    /* T_high:  heat OFF above 22°C */
        0.8,     /* alpha:   heating rate coefficient */
        0.3,     /* beta:    cooling rate coefficient */
        35.0,    /* T_hot:   ambient when heating */
        5.0      /* T_cold:  ambient when cooling */
    );

    printf("Parameters:\n");
    printf("  T_low = %.1f°C, T_high = %.1f°C (deadband = %.1f°C)\n",
           t.T_low, t.T_high, t.T_high - t.T_low);
    printf("  Heating: α = %.2f, T_ambient = %.1f°C\n",
           t.alpha, t.T_ambient_hot);
    printf("  Cooling: β = %.2f, T_ambient = %.1f°C\n",
           t.beta, t.T_ambient_cold);
    printf("  Initial: T = %.1f°C, heater = %s\n\n",
           t.temperature, t.heater_on ? "ON" : "OFF");

    double dt = 0.01;
    double sim_time = 20.0;
    double print_dt = 1.0;
    double next_print = 0.0;
    int switch_count = 0;

    printf("Simulation (dt=%.3f s):\n", dt);
    printf("  %8s  %10s  %8s  %10s\n", "Time(s)", "Temp(°C)", "Heater", "Duty");

    while (t.time < sim_time) {
        int sw = hss_thermostat_step(&t, dt);
        if (sw != 0) switch_count++;

        if (t.time >= next_print) {
            printf("  %8.2f  %10.3f  %8s  %10.3f\n",
                   t.time, t.temperature,
                   t.heater_on ? "ON " : "OFF",
                   t.duty_cycle);
            next_print += print_dt;
        }
    }

    printf("\nResults:\n");
    printf("  Switch events: %d\n", switch_count);
    printf("  Average temperature: %.3f°C\n", t.avg_temperature);
    printf("  Duty cycle: %.3f\n", t.duty_cycle);
    printf("  Temperature variation: %.3f°C (deadband %.1f°C)\n",
           fabs(t.temperature - (t.T_high + t.T_low) / 2.0),
           t.T_high - t.T_low);

    return 0;
}
