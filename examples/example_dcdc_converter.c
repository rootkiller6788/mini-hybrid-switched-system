/**
 * @file example_dcdc_converter.c
 * @brief L6 KP3: DC-DC Boost Converter — Power Electronics Switching
 *
 * Demonstrates a switched affine system modeling a DC-DC boost converter.
 *
 * Modes:
 *   ON:  L·diL/dt = Vin,         C·dvC/dt = -vC/R
 *   OFF: L·diL/dt = Vin - vC,    C·dvC/dt = iL - vC/R
 *
 * The converter alternates between ON and OFF at PWM frequency.
 * Average output: Vout ≈ Vin / (1 - D)
 */
#include "hss_simulation.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== DC-DC Boost Converter — Switched Affine System (L6) ===\n\n");

    double Vin = 12.0;
    double Vref = 24.0;
    double D = 1.0 - Vin / Vref; /* D = 0.5 for 12V→24V */

    HSS_DCDCConverter conv = hss_dcdc_init(
        Vin, 1e-3, 100e-6, 10.0, 100e3, D
    );

    printf("Parameters:\n");
    printf("  Vin = %.1f V, Vref = %.1f V, D = %.2f\n", Vin, Vref, D);
    printf("  L = %.1e H, C = %.1e F, R = %.1f Ω\n", conv.L, conv.C, conv.R);
    printf("  f_sw = %.0f kHz, T_sw = %.2f μs\n\n",
           conv.freq / 1000.0, 1e6 / conv.freq);

    double dt = 1e-7;
    double sim_time = 0.01;
    double print_interval = 0.001;
    double next_print = 0.0;

    printf("Startup transient:\n");
    printf("  %8s  %10s  %10s  %8s\n", "Time(s)", "iL(A)", "vC(V)", "State");

    while (conv.time < sim_time) {
        int sw_changed = hss_dcdc_step(&conv, dt);

        if (conv.time >= next_print) {
            printf("  %8.5f  %10.4f  %10.4f  %8s\n",
                   conv.time, conv.iL, conv.vC,
                   conv.switch_state ? "ON " : "OFF");
            next_print += print_interval;
        }
        (void)sw_changed;
    }

    printf("\nSteady-state results:\n");
    printf("  Output voltage: %.4f V (target: %.1f V)\n", conv.vC, Vref);
    printf("  Output ripple: %.4f V (%.2f%%)\n",
           conv.output_ripple,
           100.0 * conv.output_ripple / conv.vC);
    printf("  Inductor current: %.4f A\n", conv.iL);
    printf("  Ideal Vout = Vin/(1-D) = %.2f V\n",
           Vin / (1.0 - D));

    /* Efficiency check */
    double P_in = Vin * conv.iL;
    double P_out = conv.vC * conv.vC / conv.R;
    printf("  Efficiency: %.1f%%\n", 100.0 * P_out / P_in);

    return 0;
}
