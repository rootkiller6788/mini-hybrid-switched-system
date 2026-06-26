#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "switched_types.h"
#include "switched_applications.h"

/**
 * Example 1: DC-DC Boost Converter Stability Analysis (L6)
 *
 * Demonstrates how a boost converter is a switched linear system
 * with two modes (ON/OFF), and how dwell-time analysis applies.
 *
 * The boost converter topology:
 *   ON mode:  inductor charges from Vin
 *   OFF mode: inductor discharges to output through diode
 *
 * This example:
 *   1. Creates a 12V-to-24V boost converter
 *   2. Computes the A matrices for each switching mode
 *   3. Analyzes stability via spectral radius
 *   4. Computes ripple and equilibrium
 *   5. Demonstrates the switched system model
 */
int main(void) {
    printf("=== Example 1: DC-DC Boost Converter Stability ===\n\n");

    /* Boost converter: 12V -> 24V, 100 kHz, typical values */
    double Vin = 12.0, Vout = 24.0;
    double L = 100e-6;   /* 100 uH */
    double C = 220e-6;   /* 220 uF */
    double R_load = 10.0; /* 10 Ohm */
    double freq = 100e3;  /* 100 kHz */

    printf("Parameters: Vin=%.1fV, Vout=%.1fV, L=%.0fuH, C=%.0fuF, R=%.1fOhm, f=%.0fkHz\n",
           Vin, Vout, L*1e6, C*1e6, R_load, freq/1000.0);

    DCDCConverter *conv = dcdc_create(Vin, Vout, L, C, R_load, freq);
    printf("Duty cycle: D = 1 - Vin/Vout = %.4f\n", conv->duty_cycle);

    /* Get A matrices */
    SwitchedMatrix A_on = sm_create(2, 2);
    SwitchedMatrix A_off = sm_create(2, 2);
    dcdc_get_matrices(conv, &A_on, &A_off);

    printf("\nON mode A matrix:\n");
    sm_print(&A_on, stdout);

    printf("\nOFF mode A matrix:\n");
    sm_print(&A_off, stdout);

    /* Spectral radius */
    double rho_on = spectral_radius(&A_on);
    double rho_off = spectral_radius(&A_off);
    printf("\nSpectral radius: ON=%.6f, OFF=%.6f\n", rho_on, rho_off);

    /* Hurwitz check */
    bool hurwitz_on = is_hurwitz_matrix(&A_on);
    bool hurwitz_off = is_hurwitz_matrix(&A_off);
    printf("Hurwitz: ON=%s, OFF=%s\n",
           hurwitz_on ? "YES" : "NO", hurwitz_off ? "YES" : "NO");

    /* Equilibrium */
    double iL_ss, vC_ss;
    dcdc_equilibrium(conv, 0.5, &iL_ss, &vC_ss);
    printf("\nEquilibrium at D=0.5: iL_ss=%.4f A, vC_ss=%.4f V\n", iL_ss, vC_ss);

    /* Ripple */
    double ripple_iL, ripple_vC;
    dcdc_compute_ripple(conv, 0.5, &ripple_iL, &ripple_vC);
    printf("Ripple: delta_iL=%.4f A, delta_vC=%.4f V\n", ripple_iL, ripple_vC);

    /* Simulate */
    printf("\nSimulating for 1 ms...\n");
    dcdc_simulate(conv, 0.001, 1e-7);
    printf("Final state: iL=%.4f A, vC=%.4f V\n", conv->iL, conv->vC);

    /* Stability analysis */
    printf("\nStability analysis: ");
    if (dcdc_analyze_stability(conv)) {
        printf("STABLE (averaged model is stable)\n");
    } else {
        printf("UNSTABLE or INCONCLUSIVE\n");
    }

    /* Demonstrate as switched system */
    SwitchedSystem *sys = sss_create("Boost Converter", 2, 2);
    sss_add_subsystem(sys, 0, &A_on);
    sss_add_subsystem(sys, 1, &A_off);

    /* Check Lie-algebraic condition */
    LieAlgebraCondition *la = lie_check_create(2, 2);
    SwitchedMatrix *A_arr[2] = {&A_on, &A_off};
    lie_condition_solvable(la, A_arr, 2, 2);
    printf("Lie algebra: solvable=%s, commute=%s\n",
           la->is_solvable ? "YES" : "NO",
           la->pair_commute ? "YES" : "NO");

    lie_check_free(la);
    sss_free(sys);
    sm_free(&A_on);
    sm_free(&A_off);
    dcdc_free(conv);

    printf("\n=== Example 1 Complete ===\n");
    return 0;
}