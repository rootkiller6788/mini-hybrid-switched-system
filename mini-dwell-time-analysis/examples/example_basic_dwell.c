#include "dta_core.h"
#include "dta_switch_signal.h"
#include "dta_stability.h"
#include <stdio.h>
#include <stdlib.h>

/* example_basic_dwell.c - Basic dwell-time analysis example
 *
 * Demonstrates:
 *   1. Creating a 2-mode switched linear system
 *   2. Analyzing stability under a given dwell time
 *   3. Computing the minimum required dwell time
 *   4. Generating a dwell-time-constrained switching signal
 */

int main(void) {
    printf("=== Dwell-Time Analysis: Basic Example ===\n\n");

    /* Create a 2-mode, 2-state switched linear system */
    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 0, 0);

    /* Mode 1: stable focus (complex eigenvalues with negative real part) */
    double A1[] = {-0.5, 2.0, -2.0, -0.5};
    /* Mode 2: stable node (real negative eigenvalues) */
    double A2[] = {-3.0, 0.1, -0.1, -1.0};

    dta_system_set_linear_mode(sys, 0, A1, NULL, NULL);
    dta_system_set_linear_mode(sys, 1, A2, NULL, NULL);

    printf("System: 2 modes, 2 states\n");
    printf("  Mode 0: stable focus (eig: -0.5 +/- 2i)\n");
    printf("  Mode 1: stable node (eig: ~-3, ~-1)\n\n");

    /* Analyze stability for dwell time tau_d = 2.0 seconds */
    double tau_d = 2.0;
    DTA_DwellStabilityResult res = dta_analyze_dwell_stability(sys, tau_d);

    printf("Dwell time analysis (tau_d = %.1f):\n", tau_d);
    printf("  Theoretical min tau_d*: %.4f\n", res.tau_d_star);
    printf("  Coupling constant mu:   %.4f\n", res.mu);
    printf("  Decay rate:            %.6f\n", res.decay_rate);
    printf("  Verdict:               %s\n",
           res.verdict == DTA_GUES ? "GLOBALLY UNIFORMLY EXP. STABLE" :
           res.verdict == DTA_INCONCLUSIVE ? "INCONCLUSIVE" : "OTHER");

    if (res.lambda_i) free(res.lambda_i);
    if (res.common_P) free(res.common_P);

    /* Generate a dwell-time-constrained switching signal */
    printf("\nGenerating switching signal with tau_d = %.1f...\n", tau_d);
    DTA_SwitchingSignal* sig = dta_signal_constant_dwell(0.0, 10.0, 2, tau_d);

    printf("  Switches: %d\n", sig->n_switches - 1);
    printf("  Min dwell: %.4f\n", dta_signal_min_dwell(sig));

    DTA_SignalStatistics stats = dta_signal_statistics(sig, 2);
    printf("  Avg dwell: %.4f\n", stats.avg_dwell);
    printf("  Switch frequency: %.4f Hz\n", stats.switch_frequency);

    dta_signal_free(sig);
    dta_system_free(sys);

    printf("\n=== Example Complete ===\n");
    return 0;
}
