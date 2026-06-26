/* example_impulsive_consensus.c -- Impulsive consensus in multi-agent systems.
 *
 * N agents with continuous dynamics xdot_i = 0 (drift-free integrators)
 * exchange state information and apply impulsive corrections at discrete
 * communication instants.
 *
 * Consensus protocol: x_i^+ = x_i - eps * sum_{j in N_i} (x_i - x_j)
 *
 * This demonstrates how impulsive communication achieves average consensus.
 */
#include "impulsive_types.h"
#include "impulsive_control.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N_AGENTS 5
#define MAX_STEPS 200

int main(void) {
    printf("=== Impulsive Consensus Example ===\n\n");

    int n = N_AGENTS;
    double x[N_AGENTS] = {1.0, 3.0, -2.0, 5.0, -4.0};  /* initial opinions */
    double eps = 0.3;  /* coupling strength */
    /* Delta = 0.5 communication interval */

    printf("Initial states: ");
    for (int i = 0; i < n; i++) printf("%.2f ", x[i]);
    printf("\n\n");

    double avg = 0.0;
    for (int i = 0; i < n; i++) avg += x[i];
    avg /= (double)n;
    printf("Target consensus value (average): %.4f\n\n", avg);

    double max_disp = 0.0;
    for (int k = 0; k < MAX_STEPS; k++) {
        /* Compute disagreement */
        double disp = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = i+1; j < n; j++)
                disp += (x[i] - x[j]) * (x[i] - x[j]);
        disp = sqrt(disp);
        if (disp > max_disp) max_disp = disp;

        if (k % 20 == 0)
            printf("  step %3d: ", k);
        if (k % 20 == 0)
            for (int i = 0; i < n; i++) printf("%.3f ", x[i]);
        if (k % 20 == 0) printf(" | disp=%.6f\n", disp);

        /* Apply impulsive consensus update (all-to-all communication) */
        double x_new[N_AGENTS];
        for (int i = 0; i < n; i++) {
            double correction = 0.0;
            for (int j = 0; j < n; j++)
                correction += (x[i] - x[j]);
            x_new[i] = x[i] - eps * correction;
        }
        memcpy(x, x_new, (size_t)n * sizeof(double));

        if (disp < 1e-6) {
            printf("  Consensus achieved at step %d!\n", k);
            break;
        }
    }

    printf("\nFinal states: ");
    for (int i = 0; i < n; i++) printf("%.6f ", x[i]);

    double final_avg = 0.0;
    for (int i = 0; i < n; i++) final_avg += x[i];
    final_avg /= (double)n;
    printf("\nFinal average: %.6f\n", final_avg);
    printf("Max disagreement: %.6f\n", max_disp);

    double consensus_err = 0.0;
    for (int i = 0; i < n; i++)
        consensus_err += (x[i] - final_avg) * (x[i] - final_avg);
    printf("Consensus error:  %.2e\n", sqrt(consensus_err));
    return 0;
}
