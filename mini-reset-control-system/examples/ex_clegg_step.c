/* ex_clegg_step.c - Clegg Integrator Step Response Demo
 *
 * Demonstrates the Clegg integrator response to a step input.
 * The Clegg integrator resets to zero on each zero-crossing of
 * the input, producing a fundamentally different response from
 * a linear integrator.
 *
 * This example shows:
 *   - Step response of Clegg integrator vs. linear integrator
 *   - Reset events at zero-crossings of the error signal
 *   - Phase advantage of 38.15 degrees over linear integrator
 *
 * Ref: Clegg (1958), Banos & Barreiro (2012) Section 3.2
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "reset_element.h"

int main(void)
{
    printf("=== Clegg Integrator Step Response Demo ===\n");


    CleggIntegrator *ci = clegg_create();
    if (!ci) {
        printf("ERROR: Failed to create Clegg integrator\n");

        return 1;
    }

    double dt = 0.001;
    double t_end = 5.0;
    int n_steps = (int)(t_end / dt);
    double step_input = 1.0;  /* unit step */

    printf("Simulating Clegg integrator with step input = %.1f\n", step_input);

    printf("Time step = %.4f, Duration = %.1f s\n", dt, t_end);

    printf("Time(s)	Output	Resets\n");

    printf("-------	------	------\n");


    double e_prev = step_input;  /* error is the input for the integrator */
    int print_interval = (int)(0.1 / dt);
    double y = 0.0;

    for (int k = 0; k < n_steps; k++) {
        double t = k * dt;
        y = clegg_step(ci, dt, step_input, e_prev);
        e_prev = step_input;

        if (k % print_interval == 0) {
            printf("%.3f	%.4f	%d\n", t, y, clegg_num_resets(ci));

        }
    }

    printf("Final Clegg integrator output: %.4f\n", y);

    printf("Total reset events: %d\n", clegg_num_resets(ci));

    printf("Describing function phase: %.2f deg (vs. -90 deg for linear)\n",

180.0 * clegg_df_phase() / M_PI);


    /* Compare with linear integrator (no reset) */
    double y_linear = 0.0;
    for (int k = 0; k < n_steps; k++) {
        y_linear += step_input * dt;
    }
    printf("Linear integrator output (no reset): %.4f\n", y_linear);

    printf("Phase advantage of Clegg over linear: %.2f degrees\n",

90.0 + 180.0 * clegg_df_phase() / M_PI);


    clegg_free(ci);
    printf("=== Demo Complete ===\n");

    return 0;
}
