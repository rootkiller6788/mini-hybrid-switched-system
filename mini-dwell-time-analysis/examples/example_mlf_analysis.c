#include "dta_core.h"
#include "dta_mlf.h"
#include <stdio.h>
#include <stdlib.h>

/* example_mlf_analysis.c - Multiple Lyapunov Function Analysis
 *
 * Demonstrates:
 *   1. Constructing quadratic MLF from a switched system
 *   2. Computing the coupling constant mu
 *   3. Computing the required dwell time from MLF
 *   4. Evaluating Lyapunov function values
 */

int main(void) {
    printf("=== Multiple Lyapunov Function Analysis ===\n\n");

    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 0, 0);
    double A1[] = {-1.5, 0.5, -1.0, -1.0};
    double A2[] = {-2.0, 0.0, 0.5, -3.0};
    dta_system_set_linear_mode(sys, 0, A1, NULL, NULL);
    dta_system_set_linear_mode(sys, 1, A2, NULL, NULL);

    printf("Constructing quadratic MLF via Lyapunov equations...\n");
    DTA_MultipleLyapunov* mlf = dta_mlf_construct_quadratic(sys, NULL);

    if (mlf) {
        printf("  MLF constructed: %d modes, state dim %d\n", mlf->n_modes, mlf->n);
        printf("  Coupling constant mu: %.6f\n", mlf->mu);

        double tau_d_req = dta_mlf_required_dwell(mlf);
        printf("  Required dwell time tau_d*: %.6f\n", tau_d_req);

        /* Evaluate V_i at sample points */
        double x1[] = {1.0, 0.0};
        double x2[] = {0.0, 1.0};
        double x3[] = {1.0, 1.0};

        printf("\n  V_0 at (1,0): %.6f\n", dta_mlf_evaluate(mlf, 0, x1));
        printf("  V_0 at (0,1): %.6f\n", dta_mlf_evaluate(mlf, 0, x2));
        printf("  V_0 at (1,1): %.6f\n", dta_mlf_evaluate(mlf, 0, x3));
        printf("  V_1 at (1,0): %.6f\n", dta_mlf_evaluate(mlf, 1, x1));
        printf("  V_1 at (1,1): %.6f\n", dta_mlf_evaluate(mlf, 1, x3));

        /* Verify MLF conditions */
        DTA_MLFVerification v = dta_mlf_verify(mlf, sys, NULL, x1, tau_d_req);
        printf("\n  MLF Verification:\n");
        printf("    All P_i > 0:       %s\n", v.all_pd ? "yes" : "no");
        printf("    All V_dot_i < 0:   %s\n", v.all_decrease ? "yes" : "no");
        printf("    Overall valid:     %s\n", v.valid ? "yes" : "no");

        dta_mlf_free(mlf);
    }

    dta_system_free(sys);
    printf("\n=== Example Complete ===\n");
    return 0;
}
