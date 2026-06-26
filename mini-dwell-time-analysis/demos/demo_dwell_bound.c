#include "dta_core.h"
#include "dta_stability.h"
#include <stdio.h>
#include <stdlib.h>

/* demo_dwell_bound.c - Visualize dwell time stability bound
 * Outputs tau_d vs decay rate data for plotting. */

int main(void) {
    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 0, 0);
    double A1[] = {-1.0, 0.2, -0.2, -1.5};
    double A2[] = {-2.0, 0.1, -0.1, -0.5};
    dta_system_set_linear_mode(sys, 0, A1, NULL, NULL);
    dta_system_set_linear_mode(sys, 1, A2, NULL, NULL);

    printf("# tau_d  decay_rate  stable?\n");
    int i;
    for (i = 1; i <= 50; i++) {
        double tau = 0.1 * (double)i;
        DTA_DwellStabilityResult res = dta_analyze_dwell_stability(sys, tau);
        printf("%.2f  %.6f  %d\n", tau, res.decay_rate,
               res.verdict == DTA_GUES ? 1 : 0);
        if (res.lambda_i) free(res.lambda_i);
        if (res.common_P) free(res.common_P);
    }

    dta_system_free(sys);
    return 0;
}
