#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "switched_types.h"
#include "switched_stability.h"
#include "switched_lyapunov.h"
#include "switched_dwell_time.h"

int main(void) {
    printf("=== Switched Stability Overview Demo ===\n\n");

    printf("--- Demo 1: CLF for commuting matrices ---\n");
    SwitchedSystem *sys1 = sss_create("demo1", 2, 2);
    SwitchedMatrix A0 = sm_create(2, 2), A1 = sm_create(2, 2);
    sm_set(&A0, 0, 0, -1.0); sm_set(&A0, 1, 1, -2.0);
    sm_set(&A1, 0, 0, -3.0); sm_set(&A1, 1, 1, -4.0);
    sss_add_subsystem(sys1, 0, &A0);
    sss_add_subsystem(sys1, 1, &A1);
    printf("  A0 = [[-1,0],[0,-2]], A1 = [[-3,0],[0,-4]]\n");
    bool gues1 = sss_is_gues_arbitrary(sys1);
    printf("  GUES under arbitrary switching: %s\n", gues1 ? "YES" : "NO");
    sss_free(sys1); sm_free(&A0); sm_free(&A1);

    printf("\n--- Demo 2: Dwell-Time Analysis ---\n");
    SwitchedSystem *sys2 = sss_create("demo2", 2, 2);
    SwitchedMatrix B0 = sm_create(2, 2), B1 = sm_create(2, 2);
    sm_set(&B0, 0, 0, -1.0); sm_set(&B0, 0, 1, 0.5);
    sm_set(&B0, 1, 0, -0.5); sm_set(&B0, 1, 1, -2.0);
    sm_set(&B1, 0, 0, -3.0); sm_set(&B1, 0, 1, 0.2);
    sm_set(&B1, 1, 0, 0.1); sm_set(&B1, 1, 1, -4.0);
    sss_add_subsystem(sys2, 0, &B0);
    sss_add_subsystem(sys2, 1, &B1);
    DwellTimeAnalysis dta;
    sss_dwell_time_theorem(sys2, &dta);
    printf("  Stability margin: %.4f\n", dta.stability_margin);
    printf("  Minimum dwell time: %.4f\n", dta.tau_d);
    printf("  Required avg dwell: %.4f\n", dta.required_tau_a);
    sss_free(sys2); sm_free(&B0); sm_free(&B1);

    printf("\n=== Demo Complete ===\n");
    return 0;
}