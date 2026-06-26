#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "switched_types.h"
#include "switched_stability.h"
#include "switched_lyapunov.h"

int main(void) {
    printf("=== Switched Stability Benchmarks ===\n\n");
    int n_iter = 1000;
    clock_t start, end;

    SwitchedMatrix A = sm_create(3, 3);
    sm_set(&A, 0, 0, -1.0); sm_set(&A, 0, 1, 0.5); sm_set(&A, 0, 2, 0.1);
    sm_set(&A, 1, 0, -0.5); sm_set(&A, 1, 1, -2.0); sm_set(&A, 1, 2, 0.3);
    sm_set(&A, 2, 0, 0.1); sm_set(&A, 2, 1, -0.3); sm_set(&A, 2, 2, -3.0);

    /* Benchmark 1: Hurwitz check */
    start = clock();
    for (int i = 0; i < n_iter; i++) { bool h = is_hurwitz_matrix(&A); (void)h; }
    end = clock();
    printf("  Hurwitz check (3x3): %.2f us/call\n",
           (double)(end-start)/CLOCKS_PER_SEC/n_iter*1e6);

    /* Benchmark 2: Lyapunov solver */
    SwitchedMatrix Q = sm_create(3, 3);
    sm_identity(&Q);
    SwitchedMatrix P = sm_create(3, 3);
    start = clock();
    for (int i = 0; i < n_iter; i++) { sss_lyap_solve(&A, &Q, &P, 3); }
    end = clock();
    printf("  Lyapunov solve (3x3): %.2f us/call\n",
           (double)(end-start)/CLOCKS_PER_SEC/n_iter*1e6);

    /* Benchmark 3: CLF gradient descent */
    SwitchedSubsystem *sub = ssub_create(0, 3);
    sm_copy(&sub->A, &A);
    ssub_set_hurwitz_check(sub);
    SwitchedSubsystem *modes[1] = {sub};
    CommonLyapunovFunction clf;
    clf.P = sm_create(3, 3);
    start = clock();
    for (int i = 0; i < 100; i++) {
        sss_clf_gradient_descent(modes, 1, 3, &clf, 100);
    }
    end = clock();
    printf("  CLF gradient descent: %.2f ms/call\n",
           (double)(end-start)/CLOCKS_PER_SEC/100*1e3);

    printf("\n=== Benchmarks Complete ===\n");
    sm_free(&A); sm_free(&Q); sm_free(&P); sm_free(&clf.P);
    ssub_free(sub);
    return 0;
}