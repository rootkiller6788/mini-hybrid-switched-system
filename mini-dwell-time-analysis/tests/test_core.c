#include "dta_core.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

/* test_core.c - Tests for dta_core.h API
 * Tests: system lifecycle, eigenvalue computation, Hurwitz/Schur checks,
 * matrix exponential, matrix measure, Lyapunov solver, L2 gain. */

static int tests_run = 0, tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)

static double A_stable_2x2[4] = {-1.0, 1.0, -1.0, -2.0};
static double A_unstable_2x2[4] = {1.0, 0.0, 0.0, -1.0};

void test_system_create(void) {
    TEST("system create");
    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 1, 1);
    CHECK(sys != NULL, "null system");
    CHECK(sys->n_modes == 2, "wrong n_modes");
    CHECK(sys->state_dim == 2, "wrong state_dim");
    dta_system_free(sys);
    PASS();
}

void test_system_set_linear_mode(void) {
    TEST("system set linear mode");
    DTA_SwitchedSystem* sys = dta_system_create(1, 2, 1, 1);
    double B[2] = {1.0, 0.0}, C[2] = {1.0, 1.0};
    int ret = dta_system_set_linear_mode(sys, 0, A_stable_2x2, B, C);
    CHECK(ret == 0, "set_linear_mode failed");
    CHECK(sys->modes[0].A != NULL, "A not set");
    CHECK(sys->modes[0].stability == DTA_MODE_STABLE, "not detected stable");
    CHECK(sys->modes[0].max_eigenvalue_real < -0.1, "wrong max eigenvalue");
    dta_system_free(sys);
    PASS();
}

void test_eigenvalues_stable(void) {
    TEST("eigenvalues stable 2x2");
    double re[2], im[2];
    int ret = dta_eigenvalues(A_stable_2x2, 2, re, im, 200);
    CHECK(ret == 0, "eigenvalues failed");
    CHECK(re[0] < 0 && re[1] < 0, "not all negative real parts");
    PASS();
}

void test_eigenvalues_unstable(void) {
    TEST("eigenvalues unstable 2x2");
    double re[2], im[2];
    int ret = dta_eigenvalues(A_unstable_2x2, 2, re, im, 200);
    CHECK(ret == 0, "eigenvalues failed");
    CHECK(re[0] > 0 || re[1] > 0, "no positive real part detected");
    PASS();
}

void test_is_hurwitz(void) {
    TEST("is_hurwitz");
    CHECK(dta_is_hurwitz(A_stable_2x2, 2, 1e-9), "stable not detected");
    CHECK(!dta_is_hurwitz(A_unstable_2x2, 2, 1e-9), "unstable not detected");
    PASS();
}

void test_is_schur(void) {
    TEST("is_schur");
    double schur_stable[4] = {0.5, 0.1, 0.1, 0.5};
    CHECK(dta_is_schur(schur_stable, 2, 1e-9), "Schur not detected");
    PASS();
}

void test_matrix_measure(void) {
    TEST("matrix measure");
    double mu = dta_matrix_measure(A_stable_2x2, 2);
    CHECK(mu < 0, "stable matrix should have negative measure");
    PASS();
}

void test_matrix_exp(void) {
    TEST("matrix exponential");
    double result[4];
    dta_matrix_exp(A_stable_2x2, 2, 0.1, result);
    /* Check that result is not identity (i.e., computation happened) */
    CHECK(fabs(result[0] - 1.0) > 0.01 || fabs(result[3] - 1.0) > 0.01,
          "matrix exp should differ from identity for t>0");
    PASS();
}

void test_solve_lyapunov(void) {
    TEST("solve Lyapunov");
    double Q[4] = {2.0, 0.0, 0.0, 2.0};
    double P[4];
    int ret = dta_solve_lyapunov(A_stable_2x2, 2, Q, P);
    CHECK(ret == 0, "Lyapunov solve failed");
    CHECK(P[0] > 0 && P[3] > 0, "P should be positive definite");
    /* Check symmetry */
    CHECK(fabs(P[1] - P[2]) < 1e-10, "P should be symmetric");
    PASS();
}

void test_l2_gain(void) {
    TEST("L2 gain");
    double B[2] = {0.0, 1.0}, C[2] = {1.0, 0.0};
    double gain = dta_l2_gain(A_stable_2x2, B, C, 2, 1, 1);
    CHECK(isfinite(gain), "L2 gain should be finite");
    CHECK(gain > 0, "L2 gain should be positive");
    PASS();
}

void test_null_checks(void) {
    TEST("null pointer checks");
    CHECK(dta_is_hurwitz(NULL, 2, 0.0) == false, "hurwitz NULL");
    CHECK(dta_is_schur(NULL, 2, 0.0) == false, "schur NULL");
    CHECK(dta_eigenvalues(NULL, 2, NULL, NULL, 0) == -1, "eig NULL");
    CHECK(dta_matrix_measure(NULL, 0) == 0.0, "measure NULL 0-dim");
    PASS();
}

int main(void) {
    printf("=== Test: dta_core ===\n");
    test_system_create();
    test_system_set_linear_mode();
    test_eigenvalues_stable();
    test_eigenvalues_unstable();
    test_is_hurwitz();
    test_is_schur();
    test_matrix_measure();
    test_matrix_exp();
    test_solve_lyapunov();
    test_l2_gain();
    test_null_checks();
    printf("\nResults: %d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
