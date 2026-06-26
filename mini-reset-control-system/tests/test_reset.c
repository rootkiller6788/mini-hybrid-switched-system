/* test_reset.c - Tests for reset core and system */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "reset_core.h"
#include "reset_system.h"

static int test_passed = 0, test_failed = 0;
#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); test_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); test_failed++; } while(0)
#define CHECK(cond) do { if (!(cond)) { FAIL(#cond); return; } } while(0)

static void test_base_create(void) {
    TEST("base_create");
    ResetLinearBase *b = reset_base_create(2, 1, 1);
    CHECK(b != NULL);
    CHECK(b->n == 2 && b->m == 1 && b->p == 1);
    CHECK(b->A != NULL && b->B != NULL && b->C != NULL && b->D != NULL);
    reset_base_free(b);
    PASS();
}

static void test_base_clone(void) {
    TEST("base_clone");
    ResetLinearBase *b = reset_base_create(2, 1, 1);
    CHECK(b != NULL);
    b->A[0] = -1.0; b->A[3] = -2.0;
    b->B[0] = 1.0; b->B[1] = 0.0;
    b->C[0] = 1.0;
    ResetLinearBase *c = reset_base_clone(b);
    CHECK(c != NULL);
    CHECK(c->n == 2);
    CHECK(fabs(c->A[0] - (-1.0)) < 1e-10);
    CHECK(fabs(c->B[0] - 1.0) < 1e-10);
    reset_base_free(b);
    reset_base_free(c);
    PASS();
}

static void test_clegg_create(void) {
    TEST("clegg_create");
    ResetSystem *rsys = reset_clegg_create();
    CHECK(rsys != NULL);
    CHECK(rsys->nc == 1);
    CHECK(rsys->flow != NULL);
    CHECK(fabs(rsys->flow->A[0] - 0.0) < 1e-10);
    CHECK(fabs(rsys->flow->B[0] - 1.0) < 1e-10);
    CHECK(fabs(rsys->flow->C[0] - 1.0) < 1e-10);
    CHECK(fabs(rsys->ratio->rho - 0.0) < 1e-10);
    reset_sys_free(rsys);
    PASS();
}

static void test_fore_create(void) {
    TEST("fore_create");
    ResetSystem *rsys = reset_fore_create(2.0, 0.5, 0.3);
    CHECK(rsys != NULL);
    CHECK(rsys->nc == 1);
    CHECK(fabs(rsys->flow->A[0] - (-2.0)) < 1e-10);
    CHECK(fabs(rsys->flow->B[0] - 4.0) < 1e-10);
    CHECK(fabs(rsys->ratio->rho - 0.3) < 1e-10);
    reset_sys_free(rsys);
    PASS();
}

static void test_fore_create_invalid(void) {
    TEST("fore_create_invalid");
    CHECK(reset_fore_create(1.0, -0.1, 0.5) == NULL);
    CHECK(reset_fore_create(1.0, 0.5, 1.5) == NULL);
    CHECK(reset_fore_create(1.0, 0.5, -0.1) == NULL);
    PASS();
}

static void test_reset_zc_detect(void) {
    TEST("reset_zc_detect");
    ResetSystem *rsys = reset_clegg_create();
    CHECK(rsys != NULL);
    rsys->x_c[0] = 5.0;
    ResetResult r = reset_check_and_apply(rsys, 1.0, -0.5);
    CHECK(r == RESET_OK);
    CHECK(fabs(rsys->x_c[0] - 0.0) < 1e-10);
    CHECK(rsys->interval->n_resets == 1);
    /* No ZC: same sign */
    r = reset_check_and_apply(rsys, 2.0, 1.0);
    CHECK(r == RESET_NO_CROSSING);
    reset_sys_free(rsys);
    PASS();
}

static void test_reset_dwell_time(void) {
    TEST("reset_dwell_time");
    ResetSystem *rsys = reset_clegg_create();
    CHECK(rsys != NULL);
    reset_set_dwell_time(rsys, 1.0);
    rsys->x_c[0] = 3.0;
    ResetResult r = reset_check_and_apply(rsys, 1.0, -0.5);
    CHECK(r == RESET_OK);
    /* Immediate second reset blocked by dwell */
    rsys->x_c[0] = 2.0;
    r = reset_check_and_apply(rsys, 0.5, -0.3);
    CHECK(r == RESET_DWELL_BLOCKED);
    reset_sys_free(rsys);
    PASS();
}

static void test_reset_band(void) {
    TEST("reset_band");
    ResetSystem *rsys = reset_clegg_create();
    CHECK(rsys != NULL);
    reset_config_zc_trigger(rsys, RESET_TRIGGER_ZERO_CROSSING, 0.1, 1.0, 0.01);
    /* Error outside band: not armed */
    rsys->e_prev = 5.0;
    rsys->reset_armed = false; /* high error disarms */
    ResetResult r = reset_check_and_apply(rsys, -0.5, 5.0);
    CHECK(r == RESET_NO_CROSSING || r == RESET_NOT_ARMED);
    reset_sys_free(rsys);
    PASS();
}

static void test_reset_manual(void) {
    TEST("reset_manual");
    ResetSystem *rsys = reset_clegg_create();
    CHECK(rsys != NULL);
    rsys->x_c[0] = 10.0;
    ResetResult r = reset_apply_manual(rsys, 1.0, 0.0);
    CHECK(r == RESET_OK);
    CHECK(fabs(rsys->x_c[0]) < 1e-10);
    CHECK(rsys->interval->n_resets == 1);
    reset_sys_free(rsys);
    PASS();
}

static void test_hurwitz_check(void) {
    TEST("hurwitz_check");
    ResetLinearBase *b = reset_base_create(2, 1, 1);
    CHECK(b != NULL);
    /* A = [-1, 0; 0, -2] is Hurwitz */
    b->A[0] = -1.0; b->A[3] = -2.0;
    CHECK(reset_is_hurwitz(b));
    /* A = [1, 0; 0, -2] is NOT Hurwitz */
    b->A[0] = 1.0;
    CHECK(!reset_is_hurwitz(b));
    reset_base_free(b);
    PASS();
}

int main(void) {
    printf("=== Reset Control System Tests ===\n\n");
    test_base_create();
    test_base_clone();
    test_clegg_create();
    test_fore_create();
    test_fore_create_invalid();
    test_reset_zc_detect();
    test_reset_dwell_time();
    test_reset_band();
    test_reset_manual();
    test_hurwitz_check();
    printf("\n=== Results: %d passed, %d failed ===\n", test_passed, test_failed);
    return test_failed > 0 ? 1 : 0;
}
