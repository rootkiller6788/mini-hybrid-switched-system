/**
 * @file test_hss_core.c
 * @brief Tests for HSS Core (L1-L2): system lifecycle, mode configuration,
 *        transition management, invariant setup, validation.
 */

#include "hss_core.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST: %s ... ", name); \
} while(0)

#define PASS() do { \
    tests_passed++; \
    printf("PASS\n"); \
} while(0)

#define FAIL(msg) do { \
    printf("FAIL: %s\n", msg); \
} while(0)

static void test_system_create_destroy(void) {
    TEST("system_create / system_destroy");
    HSS_System *sys = hss_system_create("test_sys", 4, 3, 2);
    assert(sys != NULL);
    assert(sys->num_modes == 4);
    assert(sys->state_dim == 3);
    assert(sys->input_dim == 2);
    assert(strcmp(sys->name, "test_sys") == 0);
    assert(sys->active_mode == -1);
    hss_system_destroy(sys);
    PASS();
}

static void test_mode_dynamics_linear(void) {
    TEST("mode_set_dynamics (linear)");
    HSS_System *sys = hss_system_create("lin", 2, 2, 1);
    assert(sys != NULL);

    double A[4] = {-1, 0, 0, -2}; /* Stable diagonal */
    double B[2] = {1, 1};
    int ret = hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, B, NULL);
    assert(ret == 0);
    assert(sys->modes[0].dynamics_class == HSS_CLASS_LINEAR);
    assert(fabs(sys->modes[0].matrix.A[0] - (-1.0)) < 1e-9);

    ret = hss_mode_set_dynamics(sys, 1, HSS_CLASS_AFFINE, A, B,
                                 (double[]){0.5, 0.0});
    assert(ret == 0);
    assert(sys->modes[1].dynamics_class == HSS_CLASS_AFFINE);
    assert(sys->modes[1].matrix.has_affine);

    hss_system_destroy(sys);
    PASS();
}

static void test_add_transition(void) {
    TEST("add_transition");
    HSS_System *sys = hss_system_create("trans_test", 3, 2, 1);
    assert(sys != NULL);

    double A[4] = {-1, 0, 0, -2};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_Guard guard;
    double guard_a[2] = {1.0, 0.0};
    guard.a = guard_a;
    guard.b = 0.5;
    guard.dim = 2;
    guard.is_active = true;

    int idx = hss_add_transition(sys, 0, 1, &guard, 1,
                                  HSS_RESET_IDENTITY, NULL, NULL,
                                  "trans_01");
    assert(idx == 0);
    assert(sys->num_transitions == 1);
    assert(sys->transitions[0].src_mode == 0);
    assert(sys->transitions[0].dst_mode == 1);
    assert(sys->transitions[0].num_guards == 1);

    /* Add another transition */
    idx = hss_add_transition(sys, 1, 0, &guard, 1,
                              HSS_RESET_LINEAR,
                              (double[]){0.5, 0, 0, 0.5}, NULL,
                              "trans_10");
    assert(idx == 1);
    assert(sys->num_transitions == 2);
    assert(sys->transitions[1].reset_type == HSS_RESET_LINEAR);

    hss_system_destroy(sys);
    PASS();
}

static void test_invariant_setup(void) {
    TEST("mode_add_invariant");
    HSS_System *sys = hss_system_create("inv_test", 1, 2, 0);
    assert(sys != NULL);

    double a[2] = {0.0, 1.0};
    int ret = hss_mode_add_invariant(sys, 0, a, 10.0);
    assert(ret == 0);
    assert(sys->modes[0].num_invariants == 1);

    /* Add second invariant */
    double a2[2] = {1.0, 0.0};
    ret = hss_mode_add_invariant(sys, 0, a2, 5.0);
    assert(ret == 0);
    assert(sys->modes[0].num_invariants == 2);

    hss_system_destroy(sys);
    PASS();
}

static void test_initial_state(void) {
    TEST("set_initial_state");
    HSS_System *sys = hss_system_create("init_test", 2, 3, 1);
    assert(sys != NULL);

    double x0[3] = {1.0, -2.0, 0.5};
    int ret = hss_set_initial_state(sys, 1, x0);
    assert(ret == 0);
    assert(sys->active_mode == 1);
    assert(fabs(sys->state.data[0] - 1.0) < 1e-9);
    assert(sys->state.time == 0.0);
    assert(sys->state.jumps == 0);

    hss_system_destroy(sys);
    PASS();
}

static void test_lyapunov_setup(void) {
    TEST("mode_set_lyapunov");
    HSS_System *sys = hss_system_create("lyap_test", 1, 2, 0);
    assert(sys != NULL);

    double P[4] = {2.0, 0.0, 0.0, 3.0};
    int ret = hss_mode_set_lyapunov(sys, 0, P, 0.5);
    assert(ret == 0);
    assert(sys->modes[0].lyapunov_P != NULL);
    assert(fabs(sys->modes[0].lyapunov_decay - 0.5) < 1e-9);

    hss_system_destroy(sys);
    PASS();
}

static void test_system_validate(void) {
    TEST("system_validate");
    HSS_System *sys = hss_system_create("valid", 2, 2, 1);
    assert(sys != NULL);

    double A[4] = {-1, 0, 0, -2};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A, NULL, NULL);

    int ret = hss_system_validate(sys);
    assert(ret == 0); /* Should be valid */

    /* Invalid transition src */
    HSS_Guard g;
    double ga[2] = {1, 0};
    g.a = ga; g.b = 0; g.dim = 2; g.is_active = true;
    hss_add_transition(sys, 0, 1, &g, 1, HSS_RESET_IDENTITY, NULL, NULL, "ok");

    ret = hss_system_validate(sys);
    assert(ret == 0); /* Still valid */

    hss_system_destroy(sys);
    PASS();
}

static void test_system_print_summary(void) {
    TEST("system_print_summary");
    HSS_System *sys = hss_system_create("print_test", 3, 2, 1);
    assert(sys != NULL);

    /* Should not crash */
    hss_system_print_summary(sys, stdout);
    hss_system_print_summary(NULL, stdout); /* NULL safety */
    hss_system_print_summary(sys, NULL);    /* NULL fp → stdout */

    hss_system_destroy(sys);
    PASS();
}

static void test_null_pointer_safety(void) {
    TEST("null pointer safety");
    /* All functions should handle NULL gracefully */
    assert(hss_system_create(NULL, 1, 1, 0) == NULL);
    assert(hss_system_create("test", 0, 1, 0) == NULL);
    assert(hss_mode_set_dynamics(NULL, 0, 0, NULL, NULL, NULL) == -1);
    assert(hss_add_transition(NULL, 0, 0, NULL, 0, 0, NULL, NULL, NULL) == -1);
    assert(hss_mode_add_invariant(NULL, 0, NULL, 0) == -1);
    assert(hss_set_initial_state(NULL, 0, NULL) == -1);
    assert(hss_system_validate(NULL) == -1);
    hss_system_destroy(NULL); /* Should not crash */
    hss_system_print_summary(NULL, NULL); /* Should not crash */
    PASS();
}

int main(void) {
    printf("=== HSS Core Tests ===\n\n");

    test_system_create_destroy();
    test_mode_dynamics_linear();
    test_add_transition();
    test_invariant_setup();
    test_initial_state();
    test_lyapunov_setup();
    test_system_validate();
    test_system_print_summary();
    test_null_pointer_safety();

    printf("\n=== Results: %d/%d tests passed ===\n",
           tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
