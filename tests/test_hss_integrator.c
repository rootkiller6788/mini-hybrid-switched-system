/**
 * @file test_hss_integrator.c
 * @brief Tests for HSS Integrator (L3-L4): system composition,
 *        equivalence, abstraction transforms.
 */
#include "hss_integrator.h"
#include "hss_core.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

static int tests_run = 0, tests_passed = 0;
#define TEST(n) do { tests_run++; printf("  TEST: %s ... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

static void test_compose_series(void) {
    TEST("compose (series)");
    HSS_System *s1 = hss_system_create("s1", 2, 2, 1);
    HSS_System *s2 = hss_system_create("s2", 2, 2, 1);
    double A[4] = {-1, 0, 0, -2};
    hss_mode_set_dynamics(s1, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(s1, 1, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(s2, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(s2, 1, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_set_initial_state(s1, 0, (double[]){1,0});
    hss_set_initial_state(s2, 1, (double[]){0,1});

    HSS_Composition *comp = hss_compose(HSS_COMP_SERIES, s1, s2);
    assert(comp != NULL);
    assert(comp->product != NULL);
    assert(comp->product->state_dim >= s1->state_dim + s2->state_dim
           || comp->product->state_dim >= 4);

    hss_composition_free(comp);
    hss_system_destroy(s1);
    hss_system_destroy(s2);
    PASS();
}

static void test_compose_parallel(void) {
    TEST("compose (parallel)");
    HSS_System *s1 = hss_system_create("p1", 1, 1, 0);
    HSS_System *s2 = hss_system_create("p2", 1, 1, 0);
    double A[1] = {-1};
    hss_mode_set_dynamics(s1, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(s2, 0, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_Composition *comp = hss_compose(HSS_COMP_PARALLEL, s1, s2);
    assert(comp != NULL);
    assert(comp->product->state_dim == 2); /* 1+1 */

    hss_composition_free(comp);
    hss_system_destroy(s1);
    hss_system_destroy(s2);
    PASS();
}

static void test_equivalence(void) {
    TEST("are_equivalent");
    HSS_System *s1 = hss_system_create("eq1", 2, 2, 0);
    HSS_System *s2 = hss_system_create("eq2", 2, 2, 0);
    double A[4] = {-1, 0, 0, -2};
    hss_mode_set_dynamics(s1, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(s1, 1, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(s2, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(s2, 1, HSS_CLASS_LINEAR, A, NULL, NULL);

    assert(hss_are_equivalent(s1, s2, HSS_EQUIV_TRACE));
    assert(hss_are_equivalent(s1, s2, HSS_EQUIV_BISIMULATION));
    assert(!hss_are_equivalent(NULL, s2, HSS_EQUIV_TRACE));

    hss_system_destroy(s1);
    hss_system_destroy(s2);
    PASS();
}

static void test_stability_preservation(void) {
    TEST("verify_stability_preservation");
    HSS_System *orig = hss_system_create("orig", 1, 2, 0);
    HSS_System *trans = hss_system_create("trans", 1, 2, 0);
    double A[4] = {-1, 0, 0, -2};
    hss_mode_set_dynamics(orig, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(trans, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    double P[4] = {1, 0, 0, 1};
    hss_mode_set_lyapunov(orig, 0, P, 0.5);
    hss_mode_set_lyapunov(trans, 0, P, 0.5);

    assert(hss_verify_stability_preservation(orig, trans, 1e-6));
    assert(!hss_verify_stability_preservation(NULL, trans, 1e-6));

    hss_system_destroy(orig);
    hss_system_destroy(trans);
    PASS();
}

static void test_mode_distinguishability(void) {
    TEST("modes_distinguishable");
    HSS_System *sys = hss_system_create("dist", 2, 2, 0);
    double A1[4] = {-1, 0, 0, -2};
    double A2[4] = {-5, 0, 0, -1};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A1, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A2, NULL, NULL);

    /* Different dynamics should be distinguishable */
    bool result = hss_modes_distinguishable(sys, 0, 1, 0.01);
    assert(result); /* Different A matrices */

    /* Same mode should not be distinguishable */
    assert(!hss_modes_distinguishable(sys, 0, 0, 0.01));

    hss_system_destroy(sys);
    PASS();
}

static void test_abstraction_lift(void) {
    TEST("lift_abstraction");
    HSS_System *sys = hss_system_create("abs", 2, 2, 0);
    double A[4] = {-1, 0, 0, -2};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_System *lifted = hss_lift_abstraction(sys, HSS_LEVEL_TIMED);
    assert(lifted != NULL);
    assert(lifted->state_dim > sys->state_dim); /* Added clock */
    hss_system_destroy(lifted);

    HSS_System *lowered = hss_lower_abstraction(sys, HSS_LEVEL_DISCRETE);
    assert(lowered != NULL);
    hss_system_destroy(lowered);

    hss_system_destroy(sys);
    PASS();
}

static void test_merge_systems(void) {
    TEST("merge_systems");
    HSS_System *s1 = hss_system_create("m1", 2, 2, 0);
    HSS_System *s2 = hss_system_create("m2", 1, 2, 0);
    double A[4] = {-1, 0, 0, -2};
    hss_mode_set_dynamics(s1, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(s1, 1, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(s2, 0, HSS_CLASS_LINEAR, A, NULL, NULL);

    int ident[2] = {0, 0}; /* Identify s1.mode0 with s2.mode0 */
    HSS_System *merged = hss_merge_systems(s1, s2, ident, 1);
    assert(merged != NULL);

    hss_system_destroy(merged);
    hss_system_destroy(s1);
    hss_system_destroy(s2);
    PASS();
}

static void test_compose_feedback(void) {
    TEST("compose (feedback)");
    HSS_System *s1 = hss_system_create("fb1", 1, 1, 0);
    HSS_System *s2 = hss_system_create("fb2", 1, 1, 0);
    double A[1] = {-1};
    hss_mode_set_dynamics(s1, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(s2, 0, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_Composition *comp = hss_compose(HSS_COMP_FEEDBACK, s1, s2);
    assert(comp != NULL);

    hss_composition_free(comp);
    hss_system_destroy(s1);
    hss_system_destroy(s2);
    PASS();
}

int main(void) {
    printf("=== HSS Integrator Tests ===\n\n");
    test_compose_series();
    test_compose_parallel();
    test_compose_feedback();
    test_equivalence();
    test_stability_preservation();
    test_mode_distinguishability();
    test_abstraction_lift();
    test_merge_systems();
    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
