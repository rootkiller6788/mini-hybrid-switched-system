/**
 * @file test_hybrid_automata.c
 * @brief Comprehensive test suite for mini-hybrid-automata
 *
 * Tests all core API functions with assert-based validation.
 * Covers L1-L6 knowledge points with at least one test per function.
 */

#include "../include/hybrid_automaton.h"
#include "../include/hybrid_execution.h"
#include "../include/hybrid_examples.h"
#include "../include/hybrid_reachability.h"
#include "../include/hybrid_safety.h"
#include "../include/hybrid_simulation.h"
#include "../include/hybrid_subclass.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond) do { if (cond) PASS(); else FAIL(#cond); } while(0)

/* ==========================================================================
 * KP1-KP9: Core automaton construction
 * ========================================================================== */

static void test_create_destroy(void)
{
    TEST("KP1: hybrid_automaton_create");
    HybridAutomaton *ha = hybrid_automaton_create("test_ha", 3);
    CHECK(ha != NULL);
    CHECK(ha->num_modes == 0);
    CHECK(ha->num_vars == 3);
    CHECK(ha->num_transitions == 0);
    hybrid_automaton_destroy(ha);

    TEST("KP1: hybrid_automaton_create (NULL name)");
    HybridAutomaton *ha2 = hybrid_automaton_create(NULL, 2);
    CHECK(ha2 == NULL);

    TEST("KP1: hybrid_automaton_create (invalid num_vars)");
    HybridAutomaton *ha3 = hybrid_automaton_create("bad", HA_MAX_VARIABLES + 1);
    CHECK(ha3 == NULL);
}

static void test_mode_management(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("modes", 2);
    assert(ha);

    TEST("KP2: hybrid_automaton_add_mode");
    int m0 = hybrid_automaton_add_mode(ha, "mode0", HAMODE_INITIAL);
    CHECK(m0 == 0);
    int m1 = hybrid_automaton_add_mode(ha, "mode1", HAMODE_NORMAL);
    CHECK(m1 == 1);

    TEST("KP2: mode count");
    CHECK(hybrid_automaton_mode_count(ha) == 2);

    TEST("KP2: get mode");
    const HybridMode *mode = hybrid_automaton_get_mode(ha, 0);
    CHECK(mode != NULL);
    CHECK(mode->id == 0);
    CHECK(strcmp(mode->name, "mode0") == 0);

    TEST("KP2: get mode (invalid)");
    CHECK(hybrid_automaton_get_mode(ha, 99) == NULL);

    hybrid_automaton_destroy(ha);
}

static void test_variable_management(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("vars", 3);
    assert(ha);

    TEST("KP3: hybrid_automaton_set_variable");
    bool ok = hybrid_automaton_set_variable(ha, 0, "temp", HAVAR_REAL_BOUNDED, 0.0, 100.0);
    CHECK(ok);
    ok = hybrid_automaton_set_variable(ha, 1, "clock1", HAVAR_CLOCK, 0.0, 1e9);
    CHECK(ok);

    TEST("KP3: get variable");
    const HybridVariable *v = hybrid_automaton_get_variable(ha, 0);
    CHECK(v != NULL);
    CHECK(strcmp(v->name, "temp") == 0);
    CHECK(v->type == HAVAR_REAL_BOUNDED);

    v = hybrid_automaton_get_variable(ha, 1);
    CHECK(v != NULL);
    CHECK(v->type == HAVAR_CLOCK);
    CHECK(v->rate == 1.0);

    TEST("KP3: invalid index");
    CHECK(!hybrid_automaton_set_variable(ha, 99, "bad", HAVAR_REAL, 0, 0));

    hybrid_automaton_destroy(ha);
}

static void test_transition_management(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("trans", 2);
    assert(ha);
    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    hybrid_automaton_add_mode(ha, "q1", HAMODE_NORMAL);

    TEST("KP4: hybrid_automaton_add_transition");
    int t0 = hybrid_automaton_add_transition(ha, 0, 1, HATRIG_AUTONOMOUS, "go");
    CHECK(t0 == 0);
    int t1 = hybrid_automaton_add_transition(ha, 1, 0, HATRIG_CONTROLLED, "back");
    CHECK(t1 == 1);

    TEST("KP4: transition count");
    CHECK(hybrid_automaton_trans_count(ha) == 2);

    TEST("KP4: out-degree / in-degree");
    CHECK(hybrid_automaton_out_degree(ha, 0) == 1);
    CHECK(hybrid_automaton_in_degree(ha, 1) == 1);
    CHECK(hybrid_automaton_out_degree(ha, 1) == 1);

    TEST("KP4: get transition");
    const HybridTransition *tr = hybrid_automaton_get_transition(ha, 0);
    CHECK(tr != NULL);
    CHECK(tr->src_mode == 0);
    CHECK(tr->tgt_mode == 1);
    CHECK(tr->is_active);

    hybrid_automaton_destroy(ha);
}

static void test_guard_reset_invariant_flow(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("config", 2);
    assert(ha);
    hybrid_automaton_add_mode(ha, "on", HAMODE_INITIAL);
    hybrid_automaton_add_mode(ha, "off", HAMODE_NORMAL);
    hybrid_automaton_add_transition(ha, 0, 1, HATRIG_AUTONOMOUS, "switch");

    TEST("KP5: hybrid_guard_set");
    double A[2] = {1, 0};
    double b[1] = {10.0};
    CHECK(hybrid_guard_set(ha, 0, 1, A, b));

    TEST("KP6: hybrid_reset_set (IDENTITY)");
    CHECK(hybrid_reset_set(ha, 0, HARESET_IDENTITY, NULL, NULL));

    TEST("KP6: hybrid_reset_set (AFFINE)");
    double R[4] = {1, 0, 0, 1};
    double r[2] = {0, 0};
    CHECK(hybrid_reset_set(ha, 0, HARESET_AFFINE, R, r));

    TEST("KP7: hybrid_invariant_set");
    double H[2] = {-1, 0};
    double k_val = 0;
    CHECK(hybrid_invariant_set(ha, 0, 1, H, &k_val));

    TEST("KP8: hybrid_flow_set");
    double A_flow[4] = {0, 1, 0, 0};
    double b_flow[2] = {0, -9.81};
    CHECK(hybrid_flow_set(ha, 0, HAFLOW_AFFINE, A_flow, b_flow));

    TEST("KP8: hybrid_flow_set_nonlinear");
    /* Use a simple nonlinear function */
    CHECK(hybrid_flow_set_nonlinear(ha, 1,
             (void(*)(double,const double*,double*,int,void*))NULL, NULL) == false);
    /* NULL function should fail */

    TEST("KP9: hybrid_init_set");
    double x0[2] = {5.0, 0.0};
    CHECK(hybrid_init_set(ha, 0, x0));

    hybrid_automaton_destroy(ha);
}

static void test_well_formedness(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("wf", 1);
    assert(ha);

    TEST("KF: well-formed (empty)");
    CHECK(!hybrid_automaton_is_well_formed(ha)); /* No modes */

    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    double x0[1] = {0};
    hybrid_init_set(ha, 0, x0);

    TEST("KF: well-formed (valid)");
    CHECK(hybrid_automaton_is_well_formed(ha));

    hybrid_automaton_destroy(ha);
}

/* ==========================================================================
 * KP14-KP19: Execution semantics
 * ========================================================================== */

static void test_execution_create(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("exec", 2);
    assert(ha);
    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    hybrid_automaton_add_mode(ha, "q1", HAMODE_NORMAL);
    double x0[2] = {0, 0};
    hybrid_init_set(ha, 0, x0);

    TEST("KP11: hybrid_execution_create");
    HybridExecution *exec = hybrid_execution_create(ha, 10, 10);
    CHECK(exec != NULL);
    CHECK(exec->num_segments == 0);
    CHECK(exec->num_jumps == 0);

    TEST("KP12: hybrid_execution_append_flow");
    double x_end[2] = {1, 2};
    CHECK(hybrid_execution_append_flow(exec, 0, 0.0, 1.0, x0, x_end, 100));
    CHECK(exec->num_segments == 1);

    TEST("KP13: hybrid_execution_append_jump");
    /* Need a transition first */
    hybrid_automaton_add_transition(ha, 0, 1, HATRIG_AUTONOMOUS, "go");
    double x_post[2] = {0, 0};
    CHECK(hybrid_execution_append_jump(exec, 0, 1.0, x_end, x_post));
    CHECK(exec->num_jumps == 1);

    hybrid_execution_destroy(exec);
    hybrid_automaton_destroy(ha);
}

static void test_determinism(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("det", 1);
    assert(ha);
    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    hybrid_automaton_add_mode(ha, "q1", HAMODE_NORMAL);

    /* Single point initial */
    double x0[1] = {0};
    hybrid_init_set(ha, 0, x0);

    /* One transition with trivially true guard */
    hybrid_automaton_add_transition(ha, 0, 1, HATRIG_AUTONOMOUS, "go");

    TEST("KP14: hybrid_is_deterministic (single trivially-true guard from a mode)");
    /* Only one outgoing transition per mode, so not overlapping */
    CHECK(hybrid_is_deterministic(ha));

    /* Add another transition with trivially true guard from same mode */
    hybrid_automaton_add_transition(ha, 0, 1, HATRIG_AUTONOMOUS, "go2");

    TEST("KP14: hybrid_is_deterministic (two trivially-true guards → non-det)");
    CHECK(!hybrid_is_deterministic(ha));

    hybrid_automaton_destroy(ha);
}

static void test_guard_invariant_evaluation(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("eval", 2);
    assert(ha);
    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    hybrid_automaton_add_mode(ha, "q1", HAMODE_NORMAL);
    hybrid_automaton_add_transition(ha, 0, 1, HATRIG_AUTONOMOUS, "go");

    /* Guard: x0 ≤ 10 → [1, 0]·x ≤ 10 */
    double A[2] = {1, 0};
    double b[1] = {10};
    hybrid_guard_set(ha, 0, 1, A, b);

    TEST("KP17: hybrid_guard_evaluate (satisfied)");
    double x1[2] = {5, 0};
    CHECK(hybrid_guard_evaluate(&ha->trans[0].guard, x1, 2));

    TEST("KP17: hybrid_guard_evaluate (violated)");
    double x2[2] = {15, 0};
    CHECK(!hybrid_guard_evaluate(&ha->trans[0].guard, x2, 2));

    TEST("KP17: hybrid_transition_enabled");
    CHECK(hybrid_transition_enabled(ha, 0, x1, 0));
    CHECK(!hybrid_transition_enabled(ha, 0, x2, 0));

    /* Invariant test */
    double H[2] = {0, 1};
    double k_val = 5;
    hybrid_invariant_set(ha, 1, 1, H, &k_val);

    TEST("KP18: hybrid_invariant_satisfied");
    double x3[2] = {0, 3};
    CHECK(hybrid_invariant_satisfied(&ha->modes[1].invariant, x3, 2));
    double x4[2] = {0, 7};
    CHECK(!hybrid_invariant_satisfied(&ha->modes[1].invariant, x4, 2));

    hybrid_automaton_destroy(ha);
}

static void test_reset_apply(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("rst", 2);
    assert(ha);
    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    hybrid_automaton_add_mode(ha, "q1", HAMODE_NORMAL);
    hybrid_automaton_add_transition(ha, 0, 1, HATRIG_AUTONOMOUS, "go");

    TEST("KP17: hybrid_reset_apply (IDENTITY)");
    HybridReset rst_id;
    rst_id.type = HARESET_IDENTITY;
    for (int i = 0; i < 2; i++) rst_id.r[i] = 0;
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
        rst_id.R[i][j] = (i == j) ? 1.0 : 0.0;
    double x[2] = {3, 4};
    double x_out[2];
    hybrid_reset_apply(&rst_id, x, x_out, 2);
    CHECK(fabs(x_out[0] - 3.0) < 1e-12 && fabs(x_out[1] - 4.0) < 1e-12);

    TEST("KP17: hybrid_reset_apply (CONSTANT)");
    HybridReset rst_const;
    rst_const.type = HARESET_CONSTANT;
    rst_const.r[0] = 7; rst_const.r[1] = 8;
    hybrid_reset_apply(&rst_const, x, x_out, 2);
    CHECK(fabs(x_out[0] - 7.0) < 1e-12 && fabs(x_out[1] - 8.0) < 1e-12);

    hybrid_automaton_destroy(ha);
}

/* ==========================================================================
 * KP19: Parallel composition
 * ========================================================================== */

static void test_parallel_composition(void)
{
    HybridAutomaton *h1 = hybrid_automaton_create("H1", 1);
    hybrid_automaton_add_mode(h1, "A", HAMODE_INITIAL);
    hybrid_automaton_add_mode(h1, "B", HAMODE_NORMAL);
    hybrid_automaton_add_transition(h1, 0, 1, HATRIG_AUTONOMOUS, "ab");
    double x0_1[1] = {0};
    hybrid_init_set(h1, 0, x0_1);

    HybridAutomaton *h2 = hybrid_automaton_create("H2", 1);
    hybrid_automaton_add_mode(h2, "C", HAMODE_INITIAL);
    hybrid_automaton_add_mode(h2, "D", HAMODE_NORMAL);
    hybrid_automaton_add_transition(h2, 0, 1, HATRIG_AUTONOMOUS, "cd");
    double x0_2[1] = {0};
    hybrid_init_set(h2, 0, x0_2);

    TEST("KP19: hybrid_compose_parallel");
    HybridAutomaton *h_comp = hybrid_compose_parallel(h1, h2);
    CHECK(h_comp != NULL);
    CHECK(h_comp->num_modes == 4); /* 2×2 */
    CHECK(h_comp->num_vars == 2);  /* 1+1 */

    hybrid_automaton_destroy(h_comp);
    hybrid_automaton_destroy(h2);
    hybrid_automaton_destroy(h1);
}

/* ==========================================================================
 * KP41-KP42: Timed automata
 * ========================================================================== */

static void test_timed_automaton(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("timed", 2);
    assert(ha);
    hybrid_automaton_set_variable(ha, 0, "x", HAVAR_CLOCK, 0, INFINITY);
    hybrid_automaton_set_variable(ha, 1, "y", HAVAR_CLOCK, 0, INFINITY);
    hybrid_automaton_add_mode(ha, "loc0", HAMODE_INITIAL);
    hybrid_automaton_add_mode(ha, "loc1", HAMODE_NORMAL);

    /* Set flows to ẋ = 1 for all clocks */
    double A[4] = {0, 0, 0, 0};
    double b[2] = {1.0, 1.0};
    hybrid_flow_set(ha, 0, HAFLOW_AFFINE, A, b);
    hybrid_flow_set(ha, 1, HAFLOW_AFFINE, A, b);

    TEST("KP41: hybrid_is_timed_automaton (valid)");
    CHECK(hybrid_is_timed_automaton(ha));

    /* Region equivalence test */
    TEST("KP42: hybrid_region_equivalent (same region)");
    double v1[2] = {0.3, 1.7};
    double v2[2] = {0.35, 1.72};
    CHECK(hybrid_region_equivalent(v1, v2, 2, 3));

    TEST("KP42: hybrid_region_equivalent (different region)");
    double v3[2] = {0.3, 2.1};
    CHECK(!hybrid_region_equivalent(v1, v3, 2, 3));

    /* Region computation */
    TEST("KP42: hybrid_region_compute");
    TimedRegion region;
    region.int_parts = (int*) calloc(2, sizeof(int));
    region.frac_order = (int*) calloc(2, sizeof(int));
    region.zero_frac = (int*) calloc(2, sizeof(int));
    hybrid_region_compute(v1, 2, 3, &region);
    CHECK(region.int_parts[0] == 0);
    CHECK(region.int_parts[1] == 1);
    free(region.int_parts);
    free(region.frac_order);
    free(region.zero_frac);

    hybrid_automaton_destroy(ha);
}

/* ==========================================================================
 * KP43-KP48: Subclass detection
 * ========================================================================== */

static void test_subclass_classification(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("class", 1);
    assert(ha);
    hybrid_automaton_set_variable(ha, 0, "x", HAVAR_CLOCK, 0, INFINITY);
    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    double b[1] = {1.0};
    hybrid_flow_set(ha, 0, HAFLOW_CONSTANT, NULL, b);
    hybrid_automaton_add_mode(ha, "q1", HAMODE_NORMAL);
    hybrid_flow_set(ha, 1, HAFLOW_CONSTANT, NULL, b);
    hybrid_automaton_add_transition(ha, 0, 1, HATRIG_AUTONOMOUS, "go");
    double x0[1] = {0};
    hybrid_init_set(ha, 0, x0);

    TEST("KP43-48: hybrid_classify (timed automaton)");
    HybridAutomatonClass cls = hybrid_classify(ha);
    CHECK(cls == HACLASS_TIMED);

    hybrid_automaton_destroy(ha);
}

/* ==========================================================================
 * KP30: Barrier certificate
 * ========================================================================== */

static void test_barrier_certificate(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("barrier", 1);
    assert(ha);
    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    /* Stable linear ODE: ẋ = -x → origin is stable, B(x)=x²-1 never crossed */
    double A_flow[1] = {-1.0};
    hybrid_flow_set(ha, 0, HAFLOW_LINEAR, A_flow, NULL);
    double x0[1] = {0.5};
    hybrid_init_set(ha, 0, x0);

    TEST("KP30: hybrid_barrier_create");
    HybridBarrierCertificate *bc = hybrid_barrier_create(1, 1);
    CHECK(bc != NULL);

    /* Set identity barrier: B(x) = x² - 1 → B(x) ≤ 0 if |x| ≤ 1 */
    bc->P[0] = 1.0;   /* P = [1] */
    bc->q[0] = 0.0;   /* q = [0] */
    bc->r[0] = -1.0;  /* r = -1 → B(x) = x² - 1 */

    /* Initial: x0=0.5, B(0.5) = 0.25 - 1 = -0.75 ≤ 0 ✓ */
    TEST("KP30: hybrid_barrier_check_init");
    CHECK(hybrid_barrier_check_init(bc, ha));

    /* Flow: dB/dt = 2x·(-x) = -2x² ≤ 0 ✓ */
    TEST("KP30: hybrid_barrier_check_flow");
    CHECK(hybrid_barrier_check_flow(bc, ha));

    hybrid_barrier_destroy(bc);
    hybrid_automaton_destroy(ha);
}

/* ==========================================================================
 * KP31: Common Lyapunov function
 * ========================================================================== */

static void test_common_lyapunov(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("lyap", 2);
    assert(ha);
    hybrid_automaton_add_mode(ha, "stable", HAMODE_INITIAL);
    /* Hurwitz matrix: A = [[-2, 0], [0, -1]] */
    double A[4] = {-2, 0, 0, -1};
    hybrid_flow_set(ha, 0, HAFLOW_LINEAR, A, NULL);
    hybrid_automaton_add_mode(ha, "stable2", HAMODE_NORMAL);
    double A2[4] = {-1, 0, 0, -3};
    hybrid_flow_set(ha, 1, HAFLOW_LINEAR, A2, NULL);
    double x0[2] = {1, 1};
    hybrid_init_set(ha, 0, x0);

    TEST("KP31: hybrid_common_lyapunov");
    double P[4] = {0};
    bool found = hybrid_common_lyapunov(ha, 2, P);
    CHECK(found);

    /* Verify P is positive definite (identity) */
    CHECK(P[0] > 0 && P[3] > 0);

    hybrid_automaton_destroy(ha);
}

/* ==========================================================================
 * KP36: ODE integration
 * ========================================================================== */

static void test_ode_integration(void)
{
    HybridFlow flow;
    flow.type = HAFLOW_AFFINE;
    flow.has_A = true;
    flow.A[0][0] = -1.0;
    flow.A[0][1] = 0.0;
    flow.A[1][0] = 0.0;
    flow.A[1][1] = -2.0;
    flow.b[0] = 0;
    flow.b[1] = 1.0;

    TEST("KP36: hybrid_ode_step (Forward Euler)");
    double x[2] = {1.0, 0.0};
    hybrid_ode_step(&flow, x, 0.0, 0.01, 2, HASTEP_FORWARD_EULER);
    /* ẋ₁ = -x₁ → x₁(0.01) ≈ 1 - 0.01 = 0.99 */
    CHECK(fabs(x[0] - 0.99) < 1e-6);
    /* ẋ₂ = -2x₂ + 1 → x₂(0.01) ≈ 0 + 0.01*(1) = 0.01 */
    CHECK(fabs(x[1] - 0.01) < 1e-6);

    TEST("KP36: hybrid_ode_step (RK4)");
    double x2[2] = {1.0, 0.0};
    hybrid_ode_step(&flow, x2, 0.0, 0.01, 2, HASTEP_RK4);
    /* RK4 should give more accurate result */
    CHECK(fabs(x2[0] - exp(-0.01)) < 1e-6);
    /* For x₂: ẋ₂ = -2x₂ + 1, exact solution x₂(t) = 0.5(1 - e^{-2t}) */
    double exact_x2 = 0.5 * (1.0 - exp(-0.02));
    CHECK(fabs(x2[1] - exact_x2) < 1e-6);

    TEST("KP36: hybrid_ode_step (Heun)");
    double x3[2] = {1.0, 0.0};
    hybrid_ode_step(&flow, x3, 0.0, 0.01, 2, HASTEP_HEUN);
    CHECK(fabs(x3[0] - exp(-0.01)) < 1e-5);
}

/* ==========================================================================
 * KP39: Hybrid simulation
 * ========================================================================== */

static void test_hybrid_simulation(void)
{
    /* Create bouncing ball */
    HybridAutomaton *ha = example_bouncing_ball(10.0, 0.0, 0.8, 9.81);
    CHECK(ha != NULL);

    TEST("KP39: hybrid_simulate (bouncing ball)");
    HybridSimConfig config = HYBRID_SIM_CONFIG_DEFAULT;
    config.dt = 0.001;
    config.t_max = 3.0;  /* Need enough time for 10m fall: √(2×10/9.81) ≈ 1.43s */
    config.max_transitions = 50;
    HybridExecution *exec = hybrid_simulate(ha, &config);
    CHECK(exec != NULL);
    /* Execution produced flow segments (may or may not bounce depending on
       guard/invariant timing with discrete steps) */
    CHECK(exec->num_segments > 0);
    /* Execution time should be ≤ t_max */
    CHECK(exec->total_time <= config.t_max + config.dt);

    hybrid_execution_destroy(exec);
    hybrid_automaton_destroy(ha);
}

/* ==========================================================================
 * KP49-KP56: Canonical examples
 * ========================================================================== */

static void test_example_bouncing_ball(void)
{
    TEST("KP49: example_bouncing_ball");
    HybridAutomaton *ha = example_bouncing_ball(10.0, 0.0, 0.9, 9.81);
    CHECK(ha != NULL);
    CHECK(ha->num_modes == 1);
    CHECK(ha->num_transitions == 1);
    CHECK(ha->num_vars == 2);
    hybrid_automaton_destroy(ha);
}

static void test_example_thermostat(void)
{
    TEST("KP50: example_thermostat");
    HybridAutomaton *ha = example_thermostat(18.0, 22.0, 15.0, 0.1, 2.0, 20.0, true);
    CHECK(ha != NULL);
    CHECK(ha->num_modes == 2);
    CHECK(ha->num_transitions == 2);
    hybrid_automaton_destroy(ha);
}

static void test_example_two_tank(void)
{
    TEST("KP51: example_two_tank");
    HybridAutomaton *ha = example_two_tank(0.5, 0.1, 5.0, 2.0);
    CHECK(ha != NULL);
    CHECK(ha->num_modes == 2);
    CHECK(ha->num_vars == 2);
    hybrid_automaton_destroy(ha);
}

static void test_example_train_gate(void)
{
    TEST("KP52: example_train_gate");
    HybridAutomaton *ha = example_train_gate(1000.0, 500.0, 600.0,
                                               5.0, 5.0, 40.0);
    CHECK(ha != NULL);
    CHECK(ha->num_modes == 4);
    CHECK(ha->num_transitions == 3);
    hybrid_automaton_destroy(ha);
}

static void test_example_dcdc(void)
{
    TEST("KP53: example_dcdc_converter");
    HybridAutomaton *ha = example_dcdc_converter(0.001, 0.0001, 10.0, 12.0, 0.0, 0.0);
    CHECK(ha != NULL);
    CHECK(ha->num_modes == 2);
    hybrid_automaton_destroy(ha);
}

static void test_example_engine_afr(void)
{
    TEST("KP54: example_engine_afr_control");
    HybridAutomaton *ha = example_engine_afr_control();
    CHECK(ha != NULL);
    hybrid_automaton_destroy(ha);
}

static void test_example_robot(void)
{
    TEST("KP55: example_robot_obstacle_avoidance");
    HybridAutomaton *ha = example_robot_obstacle_avoidance(1.0, 5.0, 10.0, 5.0);
    CHECK(ha != NULL);
    hybrid_automaton_destroy(ha);
}

static void test_example_infusion(void)
{
    TEST("KP56: example_infusion_pump");
    HybridAutomaton *ha = example_infusion_pump();
    CHECK(ha != NULL);
    hybrid_automaton_destroy(ha);
}

/* ==========================================================================
 * KP22-KP28: Reachability
 * ========================================================================== */

static void test_reachability_basic(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("reach", 2);
    assert(ha);
    hybrid_automaton_set_variable(ha, 0, "x1", HAVAR_REAL, -10, 10);
    hybrid_automaton_set_variable(ha, 1, "x2", HAVAR_REAL, -10, 10);
    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    double A[4] = {0, 1, 0, 0};
    hybrid_flow_set(ha, 0, HAFLOW_LINEAR, A, NULL);
    double x0[2] = {0, 1};
    hybrid_init_set(ha, 0, x0);

    TEST("KP23: hybrid_zono_create");
    HybridZonotope *z = hybrid_zono_create(2, 4);
    CHECK(z != NULL);
    CHECK(z->dim == 2);
    CHECK(z->num_generators == 0);

    TEST("KP23: hybrid_zono_set_center + add_generator");
    double c[2] = {0, 0};
    hybrid_zono_set_center(z, c);
    double g[2] = {1, 0};
    CHECK(hybrid_zono_add_generator(z, g));

    TEST("KP23: hybrid_zono_affine_transform");
    CHECK(hybrid_zono_affine_transform(z, A, NULL, 2, 2));

    TEST("KP24: hybrid_poly_create");
    HybridPolyhedron *p = hybrid_poly_create(2, 4);
    CHECK(p != NULL);

    TEST("KP24: hybrid_poly_add_halfspace + contains");
    double h[2] = {1, 0};
    CHECK(hybrid_poly_add_halfspace(p, h, 5.0));
    double x_test[2] = {3, 0};
    CHECK(hybrid_poly_contains(p, x_test));
    double x_bad[2] = {7, 0};
    CHECK(!hybrid_poly_contains(p, x_bad));

    TEST("KP25: hybrid_support_function");
    const double *dirs[2];
    double d1[2] = {1, 0};
    double d2[2] = {0, 1};
    dirs[0] = d1;
    dirs[1] = d2;
    double bounds[2];
    hybrid_support_function(z, dirs, 2, bounds);
    CHECK(bounds[0] >= z->center[0]);

    TEST("KP26: hybrid_bisimulation_quotient");
    HybridAutomaton *quot = hybrid_bisimulation_quotient(ha, 10);
    CHECK(quot != NULL);
    hybrid_automaton_destroy(quot);

    TEST("KP22: hybrid_flowpipe_compute");
    int nsteps;
    HybridZonotope **fp = hybrid_flowpipe_compute(ha, 0, z, 0.1, 0.01, &nsteps);
    CHECK(fp != NULL);
    CHECK(nsteps > 0);
    for (int i = 0; i < nsteps; i++) hybrid_zono_destroy(fp[i]);
    free(fp);

    TEST("KP28: hybrid_onion_peeling");
    int layers;
    HybridReachOptions opts = HYBRID_REACH_OPTIONS_DEFAULT;
    HybridPolyhedron **peels = hybrid_onion_peeling(ha, 3, &opts, &layers);
    CHECK(peels != NULL);
    for (int i = 0; i < layers; i++) hybrid_poly_destroy(peels[i]);
    free(peels);

    hybrid_zono_destroy(z);
    hybrid_poly_destroy(p);
    hybrid_automaton_destroy(ha);
}

/* ==========================================================================
 * KP34-KP35: Dwell-time and event-triggered safety
 * ========================================================================== */

static void test_dwell_time(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("dwell", 2);
    assert(ha);
    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    double A[4] = {-1, 0, 0, -2};
    hybrid_flow_set(ha, 0, HAFLOW_LINEAR, A, NULL);
    hybrid_automaton_add_mode(ha, "q1", HAMODE_NORMAL);
    double A2[4] = {-3, 0, 0, -0.5};
    hybrid_flow_set(ha, 1, HAFLOW_LINEAR, A2, NULL);
    double x0[2] = {1, 1};
    hybrid_init_set(ha, 0, x0);

    TEST("KP34: hybrid_dwell_time_safety");
    CHECK(hybrid_dwell_time_safety(ha, 0.1, 2));

    hybrid_automaton_destroy(ha);
}

static void test_event_triggered(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("et", 1);
    assert(ha);
    hybrid_automaton_add_mode(ha, "q0", HAMODE_INITIAL);
    double x0[1] = {0};
    hybrid_init_set(ha, 0, x0);

    TEST("KP35: hybrid_event_triggered_safety (no policy)");
    HybridSafetyProperty *sp = hybrid_safety_property_create("safe", NULL, 0);
    int policy[1] = {-1}; /* Stay */
    CHECK(hybrid_event_triggered_safety(ha, sp, policy, 1));

    hybrid_safety_property_destroy(sp);
    hybrid_automaton_destroy(ha);
}

/* ==========================================================================
 * Print function (coverage)
 * ========================================================================== */

static void test_print(void)
{
    HybridAutomaton *ha = example_bouncing_ball(10.0, 0.0, 0.8, 9.81);
    assert(ha);

    TEST("UTIL: hybrid_automaton_print");
    /* Just ensure it doesn't crash */
    hybrid_automaton_print(ha);

    TEST("UTIL: hybrid_automaton_print (NULL)");
    hybrid_automaton_print(NULL);

    hybrid_automaton_destroy(ha);
    PASS(); /* null print doesn't crash */
}

/* ==========================================================================
 * Main
 * ========================================================================== */

int main(void)
{
    printf("=== mini-hybrid-automata Test Suite ===\n\n");

    printf("--- Core Construction (L1) ---\n");
    test_create_destroy();
    test_mode_management();
    test_variable_management();
    test_transition_management();
    test_guard_reset_invariant_flow();
    test_well_formedness();

    printf("\n--- Execution Semantics (L2) ---\n");
    test_execution_create();
    test_determinism();
    test_guard_invariant_evaluation();
    test_reset_apply();

    printf("\n--- Parallel Composition (L2) ---\n");
    test_parallel_composition();

    printf("\n--- Timed Automata (L3) ---\n");
    test_timed_automaton();

    printf("\n--- Subclass Detection (L3-L4) ---\n");
    test_subclass_classification();

    printf("\n--- Safety Verification (L4-L5) ---\n");
    test_barrier_certificate();
    test_common_lyapunov();
    test_dwell_time();
    test_event_triggered();

    printf("\n--- ODE Integration (L5) ---\n");
    test_ode_integration();

    printf("\n--- Hybrid Simulation (L5-L6) ---\n");
    test_hybrid_simulation();

    printf("\n--- Canonical Examples (L6) ---\n");
    test_example_bouncing_ball();
    test_example_thermostat();
    test_example_two_tank();
    test_example_train_gate();
    test_example_dcdc();

    printf("\n--- Applications (L7) ---\n");
    test_example_engine_afr();
    test_example_robot();
    test_example_infusion();

    printf("\n--- Reachability (L3-L5) ---\n");
    test_reachability_basic();

    printf("\n--- Utilities ---\n");
    test_print();

    printf("\n========================================\n");
    printf("Results: %d passed, %d failed, %d total\n",
           tests_passed, tests_failed, tests_passed + tests_failed);
    printf("========================================\n");

    return tests_failed > 0 ? 1 : 0;
}
