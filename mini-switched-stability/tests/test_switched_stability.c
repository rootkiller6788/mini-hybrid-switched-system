#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "switched_types.h"
#include "switched_stability.h"
#include "switched_lyapunov.h"
#include "switched_dwell_time.h"
#include "switched_applications.h"

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST %s... ", #name); } while(0)
#define PASS() do { tests_passed++; printf("PASSED\n"); } while(0)
#define FAIL(msg) do { printf("FAILED: %s\n", msg); } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)
#define CHECK_EQ(a, b, msg) do { if ((a) != (b)) { printf("FAILED: %s (got %d, expected %d)\n", msg, a, b); return; } } while(0)
#define CHECK_CLOSE(a, b, tol, msg) do { if (fabs((a)-(b)) > (tol)) { printf("FAILED: %s (got %.6f, expected %.6f)\n", msg, a, b); return; } } while(0)

/* ===================================================================
 * Test 1: Vector Operations
 * =================================================================== */
static void test_vector_ops(void) {
    TEST(vector_create);
    SwitchedVector v = sv_create(5);
    CHECK(v.n == 5, "vector dim");
    CHECK(v.data != NULL, "vector data allocation");
    sv_free(&v);
    PASS();
}

static void test_vector_set_get(void) {
    TEST(vector_set_get);
    SwitchedVector v = sv_create(3);
    sv_set(&v, 0, 1.0);
    sv_set(&v, 1, 2.0);
    sv_set(&v, 2, 3.0);
    CHECK_CLOSE(sv_get(&v, 0), 1.0, 1e-12, "v[0]");
    CHECK_CLOSE(sv_get(&v, 1), 2.0, 1e-12, "v[1]");
    CHECK_CLOSE(sv_get(&v, 2), 3.0, 1e-12, "v[2]");
    sv_free(&v);
    PASS();
}

static void test_vector_norm(void) {
    TEST(vector_norm);
    SwitchedVector v = sv_create(3);
    sv_set(&v, 0, 3.0);
    sv_set(&v, 1, 4.0);
    sv_set(&v, 2, 0.0);
    CHECK_CLOSE(sv_norm(&v), 5.0, 1e-12, "norm of [3,4,0]");
    sv_free(&v);
    PASS();
}

static void test_vector_dot(void) {
    TEST(vector_dot);
    SwitchedVector a = sv_create(3), b = sv_create(3);
    sv_set(&a, 0, 1.0); sv_set(&a, 1, 2.0); sv_set(&a, 2, 3.0);
    sv_set(&b, 0, 4.0); sv_set(&b, 1, 5.0); sv_set(&b, 2, 6.0);
    CHECK_CLOSE(sv_dot(&a, &b), 32.0, 1e-12, "dot [1,2,3].[4,5,6]");
    sv_free(&a); sv_free(&b);
    PASS();
}

/* ===================================================================
 * Test 2: Matrix Operations
 * =================================================================== */
static void test_matrix_create(void) {
    TEST(matrix_create);
    SwitchedMatrix M = sm_create(3, 3);
    CHECK(M.rows == 3, "matrix rows");
    CHECK(M.cols == 3, "matrix cols");
    CHECK(M.data != NULL, "matrix data");
    sm_free(&M);
    PASS();
}

static void test_matrix_identity(void) {
    TEST(matrix_identity);
    SwitchedMatrix M = sm_create(3, 3);
    sm_identity(&M);
    CHECK_CLOSE(sm_get(&M, 0, 0), 1.0, 1e-12, "M[0,0]");
    CHECK_CLOSE(sm_get(&M, 0, 1), 0.0, 1e-12, "M[0,1]");
    CHECK_CLOSE(sm_get(&M, 1, 1), 1.0, 1e-12, "M[1,1]");
    CHECK_CLOSE(sm_get(&M, 2, 2), 1.0, 1e-12, "M[2,2]");
    sm_free(&M);
    PASS();
}

static void test_matrix_mul(void) {
    TEST(matrix_mul);
    SwitchedMatrix A = sm_create(2, 2), B = sm_create(2, 2), C = sm_create(2, 2);
    sm_set(&A, 0, 0, 1.0); sm_set(&A, 0, 1, 2.0);
    sm_set(&A, 1, 0, 3.0); sm_set(&A, 1, 1, 4.0);
    sm_set(&B, 0, 0, 2.0); sm_set(&B, 0, 1, 0.0);
    sm_set(&B, 1, 0, 1.0); sm_set(&B, 1, 1, 2.0);
    sm_mul(&C, &A, &B);
    /* A*B = [[1*2+2*1, 1*0+2*2], [3*2+4*1, 3*0+4*2]] = [[4,4],[10,8]] */
    CHECK_CLOSE(sm_get(&C, 0, 0), 4.0, 1e-12, "C[0,0]");
    CHECK_CLOSE(sm_get(&C, 0, 1), 4.0, 1e-12, "C[0,1]");
    CHECK_CLOSE(sm_get(&C, 1, 0), 10.0, 1e-12, "C[1,0]");
    CHECK_CLOSE(sm_get(&C, 1, 1), 8.0, 1e-12, "C[1,1]");
    sm_free(&A); sm_free(&B); sm_free(&C);
    PASS();
}

static void test_matrix_det(void) {
    TEST(matrix_det2x2);
    SwitchedMatrix M = sm_create(2, 2);
    sm_set(&M, 0, 0, 1.0); sm_set(&M, 0, 1, 2.0);
    sm_set(&M, 1, 0, 3.0); sm_set(&M, 1, 1, 4.0);
    CHECK_CLOSE(sm_det_2x2(&M), -2.0, 1e-12, "det 2x2");
    sm_free(&M);
    PASS();
}

static void test_matrix_trace(void) {
    TEST(matrix_trace);
    SwitchedMatrix M = sm_create(3, 3);
    sm_set(&M, 0, 0, 1.0); sm_set(&M, 1, 1, 5.0); sm_set(&M, 2, 2, -3.0);
    CHECK_CLOSE(sm_trace(&M), 3.0, 1e-12, "trace");
    sm_free(&M);
    PASS();
}

/* ===================================================================
 * Test 3: Switching Signal
 * =================================================================== */
static void test_switching_signal(void) {
    TEST(switching_signal);
    SwitchingSignal *sig = ssig_create(10);
    CHECK(sig != NULL, "signal creation");
    CHECK(sig->n_switches == 0, "initial switches");

    ssig_record_switch(sig, 1, 1.0);
    CHECK(sig->n_switches == 1, "one switch");
    CHECK(sig->mode_sequence[1] == 1, "mode after switch");

    int mode = ssig_active_mode_at(sig, 0.5);
    CHECK(mode == 0, "mode at t=0.5");

    mode = ssig_active_mode_at(sig, 2.0);
    CHECK(mode == 1, "mode at t=2.0");

    ssig_free(sig);
    PASS();
}

/* ===================================================================
 * Test 4: Subsystem Creation
 * =================================================================== */
static void test_subsystem(void) {
    TEST(subsystem_create);
    SwitchedSubsystem *sub = ssub_create(0, 2);
    CHECK(sub != NULL, "subsystem creation");
    CHECK(sub->mode_id == 0, "mode_id");

    /* Set Hurwitz matrix A = [[-1, 0], [0, -2]] */
    sm_set(&sub->A, 0, 0, -1.0);
    sm_set(&sub->A, 1, 1, -2.0);
    ssub_set_hurwitz_check(sub);
    CHECK(sub->is_hurwitz == true, "is_hurwitz");

    ssub_free(sub);
    PASS();
}

/* ===================================================================
 * Test 5: Switched System
 * =================================================================== */
static void test_switched_system_create(void) {
    TEST(switched_system_create);
    SwitchedSystem *sys = sss_create("test_sys", 2, 2);
    CHECK(sys != NULL, "system creation");
    CHECK(sys->state_dim == 2, "state_dim");
    CHECK(sys->n_modes == 2, "n_modes");
    sss_free(sys);
    PASS();
}

/* ===================================================================
 * Test 6: Lyapunov Functions
 * =================================================================== */
static void test_lyap_solve_2x2(void) {
    TEST(lyap_solve_2x2);
    /* A = [[-1, 0], [0, -2]] (Hurwitz) */
    SwitchedMatrix A = sm_create(2, 2);
    sm_set(&A, 0, 0, -1.0);
    sm_set(&A, 1, 1, -2.0);
    SwitchedMatrix P = sm_create(2, 2);

    bool ok = sss_lyap_2x2(&A, &P);
    CHECK(ok == true, "lyap solve returned true");

    /* Check that P is positive definite */
    CHECK(sm_is_positive_definite(&P) == true, "P positive definite");

    sm_free(&A); sm_free(&P);
    PASS();
}

static void test_lyap_eval(void) {
    TEST(lyap_eval);
    SwitchedMatrix P = sm_create(2, 2);
    sm_identity(&P);
    SwitchedVector x = sv_create(2);
    sv_set(&x, 0, 3.0);
    sv_set(&x, 1, 4.0);
    double V = sss_lyap_eval(&P, &x);
    CHECK_CLOSE(V, 25.0, 1e-12, "V = x^T I x = 25");
    sm_free(&P); sv_free(&x);
    PASS();
}

/* ===================================================================
 * Test 7: Dwell-Time Analysis
 * =================================================================== */
static void test_dwell_time_computation(void) {
    TEST(dwell_time_computation);
    double lambda_0 = 1.0;
    double mu = exp(2.0); /* mu = e^2 */
    double tau_d = sdt_compute_min_dwell(lambda_0, mu);
    CHECK_CLOSE(tau_d, 2.0, 1e-10, "tau_d = ln(e^2)/1 = 2");

    double tau_a = sdt_compute_avg_dwell(lambda_0, mu);
    CHECK_CLOSE(tau_a, 2.0, 1e-10, "tau_a* = ln(e^2)/1 = 2");
    PASS();
}

static void test_dwell_check(void) {
    TEST(dwell_check);
    double times[] = {0.0, 2.0, 4.0, 6.0};
    CHECK(sdt_check_dwell(times, 3, 1.0) == true, "dwell satisfied");
    CHECK(sdt_check_dwell(times, 3, 3.0) == false, "dwell violated");
    PASS();
}

/* ===================================================================
 * Test 8: Eigenvalue Computation
 * =================================================================== */
static void test_eigenvalues_2x2(void) {
    TEST(eigenvalues_2x2);
    SwitchedMatrix A = sm_create(2, 2);
    /* A = [[0, 1], [-1, 0]] has eigenvalues +/- i */
    sm_set(&A, 0, 1, 1.0);
    sm_set(&A, 1, 0, -1.0);

    QRWorkspace *qr = qr_create(2);
    EigenvalueResult eig[2];
    int n = qr_eigenvalues(qr, &A, eig);
    CHECK(n == 2, "got 2 eigenvalues");

    /* Mag = 1 for pure imaginary */
    CHECK_CLOSE(eig[0].magnitude, 1.0, 1e-10, "magnitude of +/- i");
    CHECK_CLOSE(eig[1].magnitude, 1.0, 1e-10, "magnitude of -/+ i");

    qr_free(qr); sm_free(&A);
    PASS();
}

static void test_hurwitz_check(void) {
    TEST(hurwitz_check);
    SwitchedMatrix A = sm_create(2, 2);
    sm_set(&A, 0, 0, -1.0);
    sm_set(&A, 1, 1, -2.0);
    CHECK(is_hurwitz_matrix(&A) == true, "diag([-1,-2]) is Hurwitz");

    sm_set(&A, 0, 0, 1.0);
    CHECK(is_hurwitz_matrix(&A) == false, "diag([1,-2]) is not Hurwitz");
    sm_free(&A);
    PASS();
}

/* ===================================================================
 * Test 9: Lie Bracket
 * =================================================================== */
static void test_lie_bracket(void) {
    TEST(lie_bracket);
    SwitchedMatrix A = sm_create(2, 2), B = sm_create(2, 2);
    sm_set(&A, 0, 1, 1.0);
    sm_set(&B, 1, 0, 1.0);

    SwitchedMatrix bracket = sm_create(2, 2);
    sm_commutator(&bracket, &A, &B);

    /* [A,B] = [[0,1],[0,0]] * [[0,0],[1,0]] - [[0,0],[1,0]] * [[0,1],[0,0]]
     *       = [[1,0],[0,0]] - [[0,0],[0,1]] = [[1,0],[0,-1]] */
    CHECK_CLOSE(sm_get(&bracket, 0, 0), 1.0, 1e-12, "[A,B][0,0]");
    CHECK_CLOSE(sm_get(&bracket, 1, 1), -1.0, 1e-12, "[A,B][1,1]");

    bool commute = lie_check_commute(&A, &B);
    CHECK(commute == false, "A and B do not commute");

    sm_free(&A); sm_free(&B); sm_free(&bracket);
    PASS();
}

/* ===================================================================
 * Test 10: Lie-Algebraic Condition
 * =================================================================== */
static void test_lie_solvable(void) {
    TEST(lie_solvable);
    /* Two commuting matrices -> solvable Lie algebra */
    SwitchedMatrix *A_arr[2];
    SwitchedMatrix A0 = sm_create(2, 2), A1 = sm_create(2, 2);
    sm_set(&A0, 0, 0, -1.0);
    sm_set(&A1, 0, 0, -2.0); /* Diagonal, commute */
    A_arr[0] = &A0; A_arr[1] = &A1;

    LieAlgebraCondition *la = lie_check_create(2, 2);
    bool solv = lie_condition_solvable(la, A_arr, 2, 2);
    CHECK(solv == true, "commuting -> solvable");
    CHECK(la->pair_commute == true, "pairwise commute");
    CHECK(la->is_nilpotent == true, "nilpotent");

    lie_check_free(la);
    sm_free(&A0); sm_free(&A1);
    PASS();
}

/* ===================================================================
 * Test 11: CLF Verification
 * =================================================================== */
static void test_clf_verify(void) {
    TEST(clf_verify);
    /* Create two stable subsystems sharing the same CLF P = I */
    SwitchedSubsystem *modes[2];
    SwitchedMatrix A0 = sm_create(2, 2), A1 = sm_create(2, 2);
    sm_set(&A0, 0, 0, -1.0); sm_set(&A0, 1, 1, -2.0);
    sm_set(&A1, 0, 0, -3.0); sm_set(&A1, 1, 1, -4.0);

    SwitchedSubsystem *sub0 = ssub_create(0, 2);
    SwitchedSubsystem *sub1 = ssub_create(1, 2);
    sm_copy(&sub0->A, &A0);
    sm_copy(&sub1->A, &A1);
    modes[0] = sub0; modes[1] = sub1;

    SwitchedMatrix P = sm_create(2, 2);
    sm_identity(&P);
    bool ok = sss_clf_verify(&P, modes, 2);
    /* P = I should be a valid CLF since A_i are diagonal with negative entries */
    CHECK(ok == true, "P=I is valid CLF");

    ssub_free(sub0); ssub_free(sub1);
    sm_free(&P); sm_free(&A0); sm_free(&A1);
    PASS();
}

/* ===================================================================
 * Test 12: MLF Computation
 * =================================================================== */
static void test_mlf_compute(void) {
    TEST(mlf_compute);
    SwitchedSystem *sys = sss_create("mlf_test", 2, 2);

    SwitchedMatrix A0 = sm_create(2, 2), A1 = sm_create(2, 2);
    sm_set(&A0, 0, 0, -1.0); sm_set(&A0, 1, 1, -2.0);
    sm_set(&A1, 0, 0, -3.0); sm_set(&A1, 1, 1, -4.0);
    sss_add_subsystem(sys, 0, &A0);
    sss_add_subsystem(sys, 1, &A1);
    sm_free(&A0); sm_free(&A1);

    MultipleLyapunovFunctions mlf;
    mlf.P = (SwitchedMatrix *)malloc(2 * sizeof(SwitchedMatrix));
    mlf.min_eig = (double *)calloc(2, sizeof(double));
    mlf.max_eig = (double *)calloc(2, sizeof(double));
    mlf.is_valid = (bool *)calloc(2, sizeof(bool));
    mlf.decay_rates = (double *)calloc(2, sizeof(double));
    mlf.n_modes = 2;
    mlf.mu = 1.0;
    mlf.mlf_condition = false;
    mlf.P[0] = sm_create(2, 2); mlf.P[1] = sm_create(2, 2);

    sss_compute_mlf(sys, &mlf);
    CHECK(mlf.is_valid[0] == true, "MLF[0] valid");
    CHECK(mlf.is_valid[1] == true, "MLF[1] valid");

    double mu = sss_compute_mu(&mlf);
    CHECK(mu >= 1.0, "mu >= 1.0");

    for (int i = 0; i < 2; i++) sm_free(&mlf.P[i]);
    free(mlf.P); free(mlf.min_eig); free(mlf.max_eig);
    free(mlf.is_valid); free(mlf.decay_rates);
    sss_free(sys);
    PASS();
}

/* ===================================================================
 * Test 13: DC-DC Converter
 * =================================================================== */
static void test_dcdc_converter(void) {
    TEST(dcdc_converter);
    DCDCConverter *conv = dcdc_create(12.0, 24.0, 100e-6, 220e-6, 10.0, 100e3);
    CHECK(conv != NULL, "converter creation");
    CHECK_CLOSE(conv->duty_cycle, 0.5, 1e-6, "duty cycle = 1 - 12/24");

    SwitchedMatrix A_on = sm_create(2, 2), A_off = sm_create(2, 2);
    dcdc_get_matrices(conv, &A_on, &A_off);

    /* Check ON mode: A_on[0,0] = 0 */
    CHECK_CLOSE(sm_get(&A_on, 0, 0), 0.0, 1e-12, "A_on[0,0]");

    /* Check OFF mode: A_off[0,1] = -1/L */
    CHECK_CLOSE(sm_get(&A_off, 0, 1), -1.0 / 100e-6, 1e-6, "A_off[0,1] = -1/L");

    double iL_ss, vC_ss;
    dcdc_equilibrium(conv, 0.5, &iL_ss, &vC_ss);
    CHECK_CLOSE(vC_ss, 24.0, 0.1, "Vout for D=0.5, Vin=12");

    sm_free(&A_on); sm_free(&A_off);
    dcdc_free(conv);
    PASS();
}

/* ===================================================================
 * Test 14: Thermostat System
 * =================================================================== */
static void test_thermostat(void) {
    TEST(thermostat);
    ThermostatSystem *t = thermo_create(20.0, 1.0, 0.5, 0.3, 0.01);
    CHECK(t != NULL, "thermostat creation");

    /* Initial temp = setpoint, should be OFF */
    CHECK(t->current_mode == 0, "initial mode OFF");

    /* Lower temp below deadband -> HEATING */
    t->temp = 18.5;
    thermo_step(t, 1.0);
    CHECK(t->current_mode == 1, "mode switch to HEATING");

    thermo_free(t);
    PASS();
}

/* ===================================================================
 * Test 15: Vehicle Spacing
 * =================================================================== */
static void test_vehicle_spacing(void) {
    TEST(vehicle_spacing);
    VehicleSpacingControl *vsc = vsc_create(30.0, 25.0, 50.0, 30.0);
    CHECK(vsc != NULL, "vsc creation");

    /* Gap > safe_gap -> CRUISE */
    int mode = vsc_determine_mode(vsc);
    CHECK(mode == 1, "gap too large -> ACCEL");

    vsc_free(vsc);
    PASS();
}

/* ===================================================================
 * Test 16: Networked Control
 * =================================================================== */
static void test_networked_control(void) {
    TEST(networked_control);
    double gain = -0.5;
    NetworkedControlDropout *ncs = ncs_create(1, &gain, 0.1, 3.0);
    CHECK(ncs != NULL, "ncs creation");

    ncs_step(ncs, 0.1);
    /* After one step, state should have been updated */
    CHECK(ncs->state != NULL, "state exists");

    ncs_free(ncs);
    PASS();
}

/* ===================================================================
 * Test 17: Simulation
 * =================================================================== */
static void test_simulation(void) {
    TEST(simulation);
    SwitchedSystem *sys = sss_create("sim_test", 2, 2);

    SwitchedMatrix A0 = sm_create(2, 2), A1 = sm_create(2, 2);
    sm_set(&A0, 0, 0, -1.0); sm_set(&A0, 1, 1, -2.0);
    sm_set(&A1, 0, 0, -3.0); sm_set(&A1, 1, 1, -4.0);
    sss_add_subsystem(sys, 0, &A0);
    sss_add_subsystem(sys, 1, &A1);
    sm_free(&A0); sm_free(&A1);

    SwitchedVector x0 = sv_create(2);
    sv_set(&x0, 0, 1.0); sv_set(&x0, 1, 1.0);
    sss_set_initial_state(sys, &x0);

    int modes[] = {0, 1};
    double durs[] = {1.0, 1.0};
    SwitchSequence seq = {modes, durs, 2, true};

    SolverConfig cfg = {0.01, 1.0, 1000, 1e-8, true};
    sss_simulate(sys, &cfg, &seq);

    /* State should have converged towards origin */
    double norm = sv_norm(&sys->state);
    CHECK(norm < 10.0, "state bounded");

    sv_free(&x0);
    sss_free(sys);
    PASS();
}

/* ===================================================================
 * Main
 * =================================================================== */
int main(void) {
    printf("=== Switched Stability Test Suite ===\n\n");

    test_vector_ops();
    test_vector_set_get();
    test_vector_norm();
    test_vector_dot();
    test_matrix_create();
    test_matrix_identity();
    test_matrix_mul();
    test_matrix_det();
    test_matrix_trace();
    test_switching_signal();
    test_subsystem();
    test_switched_system_create();
    test_lyap_solve_2x2();
    test_lyap_eval();
    test_dwell_time_computation();
    test_dwell_check();
    test_eigenvalues_2x2();
    test_hurwitz_check();
    test_lie_bracket();
    test_lie_solvable();
    test_clf_verify();
    test_mlf_compute();
    test_dcdc_converter();
    test_thermostat();
    test_vehicle_spacing();
    test_networked_control();
    test_simulation();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}