/**
 * @file test_pwa.c
 * @brief Tests for Piecewise Affine System — Core L1-L4
 *
 * Tests system creation, point location, well-posedness,
 * PWQ Lyapunov evaluation, and invariant set computation.
 */

#include "pwa_defs.h"
#include "pwa_stability.h"
#include "pwa_geometry.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define TEST(name) printf("  TEST: %s ... ", name)
#define PASS() printf("PASS\n")
#define EPS 1e-10

static int test_system_create(void)
{
    TEST("System creation");
    PWASystem *sys = pwa_system_create(2, 1, 1, 4, 0, 0.1);
    assert(sys != NULL);
    assert(sys->n_state == 2);
    assert(sys->n_input == 1);
    assert(sys->n_output == 1);
    assert(sys->n_allocated == 4);
    assert(sys->n_regions == 0);
    pwa_system_destroy(sys);
    PASS();
    return 1;
}

static int test_add_dynamics_region(void)
{
    TEST("Add dynamics and region");

    PWASystem *sys = pwa_system_create(2, 1, 1, 2, 0, 0.1);
    assert(sys != NULL);

    /* Region 1: stable node */
    double A1[4] = {0.8, 0.1, 0.1, 0.8};
    double B1[2] = {1.0, 0.0};
    double f1[2] = {0.0, 0.0};
    double C1[2] = {1.0, 0.0};
    double D1[1] = {0.0};
    double g1[1] = {0.0};

    int d1 = pwa_add_dynamics(sys, A1, B1, f1, C1, D1, g1);
    assert(d1 == 0);

    /* Region 2: unstable saddle */
    double A2[4] = {0.5, 0.3, 0.3, 0.5};
    double B2[2] = {0.0, 1.0};
    double f2[2] = {0.0, 0.0};

    int d2 = pwa_add_dynamics(sys, A2, B2, f2, NULL, NULL, NULL);
    assert(d2 == 1);

    /* Region 1 constraints: x0 + x1 <= 1, -x0 - x1 <= 1, x0 >= -5, x0 <= 5, x1 >= -5, x1 <= 5 */
    double H1[18] = {1, 1, 1,  -1, -1, 0,  1, 0, 0,  -1, 0, 0,  0, 1, 0,  0, -1, 0};
    double K1[6] = {1, 1, 5, 5, 5, 5};
    int n_cons1 = 6;

    int r1 = pwa_add_region(sys, H1, K1, n_cons1, 0);
    assert(r1 == 0);

    /* Region 2 constraints: x0 + x1 >= -1 (i.e., -x0 - x1 <= 1) */
    double H2[15] = {-1, -1, 0,  1, 0, 0,  -1, 0, 0,  0, 1, 0,  0, -1, 0};
    double K2[5] = {1, 5, 5, 5, 5};
    int n_cons2 = 5;

    int r2 = pwa_add_region(sys, H2, K2, n_cons2, 1);
    assert(r2 == 1);

    assert(sys->n_regions == 2);

    /* Validate system */
    assert(pwa_system_validate(sys) == 0);

    pwa_system_destroy(sys);
    PASS();
    return 1;
}

static int test_point_location(void)
{
    TEST("Point location");

    PWASystem *sys = pwa_system_create(2, 0, 0, 2, 0, 0.1);
    assert(sys != NULL);

    double A1[4] = {0.5, 0.0, 0.0, 0.5};
    double f0[2] = {0.0, 0.0};
    int d1 = pwa_add_dynamics(sys, A1, NULL, f0, NULL, NULL, NULL);
    assert(d1 == 0);

    /* Region 0: x >= 0, y >= 0 */
    double H0[8] = {-1, 0, 0, -1};
    double K0[2] = {0, 0};
    pwa_add_region(sys, H0, K0, 2, 0);

    double A2[4] = {-0.3, 0.0, 0.0, -0.3};
    int d2 = pwa_add_dynamics(sys, A2, NULL, f0, NULL, NULL, NULL);
    assert(d2 == 1);

    /* Region 1: x < 0 or y < 0 (complement of region 0) */
    double H1[8] = {1, 0, 0, 1};
    double K1[2] = {0.1, 0.1};
    pwa_add_region(sys, H1, K1, 2, 1);

    double x_in[2] = {0.5, 0.5};
    int loc = pwa_point_location(sys, x_in, NULL);
    assert(loc == 0);

    double x_out[2] = {-0.5, -0.5};
    loc = pwa_point_location(sys, x_out, NULL);
    assert(loc == 1);

    pwa_system_destroy(sys);
    PASS();
    return 1;
}

static int test_pwq_lyapunov_create_eval(void)
{
    TEST("PWQ Lyapunov create and evaluate");

    int nr = 2, n = 2;
    double *P0 = (double*)calloc(n * n, sizeof(double));
    double *P1 = (double*)calloc(n * n, sizeof(double));
    double *q0 = (double*)calloc(n, sizeof(double));
    double *q1 = (double*)calloc(n, sizeof(double));
    double r[2] = {0.0, 0.0};

    P0[0] = 1.0; P0[3] = 1.0;  /* P0 = I */
    P1[0] = 2.0; P1[3] = 0.5;  /* P1 = diag(2, 0.5) */

    const double *P_arr[2] = {P0, P1};
    const double *q_arr[2] = {q0, q1};

    PWQLyapunov *lyap = pwa_pwq_create(nr, n, P_arr, q_arr, r);
    assert(lyap != NULL);
    assert(lyap->n_regions == 2);
    assert(lyap->n_state == 2);

    double x[2] = {1.0, 2.0};
    double V0 = pwa_pwq_evaluate(lyap, x, 0);
    assert(fabs(V0 - 5.0) < EPS);  /* V0 = x0^2 + x1^2 = 1 + 4 = 5 */

    double V1 = pwa_pwq_evaluate(lyap, x, 1);
    assert(fabs(V1 - 4.0) < EPS);  /* V1 = 2*x0^2 + 0.5*x1^2 = 2 + 2 = 4 */

    pwa_pwq_destroy(lyap);
    free(P0); free(P1); free(q0); free(q1);
    PASS();
    return 1;
}

static int test_matrix_pd_psd(void)
{
    TEST("Matrix PD/PSD check");

    /* Positive definite */
    double M_pd[4] = {2.0, 0.0, 0.0, 2.0};
    assert(pwa_matrix_is_pd(M_pd, 2) == 1);
    assert(pwa_matrix_is_psd(M_pd, 2) == 1);

    /* Positive semidefinite (singular) */
    double M_psd[4] = {1.0, 1.0, 1.0, 1.0};  /* Rank 1, eigenvalue 2 */
    assert(pwa_matrix_is_psd(M_psd, 2) == 1);

    /* Not PSD */
    double M_not[4] = {-1.0, 0.0, 0.0, 2.0};
    assert(pwa_matrix_is_psd(M_not, 2) == 0);
    assert(pwa_matrix_is_pd(M_not, 2) == 0);

    PASS();
    return 1;
}

static int test_convex_hull_2d(void)
{
    TEST("2D convex hull (Graham scan)");

    PWAPoint2D pts[6] = {
        {0, 0}, {1, 0}, {2, 1}, {1, 2}, {0, 1}, {0.5, 0.5}
    };

    PWAConvexHull2D hull;
    memset(&hull, 0, sizeof(hull));
    int n_hull = pwa_convex_hull_2d(pts, 6, &hull);
    assert(n_hull > 0);
    assert(hull.n_points >= 4);  /* Should find at least 4 hull points */
    assert(hull.area > 0.0);

    free(hull.points);
    PASS();
    return 1;
}

static int test_polyhedron_contains(void)
{
    TEST("Polyhedron contains point");

    /* Unit cube: -1 <= x_i <= 1 */
    PWAPolyhedron poly;
    poly.dim = 2;
    poly.n_halfspaces = 4;
    poly.is_bounded = 1;
    poly.is_empty = 0;
    poly.H = (double*)malloc(8 * sizeof(double));
    poly.K = (double*)malloc(4 * sizeof(double));

    /* x0 <= 1, -x0 <= 1, x1 <= 1, -x1 <= 1 */
    poly.H[0] = 1; poly.H[1] = 0; poly.K[0] = 1;
    poly.H[2] = -1; poly.H[3] = 0; poly.K[1] = 1;
    poly.H[4] = 0; poly.H[5] = 1; poly.K[2] = 1;
    poly.H[6] = 0; poly.H[7] = -1; poly.K[3] = 1;

    double x_inside[2] = {0.5, -0.3};
    assert(pwa_polyhedron_contains(&poly, x_inside) == 1);

    double x_outside[2] = {1.5, 0.0};
    assert(pwa_polyhedron_contains(&poly, x_outside) == 0);

    double x_boundary[2] = {1.0, 0.0};
    assert(pwa_polyhedron_contains(&poly, x_boundary) == 1);

    free(poly.H);
    free(poly.K);
    PASS();
    return 1;
}

static int test_polyhedron_intersection(void)
{
    TEST("Polyhedron intersection");

    /* P: x ∈ [0,2], y ∈ [0,2] */
    PWAPolyhedron P, Q, R;
    memset(&P, 0, sizeof(P)); memset(&Q, 0, sizeof(Q)); memset(&R, 0, sizeof(R));

    P.dim = 2; P.n_halfspaces = 4;
    P.H = (double*)malloc(8 * sizeof(double));
    P.K = (double*)malloc(4 * sizeof(double));
    P.H[0]=1;P.H[1]=0;P.K[0]=2; P.H[2]=-1;P.H[3]=0;P.K[1]=0;
    P.H[4]=0;P.H[5]=1;P.K[2]=2; P.H[6]=0;P.H[7]=-1;P.K[3]=0;

    /* Q: x ∈ [1,3], y ∈ [1,3] */
    Q.dim = 2; Q.n_halfspaces = 4;
    Q.H = (double*)malloc(8 * sizeof(double));
    Q.K = (double*)malloc(4 * sizeof(double));
    Q.H[0]=1;Q.H[1]=0;Q.K[0]=3; Q.H[2]=-1;Q.H[3]=0;Q.K[1]=-1;
    Q.H[4]=0;Q.H[5]=1;Q.K[2]=3; Q.H[6]=0;Q.H[7]=-1;Q.K[3]=-1;

    assert(pwa_polyhedron_intersect(&P, &Q) == 1);

    /* Compute intersection */
    pwa_polyhedron_intersection(&P, &Q, &R);
    assert(R.n_halfspaces == 8);
    assert(R.dim == 2);

    /* Check that (1.5, 1.5) is in the intersection */
    double x_mid[2] = {1.5, 1.5};
    assert(pwa_polyhedron_contains(&R, x_mid) == 1);

    /* Check that (0.5, 0.5) is NOT in the intersection */
    double x_out[2] = {0.5, 0.5};
    assert(pwa_polyhedron_contains(&R, x_out) == 0);

    free(P.H); free(P.K);
    free(Q.H); free(Q.K);
    free(R.H); free(R.K);
    PASS();
    return 1;
}

int main(void)
{
    printf("=== PWA System Tests ===\n\n");
    int passed = 0;
    passed += test_system_create();
    passed += test_add_dynamics_region();
    passed += test_point_location();
    passed += test_pwq_lyapunov_create_eval();
    passed += test_matrix_pd_psd();
    passed += test_convex_hull_2d();
    passed += test_polyhedron_contains();
    passed += test_polyhedron_intersection();
    printf("\n=== %d/8 tests passed ===\n", passed);
    return (passed == 8) ? 0 : 1;
}
