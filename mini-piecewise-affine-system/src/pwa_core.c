/**
 * @file pwa_core.c
 * @brief Piecewise Affine System Core Operations — L1+L2 Implementation
 *
 * Implements core PWA system construction, point location (region
 * membership testing), well-posedness checking, region adjacency
 * computation, and system validation.
 *
 * Knowledge coverage:
 *   L1: PWA system creation/destruction, region/dynamics management
 *   L2: Point location in polyhedral partition, well-posedness
 *       analysis (overlap detection, gap detection), adjacency
 *
 * References:
 *   Johansson, M. (2003). "Piecewise Linear Control Systems."
 *     Springer-Verlag, Chapter 2.
 *   Sontag, E. D. (1981). "Nonlinear regulation: The piecewise
 *     linear approach." IEEE TAC, 26(2):346-358.
 */

#include "pwa_defs.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/*===========================================================================
 * L1: PWA System Creation and Destruction
 *===========================================================================*/

PWASystem* pwa_system_create(int n_state, int n_input, int n_output,
                              int n_regions_max, int is_continuous, double dt)
{
    if (n_state <= 0 || n_input < 0 || n_output < 0 || n_regions_max <= 0) {
        return NULL;
    }

    PWASystem *sys = (PWASystem*)calloc(1, sizeof(PWASystem));
    if (!sys) return NULL;

    sys->n_state = n_state;
    sys->n_input = n_input;
    sys->n_output = n_output;
    sys->n_regions = 0;
    sys->n_allocated = n_regions_max;
    sys->is_continuous = is_continuous ? 1 : 0;
    sys->dt = dt;
    sys->current_region = -1;

    /* Allocate region array */
    sys->regions = (PWARegion*)calloc((size_t)n_regions_max, sizeof(PWARegion));
    if (!sys->regions) {
        free(sys);
        return NULL;
    }

    /* Allocate dynamics array */
    sys->dynamics = (PWAAffineDynamics*)calloc((size_t)n_regions_max,
                                                sizeof(PWAAffineDynamics));
    if (!sys->dynamics) {
        free(sys->regions);
        free(sys);
        return NULL;
    }

    /* Initialize regions and dynamics with safe defaults */
    for (int i = 0; i < n_regions_max; i++) {
        sys->regions[i].id = i;
        sys->regions[i].n_state = n_state;
        sys->regions[i].n_input = n_input;
        sys->regions[i].n_constraints = 0;
        sys->regions[i].H = NULL;
        sys->regions[i].K = NULL;
        sys->regions[i].dynamics_id = -1;
        sys->regions[i].is_active = 0;

        sys->dynamics[i].n_state = n_state;
        sys->dynamics[i].n_input = n_input;
        sys->dynamics[i].n_output = n_output;
        sys->dynamics[i].A = NULL;
        sys->dynamics[i].B = NULL;
        sys->dynamics[i].f = NULL;
        sys->dynamics[i].C = NULL;
        sys->dynamics[i].D = NULL;
        sys->dynamics[i].g = NULL;
    }

    /* Allocate state/input arrays */
    sys->current_state = (double*)calloc((size_t)n_state, sizeof(double));
    sys->current_input = (double*)calloc((size_t)n_input, sizeof(double));
    sys->x_min = (double*)calloc((size_t)n_state, sizeof(double));
    sys->x_max = (double*)calloc((size_t)n_state, sizeof(double));
    sys->u_min = (double*)calloc((size_t)n_input, sizeof(double));
    sys->u_max = (double*)calloc((size_t)n_input, sizeof(double));

    if (!sys->current_state || !sys->current_input ||
        !sys->x_min || !sys->x_max || !sys->u_min || !sys->u_max) {
        pwa_system_destroy(sys);
        return NULL;
    }

    return sys;
}

void pwa_system_destroy(PWASystem *sys)
{
    if (!sys) return;

    for (int i = 0; i < sys->n_allocated; i++) {
        free(sys->regions[i].H);
        free(sys->regions[i].K);
        free(sys->dynamics[i].A);
        free(sys->dynamics[i].B);
        free(sys->dynamics[i].f);
        free(sys->dynamics[i].C);
        free(sys->dynamics[i].D);
        free(sys->dynamics[i].g);
    }

    free(sys->regions);
    free(sys->dynamics);
    free(sys->current_state);
    free(sys->current_input);
    free(sys->x_min);
    free(sys->x_max);
    free(sys->u_min);
    free(sys->u_max);
    free(sys);
}

/*===========================================================================
 * L1: Dynamics and Region Management
 *===========================================================================*/

static double* pwa_alloc_matrix(int rows, int cols)
{
    if (rows <= 0 || cols <= 0) return NULL;
    return (double*)calloc((size_t)(rows * cols), sizeof(double));
}

static double* pwa_alloc_vector(int n)
{
    if (n <= 0) return NULL;
    return (double*)calloc((size_t)n, sizeof(double));
}

int pwa_add_dynamics(PWASystem *sys,
                      const double *A, const double *B, const double *f,
                      const double *C, const double *D, const double *g)
{
    if (!sys) return -1;
    if (sys->n_regions >= sys->n_allocated) return -1;

    int idx = sys->n_regions;
    int n = sys->n_state;
    int m = sys->n_input;
    int p = sys->n_output;

    PWAAffineDynamics *dyn = &sys->dynamics[idx];

    /* Allocate and copy A */
    dyn->A = pwa_alloc_matrix(n, n);
    if (A && dyn->A) {
        memcpy(dyn->A, A, (size_t)(n * n) * sizeof(double));
    }

    /* Allocate and copy B */
    dyn->B = pwa_alloc_matrix(n, m);
    if (B && dyn->B) {
        memcpy(dyn->B, B, (size_t)(n * m) * sizeof(double));
    }

    /* Allocate and copy f */
    dyn->f = pwa_alloc_vector(n);
    if (f && dyn->f) {
        memcpy(dyn->f, f, (size_t)n * sizeof(double));
    }

    /* Allocate and copy C */
    dyn->C = pwa_alloc_matrix(p, n);
    if (C && dyn->C) {
        memcpy(dyn->C, C, (size_t)(p * n) * sizeof(double));
    }

    /* Allocate and copy D */
    dyn->D = pwa_alloc_matrix(p, m);
    if (D && dyn->D) {
        memcpy(dyn->D, D, (size_t)(p * m) * sizeof(double));
    }

    /* Allocate and copy g */
    dyn->g = pwa_alloc_vector(p);
    if (g && dyn->g) {
        memcpy(dyn->g, g, (size_t)p * sizeof(double));
    }

    sys->n_regions++;
    return idx;
}

int pwa_add_region(PWASystem *sys,
                    const double *H, const double *K, int n_cons,
                    int dynamics_id)
{
    if (!sys || !H || !K || n_cons <= 0) return -1;
    if (dynamics_id < 0 || dynamics_id >= sys->n_regions) return -1;

    int idx = dynamics_id;  /* Region and dynamics share the same index */
    if (idx >= sys->n_allocated) return -1;

    PWARegion *reg = &sys->regions[idx];
    int nz = sys->n_state + sys->n_input;

    /* Free existing constraint data if re-adding */
    free(reg->H);
    free(reg->K);

    reg->n_constraints = n_cons;
    reg->H = pwa_alloc_matrix(n_cons, nz);
    reg->K = pwa_alloc_vector(n_cons);
    reg->dynamics_id = dynamics_id;
    reg->is_active = 1;

    if (reg->H) {
        memcpy(reg->H, H, (size_t)(n_cons * nz) * sizeof(double));
    }
    if (reg->K) {
        memcpy(reg->K, K, (size_t)n_cons * sizeof(double));
    }

    return idx;
}

/*===========================================================================
 * L1: System Validation
 *===========================================================================*/

int pwa_system_validate(const PWASystem *sys)
{
    if (!sys) return -1;
    if (sys->n_state <= 0) return -2;
    if (sys->n_regions <= 0) return -3;

    /* Check each active region has valid dynamics and constraints */
    for (int i = 0; i < sys->n_regions; i++) {
        const PWARegion *reg = &sys->regions[i];
        if (!reg->is_active) continue;

        if (reg->dynamics_id < 0 || reg->dynamics_id >= sys->n_regions) {
            return -4;
        }
        if (reg->n_constraints <= 0) {
            return -5;
        }
        if (!reg->H || !reg->K) {
            return -6;
        }

        const PWAAffineDynamics *dyn = &sys->dynamics[reg->dynamics_id];
        if (!dyn->A) {
            return -7;
        }
    }

    /* Check dimension consistency */
    for (int i = 0; i < sys->n_regions; i++) {
        if (sys->regions[i].n_state != sys->n_state) return -8;
        if (sys->regions[i].n_input != sys->n_input) return -9;
    }

    return 0;
}

/*===========================================================================
 * L2: Point in Region (Polyhedron Membership)
 *===========================================================================*/

int pwa_point_in_region(const PWARegion *region, const double *z)
{
    if (!region || !z || !region->H || !region->K) return 0;

    int n_cons = region->n_constraints;
    int nz = region->n_state + region->n_input;

    for (int i = 0; i < n_cons; i++) {
        double sum = 0.0;
        const double *Hi = &region->H[i * nz];
        for (int j = 0; j < nz; j++) {
            sum += Hi[j] * z[j];
        }
        /* H_i * z ≤ K_i  — allow small tolerance for floating point */
        if (sum > region->K[i] + 1e-10) {
            return 0;
        }
    }

    return 1;
}

/*===========================================================================
 * L2: Point Location in Polyhedral Partition
 *===========================================================================*/

int pwa_point_location(const PWASystem *sys, const double *x, const double *u)
{
    if (!sys || !x) return -1;

    int n_in = sys->n_input;
    int nz = sys->n_state + n_in;
    double *z = (double*)malloc((size_t)nz * sizeof(double));
    if (!z) return -1;

    /* Build extended vector z = [x; u] */
    memcpy(z, x, (size_t)sys->n_state * sizeof(double));
    if (u && n_in > 0) {
        memcpy(z + sys->n_state, u, (size_t)n_in * sizeof(double));
    } else if (n_in > 0) {
        memset(z + sys->n_state, 0, (size_t)n_in * sizeof(double));
    }

    int result = -1;
    for (int i = 0; i < sys->n_regions; i++) {
        if (!sys->regions[i].is_active) continue;
        if (pwa_point_in_region(&sys->regions[i], z)) {
            result = i;
            break;
        }
    }

    free(z);
    return result;
}

/*===========================================================================
 * L2: Well-Posedness Checking
 *===========================================================================*/

/**
 * Check if two regions overlap in their interiors.
 *
 * Two regions R_i and R_j have overlapping interiors if
 * the polyhedron R_i ∩ R_j has dimension n+m.
 *
 * Simplified check: test if their Chebyshev centers are in both.
 * For exact check, would need LP: test if ∃z s.t. H_i z ≤ K_i - ε
 * and H_j z ≤ K_j - ε.
 */
static int regions_overlap_interior(const PWARegion *r1, const PWARegion *r2)
{
    if (!r1 || !r2) return 0;
    if (!r1->is_active || !r2->is_active) return 0;
    if (r1 == r2) return 0;

    int nz = r1->n_state + r1->n_input;

    /* Simple check: sample the centroid estimate of each region
     * by averaging constraint vertices, and test membership */
    /* For each constraint of r1, compute the point that satisfies all
     * constraints as equalities except one, and check if it's in r2.
     * This is expensive, so we use a simpler heuristic: */

    /* Heuristic: compute the Chebyshev center of r1 (by solving LP)
     * and check if it's strictly inside r2. Since we don't have LP
     * solver inline here, use the average of some feasible points. */

    /* Sample the origin offset. Build a rough center estimate. */
    double *center = (double*)calloc((size_t)nz, sizeof(double));
    if (!center) return 0;

    /* Use constraint midpoints heuristic */
    int count = 0;
    for (int c = 0; c < r1->n_constraints && count < 10; c++) {
        /* For constraint c: H_c·z = K_c, find a point satisfying it
         * along with H_0·z ≤ K_0 ... (simple: project origin) */
        double norm_sq = 0.0;
        const double *Hc = &r1->H[c * nz];
        for (int j = 0; j < nz; j++) {
            norm_sq += Hc[j] * Hc[j];
        }
        if (norm_sq < 1e-12) continue;

        double alpha = r1->K[c] / norm_sq;
        for (int j = 0; j < nz; j++) {
            center[j] += alpha * Hc[j];
        }
        count++;
    }

    if (count > 0) {
        for (int j = 0; j < nz; j++) {
            center[j] /= (double)count;
        }
    }

    /* Check if center is in both regions */
    int in_r1 = pwa_point_in_region(r1, center);
    int in_r2 = pwa_point_in_region(r2, center);

    free(center);

    return (in_r1 && in_r2) ? 1 : 0;
}

PWAWellPosedStatus pwa_check_wellposed(const PWASystem *sys)
{
    if (!sys) return PWA_WELLPOSED_INVALID;
    if (sys->n_regions < 1) return PWA_WELLPOSED_INVALID;

    int has_overlap = 0;

    /* Check for region overlaps */
    for (int i = 0; i < sys->n_regions; i++) {
        if (!sys->regions[i].is_active) continue;
        for (int j = i + 1; j < sys->n_regions; j++) {
            if (!sys->regions[j].is_active) continue;
            if (regions_overlap_interior(&sys->regions[i], &sys->regions[j])) {
                has_overlap = 1;
                /* Don't break - continue scanning for diagnostic purposes,
                 * but we'll report overlap. However, some controlled overlap
                 * at boundaries is acceptable. Check strict interior: */
                break;
            }
        }
        if (has_overlap) break;
    }

    /* Check for gaps: sample a grid in the domain and check
     * if any point falls outside all regions */
    int has_gap = 0;
    int n_samples = 5;  /* 5 samples per dimension */
    int n_state = sys->n_state;
    int n_input = sys->n_input;
    int nz = n_state + n_input;

    if (nz > 0 && nz <= 4) {  /* Only do grid search for low dimensions */
        int total_samples = 1;
        for (int d = 0; d < nz; d++) total_samples *= n_samples;

        double *z = (double*)malloc((size_t)nz * sizeof(double));
        if (z) {
            for (int s = 0; s < total_samples && !has_gap; s++) {
                int tmp = s;
                for (int d = 0; d < nz; d++) {
                    int idx = tmp % n_samples;
                    tmp /= n_samples;
                    /* Interpolate between bounds */
                    double lo = -10.0, hi = 10.0;  /* Default domain */
                    if (d < n_state) {
                        lo = sys->x_min ? sys->x_min[d] : -10.0;
                        hi = sys->x_max ? sys->x_max[d] : 10.0;
                    } else {
                        int ud = d - n_state;
                        lo = sys->u_min ? sys->u_min[ud] : -10.0;
                        hi = sys->u_max ? sys->u_max[ud] : 10.0;
                    }
                    z[d] = lo + (hi - lo) * (double)idx / (double)(n_samples - 1);
                }

                int in_any = 0;
                for (int i = 0; i < sys->n_regions; i++) {
                    if (!sys->regions[i].is_active) continue;
                    if (pwa_point_in_region(&sys->regions[i], z)) {
                        in_any = 1;
                        break;
                    }
                }
                if (!in_any) {
                    has_gap = 1;
                }
            }
            free(z);
        }
    }

    if (has_overlap) return PWA_WELLPOSED_OVERLAP;
    if (has_gap) return PWA_WELLPOSED_GAPS;
    return PWA_WELLPOSED_OK;
}

/*===========================================================================
 * L2: Region Adjacency Computation
 *===========================================================================*/

/**
 * Check if two polyhedral regions share a facet.
 *
 * Two regions i and j are adjacent if their intersection is
 * an (n+m-1)-dimensional polytope — i.e., they share a facet.
 *
 * This simplified check: test if there exists a point z that
 * satisfies all constraints of both regions, with exactly one
 * constraint active for each (or shared between them).
 */
static int regions_are_adjacent(const PWARegion *r1, const PWARegion *r2)
{
    if (!r1 || !r2) return 0;
    if (!r1->is_active || !r2->is_active) return 0;

    int nz = r1->n_state + r1->n_input;
    double *z = (double*)calloc((size_t)nz, sizeof(double));
    if (!z) return 0;

    int adjacent = 0;

    /* For each possible shared constraint pair */
    for (int c1 = 0; c1 < r1->n_constraints && !adjacent; c1++) {
        for (int c2 = 0; c2 < r2->n_constraints && !adjacent; c2++) {
            /* Check if the two half-spaces are parallel and opposite:
             * H1 ≈ -α * H2 for some α > 0, and K1 ≈ -α * K2 */

            const double *h1 = &r1->H[c1 * nz];
            const double *h2 = &r2->H[c2 * nz];

            /* Compute dot product of normals */
            double dot = 0.0, n1_sq = 0.0, n2_sq = 0.0;
            for (int j = 0; j < nz; j++) {
                dot += h1[j] * h2[j];
                n1_sq += h1[j] * h1[j];
                n2_sq += h2[j] * h2[j];
            }

            /* Check if normals are nearly opposite (cos angle ≈ -1) */
            if (n1_sq < 1e-12 || n2_sq < 1e-12) continue;

            double cos_angle = dot / sqrt(n1_sq * n2_sq);

            if (cos_angle < -0.9) {
                /* Potential shared facet. Check if the two half-spaces
                 * define the same hyperplane with opposite orientation. */
                double n1 = sqrt(n1_sq);
                double n2 = sqrt(n2_sq);

                /* Normalize and check if -h1/n1 ≈ h2/n2 */
                int aligned = 1;
                for (int j = 0; j < nz; j++) {
                    double diff = (-h1[j] / n1) - (h2[j] / n2);
                    if (fabs(diff) > 0.05) {
                        aligned = 0;
                        break;
                    }
                }
                if (aligned) {
                    /* Check K values: K1/n1 + K2/n2 ≈ 0 for opposite half-spaces */
                    double k_sum = r1->K[c1] / n1 + r2->K[c2] / n2;
                    if (fabs(k_sum) < 0.1) {
                        adjacent = 1;
                        break;
                    }
                }
            }
        }
    }

    free(z);
    return adjacent;
}

int pwa_compute_adjacency(const PWASystem *sys, int *adjacency)
{
    if (!sys || !adjacency) return -1;

    int n = sys->n_regions;
    memset(adjacency, 0, (size_t)(n * n) * sizeof(int));

    for (int i = 0; i < n; i++) {
        if (!sys->regions[i].is_active) continue;
        for (int j = i + 1; j < n; j++) {
            if (!sys->regions[j].is_active) continue;
            if (regions_are_adjacent(&sys->regions[i], &sys->regions[j])) {
                adjacency[i * n + j] = 1;
                adjacency[j * n + i] = 1;
            }
        }
    }

    return 0;
}

/*===========================================================================
 * L2: System Printing
 *===========================================================================*/

void pwa_system_print(const PWASystem *sys, int verbose)
{
    if (!sys) {
        printf("PWASystem: NULL\n");
        return;
    }

    printf("=== PWA System ===\n");
    printf("Dimensions: n_state=%d, n_input=%d, n_output=%d\n",
           sys->n_state, sys->n_input, sys->n_output);
    printf("Type: %s\n", sys->is_continuous ? "Continuous-time" : "Discrete-time");
    if (!sys->is_continuous) {
        printf("Sampling time: dt=%.6f\n", sys->dt);
    }
    printf("Regions: %d active / %d allocated\n",
           sys->n_regions, sys->n_allocated);
    printf("Current region: %d\n", sys->current_region);

    if (!verbose) return;

    for (int i = 0; i < sys->n_regions; i++) {
        const PWARegion *reg = &sys->regions[i];
        printf("\n--- Region %d ---\n", i);
        printf("  Active: %s\n", reg->is_active ? "yes" : "no");
        printf("  Constraints: %d\n", reg->n_constraints);
        printf("  Dynamics ID: %d\n", reg->dynamics_id);

        /* Print constraints */
        if (reg->H && reg->K && reg->n_constraints > 0) {
            int nz = reg->n_state + reg->n_input;
            for (int c = 0; c < reg->n_constraints; c++) {
                printf("  C%02d: [", c);
                for (int j = 0; j < nz; j++) {
                    printf(" %+.3f", reg->H[c * nz + j]);
                }
                printf(" ]·z ≤ %+.3f\n", reg->K[c]);
            }
        }

        /* Print dynamics */
        if (reg->dynamics_id >= 0 && reg->dynamics_id < sys->n_regions) {
            const PWAAffineDynamics *dyn = &sys->dynamics[reg->dynamics_id];
            printf("  Dynamics: A[%dx%d] f[%d]\n",
                   dyn->n_state, dyn->n_state, dyn->n_state);
            if (dyn->A && verbose >= 2) {
                for (int r = 0; r < dyn->n_state; r++) {
                    printf("    A[%d]:", r);
                    for (int c = 0; c < dyn->n_state; c++) {
                        printf(" %+.3f", dyn->A[r * dyn->n_state + c]);
                    }
                    printf("\n");
                }
            }
            if (dyn->f) {
                printf("  f: ");
                for (int r = 0; r < dyn->n_state; r++) {
                    printf(" %+.3f", dyn->f[r]);
                }
                printf("\n");
            }
        }
    }
}
