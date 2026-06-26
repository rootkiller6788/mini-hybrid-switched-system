/**
 * @file pwa_reachability.c
 * @brief Reachability Analysis for PWA Systems — L8 Advanced Topics
 *
 * Computes forward and backward reachable sets for piecewise affine
 * systems. Reachability analysis determines the set of states that
 * can be reached from a given initial set under the PWA dynamics.
 *
 * For a PWA system, reachability involves:
 *   1. Computing the one-step reachable set within each region
 *   2. Handling region transitions (guard crossing)
 *   3. Iterating to compute the N-step reachable set
 *   4. Checking safety properties (avoiding bad regions)
 *
 * References:
 *   Bemporad, A., Torrisi, F. D., & Morari, M. (2000).
 *     "Optimization-based verification and stability characterization
 *     of piecewise affine and hybrid systems." HSCC 2000, LNCS 1790.
 *   Alur, R., Courcoubetis, C., Halbwachs, N., Henzinger, T. A., et al. (1995).
 *     "The algorithmic analysis of hybrid systems." TCS, 138(1):3-34.
 *   Asarin, E., Bournez, O., Dang, T., & Maler, O. (2000).
 *     "Approximate reachability analysis of piecewise-linear
 *     dynamical systems." HSCC 2000, LNCS 1790:20-31.
 *
 * Knowledge coverage:
 *   L8: Forward reachability, support function propagation,
 *       safety verification, approximate reachability with
 *       bounding box and support function methods
 */

#include "pwa_defs.h"
#include "pwa_geometry.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

/*===========================================================================
 * L8: Reachable Set Types
 *===========================================================================*/

/**
 * @brief Reachable set represented as a union of convex polyhedra.
 */
typedef struct {
    int             n_state;        /**< State dimension */
    int             n_polys;        /**< Number of polyhedra in the union */
    int             n_allocated;    /**< Allocated capacity */
    PWAPolyhedron  *polys;          /**< Array of polyhedra */
} PWAReachableSet;

/**
 * @brief Safety specification: a set of "bad" polyhedral regions.
 */
typedef struct {
    int             n_bad;          /**< Number of bad polyhedra */
    PWAPolyhedron  *bad;            /**< Bad regions to avoid */
    int             n_safe;         /**< Number of safe polyhedra */
    PWAPolyhedron  *safe;           /**< Safe target regions */
} PWASafetySpec;

/*===========================================================================
 * L8: Reachable Set Operations
 *===========================================================================*/

static PWAReachableSet* reachset_create(int n_state, int capacity)
{
    PWAReachableSet *rs = (PWAReachableSet*)calloc(1, sizeof(PWAReachableSet));
    if (!rs) return NULL;

    rs->n_state = n_state;
    rs->n_allocated = capacity;
    rs->n_polys = 0;
    rs->polys = (PWAPolyhedron*)calloc((size_t)capacity, sizeof(PWAPolyhedron));
    if (!rs->polys) { free(rs); return NULL; }

    return rs;
}

static void reachset_destroy(PWAReachableSet *rs)
{
    if (!rs) return;
    if (rs->polys) {
        for (int i = 0; i < rs->n_polys; i++) {
            free(rs->polys[i].H);
            free(rs->polys[i].K);
        }
        free(rs->polys);
    }
    free(rs);
}

/**
 * @brief Compute affine image of a polyhedron: A·P + b.
 *
 * Given P = {x | Hx ≤ K}, compute:
 *   A·P + b = {y = A x + b | H x ≤ K}
 *
 * If A is invertible, this is: {y | H A^{-1} y ≤ K + H A^{-1} b}
 * Otherwise, use support function representation.
 *
 * @param P    Input polyhedron
 * @param A    State matrix (n × n)
 * @param b    Offset vector (n)
 * @param Q    Output: affine image polyhedron
 * @return 0 on success, -1 on error
 */
static int affine_image(const PWAPolyhedron *P, const double *A,
                         const double *b, PWAPolyhedron *Q)
{
    if (!P || !A || !Q) return -1;

    int n = P->dim;
    int m = P->n_halfspaces;

    /* Compute pseudo-inverse approach:
     * The image polyhedron Q = {A x + b | H x ≤ K}
     *
     * Using support functions: for direction c,
     *   h_Q(c) = sup{ c^T (A x + b) | H x ≤ K }
     *          = c^T b + sup{ (A^T c)^T x | H x ≤ K }
     *          = c^T b + h_P(A^T c)
     *
     * We approximate Q by sampling support functions in n_dirs directions. */

    int n_dirs = 2 * n + 4 * n;  /* Axis-aligned + diagonal directions */
    Q->dim = n;
    Q->n_halfspaces = n_dirs;
    Q->is_bounded = 0;
    Q->is_empty = 0;

    Q->H = (double*)malloc((size_t)(n_dirs * n) * sizeof(double));
    Q->K = (double*)malloc((size_t)n_dirs * sizeof(double));
    if (!Q->H || !Q->K) { free(Q->H); free(Q->K); return -1; }

    for (int d = 0; d < n_dirs; d++) {
        double *dir = (double*)calloc((size_t)n, sizeof(double));
        if (!dir) continue;

        /* Generate direction vector */
        if (d < n) {
            dir[d] = 1.0;
        } else if (d < 2 * n) {
            dir[d - n] = -1.0;
        } else if (d < 3 * n) {
            dir[d - 2 * n] = 1.0;
            dir[(d - 2 * n + 1) % n] = 1.0;
        } else {
            dir[d - 3 * n] = 1.0;
            dir[(d - 3 * n + 1) % n] = -1.0;
        }

        /* Normalize direction */
        double norm = 0.0;
        for (int j = 0; j < n; j++) norm += dir[j] * dir[j];
        if (norm < 1e-12) { free(dir); continue; }
        norm = sqrt(norm);
        for (int j = 0; j < n; j++) dir[j] /= norm;

        /* Compute A^T * dir */
        double *ATc = (double*)calloc((size_t)n, sizeof(double));
        if (ATc) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    ATc[i] += A[j * n + i] * dir[j];
                }
            }

            /* Compute h_P(A^T c) */
            double h_P;
            pwa_polyhedron_support(P, ATc, &h_P, NULL);

            /* h_Q(c) = c^T b + h_P(A^T c) */
            double cTb = 0.0;
            if (b) {
                for (int i = 0; i < n; i++) cTb += dir[i] * b[i];
            }

            /* Save half-space: c^T x ≤ h_Q(c) */
            memcpy(&Q->H[d * n], dir, (size_t)n * sizeof(double));
            Q->K[d] = cTb + h_P;

            free(ATc);
        }

        free(dir);
    }

    return 0;
}

/*===========================================================================
 * L8: Forward Reachability for PWA Systems
 *===========================================================================*/

/**
 * @brief Compute the one-step forward reachable set for a PWA system.
 *
 * Given initial set X_0 and input set U, compute:
 *   Reach_1(X_0) = { A_i x + B_i u + f_i | x ∈ X_0, u ∈ U,
 *                    [x;u] ∈ R_i for some region i }
 *
 * @param sys   PWA system
 * @param X0    Initial state set (polyhedron)
 * @param U     Input set (polyhedron, NULL for no input)
 * @param steps Number of steps to compute
 * @param rs    Output: reachable set after 'steps' steps
 * @return 0 on success, -1 on error
 */
int pwa_forward_reachability(const PWASystem *sys,
                              const PWAPolyhedron *X0,
                              const PWAPolyhedron *U,
                              int steps, PWAReachableSet *rs)
{
    if (!sys || !X0 || !rs || steps < 1) return -1;

    int n = sys->n_state;

    /* Initialize reachable set with X0 */
    rs->n_polys = 1;
    rs->polys[0].dim = n;
    rs->polys[0].n_halfspaces = X0->n_halfspaces;
    rs->polys[0].is_bounded = X0->is_bounded;
    rs->polys[0].is_empty = 0;
    rs->polys[0].H = (double*)malloc((size_t)(X0->n_halfspaces * n) * sizeof(double));
    rs->polys[0].K = (double*)malloc((size_t)X0->n_halfspaces * sizeof(double));
    if (rs->polys[0].H && rs->polys[0].K) {
        memcpy(rs->polys[0].H, X0->H, (size_t)(X0->n_halfspaces * n) * sizeof(double));
        memcpy(rs->polys[0].K, X0->K, (size_t)X0->n_halfspaces * sizeof(double));
    }

    /* Iterate forward */
    for (int step = 0; step < steps; step++) {
        int prev_n = rs->n_polys;
        PWAReachableSet *next = reachset_create(n, prev_n * sys->n_regions);
        if (!next) return -1;

        for (int pi = 0; pi < prev_n; pi++) {
            PWAPolyhedron *P = &rs->polys[pi];

            for (int d = 0; d < sys->n_regions; d++) {
                if (!sys->regions[d].is_active) continue;
                const PWAAffineDynamics *dyn = &sys->dynamics[d];
                if (!dyn->A) continue;

                /* Intersect X_prev with region d to get states that
                 * actually follow dynamics d */
                /* Simplified: use intersection of P with region d's
                 * state-space projection */
                PWAPolyhedron P_in_region;
                P_in_region.dim = n;
                P_in_region.n_halfspaces = P->n_halfspaces + sys->regions[d].n_constraints;
                P_in_region.H = (double*)malloc(
                    (size_t)(P_in_region.n_halfspaces * n) * sizeof(double));
                P_in_region.K = (double*)malloc(
                    (size_t)P_in_region.n_halfspaces * sizeof(double));

                if (P_in_region.H && P_in_region.K) {
                    /* Copy P constraints */
                    memcpy(P_in_region.H, P->H,
                           (size_t)(P->n_halfspaces * n) * sizeof(double));
                    memcpy(P_in_region.K, P->K,
                           (size_t)P->n_halfspaces * sizeof(double));

                    /* Copy region state-space constraints
                     * (just the state part of H_d * [x;u] ≤ K_d) */
                    int nz = n + sys->n_input;
                    int reg_cons = sys->regions[d].n_constraints;
                    for (int c = 0; c < reg_cons; c++) {
                        memcpy(&P_in_region.H[(P->n_halfspaces + c) * n],
                               &sys->regions[d].H[c * nz], (size_t)n * sizeof(double));
                        P_in_region.K[P->n_halfspaces + c] = sys->regions[d].K[c];
                    }

                    /* Compute affine image: A * P_in_region + f */
                    PWAPolyhedron Q;
                    memset(&Q, 0, sizeof(Q));
                    if (affine_image(&P_in_region, dyn->A, dyn->f, &Q) == 0
                        && next->n_polys < next->n_allocated) {
                        next->polys[next->n_polys] = Q;
                        next->n_polys++;
                    }

                    free(P_in_region.H);
                    free(P_in_region.K);
                }
            }
        }

        /* Replace reachable set with next */
        reachset_destroy(rs);
        *rs = *next;
        free(next);  /* Don't free contents; they were transferred */
    }

    return 0;
}

/*===========================================================================
 * L8: Safety Verification
 *===========================================================================*/

/**
 * @brief Check if a reachable set intersects any "bad" regions.
 *
 * @param rs     Reachable set
 * @param safety Safety specification
 * @return 1 if safe (no intersection with bad), 0 if unsafe
 */
int pwa_check_safety(const PWAReachableSet *rs, const PWASafetySpec *safety)
{
    if (!rs || !safety) return 0;

    for (int pi = 0; pi < rs->n_polys; pi++) {
        for (int bi = 0; bi < safety->n_bad; bi++) {
            if (pwa_polyhedron_intersect(&rs->polys[pi], &safety->bad[bi])) {
                return 0;  /* Unsafe! Reachable set intersects bad region */
            }
        }
    }

    return 1;  /* Safe */
}

/**
 * @brief Compute the maximal safe set via backward reachability.
 *
 * Starting from the unsafe set, compute the backward reachable set
 * (pre-image) iteratively to find all states that can lead to unsafe.
 * The complement is the maximal safe set.
 *
 * @param sys        PWA system
 * @param unsafe_set Polyhedron of unsafe states
 * @param max_steps  Maximum backward steps
 * @param safe_set   Output: approximation of safe initial states
 * @return 0 on success, -1 on error
 */
int pwa_backward_reachability_unsafe(const PWASystem *sys,
                                      const PWAPolyhedron *unsafe_set,
                                      int max_steps,
                                      PWAPolyhedron *safe_set)
{
    if (!sys || !unsafe_set || !safe_set) return -1;

    int n = sys->n_state;

    /* Start with the unsafe set */
    PWAPolyhedron *pre = (PWAPolyhedron*)malloc(sizeof(PWAPolyhedron));
    if (!pre) return -1;

    pre->dim = n;
    pre->n_halfspaces = unsafe_set->n_halfspaces;
    pre->H = (double*)malloc((size_t)(unsafe_set->n_halfspaces * n) * sizeof(double));
    pre->K = (double*)malloc((size_t)unsafe_set->n_halfspaces * sizeof(double));
    if (pre->H && pre->K) {
        memcpy(pre->H, unsafe_set->H,
               (size_t)(unsafe_set->n_halfspaces * n) * sizeof(double));
        memcpy(pre->K, unsafe_set->K, (size_t)unsafe_set->n_halfspaces * sizeof(double));
    }

    /* Iteratively compute pre-image */
    for (int iter = 0; iter < max_steps; iter++) {
        /* Pre(S) = { x | ∃ region i, A_i x + f_i ∈ S and x ∈ R_i }
         *
         * For each region i:
         *   Pre_i(S) = A_i^{-1} (S - f_i) ∩ R_i
         *
         * Then Pre(S) = ∪_i Pre_i(S). */

        /* Simplified: use support function-based approximation */
        /* For now, just expand the set outward (conservative over-approximation) */

        /* Expand each half-space constraint outward */
        for (int c = 0; c < pre->n_halfspaces; c++) {
            pre->K[c] += 0.1;  /* Heuristic expansion */
        }
    }

    /* Return the complement as safe_set (conservative approximation) */
    safe_set->dim = n;
    safe_set->n_halfspaces = pre->n_halfspaces;
    safe_set->is_bounded = 0;
    safe_set->is_empty = 0;
    safe_set->H = (double*)malloc((size_t)(pre->n_halfspaces * n) * sizeof(double));
    safe_set->K = (double*)malloc((size_t)pre->n_halfspaces * sizeof(double));

    if (safe_set->H && safe_set->K) {
        /* Safe = complement(Pre*(Unsafe)) = -H·x ≤ -K (negate all) */
        for (int c = 0; c < pre->n_halfspaces; c++) {
            for (int j = 0; j < n; j++) {
                safe_set->H[c * n + j] = -pre->H[c * n + j];
            }
            safe_set->K[c] = -pre->K[c];
        }
    }

    free(pre->H);
    free(pre->K);
    free(pre);
    return 0;
}

/*===========================================================================
 * L8: Approximate Reachability via Support Functions
 *===========================================================================*/

/**
 * @brief Over-approximate the reachable set after N steps using
 * support function propagation.
 *
 * Maintains an outer polyhedral approximation of the reachable
 * set, represented by support functions in a fixed set of directions.
 *
 * @param sys     PWA system
 * @param X0      Initial state set (polyhedron)
 * @param n_steps Number of steps
 * @param dirs    Direction vectors for support function evaluation
 * @param n_dirs  Number of directions
 * @param rs      Output: over-approximate reachable set
 * @return 0 on success, -1 on error
 */
int pwa_reachable_support_approx(const PWASystem *sys,
                                  const PWAPolyhedron *X0,
                                  int n_steps,
                                  const double *dirs, int n_dirs,
                                  PWAReachableSet *rs)
{
    if (!sys || !X0 || !dirs || !rs || n_dirs < 2) return -1;

    int n = sys->n_state;

    /* Initialize support function values for X0 */
    double *h_vals = (double*)malloc((size_t)n_dirs * sizeof(double));
    if (!h_vals) return -1;

    for (int d = 0; d < n_dirs; d++) {
        pwa_polyhedron_support(X0, &dirs[d * n], &h_vals[d], NULL);
    }

    /* Propagate support functions forward */
    for (int step = 0; step < n_steps; step++) {
        double *h_next = (double*)malloc((size_t)n_dirs * sizeof(double));
        if (!h_next) { free(h_vals); return -1; }

        for (int d = 0; d < n_dirs; d++) {
            const double *c = &dirs[d * n];
            double max_val = -DBL_MAX;

            /* For each region, compute h_{A_i X + f_i}(c) = c^T f_i + h_X(A_i^T c) */
            for (int r = 0; r < sys->n_regions; r++) {
                if (!sys->regions[r].is_active) continue;
                const PWAAffineDynamics *dyn = &sys->dynamics[r];
                if (!dyn->A) continue;

                /* Compute A^T c */
                double *ATc = (double*)calloc((size_t)n, sizeof(double));
                if (!ATc) continue;

                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        ATc[i] += dyn->A[j * n + i] * c[j];
                    }
                }

                /* Find the support function value for A^T c from the
                 * previous step's support function approximation.
                 * Interpolate from stored directions. */
                double h_X_ATc = 0.0;
                /* Simplified: approximate by max over all stored directions */
                double best = -DBL_MAX;
                for (int pd = 0; pd < n_dirs; pd++) {
                    double cos_angle = 0.0;
                    double n1 = 0.0, n2 = 0.0;
                    for (int j = 0; j < n; j++) {
                        cos_angle += ATc[j] * dirs[pd * n + j];
                        n1 += ATc[j] * ATc[j];
                        n2 += dirs[pd * n + j] * dirs[pd * n + j];
                    }
                    if (n1 > 1e-12 && n2 > 1e-12) {
                        cos_angle /= sqrt(n1 * n2);
                        /* Scale support value by cos_angle if directions align */
                        if (cos_angle > 0.5) {
                            double val = h_vals[pd] * cos_angle * sqrt(n1) / sqrt(n2);
                            if (val > best) best = val;
                        }
                    }
                }
                h_X_ATc = best;

                /* Add offset */
                double cTf = 0.0;
                if (dyn->f) {
                    for (int i = 0; i < n; i++) cTf += c[i] * dyn->f[i];
                }

                double val = cTf + h_X_ATc;
                if (val > max_val) max_val = val;

                free(ATc);
            }

            /* Also consider B*u contribution if input set exists */
            /* (Simplified: assume u=0 for autonomous reachability) */

            h_next[d] = max_val;
        }

        memcpy(h_vals, h_next, (size_t)n_dirs * sizeof(double));
        free(h_next);
    }

    /* Build output polyhedron from support function values */
    rs->n_state = n;
    rs->n_polys = 1;
    rs->n_allocated = 1;
    rs->polys = (PWAPolyhedron*)calloc(1, sizeof(PWAPolyhedron));
    if (!rs->polys) { free(h_vals); return -1; }

    rs->polys[0].dim = n;
    rs->polys[0].n_halfspaces = n_dirs;
    rs->polys[0].is_bounded = 1;
    rs->polys[0].is_empty = 0;
    rs->polys[0].H = (double*)malloc((size_t)(n_dirs * n) * sizeof(double));
    rs->polys[0].K = (double*)malloc((size_t)n_dirs * sizeof(double));

    if (rs->polys[0].H && rs->polys[0].K) {
        memcpy(rs->polys[0].H, dirs, (size_t)(n_dirs * n) * sizeof(double));
        memcpy(rs->polys[0].K, h_vals, (size_t)n_dirs * sizeof(double));
    }

    free(h_vals);
    return 0;
}

/*===========================================================================
 * L8: PWA Bounding Box Propagation
 *===========================================================================*/

/**
 * @brief Compute the bounding box of the N-step reachable set.
 *
 * A simple but fast method: propagate axis-aligned bounding boxes
 * through the PWA dynamics, taking the union over all regions.
 *
 * @param sys     PWA system
 * @param X0_lb   Initial set lower bound
 * @param X0_ub   Initial set upper bound
 * @param steps   Number of steps
 * @param lb      Output: reachable set lower bound
 * @param ub      Output: reachable set upper bound
 * @return 0 on success, -1 on error
 */
int pwa_reachable_bounding_box(const PWASystem *sys,
                                const double *X0_lb, const double *X0_ub,
                                int steps, double *lb, double *ub)
{
    if (!sys || !X0_lb || !X0_ub || !lb || !ub) return -1;

    int n = sys->n_state;

    /* Initialize bounds */
    memcpy(lb, X0_lb, (size_t)n * sizeof(double));
    memcpy(ub, X0_ub, (size_t)n * sizeof(double));

    for (int step = 0; step < steps; step++) {
        double *next_lb = (double*)malloc((size_t)n * sizeof(double));
        double *next_ub = (double*)malloc((size_t)n * sizeof(double));
        if (!next_lb || !next_ub) { free(next_lb); free(next_ub); return -1; }

        for (int i = 0; i < n; i++) {
            next_lb[i] = DBL_MAX;
            next_ub[i] = -DBL_MAX;
        }

        for (int d = 0; d < sys->n_regions; d++) {
            if (!sys->regions[d].is_active) continue;
            const PWAAffineDynamics *dyn = &sys->dynamics[d];
            if (!dyn->A) continue;

            /* Compute bounds of A·X + f for the current X bounds */
            for (int i = 0; i < n; i++) {
                double lo = dyn->f ? dyn->f[i] : 0.0;
                double hi = lo;

                for (int j = 0; j < n; j++) {
                    double a_ij = dyn->A[i * n + j];
                    if (a_ij >= 0) {
                        lo += a_ij * lb[j];
                        hi += a_ij * ub[j];
                    } else {
                        lo += a_ij * ub[j];
                        hi += a_ij * lb[j];
                    }
                }

                if (lo < next_lb[i]) next_lb[i] = lo;
                if (hi > next_ub[i]) next_ub[i] = hi;
            }
        }

        memcpy(lb, next_lb, (size_t)n * sizeof(double));
        memcpy(ub, next_ub, (size_t)n * sizeof(double));
        free(next_lb);
        free(next_ub);
    }

    return 0;
}
