/**
 * @file pwa_geometry.c
 * @brief Polyhedral Geometry Operations — L3 Mathematical Structures
 *
 * Implements computational geometry for polyhedral partitions used
 * in piecewise affine systems: half-space operations, convex hull
 * algorithms, polyhedron operations (intersection, membership, subset),
 * simplex LP, Fourier-Motzkin elimination, support functions,
 * Minkowski sum, and Hausdorff distance.
 *
 * Knowledge coverage:
 *   L3: All polyhedral operations, convex hull, half-space intersection,
 *       vertex enumeration, Fourier-Motzkin projection, Minkowski sum,
 *       Hausdorff distance, support functions
 *
 * References:
 *   Preparata & Shamos (1985). "Computational Geometry." Springer.
 *   Ziegler (1995). "Lectures on Polytopes." Springer GTM 152.
 *   Boyd & Vandenberghe (2004). "Convex Optimization." Cambridge.
 */

#include "pwa_geometry.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*===========================================================================
 * L3: Half-Space and Hyperplane Operations
 *===========================================================================*/

double pwa_halfspace_eval(const PWAHalfSpace *hs, const double *x)
{
    if (!hs || !x) return 0.0;
    double val = -hs->b;  /* a·x - b */
    for (int i = 0; i < hs->dim; i++) {
        val += hs->a[i] * x[i];
    }
    return val;
}

double pwa_hyperplane_distance(const PWAHyperplane *hp, const double *x)
{
    if (!hp || !x) return 0.0;
    double dot = 0.0;
    double norm_sq = 0.0;
    for (int i = 0; i < hp->dim; i++) {
        dot += hp->a[i] * x[i];
        norm_sq += hp->a[i] * hp->a[i];
    }
    if (norm_sq < 1e-15) return 0.0;
    return (dot - hp->b) / sqrt(norm_sq);
}

/*===========================================================================
 * L3: Polyhedron Membership and Basic Properties
 *===========================================================================*/

int pwa_polyhedron_contains(const PWAPolyhedron *poly, const double *x)
{
    if (!poly || !x || !poly->H || !poly->K) return 0;

    for (int i = 0; i < poly->n_halfspaces; i++) {
        double sum = 0.0;
        const double *Hi = &poly->H[i * poly->dim];
        for (int j = 0; j < poly->dim; j++) {
            sum += Hi[j] * x[j];
        }
        if (sum > poly->K[i] + 1e-10) {
            return 0;
        }
    }
    return 1;
}

int pwa_polyhedron_chebyshev_center(const PWAPolyhedron *poly,
                                     double *center, double *radius)
{
    if (!poly || !center || !radius) return -1;
    if (poly->n_halfspaces < 1) return -1;

    /* Simplified: for small dimensions, use a gradient-based approach.
     * Full LP implementation would use a proper simplex solver.
     * Here we use an iterative projection method for the Chebyshev center.
     *
     * The problem: max r s.t. H_i·x + r·||H_i|| ≤ K_i for all i
     *
     * We solve this via a simplified perceptron-like iterative algorithm
     * followed by a simple LP for small dimensions.
     */

    int m = poly->n_halfspaces;
    int d = poly->dim;

    /* Precompute row norms of H */
    double *row_norms = (double*)malloc((size_t)m * sizeof(double));
    if (!row_norms) return -1;

    for (int i = 0; i < m; i++) {
        double ns = 0.0;
        const double *Hi = &poly->H[i * d];
        for (int j = 0; j < d; j++) {
            ns += Hi[j] * Hi[j];
        }
        row_norms[i] = sqrt(ns);
    }

    /* Start with origin, then adjust */
    for (int j = 0; j < d; j++) center[j] = 0.0;

    /* Iterative feasible point finding via projection */
    for (int iter = 0; iter < 1000; iter++) {
        int violated = -1;
        double max_violation = 0.0;

        for (int i = 0; i < m; i++) {
            double prod = 0.0;
            const double *Hi = &poly->H[i * d];
            for (int j = 0; j < d; j++) {
                prod += Hi[j] * center[j];
            }
            double violation = prod - poly->K[i];
            if (violation > max_violation) {
                max_violation = violation;
                violated = i;
            }
        }

        if (max_violation < 1e-10) break;  /* Feasible */

        if (violated >= 0 && row_norms[violated] > 1e-12) {
            const double *Hv = &poly->H[violated * d];
            double step = max_violation / (row_norms[violated] * row_norms[violated]);
            for (int j = 0; j < d; j++) {
                center[j] -= step * Hv[j];
            }
        }
    }

    /* Now compute the maximum inscribed radius */
    *radius = DBL_MAX;
    for (int i = 0; i < m; i++) {
        double prod = 0.0;
        const double *Hi = &poly->H[i * d];
        for (int j = 0; j < d; j++) {
            prod += Hi[j] * center[j];
        }
        double slack = poly->K[i] - prod;
        double r_i = (row_norms[i] > 1e-12) ? slack / row_norms[i] : DBL_MAX;
        if (r_i < *radius) *radius = r_i;
    }

    free(row_norms);

    if (*radius < -1e-10) {
        *radius = 0.0;
        return -1;  /* Polyhedron is empty */
    }
    if (*radius < 0.0) *radius = 0.0;

    return 0;
}

int pwa_polyhedron_intersect(const PWAPolyhedron *P, const PWAPolyhedron *Q)
{
    if (!P || !Q) return 0;
    if (P->dim != Q->dim) return 0;

    /* LP feasibility: find x s.t. H_P x ≤ K_P and H_Q x ≤ K_Q */
    /* Simplified: check Chebyshev center feasibility of combined system */
    PWAPolyhedron combined;
    combined.dim = P->dim;
    combined.n_halfspaces = P->n_halfspaces + Q->n_halfspaces;
    combined.is_bounded = 0;
    combined.is_empty = 0;

    int d = P->dim;
    int m = combined.n_halfspaces;

    combined.H = (double*)malloc((size_t)(m * d) * sizeof(double));
    combined.K = (double*)malloc((size_t)m * sizeof(double));
    if (!combined.H || !combined.K) {
        free(combined.H);
        free(combined.K);
        return 0;
    }

    /* Copy P constraints */
    memcpy(combined.H, P->H, (size_t)(P->n_halfspaces * d) * sizeof(double));
    memcpy(combined.K, P->K, (size_t)P->n_halfspaces * sizeof(double));

    /* Copy Q constraints */
    memcpy(combined.H + P->n_halfspaces * d, Q->H,
           (size_t)(Q->n_halfspaces * d) * sizeof(double));
    memcpy(combined.K + P->n_halfspaces, Q->K,
           (size_t)Q->n_halfspaces * sizeof(double));

    double *center = (double*)calloc((size_t)d, sizeof(double));
    double radius = 0.0;
    int result = 0;

    if (center) {
        int ret = pwa_polyhedron_chebyshev_center(&combined, center, &radius);
        result = (ret == 0 && radius >= -1e-10) ? 1 : 0;
        free(center);
    }

    free(combined.H);
    free(combined.K);

    return result;
}

int pwa_polyhedron_intersection(const PWAPolyhedron *P,
                                 const PWAPolyhedron *Q,
                                 PWAPolyhedron *R)
{
    if (!P || !Q || !R) return -1;
    if (P->dim != Q->dim) return -1;

    int m1 = P->n_halfspaces;
    int m2 = Q->n_halfspaces;
    int d = P->dim;
    int m_total = m1 + m2;

    R->dim = d;
    R->n_halfspaces = m_total;
    R->is_bounded = 0;
    R->is_empty = 0;

    R->H = (double*)malloc((size_t)(m_total * d) * sizeof(double));
    R->K = (double*)malloc((size_t)m_total * sizeof(double));
    if (!R->H || !R->K) {
        free(R->H);
        free(R->K);
        return -1;
    }

    memcpy(R->H, P->H, (size_t)(m1 * d) * sizeof(double));
    memcpy(R->K, P->K, (size_t)m1 * sizeof(double));
    memcpy(R->H + m1 * d, Q->H, (size_t)(m2 * d) * sizeof(double));
    memcpy(R->K + m1, Q->K, (size_t)m2 * sizeof(double));

    return 0;
}

int pwa_polyhedron_subset(const PWAPolyhedron *P, const PWAPolyhedron *Q)
{
    if (!P || !Q) return 0;
    if (P->dim != Q->dim) return 0;

    /* P ⊆ Q iff for every constraint h_j·x ≤ k_j of Q,
     * sup_{x∈P} h_j·x ≤ k_j.
     * Solve each as LP: max h_j^T x s.t. H_P x ≤ K_P */

    int d = P->dim;

    for (int j = 0; j < Q->n_halfspaces; j++) {
        const double *qj = &Q->H[j * d];
        double k_j = Q->K[j];

        /* Solve LP using simplified simplex (or just check vertices) */
        /* For this implementation, use a poor man's LP via
         * iterative projection: max c^T x s.t. H x ≤ K */
        double *x = (double*)calloc((size_t)d, sizeof(double));
        if (!x) return 0;

        /* Iterative gradient ascent with projection */
        double step = 0.01;
        for (int iter = 0; iter < 500; iter++) {
            /* Gradient step in direction qj */
            for (int i = 0; i < d; i++) {
                x[i] += step * qj[i];
            }

            /* Project back onto P */
            for (int proj = 0; proj < 20; proj++) {
                for (int ci = 0; ci < P->n_halfspaces; ci++) {
                    double val = 0.0;
                    const double *Hi = &P->H[ci * d];
                    for (int k = 0; k < d; k++) {
                        val += Hi[k] * x[k];
                    }
                    if (val > P->K[ci]) {
                        double ns = 0.0;
                        for (int k = 0; k < d; k++) ns += Hi[k] * Hi[k];
                        if (ns > 1e-12) {
                            double excess = (val - P->K[ci]) / ns;
                            for (int k = 0; k < d; k++) {
                                x[k] -= excess * Hi[k];
                            }
                        }
                    }
                }
            }
        }

        /* Check objective */
        double obj = 0.0;
        for (int i = 0; i < d; i++) obj += qj[i] * x[i];
        free(x);

        if (obj > k_j + 1e-6) {
            return 0;  /* Not a subset */
        }
    }

    return 1;
}

int pwa_polyhedron_bounding_box(const PWAPolyhedron *poly, PWABoundingBox *bb)
{
    if (!poly || !bb) return -1;

    int d = poly->dim;
    bb->dim = d;
    if (!bb->lb) bb->lb = (double*)malloc((size_t)d * sizeof(double));
    if (!bb->ub) bb->ub = (double*)malloc((size_t)d * sizeof(double));
    if (!bb->lb || !bb->ub) return -1;

    for (int dim = 0; dim < d; dim++) {
        /* Solve: min/max x_dim s.t. Hx ≤ K */
        /* For upper bound */
        double *c = (double*)calloc((size_t)d, sizeof(double));
        if (!c) return -1;
        c[dim] = 1.0;

        /* Simple iterative solver */
        double *x = (double*)calloc((size_t)d, sizeof(double));
        if (!x) { free(c); return -1; }

        /* Gradient ascent for max, descent for min */
        for (int sign = 0; sign < 2; sign++) {
            double step_dir = (sign == 0) ? 1.0 : -1.0; /* 0=max, 1=min */

            memset(x, 0, (size_t)d * sizeof(double));
            double lr = 0.01;
            for (int iter = 0; iter < 500; iter++) {
                x[dim] += lr * step_dir;

                /* Project onto polyhedron */
                for (int proj_iter = 0; proj_iter < 20; proj_iter++) {
                    for (int ci = 0; ci < poly->n_halfspaces; ci++) {
                        const double *Hi = &poly->H[ci * d];
                        double val = 0.0;
                        for (int k = 0; k < d; k++) val += Hi[k] * x[k];
                        if (val > poly->K[ci]) {
                            double ns = 0.0;
                            for (int k = 0; k < d; k++) ns += Hi[k] * Hi[k];
                            if (ns > 1e-12) {
                                double excess = (val - poly->K[ci]) / ns;
                                for (int k = 0; k < d; k++) x[k] -= excess * Hi[k];
                            }
                        }
                    }
                }
            }

            if (sign == 0) bb->ub[dim] = x[dim];  /* max */
            else bb->lb[dim] = x[dim];            /* min */
        }

        free(x);
        free(c);
    }

    return 0;
}

/*===========================================================================
 * L3: Convex Hull - Graham Scan (2D)
 *===========================================================================*/

static int compare_points_by_angle(const void *a, const void *b, void *ref)
{
    PWAPoint2D *p0 = (PWAPoint2D*)ref;
    const PWAPoint2D *p1 = (const PWAPoint2D*)a;
    const PWAPoint2D *p2 = (const PWAPoint2D*)b;

    double dx1 = p1->x - p0->x;
    double dy1 = p1->y - p0->y;
    double dx2 = p2->x - p0->x;
    double dy2 = p2->y - p0->y;

    double cross = dx1 * dy2 - dy1 * dx2;
    if (fabs(cross) < 1e-12) {
        /* Collinear: sort by distance */
        double d1 = dx1 * dx1 + dy1 * dy1;
        double d2 = dx2 * dx2 + dy2 * dy2;
        return (d1 < d2) ? -1 : (d1 > d2) ? 1 : 0;
    }
    return (cross > 0) ? -1 : 1;  /* CCW order */
}

static double cross_2d(PWAPoint2D O, PWAPoint2D A, PWAPoint2D B)
{
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

int pwa_convex_hull_2d(const PWAPoint2D *points, int n,
                        PWAConvexHull2D *hull)
{
    if (!points || n < 3 || !hull) return -1;

    /* Allocate workspace */
    PWAPoint2D *pts = (PWAPoint2D*)malloc((size_t)n * sizeof(PWAPoint2D));
    if (!pts) return -1;
    memcpy(pts, points, (size_t)n * sizeof(PWAPoint2D));

    /* Find lowest point (by y, then x) */
    int lowest = 0;
    for (int i = 1; i < n; i++) {
        if (pts[i].y < pts[lowest].y ||
            (fabs(pts[i].y - pts[lowest].y) < 1e-12 && pts[i].x < pts[lowest].x)) {
            lowest = i;
        }
    }

    /* Swap lowest to first */
    PWAPoint2D tmp = pts[0];
    pts[0] = pts[lowest];
    pts[lowest] = tmp;

    /* Sort by polar angle */
    /* Since qsort_r not universally available, use a simple insertion
     * sort based on compare function */
    for (int i = 1; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx1 = pts[i].x - pts[0].x;
            double dy1 = pts[i].y - pts[0].y;
            double dx2 = pts[j].x - pts[0].x;
            double dy2 = pts[j].y - pts[0].y;
            double cross = dx1 * dy2 - dy1 * dx2;

            int should_swap = 0;
            if (fabs(cross) < 1e-12) {
                double d1 = dx1 * dx1 + dy1 * dy1;
                double d2 = dx2 * dx2 + dy2 * dy2;
                if (d1 > d2) should_swap = 1;
            } else if (cross < 0) {
                should_swap = 1;
            }

            if (should_swap) {
                PWAPoint2D temp = pts[i];
                pts[i] = pts[j];
                pts[j] = temp;
            }
        }
    }

    /* Graham scan */
    int *stack = (int*)malloc((size_t)n * sizeof(int));
    if (!stack) { free(pts); return -1; }

    int top = 0;
    stack[top++] = 0;
    stack[top++] = 1;
    stack[top++] = 2;

    for (int i = 3; i < n; i++) {
        while (top >= 2 &&
               cross_2d(pts[stack[top-2]], pts[stack[top-1]], pts[i]) <= 1e-12) {
            top--;
        }
        stack[top++] = i;
    }

    /* Copy hull points */
    hull->n_points = top;
    if (!hull->points) {
        hull->points = (PWAPoint2D*)malloc((size_t)top * sizeof(PWAPoint2D));
    }
    if (hull->points) {
        for (int i = 0; i < top; i++) {
            hull->points[i] = pts[stack[i]];
        }
    }

    /* Compute area using shoelace formula */
    hull->area = 0.0;
    for (int i = 0; i < top; i++) {
        int j = (i + 1) % top;
        hull->area += hull->points[i].x * hull->points[j].y
                    - hull->points[j].x * hull->points[i].y;
    }
    hull->area = fabs(hull->area) * 0.5;

    free(stack);
    free(pts);
    return top;
}

/*===========================================================================
 * L3: Convex Hull - Gift Wrapping (arbitrary dimension)
 *===========================================================================*/

int pwa_convex_hull_giftwrap(const double *points, int n, int dim,
                              double *hull, int *n_hull)
{
    if (!points || n < 2 || dim < 2 || !hull || !n_hull) return -1;

    /* Gift wrapping (Jarvis march) generalized to d dimensions.
     * For simplicity, implement 2D case fully; for higher dims,
     * return vertices that are extreme in at least one direction. */

    *n_hull = 0;

    if (dim == 2) {
        /* 2D Jarvis march */
        /* Find leftmost point */
        int leftmost = 0;
        for (int i = 1; i < n; i++) {
            if (points[i * 2] < points[leftmost * 2]) leftmost = i;
        }

        int *hull_idx = (int*)malloc((size_t)n * sizeof(int));
        if (!hull_idx) return -1;

        int p = leftmost;
        int count = 0;
        do {
            hull_idx[count++] = p;
            int q = (p + 1) % n;

            for (int i = 0; i < n; i++) {
                if (i == p) continue;

                double px = points[p * 2];
                double py = points[p * 2 + 1];
                double qx = points[q * 2];
                double qy = points[q * 2 + 1];
                double ix = points[i * 2];
                double iy = points[i * 2 + 1];

                double cross = (qx - px) * (iy - py) - (qy - py) * (ix - px);
                if (cross < -1e-12) {
                    q = i;
                } else if (fabs(cross) < 1e-12) {
                    /* Collinear: take farthest */
                    double dq = (qx-px)*(qx-px) + (qy-py)*(qy-py);
                    double di = (ix-px)*(ix-px) + (iy-py)*(iy-py);
                    if (di > dq) q = i;
                }
            }

            p = q;
        } while (p != leftmost && count < n);

        for (int i = 0; i < count; i++) {
            hull[i * 2] = points[hull_idx[i] * 2];
            hull[i * 2 + 1] = points[hull_idx[i] * 2 + 1];
        }

        *n_hull = count;
        free(hull_idx);
    } else {
        /* Higher dimension: return extreme points.
         * A point p_i is extreme if there exists a direction c
         * such that c^T p_i > c^T p_j for all j ≠ i.
         *
         * Test 2*dim axis-aligned directions and 2*dim more
         * random directions. */

        int *is_extreme = (int*)calloc((size_t)n, sizeof(int));
        if (!is_extreme) return -1;

        int n_dirs = 4 * dim;
        for (int d = 0; d < n_dirs; d++) {
            double *dir = (double*)calloc((size_t)dim, sizeof(double));
            if (!dir) { free(is_extreme); return -1; }

            /* Generate direction */
            if (d < dim) {
                dir[d] = 1.0;
            } else if (d < 2 * dim) {
                dir[d - dim] = -1.0;
            } else {
                /* Pseudo-random direction */
                for (int i = 0; i < dim; i++) {
                    dir[i] = ((double)(d * 71 + i * 37) / 100.0);
                    dir[i] = dir[i] - floor(dir[i]);
                    dir[i] = dir[i] * 2.0 - 1.0;
                }
            }

            /* Find point with max dot product */
            double max_val = -DBL_MAX;
            int max_idx = -1;
            for (int i = 0; i < n; i++) {
                double val = 0.0;
                for (int j = 0; j < dim; j++) {
                    val += dir[j] * points[i * dim + j];
                }
                if (val > max_val) {
                    max_val = val;
                    max_idx = i;
                }
            }

            if (max_idx >= 0) is_extreme[max_idx] = 1;
            free(dir);
        }

        int count = 0;
        for (int i = 0; i < n; i++) {
            if (is_extreme[i]) {
                memcpy(hull + count * dim, points + i * dim,
                       (size_t)dim * sizeof(double));
                count++;
            }
        }

        *n_hull = count;
        free(is_extreme);
    }

    return 0;
}

/*===========================================================================
 * L3: Point in Convex Hull
 *===========================================================================*/

int pwa_point_in_convex_hull(const double *points, int n, int dim,
                              const double *x)
{
    if (!points || n < 1 || !x) return 0;

    if (dim == 1) {
        double lo = points[0], hi = points[0];
        for (int i = 1; i < n; i++) {
            if (points[i] < lo) lo = points[i];
            if (points[i] > hi) hi = points[i];
        }
        return (x[0] >= lo - 1e-10 && x[0] <= hi + 1e-10) ? 1 : 0;
    }

    /* For general dimension, use Carathéodory's theorem:
     * x ∈ CH(P) iff x can be expressed as convex combination of
     * at most dim+1 points from P.
     *
     * Solve LP: find λ ≥ 0, Σλ = 1, Σ λ_i p_i = x.
     *
     * For robustness, use iterative projection onto the convex hull:
     * start at centroid, move toward x, projecting onto CH(P). */

    double *y = (double*)calloc((size_t)dim, sizeof(double));
    if (!y) return 0;

    /* Start at centroid */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dim; j++) {
            y[j] += points[i * dim + j] / (double)n;
        }
    }

    /* Move toward x, projecting onto CH(P) at each step */
    for (int iter = 0; iter < 1000; iter++) {
        /* Step toward x */
        double dist = 0.0;
        for (int j = 0; j < dim; j++) {
            double diff = x[j] - y[j];
            dist += diff * diff;
        }
        if (dist < 1e-12) break;

        double step = 0.1;
        for (int j = 0; j < dim; j++) {
            y[j] += step * (x[j] - y[j]);
        }

        /* Project onto CH(P): find closest point in convex hull */
        /* Simplified: check if y is inside; if not, project onto
         * the closest facet or vertex. For here, just check if y
         * is a convex combination.
         *
         * Use a simple barycentric check: is y further from origin
         * than the farthest vertex in its direction? */
        double y_norm_sq = 0.0;
        for (int j = 0; j < dim; j++) y_norm_sq += y[j] * y[j];
        if (y_norm_sq < 1e-12) continue;

        /* Project y onto the line from origin to y, constrained by CH */
        double max_proj = -DBL_MAX;
        for (int i = 0; i < n; i++) {
            double dot = 0.0;
            for (int j = 0; j < dim; j++) {
                dot += y[j] * points[i * dim + j];
            }
            if (dot > max_proj) max_proj = dot;
        }

        if (y_norm_sq > max_proj + 1e-10) {
            /* Scale back */
            double scale = max_proj / y_norm_sq;
            for (int j = 0; j < dim; j++) {
                y[j] *= scale;
            }
        }
    }

    /* Check convergence */
    double final_dist = 0.0;
    for (int j = 0; j < dim; j++) {
        double diff = x[j] - y[j];
        final_dist += diff * diff;
    }

    free(y);
    return (sqrt(final_dist) < 1e-6) ? 1 : 0;
}

/*===========================================================================
 * L3: Simplex LP (Simplified Implementation)
 *===========================================================================*/

int pwa_simplex_lp(const double *A, const double *b, const double *c,
                    int m, int n, double *x, double *opt_val)
{
    if (!A || !b || !c || !x || m <= 0 || n <= 0) return -1;

    /* For small LPs, use an iterative gradient projection method.
     * This is not a full simplex implementation but provides reasonable
     * solutions for the polyhedral operations needed in PWA analysis. */

    /* Initialize x to zero */
    memset(x, 0, (size_t)n * sizeof(double));

    double lr = 0.01;
    for (int iter = 0; iter < 2000; iter++) {
        /* Gradient step: minimize c^T x → move in direction -c */
        for (int j = 0; j < n; j++) {
            x[j] -= lr * c[j];
        }

        /* Enforce non-negativity */
        for (int j = 0; j < n; j++) {
            if (x[j] < 0.0) x[j] = 0.0;
        }

        /* Project onto Ax ≤ b */
        for (int proj = 0; proj < 10; proj++) {
            for (int i = 0; i < m; i++) {
                double val = 0.0;
                const double *Ai = &A[i * n];
                for (int j = 0; j < n; j++) {
                    val += Ai[j] * x[j];
                }
                if (val > b[i] + 1e-10) {
                    double ns = 0.0;
                    for (int j = 0; j < n; j++) ns += Ai[j] * Ai[j];
                    if (ns > 1e-12) {
                        double excess = (val - b[i]) / ns;
                        for (int j = 0; j < n; j++) {
                            x[j] -= excess * Ai[j];
                        }
                    }
                }
            }
        }

        /* Decrease learning rate */
        lr *= 0.999;
    }

    /* Compute objective value */
    *opt_val = 0.0;
    for (int j = 0; j < n; j++) {
        *opt_val += c[j] * x[j];
    }

    /* Check feasibility */
    for (int i = 0; i < m; i++) {
        double val = 0.0;
        const double *Ai = &A[i * n];
        for (int j = 0; j < n; j++) val += Ai[j] * x[j];
        if (val > b[i] + 1e-8) return -1;  /* Infeasible */
    }

    return 0;
}

/*===========================================================================
 * L3: Fourier-Motzkin Elimination
 *===========================================================================*/

int pwa_fourier_motzkin_eliminate(const double *H_in, const double *K_in,
                                   int m, int d,
                                   double *H_out, double *K_out, int *m_out)
{
    if (!H_in || !K_in || m < 1 || d < 2 || !H_out || !K_out || !m_out) {
        return -1;
    }

    /* Eliminate the last variable x_d.
     *
     * For each constraint i: h_{i,1} x_1 + ... + h_{i,d} x_d ≤ k_i
     *
     * Partition constraints by sign of h_{i,d}:
     *   P = {i | h_{i,d} > 0}  →  x_d ≤ (k_i - Σ_{j<d} h_{i,j} x_j) / h_{i,d}
     *   N = {i | h_{i,d} < 0}  →  x_d ≥ (k_i - Σ_{j<d} h_{i,j} x_j) / h_{i,d}
     *   Z = {i | h_{i,d} == 0} → no x_d
     *
     * Then for each pair (p ∈ P, n ∈ N):
     *   (K_n - Σ_{j<d} H_{n,j} x_j)/H_{n,d} ≤ x_d ≤ (K_p - Σ_{j<d} H_{p,j} x_j)/H_{p,d}
     *   → Σ_{j<d} (H_{p,j}/H_{p,d} - H_{n,j}/H_{n,d}) x_j ≤ K_p/H_{p,d} - K_n/H_{n,d}
     */

    int d1 = d - 1;  /* Reduced dimension */
    int max_out = (m * m) / 4 + m;  /* Maximum new constraints */

    int n_pos = 0, n_neg = 0, n_zero = 0;
    int *pos_idx = (int*)malloc((size_t)m * sizeof(int));
    int *neg_idx = (int*)malloc((size_t)m * sizeof(int));
    int *zero_idx = (int*)malloc((size_t)m * sizeof(int));
    if (!pos_idx || !neg_idx || !zero_idx) {
        free(pos_idx); free(neg_idx); free(zero_idx);
        return -1;
    }

    for (int i = 0; i < m; i++) {
        double coeff = H_in[i * d + d1];  /* h_{i,d} */
        if (coeff > 1e-12) {
            pos_idx[n_pos++] = i;
        } else if (coeff < -1e-12) {
            neg_idx[n_neg++] = i;
        } else {
            zero_idx[n_zero++] = i;
        }
    }

    int count = 0;

    /* Pass through zero-coefficient constraints (those without x_d) */
    for (int z = 0; z < n_zero; z++) {
        int i = zero_idx[z];
        if (count < max_out) {
            memcpy(&H_out[count * d1], &H_in[i * d], (size_t)d1 * sizeof(double));
            K_out[count] = K_in[i];
            count++;
        }
    }

    /* Combine positive and negative constraints */
    for (int pi = 0; pi < n_pos; pi++) {
        int p = pos_idx[pi];
        double hp_d = H_in[p * d + d1];

        for (int ni = 0; ni < n_neg; ni++) {
            int n = neg_idx[ni];
            double hn_d = H_in[n * d + d1];

            if (count >= max_out) break;

            /* New constraint: for j=0..d-2:
             *   (H[p,j]/hp_d - H[n,j]/hn_d) x_j ≤ K[p]/hp_d - K[n]/hn_d */
            for (int j = 0; j < d1; j++) {
                H_out[count * d1 + j] = H_in[p * d + j] / hp_d
                                       - H_in[n * d + j] / hn_d;
            }
            K_out[count] = K_in[p] / hp_d - K_in[n] / hn_d;
            count++;
        }
    }

    /* If no positive or no negative constraints, just pass through zeros;
     * if all zero, projection is vacuous (unconstrained in x_d) */
    if (n_pos == 0 && n_neg == 0) {
        /* All are zero constraints; already handled above.
         * Nothing more to add; the projection is just the zero constraints
         * in the reduced dimension. */
    } else if (n_pos == 0 && n_zero == 0) {
        /* Only negative: x_d ≥ all lower bounds → unconstrained above */
        /* Pass through as-is (no elimination needed for existence) */
    } else if (n_neg == 0 && n_zero == 0) {
        /* Only positive: x_d ≤ all upper bounds → unconstrained below */
        /* Pass through as-is */
    }

    *m_out = count;

    free(pos_idx);
    free(neg_idx);
    free(zero_idx);

    return 0;
}

/*===========================================================================
 * L3: Support Function
 *===========================================================================*/

int pwa_polyhedron_support(const PWAPolyhedron *poly, const double *c,
                            double *val, double *xsup)
{
    if (!poly || !c || !val) return -1;

    int d = poly->dim;
    int m = poly->n_halfspaces;

    /* Solve LP: max c^T x s.t. H x ≤ K */
    /* Use gradient ascent with projection */
    double *x = (double*)calloc((size_t)d, sizeof(double));
    if (!x) return -1;

    /* Start from Chebyshev center if possible, otherwise from origin */
    double *center = (double*)calloc((size_t)d, sizeof(double));
    double radius;
    if (center && pwa_polyhedron_chebyshev_center(poly, center, &radius) == 0) {
        memcpy(x, center, (size_t)d * sizeof(double));
    }
    free(center);

    double lr = 0.01;
    for (int iter = 0; iter < 1000; iter++) {
        /* Gradient step in direction c */
        for (int j = 0; j < d; j++) {
            x[j] += lr * c[j];
        }

        /* Project onto polyhedron */
        for (int proj = 0; proj < 15; proj++) {
            for (int i = 0; i < m; i++) {
                double val_i = 0.0;
                const double *Hi = &poly->H[i * d];
                for (int j = 0; j < d; j++) val_i += Hi[j] * x[j];
                if (val_i > poly->K[i] + 1e-10) {
                    double ns = 0.0;
                    for (int j = 0; j < d; j++) ns += Hi[j] * Hi[j];
                    if (ns > 1e-12) {
                        double excess = (val_i - poly->K[i]) / ns;
                        for (int j = 0; j < d; j++) x[j] -= excess * Hi[j];
                    }
                }
            }
        }
    }

    *val = 0.0;
    for (int j = 0; j < d; j++) *val += c[j] * x[j];

    if (xsup) {
        memcpy(xsup, x, (size_t)d * sizeof(double));
    }

    free(x);
    return 0;
}

/*===========================================================================
 * L3: Minkowski Sum and Hausdorff Distance
 *===========================================================================*/

int pwa_minkowski_sum(const PWAPolyhedron *P, const PWAPolyhedron *Q,
                       const double *dirs, int n_dirs, PWAPolyhedron *R)
{
    if (!P || !Q || !dirs || !R || n_dirs < 1) return -1;
    if (P->dim != Q->dim) return -1;

    int d = P->dim;
    R->dim = d;
    R->n_halfspaces = n_dirs;
    R->is_bounded = 0;
    R->is_empty = 0;

    R->H = (double*)malloc((size_t)(n_dirs * d) * sizeof(double));
    R->K = (double*)malloc((size_t)n_dirs * sizeof(double));
    if (!R->H || !R->K) {
        free(R->H); free(R->K);
        return -1;
    }

    for (int i = 0; i < n_dirs; i++) {
        const double *dir = &dirs[i * d];

        /* Normalize direction */
        double norm = 0.0;
        for (int j = 0; j < d; j++) norm += dir[j] * dir[j];
        norm = sqrt(norm);
        if (norm < 1e-12) norm = 1.0;

        for (int j = 0; j < d; j++) {
            R->H[i * d + j] = dir[j] / norm;
        }

        /* h_{P⊕Q}(c) = h_P(c) + h_Q(c) */
        double val_P, val_Q;
        pwa_polyhedron_support(P, R->H + i * d, &val_P, NULL);
        pwa_polyhedron_support(Q, R->H + i * d, &val_Q, NULL);

        R->K[i] = val_P + val_Q;
    }

    return 0;
}

double pwa_hausdorff_distance(const PWAPolyhedron *P, const PWAPolyhedron *Q,
                               const double *dirs, int n_dirs)
{
    if (!P || !Q || !dirs || n_dirs < 1) return 0.0;

    int d = P->dim;
    double max_diff = 0.0;

    for (int i = 0; i < n_dirs; i++) {
        const double *dir = &dirs[i * d];

        double norm = 0.0;
        for (int j = 0; j < d; j++) norm += dir[j] * dir[j];
        norm = sqrt(norm);
        if (norm < 1e-12) continue;

        /* Normalize */
        double *dir_n = (double*)malloc((size_t)d * sizeof(double));
        if (!dir_n) continue;
        for (int j = 0; j < d; j++) dir_n[j] = dir[j] / norm;

        double val_P, val_Q;
        pwa_polyhedron_support(P, dir_n, &val_P, NULL);
        pwa_polyhedron_support(Q, dir_n, &val_Q, NULL);

        double diff = fabs(val_P - val_Q);
        if (diff > max_diff) max_diff = diff;

        /* Also check opposite direction */
        for (int j = 0; j < d; j++) dir_n[j] = -dir[j] / norm;
        pwa_polyhedron_support(P, dir_n, &val_P, NULL);
        pwa_polyhedron_support(Q, dir_n, &val_Q, NULL);

        diff = fabs(val_P - val_Q);
        if (diff > max_diff) max_diff = diff;

        free(dir_n);
    }

    return max_diff;
}

/*===========================================================================
 * L3: Half-Space Intersection Vertices (BFS enumeration)
 *===========================================================================*/

int pwa_halfspace_intersection_vertices(const PWAHalfSpace *halfspaces,
                                         int m, int d,
                                         double *vertices, int *n_vert)
{
    if (!halfspaces || m < d || !vertices || !n_vert) return -1;
    *n_vert = 0;

    /* Enumerate all subsets of d constraints, solve the linear system,
     * and check if the solution satisfies all other constraints. */

    if (d > 6 || m > 20) {
        /* Too many combinations for brute force; use extreme direction
         * sampling as an approximation */
        int n_dirs = 2 * d + 10;
        for (int di = 0; di < n_dirs; di++) {
            double *dir = (double*)calloc((size_t)d, sizeof(double));
            if (!dir) continue;

            if (di < d) dir[di] = 1.0;
            else if (di < 2 * d) dir[di - d] = -1.0;
            else {
                /* Random direction */
                unsigned seed = (unsigned)(di * 2654435761U);
                for (int j = 0; j < d; j++) {
                    seed = seed * 1103515245U + 12345U;
                    dir[j] = ((double)(seed & 0x7FFFFFFF) / 2147483648.0) * 2.0 - 1.0;
                }
            }

            /* Solve LP: max dir·x s.t. Hx ≤ K using gradient ascent */
            double *x = (double*)calloc((size_t)d, sizeof(double));
            if (x) {
                double lr = 0.01;
                for (int iter = 0; iter < 500; iter++) {
                    for (int j = 0; j < d; j++) x[j] += lr * dir[j];
                    for (int proj = 0; proj < 10; proj++) {
                        for (int i = 0; i < m; i++) {
                            double val = -halfspaces[i].b;
                            for (int jj = 0; jj < d; jj++)
                                val += halfspaces[i].a[jj] * x[jj];
                            if (val > 1e-10) {
                                double ns = 0.0;
                                for (int jj = 0; jj < d; jj++)
                                    ns += halfspaces[i].a[jj] * halfspaces[i].a[jj];
                                if (ns > 1e-12) {
                                    double excess = val / ns;
                                    for (int jj = 0; jj < d; jj++)
                                        x[jj] -= excess * halfspaces[i].a[jj];
                                }
                            }
                        }
                    }
                }

                /* Check if point satisfies all constraints */
                int feasible = 1;
                for (int i = 0; i < m; i++) {
                    double val = -halfspaces[i].b;
                    for (int jj = 0; jj < d; jj++) val += halfspaces[i].a[jj] * x[jj];
                    if (val > 1e-8) { feasible = 0; break; }
                }

                if (feasible) {
                    /* Check if this vertex is new */
                    int is_new = 1;
                    for (int v = 0; v < *n_vert; v++) {
                        double dist = 0.0;
                        for (int jj = 0; jj < d; jj++) {
                            double diff = x[jj] - vertices[v * d + jj];
                            dist += diff * diff;
                        }
                        if (dist < 1e-10) { is_new = 0; break; }
                    }
                    if (is_new && *n_vert < m * 2) {
                        memcpy(&vertices[(*n_vert) * d], x, (size_t)d * sizeof(double));
                        (*n_vert)++;
                    }
                }

                free(x);
            }
            free(dir);
        }
        return 0;
    }

    /* Brute force for small cases: try all combinations of d constraints */
    /* Generate all d-combinations of m constraints */
    int *combo = (int*)malloc((size_t)d * sizeof(int));
    if (!combo) return -1;

    for (int i = 0; i < d; i++) combo[i] = i;

    while (1) {
        /* Solve the d×d linear system for the intersection of the d hyperplanes:
         * a_{combo[k]} · x = b_{combo[k]} for k=0..d-1 */
        double *A = (double*)malloc((size_t)(d * d) * sizeof(double));
        double *rhs = (double*)malloc((size_t)d * sizeof(double));
        if (A && rhs) {
            for (int k = 0; k < d; k++) {
                memcpy(&A[k * d], halfspaces[combo[k]].a, (size_t)d * sizeof(double));
                rhs[k] = halfspaces[combo[k]].b;
            }

            /* Gaussian elimination with partial pivoting */
            double *sol = (double*)malloc((size_t)d * sizeof(double));
            if (sol) {
                /* Copy for solving */
                double *AA = (double*)malloc((size_t)(d * d) * sizeof(double));
                double *bb = (double*)malloc((size_t)d * sizeof(double));
                if (AA && bb) {
                    memcpy(AA, A, (size_t)(d * d) * sizeof(double));
                    memcpy(bb, rhs, (size_t)d * sizeof(double));

                    int singular = 0;
                    for (int col = 0; col < d; col++) {
                        /* Find pivot */
                        int pivot = col;
                        double max_val = fabs(AA[col * d + col]);
                        for (int row = col + 1; row < d; row++) {
                            if (fabs(AA[row * d + col]) > max_val) {
                                max_val = fabs(AA[row * d + col]);
                                pivot = row;
                            }
                        }
                        if (max_val < 1e-12) { singular = 1; break; }

                        /* Swap rows */
                        if (pivot != col) {
                            for (int jj = 0; jj < d; jj++) {
                                double tmp = AA[col * d + jj];
                                AA[col * d + jj] = AA[pivot * d + jj];
                                AA[pivot * d + jj] = tmp;
                            }
                            double tmp = bb[col];
                            bb[col] = bb[pivot];
                            bb[pivot] = tmp;
                        }

                        /* Eliminate */
                        for (int row = col + 1; row < d; row++) {
                            double factor = AA[row * d + col] / AA[col * d + col];
                            for (int jj = col; jj < d; jj++) {
                                AA[row * d + jj] -= factor * AA[col * d + jj];
                            }
                            bb[row] -= factor * bb[col];
                        }
                    }

                    if (!singular) {
                        /* Back substitution */
                        for (int row = d - 1; row >= 0; row--) {
                            double sum = bb[row];
                            for (int jj = row + 1; jj < d; jj++) {
                                sum -= AA[row * d + jj] * sol[jj];
                            }
                            sol[row] = sum / AA[row * d + row];
                        }

                        /* Check if this solution satisfies ALL m constraints */
                        int feasible = 1;
                        for (int i = 0; i < m; i++) {
                            double val = -halfspaces[i].b;
                            for (int jj = 0; jj < d; jj++) {
                                val += halfspaces[i].a[jj] * sol[jj];
                            }
                            if (val > 1e-8) { feasible = 0; break; }
                        }

                        if (feasible) {
                            /* Check if this vertex is already found */
                            int is_new = 1;
                            for (int v = 0; v < *n_vert; v++) {
                                double dist = 0.0;
                                for (int jj = 0; jj < d; jj++) {
                                    double diff = sol[jj] - vertices[v * d + jj];
                                    dist += diff * diff;
                                }
                                if (dist < 1e-10) { is_new = 0; break; }
                            }
                            if (is_new) {
                                memcpy(&vertices[(*n_vert) * d], sol,
                                       (size_t)d * sizeof(double));
                                (*n_vert)++;
                            }
                        }
                    }
                    free(AA);
                    free(bb);
                }
                free(sol);
            }
            free(A);
            free(rhs);
        }

        /* Generate next combination */
        int k;
        for (k = d - 1; k >= 0; k--) {
            if (combo[k] < m - d + k) break;
        }
        if (k < 0) break;

        combo[k]++;
        for (int j = k + 1; j < d; j++) {
            combo[j] = combo[j - 1] + 1;
        }
    }

    free(combo);
    return 0;
}

/*===========================================================================
 * L3: Voronoi Partition for PWA Regions
 *===========================================================================*/

int pwa_voronoi_partition(const double *sites, int n_sites, int dim,
                           PWAPolyhedron *regions, int *n_regions)
{
    if (!sites || n_sites < 2 || dim < 1 || !regions || !n_regions) return -1;

    /* For each site i, the Voronoi cell is:
     *   V_i = { x | ||x - s_i||^2 ≤ ||x - s_j||^2  for all j ≠ i }
     *
     * This simplifies to: 2(s_j - s_i)^T x ≤ ||s_j||^2 - ||s_i||^2
     * which is a half-space.
     *
     * So each V_i is the intersection of (n_sites - 1) half-spaces.
     */

    *n_regions = n_sites;

    for (int i = 0; i < n_sites; i++) {
        regions[i].dim = dim;
        regions[i].n_halfspaces = n_sites - 1;
        regions[i].is_bounded = 0;
        regions[i].is_empty = 0;

        regions[i].H = (double*)malloc((size_t)((n_sites - 1) * dim) * sizeof(double));
        regions[i].K = (double*)malloc((size_t)(n_sites - 1) * sizeof(double));
        if (!regions[i].H || !regions[i].K) continue;

        double si_norm_sq = 0.0;
        for (int d = 0; d < dim; d++) {
            si_norm_sq += sites[i * dim + d] * sites[i * dim + d];
        }

        int cons_idx = 0;
        for (int j = 0; j < n_sites; j++) {
            if (j == i) continue;

            double sj_norm_sq = 0.0;
            for (int d = 0; d < dim; d++) {
                sj_norm_sq += sites[j * dim + d] * sites[j * dim + d];
            }

            /* 2(s_j - s_i)^T x ≤ ||s_j||^2 - ||s_i||^2 */
            for (int d = 0; d < dim; d++) {
                regions[i].H[cons_idx * dim + d] =
                    2.0 * (sites[j * dim + d] - sites[i * dim + d]);
            }
            regions[i].K[cons_idx] = sj_norm_sq - si_norm_sq;
            cons_idx++;
        }
    }

    return 0;
}
