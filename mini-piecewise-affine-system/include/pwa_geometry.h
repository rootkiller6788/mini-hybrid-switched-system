/**
 * @file pwa_geometry.h
 * @brief Polyhedral Geometry Operations — L3 Mathematical Structures
 *
 * Implements computational geometry primitives for polyhedral
 * regions used in piecewise affine systems. A polyhedron P ⊂ R^d
 * is the intersection of finitely many closed half-spaces:
 *
 *   P = { x ∈ R^d | H x ≤ K }
 *
 * where H ∈ R^{m×d} and K ∈ R^m.
 *
 * A polytope is a bounded polyhedron. All regions in a PWA system
 * must be polytopes (bounded) to ensure well-posedness.
 *
 * References:
 *   Ziegler, G. M. (1995). "Lectures on Polytopes." Springer GTM 152.
 *   Grünbaum, B. (2003). "Convex Polytopes." Springer GTM 221, 2nd ed.
 *   Fukuda, K. (2004). "Frequently Asked Questions in Polyhedral Computation."
 *   Preparata, F. P. & Shamos, M. I. (1985). "Computational Geometry."
 *     Springer-Verlag. Chapters 3, 7.
 *
 * Knowledge coverage:
 *   L3: Polyhedron, half-space, vertex, face, convex hull,
 *       hyperplane, facet enumeration, Minkowski sum
 *
 * Nine-school course alignment:
 *   MIT 6.838 Computational Geometry — Lec 3 Convex Hulls
 *   Stanford CS368 Computational Geometry — Lec 4 Polyhedra
 *   Berkeley CS274 Computational Geometry — Ch 5 Convex Polytopes
 *   CMU 15-451 Algorithms — Lec 17 Convex Hull
 *   ETH 263-5210 Computational Geometry — Lec 6 Polyhedra
 */

#ifndef PWA_GEOMETRY_H
#define PWA_GEOMETRY_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
 * L3: Core Geometry Types
 *===========================================================================*/

/**
 * @brief A half-space in R^d: { x ∈ R^d | a·x ≤ b }
 *
 * The vector a is the outward normal, b is the offset.
 * A polyhedron is the intersection of m half-spaces.
 */
typedef struct {
    int     dim;       /**< Ambient dimension d */
    double *a;         /**< Normal vector, length d */
    double  b;         /**< Right-hand side offset */
} PWAHalfSpace;

/**
 * @brief A hyperplane in R^d: { x ∈ R^d | a·x = b }
 *
 * Divides R^d into two half-spaces: a·x ≤ b and a·x ≥ b.
 */
typedef struct {
    int     dim;       /**< Ambient dimension d */
    double *a;         /**< Normal vector, length d */
    double  b;         /**< Right-hand side offset */
} PWAHyperplane;

/**
 * @brief A convex polyhedron as intersection of half-spaces (H-representation).
 *
 * P = { x ∈ R^d | H x ≤ K } where H ∈ R^{m×d}, K ∈ R^m.
 */
typedef struct {
    int      dim;          /**< Ambient dimension d */
    int      n_halfspaces; /**< Number of half-space constraints m */
    double  *H;            /**< Constraint matrix, m × d, row-major */
    double  *K;            /**< Constraint RHS, length m */
    int      is_bounded;   /**< 1 if polyhedron is bounded (polytope) */
    int      is_empty;     /**< 1 if polyhedron is empty */
} PWAPolyhedron;

/**
 * @brief A vertex of a polytope in R^d.
 */
typedef struct {
    int     dim;       /**< Ambient dimension d */
    double *x;         /**< Coordinates, length d */
    int     is_extreme; /**< 1 if vertex is extreme point of polytope */
} PWAVertex;

/**
 * @brief A face of a polytope.
 *
 * A k-face is a k-dimensional face. Vertices are 0-faces,
 * edges are 1-faces, facets are (d-1)-faces.
 */
typedef struct {
    int      dim;          /**< Dimension of the face k */
    int      n_vertices;   /**< Number of vertices defining this face */
    int     *vertex_ids;   /**< Indices of vertices in this face */
    double  *normal;       /**< Outward normal (for facets, length d) */
    double   offset;       /**< Offset for facet hyperplane */
} PWAFace;

/**
 * @brief A full polytope with both V-representation and H-representation.
 *
 * V-representation: convex hull of vertices
 * H-representation: intersection of half-spaces
 */
typedef struct {
    int          dim;           /**< Ambient dimension d */
    int          n_vertices;    /**< Number of vertices */
    PWAVertex   *vertices;      /**< Vertex array */
    int          n_facets;      /**< Number of facets */
    PWAFace     *facets;        /**< Facet array */
    double      *centroid;      /**< Centroid of the polytope (length d) */
    double       volume;        /**< Volume of the polytope */
    PWAPolyhedron h_rep;        /**< Half-space representation */
} PWAPolytope;

/**
 * @brief Point in 2D for convex hull computations.
 */
typedef struct {
    double x;
    double y;
} PWAPoint2D;

/**
 * @brief Convex hull in 2D.
 */
typedef struct {
    int          n_points;      /**< Number of hull points */
    PWAPoint2D  *points;        /**< Hull points in CCW order */
    double       area;          /**< Hull area */
} PWAConvexHull2D;

/**
 * @brief Bounding box (axis-aligned) in R^d.
 */
typedef struct {
    int     dim;       /**< Dimension d */
    double *lb;        /**< Lower bound, length d */
    double *ub;        /**< Upper bound, length d */
} PWABoundingBox;

/*===========================================================================
 * L3: Half-Space Operations
 *===========================================================================*/

/**
 * @brief Evaluate half-space constraint: computes a·x - b.
 *
 * @param hs Half-space
 * @param x  Point in R^d
 * @return a·x - b (≤ 0 means inside, = 0 means on boundary)
 *
 * Complexity: O(d)
 */
double pwa_halfspace_eval(const PWAHalfSpace *hs, const double *x);

/**
 * @brief Signed distance from a point to a hyperplane.
 *
 * d(x, H) = (a·x - b) / ||a||
 *
 * @param hp Hyperplane
 * @param x  Point in R^d
 * @return Signed distance (positive = in direction of normal)
 *
 * Complexity: O(d)
 */
double pwa_hyperplane_distance(const PWAHyperplane *hp, const double *x);

/**
 * @brief Check if a point satisfies all polyhedron constraints.
 *
 * @param poly Polyhedron
 * @param x    Point in R^d
 * @return 1 if x ∈ poly (with tolerance 1e-10), 0 otherwise
 *
 * Complexity: O(m * d)
 */
int pwa_polyhedron_contains(const PWAPolyhedron *poly, const double *x);

/**
 * @brief Compute the Chebyshev center of a polyhedron.
 *
 * The Chebyshev center is the center of the largest inscribed ball.
 * Solved via linear programming: max r s.t. H_i·x + r·||H_i|| ≤ K_i
 *
 * @param poly    Polyhedron
 * @param center  Output: Chebyshev center (length d)
 * @param radius  Output: radius of largest inscribed ball
 * @return 0 on success, -1 if polyhedron is empty
 *
 * Complexity: O(m * d) iterations (simplex-like LP)
 */
int pwa_polyhedron_chebyshev_center(const PWAPolyhedron *poly,
                                     double *center, double *radius);

/**
 * @brief Check if two polyhedra intersect.
 *
 * Uses linear programming feasibility check:
 *   ∃ x s.t. H_1 x ≤ K_1 ∧ H_2 x ≤ K_2
 *
 * @param P First polyhedron
 * @param Q Second polyhedron
 * @return 1 if intersection is non-empty, 0 otherwise
 *
 * Complexity: O((m1+m2) * d) iterations
 */
int pwa_polyhedron_intersect(const PWAPolyhedron *P, const PWAPolyhedron *Q);

/**
 * @brief Intersection of two polyhedra.
 *
 * Result: R = P ∩ Q = { x | [H_P; H_Q] x ≤ [K_P; K_Q] }
 *
 * @param P    First polyhedron
 * @param Q    Second polyhedron
 * @param R    Output: P ∩ Q (caller must allocate)
 * @return 0 on success, -1 on error
 */
int pwa_polyhedron_intersection(const PWAPolyhedron *P,
                                 const PWAPolyhedron *Q,
                                 PWAPolyhedron *R);

/**
 * @brief Check if one polyhedron is a subset of another.
 *
 * P ⊆ Q iff every constraint of Q is redundant for the system (P, Q_constraints).
 * Solved by checking: sup_{x∈P} h_j·x ≤ k_j for each constraint j of Q.
 *
 * @param P Inner polyhedron
 * @param Q Outer polyhedron
 * @return 1 if P ⊆ Q, 0 otherwise
 */
int pwa_polyhedron_subset(const PWAPolyhedron *P, const PWAPolyhedron *Q);

/**
 * @brief Compute the bounding box of a polyhedron.
 *
 * For each dimension i, solves:
 *   lb_i = min{x_i | Hx ≤ K}
 *   ub_i = max{x_i | Hx ≤ K}
 *
 * @param poly Polyhedron
 * @param bb   Output: bounding box (caller must allocate)
 * @return 0 on success, -1 if polyhedron is empty
 */
int pwa_polyhedron_bounding_box(const PWAPolyhedron *poly, PWABoundingBox *bb);

/*===========================================================================
 * L3: Convex Hull Operations
 *===========================================================================*/

/**
 * @brief Compute the convex hull of a set of 2D points (Graham scan).
 *
 * @param points      Input points
 * @param n           Number of input points
 * @param hull        Output: convex hull (caller must allocate hull->points
 *                    with at least n entries)
 * @return Number of hull points (CCW order), -1 on error
 *
 * Theorem: Graham scan computes CH(P) in O(n log n) time.
 * Reference: Graham, R. L. (1972). "An Efficient Algorithm for
 *   Determining the Convex Hull of a Finite Planar Set."
 *   Info. Proc. Lett., 1:132-133.
 */
int pwa_convex_hull_2d(const PWAPoint2D *points, int n,
                        PWAConvexHull2D *hull);

/**
 * @brief Compute convex hull using the Gift Wrapping (Jarvis march) algorithm.
 *
 * Works in any dimension but practical for 2D/3D.
 * For d=2: O(n·h) where h is number of hull points.
 *
 * @param points    Input points
 * @param n         Number of points
 * @param dim       Ambient dimension
 * @param hull      Output: hull points (length n)
 * @param n_hull    Output: number of hull vertices
 * @return 0 on success, -1 on error
 *
 * Reference: Jarvis, R. A. (1973). "On the identification of the convex
 *   hull of a finite set of points in the plane." Info. Proc. Lett., 2:18-21.
 */
int pwa_convex_hull_giftwrap(const double *points, int n, int dim,
                              double *hull, int *n_hull);

/**
 * @brief Check if a point is inside the convex hull of a set of points.
 *
 * Uses linear programming: point ∈ CH(P) iff ∃ λ ≥ 0, Σ λ_i = 1, x = Σ λ_i p_i
 *
 * @param points    Point set (n × dim, row-major)
 * @param n         Number of points
 * @param dim       Ambient dimension
 * @param x         Query point (length dim)
 * @return 1 if x ∈ CH(points), 0 otherwise
 */
int pwa_point_in_convex_hull(const double *points, int n, int dim,
                              const double *x);

/*===========================================================================
 * L3: Linear Programming (Simplex) for Polyhedral Operations
 *===========================================================================*/

/**
 * @brief Simplex method for LP: min c^T x s.t. A x ≤ b, x ≥ 0
 *
 * Standard form: min c^T x s.t. A x ≤ b, x ≥ 0
 *
 * @param A       Constraint matrix (m × n)
 * @param b       Constraint RHS (m)
 * @param c       Cost vector (n)
 * @param m       Number of constraints
 * @param n       Number of variables
 * @param x       Output: optimal solution (n)
 * @param opt_val Output: optimal value
 * @return 0=optimal, 1=unbounded, -1=infeasible
 *
 * Complexity: worst-case exponential, average O(m^3 n)
 *
 * Theorem (Dantzig 1947): The simplex method terminates at an optimal
 * basic feasible solution or detects unboundedness/infeasibility.
 */
int pwa_simplex_lp(const double *A, const double *b, const double *c,
                    int m, int n, double *x, double *opt_val);

/**
 * @brief Fourier-Motzkin elimination for polyhedral projection.
 *
 * Eliminates one variable from a polyhedron Hx ≤ K, producing
 * the projection onto the remaining coordinates.
 *
 * Given P = { (x_1, ..., x_n) | H·[x_1;...;x_n] ≤ K },
 * computes proj_{x_1,...,x_{n-1}}(P).
 *
 * @param H_in       Input constraint matrix (m × d)
 * @param K_in       Input constraint RHS (m)
 * @param m          Number of constraints
 * @param d          Dimension
 * @param H_out      Output constraint matrix (max m^2/4 × (d-1))
 * @param K_out      Output constraint RHS
 * @param m_out      Output: number of resulting constraints
 * @return 0 on success, -1 on memory error
 *
 * Complexity: O(m^2) constraints produced per elimination step.
 * Reference: Fourier (1826), Motzkin (1936).
 */
int pwa_fourier_motzkin_eliminate(const double *H_in, const double *K_in,
                                   int m, int d,
                                   double *H_out, double *K_out, int *m_out);

/**
 * @brief Compute the support function of a polyhedron.
 *
 * h_P(c) = sup{ c^T x | H x ≤ K }
 *
 * Solved via LP.
 *
 * @param poly Polyhedron
 * @param c    Direction vector (length d)
 * @param val  Output: support function value
 * @param xsup Output: supporting point (length d, may be NULL)
 * @return 0 on success, -1 if unbounded or empty
 *
 * Complexity: O(m * d) per iteration (LP)
 */
int pwa_polyhedron_support(const PWAPolyhedron *poly, const double *c,
                            double *val, double *xsup);

/*===========================================================================
 * L3: Minkowski Operations
 *===========================================================================*/

/**
 * @brief Minkowski sum of two polytopes (approximate via support functions).
 *
 * P ⊕ Q = { p + q | p ∈ P, q ∈ Q }
 *
 * For polytopes: h_{P⊕Q}(c) = h_P(c) + h_Q(c)
 *
 * @param P    First polytope
 * @param Q    Second polytope
 * @param dirs Direction vectors for outer approximation
 * @param n_dirs Number of directions
 * @param R    Output: approximate Minkowski sum (H-representation)
 * @return 0 on success, -1 on error
 */
int pwa_minkowski_sum(const PWAPolyhedron *P, const PWAPolyhedron *Q,
                       const double *dirs, int n_dirs, PWAPolyhedron *R);

/**
 * @brief Hausdorff distance between two polytopes (approximate).
 *
 * d_H(P, Q) = max{ sup_{p∈P} inf_{q∈Q} ||p-q||,
 *                   sup_{q∈Q} inf_{p∈P} ||p-q|| }
 *
 * @param P    First polytope
 * @param Q    Second polytope
 * @param dirs Directions for approximation
 * @param n_dirs Number of directions
 * @return Approximate Hausdorff distance
 */
double pwa_hausdorff_distance(const PWAPolyhedron *P, const PWAPolyhedron *Q,
                               const double *dirs, int n_dirs);

/*===========================================================================
 * L3: Half-space intersection (dual of convex hull)
 *===========================================================================*/

/**
 * @brief Compute the intersection of half-spaces and extract extreme points.
 *
 * Given m half-spaces a_i·x ≤ b_i, compute all vertices of the resulting
 * polytope by enumerating all (d choose d) vertex candidates from solving
 * subsystems of d equations a_i·x = b_i.
 *
 * @param halfspaces Array of m half-spaces
 * @param m          Number of half-spaces
 * @param d          Dimension
 * @param vertices   Output: array of extreme points
 * @param n_vert     Output: number of vertices found
 * @return 0 on success, -1 on error
 *
 * Complexity: O(m^d * d^3) — the McMullen upper bound theorem limits
 * the number of vertices to O(m^{⌊d/2⌋}).
 */
int pwa_halfspace_intersection_vertices(const PWAHalfSpace *halfspaces,
                                         int m, int d,
                                         double *vertices, int *n_vert);

/**
 * @brief Delaunay-like triangulation of a 2D point set into polyhedral regions.
 *
 * Partitions the convex hull of a point set into regions where each
 * region is closer to its site than to any other (Voronoi-like for PWA).
 *
 * @param sites     Site points (n_sites × dim)
 * @param n_sites   Number of sites
 * @param dim       Dimension
 * @param regions   Output: partition regions as H-representation
 * @param n_regions Output: number of regions
 * @return 0 on success, -1 on error
 */
int pwa_voronoi_partition(const double *sites, int n_sites, int dim,
                           PWAPolyhedron *regions, int *n_regions);

#ifdef __cplusplus
}
#endif

#endif /* PWA_GEOMETRY_H */
