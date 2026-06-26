/**
 * @file pwa_identification.h
 * @brief Data-Driven PWA System Identification — L5 Algorithms
 *
 * Identifies piecewise affine models from input-output data.
 * The identification problem is: given data {(x_k, u_k, y_k, x_{k+1})},
 * find:
 *   1. A partition of the regressor space into polyhedral regions
 *   2. Affine dynamics parameters (A_i, B_i, f_i) for each region
 *
 * The approach combines clustering (to find regions in the data),
 * classification (to find separating hyperplanes), and least-squares
 * regression (to estimate affine parameters within each region).
 *
 * This is a mixed-integer/continuous estimation problem. Common methods:
 *   - Clustering-based: cluster data, fit hyperplanes between clusters,
 *     then fit affine models within clusters
 *   - Bounded-error: use set-membership identification
 *   - Mixed-integer programming: MINLP formulation
 *   - Bayesian: EM-like alternating between assignment and estimation
 *
 * References:
 *   Ferrari-Trecate, G., Muselli, M., Liberati, D., & Morari, M. (2003).
 *     "A clustering technique for the identification of piecewise affine
 *     systems." Automatica, 39(2):205-217.
 *   Bemporad, A., Garulli, A., Paoletti, S., & Vicino, A. (2005).
 *     "A bounded-error approach to piecewise affine system identification."
 *     IEEE TAC, 50(10):1567-1580.
 *   Roll, J., Bemporad, A., & Ljung, L. (2004).
 *     "Identification of piecewise affine systems via mixed-integer
 *     programming." Automatica, 40(1):37-50.
 *
 * Knowledge coverage:
 *   L5: K-means clustering, SVM-like hyperplane separation,
 *       least-squares regression, model validation
 *
 * Nine-school course alignment:
 *   MIT 6.241J — Lec 22 System identification
 *   Stanford AA203 — Lec 18 Learning for control
 *   Berkeley EE221A — Lec 24 Identification
 *   CMU 24-677 — Lec 22 Data-driven methods
 *   ETH 227-0216 — Lec 20-24 System identification
 */

#ifndef PWA_IDENTIFICATION_H
#define PWA_IDENTIFICATION_H

#include "pwa_defs.h"
#include "pwa_geometry.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
 * L5: Identification Types
 *===========================================================================*/

/**
 * @brief Data point for PWA identification.
 *
 * A regressor is z_k = [x_k; u_k] and the output is either
 * x_{k+1} - x_k (CT derivative approximation) or x_{k+1} (DT).
 */
typedef struct {
    int      index;       /**< Data point index */
    int      n_state;     /**< State dimension */
    int      n_input;     /**< Input dimension */
    double  *z;           /**< Regressor [x; u] of length n_state + n_input */
    double  *y;           /**< Output (x_{k+1} or dx/dt) of length n_state */
    int      cluster_id;  /**< Assigned cluster (-1 if unassigned) */
    double   weight;      /**< Data point weight */
} PWADataPoint;

/**
 * @brief Dataset for PWA identification.
 */
typedef struct {
    int          n_points;      /**< Number of data points */
    int          n_state;       /**< State dimension */
    int          n_input;       /**< Input dimension */
    int          n_features;    /**< Feature dimension = n_state + n_input */
    PWADataPoint *points;       /**< Data points array */
    double      *stats_mean;    /**< Feature means */
    double      *stats_std;     /**< Feature standard deviations */
} PWADataset;

/**
 * @brief Cluster in the regressor space.
 */
typedef struct {
    int          id;            /**< Cluster ID */
    int          n_points;      /**< Number of points in cluster */
    int         *point_ids;     /**< Indices of points in this cluster */
    double      *centroid;      /**< Cluster centroid (n_features) */
    double      *covariance;    /**< Cluster covariance (n_features × n_features) */
    double       within_ss;     /**< Within-cluster sum of squares */
} PWACluster;

/**
 * @brief Clustering result from K-means or similar.
 */
typedef struct {
    int          n_clusters;    /**< Number of clusters */
    int          n_features;    /**< Feature dimension */
    PWACluster  *clusters;      /**< Cluster array */
    double       total_wss;     /**< Total within-cluster sum of squares */
    double       silhouette;    /**< Average silhouette score */
} PWAClustering;

/**
 * @brief Regression result for a single region.
 *
 * Solves: Y = X * Θ^T  where Θ = [A | B | f]
 * via least squares: Θ = (X^T X)^{-1} X^T Y
 */
typedef struct {
    int          n_points;      /**< Number of points used */
    double      *A;             /**< Estimated A matrix (n_state × n_state) */
    double      *B;             /**< Estimated B matrix (n_state × n_input) */
    double      *f;             /**< Estimated f vector (n_state) */
    double       rss;           /**< Residual sum of squares */
    double       r_squared;     /**< R² goodness of fit */
    double      *residuals;     /**< Residuals per data point (n_points × n_state) */
} PWARegression;

/**
 * @brief PWA identification result.
 */
typedef struct {
    int          n_regions;     /**< Number of identified regions */
    int          n_state;       /**< State dimension */
    int          n_input;       /**< Input dimension */
    PWAClustering clustering;   /**< Clustering result */
    double      **H;            /**< Region boundaries: H[i] * z ≤ K[i] */
    double      **K;            /**< Region boundary RHS */
    PWARegression *regressions; /**< Regression for each region */
    PWASystem   *pwa_system;    /**< Resulting PWA system */
    double       total_error;   /**< Total fitting error */
} PWAIdentification;

/*===========================================================================
 * L5: Dataset Operations
 *===========================================================================*/

/**
 * @brief Create a PWA identification dataset.
 *
 * @param n_points Maximum number of data points
 * @param n_state  State dimension
 * @param n_input  Input dimension
 * @return New dataset, or NULL on error
 */
PWADataset* pwa_dataset_create(int n_points, int n_state, int n_input);

/**
 * @brief Destroy a dataset and free all memory.
 *
 * @param ds Dataset to destroy
 */
void pwa_dataset_destroy(PWADataset *ds);

/**
 * @brief Add a data point to the dataset.
 *
 * @param ds Dataset
 * @param z  Regressor vector [x; u] (n_state + n_input)
 * @param y  Output vector x_{k+1} (n_state)
 * @param weight Data point weight
 * @return Data point index, -1 on error
 */
int pwa_dataset_add_point(PWADataset *ds,
                           const double *z, const double *y, double weight);

/**
 * @brief Normalize dataset to zero mean, unit variance.
 *
 * Stores normalization parameters in ds->stats_mean, ds->stats_std.
 *
 * @param ds Dataset to normalize (modified in place)
 * @return 0 on success, -1 on error
 */
int pwa_dataset_normalize(PWADataset *ds);

/**
 * @brief Generate synthetic PWA dataset from a known PWA system.
 *
 * Simulates the system with random initial conditions and inputs
 * to generate identification data.
 *
 * @param sys       Source PWA system
 * @param n_traj    Number of trajectories to generate
 * @param steps_per Number of steps per trajectory
 * @param noise_std Standard deviation of measurement noise
 * @param ds        Output: dataset (must be pre-allocated)
 * @return 0 on success, -1 on error
 */
int pwa_dataset_generate(const PWASystem *sys, int n_traj, int steps_per,
                          double noise_std, PWADataset *ds);

/*===========================================================================
 * L5: K-Means Clustering
 *===========================================================================*/

/**
 * @brief K-means clustering of regression data.
 *
 * Minimizes the within-cluster sum of squares:
 *   WSS = Σ_i Σ_{z∈C_i} ||z - μ_i||^2
 *
 * Standard Lloyd's algorithm with random initialization.
 * Repeats n_restarts with different initial centroids and keeps
 * the best result.
 *
 * @param ds          Dataset
 * @param k           Number of clusters
 * @param max_iters   Maximum iterations per run
 * @param n_restarts  Number of random restarts
 * @param tol         Convergence tolerance
 * @param clustering  Output: clustering result
 * @return 0 on success, -1 on error
 *
 * Complexity: O(n_restarts · max_iters · k · n_points · n_features)
 *
 * Reference: Lloyd, S. P. (1982). "Least squares quantization in PCM."
 *   IEEE Trans. Info. Theory, 28(2):129-137.
 *   (Originally published as Bell Labs Tech. Note, 1957.)
 */
int pwa_kmeans_cluster(const PWADataset *ds, int k, int max_iters,
                        int n_restarts, double tol, PWAClustering *clustering);

/**
 * @brief Run K-means with K-means++ initialization.
 *
 * K-means++ selects initial centroids with probability proportional
 * to squared distance from nearest existing centroid, giving
 * O(log k) approximation guarantee in expectation.
 *
 * @param ds          Dataset
 * @param k           Number of clusters
 * @param max_iters   Maximum iterations
 * @param tol         Convergence tolerance
 * @param clustering  Output: clustering result
 * @return 0 on success, -1 on error
 *
 * Reference: Arthur, D. & Vassilvitskii, S. (2007).
 *   "K-means++: The Advantages of Careful Seeding." SODA 2007.
 */
int pwa_kmeans_plusplus(const PWADataset *ds, int k, int max_iters,
                         double tol, PWAClustering *clustering);

/**
 * @brief Compute the silhouette score for a clustering.
 *
 * For each point i: s(i) = (b(i) - a(i)) / max{a(i), b(i)}
 * where a(i) = mean distance to own cluster,
 *       b(i) = min mean distance to other clusters.
 *
 * @param ds          Dataset
 * @param clustering  Clustering to evaluate
 * @return Average silhouette score in [-1, 1]
 */
double pwa_silhouette_score(const PWADataset *ds, const PWAClustering *clustering);

/*===========================================================================
 * L5: Region Boundary Estimation
 *===========================================================================*/

/**
 * @brief Estimate separating hyperplanes between clusters using SVM-like approach.
 *
 * For each pair of clusters (i, j), finds a hyperplane w·z + b = 0
 * that separates the points with maximum margin.
 *
 * Uses a simplified perceptron-like linear classifier.
 *
 * @param ds          Dataset
 * @param clustering  Clustering result
 * @param H           Output: constraint matrices (n_clusters × max_cons × n_features)
 * @param K           Output: constraint RHS (n_clusters × max_cons)
 * @param n_cons_arr  Output: number of constraints per cluster
 * @return 0 on success, -1 on error
 *
 * Complexity: O(k^2 · n_points · n_features)
 */
int pwa_estimate_region_boundaries(const PWADataset *ds,
                                    const PWAClustering *clustering,
                                    double **H, double **K, int *n_cons_arr);

/**
 * @brief Assign a polyhedral region to each cluster.
 *
 * Region i is defined by the intersection of half-spaces
 * that separate cluster i from all other clusters j ≠ i:
 *
 *   R_i = { z | w_{ij}·z + b_{ij} ≤ 0 for all j ≠ i }
 *
 * @param clustering  Clustering
 * @param n_features  Feature dimension
 * @param H           Full constraint matrices per region
 * @param K           Full constraint RHS per region
 * @param n_cons      Number of constraints per region
 * @return 0 on success, -1 on error
 */
int pwa_build_polyhedral_regions(const PWAClustering *clustering,
                                  int n_features,
                                  double **H, double **K, int *n_cons);

/**
 * @brief Ensure region partition covers the dataset domain.
 *
 * Checks for gaps and overlaps, applies refinement if needed.
 *
 * @param regions    Region polyhedra
 * @param n_regions  Number of regions
 * @param H          Constraint matrices
 * @param K          Constraint RHS
 * @param n_cons     Constraints per region
 * @return 0 on success, -1 on error
 */
int pwa_refine_partition(int n_regions, double **H, double **K, int *n_cons);

/*===========================================================================
 * L5: Affine Regression per Region
 *===========================================================================*/

/**
 * @brief Fit affine model to data points in a region using least squares.
 *
 * For DT-PWA: x_{k+1} = A x_k + B u_k + f
 * Solves: min_{A,B,f} Σ ||x_{k+1} - A x_k - B u_k - f||^2
 *
 * Formulation: Y = Θ · [X; 1] where Θ = [A B f]
 * Solution: Θ = Y · [X; 1]^T · ([X; 1] · [X; 1]^T)^{-1}
 *
 * @param ds            Dataset
 * @param point_ids     Indices of points in region
 * @param n_region_pts  Number of points in region
 * @param n_state       State dimension
 * @param n_input       Input dimension
 * @param reg           Output: regression result
 * @return 0 on success, -1 on error
 *
 * Complexity: O(n_region_pts · (n_state + n_input)^2 + (n_state + n_input)^3)
 */
int pwa_fit_affine_model(const PWADataset *ds,
                          const int *point_ids, int n_region_pts,
                          int n_state, int n_input, PWARegression *reg);

/*===========================================================================
 * L5: Full Identification Pipeline
 *===========================================================================*/

/**
 * @brief Run the full PWA identification pipeline.
 *
 * Steps:
 *   1. Cluster data using K-means++
 *   2. Estimate region boundaries (separating hyperplanes)
 *   3. Build polyhedral regions
 *   4. Fit affine models per region via least squares
 *   5. Assemble PWA system from identified components
 *
 * @param ds            Dataset (normalized recommended)
 * @param n_regions     Desired number of regions
 * @param max_iters     Max K-means iterations
 * @param ident         Output: identification result
 * @return 0 on success, -1 on error
 *
 * Reference: Ferrari-Trecate et al. (2003), Algorithm 1.
 */
int pwa_identify(const PWADataset *ds, int n_regions, int max_iters,
                  PWAIdentification *ident);

/**
 * @brief Validate identified PWA model against test data.
 *
 * Computes:
 *   - RMS prediction error
 *   - Maximum absolute error
 *   - Variance accounted for (VAF): 1 - var(e) / var(y)
 *
 * @param ident     Identification result
 * @param test_ds   Test dataset
 * @param rms_error Output: root mean square error
 * @param vaf       Output: variance accounted for [0, 1]
 * @return 0 on success, -1 on error
 */
int pwa_validate_model(const PWAIdentification *ident,
                        const PWADataset *test_ds,
                        double *rms_error, double *vaf);

/**
 * @brief Determine the optimal number of regions using gap statistic.
 *
 * Tests k = k_min .. k_max and selects the k with the largest
 * gap between within-cluster dispersion and its expectation
 * under a null reference distribution.
 *
 * @param ds       Dataset
 * @param k_min    Minimum k to test
 * @param k_max    Maximum k to test
 * @param n_boot   Number of bootstrap reference samples
 * @param best_k   Output: optimal number of regions
 * @return 0 on success, -1 on error
 *
 * Reference: Tibshirani, R., Walther, G., & Hastie, T. (2001).
 *   "Estimating the number of clusters in a data set via the gap statistic."
 *   JRSS B, 63(2):411-423.
 */
int pwa_select_num_regions(const PWADataset *ds, int k_min, int k_max,
                            int n_boot, int *best_k);

/**
 * @brief Destroy a PWAIdentification result.
 *
 * @param ident Identification to destroy
 */
void pwa_identification_destroy(PWAIdentification *ident);

#ifdef __cplusplus
}
#endif

#endif /* PWA_IDENTIFICATION_H */
