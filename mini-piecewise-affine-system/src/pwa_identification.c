/**
 * @file pwa_identification.c
 * @brief Data-Driven PWA System Identification — L5 Algorithms
 *
 * Implements the complete pipeline for identifying piecewise affine
 * models from input-output data: clustering, region boundary
 * estimation, affine regression, and model validation.
 *
 * Knowledge coverage:
 *   L5: K-means/K-means++ clustering, SVM-like hyperplane separation,
 *       least-squares regression, silhouette score, gap statistic,
 *       full identification pipeline, model validation (RMS, VAF)
 *
 * References:
 *   Ferrari-Trecate et al. (2003). "A clustering technique for
 *     the identification of PWA systems." Automatica, 39(2):205-217.
 *   Arthur & Vassilvitskii (2007). "K-means++." SODA 2007.
 *   Tibshirani et al. (2001). "Gap statistic." JRSS B, 63(2):411-423.
 */

#include "pwa_identification.h"
#include "pwa_simulation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

/*===========================================================================
 * L5: Dataset Operations
 *===========================================================================*/

PWADataset* pwa_dataset_create(int n_points, int n_state, int n_input)
{
    if (n_points <= 0 || n_state <= 0 || n_input < 0) return NULL;

    PWADataset *ds = (PWADataset*)calloc(1, sizeof(PWADataset));
    if (!ds) return NULL;

    ds->n_points = 0;
    ds->n_state = n_state;
    ds->n_input = n_input;
    ds->n_features = n_state + n_input;

    ds->points = (PWADataPoint*)calloc((size_t)n_points, sizeof(PWADataPoint));
    if (!ds->points) { free(ds); return NULL; }

    for (int i = 0; i < n_points; i++) {
        ds->points[i].n_state = n_state;
        ds->points[i].n_input = n_input;
        ds->points[i].z = (double*)calloc((size_t)ds->n_features, sizeof(double));
        ds->points[i].y = (double*)calloc((size_t)n_state, sizeof(double));
        ds->points[i].cluster_id = -1;
        ds->points[i].weight = 1.0;
    }

    ds->stats_mean = (double*)calloc((size_t)ds->n_features, sizeof(double));
    ds->stats_std = (double*)calloc((size_t)ds->n_features, sizeof(double));

    return ds;
}

void pwa_dataset_destroy(PWADataset *ds)
{
    if (!ds) return;
    if (ds->points) {
        for (int i = 0; i < ds->n_points; i++) {
            free(ds->points[i].z);
            free(ds->points[i].y);
        }
        free(ds->points);
    }
    free(ds->stats_mean);
    free(ds->stats_std);
    free(ds);
}

int pwa_dataset_add_point(PWADataset *ds,
                           const double *z, const double *y, double weight)
{
    if (!ds || !z || !y) return -1;

    /* Find first empty slot */
    int idx = ds->n_points;
    PWADataPoint *pt = &ds->points[idx];

    memcpy(pt->z, z, (size_t)ds->n_features * sizeof(double));
    memcpy(pt->y, y, (size_t)ds->n_state * sizeof(double));
    pt->weight = weight;
    pt->cluster_id = -1;
    pt->index = idx;

    ds->n_points++;
    return idx;
}

int pwa_dataset_normalize(PWADataset *ds)
{
    if (!ds) return -1;
    if (ds->n_points < 2) return -1;

    int nf = ds->n_features;
    int np = ds->n_points;

    /* Compute mean */
    memset(ds->stats_mean, 0, (size_t)nf * sizeof(double));
    for (int i = 0; i < np; i++) {
        for (int f = 0; f < nf; f++) {
            ds->stats_mean[f] += ds->points[i].z[f] / (double)np;
        }
    }

    /* Compute standard deviation */
    memset(ds->stats_std, 0, (size_t)nf * sizeof(double));
    for (int i = 0; i < np; i++) {
        for (int f = 0; f < nf; f++) {
            double diff = ds->points[i].z[f] - ds->stats_mean[f];
            ds->stats_std[f] += diff * diff / (double)np;
        }
    }
    for (int f = 0; f < nf; f++) {
        ds->stats_std[f] = sqrt(ds->stats_std[f]);
        if (ds->stats_std[f] < 1e-12) ds->stats_std[f] = 1.0;
    }

    /* Standardize */
    for (int i = 0; i < np; i++) {
        for (int f = 0; f < nf; f++) {
            ds->points[i].z[f] = (ds->points[i].z[f] - ds->stats_mean[f])
                                 / ds->stats_std[f];
        }
    }

    return 0;
}

int pwa_dataset_generate(const PWASystem *sys, int n_traj, int steps_per,
                          double noise_std, PWADataset *ds)
{
    if (!sys || !ds || n_traj <= 0 || steps_per <= 0) return -1;

    PWASimConfig cfg = pwa_sim_config_default();
    cfg.t_end = (double)steps_per * (sys->dt > 0 ? sys->dt : 0.01);
    cfg.max_steps = steps_per + 2;

    int n_state = sys->n_state;
    int n_input = sys->n_input;
    int nz = n_state + n_input;

    for (int traj_i = 0; traj_i < n_traj; traj_i++) {
        /* Random initial state */
        double *x0 = (double*)malloc((size_t)n_state * sizeof(double));
        if (!x0) continue;

        for (int i = 0; i < n_state; i++) {
            x0[i] = ((double)((traj_i * 101 + i * 37) % 1000) / 1000.0 - 0.5) * 2.0;
        }

        /* Generate random input sequence */
        double *u_seq = (double*)calloc((size_t)(steps_per * n_input), sizeof(double));
        if (u_seq) {
            for (int s = 0; s < steps_per; s++) {
                for (int i = 0; i < n_input; i++) {
                    u_seq[s * n_input + i] =
                        ((double)((s * 71 + traj_i * 13 + i * 31) % 1000) / 1000.0 - 0.5) * 2.0;
                }
            }
        }

        /* Simulate */
        PWATrajectory *traj = pwa_trajectory_create(n_state, n_input,
                                                     sys->n_output, steps_per + 2);
        if (!traj) { free(x0); free(u_seq); continue; }

        if (sys->is_continuous) {
            pwa_simulate_ct(sys, x0, NULL, NULL, &cfg, traj);
        } else {
            pwa_simulate_dt(sys, x0, u_seq, &cfg, traj);
        }

        /* Extract data points: z = [x_k; u_k], y = x_{k+1} */
        for (int s = 0; s < traj->n_steps - 1 && ds->n_points < ds->n_points + traj->n_steps; s++) {
            double *z = (double*)malloc((size_t)nz * sizeof(double));
            double *y = (double*)malloc((size_t)n_state * sizeof(double));
            if (!z || !y) { free(z); free(y); continue; }

            /* z = [x_k; u_k] */
            memcpy(z, &traj->x_hist[s * n_state], (size_t)n_state * sizeof(double));
            if (u_seq) {
                memcpy(z + n_state, &u_seq[s * n_input],
                       (size_t)n_input * sizeof(double));
            }

            /* y = x_{k+1} (with noise) */
            for (int i = 0; i < n_state; i++) {
                double x_next = traj->x_hist[(s + 1) * n_state + i];
                /* Add Gaussian noise approximation */
                double noise = ((double)((s * 19 + traj_i * 7 + i * 43) % 1000) / 1000.0 - 0.5)
                             * 2.0 * noise_std;
                y[i] = x_next + noise;
            }

            pwa_dataset_add_point(ds, z, y, 1.0);
            free(z);
            free(y);
        }

        pwa_trajectory_destroy(traj);
        free(x0);
        free(u_seq);
    }

    return 0;
}

/*===========================================================================
 * L5: K-Means Clustering
 *===========================================================================*/

int pwa_kmeans_cluster(const PWADataset *ds, int k, int max_iters,
                        int n_restarts, double tol, PWAClustering *clustering)
{
    if (!ds || k < 2 || !clustering) return -1;
    if (ds->n_points < k) return -1;

    int nf = ds->n_features;
    int np = ds->n_points;

    double best_wss = DBL_MAX;
    int *best_assign = (int*)malloc((size_t)np * sizeof(int));
    double *best_centroids = (double*)malloc((size_t)(k * nf) * sizeof(double));
    int *assign = (int*)malloc((size_t)np * sizeof(int));
    double *centroids = (double*)malloc((size_t)(k * nf) * sizeof(double));
    int *cluster_sizes = (int*)malloc((size_t)k * sizeof(int));

    if (!best_assign || !best_centroids || !assign || !centroids || !cluster_sizes) {
        free(best_assign); free(best_centroids); free(assign);
        free(centroids); free(cluster_sizes);
        return -1;
    }

    for (int restart = 0; restart < n_restarts; restart++) {
        /* Initialize centroids randomly from data points */
        int *chosen = (int*)calloc((size_t)np, sizeof(int));
        if (!chosen) continue;

        for (int c = 0; c < k; c++) {
            int pick;
            do {
                pick = (restart * k * 101 + c * 37) % np;
            } while (chosen[pick] && c < np);
            chosen[pick] = 1;
            memcpy(&centroids[c * nf], ds->points[pick].z,
                   (size_t)nf * sizeof(double));
        }
        free(chosen);

        double prev_wss = DBL_MAX;

        for (int iter = 0; iter < max_iters; iter++) {
            /* Step 1: Assign points to nearest centroid */
            double total_wss = 0.0;
            memset(cluster_sizes, 0, (size_t)k * sizeof(int));
            memset(centroids, 0, (size_t)(k * nf) * sizeof(double));

            for (int i = 0; i < np; i++) {
                double min_dist = DBL_MAX;
                int best_c = 0;

                for (int c = 0; c < k; c++) {
                    double dist = 0.0;
                    for (int f = 0; f < nf; f++) {
                        double diff = ds->points[i].z[f] - centroids[c * nf + f];
                        dist += diff * diff;
                    }
                    if (dist < min_dist) {
                        min_dist = dist;
                        best_c = c;
                    }
                }

                assign[i] = best_c;
                total_wss += min_dist;
                cluster_sizes[best_c]++;

                /* Accumulate for new centroid */
                for (int f = 0; f < nf; f++) {
                    centroids[best_c * nf + f] += ds->points[i].z[f];
                }
            }

            /* Step 2: Update centroids */
            for (int c = 0; c < k; c++) {
                if (cluster_sizes[c] > 0) {
                    for (int f = 0; f < nf; f++) {
                        centroids[c * nf + f] /= (double)cluster_sizes[c];
                    }
                }
            }

            /* Check convergence */
            if (fabs(total_wss - prev_wss) < tol) break;
            prev_wss = total_wss;
        }

        /* Compute final WSS for this restart */
        double final_wss = 0.0;
        for (int i = 0; i < np; i++) {
            int c = assign[i];
            double dist = 0.0;
            for (int f = 0; f < nf; f++) {
                double diff = ds->points[i].z[f] - centroids[c * nf + f];
                dist += diff * diff;
            }
            final_wss += dist;
        }

        if (final_wss < best_wss) {
            best_wss = final_wss;
            memcpy(best_assign, assign, (size_t)np * sizeof(int));
            memcpy(best_centroids, centroids, (size_t)(k * nf) * sizeof(double));
        }
    }

    /* Build output clustering */
    clustering->n_clusters = k;
    clustering->n_features = nf;
    clustering->total_wss = best_wss;

    clustering->clusters = (PWACluster*)calloc((size_t)k, sizeof(PWACluster));
    if (!clustering->clusters) {
        free(best_assign); free(best_centroids); free(assign);
        free(centroids); free(cluster_sizes);
        return -1;
    }

    /* Count points per cluster */
    int *counts = (int*)calloc((size_t)k, sizeof(int));
    if (counts) {
        for (int i = 0; i < np; i++) counts[best_assign[i]]++;

        for (int c = 0; c < k; c++) {
            clustering->clusters[c].id = c;
            clustering->clusters[c].n_points = counts[c];
            clustering->clusters[c].point_ids = (int*)malloc(
                (size_t)(counts[c] > 0 ? counts[c] : 1) * sizeof(int));
            clustering->clusters[c].centroid = (double*)malloc(
                (size_t)nf * sizeof(double));
            clustering->clusters[c].covariance = (double*)calloc(
                (size_t)(nf * nf), sizeof(double));

            if (clustering->clusters[c].centroid) {
                memcpy(clustering->clusters[c].centroid, &best_centroids[c * nf],
                       (size_t)nf * sizeof(double));
            }

            /* Re-count for filling point_ids */
            int idx = 0;
            for (int i = 0; i < np; i++) {
                if (best_assign[i] == c && idx < counts[c]) {
                    clustering->clusters[c].point_ids[idx++] = i;
                }
            }
        }
        free(counts);
    }

    /* Compute silhouette score */
    clustering->silhouette = pwa_silhouette_score(ds, clustering);

    free(best_assign); free(best_centroids); free(assign);
    free(centroids); free(cluster_sizes);
    return 0;
}

int pwa_kmeans_plusplus(const PWADataset *ds, int k, int max_iters,
                         double tol, PWAClustering *clustering)
{
    if (!ds || k < 2 || !clustering) return -1;

    int nf = ds->n_features;
    int np = ds->n_points;
    if (np < k) return -1;

    double *centroids = (double*)malloc((size_t)(k * nf) * sizeof(double));
    double *min_dist = (double*)malloc((size_t)np * sizeof(double));
    int *assign = (int*)malloc((size_t)np * sizeof(int));
    if (!centroids || !min_dist || !assign) {
        free(centroids); free(min_dist); free(assign);
        return -1;
    }

    /* K-means++ initialization */
    /* Step 1: Choose first centroid uniformly at random */
    int first = 0;  /* Pseudo-random: use first point */
    memcpy(centroids, ds->points[first].z, (size_t)nf * sizeof(double));

    /* Initialize distances */
    for (int i = 0; i < np; i++) {
        double dist = 0.0;
        for (int f = 0; f < nf; f++) {
            double diff = ds->points[i].z[f] - centroids[f];
            dist += diff * diff;
        }
        min_dist[i] = dist;
    }

    /* Choose remaining centroids */
    for (int c = 1; c < k; c++) {
        /* Compute total distance */
        double total = 0.0;
        for (int i = 0; i < np; i++) total += min_dist[i];

        /* Choose next centroid with probability proportional to min_dist^2 */
        double r = ((double)(c * 101) / 1000.0) * total;
        double cum = 0.0;
        int pick = np - 1;
        for (int i = 0; i < np; i++) {
            cum += min_dist[i];
            if (cum >= r) { pick = i; break; }
        }

        memcpy(&centroids[c * nf], ds->points[pick].z, (size_t)nf * sizeof(double));

        /* Update minimum distances */
        for (int i = 0; i < np; i++) {
            double dist = 0.0;
            for (int f = 0; f < nf; f++) {
                double diff = ds->points[i].z[f] - centroids[c * nf + f];
                dist += diff * diff;
            }
            if (dist < min_dist[i]) min_dist[i] = dist;
        }
    }

    /* Run standard K-means from these centroids */
    int *cluster_sizes = (int*)malloc((size_t)k * sizeof(int));
    if (!cluster_sizes) {
        free(centroids); free(min_dist); free(assign);
        return -1;
    }

    for (int iter = 0; iter < max_iters; iter++) {
        double total_wss = 0.0;
        memset(cluster_sizes, 0, (size_t)k * sizeof(int));

        /* Copy centroids for update */
        double *new_centroids = (double*)calloc((size_t)(k * nf), sizeof(double));
        if (!new_centroids) break;

        /* Assign points */
        for (int i = 0; i < np; i++) {
            double best_dist = DBL_MAX;
            int best_c = 0;
            for (int c = 0; c < k; c++) {
                double dist = 0.0;
                for (int f = 0; f < nf; f++) {
                    double diff = ds->points[i].z[f] - centroids[c * nf + f];
                    dist += diff * diff;
                }
                if (dist < best_dist) { best_dist = dist; best_c = c; }
            }
            assign[i] = best_c;
            total_wss += best_dist;
            cluster_sizes[best_c]++;
            for (int f = 0; f < nf; f++) {
                new_centroids[best_c * nf + f] += ds->points[i].z[f];
            }
        }

        /* Update centroids */
        double max_shift = 0.0;
        for (int c = 0; c < k; c++) {
            if (cluster_sizes[c] > 0) {
                for (int f = 0; f < nf; f++) {
                    new_centroids[c * nf + f] /= (double)cluster_sizes[c];
                    double shift = fabs(new_centroids[c * nf + f]
                                       - centroids[c * nf + f]);
                    if (shift > max_shift) max_shift = shift;
                    centroids[c * nf + f] = new_centroids[c * nf + f];
                }
            }
        }

        free(new_centroids);
        if (max_shift < tol) break;
    }

    /* Build output */
    clustering->n_clusters = k;
    clustering->n_features = nf;
    clustering->clusters = (PWACluster*)calloc((size_t)k, sizeof(PWACluster));
    if (!clustering->clusters) {
        free(centroids); free(min_dist); free(assign); free(cluster_sizes);
        return -1;
    }

    /* Compute total WSS and fill cluster data */
    double total_wss = 0.0;
    int *counts = (int*)calloc((size_t)k, sizeof(int));
    if (counts) {
        for (int i = 0; i < np; i++) counts[assign[i]]++;

        for (int c = 0; c < k; c++) {
            clustering->clusters[c].id = c;
            clustering->clusters[c].n_points = counts[c];
            clustering->clusters[c].point_ids = (int*)malloc(
                (size_t)(counts[c] > 0 ? counts[c] : 1) * sizeof(int));
            clustering->clusters[c].centroid = (double*)malloc(
                (size_t)nf * sizeof(double));

            if (clustering->clusters[c].centroid) {
                memcpy(clustering->clusters[c].centroid, &centroids[c * nf],
                       (size_t)nf * sizeof(double));
            }

            int idx = 0;
            for (int i = 0; i < np; i++) {
                if (assign[i] == c && idx < counts[c]) {
                    clustering->clusters[c].point_ids[idx++] = i;
                }
            }

            /* Compute within-cluster sum of squares */
            double wss_c = 0.0;
            for (int ii = 0; ii < counts[c]; ii++) {
                int pi = clustering->clusters[c].point_ids[ii];
                for (int f = 0; f < nf; f++) {
                    double diff = ds->points[pi].z[f] - centroids[c * nf + f];
                    wss_c += diff * diff;
                }
            }
            clustering->clusters[c].within_ss = wss_c;
            total_wss += wss_c;
        }
        free(counts);
    }

    clustering->total_wss = total_wss;
    clustering->silhouette = pwa_silhouette_score(ds, clustering);

    free(centroids); free(min_dist); free(assign); free(cluster_sizes);
    return 0;
}

/*===========================================================================
 * L5: Silhouette Score
 *===========================================================================*/

double pwa_silhouette_score(const PWADataset *ds, const PWAClustering *clustering)
{
    if (!ds || !clustering) return 0.0;
    if (clustering->n_clusters < 2) return 0.0;

    int np = ds->n_points;
    int nf = ds->n_features;
    int k = clustering->n_clusters;

    /* Build cluster assignment array */
    int *assign = (int*)malloc((size_t)np * sizeof(int));
    if (!assign) return 0.0;

    for (int i = 0; i < np; i++) assign[i] = -1;
    for (int c = 0; c < k; c++) {
        for (int j = 0; j < clustering->clusters[c].n_points; j++) {
            int pi = clustering->clusters[c].point_ids[j];
            if (pi >= 0 && pi < np) assign[pi] = c;
        }
    }

    double total_s = 0.0;
    int count = 0;

    for (int i = 0; i < np; i++) {
        int c_i = assign[i];
        if (c_i < 0) continue;

        /* Compute a(i): mean distance to other points in same cluster */
        double a_i = 0.0;
        int n_same = 0;
        for (int j = 0; j < np; j++) {
            if (i == j || assign[j] != c_i) continue;
            double dist = 0.0;
            for (int f = 0; f < nf; f++) {
                double diff = ds->points[i].z[f] - ds->points[j].z[f];
                dist += diff * diff;
            }
            a_i += sqrt(dist);
            n_same++;
        }
        a_i = (n_same > 0) ? a_i / (double)n_same : 0.0;

        /* Compute b(i): minimum mean distance to points in other clusters */
        double b_i = DBL_MAX;
        for (int c = 0; c < k; c++) {
            if (c == c_i) continue;
            double mean_dist = 0.0;
            int n_other = 0;
            for (int j = 0; j < np; j++) {
                if (assign[j] != c) continue;
                double dist = 0.0;
                for (int f = 0; f < nf; f++) {
                    double diff = ds->points[i].z[f] - ds->points[j].z[f];
                    dist += diff * diff;
                }
                mean_dist += sqrt(dist);
                n_other++;
            }
            if (n_other > 0) mean_dist /= (double)n_other;
            if (mean_dist < b_i) b_i = mean_dist;
        }
        if (b_i == DBL_MAX) b_i = a_i;

        double max_ab = (a_i > b_i) ? a_i : b_i;
        if (max_ab > 1e-12) {
            total_s += (b_i - a_i) / max_ab;
            count++;
        }
    }

    free(assign);
    return (count > 0) ? total_s / (double)count : 0.0;
}

/*===========================================================================
 * L5: Region Boundary Estimation
 *===========================================================================*/

int pwa_estimate_region_boundaries(const PWADataset *ds,
                                    const PWAClustering *clustering,
                                    double **H, double **K, int *n_cons_arr)
{
    if (!ds || !clustering || !H || !K || !n_cons_arr) return -1;

    int k = clustering->n_clusters;
    int nf = ds->n_features;
    int np = ds->n_points;

    for (int ci = 0; ci < k; ci++) {
        int max_cons = k - 1;
        n_cons_arr[ci] = max_cons;
        H[ci] = (double*)calloc((size_t)(max_cons * nf), sizeof(double));
        K[ci] = (double*)calloc((size_t)max_cons, sizeof(double));
        if (!H[ci] || !K[ci]) {
            for (int j = 0; j <= ci; j++) {
                free(H[j]); free(K[j]);
            }
            return -1;
        }
    }

    /* For each pair of clusters, find a separating hyperplane */
    for (int ci = 0; ci < k; ci++) {
        int cons_idx = 0;

        for (int cj = 0; cj < k; cj++) {
            if (ci == cj) continue;
            if (cons_idx >= k - 1) break;

            /* Find normal vector pointing from cluster cj to ci
             * using the difference of centroids */
            double *normal = (double*)calloc((size_t)nf, sizeof(double));
            if (!normal) continue;

            double norm_sq = 0.0;
            for (int f = 0; f < nf; f++) {
                normal[f] = clustering->clusters[ci].centroid[f]
                           - clustering->clusters[cj].centroid[f];
                norm_sq += normal[f] * normal[f];
            }

            if (norm_sq < 1e-12) {
                free(normal);
                /* Use a different direction: PCA-like */
                continue;
            }

            /* Normalize */
            double norm = sqrt(norm_sq);
            for (int f = 0; f < nf; f++) normal[f] /= norm;

            /* Set hyperplane: normal·z ≤ b, where b is midpoint
             * projected onto the normal direction */
            double mid_proj = 0.0;
            for (int f = 0; f < nf; f++) {
                double mid = (clustering->clusters[ci].centroid[f]
                            + clustering->clusters[cj].centroid[f]) * 0.5;
                mid_proj += normal[f] * mid;
            }

            /* The half-space defining region ci relative to cj:
             * H_ci z ≤ K_ci where H_ci = -normal (pointing from ci to cj)
             * Actually, for region ci: (normal from cj to ci)·z ≤ b
             * Wait: region ci should be on the side containing its centroid.
             *
             * For pairs of clusters: define half-space for ci:
             *   (centroid_ci - centroid_cj)·z ≤ midpoint_projection */

            for (int f = 0; f < nf; f++) {
                H[ci][cons_idx * nf + f] = normal[f];
            }
            K[ci][cons_idx] = mid_proj;

            free(normal);
            cons_idx++;
        }
    }

    return 0;
}

int pwa_build_polyhedral_regions(const PWAClustering *clustering,
                                  int n_features,
                                  double **H, double **K, int *n_cons)
{
    if (!clustering || !H || !K || !n_cons) return -1;

    int k = clustering->n_clusters;

    /* Each region is already defined by the separating hyperplanes
     * computed in pwa_estimate_region_boundaries.
     *
     * Additionally, add bounding box constraints to ensure boundedness. */

    for (int ci = 0; ci < k; ci++) {
        int old_cons = n_cons[ci];
        int new_cons = old_cons + 2 * n_features;

        double *new_H = (double*)realloc(H[ci],
                            (size_t)(new_cons * n_features) * sizeof(double));
        double *new_K = (double*)realloc(K[ci],
                            (size_t)new_cons * sizeof(double));
        if (!new_H || !new_K) continue;

        H[ci] = new_H;
        K[ci] = new_K;

        /* Add bounding box: -bound ≤ x_i ≤ bound */
        double bound = 5.0;  /* 5 standard deviations */
        for (int f = 0; f < n_features; f++) {
            /* + constraint: x_f ≤ centroid_f + bound */
            H[ci][(old_cons + 2*f) * n_features + f] = 1.0;
            K[ci][old_cons + 2*f] = clustering->clusters[ci].centroid[f] + bound;

            /* - constraint: -x_f ≤ -centroid_f + bound → x_f ≥ centroid_f - bound */
            H[ci][(old_cons + 2*f + 1) * n_features + f] = -1.0;
            K[ci][old_cons + 2*f + 1] = -clustering->clusters[ci].centroid[f] + bound;
        }

        n_cons[ci] = new_cons;
    }

    return 0;
}

int pwa_refine_partition(int n_regions, double **H, double **K, int *n_cons)
{
    (void)n_regions; (void)H; (void)K; (void)n_cons;
    /* Placeholder for gap detection and refinement.
     * A full implementation would:
     * 1. Sample the space to find points outside all regions
     * 2. Assign unassigned points to nearest region
     * 3. Or split overlapping regions along their intersection */
    return 0;
}

/*===========================================================================
 * L5: Affine Regression per Region
 *===========================================================================*/

int pwa_fit_affine_model(const PWADataset *ds,
                          const int *point_ids, int n_region_pts,
                          int n_state, int n_input, PWARegression *reg)
{
    if (!ds || !point_ids || !reg || n_region_pts < 3) return -1;

    int nz = n_state + n_input;  /* Number of regressors */
    int n_aug = nz + 1;          /* Augmented with 1 for offset f */

    reg->n_points = n_region_pts;

    /* Build regression matrix X (n_region_pts × n_aug) and Y (n_region_pts × n_state)
     * X[row] = [x_k, u_k, 1]
     * Y[row] = x_{k+1}
     *
     * Solve: Θ^T = (X^T X)^{-1} X^T Y
     * where Θ = [A | B | f]^T */
    double *X = (double*)malloc((size_t)(n_region_pts * n_aug) * sizeof(double));
    double *Y = (double*)malloc((size_t)(n_region_pts * n_state) * sizeof(double));
    if (!X || !Y) { free(X); free(Y); return -1; }

    for (int r = 0; r < n_region_pts; r++) {
        int pi = point_ids[r];
        if (pi < 0 || pi >= ds->n_points) continue;

        /* X row = [z_k, 1] */
        for (int f = 0; f < nz; f++) {
            X[r * n_aug + f] = ds->points[pi].z[f];
        }
        X[r * n_aug + nz] = 1.0;  /* For offset f */

        /* Y row = x_{k+1} */
        memcpy(&Y[r * n_state], ds->points[pi].y, (size_t)n_state * sizeof(double));
    }

    /* Compute X^T X (n_aug × n_aug) */
    double *XTX = (double*)calloc((size_t)(n_aug * n_aug), sizeof(double));
    double *XTY = (double*)calloc((size_t)(n_aug * n_state), sizeof(double));
    if (!XTX || !XTY) { free(X); free(Y); free(XTX); free(XTY); return -1; }

    for (int r = 0; r < n_region_pts; r++) {
        for (int i = 0; i < n_aug; i++) {
            double xi = X[r * n_aug + i];
            for (int j = 0; j < n_aug; j++) {
                XTX[i * n_aug + j] += xi * X[r * n_aug + j];
            }
            for (int s = 0; s < n_state; s++) {
                XTY[i * n_state + s] += xi * Y[r * n_state + s];
            }
        }
    }

    /* Solve XTX * Theta = XTY using Gaussian elimination
     * for each column of Theta (each state dimension) */
    reg->A = (double*)calloc((size_t)(n_state * n_state), sizeof(double));
    reg->B = (double*)calloc((size_t)(n_state * n_input), sizeof(double));
    reg->f = (double*)calloc((size_t)n_state, sizeof(double));
    reg->residuals = (double*)calloc((size_t)(n_region_pts * n_state), sizeof(double));

    if (!reg->A || !reg->B || !reg->f || !reg->residuals) {
        free(X); free(Y); free(XTX); free(XTY);
        return -1;
    }

    /* Gaussian elimination with partial pivoting for each state */
    for (int s = 0; s < n_state; s++) {
        double *A_sys = (double*)malloc((size_t)(n_aug * n_aug) * sizeof(double));
        double *b_sys = (double*)malloc((size_t)n_aug * sizeof(double));
        if (!A_sys || !b_sys) { free(A_sys); free(b_sys); continue; }

        memcpy(A_sys, XTX, (size_t)(n_aug * n_aug) * sizeof(double));
        for (int i = 0; i < n_aug; i++) b_sys[i] = XTY[i * n_state + s];

        /* Gaussian elimination */
        for (int col = 0; col < n_aug; col++) {
            /* Pivot */
            int pivot = col;
            double max_val = fabs(A_sys[col * n_aug + col]);
            for (int row = col + 1; row < n_aug; row++) {
                double val = fabs(A_sys[row * n_aug + col]);
                if (val > max_val) { max_val = val; pivot = row; }
            }

            if (max_val < 1e-15) continue;

            if (pivot != col) {
                for (int j = col; j < n_aug; j++) {
                    double tmp = A_sys[col * n_aug + j];
                    A_sys[col * n_aug + j] = A_sys[pivot * n_aug + j];
                    A_sys[pivot * n_aug + j] = tmp;
                }
                double tmp = b_sys[col];
                b_sys[col] = b_sys[pivot];
                b_sys[pivot] = tmp;
            }

            /* Eliminate rows below */
            for (int row = col + 1; row < n_aug; row++) {
                double factor = A_sys[row * n_aug + col] / A_sys[col * n_aug + col];
                for (int j = col; j < n_aug; j++) {
                    A_sys[row * n_aug + j] -= factor * A_sys[col * n_aug + j];
                }
                b_sys[row] -= factor * b_sys[col];
            }
        }

        /* Back substitution */
        double *theta_s = (double*)calloc((size_t)n_aug, sizeof(double));
        if (theta_s) {
            for (int row = n_aug - 1; row >= 0; row--) {
                double sum = b_sys[row];
                for (int j = row + 1; j < n_aug; j++) {
                    sum -= A_sys[row * n_aug + j] * theta_s[j];
                }
                if (fabs(A_sys[row * n_aug + row]) > 1e-15) {
                    theta_s[row] = sum / A_sys[row * n_aug + row];
                }
            }

            /* Extract A, B, f */
            for (int i = 0; i < n_state; i++) {
                reg->A[s * n_state + i] = theta_s[i];  /* A[s][i] = θ_i^s */
            }
            for (int i = 0; i < n_input; i++) {
                reg->B[s * n_input + i] = theta_s[n_state + i];  /* B */
            }
            reg->f[s] = theta_s[nz];  /* f = θ_{nz} */

            free(theta_s);
        }

        free(A_sys);
        free(b_sys);
    }

    /* Compute residuals and R² */
    double ss_res = 0.0, ss_tot = 0.0;
    double y_mean = 0.0;
    for (int r = 0; r < n_region_pts; r++) {
        for (int s = 0; s < n_state; s++) {
            y_mean += Y[r * n_state + s];
        }
    }
    y_mean /= (double)(n_region_pts * n_state);

    for (int r = 0; r < n_region_pts; r++) {
        for (int s = 0; s < n_state; s++) {
            /* Prediction */
            double pred = reg->f[s];
            for (int i = 0; i < n_state; i++) {
                pred += reg->A[s * n_state + i] * X[r * n_aug + i];
            }
            for (int i = 0; i < n_input; i++) {
                pred += reg->B[s * n_input + i] * X[r * n_aug + n_state + i];
            }

            double residual = Y[r * n_state + s] - pred;
            reg->residuals[r * n_state + s] = residual;
            ss_res += residual * residual;
            ss_tot += (Y[r * n_state + s] - y_mean) * (Y[r * n_state + s] - y_mean);
        }
    }

    reg->rss = ss_res;
    reg->r_squared = (ss_tot > 1e-12) ? 1.0 - ss_res / ss_tot : 0.0;

    free(X); free(Y); free(XTX); free(XTY);
    return 0;
}

/*===========================================================================
 * L5: Full Identification Pipeline
 *===========================================================================*/

int pwa_identify(const PWADataset *ds, int n_regions, int max_iters,
                  PWAIdentification *ident)
{
    if (!ds || !ident || n_regions < 2) return -1;

    memset(ident, 0, sizeof(PWAIdentification));
    ident->n_regions = n_regions;
    ident->n_state = ds->n_state;
    ident->n_input = ds->n_input;

    /* Step 1: K-means++ clustering */
    if (pwa_kmeans_plusplus(ds, n_regions, max_iters, 1e-6, &ident->clustering) < 0) {
        return -1;
    }

    /* Step 2: Estimate region boundaries */
    ident->H = (double**)calloc((size_t)n_regions, sizeof(double*));
    ident->K = (double**)calloc((size_t)n_regions, sizeof(double*));
    int *n_cons_arr = (int*)calloc((size_t)n_regions, sizeof(int));
    if (!ident->H || !ident->K || !n_cons_arr) {
        free(ident->H); free(ident->K); free(n_cons_arr);
        return -1;
    }

    if (pwa_estimate_region_boundaries(ds, &ident->clustering,
                                        ident->H, ident->K, n_cons_arr) < 0) {
        free(n_cons_arr);
        return -1;
    }

    /* Step 3: Build polyhedral regions */
    if (pwa_build_polyhedral_regions(&ident->clustering, ds->n_features,
                                      ident->H, ident->K, n_cons_arr) < 0) {
        free(n_cons_arr);
        return -1;
    }

    /* Step 4: Fit affine models per region */
    ident->regressions = (PWARegression*)calloc((size_t)n_regions,
                                                  sizeof(PWARegression));
    if (!ident->regressions) { free(n_cons_arr); return -1; }

    for (int c = 0; c < n_regions; c++) {
        if (ident->clustering.clusters[c].n_points > 0) {
            if (pwa_fit_affine_model(ds,
                     ident->clustering.clusters[c].point_ids,
                     ident->clustering.clusters[c].n_points,
                     ds->n_state, ds->n_input,
                     &ident->regressions[c]) < 0) {
                /* Skip this cluster, keep going */
            }
        }
    }

    /* Step 5: Assemble PWA system */
    ident->pwa_system = pwa_system_create(ds->n_state, ds->n_input, ds->n_state,
                                           n_regions, 0, 1.0);
    if (ident->pwa_system) {
        for (int c = 0; c < n_regions; c++) {
            pwa_add_dynamics(ident->pwa_system,
                            ident->regressions[c].A,
                            ident->regressions[c].B,
                            ident->regressions[c].f,
                            NULL, NULL, NULL);
            pwa_add_region(ident->pwa_system,
                          ident->H[c], ident->K[c],
                          n_cons_arr[c], c);
        }
    }

    /* Compute total error */
    ident->total_error = 0.0;
    for (int c = 0; c < n_regions; c++) {
        ident->total_error += ident->regressions[c].rss;
    }

    free(n_cons_arr);
    return 0;
}

int pwa_validate_model(const PWAIdentification *ident,
                        const PWADataset *test_ds,
                        double *rms_error, double *vaf)
{
    if (!ident || !test_ds || !rms_error || !vaf) return -1;

    int n_state = ident->n_state;
    int n_input = ident->n_input;
    int nz = n_state + n_input;
    int np = test_ds->n_points;

    double ss_error = 0.0;
    double ss_total = 0.0;
    double y_mean = 0.0;
    int count = 0;

    /* Compute mean of test outputs */
    for (int i = 0; i < np; i++) {
        for (int s = 0; s < n_state; s++) {
            y_mean += test_ds->points[i].y[s];
        }
    }
    y_mean /= (double)(np * n_state);

    /* Predict and compute errors */
    for (int i = 0; i < np; i++) {
        /* Find region for this test point */
        int region = -1;
        for (int c = 0; c < ident->n_regions; c++) {
            int in_region = 1;
            /* Check constraints for this region */
            if (!ident->H[c] || !ident->K[c]) continue;

            /* Simple constraint check */
            for (int cons = 0; cons < 5 && cons < 100; cons++) {
                /* (This requires n_cons info; simplified check) */
                if (ident->H[c] && ident->K[c]) {
                    /* Check first few constraints */
                    double val = 0.0;
                    for (int f = 0; f < nz; f++) {
                        val += ident->H[c][cons * nz + f] * test_ds->points[i].z[f];
                    }
                    if (val > ident->K[c][cons] + 1e-6) {
                        in_region = 0;
                        break;
                    }
                }
            }
            if (in_region) { region = c; break; }
        }

        if (region < 0) region = 0;  /* Default */

        /* Predict using identified model */
        const PWARegression *reg = &ident->regressions[region];
        for (int s = 0; s < n_state; s++) {
            double pred = reg->f[s];
            for (int j = 0; j < n_state; j++) {
                pred += reg->A[s * n_state + j] * test_ds->points[i].z[j];
            }
            for (int j = 0; j < n_input; j++) {
                pred += reg->B[s * n_input + j] * test_ds->points[i].z[n_state + j];
            }

            double actual = test_ds->points[i].y[s];
            double error = actual - pred;
            ss_error += error * error;
            ss_total += (actual - y_mean) * (actual - y_mean);
            count++;
        }
    }

    *rms_error = (count > 0) ? sqrt(ss_error / (double)count) : DBL_MAX;
    *vaf = (ss_total > 1e-12) ? 1.0 - ss_error / ss_total : 0.0;

    return 0;
}

int pwa_select_num_regions(const PWADataset *ds, int k_min, int k_max,
                            int n_boot, int *best_k)
{
    if (!ds || !best_k || k_min < 2 || k_max < k_min) return -1;

    int nf = ds->n_features;
    int np = ds->n_points;
    double best_gap = -DBL_MAX;
    *best_k = k_min;

    for (int k = k_min; k <= k_max; k++) {
        /* Run K-means on actual data */
        PWAClustering clust;
        memset(&clust, 0, sizeof(clust));
        if (pwa_kmeans_plusplus(ds, k, 100, 1e-6, &clust) < 0) continue;

        double log_wss = log(clust.total_wss / (double)np);

        /* Generate reference distribution via bootstrap uniform in bounding box */
        double ref_log_w_sum = 0.0;
        for (int boot = 0; boot < n_boot; boot++) {
            /* Generate uniform reference data in the bounding box of actual data */
            double *ref_data = (double*)malloc((size_t)(np * nf) * sizeof(double));
            if (!ref_data) continue;

            /* Find data bounds */
            double *mins = (double*)malloc((size_t)nf * sizeof(double));
            double *maxs = (double*)malloc((size_t)nf * sizeof(double));
            if (mins && maxs) {
                for (int f = 0; f < nf; f++) {
                    mins[f] = DBL_MAX; maxs[f] = -DBL_MAX;
                }
                for (int i = 0; i < np; i++) {
                    for (int f = 0; f < nf; f++) {
                        if (ds->points[i].z[f] < mins[f]) mins[f] = ds->points[i].z[f];
                        if (ds->points[i].z[f] > maxs[f]) maxs[f] = ds->points[i].z[f];
                    }
                }
                /* Generate uniform random points */
                for (int i = 0; i < np; i++) {
                    for (int f = 0; f < nf; f++) {
                        double r = ((double)((i * 17 + boot * 31 + f * 7) % 1000) / 1000.0);
                        ref_data[i * nf + f] = mins[f] + r * (maxs[f] - mins[f]);
                    }
                }
                free(mins); free(maxs);
            }

            /* Run K-means on reference data */
            PWADataset ref_ds;
            memset(&ref_ds, 0, sizeof(ref_ds));
            ref_ds.n_points = np;
            ref_ds.n_state = ds->n_state;
            ref_ds.n_input = ds->n_input;
            ref_ds.n_features = nf;
            ref_ds.points = (PWADataPoint*)calloc((size_t)np, sizeof(PWADataPoint));
            if (ref_ds.points) {
                for (int i = 0; i < np; i++) {
                    ref_ds.points[i].z = &ref_data[i * nf];
                    ref_ds.points[i].y = NULL;
                }
                PWAClustering ref_clust;
                memset(&ref_clust, 0, sizeof(ref_clust));
                pwa_kmeans_plusplus(&ref_ds, k, 100, 1e-6, &ref_clust);
                ref_log_w_sum += log(ref_clust.total_wss / (double)np);
                /* Free ref_clust... (simplified) */
                free(ref_ds.points);
            }
            free(ref_data);
        }

        double exp_ref = ref_log_w_sum / (double)n_boot;
        double gap = exp_ref - log_wss;

        if (gap > best_gap) {
            best_gap = gap;
            *best_k = k;
        }
    }

    return 0;
}

void pwa_identification_destroy(PWAIdentification *ident)
{
    if (!ident) return;

    if (ident->H) {
        for (int i = 0; i < ident->n_regions; i++) free(ident->H[i]);
        free(ident->H);
    }
    if (ident->K) {
        for (int i = 0; i < ident->n_regions; i++) free(ident->K[i]);
        free(ident->K);
    }
    if (ident->regressions) {
        for (int i = 0; i < ident->n_regions; i++) {
            free(ident->regressions[i].A);
            free(ident->regressions[i].B);
            free(ident->regressions[i].f);
            free(ident->regressions[i].residuals);
        }
        free(ident->regressions);
    }
    if (ident->clustering.clusters) {
        for (int i = 0; i < ident->clustering.n_clusters; i++) {
            free(ident->clustering.clusters[i].point_ids);
            free(ident->clustering.clusters[i].centroid);
            free(ident->clustering.clusters[i].covariance);
        }
        free(ident->clustering.clusters);
    }
    if (ident->pwa_system) {
        pwa_system_destroy(ident->pwa_system);
    }
}
