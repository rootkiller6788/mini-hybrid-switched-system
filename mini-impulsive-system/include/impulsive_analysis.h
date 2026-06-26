#ifndef IMPULSIVE_ANALYSIS_H
#define IMPULSIVE_ANALYSIS_H
/*
 * impulsive_analysis.h -- System analysis for impulsive dynamical systems
 *
 * Analysis tools: dwell-time, invariant sets, basin of attraction,
 * robustness, frequency-domain, monodromy/Floquet theory.
 *
 * References:
 *   Goebel, Sanfelice, Teel (2012) "Hybrid Dynamical Systems"
 *   Haddad et al. (2006) "Impulsive and Hybrid Dynamical Systems"
 *   Bainov & Simeonov (1989) "Systems with Impulse Effect"
 */
#include "impulsive_types.h"
#include "impulsive_lyapunov.h"
#include "impulsive_solver.h"

int imp_analysis_dwell_stats(const ImpTimeSeq *seq, ImpDwellStats *stats);
bool imp_analysis_dwell_satisfies_bound(const ImpTimeSeq *seq, double min_dwell);
double imp_analysis_compute_required_dwell(double c_flow, double d_jump, double safety);
int imp_analysis_maximum_dwell_interval(const ImpTimeSeq *seq, double *start, double *end, double *max_len);

typedef struct {
    double *center; double *shape; double radius_sq; int n;
} ImpInvariantSet;

ImpInvariantSet* imp_analysis_invariant_set_create(int n);
void imp_analysis_invariant_set_free(ImpInvariantSet *inv);
int imp_analysis_invariant_set_compute(const ImpSystem *sys, ImpLyapunovFn *lyap,
                                        double alpha, double rho, ImpInvariantSet *inv);
bool imp_analysis_is_in_set(const ImpInvariantSet *inv, const double *x);

typedef struct {
    double *grid; int *labels; double *lb; double *ub; int N, n; double volume;
} ImpBasinAttraction;

ImpBasinAttraction* imp_analysis_basin_create(const double *lb, const double *ub, int N, int n);
void imp_analysis_basin_free(ImpBasinAttraction *basin);
int imp_analysis_basin_estimate(ImpBasinAttraction *basin, const ImpSystem *sys,
                                 const ImpSolverConfig *cfg, double max_T);
double imp_analysis_basin_volume(ImpBasinAttraction *basin);

typedef struct {
    double *delta_t; double *delta_x;
    double sensitivity_t, sensitivity_x; int nimp, n;
} ImpRobustness;

ImpRobustness* imp_analysis_robustness_create(int nimp, int n);
void imp_analysis_robustness_free(ImpRobustness *rob);
int imp_analysis_robustness_compute(const ImpSystem *sys, const double *x0,
                                     const ImpSolverConfig *cfg, ImpRobustness *rob);

typedef struct {
    double *freq; double *mag; double *phase; int Nf; double f_min, f_max;
} ImpFreqResponse;

ImpFreqResponse* imp_analysis_freq_response_create(int Nf, double f_min, double f_max);
void imp_analysis_freq_response_free(ImpFreqResponse *fr);
int imp_analysis_freq_response_compute(const ImpSystem *sys, const double *x0,
                                        int n, ImpFreqResponse *fr);

int imp_analysis_monodromy_matrix(const ImpSystem *sys, const ImpSolverConfig *cfg,
                                   double *Phi, int n);
double imp_analysis_floquet_multipliers(const double *Phi, int n, double *multipliers);
double imp_analysis_impulse_jitter_margin(const ImpSystem *sys, const double *x0,
                                            int n, const ImpSolverConfig *cfg);

int imp_analysis_comparison_check(const ImpSystem *sys, ImpLyapunovFn *lyap,
                                   const ImpComparisonSys *comp, double *v_bound);

double imp_analysis_energy(const double *x, int n);
double imp_analysis_total_impulse_energy(const ImpSolution *sol);
double imp_analysis_L2_gain(const ImpSolution *sol_input, const ImpSolution *sol_output);

#endif
