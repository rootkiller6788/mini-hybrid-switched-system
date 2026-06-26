#ifndef IMPULSIVE_SYNC_H
#define IMPULSIVE_SYNC_H
/*
 * impulsive_sync.h -- Impulsive synchronization of chaotic systems
 *
 * Master-slave synchronization via impulsive control:
 *
 *   Master:  dx/dt = f(t, x)        (autonomous chaotic system)
 *   Slave:   dy/dt = f(t, y), t != tau_k
 *            Delta y = -B * (y - x), t = tau_k
 *
 * The slave synchronizes to the master through discrete-time
 * impulsive corrections. Error dynamics:
 *
 *   e = y - x
 *   de/dt = f(t, y) - f(t, x) = Df(xi) * e  (linearized)
 *   Delta e = (I - B) * e  at t = tau_k
 *
 * Synchronization condition (Yang 2001):
 *   If the largest Lyapunov exponent of the flow is lambda_max,
 *   and impulses occur with period Delta, then synchronization
 *   requires:  rho = ||I - B|| < exp(-lambda_max * Delta)
 *
 * Applications:
 *   - Secure communication (chaotic masking)
 *   - Spread-spectrum communication
 *   - Biological rhythm synchronization (circadian clocks)
 *   - Multi-agent coordination with intermittent communication
 *
 * References:
 *   Yang & Chua (1997) "Impulsive stabilization for control and
 *     synchronization of chaotic systems"
 *   Yang (2001) "Impulsive Control Theory", Chapter 5
 *   Stojanovski, Kocarev, Parlitz (1996) "Driving and synchronizing
 *     by chaotic impulses"
 */

#include "impulsive_types.h"
#include "impulsive_solver.h"

/* -- Synchronization Configuration -- */

typedef enum {
    IMP_SYNC_LINEAR_ERROR    = 0,  /* e = y - x (linear error) */
    IMP_SYNC_COMPLETE        = 1,  /* full state synchronization */
    IMP_SYNC_PROJECTIVE      = 2,  /* y -> alpha * x (scaled sync) */
    IMP_SYNC_PHASE           = 3,  /* phase synchronization only */
    IMP_SYNC_GENERALIZED     = 4,  /* generalized synchro via function H */
    IMP_SYNC_LAG             = 5   /* lag synchronization y(t) = x(t-tau) */
} ImpSyncType;

typedef struct {
    ImpSyncType  sync_type;
    double       sync_threshold;    /* ||e|| < threshold => synchronized */
    double       coupling_gain;     /* global coupling strength */
    int          sync_vars_mask;    /* bitmask: which states to sync */
    bool         use_adaptive_gain; /* adapt coupling online */
    double       gain_adapt_rate;
    double       min_gain;
    double       max_gain;
} ImpSyncConfig;

ImpSyncConfig imp_sync_config_default(void);

/* -- Synchronization Error Metrics -- */

typedef struct {
    double  error_norm;        /* ||e|| at current time */
    double  max_error;         /* max ||e|| over simulation */
    double  mean_error;        /* time-averaged ||e|| */
    double  sync_time;         /* first time ||e|| < threshold */
    bool    is_synchronized;   /* true if currently synched */
} ImpSyncMetrics;

/* -- Core Synchronization Functions -- */

int imp_sync_error_compute(const double *x_master,
                            const double *x_slave,
                            int n, ImpSyncType type,
                            double *error);

int imp_sync_master_slave_simulate(const ImpSystem *master,
                                    const ImpSystem *slave,
                                    const ImpJumpMap *sync_jumps,
                                    void **sync_ctxs,
                                    int num_sync_jumps,
                                    const ImpSyncConfig *cfg,
                                    const ImpSolverConfig *sol_cfg,
                                    ImpSolution *sol_master,
                                    ImpSolution *sol_slave,
                                    ImpSyncMetrics *metrics);

/* -- Impulsive Synchronization Jump Maps -- */

typedef struct {
    double *B;       /* [n*n] linear error feedback gain */
    int     n;
} ImpSyncJumpLinear;

ImpSyncJumpLinear* imp_sync_jump_linear_create(const double *B, int n);
void imp_sync_jump_linear_free(ImpSyncJumpLinear *sjl);
int imp_sync_jump_linear_eval(double t_k, const double *e_before, int n,
                               double *e_after, void *ctx);

typedef struct {
    double  gain;      /* scalar gain */
    double *mask;      /* [n] which components to correct */
    int     n;
} ImpSyncJumpScalar;

ImpSyncJumpScalar* imp_sync_jump_scalar_create(double gain, int n);
void imp_sync_jump_scalar_free(ImpSyncJumpScalar *sjs);
int imp_sync_jump_scalar_eval(double t_k, const double *e_before, int n,
                               double *e_after, void *ctx);

/* -- Synchronization Stability Analysis -- */

double imp_sync_estimate_largest_lyapunov(ImpVectorField f, void *ctx,
                                           const double *x0, int n,
                                           double T, double dt);

bool imp_sync_check_synchronization_condition(double lambda_max,
                                                double impulse_period,
                                                double max_eigenvalue_B);

int imp_sync_compute_error_derivative(ImpVectorField f_master, void *ctx_m,
                                       ImpVectorField f_slave, void *ctx_s,
                                       double t, const double *x_m,
                                       const double *x_s, int n,
                                       double *dedt);

/* -- Applications: Chaotic Masking -- */

typedef struct {
    double *message;     /* [msg_len] original signal */
    double *encrypted;   /* [msg_len] encrypted/chaos-masked signal */
    double *decrypted;   /* [msg_len] recovered signal */
    int     msg_len;
    double  sample_dt;
} ImpSyncChaoticMasking;

ImpSyncChaoticMasking* imp_sync_chaotic_masking_create(int msg_len, double dt);
void imp_sync_chaotic_masking_free(ImpSyncChaoticMasking *cm);
int imp_sync_chaotic_masking_encode(ImpSyncChaoticMasking *cm,
                                     const double *chaos_signal);
int imp_sync_chaotic_masking_decode(ImpSyncChaoticMasking *cm,
                                     const double *chaos_signal);
double imp_sync_chaotic_masking_snr(const ImpSyncChaoticMasking *cm);

#endif /* IMPULSIVE_SYNC_H */
