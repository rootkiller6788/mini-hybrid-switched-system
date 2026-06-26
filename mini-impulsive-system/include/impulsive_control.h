#ifndef IMPULSIVE_CONTROL_H
#define IMPULSIVE_CONTROL_H
/*
 * impulsive_control.h -- Impulsive control design methods
 *
 * Impulsive control applies control action only at discrete instants:
 *   dx/dt = f(t, x),            t != tau_k  (uncontrolled flow)
 *   Delta x = B_k * u_k(x),     t = tau_k   (impulsive control)
 *
 * Advantages: reduced control effort, suitable for intermittent
 * actuation, natural for sampled-data and networked control.
 *
 * Design methods: Impulsive LQR, pole placement, Lyapunov-based,
 * model predictive, adaptive.
 *
 * References:
 *   Yang (2001) "Impulsive Control Theory", Chapters 3-4
 *   Sun, Naghshtabrizi, Packard (2013)
 *   Haddad et al. (2006), Chapters 6-8
 */
#include "impulsive_types.h"

typedef struct {
    double *K; double *P; double cost;
    int n, m; bool solved;
} ImpLQR;

ImpLQR* imp_lqr_create(const double *A, const double *B,
                        const double *Q, const double *R,
                        int n, int m, double dt);
void imp_lqr_free(ImpLQR *lqr);
int imp_lqr_compute_gain(ImpLQR *lqr);
int imp_lqr_apply(const ImpLQR *lqr, const double *x, double *u, int n, int m);
double imp_lqr_cost_to_go(const ImpLQR *lqr, const double *x, int n);

typedef struct {
    double *K; double *poles; int n, m;
} ImpPolePlace;

ImpPolePlace* imp_pole_place_create(int n, int m, const double *desired_poles);
void imp_pole_place_free(ImpPolePlace *pp);
int imp_pole_place_design(ImpPolePlace *pp, const double *A, const double *B);
int imp_pole_place_apply(const ImpPolePlace *pp, const double *x, double *u);

typedef struct {
    ImpLyapunovFn *lyap; double rho_des; double *K; int n, m;
} ImpLyapControl;

ImpLyapControl* imp_lyap_control_create(ImpLyapunovFn *lyap, double rho, int n, int m);
void imp_lyap_control_free(ImpLyapControl *lc);
int imp_lyap_control_design(ImpLyapControl *lc, const double *A, const double *B);
int imp_lyap_control_apply(ImpLyapControl *lc, const double *x, double *u);

typedef struct {
    double period; double *K; double *x_ref; int n, m;
} ImpPeriodicControl;

ImpPeriodicControl* imp_periodic_control_create(double period, const double *K, int n, int m);
void imp_periodic_control_free(ImpPeriodicControl *pc);
int imp_periodic_control_update(ImpPeriodicControl *pc, double t, const double *x, double *u);

typedef struct {
    double threshold; double *K; double *x_last; double t_last;
    bool triggered; int n, m;
} ImpEventControl;

ImpEventControl* imp_event_control_create(double threshold, const double *K, int n, int m);
void imp_event_control_free(ImpEventControl *ec);
bool imp_event_control_check(ImpEventControl *ec, double t, const double *x);
int imp_event_control_trigger(ImpEventControl *ec, double t, const double *x, double *u);

typedef struct {
    double *K; double *K_prev; double gamma; double *e; int n, m;
} ImpAdaptiveControl;

ImpAdaptiveControl* imp_adaptive_control_create(double gamma, const double *K0, int n, int m);
void imp_adaptive_control_free(ImpAdaptiveControl *ac);
int imp_adaptive_control_adapt(ImpAdaptiveControl *ac, const double *x, const double *x_ref);
int imp_adaptive_control_apply(ImpAdaptiveControl *ac, const double *x, double *u);

double imp_control_effort_L2(const double *u, int m);
double imp_control_effort_peak(const double *u_seq, int num_steps, int m);
double imp_control_effort_total_variation(const double *u_seq, int num_steps, int m);

#endif
