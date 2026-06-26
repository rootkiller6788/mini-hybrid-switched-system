#ifndef DTA_AVERAGE_H
#define DTA_AVERAGE_H
#include "dta_core.h"
#include "dta_switch_signal.h"

/* ==============================================================
 * dta_average.h - Average Dwell-Time Theory
 *
 * Average dwell time (ADT) relaxes the constant dwell time
 * requirement. A switching signal σ has average dwell time τ_a
 * if there exists N0 > 0 such that:
 *
 *   N_σ(T, t) ≤ N0 + (T - t) / τ_a   for all T ≥ t ≥ 0
 *
 * where N_σ(T,t) = number of switches in interval [t, T].
 *
 * N0 is the "chatter bound" — allows for occasional fast switching.
 *
 * Theorem (Hespanha & Morse 1999):
 *   If all subsystems are stable with
 *     ||e^{A_i t}|| ≤ c e^{-λ_i t} for some c≥1, λ_i>0,
 *   then the switched system is GAS for any switching signal with
 *     τ_a > ln(c) / λ,  where λ = min_i λ_i.
 *
 * Mode-dependent ADT (Zhao & Hill 2008):
 *   τ_{a,i} > ln(μ) / λ_i  for each mode i.
 *
 * References:
 *   Hespanha & Morse (1999) "Stability of switched systems with
 *     average dwell-time", IEEE CDC
 *   Hespanha (2004) IEEE TAC 49(4):470-482
 *   Zhao & Hill (2008) Automatica 44(7):1809-1816
 * ============================================================== */

/* --- ADT configuration --- */
typedef struct {
    double tau_a;              /* Average dwell time */
    double N0;                 /* Chatter bound */
    bool mode_dependent;       /* Use mode-dependent ADT */
    double* tau_a_i;           /* Per-mode average dwell times */
    int n_modes;
} DTA_ADTConfig;

/* --- ADT stability result --- */
typedef struct {
    DTA_StabilityVerdict verdict;
    double tau_a_min;          /* Computed minimum ADT */
    double actual_tau_a;       /* Actual ADT of the signal */
    double N0;
    double decay_rate;         /* Guaranteed α */
    bool mode_dependent_ok;
} DTA_ADTResult;

/* --- API --- */

/** Compute the minimum average dwell time from system data.
 *  τ_a* = max_i { ln(c_i) / λ_i } where c_i = cond(P_i)^{1/2} */
double dta_adt_compute_min(const DTA_SwitchedSystem* sys);

/** Compute mode-dependent minimum ADT:
 *  τ_{a,i}* = ln(μ_i) / λ_i  for each mode i */
void dta_adt_compute_mode_dependent(const DTA_SwitchedSystem* sys,
                                     double* tau_a_i_out);

/** Check if a given switching signal satisfies the ADT condition
 *  N_σ(T,t) ≤ N0 + (T-t)/τ_a */
bool dta_adt_signal_satisfies(const DTA_SwitchingSignal* sig,
                               double tau_a, double N0);

/** Stability analysis under average dwell time */
DTA_ADTResult dta_adt_analyze(const DTA_SwitchedSystem* sys,
                               const DTA_SwitchingSignal* sig,
                               double N0);

/** Compute the chatter bound N0 from a switching signal.
 *  N0 = max_{T≥t≥0} { N_σ(T,t) - (T-t)/τ_a } */
double dta_adt_compute_N0(const DTA_SwitchingSignal* sig, double tau_a);

/** Construct a switching signal with guaranteed average dwell time τ_a.
 *  Inserts fast-switch bursts compatible with the chatter bound N0. */
DTA_SwitchingSignal* dta_adt_generate_signal(double t_start, double t_end,
    int n_modes, double tau_a, double N0);

/** Rate of convergence bound for ADT switching:
 *  ||x(t)|| ≤ c e^{-α t} ||x(0)|| where α > 0 if τ_a > τ_a* */
double dta_adt_convergence_rate(const DTA_SwitchedSystem* sys,
                                 double tau_a, double N0);

/** Compare constant dwell time vs ADT: given a signal, compute both */
void dta_adt_compare(const DTA_SwitchingSignal* sig,
                      double* out_const_dwell, double* out_avg_dwell);

/** Print ADT analysis report */
void dta_adt_report(const DTA_ADTResult* result);

#endif /* DTA_AVERAGE_H */
