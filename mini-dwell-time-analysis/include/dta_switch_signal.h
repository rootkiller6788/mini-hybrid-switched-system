#ifndef DTA_SWITCH_SIGNAL_H
#define DTA_SWITCH_SIGNAL_H
#include "dta_core.h"

/* ==============================================================
 * dta_switch_signal.h - Switching Signal Generation & Analysis
 *
 * A switching signal σ: [0, ∞) → Q = {1, 2, ..., m} is a
 * piecewise-constant, right-continuous function.
 *
 * Three classes of switching signals:
 * 1. Arbitrary switching — no constraints (worst case)
 * 2. Dwell-time constrained — minimum τ_d between switches
 * 3. Average dwell-time — bound on switching frequency
 *
 * References:
 *   Liberzon (2003) Ch. 3 — Stability under constrained switching
 *   Hespanha (2004) "Uniform stability of switched linear systems:
 *     extensions of LaSalle's invariance principle", IEEE TAC
 * ============================================================== */

/* --- Switching signal generation strategy --- */
typedef enum {
    DTA_SIG_PERIODIC = 0,
    DTA_SIG_RANDOM_UNIFORM = 1,
    DTA_SIG_DWELL_CONSTRAINED = 2,
    DTA_SIG_AVG_DWELL = 3,
    DTA_SIG_STATE_DEPENDENT = 4,
    DTA_SIG_TIME_VARYING = 5,
    DTA_SIG_HYSTERESIS = 6
} DTA_SignalType;

/* --- Switching signal statistics --- */
typedef struct {
    double min_dwell;
    double max_dwell;
    double avg_dwell;
    int total_switches;
    double duration;
    double switch_frequency;
    int* mode_visits;
    double* mode_durations;
    int n_modes;
} DTA_SignalStatistics;

/* --- API --- */

/** Create an empty switching signal with given capacity */
DTA_SwitchingSignal* dta_signal_create(int capacity);
void dta_signal_free(DTA_SwitchingSignal* sig);

/** Add a (time, mode) pair to the signal */
int dta_signal_append(DTA_SwitchingSignal* sig, double t, int mode);

/** Generate a periodic switching signal with period T_per */
DTA_SwitchingSignal* dta_signal_periodic(double t_start, double t_end,
    const int* mode_pattern, int pattern_len, double T_per);

/** Generate a random switching signal respecting minimum dwell time τ_d.
 *  Uses rand() with 'seed'/srand; call srand() before this. */
DTA_SwitchingSignal* dta_signal_random_dwell(double t_start, double t_end,
    int n_modes, double tau_d);

/** Generate a switching signal with constant dwell time τ_d */
DTA_SwitchingSignal* dta_signal_constant_dwell(double t_start, double t_end,
    int n_modes, double tau_d);

/** Generate a switching signal with average dwell time τ_a */
DTA_SwitchingSignal* dta_signal_average_dwell(double t_start, double t_end,
    int n_modes, double tau_a, double N0);

/** Compute statistics of a switching signal */
DTA_SignalStatistics dta_signal_statistics(const DTA_SwitchingSignal* sig,
                                            int n_modes);

/** Compute the minimum dwell time of an existing signal */
double dta_signal_min_dwell(const DTA_SwitchingSignal* sig);

/** Compute the average dwell time over [0, s]:
 *  τ_a(s) = s / N_σ(0, s) */
double dta_signal_avg_dwell_at(const DTA_SwitchingSignal* sig, double s);

/** Count switches in interval [t1, t2] */
int dta_signal_count_switches(const DTA_SwitchingSignal* sig,
                               double t1, double t2);

/** Get the active mode at time t */
int dta_signal_active_mode(const DTA_SwitchingSignal* sig, double t);

/** Validate signal respects minimum dwell time τ_d_min */
bool dta_signal_validate_dwell(const DTA_SwitchingSignal* sig,
                                double tau_d_min);

/** Validate signal respects average dwell time (N0, tau_a) */
bool dta_signal_validate_avg_dwell(const DTA_SwitchingSignal* sig,
                                    double tau_a, double N0);

/** Check if signal satisfies: N_σ(T,t) ≤ N0 + (T-t)/τ_a */
bool dta_signal_check_adt_bound(const DTA_SwitchingSignal* sig,
                                 double tau_a, double N0, double T, double t);

/** Merge two switching signals sequentially */
DTA_SwitchingSignal* dta_signal_concat(const DTA_SwitchingSignal* sig1,
                                        const DTA_SwitchingSignal* sig2);

/** Extract a sub-signal on [t1, t2] */
DTA_SwitchingSignal* dta_signal_slice(const DTA_SwitchingSignal* sig,
                                       double t1, double t2);

/** Print signal to stdout for debugging */
void dta_signal_print(const DTA_SwitchingSignal* sig);

#endif /* DTA_SWITCH_SIGNAL_H */
