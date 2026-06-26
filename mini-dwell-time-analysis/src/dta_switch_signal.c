#include "dta_core.h"
#include "dta_switch_signal.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * dta_switch_signal.c - Switching Signal Generation and Analysis
 *
 * Implements generation of periodic, random (dwell-constrained),
 * constant-dwell, and average-dwell switching signals.
 * Also provides signal statistics (min/avg/max dwell, frequency),
 * validation against dwell time constraints, and the ADT bound
 * check N_sigma(T,t) <= N0 + (T-t)/tau_a.
 *
 * Each function implements an independent knowledge point.
 * References:
 *   Liberzon (2003) "Switching in Systems and Control" Ch.3
 *   Hespanha & Morse (1999) IEEE CDC
 * ============================================================== */

static void* safe_alloc(size_t sz) {
    void* p = malloc(sz);
    if (!p) { fprintf(stderr, "DTA signal: alloc fail\n"); exit(1); }
    return p;
}

/* --- Signal lifecycle --- */

DTA_SwitchingSignal* dta_signal_create(int capacity) {
    DTA_SwitchingSignal* sig = safe_alloc(sizeof(DTA_SwitchingSignal));
    sig->capacity = capacity > 0 ? capacity : 256;
    sig->n_switches = 0;
    sig->t_start = 0.0;
    sig->t_end = 0.0;
    sig->switch_times = safe_alloc((size_t)sig->capacity * sizeof(double));
    sig->mode_sequence = safe_alloc((size_t)sig->capacity * sizeof(int));
    return sig;
}

void dta_signal_free(DTA_SwitchingSignal* sig) {
    if (!sig) return;
    free(sig->switch_times);
    free(sig->mode_sequence);
    free(sig);
}

int dta_signal_append(DTA_SwitchingSignal* sig, double t, int mode) {
    if (!sig) return -1;
    if (sig->n_switches >= sig->capacity) {
        sig->capacity *= 2;
        sig->switch_times = realloc(sig->switch_times,
            (size_t)sig->capacity * sizeof(double));
        sig->mode_sequence = realloc(sig->mode_sequence,
            (size_t)sig->capacity * sizeof(int));
        if (!sig->switch_times || !sig->mode_sequence) return -1;
    }
    sig->switch_times[sig->n_switches] = t;
    sig->mode_sequence[sig->n_switches] = mode;
    sig->n_switches++;
    sig->t_end = t;
    if (sig->n_switches == 1) sig->t_start = t;
    return 0;
}

/* --- Periodic switching signal --- */

DTA_SwitchingSignal* dta_signal_periodic(double t_start, double t_end,
    const int* mode_pattern, int pattern_len, double T_per) {
    if (!mode_pattern || pattern_len <= 0 || T_per <= 0 || t_end <= t_start)
        return NULL;
    DTA_SwitchingSignal* sig = dta_signal_create(
        (int)((t_end - t_start) / T_per) * pattern_len + 2);
    double t = t_start;
    int idx = 0;
    while (t < t_end) {
        dta_signal_append(sig, t, mode_pattern[idx % pattern_len]);
        t += T_per;
        idx++;
    }
    dta_signal_append(sig, t_end, mode_pattern[(idx - 1) % pattern_len]);
    return sig;
}

/* --- Random switching with dwell time constraint --- */

DTA_SwitchingSignal* dta_signal_random_dwell(double t_start, double t_end,
    int n_modes, double tau_d) {
    if (n_modes <= 0 || tau_d <= 0 || t_end <= t_start) return NULL;
    DTA_SwitchingSignal* sig = dta_signal_create(
        (int)((t_end - t_start) / tau_d) + 2);
    double t = t_start;
    int current_mode = rand() % n_modes;
    dta_signal_append(sig, t, current_mode);
    while (t + tau_d < t_end) {
        t += tau_d + ((double)rand() / RAND_MAX) * tau_d * 0.5;
        if (t >= t_end) break;
        int next_mode;
        do { next_mode = rand() % n_modes; } while (next_mode == current_mode);
        current_mode = next_mode;
        dta_signal_append(sig, t, current_mode);
    }
    return sig;
}

/* --- Constant dwell time signal: cycles through modes --- */

DTA_SwitchingSignal* dta_signal_constant_dwell(double t_start, double t_end,
    int n_modes, double tau_d) {
    if (n_modes <= 0 || tau_d <= 0 || t_end <= t_start) return NULL;
    DTA_SwitchingSignal* sig = dta_signal_create(
        (int)((t_end - t_start) / tau_d) + 2);
    double t = t_start;
    int mode = 0;
    dta_signal_append(sig, t, mode);
    while (t + tau_d < t_end) {
        t += tau_d;
        mode = (mode + 1) % n_modes;
        dta_signal_append(sig, t, mode);
    }
    return sig;
}

/* --- Average dwell time signal ---
 * Satisfies: N_sigma(T,t) <= N0 + (T-t)/tau_a
 * Strategy: alternate between fast and slow segments. */

DTA_SwitchingSignal* dta_signal_average_dwell(double t_start, double t_end,
    int n_modes, double tau_a, double N0) {
    if (n_modes <= 0 || tau_a <= 0 || t_end <= t_start) return NULL;
    DTA_SwitchingSignal* sig = dta_signal_create((int)((t_end - t_start) / (tau_a * 0.5)) + 10);
    double t = t_start;
    int mode = 0;
    dta_signal_append(sig, t, mode);
    while (t < t_end) {
        double segment_len = tau_a * (1.0 + 0.5 * ((double)rand() / RAND_MAX));
        t += segment_len;
        if (t >= t_end) break;
        mode = (mode + 1) % n_modes;
        dta_signal_append(sig, t, mode);
        if (N0 > 0.5 && ((double)rand() / RAND_MAX) < 0.2) {
            int burst_count = (int)(N0 * ((double)rand() / RAND_MAX));
            int bi;
            for (bi = 0; bi < burst_count && t < t_end; bi++) {
                t += tau_a * 0.1;
                if (t >= t_end) break;
                mode = (mode + 1) % n_modes;
                dta_signal_append(sig, t, mode);
            }
        }
    }
    return sig;
}

/* --- Signal statistics --- */

DTA_SignalStatistics dta_signal_statistics(const DTA_SwitchingSignal* sig,
                                            int n_modes) {
    DTA_SignalStatistics stats;
    memset(&stats, 0, sizeof(stats));
    if (!sig || sig->n_switches < 2) return stats;
    stats.n_modes = n_modes;
    stats.total_switches = sig->n_switches - 1;
    stats.duration = sig->switch_times[sig->n_switches-1] - sig->switch_times[0];
    stats.min_dwell = INFINITY;
    stats.max_dwell = 0.0;
    double sum_dwell = 0.0;
    int i;
    for (i = 1; i < sig->n_switches; i++) {
        double dwell = sig->switch_times[i] - sig->switch_times[i-1];
        if (dwell < stats.min_dwell) stats.min_dwell = dwell;
        if (dwell > stats.max_dwell) stats.max_dwell = dwell;
        sum_dwell += dwell;
    }
    stats.avg_dwell = sum_dwell / (double)stats.total_switches;
    stats.switch_frequency = (stats.duration > 0) ?
        (double)stats.total_switches / stats.duration : 0.0;
    stats.mode_visits = safe_alloc((size_t)n_modes * sizeof(int));
    stats.mode_durations = safe_alloc((size_t)n_modes * sizeof(double));
    for (i = 0; i < n_modes; i++) {
        stats.mode_visits[i] = 0;
        stats.mode_durations[i] = 0.0;
    }
    for (i = 0; i < sig->n_switches; i++) {
        int m = sig->mode_sequence[i];
        if (m >= 0 && m < n_modes) {
            stats.mode_visits[m]++;
            if (i < sig->n_switches - 1)
                stats.mode_durations[m] += sig->switch_times[i+1] - sig->switch_times[i];
        }
    }
    return stats;
}

double dta_signal_min_dwell(const DTA_SwitchingSignal* sig) {
    if (!sig || sig->n_switches < 2) return INFINITY;
    double min_d = INFINITY;
    int i;
    for (i = 1; i < sig->n_switches; i++) {
        double d = sig->switch_times[i] - sig->switch_times[i-1];
        if (d < min_d) min_d = d;
    }
    return min_d;
}

double dta_signal_avg_dwell_at(const DTA_SwitchingSignal* sig, double s) {
    if (!sig || s <= 0) return INFINITY;
    int N = dta_signal_count_switches(sig, 0.0, s);
    return (N > 0) ? s / (double)N : INFINITY;
}

int dta_signal_count_switches(const DTA_SwitchingSignal* sig,
                               double t1, double t2) {
    if (!sig || sig->n_switches < 2 || t2 <= t1) return 0;
    int count = 0, i;
    for (i = 1; i < sig->n_switches; i++) {
        double sw = sig->switch_times[i];
        if (sw > t1 && sw <= t2) count++;
    }
    return count;
}

int dta_signal_active_mode(const DTA_SwitchingSignal* sig, double t) {
    if (!sig || sig->n_switches == 0) return -1;
    if (t < sig->switch_times[0]) return sig->mode_sequence[0];
    int i;
    for (i = sig->n_switches - 1; i >= 0; i--) {
        if (sig->switch_times[i] <= t)
            return sig->mode_sequence[i];
    }
    return sig->mode_sequence[sig->n_switches - 1];
}

bool dta_signal_validate_dwell(const DTA_SwitchingSignal* sig,
                                double tau_d_min) {
    if (!sig) return false;
    return dta_signal_min_dwell(sig) >= tau_d_min;
}

bool dta_signal_validate_avg_dwell(const DTA_SwitchingSignal* sig,
                                    double tau_a, double N0) {
    if (!sig || sig->n_switches < 2) return false;
    int i, j;
    for (i = 0; i < sig->n_switches; i++) {
        for (j = i + 1; j < sig->n_switches; j++) {
            double T = sig->switch_times[j];
            double t = sig->switch_times[i];
            if (T <= t) continue;
            int N = dta_signal_count_switches(sig, t, T);
            if ((double)N > N0 + (T - t) / tau_a)
                return false;
        }
    }
    return true;
}

bool dta_signal_check_adt_bound(const DTA_SwitchingSignal* sig,
                                 double tau_a, double N0,
                                 double T, double t) {
    if (!sig) return false;
    int N = dta_signal_count_switches(sig, t, T);
    return (double)N <= N0 + (T - t) / tau_a;
}

DTA_SwitchingSignal* dta_signal_concat(const DTA_SwitchingSignal* sig1,
                                        const DTA_SwitchingSignal* sig2) {
    if (!sig1 || !sig2) return NULL;
    int total = (sig1 ? sig1->n_switches : 0) + (sig2 ? sig2->n_switches : 0);
    DTA_SwitchingSignal* result = dta_signal_create(total);
    int i;
    for (i = 0; i < sig1->n_switches; i++)
        dta_signal_append(result, sig1->switch_times[i], sig1->mode_sequence[i]);
    double offset = sig1->switch_times[sig1->n_switches - 1];
    for (i = 0; i < sig2->n_switches; i++)
        dta_signal_append(result, sig2->switch_times[i] + offset,
                          sig2->mode_sequence[i]);
    return result;
}

DTA_SwitchingSignal* dta_signal_slice(const DTA_SwitchingSignal* sig,
                                       double t1, double t2) {
    if (!sig || t2 <= t1) return NULL;
    DTA_SwitchingSignal* result = dta_signal_create(sig->n_switches + 2);
    int active = dta_signal_active_mode(sig, t1);
    dta_signal_append(result, t1, active);
    int i;
    for (i = 0; i < sig->n_switches; i++) {
        double ts = sig->switch_times[i];
        if (ts > t1 && ts < t2)
            dta_signal_append(result, ts, sig->mode_sequence[i]);
    }
    active = dta_signal_active_mode(sig, t2);
    dta_signal_append(result, t2, active);
    return result;
}

void dta_signal_print(const DTA_SwitchingSignal* sig) {
    if (!sig) { printf("NULL signal\n"); return; }
    printf("Switching signal: %d switches, [%.4f, %.4f]\n",
           sig->n_switches - 1, sig->t_start, sig->t_end);
    int i;
    for (i = 0; i < sig->n_switches && i < 20; i++)
        printf("  t=%.4f -> mode %d\n", sig->switch_times[i], sig->mode_sequence[i]);
    if (sig->n_switches > 20) printf("  ... (%d more)\n", sig->n_switches - 20);
}
