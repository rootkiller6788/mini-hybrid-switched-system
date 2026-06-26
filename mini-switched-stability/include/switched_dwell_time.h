#ifndef SWITCHED_DWELL_TIME_H
#define SWITCHED_DWELL_TIME_H

#include "switched_types.h"

/* ============================================================================
 * Dwell-Time Analysis for Switched Systems (L4, L5)
 *
 * Implements:
 *   - Minimum dwell time computation (Morse 1996)
 *   - Average dwell time computation (Hespanha & Morse 1999)
 *   - Multiple Lyapunov function dwell-time analysis
 *   - Switching signal dwell-time verification
 * ============================================================================ */

/**
 * Compute minimum dwell time tau_d from system parameters.
 *
 * Given stability margin lambda_0 > 0 and Lyapunov mismatch mu >= 1:
 *   tau_d = ln(mu) / lambda_0
 *
 * If the switching signal satisfies t_{k+1} - t_k >= tau_d for all k,
 * the switched system is GUES with decay rate:
 *   ||x(t)|| <= K * exp(-(lambda_0 - ln(mu)/tau_d) * t) * ||x(0)||
 *
 * Reference: Liberzon (2003), Section 3.2.1, Theorem 3.1.
 *
 * @param lambda_0  Stability margin (max decay rate)
 * @param mu        Lyapunov function mismatch parameter
 * @return          Minimum dwell time tau_d
 */
double sdt_compute_min_dwell(double lambda_0, double mu);

/**
 * Compute average dwell time tau_a* from system parameters.
 *
 * Given lambda_0 > 0, mu >= 1:
 *   tau_a* = ln(mu) / lambda_0
 *
 * For any tau_a > tau_a*, the switched system is GUES with:
 *   ||x(t)|| <= K * N_0 * exp(-(lambda_0 - ln(mu)/tau_a) * t) * ||x(0)||
 *
 * Reference: Hespanha & Morse (1999), Liberzon (2003) Theorem 3.2.
 *
 * @param lambda_0  Stability margin
 * @param mu        Lyapunov mismatch parameter
 * @return          Critical average dwell time tau_a*
 */
double sdt_compute_avg_dwell(double lambda_0, double mu);

/**
 * Compute the actual average dwell time from a switching signal.
 *
 * tau_avg = (t_N - t_0) / N_s  where N_s = number of switches.
 *
 * @param signal  Switching signal history
 * @return        Actual average dwell time
 */
double sdt_actual_avg_dwell(const SwitchingSignal *signal);

/**
 * Compute the chatter bound N_0 for average dwell time analysis.
 *
 * N_0 is the maximum number of switches that can occur in
 * an interval of length tau_a without violating stability.
 *
 * Typically N_0 = 0 or 1 (0 = no fast switching allowed).
 */
int sdt_chatter_bound(const SwitchingSignal *signal);

/**
 * Verify that a switching signal satisfies the minimum dwell time constraint.
 *
 * Checks: t_{k+1} - t_k >= tau_d for all k = 0, ..., N-1.
 *
 * @param switch_times  Array of switch instants (size n_switches+1)
 * @param n_switches    Number of switch intervals
 * @param tau_d         Required minimum dwell time
 * @return              true if constraint satisfied
 */
bool sdt_check_dwell(const double *switch_times, int n_switches, double tau_d);

/**
 * Verify the average dwell time condition.
 *
 * For all T > 0:
 *   N_sigma(T, 0) <= N_0 + T / tau_a
 *
 * where N_sigma(T, 0) counts switches in [0, T].
 *
 * @param signal  Switching signal
 * @param tau_a   Average dwell time parameter
 * @param N0      Chatter bound
 * @return        true if average dwell time condition holds for all T
 */
bool sdt_check_avg_dwell(const SwitchingSignal *signal, double tau_a, int N0);

/**
 * Perform full dwell-time analysis on a switched system.
 *
 * This function:
 *   1. Computes stability margin lambda_0 for all subsystems
 *   2. Solves Lyapunov equations for each Hurwitz A_i
 *   3. Computes mismatch parameter mu
 *   4. Computes required tau_d and tau_a*
 *   5. Evaluates whether the actual switching signal is slow enough
 *
 * @param sys  Switched system
 * @param dta  Output: dwell time analysis results
 */
void sdt_full_analysis(SwitchedSystem *sys, DwellTimeAnalysis *dta);

/**
 * Compute stability margin lambda_0 = min_i [-max_j Re(lambda_j(A_i))].
 *
 * For each subsystem A_i:
 *   alpha_i = -max_j Re(lambda_j(A_i))
 * Then: lambda_0 = min_i alpha_i
 *
 * lambda_0 > 0 iff all subsystems are Hurwitz.
 *
 * @param sys  Switched system
 * @return     Stability margin (positive if all modes stable)
 */
double sdt_stability_margin(const SwitchedSystem *sys);

/**
 * Compute Lyapunov mismatch parameter mu.
 *
 * mu = max_{i,j} lambda_max(P_j) / lambda_min(P_i)
 *
 * This bounds the ratio between different Lyapunov function values
 * at switching instants. mu = 1 means perfect matching (ideal case).
 *
 * @param mlf  Multiple Lyapunov functions
 * @return     Mismatch parameter mu >= 1
 */
double sdt_lyap_mismatch(const MultipleLyapunovFunctions *mlf);

/**
 * Evaluate the dwell-time stability certificate.
 *
 * Given tau_d, lambda_0, mu, and the actual switching signal:
 *   - If tau_actual >= tau_d for all intervals -> stable
 *   - If tau_avg > tau_a* -> stable (average dwell time)
 *   - Otherwise: may be unstable under this switching
 *
 * @param dta     Dwell time analysis structure
 * @param signal  Actual switching signal used
 * @return        Stability classification
 */
SwitchedStabilityType sdt_certify(const DwellTimeAnalysis *dta, const SwitchingSignal *signal);

/**
 * Find the largest tau_d that guarantees stability.
 * This is the "slow switching" stability margin.
 *
 * Binary search over tau_d in [0, T_max] to find the boundary
 * at which stability transitions from unstable to stable.
 *
 * @param sys    Switched system
 * @param T_max  Maximum time to consider for search
 * @return       Critical dwell time tau_d_crit
 */
double sdt_critical_dwell_time(SwitchedSystem *sys, double T_max);

/**
 * Scale dwell time to account for different initial conditions.
 *
 * Given a system with stability margin lambda_0 and initial
 * condition magnitude ||x_0||, compute the adjusted dwell time
 * that guarantees ||x(T)|| <= eps * ||x_0|| for a given tolerance eps.
 *
 * @param lambda_0   Stability margin
 * @param eps        Desired attenuation factor
 * @param mu         Lyapunov mismatch parameter
 * @return           Adjusted minimum dwell time
 */
double sdt_dwell_for_tolerance(double lambda_0, double eps, double mu);

#endif /* SWITCHED_DWELL_TIME_H */
