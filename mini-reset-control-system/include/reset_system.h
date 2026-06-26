/* reset_system.h - Reset Control System Composition and Analysis
 *
 * Functions for composing reset control elements into feedback loops,
 * analyzing closed-loop stability, computing describing functions,
 * and evaluating hybrid system trajectories.
 *
 * Knowledge coverage:
 *   L2: Reset feedback composition, closed-loop reset system
 *   L3: Interconnection of hybrid systems, series/parallel/feedback
 *   L4: Circle criterion for reset systems, passivity-based stability
 *   L5: Numerical frequency-response computation, H2/Hinf norms
 *
 * Ref: [BB12] Ch.4-6; [NZT08] on stability; Nešić, Teel, Zaccarian (2011)
 *      "Stability and performance of SISO control systems with
 *       First-Order Reset Elements", IEEE TAC
 */

#ifndef RESET_SYSTEM_H
#define RESET_SYSTEM_H

#include "reset_core.h"
#include "reset_element.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L3: System Composition Structures
 * ================================================================ */

/** ResetFeedbackLoop - standard unity-feedback reset control system.
 *
 *  Architecture:
 *    r(t) --->(+)--->[ Reset Ctrl C ]--->[ Plant P ]---+---> y(t)
 *             ^-                                      |
 *             |_______________________________________|
 *
 *  The reset controller C has internal reset logic (zero-crossing
 *  of error e = r - y triggers state reset in C).
 *
 *  The plant P is a standard linear time-invariant system.
 *
 *  Ref: [BB12] Section 4.1, Figure 4.1
 */
typedef struct {
    ResetSystem     *controller;  /**< reset controller (Clegg/FORE/PID) */
    ResetLinearBase *plant;       /**< LTI plant (can be any order)     */
    double           r;           /**< current reference input          */
    double           e;           /**< current tracking error           */
    double           y;           /**< current plant output             */
    double           u;           /**< current control signal           */
    double          *xp;          /**< plant state vector               */
    double          *xc;          /**< controller state (alias to rsys) */
    int              np;          /**< plant state dimension            */
    int              nc;          /**< controller state dimension       */
    double           dt;          /**< simulation time step             */
    bool             is_initialized;
} ResetFeedbackLoop;

/** ResetSeries - series connection of reset elements.
 *  G(s) = G2_reset(s) * G1(s)  (G1 first, then G2).
 *  Useful for designing higher-order reset compensators. */
typedef struct {
    ResetSystem     *elem1;       /**< first reset element in series   */
    ResetSystem     *elem2;       /**< second reset element in series  */
    double           u;           /**< input to first element          */
    double           y1;          /**< output of first element         */
    double           y;           /**< output of second element        */
    double           e_prev1;     /**< previous error for elem1        */
    double           e_prev2;     /**< previous error for elem2        */
} ResetSeries;

/** ResetParallel - parallel connection of reset elements.
 *  G(s) = G1_reset(s) + G2(s). */
typedef struct {
    ResetSystem     *elem1;       /**< first reset element             */
    ResetLinearBase *elem2;       /**< second (linear) element         */
    double           u;           /**< common input                    */
    double           y1;          /**< output of first element         */
    double           y2;          /**< output of second element        */
    double           y;           /**< total output                    */
    double           e_prev1;     /**< previous input for elem1        */
} ResetParallel;

/* ================================================================
 * L4: Stability Analysis Functions
 * ================================================================ */

/** Compute the eigenvalues of the base linear system matrix A.
 *  Uses QR algorithm via Francis double-shift method.
 *  eigen_real, eigen_imag: pre-allocated arrays of length n.
 *  Complexity: O(n^3).
 *  Returns: number of eigenvalues found (should equal n). */
int reset_eigenvalues(const ResetLinearBase *base,
                       double *eigen_real, double *eigen_imag);

/** Check if the base linear system is Hurwitz stable
 *  (all eigenvalues of A have strictly negative real part).
 *  Complexity: O(n^3). Returns true if Hurwitz stable. */
bool reset_is_hurwitz(const ResetLinearBase *base);

/** Compute the H2 norm of the linear base system.
 *  ||G||_2^2 = trace(C * P * C^T) where P solves
 *  A*P + P*A^T + B*B^T = 0 (controllability Gramian).
 *  Complexity: O(n^3) for Lyapunov equation.
 *  Returns infinity if A is not Hurwitz. */
double reset_h2_norm(const ResetLinearBase *base);

/** Compute the H-infinity norm via bisection on Hamiltonian matrix.
 *  ||G||_inf = sup_omega sigma_max(G(j omega)).
 *  Complexity: O(n^3 * iterations).
 *  max_iter: maximum bisection iterations.
 *  tol: convergence tolerance for the bisection. */
double reset_hinf_norm(const ResetLinearBase *base, int max_iter, double tol);

/** Compute the closed-loop A-matrix for unity feedback:
 *  A_cl = [Ap - Bp*Dc*Cp,  Bp*Cc  ;
 *          -Bc*Cp,          Ac     ]
 *  where (Ac,Bc,Cc,Dc) = controller, (Ap,Bp,Cp,Dp) = plant.
 *  A_cl is allocated and returned; caller must free() it.
 *  Complexity: O((np+nc)^3) for matrix multiplication. */
double* reset_closed_loop_A(const ResetLinearBase *plant,
                             const ResetLinearBase *controller,
                             int *n_cl);

/** Check closed-loop stability of the reset feedback system
 *  using the H_beta condition (Banos & Barreiro, 2012).
 *
 *  The reset system with FORE having reset ratio rho is
 *  asymptotically stable if:
 *    1. The base linear closed-loop is stable
 *    2. There exists beta > 0 and P > 0 such that
 *       H_beta(P) < 0  (LMI condition)
 *
 *  This function checks condition (1) only (linear base stability).
 *  For the full LMI-based condition, use reset_lyapunov functions.
 *
 *  Returns true if base closed-loop is Hurwitz stable. */
bool reset_check_hbeta_stability(const ResetLinearBase *plant,
                                  const ResetSystem *controller);

/* ================================================================
 * L3/L5: System Evaluation Functions
 * ================================================================ */

/* ---- Feedback Loop Simulation ---- */

/** Initialize a reset feedback loop structure.
 *  Allocates plant state vector and sets initial conditions to zero. */
ResetFeedbackLoop* reset_feedback_create(const ResetSystem *controller,
                                          const ResetLinearBase *plant);

/** Free a reset feedback loop and all internal allocations. */
void reset_feedback_free(ResetFeedbackLoop *loop);

/** Execute one simulation step of the reset feedback loop.
 *  dt: time step duration.
 *  r: reference input at current time.
 *  Handles: error computation, controller step (with ZC reset),
 *  plant state integration (Euler), and output computation.
 *  Returns the plant output y. */
double reset_feedback_step(ResetFeedbackLoop *loop, double dt, double r);

/** Run the reset feedback loop for multiple steps.
 *  t_span: array of time points [n_steps].
 *  r_vals: reference input at each time point [n_steps].
 *  y_out: output array to be filled [n_steps], pre-allocated.
 *  n_steps: number of simulation steps.
 *  Returns number of reset events that occurred during simulation. */
int reset_feedback_simulate(ResetFeedbackLoop *loop,
                             const double *t_span,
                             const double *r_vals,
                             double *y_out,
                             int n_steps);

/* ---- Frequency Response of Reset Systems ---- */

/** Compute the closed-loop frequency response magnitude
 *  at a given frequency omega (rad/s) for the base linear system
 *  (ignoring resets).
 *  Complexity: O(n^3) for matrix inversion. */
double reset_closed_loop_freqresp(const ResetLinearBase *plant,
                                   const ResetSystem *controller,
                                   double omega);

/** Compute the sensitivity function S(jw) = 1/(1 + P(jw)*C(jw))
 *  at frequency omega for the base linear system.
 *  Returns complex magnitude of sensitivity.*/
double reset_sensitivity(const ResetLinearBase *plant,
                          const ResetSystem *controller,
                          double omega);

/** Compute the complementary sensitivity T(jw) = P(jw)*C(jw)/(1 + P(jw)*C(jw))
 *  at frequency omega for the base linear system.
 *  Returns complex magnitude of complementary sensitivity. */
double reset_complementary_sensitivity(const ResetLinearBase *plant,
                                        const ResetSystem *controller,
                                        double omega);

/* ---- Reset System Interconnection ---- */

/** Create a reset system representing the series connection
 *  of a linear system followed by a reset controller.
 *  G_series = C_reset(s) * P(s).
 *  The reset condition of the original controller is preserved. */
ResetSystem* reset_series_connection(const ResetLinearBase *sys1,
                                      const ResetSystem *sys2);

/** Create a reset system representing the feedback connection
 *  of a plant and reset controller.
 *  G_fb = P(s) / (1 + P(s) * C_reset(s)) (base linear closed-loop).
 *  Only the linear base of the closed-loop is returned;
 *  reset conditions must be handled separately. */
ResetLinearBase* reset_feedback_connection(const ResetLinearBase *plant,
                                            const ResetSystem *controller);

#ifdef __cplusplus
}
#endif

#endif /* RESET_SYSTEM_H */