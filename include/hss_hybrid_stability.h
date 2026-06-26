/**
 * @file hss_hybrid_stability.h
 * @brief Unified Stability Theory for Hybrid Switched Systems (L4-L7)
 *
 * This module provides a unified framework for analyzing stability
 * across all hybrid switched system formalisms. It bridges the
 * stability theories from:
 *
 *   - Switched Systems: CLF, MLF (Liberzon, Branicky)
 *   - Dwell-Time Analysis: constant & average dwell-time (Morse, Hespanha)
 *   - Hybrid Automata: LaSalle invariance, barrier certificates
 *   - Impulsive Systems: Lyapunov with jumps (Lakshmikantham)
 *   - Reset Control: Clegg Integrator, FORE stability (Beker/Hollot/Chait)
 *   - PWA Systems: Piecewise quadratic Lyapunov (Johansson/Rantzer)
 *   - Event-Triggered: ISS, Zeno-free stability (Tabuada)
 *   - Supervisory: Language-based stability (Ramadge/Wonham)
 *
 * Knowledge Points:
 *   L4 KP8: Unified stability theorem — common sufficient conditions
 *   L4 KP9: Comparison principle for hybrid systems
 *   L4 KP10: Converse Lyapunov theorem for switched systems
 *   L5 KP13: Cross-framework stability verification
 *   L5 KP14: Stability margin computation
 *   L5 KP15: Robust stability under perturbations
 *   L7 KP1: Power electronics stability (DC-DC converters)
 *   L7 KP2: Automotive control stability (engine, transmission)
 *   L7 KP3: HVAC control stability (building thermal)
 *
 * Course Mapping:
 *   MIT 6.241 — Nonlinear & hybrid stability
 *   Stanford CS359 — Verification of hybrid systems
 *   Berkeley EECS 291E — CPS stability
 *   CMU 15-424 — Stability proofs for CPS
 *   Cambridge Part II — Control theory
 *   Georgia Tech CS 6290 — Advanced architecture + control
 */

#ifndef HSS_HYBRID_STABILITY_H
#define HSS_HYBRID_STABILITY_H

#include "hss_core.h"

/* ============================================================================
 * L4 KP8: Unified Stability Framework
 * ============================================================================ */

/**
 * @brief Unified hybrid stability result (L4 KP8)
 *
 * Combines multiple stability analysis methods into a single framework.
 * The unified stability theorem states:
 *
 * A hybrid switched system is GUAS if there exists a family of
 * Lyapunov-like functions {V_q}_{q∈Q} and constants α₁, α₂, μ ≥ 1 such that:
 *   1. α₁||x||² ≤ V_q(x) ≤ α₂||x||²        (boundedness)
 *   2. V̇_q(x) ≤ -λ_q V_q(x) during flow    (decrease in mode)
 *   3. V_{q'}(x⁺) ≤ μ V_q(x⁻) at switches   (bounded jump)
 *   4. τ_d > ln(μ)/min_q(λ_q)               (sufficient dwell time)
 *
 * This unifies CLF (μ=1), MLF, dwell-time, and average dwell-time.
 */
typedef struct {
    /* Lyapunov function properties */
    double   alpha_1;          /**< Lower bound coefficient α₁           */
    double   alpha_2;          /**< Upper bound coefficient α₂           */
    double   lambda_min;       /**< Minimum decay rate across modes      */
    double   mu_max;           /**< Maximum overshoot at switches        */
    double   condition_number; /**< κ = α₂/α₁ (should be moderate)      */

    /* Dwell-time properties */
    double   tau_d;            /**< Actual/required dwell time           */
    double   tau_a;            /**< Average dwell time                   */
    bool     dwell_sufficient; /**< True if dwell condition satisfied    */

    /* Stability verdict */
    HSS_StabilityConcept stability_type; /**< Stability type achieved    */
    bool     is_stable;        /**< Overall stability                    */
    double   stability_margin; /**< Distance from instability boundary   */
    double   decay_rate;       /**< Guaranteed decay rate                */
    double   overshoot_bound;  /**< Maximum possible overshoot           */

    /* Cross-validation */
    bool     clf_verified;     /**< CLF condition verified               */
    bool     mlf_verified;     /**< MLF condition verified               */
    bool     dwell_verified;   /**< Dwell-time condition verified        */
    bool     barrier_verified; /**< Safety barrier verified              */

    char     diagnosis[256];   /**< Human-readable analysis summary      */
} HSS_UnifiedStability;

/**
 * @brief Comparison function for hybrid systems (L4 KP9)
 *
 * The comparison principle: if V̇ ≤ ω(V, t) for a comparison
 * function ω, and the comparison system u̇ = ω(u, t) is stable,
 * then the original system is stable.
 *
 * For hybrid systems, this extends to both continuous flow
 * and discrete jump comparisons.
 */
typedef struct {
    double   comparison_state;  /**< State u of comparison system        */
    double   (*omega_flow)(double u, double t, void *params);
                                /**< Comparison function during flow     */
    double   (*omega_jump)(double u, void *params);
                                /**< Comparison function at jumps        */
    void    *params;            /**< Parameters for comparison functions */
    bool     is_comparison_stable; /**< True if comparison system stable */
    double   bound_V;           /**< Upper bound V(x,t) ≤ u(t)          */
} HSS_ComparisonSystem;

/**
 * @brief Converse Lyapunov theorem result (L4 KP10)
 *
 * For GUES switched systems, a Lyapunov function always exists.
 * This constructs the converse Lyapunov function:
 *   V(x) = sup_{t≥0} ||φ(t, x)|| · e^{λ t}
 * which is Lipschitz, positive definite, and decreasing.
 */
typedef struct {
    double  *V_values;         /**< Converse LF values at sample points  */
    int      num_samples;      /**< Number of sample points              */
    bool     is_constructed;   /**< True if LF was constructed           */
    double   lambda;           /**< Exponential rate used                */
    double   lipschitz_const;  /**< Estimated Lipschitz constant          */
    double   max_overshoot;    /**< Maximum overshoot bound              */
} HSS_ConverseLyapunov;

/* ============================================================================
 * L5 KP13-KP15: Advanced Stability Algorithms
 * ============================================================================ */

/**
 * @brief Cross-framework stability result (L5 KP13)
 *
 * Verifies stability using multiple frameworks and consolidates results.
 * If one framework proves stability, the system is stable.
 * Consensus among frameworks increases confidence.
 */
typedef struct {
    HSS_UnifiedStability unified;     /**< Unified stability result       */

    /* Per-framework results */
    int     frameworks_checked;        /**< Number of frameworks checked  */
    int     frameworks_stable;         /**< Number proving stability      */
    char    **framework_names;        /**< Names of frameworks used       */
    bool    *framework_stable;        /**< Per-framework stability flag   */
    double  *framework_margin;        /**< Per-framework stability margin */

    /* Consensus */
    bool    consensus_stable;         /**< Majority vote: stable          */
    double  confidence;               /**< Confidence score [0,1]         */
} HSS_CrossFrameworkStability;

/**
 * @brief Stability margin computation result (L5 KP14)
 *
 * Computes the distance to instability for a switched system.
 * The stability margin is the smallest perturbation that
 * destabilizes the system.
 *
 * For linear modes: margin = min_q (-α(A_q)) under switching.
 * For nonlinear: margin from the LF decrease inequality.
 */
typedef struct {
    double   margin;            /**< Stability margin (positive = stable) */
    double   perturbation_norm; /**< Norm of destabilizing perturbation   */
    int      critical_mode;     /**< Mode with smallest margin            */
    double   critical_eigenvalue;/**< Most dangerous eigenvalue            */
    bool     robust;            /**< True if margin > tolerance           */
    double   robust_tolerance;  /**< Tolerance for robust margin          */
} HSS_StabilityMargin;

/**
 * @brief Robust stability under perturbation (L5 KP15)
 *
 * Analyzes stability when system matrices are perturbed:
 *   A_q → A_q + ΔA_q with ||ΔA_q|| ≤ δ
 *
 * Uses Lyapunov-based robustness: perturbation tolerated while
 * A_qᵀP + PA_q + 2δ||P||I < 0.
 */
typedef struct {
    double   max_perturbation;   /**< Max tolerable perturbation norm δ   */
    double  *perturbation_direction; /**< Worst-case direction (n×n)     */
    int      n;                  /**< State dimension                      */
    bool     has_stable_perturbation; /**< True if some perturbation OK  */
    double   degradation_rate;   /**< How fast margin degrades with δ    */
} HSS_RobustStability;

/* ============================================================================
 * L7 KP1-KP3: Application-Specific Stability
 * ============================================================================ */

/**
 * @brief Power electronics stability (L7 KP1) — DC-DC converter
 *
 * Verifies stability of a buck/boost converter under PWM switching.
 * Uses averaged model and switched Lyapunov analysis.
 */
typedef struct {
    double   L;                 /**< Inductance                          */
    double   C;                 /**< Capacitance                          */
    double   R;                 /**< Load resistance                      */
    double   Vin;               /**< Input voltage                        */
    double   Vref;              /**< Reference output voltage              */
    double   duty_cycle;        /**< Steady-state duty cycle              */
    bool     is_stable;         /**< Stability verdict                    */
    double   phase_margin;      /**< Phase margin from averaged model     */
    double   gain_margin;       /**< Gain margin from averaged model      */
    double   settling_time;     /**< Estimated settling time              */
    double   overshoot;         /**< Estimated overshoot percentage       */
} HSS_PowerElectronicsStability;

/**
 * @brief Automotive stability (L7 KP2) — Engine idle speed control
 *
 * Hybrid system: engine operating modes (idle, cruise, WOT)
 * with different dynamics in each mode.
 * Verifies stability across mode transitions.
 */
typedef struct {
    double   rpm;               /**< Current engine speed (RPM)           */
    double   rpm_target;        /**< Target idle speed (RPM)              */
    double   throttle_angle;    /**< Throttle position (%)                */
    double   load_torque;       /**< External load torque (Nm)            */
    int      gear;              /**< Current gear (-1=reverse, 0=neutral) */
    bool     is_idle_stable;    /**< Stable idle speed regulation         */
    double   rpm_overshoot;     /**< Max RPM overshoot at mode change     */
    double   recovery_time;     /**< Time to stabilize after disturbance  */
    double   disturbance_margin;/**< Max load torque before stall          */
} HSS_AutomotiveStability;

/**
 * @brief HVAC stability (L7 KP3) — Building thermal control
 *
 * Multi-zone thermal system with discrete ON/OFF heating/cooling.
 * Verifies temperature regulation stability with switching.
 */
typedef struct {
    double   zone_temp;         /**< Current zone temperature (°C)        */
    double   target_temp;       /**< Target temperature (°C)              */
    double   ambient_temp;      /**< Outside ambient temperature (°C)     */
    double   thermal_capacity;  /**< Building thermal mass (kJ/K)         */
    double   heating_power;     /**< Heater power (kW)                    */
    double   cooling_power;     /**< Cooler power (kW)                    */
    double   heat_loss_coeff;   /**< Thermal loss coefficient (kW/K)      */
    bool     heating_on;        /**< Heater state                         */
    bool     cooling_on;        /**< Cooler state                         */
    bool     is_temperature_stable; /**< Temperature within deadband      */
    double   temp_variation;    /**< Peak-to-peak temperature variation   */
    double   energy_consumption;/**< Energy consumption estimate           */
} HSS_HVACStability;

/* ============================================================================
 * L4-L7 API: Stability Functions
 * ============================================================================ */

/* ---- L4: Unified Stability Theorems ---- */

/**
 * @brief Compute unified stability analysis (L4 KP8)
 *
 * Applies all applicable stability theorems and returns
 * a consolidated stability verdict.
 *
 * Theory (unified): Combines CLF, MLF, dwell-time, and barrier
 * conditions into a single stability evaluation.
 *
 * @param sys HSS system
 * @param signal Switching signal (NULL for arbitrary switching)
 * @return Unified stability result
 */
HSS_UnifiedStability hss_unified_stability_analysis(
    const HSS_System *sys, const HSS_SwitchingSignal *signal);

/**
 * @brief Set up comparison system (L4 KP9)
 *
 * Creates a scalar comparison system u̇ = ω(u, t) that bounds
 * the Lyapunov function decay of the original hybrid system.
 *
 * @param sys HSS system
 * @param V_max Maximum initial Lyapunov value
 * @return Initialized comparison system
 */
HSS_ComparisonSystem hss_comparison_system_setup(
    const HSS_System *sys, double V_max);

/**
 * @brief Construct converse Lyapunov function (L4 KP10)
 *
 * For a GUES system, constructs the converse Lyapunov function
 * via simulation and supremal bound.
 *
 * Theorem (Mancilla-Aguilar & Garcia 2001): For switched systems
 * with GUES, there exists a smooth Lyapunov function.
 *
 * @param sys HSS system (must be GUES)
 * @param trace Execution trace from simulation
 * @param lambda Exponential rate
 * @return Converse Lyapunov result
 */
HSS_ConverseLyapunov hss_converse_lyapunov(
    const HSS_System *sys, const HSS_ExecutionTrace *trace,
    double lambda);

/* ---- L5: Advanced Stability Algorithms ---- */

/**
 * @brief Cross-framework stability verification (L5 KP13)
 *
 * Verifies stability using multiple theoretical frameworks
 * and provides a consensus assessment.
 *
 * @param sys HSS system
 * @param signal Switching signal
 * @return Cross-framework result
 */
HSS_CrossFrameworkStability hss_cross_framework_verify(
    const HSS_System *sys, const HSS_SwitchingSignal *signal);

/**
 * @brief Compute stability margin (L5 KP14)
 *
 * Finds the smallest perturbation that would destabilize the system.
 *
 * Algorithm: For each mode, compute the spectral abscissa α(A_q).
 * The stability margin is min_q (-α(A_q)) for Hurwitz A_q.
 * For nonlinear modes, uses the LF decay rate.
 *
 * @param sys HSS system
 * @return Stability margin result
 */
HSS_StabilityMargin hss_compute_stability_margin(const HSS_System *sys);

/**
 * @brief Analyze robust stability (L5 KP15)
 *
 * Determines the maximum tolerable perturbation that preserves
 * stability, using Lyapunov-based robustness analysis.
 *
 * Theorem (Khalil 2002, §9.3): If V̇ ≤ -λ||x||² and the perturbation
 * satisfies ||Δ|| ≤ ε for sufficiently small ε, then stability
 * is preserved.
 *
 * @param sys HSS system
 * @param perturbation_bound Maximum perturbation norm to consider
 * @param num_samples Number of perturbation samples
 * @return Robust stability result
 */
HSS_RobustStability hss_robust_stability_analysis(
    const HSS_System *sys, double perturbation_bound, int num_samples);

/* ---- L7: Application-Specific Stability ---- */

/**
 * @brief Analyze power electronics converter stability (L7 KP1)
 *
 * Evaluates stability of a DC-DC converter under PWM switching.
 *
 * @param L Inductance (H)
 * @param C Capacitance (F)
 * @param R Load resistance (Ω)
 * @param Vin Input voltage (V)
 * @param Vref Reference voltage (V)
 * @param duty_cycle PWM duty cycle
 * @return Power electronics stability result
 */
HSS_PowerElectronicsStability hss_power_electronics_stability(
    double L, double C, double R, double Vin, double Vref,
    double duty_cycle);

/**
 * @brief Analyze automotive engine stability (L7 KP2)
 *
 * Evaluates stability of engine idle speed control across modes.
 *
 * @param rpm_current Current engine speed
 * @param rpm_target Target idle speed
 * @param load_torque Current load torque
 * @param gear Current gear
 * @return Automotive stability result
 */
HSS_AutomotiveStability hss_automotive_stability(
    double rpm_current, double rpm_target,
    double load_torque, int gear);

/**
 * @brief Analyze HVAC building thermal stability (L7 KP3)
 *
 * Evaluates temperature regulation stability with switching.
 *
 * @param zone_temp Current temperature
 * @param target_temp Target temperature
 * @param ambient_temp Outdoor temperature
 * @param heat_capacity Building thermal mass
 * @param heat_power Heater capacity
 * @param cool_power Cooler capacity
 * @return HVAC stability result
 */
HSS_HVACStability hss_hvac_stability(
    double zone_temp, double target_temp, double ambient_temp,
    double heat_capacity, double heat_power, double cool_power);

/**
 * @brief Print unified stability report
 * @param result Stability result
 * @param fp Output stream
 */
void hss_unified_stability_report(const HSS_UnifiedStability *result,
                                   FILE *fp);

/**
 * @brief Print cross-framework stability report
 * @param result Cross-framework result
 * @param fp Output stream
 */
void hss_cross_framework_report(const HSS_CrossFrameworkStability *result,
                                 FILE *fp);

#endif /* HSS_HYBRID_STABILITY_H */
