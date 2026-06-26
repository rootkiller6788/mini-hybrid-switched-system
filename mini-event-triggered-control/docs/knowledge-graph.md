# Knowledge Graph — Event-Triggered Control

## L1: Definitions (Complete ✅)

| # | Definition | C Representation | Lean Representation |
|---|-----------|-----------------|---------------------|
| 1.1 | Event-Triggered Control System | `ETCSystem` struct | — |
| 1.2 | Trigger Function Γ(x,e) | `etc_trigger_static/quadratic/...` | `trigger_static/quadratic/...` |
| 1.3 | Measurement Error e(t) = x̂(t_k) − x(t) | `ETCVector e` in ETCSystem | — |
| 1.4 | Minimum Inter-Event Time (MIET) | `etc_inter_event_time_bound()` | `miet_bound` |
| 1.5 | Zeno Behavior | `ETC_REGIME_ZENO` enum | `ETCRegime.zeno` |
| 1.6 | ISS Lyapunov Function | `ETCLyapunovFunction` struct | `ISSLyapunovData` |
| 1.7 | Comparison Functions (K, K∞, KL) | `ETCComparisonClass` enum | `ComparisonClass.K/Kinf/KL` |
| 1.8 | Trigger Type Classification | `ETCTriggerType` enum | `TriggerType` inductive type |

## L2: Core Concepts (Complete ✅)

| # | Concept | Implementation |
|---|---------|---------------|
| 2.1 | Event-triggered vs time-triggered sampling | `ETCSamplingType` (Riemann/Lebesgue) |
| 2.2 | State-dependent triggering | `etc_trigger_static`, `etc_trigger_quadratic` |
| 2.3 | Dynamic triggering with internal variable η | `etc_trigger_dynamic`, `etc_system_step_dynamic_trigger` |
| 2.4 | Self-triggered control (no monitoring) | `etc_self.c` entire module |
| 2.5 | Periodic event-triggered control (PETC) | `etc_periodic.c` entire module |
| 2.6 | Zero-order hold control update | `etc_system_compute_control` |
| 2.7 | Communication reduction | `etc_communication_ratio` |
| 2.8 | Trigger margin monitoring | `etc_adaptive_step` |

## L3: Mathematical Structures (Complete ✅)

| # | Structure | Implementation |
|---|-----------|---------------|
| 3.1 | State-space (A, B, K) for ETC | `ETCSystem.A/B/K/Acl` matrices |
| 3.2 | Quadratic Lyapunov function V(x)=xᵀPx | `ETCLyapunovFunction` with P matrix, eigenvalues |
| 3.3 | Continuous-time Lyapunov equation AᵀP+PA=−Q | `etc_solve_lyapunov` |
| 3.4 | Matrix exponential e^{At} | `etc_matrix_exponential` (scaling+squaring) |
| 3.5 | Flow map Φ_t(x₀,u) | `etc_flow_map` |
| 3.6 | ISS gain characterization (α₁,α₂,α₃,γ) | `etc_verify_iss_lyapunov` |
| 3.7 | K-class function verification | `etc_verify_kclass` |
| 3.8 | KL-class exponential decay | `etc_kl_function` |

## L4: Fundamental Laws/Theorems (Complete ✅)

| # | Theorem | Code |
|---|---------|------|
| 4.1 | **Tabuada (2007) Theorem**: MIET lower bound τ_min = σ / (||Acl|| + σ||BK||) | `etc_compute_iet_lower_bound` |
| 4.2 | **Stability condition**: σ < σ_max = λ_min(Q)/(2||PBK||) | `etc_is_sigma_stabilizing`, `etc_compute_sigma_max` |
| 4.3 | **ISS Lyapunov characterization**: V̇ ≤ −α₃(|x|) + γ(|e|) | `etc_verify_iss_lyapunov` |
| 4.4 | **Zeno-freeness condition**: strictly positive MIET ⇒ Zeno-free | `etc_check_zeno_free` |
| 4.5 | **PETC stability** (Heemels 2013): h ≤ h_max preserves stability | `etc_petc_max_sampling_period` |
| 4.6 | **STC conservatism**: τ_k ≤ τ_ETC(x_k) | `etc_stc_worst_case_interval` |
| 4.7 | **Dynamic trigger positivity**: η(t) ≥ 0 for all t | `etc_system_step_dynamic_trigger` (clamped) |
| 4.8 | **MIET positivity theorem** | `miet_positive` in Lean |
| 4.9 | **Trigger type exhaustiveness** | `trigger_types_exhaustive` in Lean |
| 4.10 | **Regime classification exhaustiveness** | `regimes_exhaustive` in Lean |

## L5: Algorithms/Methods (Complete ✅)

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 5.1 | RK4 integration with event detection | `etc_system_step_rk4` |
| 5.2 | Event-triggered simulation loop | `etc_system_simulate`, `etc_simulate_full` |
| 5.3 | Trigger zero-crossing detection | `etc_detect_event` |
| 5.4 | Lyapunov equation iterative solver | `etc_solve_lyapunov` |
| 5.5 | Matrix exponential (scaling+squaring) | `etc_matrix_exponential` |
| 5.6 | Flow map computation (augmented system) | `etc_flow_map` |
| 5.7 | Power iteration for spectral radius | `etc_matrix_spectral_radius` |
| 5.8 | Adaptive step-size control | `etc_adaptive_step` |
| 5.9 | PETC simulation with periodic checks | `etc_petc_simulate` |
| 5.10 | STC next-interval bisection search | `etc_stc_next_interval` |
| 5.11 | Threshold design from performance specs | `etc_design_threshold` |
| 5.12 | Gershgorin eigenvalue estimation | `etc_matrix_is_positive_definite` |

## L6: Canonical Problems (Complete ✅)

| # | Problem | Example |
|---|---------|---------|
| 6.1 | First-order stabilization with ETC | `example1_first_order.c` |
| 6.2 | Double integrator with PD control | `example2_double_integrator.c` |
| 6.3 | Networked control: ETC vs PETC vs STC | `example3_networked_control.c` |

## L7: Applications (Partial+ ✅)

| # | Application | Implementation |
|---|-------------|---------------|
| 7.1 | Networked embedded control with bandwidth constraints | `example3_networked_control.c` (Tesla-like wireless sensor/actuator) |
| 7.2 | Battery-aware IoT sensor nodes | Communication ratio analysis in simulation engine |
| 7.3 | Wireless sensor/actuator networks | PETC/STC comparison for resource-constrained nodes |

## L8: Advanced Topics (Partial+ ✅)

| # | Topic | Implementation |
|---|-------|---------------|
| 8.1 | Dynamic event-triggered control (Girard 2015) | `etc_trigger_dynamic`, dynamic η evolution |
| 8.2 | Self-triggered control with state prediction | `etc_stc_next_interval`, `etc_stc_simulate` |
| 8.3 | Robust STC under bounded disturbances | `etc_stc_robust_interval` |
| 8.4 | Adaptive STC with online τ_min adjustment | `etc_stc_adapt` |
| 8.5 | Mixed trigger with time-decaying threshold | `etc_trigger_mixed` |
| 8.6 | L₂-gain bound for ETC under disturbances | `etc_compute_l2_gain` |

## L9: Research Frontiers (Partial ✅ — documented)

| # | Frontier | Status |
|---|----------|--------|
| 9.1 | Learning-based event-triggering (RL for σ selection) | Documented in course-tree.md |
| 9.2 | Security of event-triggered systems (DoS attacks) | Documented in gap-report.md |
| 9.3 | Decentralized event-triggered multi-agent consensus | Documented in course-alignment.md |
| 9.4 | Event-triggered MPC | Documented |
| 9.5 | Quantum event-triggered measurement | Documented |
