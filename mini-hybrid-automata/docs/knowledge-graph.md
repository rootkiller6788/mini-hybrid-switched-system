# Knowledge Graph — mini-hybrid-automata

## L1: Definitions (Complete)

| # | Knowledge Point | C Implementation | Header |
|---|----------------|-----------------|--------|
| KP1 | Hybrid Automaton H = (Q,X,Init,Flow,Inv,E,G,R) | `hybrid_automaton_create` | hybrid_automaton.h |
| KP2 | Discrete mode/location q ∈ Q | `hybrid_automaton_add_mode` | hybrid_automaton.h |
| KP3 | Continuous variable x_i ∈ X | `hybrid_automaton_set_variable` | hybrid_automaton.h |
| KP4 | Discrete transition e ∈ E | `hybrid_automaton_add_transition` | hybrid_automaton.h |
| KP5 | Guard condition G(e) ⊆ ℝⁿ | `hybrid_guard_set` | hybrid_automaton.h |
| KP6 | Reset map R(e): ℝⁿ → ℝⁿ | `hybrid_reset_set` | hybrid_automaton.h |
| KP7 | Invariant condition Inv(q) ⊆ ℝⁿ | `hybrid_invariant_set` | hybrid_automaton.h |
| KP8 | Continuous flow f_q: ℝⁿ → ℝⁿ | `hybrid_flow_set` / `hybrid_flow_set_nonlinear` | hybrid_automaton.h |
| KP9 | Initial condition Init ⊆ Q × X | `hybrid_init_set` / `hybrid_init_set_rect` | hybrid_automaton.h |

## L2: Core Concepts (Complete)

| # | Knowledge Point | C Implementation | Header |
|---|----------------|-----------------|--------|
| KP10 | Hybrid time domain | `HybridTimeSet` struct | hybrid_execution.h |
| KP11 | Execution/trajectory | `hybrid_execution_create` | hybrid_execution.h |
| KP12 | Continuous flow segment | `hybrid_execution_append_flow` | hybrid_execution.h |
| KP13 | Discrete jump point | `hybrid_execution_append_jump` | hybrid_execution.h |
| KP14 | Determinism check | `hybrid_is_deterministic` | hybrid_execution.h |
| KP15 | Non-blocking check | `hybrid_is_nonblocking` | hybrid_execution.h |
| KP16 | Zeno execution detection | `hybrid_execution_is_zeno` | hybrid_execution.h |
| KP17 | Transition enabled check | `hybrid_transition_enabled` | hybrid_execution.h |
| KP18 | Invariant satisfaction | `hybrid_invariant_satisfied` | hybrid_execution.h |
| KP19 | Parallel composition H₁‖H₂ | `hybrid_compose_parallel` | hybrid_execution.h |

## L3: Mathematical Structures (Complete)

| # | Knowledge Point | C Implementation | Header |
|---|----------------|-----------------|--------|
| KP20 | Forward reachable set | `hybrid_reachable_forward` | hybrid_reachability.h |
| KP21 | Backward reachable set | `hybrid_reachable_backward` | hybrid_reachability.h |
| KP22 | Flowpipe computation | `hybrid_flowpipe_compute` | hybrid_reachability.h |
| KP23 | Zonotope representation | `hybrid_zono_create` / `*` | hybrid_reachability.h |
| KP24 | Polyhedron H-representation | `hybrid_poly_create` / `*` | hybrid_reachability.h |
| KP25 | Support function method | `hybrid_support_function` | hybrid_reachability.h |
| KP41 | Timed automaton | `hybrid_is_timed_automaton` | hybrid_subclass.h |
| KP42 | Region equivalence | `hybrid_region_compute` | hybrid_subclass.h |
| KP43 | Rectangular HA | `hybrid_is_rectangular` | hybrid_subclass.h |
| KP44 | Initialized rectangular HA | `hybrid_is_initialized_rectangular` | hybrid_subclass.h |
| KP45 | Linear HA (LHA) | `hybrid_is_linear_hybrid` | hybrid_subclass.h |
| KP46 | Piecewise affine (PWA) | `hybrid_is_pwa` | hybrid_subclass.h |
| KP47 | Mixed Logical Dynamical | `hybrid_to_mld` | hybrid_subclass.h |
| KP48 | O-minimal HA | `hybrid_is_ominimal` | hybrid_subclass.h |

## L4: Fundamental Theorems (Complete)

| # | Knowledge Point | C Implementation | Header |
|---|----------------|-----------------|--------|
| KP30 | Barrier certificate theorem | `hybrid_barrier_verify` | hybrid_safety.h |
| KP31 | Common Lyapunov function | `hybrid_common_lyapunov` | hybrid_safety.h |
| KP32 | Multiple Lyapunov functions | `hybrid_multiple_lyapunov` | hybrid_safety.h |
| KP33 | Inductive invariant | `hybrid_inductive_invariant_check` | hybrid_safety.h |
| KP34 | Dwell-time safety theorem | `hybrid_dwell_time_safety` | hybrid_safety.h |
| -- | Undecidability of reachability | Documented in code comments | hybrid_reachability.h |
| -- | Decidability for timed automata | Documented in code comments | hybrid_subclass.h |
| -- | Decidability for init. rectangular HA | Documented in code comments | hybrid_subclass.h |

## L5: Algorithms/Methods (Complete)

| # | Knowledge Point | C Implementation | Header |
|---|----------------|-----------------|--------|
| KP26 | Bisimulation quotient | `hybrid_bisimulation_quotient` | hybrid_reachability.h |
| KP27 | CEGAR loop | `hybrid_cegar` | hybrid_reachability.h |
| KP28 | Onion-peeling reachability | `hybrid_onion_peeling` | hybrid_reachability.h |
| KP36 | ODE integration (Euler, Heun, RK4) | `hybrid_ode_step` | hybrid_simulation.h |
| KP37 | Guard/invariant event detection | `hybrid_event_detect` | hybrid_simulation.h |
| KP38 | Precise event location | `hybrid_event_locate` | hybrid_simulation.h |
| KP39 | Hybrid simulation main loop | `hybrid_simulate` | hybrid_simulation.h |
| KP40 | Nondeterministic simulation | `hybrid_simulate_nondet` | hybrid_simulation.h |

## L6: Canonical Problems (Complete)

| # | Knowledge Point | C Implementation |
|---|----------------|-----------------|
| KP49 | Bouncing ball | `example_bouncing_ball` |
| KP50 | Thermostat | `example_thermostat` |
| KP51 | Two-tank system | `example_two_tank` |
| KP52 | Train-gate controller | `example_train_gate` |
| KP53 | DC-DC converter | `example_dcdc_converter` |

## L7: Applications (Complete — 3 applications)

| # | Application | C Implementation | Domain |
|---|------------|-----------------|--------|
| KP54 | Engine air-fuel ratio control | `example_engine_afr_control` | Automotive |
| KP55 | Robot obstacle avoidance | `example_robot_obstacle_avoidance` | Robotics |
| KP56 | Medical infusion pump | `example_infusion_pump` | Medical devices |

## L8: Advanced Topics (Partial)

| Topic | Status | Notes |
|-------|--------|-------|
| Stochastic hybrid automata | Documented | Not implemented (requires probability libraries) |
| Compositional verification | Partial | Parallel composition implemented |
| CEGAR with predicate abstraction | Partial | CEGAR loop with coarse abstraction |

## L9: Research Frontiers (Partial)

| Topic | Status | Notes |
|-------|--------|-------|
| Symbolic reachability | Documented | Zonotope/polyhedral abstractions as building blocks |
| Learning-based hybrid control | Documented | Architecture supports it |
| Meta-complexity of HA verification | Documented | Undecidability results referenced |
