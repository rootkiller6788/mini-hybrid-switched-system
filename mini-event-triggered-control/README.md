# Event-Triggered Control (ETC)

## Core Concept

Event-triggered control replaces periodic sampling with state-dependent triggering:
control updates occur only when a **triggering condition** Γ(x, e) ≥ 0 is violated.
This reduces communication and computation while preserving stability.

```
Time-triggered:    |--*--*--*--*--*--*--*--*--*--|  (fixed period T_s)
Event-triggered:   |--*-----*--------*---*-----*--|  (state-dependent)
```

## Key Equations

### Trigger Function (Tabuada, 2007)
```
Γ(x, e) = |e| − σ|x|        (static, relative threshold)
Γ(x, e) = |e|² − σ|x|²      (quadratic, energy-based)
Γ(x, e) = η + θ(σ|x|²−|e|²) (dynamic, Girard 2015)
```
where e(t) = x(t_k) − x(t) is the measurement error.

### Minimum Inter-Event Time
```
τ_min = σ / (||Acl|| + σ||BK||)   (Tabuada 2007, Lemma 4)
```
A strictly positive τ_min guarantees Zeno-freeness.

### Stability Condition
```
σ < σ_max = λ_min(Q) / (2 ||PBK||)
```
Any σ ∈ (0, σ_max) ensures asymptotic stability of the ETC system.

### ISS Lyapunov Characterization
```
α₁(|x|) ≤ V(x) ≤ α₂(|x|)
V̇(x) ≤ −α₃(|x|) + γ(|e|)
```
With V(x) = xᵀPx, α₁(s) = λ_min(P)·s², α₃(s) = (λ_min(Q) − ||PBK||)·s².

## Core Definitions

| Definition | Type | Description |
|-----------|------|-------------|
| ETCSystem | struct | Complete event-triggered control system (A,B,K, state, trigger, Lyapunov) |
| ETCVector / ETCMatrix | struct | N-dimensional vector and matrix types |
| ETCEvent / ETCEventHistory | struct | Event record and history with statistics |
| ETCLyapunovFunction | struct | Quadratic Lyapunov function V(x) = xᵀPx |
| ETCTriggerType | enum | Static, dynamic, mixed, absolute, self-triggered |
| ETCRegime | enum | Asymptotically stable, ISS practical, Zeno, minimal interval |

## Core Theorems

| Theorem | Formula | Source |
|---------|---------|--------|
| MIET lower bound | τ_min = σ / (||Acl|| + σ||BK||) | Tabuada (2007) Lemma 4 |
| Stability threshold | σ_max = λ_min(Q) / (2||PBK||) | Tabuada (2007) Theorem 1 |
| ISS dissipation | V̇ ≤ −(λ_min(Q)−||PBK||)|x|² + ||PBK|||e|² | Sontag (2008) |
| PETC h_max | h_max = (1/||Acl||) ln(1 + σ||Acl||/(||Acl||+σ||BK||)) | Heemels (2013) Thm III.2 |
| STC worst-case | τ_wc = σ / (||Acl|| + σ||BK||) | Mazo (2010) |
| Dynamic η positivity | η(t) ≥ 0 for all t | Girard (2015) Lemma 1 |

## Core Algorithms

| Algorithm | Complexity | Implementation |
|-----------|-----------|----------------|
| RK4 integration with event detection | O(n²) per step | `etc_system_step_rk4` |
| Trigger zero-crossing (linear interp) | O(n) | `etc_detect_event` |
| Matrix exponential (scaling+squaring) | O(n³·log(||A||t)) | `etc_matrix_exponential` |
| Lyapunov equation solver (iterative) | O(n³·iters) | `etc_solve_lyapunov` |
| Flow map (augmented system) | O((n+1)³) | `etc_flow_map` |
| PETC simulation | O(T/dt · n²) | `etc_petc_simulate` |
| STC next-interval bisection | O(n³·log(1/ε)) | `etc_stc_next_interval` |

## Classic Problems

| Problem | Example File | Description |
|---------|-------------|-------------|
| First-order ETC | `example1_first_order.c` | σ vs event-count tradeoff, communication savings |
| Double integrator | `example2_double_integrator.c` | CLF, ISS verification, event timeline visualization |
| Networked control | `example3_networked_control.c` | ETC vs PETC vs STC comparison, bandwidth analysis |

## Nine-School Course Mapping

| School | Course | Topic Link |
|--------|--------|------------|
| MIT | 6.241J Dynamic Systems | Lyapunov stability for ETC |
| MIT | 16.323 Optimal Ctrl | Optimal event scheduling |
| Stanford | AA203 Optimal Ctrl | Optimization-based triggering |
| Stanford | EE363 Convex Opt | LMI conditions for ETC |
| Berkeley | EE221A Linear Systems | Lyapunov equation, Hurwitz analysis |
| Berkeley | EE222 Nonlinear Systems | ISS, comparison functions |
| CMU | 24-677 Nonlinear Ctrl | Dynamic triggering, nonlinear ETC |
| Princeton | MAE 546 Optimal Ctrl | State-feedback design |
| Cambridge | 4F3 Nonlinear Ctrl | Zeno analysis |

## Module Status: COMPLETE ✅

- **L1** Definitions: Complete — 6 typedef structs, 5 enums, 8 core definitions with Lean types
- **L2** Core Concepts: Complete — 6 headers, 6 source files, 8 concepts implemented
- **L3** Math Structures: Complete — Matrix/vector library, Lyapunov eq solver, matrix exponential, flow map
- **L4** Fundamental Laws: Complete — 10 theorems (Tabuada MIET, ISS, Zeno-freeness), 62 assert tests, 15+ Lean theorems
- **L5** Algorithms: Complete — 12 algorithms (RK4, event detection, matrix exp, Lyapunov solver, PETC/STC sim)
- **L6** Canonical Problems: Complete — 3 examples (>30 lines, main+printf)
- **L7** Applications: Partial+ — Networked control, IoT sensor nodes, bandwidth-constrained systems
- **L8** Advanced Topics: Partial+ — Dynamic ETC, STC, robust STC, mixed triggers, L₂-gain analysis
- **L9** Research Frontiers: Partial — 5 frontiers documented (learning-based ETC, security, decentralized, MPC, quantum)

## Build & Test

```bash
make          # Build static library libevent_triggered_control.a
make test     # Build and run test suite (62 asserts)
make examples # Build all 3 examples
make demo     # Build and run demo
make bench    # Build and run benchmarks
make clean    # Clean build artifacts
```

## Directory Structure

```
mini-event-triggered-control/
  Makefile                    # Build system
  README.md                   # This file
  include/
    etc_core.h                # Core types, ETCSystem, vectors, matrices (325 lines)
    etc_trigger.h             # Trigger functions, threshold design (185 lines)
    etc_stability.h           # ISS, Lyapunov, Zeno, L2-gain analysis (185 lines)
    etc_dynamics.h            # Simulation engine, flow maps, IET analysis (201 lines)
    etc_periodic.h            # Periodic ETC (PETC) (108 lines)
    etc_self.h                # Self-triggered control (STC) (149 lines)
  src/
    etc_core.c                # Core implementation: vectors, matrices, ETCSystem (673 lines)
    etc_trigger.c             # Trigger evaluation, sigma_max, threshold design (197 lines)
    etc_stability.c           # Lyapunov solver, ISS, Zeno, L2-gain, practical stability (425 lines)
    etc_dynamics.c            # Simulation, matrix exp, flow map, adaptive steps (475 lines)
    etc_periodic.c            # PETC implementation with stability guarantees (189 lines)
    etc_self.c                # STC with prediction, robustness, adaptation (299 lines)
    event_triggered.lean      # Lean 4: trigger types, regimes, MIET, ISS theorems (254 lines)
  tests/
    test_event_triggered.c    # 62 asserts covering all APIs
  examples/
    example1_first_order.c          # σ sweep, communication savings
    example2_double_integrator.c    # CLF, ISS, event timeline
    example3_networked_control.c    # ETC vs PETC vs STC comparison
  demos/
    demo_etc_overview.c       # Interactive concept demonstration
  benches/
    bench_trigger.c           # Trigger evaluation, simulation throughput
  docs/
    knowledge-graph.md        # L1-L9 knowledge coverage table
    coverage-report.md        # Per-level assessment
    gap-report.md             # Missing knowledge items
    course-alignment.md       # Nine-school curriculum mapping
    course-tree.md            # Prerequisites and dependency graph
```

## Quality Metrics

| Metric | Value | Requirement |
|--------|-------|-------------|
| include/ .h files | 6 | ≥ 4 |
| src/ .c files | 6 | ≥ 4 |
| src/ .lean files | 1 | ≥ 1 |
| include/ + src/ lines | 3,411 | ≥ 3,000 |
| typedef struct count | 16 | ≥ 5 |
| Exported functions | 70+ | ≥ 20 |
| Test asserts | 62 | ≥ 15 |
| Examples | 3 | ≥ 3 |
| Docs | 5 | ≥ 5 |
| Lean theorems | 15+ | ≥ 1 |
| make compiles | YES | required |
| make test passes | YES | required |

## Key References

- Tabuada, P. (2007). Event-triggered real-time scheduling of stabilizing control tasks. *IEEE TAC*, 52(9), 1680–1685.
- Heemels, W.P.M.H., Johansson, K.H., & Tabuada, P. (2012). An introduction to event-triggered and self-triggered control. *IEEE CDC*.
- Girard, A. (2015). Dynamic triggering mechanisms for event-triggered control. *IEEE TAC*, 60(7), 1992–1997.
- Mazo, M. Jr., Anta, A., & Tabuada, P. (2010). On self-triggered control for linear systems. *IEEE TAC*, 55(8).
- Heemels, W.P.M.H., Donkers, M.C.F., & Teel, A.R. (2013). Periodic event-triggered control. *IEEE TAC*, 58(4), 847–861.
- Sontag, E.D. (2008). Input to state stability. In *Encyclopedia of Complexity and Systems Science*. Springer.
- Khalil, H.K. (2002). *Nonlinear Systems*, 3rd ed. Prentice Hall.
- Astrom, K.J. & Bernhardsson, B. (2002). Comparison of Riemann and Lebesgue sampling. *IEEE TAC*.
