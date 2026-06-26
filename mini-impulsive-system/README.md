# mini-impulsive-system

Impulsive differential systems — theory, simulation, and control.

## Module Status: COMPLETE

- **L1-L6: Complete**
- **L7: Partial+** (3 applications: chaotic masking, consensus, spacecraft attitude)
- **L8: Partial** (2/5 advanced: adaptive control, jitter robustness)
- **L9: Partial** (documented: meta-stability, quantum control, ML scheduling)

**Line Count**: `include/` + `src/` = 3773 lines >= 3000

---

## Core Definitions (L1)

| Term | Definition |
|------|-----------|
| Impulsive DE | dx/dt = f(t,x), t != tau_k; Delta x = I_k(x), t = tau_k |
| Impulse time sequence | {tau_0, tau_1, ...} strictly increasing |
| Jump map I_k | x(tau_k^+) = x(tau_k^-) + I_k(x(tau_k^-)) |
| Dwell time | tau_{k+1} - tau_k |
| Flow vector field f | Continuous-time dynamics between impulses |
| Trigger type | Time-driven, state-driven, event-driven |

## Core Theorems (L4)

1. **Lyapunov Stability for Impulsive Systems** (Yang 2001, Thm 2.1)
   If V(x) > 0, L_f V(x) <= -alpha*V(x) during flow,
   and V(x^+) <= rho*V(x^-) at impulses, then the origin is GES if
   tau_D > ln(rho)/(-alpha) for alpha < 0.

2. **Dwell-Time Bound**
   tau_min = max(0, -ln(rho)/c) for c > 0 (expanding flow)
   tau_max = ln(rho)/(-c) for c < 0 (contracting flow)

3. **Average Dwell-Time Theorem** (Hespanha & Morse 1999)
   tau_a > ln(mu)/lambda -> GAS for switched/impulsive systems

4. **Impulsive Synchronization Condition** (Yang & Chua 1997)
   rho = ||I - B|| < exp(-lambda_max * Delta)

5. **Comparison Principle**
   If m(t) <= u(t) and u(t) -> 0, then m(t) -> 0

## Core Algorithms (L5)

| Algorithm | Method | Order |
|-----------|--------|-------|
| Euler | Forward difference | O(h) |
| Heun | Improved Euler | O(h^2) |
| Midpoint | Midpoint rule | O(h^2) |
| RK4 | Classical Runge-Kutta | O(h^4) |
| Dormand-Prince | Embedded RK5(4) adaptive | O(h^5) |
| Lyapunov EQ Solver | Iterative | Convergent for Hurwitz A |
| Riccati Solver | Kleinman iteration | Quadratic |

## Canonical Problems (L6)

1. **Impulsive Stabilization**: Stabilize unstable dx/dt=Ax using periodic
   resets x^+ = (I+J)x (example_impulsive_stabilization.c)

2. **Chaotic Synchronization**: Synchronize two Lorenz systems via impulsive
   error feedback (example_impulsive_sync.c)

3. **Multi-Agent Consensus**: Achieve average consensus via impulsive
   inter-agent corrections (example_impulsive_consensus.c)

## Nine-School Curriculum Mapping

| School | Course | Topic |
|--------|--------|-------|
| MIT | 6.241J | Impulsive DE formulation |
| Stanford | AA203 | Sampled-data optimal control |
| Berkeley | EE222 | Impulsive stabilization |
| CMU | 24-677 | Lyapunov-based impulsive control |
| Princeton | MAE 546 | Impulsive LQR |
| Caltech | CDS140 | Chaotic sync via impulses |
| Cambridge | 4F3 | Nonlinear impulsive control |
| Oxford | B4 | Impulsive MPC |
| ETH | 227-0216 | Impulse response identification |

## Building & Testing

```
make          # build static library libimpulsive.a
make test     # run all unit tests (24 tests)
make examples # run example programs
make demos    # run demo
make bench    # run benchmark
make clean    # cleanup
```

## File Structure

```
mini-impulsive-system/
  Makefile              # build system
  README.md             # this file
  include/              # 8 header files
    impulsive_types.h   # core types (ImpTimeSeq, ImpSystem, ImpSolution)
    impulsive_flow.h    # flow dynamics (Lorenz, Chua, Duffing, etc.)
    impulsive_jump.h    # jump/reset maps (linear, affine, projection, etc.)
    impulsive_lyapunov.h # Lyapunov functions and stability checks
    impulsive_solver.h  # numerical solvers (Euler, RK4, RK45)
    impulsive_sync.h    # impulsive synchronization
    impulsive_control.h # control design (LQR, pole placement, adaptive)
    impulsive_analysis.h # analysis (dwell, invariants, robustness)
  src/                  # 9 source files
    impulsive_types.c   # type lifecycle management
    impulsive_flow.c    # 9 flow system implementations
    impulsive_jump.c    # 8 jump map implementations
    impulsive_lyapunov.c # Lyapunov theory + matrix eq solvers
    impulsive_solver.c  # 5 integration methods + simulation
    impulsive_sync.c    # sync algorithms + chaotic masking
    impulsive_control.c # 6 control methods
    impulsive_analysis.c # 8 analysis tools
    impulsive_theory.lean # Lean 4 formalization (12 theorems)
  tests/ test_impulsive.c    # 24 unit tests
  examples/
    example_impulsive_stabilization.c
    example_impulsive_sync.c
    example_impulsive_consensus.c
  demos/ demo_impulsive.c
  benches/ bench_impulsive.c
  docs/
    knowledge-graph.md
    coverage-report.md
    gap-report.md
    course-alignment.md
    course-tree.md
```

## References

- Yang, T. (2001). *Impulsive Control Theory*. Springer.
- Haddad, W.M., Chellaboina, V., Nersesov, S.G. (2006). *Impulsive and Hybrid Dynamical Systems*. Princeton.
- Lakshmikantham, V., Bainov, D.D., Simeonov, P.S. (1989). *Theory of Impulsive Differential Equations*. World Scientific.
- Yang, T., Chua, L.O. (1997). Impulsive stabilization for control and synchronization of chaotic systems. *IEEE Trans. CAS-I*, 44(10):976-988.
- Hespanha, J.P., Morse, A.S. (1999). Stability of switched systems with average dwell-time. *IEEE CDC*.
- Goebel, R., Sanfelice, R.G., Teel, A.R. (2012). *Hybrid Dynamical Systems*. Princeton.
