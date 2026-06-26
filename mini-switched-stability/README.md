# Mini-Switched-Stability

Switched system stability analysis — Common/Multiple Lyapunov Functions,
Dwell-Time, and Lie-Algebraic stability conditions for hybrid switched systems.

## Core Concept

A **switched system** has dynamics:
```
dx/dt = f_{sigma(t)}(x, t)
```
where sigma(t): [0, infinity) -> P = {1, ..., p} is the switching signal
selecting which subsystem governs the evolution.

**Key question**: Under what conditions on sigma(t) is the system stable?

## Key Theorems

### Common Lyapunov Function (Liberzon 2003, Theorem 2.1)
```
If there exists P = P^T > 0 such that:
  A_i^T P + P A_i < 0  for all i in P
then the system is GUES under arbitrary switching.
```

### Multiple Lyapunov Functions (Branicky 1998)
```
For each mode i: V_i(x) = x^T P_i x
Stability if: dV_i/dt < 0 when sigma = i
AND V_j(x(t_k)) <= V_i(x(t_k)) at switch times
```

### Dwell-Time Theorem (Morse 1996)
```
If all A_i are Hurwitz with margin lambda_0 > 0,
there exists tau_d = ln(mu) / lambda_0 such that
t_{k+1} - t_k >= tau_d => GUES
```

### Average Dwell-Time (Hespanha & Morse 1999)
```
If tau_a > ln(mu) / lambda_0, then:
N_sigma(T, t) <= N_0 + (T - t) / tau_a => GUES
```

### Lie-Algebraic Condition (Liberzon, Hespanha & Morse 1999)
```
If Lie{A_1, ..., A_p} is solvable => CLF exists => GUES
Special case: pairwise commuting A_i => trivial CLF
```

## Nine-Layer Knowledge Coverage

| Level | Name | Status | Lines/Count |
|-------|------|--------|-------------|
| **L1** | Definitions | Complete | 15+ typedefs, 4 enums |
| **L2** | Core Concepts | Complete | CLF, MLF, dwell-time, Lie algebra |
| **L3** | Math Structures | Complete | Vector, Matrix, Kronecker, LMI |
| **L4** | Fundamental Laws | Complete | CLF Theorem, MLF Theorem, Dwell-Time, Avg Dwell-Time |
| **L5** | Algorithms | Complete | QR eigenvalues, Lyapunov eq solver, CLF gradient descent, MLF computation |
| **L6** | Canonical Problems | Complete | DC-DC converter, thermostat, vehicle spacing |
| **L7** | Applications | Complete | Power electronics, HVAC, NCS with dropouts |
| **L8** | Advanced Topics | Partial | Lie-algebraic condition, dwell-time certification |
| **L9** | Research Frontiers | Partial | Documented only |

## Core Definitions (L1)

- **SwitchedSystem**: Complete switched system with multiple LTI subsystems
- **SwitchingSignal**: sigma(t) with 8 signal types (arbitrary, periodic, dwell-time, etc.)
- **CommonLyapunovFunction**: Single P = P^T > 0 for all subsystems
- **MultipleLyapunovFunctions**: Per-mode P_i with switching attenuation mu
- **DwellTimeAnalysis**: tau_d, tau_a, stability margin, certification
- **LieAlgebraCondition**: Solvability, nilpotency, commutativity checks

## Core Algorithms (L5)

| Algorithm | Complexity | Description |
|-----------|-----------|-------------|
| QR Eigenvalue | O(n^3) | Wilkinson shift QR with Householder reflections |
| Lyapunov Solver | O(n^6) | Kronecker product + Gaussian elimination |
| CLF Gradient Descent | O(max_iter * p * n^3) | LMI feasibility via projected gradient |
| MLF Computation | O(p * n^6) | Per-mode Lyapunov equation + mu parameter |
| Dwell-Time Analysis | O(p * n^3) | lambda_0, mu, tau_d, tau_a* computation |
| Lie-Algebraic Check | O(p^2 * n^3) | Commutator [A_i, A_j] and solvability |

## Course Alignment

| School | Course | Topic |
|--------|--------|-------|
| MIT | 6.241J Dynamic Systems | Lyapunov stability, state-space |
| Stanford | AA203 Optimal Control | Switched LQR, hybrid optimal control |
| Berkeley | EE222 Nonlinear Systems | Multiple Lyapunov, dwell-time |
| Cambridge | 4F2 Robust Control | Switching stability, dwell-time |
| ETH | 227-0216 Sys Identification | Switched system modeling |
| Caltech | CDS140 Nonlinear Dynamics | Hybrid/switched dynamics |

## Build & Test

```bash
make          # Build static library
make test     # Run 27 test assertions
make examples # Build 3 examples
make demo     # Run all examples
make clean    # Clean build artifacts
```

## Directory Structure

```
mini-switched-stability/
  Makefile                  # make test/demo/clean targets
  README.md                 # This file
  include/
    switched_types.h        # Core types, 25+ structs/enums (433 lines)
    switched_stability.h    # Stability analysis API (262 lines)
    switched_lyapunov.h     # Lyapunov computation API (182 lines)
    switched_dwell_time.h   # Dwell-time analysis API (181 lines)
    switched_applications.h # DC-DC, thermostat, vehicle, NCS (238 lines)
  src/
    switched_types.c        # Vector, matrix, Lie bracket, QR algorithm (856 lines)
    switched_stability.c    # CLF/MLF theorems, dwell-time, simulation (759 lines)
    switched_lyapunov.c     # Eigenvalue bounds, Lyapunov equation solvers (503 lines)
    switched_dwell_time.c   # Dwell-time computation, certification (245 lines)
    switched_applications.c # DC-DC converter, thermostat, VSC, NCS (350 lines)
    switched_stability.lean # Lean 4 formal verification: 7 theorems
  tests/
    test_switched_stability.c  # 27 test assertions
  examples/
    example1_dcdc_converter.c  # Boost converter stability
    example2_thermostat.c      # Bang-bang thermostat control
    example3_networked_control.c # NCS with packet dropouts
  docs/
    knowledge-graph.md      # Nine-layer knowledge graph
    coverage-report.md      # Coverage evaluation
    gap-report.md           # Knowledge gaps
    course-alignment.md     # Nine-school curriculum mapping
    course-tree.md          # Prerequisite dependency tree
```

## Key References

- Liberzon, D. (2003). *Switching in Systems and Control.* Birkhauser.
- Branicky, M.S. (1998). Multiple Lyapunov functions. *IEEE TAC*, 43(4), 475-482.
- Hespanha, J.P. & Morse, A.S. (1999). Stability of switched systems with average dwell-time. *IEEE CDC*.
- Morse, A.S. (1996). Supervisory control of families of linear set-point controllers. *IEEE TAC*, 41(10), 1413-1431.
- Liberzon, D., Hespanha, J.P. & Morse, A.S. (1999). Stability of switched systems: a Lie-algebraic condition. *Systems & Control Letters*, 37(3), 117-122.
- Sun, Z. & Ge, S.S. (2011). *Stability Theory of Switched Dynamical Systems.* Springer.
- Boyd, S. et al. (1994). *Linear Matrix Inequalities in System and Control Theory.* SIAM.
- Golub, G.H. & Van Loan, C.F. (2013). *Matrix Computations.* 4th ed. Johns Hopkins.

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Complete (3 applications: power electronics, HVAC, NCS)
- **L8**: Partial (lie-algebraic condition, dwell-time certification implemented)
- **L9**: Partial (documented, not implemented)