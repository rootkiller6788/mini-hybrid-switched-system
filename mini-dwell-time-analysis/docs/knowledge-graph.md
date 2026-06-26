# Knowledge Graph: Dwell-Time Analysis

## L1 Definitions
- Dwell time τ_d: minimum time between consecutive switches
- Average dwell time τ_a: asymptotic switching frequency bound
- Mode-dependent dwell time τ_{d,i}: per-mode minimum dwell
- Persistent dwell time: τ_d with bounded violations
- Switched system: ẋ = f_σ(t)(x)
- Switching signal σ(t): [0,∞) → {1,...,m}
- Multiple Lyapunov functions: {V_i(x)} for each mode
- Coupling constant μ: V_i(x) ≤ μ V_j(x)
- Chatter bound N0: ADT slack variable
- Hurwitz/Schur stability of individual modes

## L2 Core Concepts
- Stability under arbitrary switching
- Slow switching (large τ_d → stability)
- Fast switching & chattering
- Common Lyapunov function vs multiple Lyapunov functions
- Hysteresis switching to prevent chattering
- State-dependent vs time-dependent switching
- Exponential convergence rate α(τ_d)
- L2 gain of switched systems

## L3 Mathematical Structures
- State-space triple (A_i, B_i, C_i) per mode
- Switching signal as piecewise-constant function
- Dwell-time inequality: N_σ(T,t) ≤ N0 + (T-t)/τ_a
- Lyapunov-Metzler matrix Π = (π_ij)
- Lyapunov equation: A_i^T P_i + P_i A_i = -Q_i
- Matrix exponential for transition: Φ(t,t0) = e^{A_i(t-t0)}
- LMI formulation: P_i > 0, A_i^T P_i + P_i A_i + 2λ_i P_i < 0

## L4 Fundamental Theorems
- Dwell-Time Stability Theorem (Liberzon 2003, Thm 3.1): τ_d > ln(μ)/(2λ) ⇒ GUES
- Average Dwell-Time Theorem (Hespanha & Morse 1999): τ_a > ln(μ)/λ ⇒ GUES
- Multiple Lyapunov Function Theorem (Branicky 1998)
- Common Lyapunov Function Theorem: ∃P: A_i^T P + P A_i < 0 ⇒ GAS under arbitrary switching
- Lyapunov-Metzler Theorem (Geromel & Colaneri 2006)
- Mode-Dependent ADT Theorem (Zhao & Hill 2008)

## L5 Algorithms/Methods
- QR algorithm for eigenvalue computation
- Bartels-Stewart Lyapunov equation solver
- Scaling+squaring matrix exponential
- Bisection search for minimum dwell time
- Gradient descent for common Lyapunov function
- Primal-dual interior-point LMI solver
- MLF construction via Lyapunov equations
- RK4 simulation with switching events
- Hysteresis band design algorithm

## L6 Canonical Problems
- Two-mode switched linear system stability
- Boost converter (DC-DC) switching analysis
- Thermostat on-off control with dwell time
- LQR with switching controllers
- Networked control with packet dropout (modeled as switching)
- Chattering detection and prevention

## L7 Applications
- Power electronics: DC-DC converter switching control
- Automotive: engine idle speed control with mode switching
- Aerospace: flight control mode transitions (NASA/GPS context)
- Process control: multi-mode chemical reactor switching
- Smart grid: power converter switching with dwell time
- Robotics: gait switching in legged robots

## L8 Advanced Topics
- Time-varying dwell time for non-autonomous switching
- Stochastic switching with Markov chain σ(t)
- Lyapunov-based model predictive control with switching
- Balanced realization for switched model reduction
- State-dependent dwell time adaptation
- Fuzzy-logic switching with dwell time guarantees

## L9 Research Frontiers
- Data-driven dwell time estimation from trajectories
- Safe reinforcement learning with dwell-time constraints
- Quantum control switching dwell times
- Information-theoretic bounds on switching frequency
- Meta-complexity of switching signal design
