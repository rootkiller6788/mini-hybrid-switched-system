# Knowledge Graph - Reset Control Systems

## L1: Definitions (Complete)
- Reset Controller, Reset Condition, Reset Map, After-Reset State
- Reset Trigger Types: ZC, rising, falling, periodic, state-condition, time-dep
- Reset Map Types: zero, ratio, matrix, setpoint, function
- Clegg Integrator: integrator with full reset on input zero-crossing
- FORE: First-Order Reset Element with proportional reset
- SORE: Second-Order Reset Element with dual reset ratios
- Reset PID: PID with reset applied to integrator term
- Reset Lead-Lag: lead/lag compensator with reset capability
- ResetMemory: snapshot of before/after reset states

## L2: Core Concepts (Complete)
- Reset Times and Reset Interval tracking
- Reset Band with arming/disarming hysteresis
- Reset Ratio rho in [0,1): fraction of state surviving after reset
- Dwell Time: minimum time between resets for Zeno prevention
- Zeno Behavior: infinitely many resets in finite time
- Hybrid Time Domain (t,j): continuous time + jump count
- Describing Function: fundamental-component frequency-domain analysis
- Phase Advantage: Clegg integrator provides 38.15 deg lead over linear
- H_beta Condition: Lyapunov-based stability condition

## L3: Mathematical Structures (Complete)
- ResetLinearBase: (A,B,C,D) state-space with dimensions (n,m,p)
- ResetJumpMap: (Ar,Br) reset matrix and input vector
- ResetCondition: trigger type + threshold + dwell time
- ResetSystem: complete hybrid automaton bundle
- Flow Set C and Jump Set D partition of state space
- Matrix Algebra: multiply, invert, transpose, determinant, trace
- Transfer Function: G(s) = C(sI-A)^-1 B + D

## L4: Fundamental Theorems (Complete)
- Lyapunov Stability for Reset: V(x)>0, Vdot<0 on flow, DeltaV<=0 on jump
- H_beta Theorem: reset system stable if base closed-loop is Hurwitz [BB12, Thm 5.1]
- Passivity Theorem: reset preserves passivity if ||Ar|| <= 1
- Circle Criterion: frequency-domain sufficient stability condition
- Dwell-Time Stability: tau_d > ln(beta)/(2*alpha) guarantees stability
- Hurwitz Criterion: A is Hurwitz iff all eigenvalues have Re < 0

## L5: Algorithms/Methods (Complete)
- Zero-Crossing Detection with bisection refinement
- QR Algorithm: Francis double-shift for eigenvalues
- Hessenberg Reduction: Householder-based for QR
- Bartels-Stewart Algorithm: Lyapunov equation solver
- Power Iteration: dominant eigenvalue computation
- Forward Euler, Heun, RK4: numerical integration
- Hybrid Simulation: flow integration + ZC detection + jump application
- Frequency Grid Search: H-infinity norm computation
- Newton Iteration for ARE: algebraic Riccati equation

## L6: Canonical Problems (Complete)
- Clegg Integrator Step Response (ex_clegg_step.c)
- FORE Response at Multiple Reset Ratios (ex_fore_response.c)
- Reset PID Step Response with Component Decomposition (ex_reset_pid.c)
- Reset Feedback Loop: unity feedback with reset controller (ex_reset_feedback.c)

## L7: Applications (Partial+)
- Servomechanism Phase Compensation
- Industrial Process Control with Reset PID
- Vibration Control using Reset Lead-Lag
- Precision Motion Systems with SORE

## L8: Advanced Topics (Partial+)
- Lyapunov Function Search Algorithm (lyap_search)
- Passivity-Based Reset Analysis
- Time-Varying Reset Ratio
- Reset System Interconnection (series/parallel/feedback)

## L9: Research Frontiers (Partial)
- Reset-Based Event-Triggered Control
- Stochastic Reset Control
- Learning-Based Reset Optimization
- Multi-Agent Reset Coordination
