# Knowledge Graph — mini-impulsive-system

## L1: Definitions (Complete)
- Impulsive Differential Equation (IDE): dx/dt=f(t,x), t!=tau_k; Delta x=I_k(x), t=tau_k
- Impulse time sequence {tau_k}
- Jump/reset map I_k(x)
- Dwell time: tau_{k+1} - tau_k
- Trigger types: time-driven, state-driven, event-driven
- Flow vector field f(t,x)
- C struct: ImpTimeSeq, ImpSystem, ImpSolution

## L2: Core Concepts (Complete)
- Stability of impulsive systems
- Comparison principle for impulsive ODEs
- Dwell-time condition for stability
- Average dwell-time (Hespanha & Morse 1999)
- Impulsive synchronization (master-slave)
- Impulsive consensus in multi-agent systems
- Event-triggered impulsive control
- C: imp_solver_simulate, imp_sync_*, imp_control_*

## L3: Mathematical Structures (Complete)
- State space R^n with matrix operations
- Linear flow: dx/dt = A*x
- Nonlinear flows: Lorenz, Chua, Duffing, Rossler, VdP, FHN, pendulum
- Jump maps: linear, affine, projection, hard-reset, nonlinear
- Lyapunov functions: quadratic V=x^TPx, norm-based V=||x||_p^q
- C struct: ImpFlowLinear, ImpFlowLorenz, ImpFlowChua, etc.

## L4: Fundamental Laws (Complete)
- Lyapunov stability theorem for impulsive systems
- Dwell-time bound: tau_D > -ln(rho)/alpha
- Comparison principle for scalar impulsive ODEs
- Leibnitz-Newton formula for impulsive integrals
- Floquet theory for periodic impulsive systems
- C: imp_lyap_check_stability, imp_lyap_dwell_time_bound
- Lean: dwell_time_stability, lyapunov_jump_condition

## L5: Algorithms/Methods (Complete)
- Euler method for IDE (1st order)
- RK4 method (4th order)
- Dormand-Prince RK45 (adaptive, 5th order)
- Lyapunov equation solver (iterative)
- Algebraic Riccati equation solver (Kleinman iteration)
- Binary search for impulse interval lookup
- C: imp_solver_euler_step, imp_solver_rk4_step, imp_solver_adaptive_step

## L6: Canonical Problems (Complete)
- Impulsive stabilization of unstable linear systems
- Chaotic synchronization (Lorenz master-slave)
- Multi-agent consensus with impulsive communication
- Event-triggered stabilization
- C: examples/example_impulsive_stabilization, sync, consensus

## L7: Applications (Partial+)
- Chaotic secure communication (chaotic masking)
- Multi-agent coordination (consensus protocol)
- Biological rhythm synchronization
- Spacecraft attitude control with impulsive thrusters
- C: imp_sync_chaotic_masking_*, consensus_step
- Keywords: NASA, Quadrotor, smart grid

## L8: Advanced Topics (Partial+)
- Meta-stability in impulsive networks
- Time-varying impulse sequences
- Adaptive impulsive control with online gain tuning
- Robustness to impulse timing jitter
- Stochastic impulsive systems
- C: imp_adaptive_control_*, imp_analysis_impulse_jitter_margin

## L9: Research Frontiers (Partial)
- Meta-complexity of impulse time optimization
- Quantum impulsive control
- Machine learning for optimal impulse scheduling
- Documented in course-tree.md and this file
