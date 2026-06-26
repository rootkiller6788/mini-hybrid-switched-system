# Course Tree: Dwell-Time Analysis

## Prerequisites
- **Linear algebra**: Eigenvalues, matrix exponential, Lyapunov equation
- **Linear system theory**: Stability, Hurwitz matrices, L2 gain
- **Nonlinear systems**: Lyapunov stability theory, invariance principles
- **Convex optimization**: LMI formulation, semidefinite programming

## Dependency Tree
```
Linear System Theory
    |
    v
Lyapunov Stability Theory (LDM)
    |
    +---> Common Lyapunov Function
    |
    +---> Multiple Lyapunov Functions
    |         |
    |         v
    |     Dwell-Time Analysis <--- This module
    |         |
    |         +---> Average Dwell Time
    |         |
    |         +---> Mode-Dependent Dwell Time
    |         |
    |         +---> State-Dependent Switching
    |
    v
Hybrid/Switched System Stability
    |
    +---> Event-Triggered Control
    |
    +---> Hybrid Automata
    |
    +---> Reset Control Systems
```

## Downstream Dependencies
- Event-triggered control (uses dwell time to guarantee Zeno-free behavior)
- Hybrid automata verification (dwell time as safety constraint)
- Model predictive control with switching (dwell time for recursive feasibility)
- Networked control (dwell time from maximum allowable transfer interval)

## L9 Research Frontiers
- Data-driven estimation of dwell times from trajectory data
- Safe RL with switching constraints (dwell time as safety shield)
- Quantum control: minimum dwell times for adiabatic transitions
