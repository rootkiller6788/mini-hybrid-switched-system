# Course Tree — mini-impulsive-system

## Prerequisites (what this module depends on)
```
mini-general-system-theory
  └── State-space representation
mini-system-dynamics
  └── ODE theory, numerical integration
mini-stability-theory
  └── Lyapunov stability, dwell-time conditions
mini-linear-system-theory
  └── Matrix algebra, Riccati equations
mini-optimal-control
  └── LQR, discrete-time Riccati
mini-nonlinear-dynamics-chaos
  └── Chaotic systems (Lorenz, Rossler, Chua)
```

## Dependents (modules that depend on this one)
```
mini-impulsive-system  <-- this module
  ├── mini-event-triggered-control (sibling)
  ├── mini-hybrid-automata (sibling)
  ├── mini-switched-stability (sibling)
  ├── mini-networked-control
  │   └── Impulsive communication protocols
  ├── mini-multi-agent-coordination
  │   └── Impulsive consensus
  └── mini-quantum-control
      └── Quantum impulsive measurements
```

## Knowledge Flow
```
ODE Theory -> Impulsive DE -> Stability -> Dwell-Time -> Sync -> Applications
    |              |              |           |           |
    v              v              v           v           v
  Euler/RK4    TimeSeq        Lyapunov    Condition    Chaotic
                JumpMap       Bound       Adaptive     Masking
```

## Research Frontiers (L9)
- Meta-stability in impulsive networks
- Quantum Zeno effect as impulsive measurement
- ML for optimal impulse scheduling
- Network-coupled impulsive systems
