# Course Tree — mini-piecewise-affine-system

## Prerequisites
```
                            ┌─────────────────────┐
                            │ mini-piecewise-affine │
                            │      -system         │
                            └──────────┬──────────┘
                                       │
        ┌──────────────┬───────────────┼───────────────┬──────────────┐
        │              │               │               │              │
   ┌────▼────┐   ┌─────▼─────┐  ┌──────▼──────┐  ┌─────▼─────┐ ┌─────▼─────┐
   │ Linear  │   │ Polyhedral│  │  Lyapunov   │  │  Hybrid   │ │   ODE     │
   │ Algebra │   │ Geometry  │  │  Stability  │  │ Automata  │ │ Solvers   │
   └────┬────┘   └─────┬─────┘  └──────┬──────┘  └─────┬─────┘ └─────┬─────┘
        │              │               │               │              │
   ┌────▼────┐   ┌─────▼─────┐  ┌──────▼──────┐  ┌─────▼─────┐ ┌─────▼─────┐
   │ Matrix  │   │ Convex    │  │ Quadratic   │  │  State    │ │ Numerical │
   │ Ops     │   │ Polytopes │  │ Forms      │  │ Machines  │ │ Methods   │
   └─────────┘   └───────────┘  └─────────────┘  └───────────┘ └───────────┘
```

## Downstream Dependencies
- mini-hybrid-automata (equivalent representation)
- mini-switched-stability (PWQ Lyapunov functions)
- mini-dwell-time-analysis (mode-dependent dwell time)
- mini-reset-control-system (impulsive/PWA with resets)
- mini-event-triggered-control (boundary crossing detection)
- mini-supervisory-control (discrete supervisor synthesis)
- mini-impulsive-system (PWA with state jumps)

## Knowledge Dependencies
1. **Linear System Theory**: A,B,C,D matrices, state-space
2. **Convex Geometry**: Polyhedra, half-spaces, convex hull
3. **Lyapunov Theory**: Quadratic forms, LMIs, stability
4. **Hybrid Systems**: Modes, guards, transitions
5. **Numerical Analysis**: ODE solvers, root-finding
6. **Optimization**: LP, QP, SDP (for stability/MPC)
7. **Machine Learning**: Clustering, regression (for identification)
