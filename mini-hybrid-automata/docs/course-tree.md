# Course Tree — mini-hybrid-automata

## Prerequisites

```
mini-hybrid-automata
├── mini-complexity-foundations (decidability concepts)
├── mini-circuit-complexity (guards as boolean functions)
├── Calculus (ODE theory, Lyapunov stability)
├── Linear Algebra (matrix exponential, polyhedra)
└── Graph Theory (discrete mode graph, adjacency)
```

## Dependency Graph Within Module

```
L1: Definitions (KP1-9)
 ├── L2: Execution Semantics (KP10-19)
 │    ├── L5: Simulation (KP36-40)
 │    ├── L3: Reachability (KP20-28)
 │    │    ├── L5: CEGAR (KP27)
 │    │    └── L3: Zonotopes/Polyhedra (KP23-25)
 │    └── L4: Safety (KP30-35)
 ├── L3: Subclasses (KP41-48)
 │    └── L4: Decidability Theorems
 └── L6: Canonical Problems (KP49-56)
      └── L7: Applications (KP54-56)
```

## Build Order

1. hybrid_automaton.h/c — Core data structures (L1)
2. hybrid_execution.h/c — Execution semantics (L2)
3. hybrid_reachability.h/c — Reachability sets + operations (L3, L5)
4. hybrid_simulation.h/c — ODE integration + event detection (L5)
5. hybrid_safety.h/c — Barrier certs, Lyapunov, invariants (L4)
6. hybrid_subclass.h/c — Timed, rectangular, LHA, PWA (L3-L4)
7. hybrid_examples.h/c — Canonical models (L6-L7)
