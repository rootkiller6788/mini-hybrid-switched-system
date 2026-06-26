# Course Tree — Event-Triggered Control

## Prerequisites

```
Linear System Theory (A,B,C,D, state-space)
    │
    ├── Lyapunov Stability Theory (V(x), V̇ < 0)
    │       │
    │       └── Input-to-State Stability (ISS)
    │               │
    ├── State Feedback Control (u = Kx, pole placement)
    │       │
    │       └── LQR Design (optimal K)
    │               │
    └── Sampled-Data Control (periodic sampling, ZOH)
            │
            └── ★ EVENT-TRIGGERED CONTROL ★
                    │
                    ├── Static ETC (Tabuada 2007)
                    │   ├── Trigger: |e| > σ|x|
                    │   ├── MIET: τ_min = σ/(||Acl|| + σ||BK||)
                    │   └── Stability: σ < λ_min(Q)/2||PBK||
                    │
                    ├── Dynamic ETC (Girard 2015)
                    │   ├── Internal variable η ≥ 0
                    │   ├── η̇ = −βη + σ|x|² − |e|²
                    │   └── Fewer events than static
                    │
                    ├── Periodic ETC (Heemels 2013)
                    │   ├── Check Γ at t = k·h
                    │   ├── h_max from stability condition
                    │   └── Zeno-free by construction
                    │
                    ├── Self-Triggered Control (Mazo 2010)
                    │   ├── τ_k = f(x(t_k)) computed at event
                    │   ├── No continuous monitoring
                    │   └── Conservative vs ETC
                    │
                    └── Advanced Topics
                        ├── Output-feedback ETC
                        ├── Decentralized ETC
                        ├── Event-triggered MPC
                        ├── Robust ETC (disturbances)
                        ├── Stochastic ETC
                        └── L9: Learning-based ETC
```

## Dependency Graph (within this module)

```
etc_core.h ──────┬── etc_trigger.h ──── etc_stability.h
                 │
                 ├── etc_dynamics.h ─── etc_periodic.h
                 │                    └── etc_self.h
                 │
                 └── (all headers depend on etc_core.h types)

src/etc_core.c           ─── Core data structures, vector/matrix math
src/etc_trigger.c        ─── Trigger function evaluations
src/etc_stability.c      ─── Lyapunov analysis, ISS, Zeno detection
src/etc_dynamics.c       ─── Simulation engine, flow maps
src/etc_periodic.c       ─── PETC implementation
src/etc_self.c           ─── STC implementation
src/event_triggered.lean ─── Formal verification (independent)
```

## Knowledge Dependencies (Across mini-complex-control-theory)

| This module depends on | From module |
|------------------------|-------------|
| State-space representation (A,B,C,D) | 3. mini-linear-system-theory |
| Lyapunov stability theory | 6. mini-stability-theory |
| State feedback control design | 3. mini-linear-system-theory |
| Nonlinear system concepts (ISS) | 15. mini-nonlinear-system-theory |
| Sampled-data / hybrid systems | 19. mini-hybrid-switched-system |
