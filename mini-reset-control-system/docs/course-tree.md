# Course Dependency Tree - Reset Control Systems

## Prerequisites
```
Linear System Theory
  +-- Nonlinear Dynamics -> Hybrid Systems -> [RESET CONTROL]
  +-- Stability Theory (Lyapunov) -> [RESET CONTROL]
  +-- Optimal Control (Riccati) -> [RESET CONTROL]
```

## Internal Dependencies
```
reset_core.h/c
  |-- reset_element.h/c
  |-- reset_system.h/c
  |-- reset_simulation.h/c
  |-- reset_math.h/c
  |-- reset_lyapunov.h/c
```

## Postrequisites
```
[RESET CONTROL] -> Event-Triggered Control
                -> Switched Systems Stability
                -> Hybrid Optimal Control
                -> Learning-Based Control
```
