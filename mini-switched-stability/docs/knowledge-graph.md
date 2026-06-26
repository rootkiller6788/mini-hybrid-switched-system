# Knowledge Graph — Switched Stability

## L1: Definitions (Complete)
- SwitchedSystem: dx/dt = A_{sigma(t)} x(t)
- SwitchingSignal: sigma: [0,inf) -> P
- Common Lyapunov Function: V(x) = x^T P x, same P for all modes
- Multiple Lyapunov Functions: V_i(x) = x^T P_i x per mode
- Dwell time: tau_d = min interval between switches
- Average dwell time: tau_a = average interval
- Lie algebra: commutator [A_i, A_j] = A_i A_j - A_j A_i

## L2: Core Concepts (Complete)
- Arbitrary switching stability (CLF)
- Constrained switching stability (MLF)
- Slow switching (dwell time)
- Fast switching instability
- Switching attenuation mu

## L3: Mathematical Structures (Complete)
- Kronecker product for Lyapunov equation
- LMI feasibility problem
- QR eigenvalue algorithm
- Gershgorin circle theorem for eigenvalue bounds
- Sylvester criterion for positive definiteness

## L4: Fundamental Theorems (Complete)
- CLF Theorem (Liberzon 2003, Theorem 2.1)
- MLF Theorem (Branicky 1998, Theorem 1)
- Dwell-Time Theorem (Morse 1996)
- Average Dwell-Time Theorem (Hespanha & Morse 1999)
- Lie-Algebraic Stability Condition (Liberzon, Hespanha & Morse 1999)

## L5: Algorithms (Complete)
- QR eigenvalue decomposition with Wilkinson shift
- Kronecker product Lyapunov equation solver
- CLF gradient descent on LMI cone
- MLF per-mode Lyapunov equation
- Dwell-time computation and certification

## L6: Canonical Problems (Complete)
- DC-DC boost converter switching stability
- Thermostat bang-bang control
- Vehicle spacing with mode switching
- Networked control with packet dropouts

## L7: Applications (Complete)
- Power electronics (DC-DC converters)
- HVAC/building control (thermostat)
- Automotive (adaptive cruise control)
- Networked control systems (NCS/IoT)

## L8: Advanced Topics (Partial)
- Lie-algebraic stability condition (implemented)
- Dwell-time certification (implemented)

## L9: Research Frontiers (Partial)
- Data-driven stability of switched systems (documented)
- Switched nonlinear systems (documented)
- Safe switching in learning-based control (documented)