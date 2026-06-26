# Gap Report — mini-piecewise-affine-system

## Missing Items

### L7 Applications (2 of 2+ met)
- [x] Actuator saturation control
- [x] DC-DC power converter
- [ ] HVAC on/off control (documented)
- [ ] Robot contact dynamics (documented)

### L8 Advanced Topics (5 of 1+ met)
- [x] PWA-MPC (enumeration-based)
- [x] Explicit PWA control law
- [x] Forward reachability (support function)
- [x] Backward reachability (safety)
- [x] Mode-dependent LQR
- [ ] Binary search tree acceleration for point location
- [ ] SDP solver integration for LMI-based stability

### L9 Research Frontiers (documented)
- [ ] Full SDP-based PWQ Lyapunov computation (requires CSDP/SDPA)
- [ ] Neural PWA model equivalence (ReLU → PWA)
- [ ] Formal verification via reachability (requires dReal/SpaceEx)
- [ ] Stochastic PWA systems
- [ ] Online PWA identification
- [ ] PWA model reduction

## Priority Queue
1. HIGH: SDP solver interface for rigorous PWQ Lyapunov synthesis
2. MEDIUM: More L7 application examples (HVAC, robotics)
3. LOW: Neural network PWA conversion utilities
