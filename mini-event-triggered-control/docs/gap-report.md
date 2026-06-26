# Gap Report — Event-Triggered Control

## Current Gaps

### L4: Missing Formal Proofs (Low Priority)
- **Tabuada Theorem formal proof**: The MIET lower bound is implemented computationally and stated in Lean as an axiom. A full constructive proof would require real analysis in Lean.
- **Zeno-freeness formal proof**: Currently relies on computational checks and axiomatized Lean theorems.

### L7: Missing Applications (Medium Priority)
- **Smart grid load control with ETC**: IEEE 1547-compliant distributed energy resource coordination
- **Automotive CAN bus**: Event-triggered message scheduling for in-vehicle networks (Toyota, ISO 11898)
- **Quadcopter formation**: Multi-agent event-triggered consensus for drone swarms

### L8: Missing Advanced Topics (Medium Priority)
- **Decentralized event-triggered control**: Node-local triggering without global state knowledge
- **Event-triggered output feedback**: Observer-based ETC when full state is not measurable
- **Event-triggered MPC**: Combine model predictive control with event-triggered optimization

### L9: Research Frontiers (Low Priority)
- **Learning-based threshold adaptation**: Use reinforcement learning to select σ online
- **Security-aware ETC**: Detect and mitigate Denial-of-Service attacks on event channels
- **Quantum event-triggered measurement**: Exploit quantum Zeno effect for measurement scheduling

## Prioritized Action Items

| Priority | Gap | Effort | Impact |
|----------|-----|--------|--------|
| 1 | Output-feedback ETC (L8) | Medium | High — enables sensor-less operation |
| 2 | Smart grid application (L7) | Medium | High — real-world relevance |
| 3 | Decentralized ETC (L8) | Large | High — multi-agent systems |
| 4 | Tabuada formal proof (L4) | Large | Low — already computationally verified |
| 5 | Learning-based ETC (L9) | Large | Low — research frontier |
