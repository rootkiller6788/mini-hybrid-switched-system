# Knowledge Graph - mini-supervisory-control

## L1: Definitions (Complete)
- DES automaton G = (Q, Sigma, delta, q0, Qm) as typedef struct des_automaton_t
- Controllable events Sigma_c, Uncontrollable Sigma_u (des_event_type_t)
- Supervisor S: L(G) -> Gamma (supervisor_t)
- Hybrid automaton H = (Q, X, f, Init, Inv, E, G, R) (hybrid_automaton_t)
- Guard, Invariant, Reset map (hyb_guard_t, hyb_invariant_t, hyb_transition_t)

## L2: Core Concepts (Complete)
- Controllability (Ramadge-Wonham): K prefix-closed, K subset L(G)
- Nonblocking: all reachable states can reach marked states
- Prefix-closure, natural projection, synchronous composition
- Modular/Decentralized/Hierarchical supervision
- DES abstraction of hybrid systems

## L3: Mathematical Structures (Complete)
- Finite automaton with adjacency matrix delta[Q][Sigma]
- Closed language L(G) and marked language Lm(G)
- Supervisor as state-feedback control pattern map
- Product automaton G1 || G2
- Linear continuous dynamics dx/dt = A*x + B*u per mode

## L4: Fundamental Laws (Complete)
- Controllability Theorem: K controllable <=> supervisor exists
- Nonblocking Controllability Theorem
- Existence of supremal controllable sublanguage supC(K,G)
- Observer property for partial observation

## L5: Algorithms/Methods (Complete)
- Controllability check (BFS on G x K product)
- supC(K,G) computation (fixpoint iteration)
- supCNB(K,G) controllable + nonblocking
- DFA minimization (Hopcroft refinement)
- Hybrid supervisor synthesis via DES abstraction

## L6: Canonical Problems (Complete)
- Transfer Line (manufacturing with buffer)
- Small Factory (two machines, one robot)
- Switched Power Converter (buck/boost safety)

## L7: Applications (Partial+)
- Power electronics (Tesla, DC-DC converters)
- Smart grid supervisory control (IEEE 1547)
- Building HVAC (Siemens, Johnson Controls)
- Aircraft power management (Boeing 787, F-35)

## L8: Advanced Topics (Partial+)
- Decentralized supervision with conjunctive fusion
- Partial observation and observer property
- Hybrid reachability analysis for safety

## L9: Research Frontiers (Partial)
- Learning-based supervisor synthesis (PAC learning)
- CEGIS for DES specifications
- Stochastic DES supervisors
