# mini-supervisory-control

Ramadge-Wonham Supervisory Control Theory for Hybrid Switched Systems

## Module Status: COMPLETE

- L1-L6: Complete
- L7: Partial+ (2+ applications)
- L8: Partial+ (decentralized, observer, minimization)
- L9: Partial (documented, not implemented)

## Core Definitions

| Definition | Type | Description |
|-----------|------|-------------|
| Automaton G | des_automaton_t | G = (Q, Sigma, delta, q0, Qm) |
| Controllable Event | DES_EVENT_CONTROLLABLE | Can be disabled by supervisor |
| Uncontrollable Event | DES_EVENT_UNCONTROLLABLE | Always enabled |
| Supervisor S | supervisor_t | S: L(G) -> Gamma, Gamma = {gamma | Sigma_u subset gamma} |
| Hybrid Automaton | hybrid_automaton_t | Modes + continuous dynamics + guards |

## Core Theorems

1. **Controllability Theorem** (Ramadge-Wonham 1987):
   K controllable wrt G AND K prefix-closed AND K subset L(G)
   ⇔ EXISTS supervisor S: L(S/G) = K

2. **Nonblocking Controllability**: supCNB(K,G) = maximal sublanguage
   of K that is both controllable and nonblocking

3. **Supremal Controllable Sublanguage**: supC(K,G) always exists
   and is computable via fixpoint iteration

## Core Algorithms

| Algorithm | File | Complexity |
|-----------|------|------------|
| Controllability Check | controllability.c | O(|Q|·|Sigma|) |
| supC(K,G) Computation | synthesis.c | O(|Q|^2·|Sigma|) |
| supCNB(K,G) | synthesis.c | O(|Q|^2·|Sigma|·iter) |
| DFA Minimization | synthesis.c | O(|Q|^2·|Sigma|) |
| DES Abstraction | hybrid_supervisor.c | O(|modes|·|trans|) |
| Hybrid Simulation | hybrid_supervisor.c | O(T/dt·|modes|) |

## Canonical Problems

1. Transfer Line (manufacturing buffer control)
2. Small Factory (two machines, one robot)
3. Switched Power Converter (buck/boost safety)

## Course Mapping

MIT 6.241J · Stanford AA274 · Berkeley EE222 · CMU 24-654/24-677 ·
Princeton MAE 546 · Caltech CDS110 · Cambridge 4F3 · Oxford C20 · ETH 227-0220

## Build & Test

    make test        # Run 20 tests (all pass)
    make examples    # Run 3 end-to-end examples
    make demo        # Interactive demo
    make bench       # Performance benchmark
    make lines       # Line count

## File Structure

    include/  5 header files (624 lines)
    src/      5 C source files + 1 Lean file (2570+287 lines)
    tests/    1 test file (20 test cases)
    examples/ 3 example programs
    docs/     5 knowledge documents
    demos/    1 interactive demo
    benches/  1 benchmark

include/ + src/ C code total: 3189 lines
