# Gap Report — mini-hybrid-automata

## Identified Gaps

### L8 Advanced Topics

| Gap | Priority | Rationale |
|-----|----------|-----------|
| Stochastic hybrid automata | Medium | Requires probability distribution libraries; conceptually mapped but not implemented |
| Full CEGAR with predicate refinement | Low | Basic CEGAR loop implemented; predicate discovery would require SMT/SAT solver integration |

### L9 Research Frontiers

| Gap | Priority | Rationale |
|-----|----------|-----------|
| Learning-based abstraction | Low | Active research area; foundational reachability primitives provided |
| Quantum hybrid control | Low | Out of scope for this module |

## Resolved Gaps

- **Barrier certificate synthesis**: Implemented as numerical validation of all four barrier conditions (I-IV)
- **Common Lyapunov function**: Implemented with Hurwitz matrix check for switched linear systems
- **Multiple Lyapunov functions**: Implemented with transition non-increase verification
- **Region equivalence for timed automata**: Implemented with full fractional ordering and successor computation
- **Zonotope reachability**: Implemented with affine transform, Minkowski sum, and support function methods

## Recommendations

1. **L8**: Integrate with an SMT solver (e.g., Z3) for full CEGAR predicate discovery
2. **L8**: Add stochastic simulation wrapper using rejection sampling
3. **L9**: Add connection to learning-based control via reward function over hybrid traces
