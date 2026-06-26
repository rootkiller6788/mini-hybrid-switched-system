# Gap Report - mini-supervisory-control

## Priority 1: Essential
- None. All L1-L6 requirements met.

## Priority 2: Valuable
1. Full controllability proof in Lean (non-trivial)
2. Observer-based synthesis for partial observation
3. Timed DES extension

## Priority 3: Future
1. PAC learning of specifications (L9)
2. Stochastic DES with MDPs
3. Runtime verification integration

## Fixed Since Last Audit
- Controllability check: event matching by label (not index)
- Trim function: correct state 0 reuse
- Livelock detection: excludes trivial SCCs
- Supervisor observe(): enables check before state update
