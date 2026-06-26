# Mini Hybrid Switched System

A collection of **from-scratch, zero-dependency C implementations** of hybrid and switched dynamical system theory. Hybrid systems combine continuous dynamics (differential equations) with discrete mode transitions (automata, guards, resets), modeling phenomena from bouncing balls and thermostats to advanced power electronics, networked control, and autonomous systems.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|-----------|--------|-------------|
| [mini-dwell-time-analysis](mini-dwell-time-analysis/) | Dwell-time stability, average dwell time, LMI-based computation, common/multiple Lyapunov functions, slow switching | MIT 6.241J, Stanford AA203 |
| [mini-event-triggered-control](mini-event-triggered-control/) | Event-triggered sampling, self-triggered control, Zeno detection, inter-event time analysis, relative/absolute thresholds | MIT 6.241J, UC Santa Barbara (Tabuada) |
| [mini-hybrid-automata](mini-hybrid-automata/) | Hybrid automaton modeling (Henzinger), forward/backward reachability, barrier certificates, hybrid simulation, Zeno analysis | MIT 6.841, Stanford CS359, Berkeley EECS 291E |
| [mini-impulsive-system](mini-impulsive-system/) | Impulsive dynamics, jump maps, B-equivalence method, impulsive stabilization, Lyapunov-based impulsive control | MIT 6.241J, Caltech CDS140 |
| [mini-piecewise-affine-system](mini-piecewise-affine-system/) | PWA system identification, polyhedral geometry, PWA simulation/reachability, explicit MPC | ETH 227-0216, MIT 6.832 |
| [mini-reset-control-system](mini-reset-control-system/) | Clegg integrator, reset band, first-order reset element (FORE), Bode integral constraint overcoming, hybrid loop analysis | Cambridge 4F2, MIT 6.241J |
| [mini-supervisory-control](mini-supervisory-control/) | DES automata, Ramadge-Wonham supervisory theory, supremal controllable sublanguage, nonblocking, modular supervision | MIT 6.241J, Cambridge 4F3 |
| [mini-switched-stability](mini-switched-stability/) | Common Lyapunov functions, multiple Lyapunov functions (Branicky), Lie algebra conditions, slow switching, Lyapunov-Metzler inequalities | Berkeley EECS 291E, MIT 6.241J |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes and knowledge layer documentation
- **Practical demos** — bouncing ball, thermostat, DC-DC converter, vehicle platoon, manufacturing cell, and more

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-dwell-time-analysis
make all    # build everything
make test   # run tests
```

Or build the top-level unified framework:

```bash
make          # build library, tests, and examples
make test     # run all top-level tests
make test-all # run top-level + all sub-module tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-hybrid-switched-system/
├── include/                          # Top-level unified HSS framework headers
├── src/                              # Top-level unified HSS framework source
├── demos/                            # Demonstration programs
├── examples/                         # Usage examples
├── tests/                            # Top-level test suites (48 tests)
├── benches/                          # Performance benchmarks
├── docs/                             # Course-alignment and theory documentation
├── mini-dwell-time-analysis/         # Dwell-time stability analysis & LMI computation
├── mini-event-triggered-control/     # Event-triggered & self-triggered control
├── mini-hybrid-automata/             # Hybrid automata modeling & reachability
├── mini-impulsive-system/            # Impulsive dynamics & stabilization
├── mini-piecewise-affine-system/     # Piecewise affine system identification
├── mini-reset-control-system/        # Reset control (Clegg integrator, FORE)
├── mini-supervisory-control/         # Supervisory control (Ramadge-Wonham)
├── mini-switched-stability/          # Switched stability (CLF, MLF, Lie algebra)
└── Makefile                          # Top-level build system
```

## License

MIT
