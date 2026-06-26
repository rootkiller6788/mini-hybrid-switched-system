# Course Alignment — Event-Triggered Control

## Nine-School Curriculum Mapping

| School | Course | Topics | Our Implementation |
|--------|--------|--------|--------------------|
| **MIT** | 6.241J Dynamic Systems | Lyapunov stability, state feedback | `etc_stability.c` — CLF computation, ISS verification |
| **MIT** | 16.323 Optimal Ctrl | Optimal event scheduling | Communication ratio analysis, threshold design |
| **Stanford** | AA203 Optimal Ctrl | Optimization-based control | `etc_design_threshold` — performance-driven σ selection |
| **Stanford** | EE363 Convex Opt | LMI-based stability conditions | L₂-gain analysis via Lyapunov inequalities |
| **Berkeley** | EE221A Linear Systems | Lyapunov equation, Hurwitz stability | `etc_solve_lyapunov`, `etc_compute_clf` |
| **Berkeley** | EE222 Nonlinear Systems | ISS, comparison functions | K/KL-class verification, ISS Lyapunov analysis |
| **CMU** | 24-677 Nonlinear Ctrl | Nonlinear event-triggered control | Dynamic triggering, mixed thresholds |
| **Princeton** | MAE 546 Optimal Ctrl | Optimal feedback design | State-feedback gain design for ETC |
| **Cambridge** | 4F3 Nonlinear Ctrl | Nonlinear stability analysis | Zeno analysis, ultimate boundedness |
| **Cambridge** | 4F2 Robust Ctrl | Robustness to disturbances | Robust STC, L₂-gain analysis |
| **Oxford** | B4 Predictive Ctrl | Model predictive control | Event-triggered MPC (documented in gaps) |
| **ETH** | 227-0216 Sys Identification | System modeling | State-space representation for ETC design |
| **ETH** | 227-0220 Model Reduction | Reduced-order ETC | (future work) |

## Topic-to-Course Matrix

| Topic | Primary Courses |
|-------|----------------|
| Event-triggered sampling | MIT 6.241J, Berkeley EE222 |
| Lyapunov-based triggering | MIT 6.241J, CMU 24-677 |
| Minimum inter-event time | MIT 6.241J, Cambridge 4F3 |
| Zeno behavior analysis | Berkeley EE222, Cambridge 4F3 |
| Periodic ETC | ETH 227-0216, Stanford AA203 |
| Self-triggered control | Princeton MAE 546, Oxford B4 |
| Dynamic triggering | CMU 24-677, Cambridge 4F3 |
| ISS characterization | Berkeley EE222, MIT 6.241J |
| L₂-gain analysis | Cambridge 4F2, Stanford EE363 |
| Networked control | ETH 227-0216, Oxford B4 |

## Reference Textbooks

| Textbook | Chapter | Topic |
|----------|---------|-------|
| Khalil (2002) — Nonlinear Systems | Ch. 4 | Lyapunov stability |
| Khalil (2002) | Ch. 9 | Input-to-state stability |
| Astrom & Wittenmark (2013) | Ch. 7 | Adaptive and event-based sampling |
| Tabuada (2007) — IEEE TAC | Full | Foundational ETC paper |
| Heemels et al. (2012) — IEEE TAC | Full | Introduction to ETC and STC |
| Girard (2015) — IEEE TAC | Full | Dynamic triggering mechanisms |
