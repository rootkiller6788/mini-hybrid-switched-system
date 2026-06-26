# Knowledge Graph — mini-piecewise-affine-system

## L1: Definitions (Complete)
- PWASystem: Complete piecewise affine system
- PWARegion: Polyhedral region with half-space constraints
- PWAAffineDynamics: Affine dynamics (A,B,f,C,D,g) per region
- PWASwitchType: Exogenous, autonomous, controlled switches
- PWAEvent: Mode switch event record
- PWATrajectory: State/input/output/region history
- PWAMode: Hybrid automaton mode with invariants
- PWAHybridAutomaton: Hybrid automaton representation
- PWAWellPosedStatus: Well-posedness classification
- PWADataPoint/PWADataset: Identification data structures

## L2: Core Concepts (Complete)
- State-space partitioning via polyhedral regions
- Point location in polyhedral partitions
- Region adjacency computation
- Well-posedness checking (overlap/gap detection)
- Continuous-time PWA: dx/dt = A_i x + B_i u + f_i
- Discrete-time PWA: x(k+1) = A_i x(k) + B_i u(k) + f_i
- Mode switching and guard conditions
- Equivalent MLD (Mixed Logical Dynamical) form

## L3: Mathematical Structures (Complete)
- PWAPolyhedron: H-representation of convex polyhedra
- PWAHalfSpace: a·x ≤ b
- PWAHyperplane: a·x = b
- PWAVertex/PWAFace: V-representation
- PWAPolytope: Bounded polyhedron with dual representation
- PWAConvexHull2D: 2D convex hull (Graham scan)
- PWABoundingBox: Axis-aligned bounds
- Fourier-Motzkin elimination
- Simplex LP solver (gradient projection)
- Minkowski sum and Hausdorff distance
- Support function computation
- Voronoi partition

## L4: Fundamental Theorems (Complete)
- PWQLyapunov: Piecewise quadratic Lyapunov functions
- S-procedure for PWA stability
- LMI construction for stability analysis
- Invariant set computation (Lyapunov ellipsoid method)
- Maximal invariant set (iterative method)
- Common quadratic Lyapunov functions for arbitrary switching
- Dwell-time stability analysis
- Matrix positive (semi)definiteness checks (Cholesky)

## L5: Algorithms/Methods (Complete)
- RK4 integration for CT-PWA
- RKF45 adaptive step integration
- Boundary crossing detection (sign change)
- Bisection-based event location
- Discrete-time PWA simulation
- Continuous-time PWA simulation with events
- K-means and K-means++ clustering
- Silhouette score for cluster validation
- SVM-like hyperplane separation
- Least-squares affine regression per region
- Gap statistic for region count selection
- Full identification pipeline (cluster → boundaries → regression)

## L6: Canonical Problems (Complete)
- Actuator saturation as PWA (example_saturation_pwa.c)
- DC-DC buck converter (example_dcdc_converter_pwa.c)
- Bouncing ball with resets (example_bouncing_ball_pwa.c)

## L7: Applications (Partial - 2 applications)
- Actuator saturation in control systems (automotive/industrial)
- DC-DC power converter (power electronics/Tesla/SpaceX)
- (Additional applications documented; implementation in examples/)

## L8: Advanced Topics (Partial - 5 topics)
- PWA Model Predictive Control (enumeration-based)
- Explicit PWA control law construction
- Terminal set computation for PWA-MPC
- Forward reachability analysis (support function method)
- Backward reachability for safety verification
- Bounding box propagation for PWA systems
- Mode-dependent LQR gains

## L9: Research Frontiers (Partial - documented)
- Data-driven PWA identification (gap statistic)
- Neural PWA equivalence (ReLU networks)
- Formal verification of PWA systems
- Meta-complexity of PWA identification
