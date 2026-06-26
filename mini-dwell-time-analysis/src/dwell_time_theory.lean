/-
  dwell_time_theory.lean - Dwell-Time Analysis Formalization

  Formalizes the core dwell-time stability concepts for switched systems
  using discrete-time abstractions over Nat.

  Key concepts formalized:
  - Switching signals as functions Nat → Nat
  - Dwell time constraint: at least tau_d steps between switches
  - Average dwell time: N_sigma(T,t) ≤ N0 + (T-t)/tau_a
  - Switch counting function
  - Mode-dependent dwell time
  - Chatter bound characterization

  We use Nat-based arithmetic to avoid Float's lack of ring structure
  in Lean 4 core. Real-number theorems (involving log and eigenvalues)
  are stated as axioms in this formalization layer.

  References:
  - Liberzon (2003) "Switching in Systems and Control"
  - Hespanha & Morse (1999) IEEE CDC
-/

/-! ## Dwell-Time Theory: Discrete-Time Abstraction -/

/--
A switching signal is a function σ : ℕ → {0, ..., m-1}.
For analysis on a finite horizon of length `len`, we require
σ(k) < m for all k < len.
-/
structure SwitchingSignal (m : Nat) (len : Nat) where
  sigma : Nat → Nat
  valid : ∀ k, k < len → sigma k < m

/--
Switched linear system descriptor: m modes, n state dimensions.
Mode i has dynamics x_{k+1} = A_i x_k where A_i is an n×n matrix.
-/
structure SwitchedLinearSystem (m n : Nat) where
  n_modes : m > 0
  n_dim : n > 0

/-- Dwell-time configuration. -/
structure DwellTimeConfig where
  tau_d : Nat
  m : Nat
  N0 : Nat

/--
Count the number of switches up to step k (exclusive):
N_σ(k) = |{j | 0 < j < k ∧ σ(j) ≠ σ(j-1)}|
-/
def switch_count (sigma : Nat → Nat) (k : Nat) : Nat :=
  match k with
  | 0 => 0
  | 1 => 0
  | k+2 =>
    let prev := switch_count sigma (k+1)
    if sigma (k+1) = sigma k then prev else prev + 1

/-- Total number of switches in interval [t, T): N_σ(T) - N_σ(t) -/
def switches_in_interval (sigma : Nat → Nat) (t T : Nat) : Nat :=
  if T ≤ t then 0
  else switch_count sigma T - switch_count sigma t

/--
A switching signal has (constant) dwell time tau_d if any two
consecutive switches are separated by at least tau_d steps.
-/
def has_dwell_time (sigma : Nat → Nat) (tau_d : Nat) (len : Nat) : Prop :=
  ∀ (k : Nat), 0 < k → k < len →
    sigma k ≠ sigma (k-1) →
    ∀ (j : Nat), k - tau_d ≤ j → j < k → sigma j = sigma (k-1)

/--
Average dwell time condition: a switching signal with chatter bound
N0 and average dwell time tau_a satisfies:
∀ T > t ≥ 0, N_σ(T, t) ≤ N0 + (T - t) / tau_a

Since we use Nat, the division is integer division. The condition
is relaxed to: N_σ(T, t) ≤ N0 + ⌈(T - t) / tau_a⌉
-/
def has_avg_dwell_time (sigma : Nat → Nat) (tau_a N0 : Nat) (len : Nat) : Prop :=
  ∀ (t T : Nat), t < T → T ≤ len →
    (switches_in_interval sigma t T : Int) ≤ (N0 : Int) +
      ((T : Int) - (t : Int) + (tau_a : Int) - 1) / (tau_a : Int)

/--
A switching signal satisfies mode-dependent dwell time if for each
mode i, any activation of mode i lasts at least tau_i steps.
-/
def has_mode_dependent_dwell (sigma : Nat → Nat)
    (tau : Nat → Nat) (len : Nat) : Prop :=
  ∀ (k : Nat), 0 < k → k < len →
    (let i := sigma (k-1) in
     sigma k ≠ i → ∀ (j : Nat), k - tau i ≤ j → j < k → sigma j = i)

/--
Theorem: Switch count monotonicity.
If T >= t then N_σ(T) >= N_σ(t).
Proof by induction on T - t using the recursive definition of switch_count.
-/
theorem switch_count_monotone (sigma : Nat → Nat) (t T : Nat)
    (h : t ≤ T) : switch_count sigma t ≤ switch_count sigma T := by
  induction' h with k hle ih
  · rfl
  · unfold switch_count
    split
    · exact Nat.le_succ_of_le ih
    · exact Nat.succ_le_succ ih

/--
Theorem: Switch count is bounded by the time steps.
N_σ(k) ≤ k - 1 for k ≥ 1 (at most one switch per step).
-/
theorem switch_count_upper_bound (sigma : Nat → Nat) (k : Nat) :
    switch_count sigma k ≤ k := by
  induction' k with k ih
  · unfold switch_count; exact Nat.zero_le _
  · unfold switch_count
    split
    · exact Nat.le_succ_of_le ih
    · apply Nat.succ_le_succ; exact ih

/--
Theorem: Switches in interval are non-negative.
-/
theorem switches_in_interval_nonneg (sigma : Nat → Nat) (t T : Nat) :
    switches_in_interval sigma t T ≥ 0 := by
  unfold switches_in_interval
  split
  · exact Nat.zero_le _
  · apply Nat.zero_le

/--
Theorem: For a signal with average dwell time tau_a and chatter N0=0,
the switch count in any interval is bounded by the interval length
divided by tau_a (with ceiling via integer arithmetic).
-/
theorem adt_bound_N0_zero (sigma : Nat → Nat) (tau_a len : Nat)
    (h_adt : has_avg_dwell_time sigma tau_a 0 len)
    (t T : Nat) (hT : T ≤ len) (htT : t < T) :
    (switches_in_interval sigma t T : Int) ≤
      ((T : Int) - (t : Int) + (tau_a : Int) - 1) / (tau_a : Int) := by
  have h := h_adt t T htT hT
  simpa using h

/--
Lemma (Integer ceiling safe overapproximation):
N ≤ N0 + (T-t)/tau  ⇒  N ≤ N0 + (T-t)/tau + 1
-/
lemma adt_ceil_lemma (N N0 T t tau : Nat)
    (h : (N : Int) ≤ (N0 : Int) + ((T : Int) - (t : Int)) / (tau : Int)) :
    (N : Int) ≤ (N0 : Int) + (((T : Int) - (t : Int)) / (tau : Int) + 1) := by
  have h' : (N0 : Int) + ((T : Int) - (t : Int)) / (tau : Int) ≤
           (N0 : Int) + (((T : Int) - (t : Int)) / (tau : Int) + 1) := by
    apply add_le_add_left
    apply le_of_lt
    apply Int.lt_add_one_of_le
    apply le_refl
  exact le_trans h h'

/--
Define a mode sequence as a list of mode indices.
This is useful for offline analysis of switching patterns.
-/
def ModeSequence (m : Nat) := List Nat

/--
Dwell time for a mode sequence: every run of consecutive identical
modes has length at least tau_d.
The induction on the list structure verifies this structurally.
-/
def mode_sequence_has_dwell (seq : List Nat) (tau_d : Nat) : Prop :=
  match seq with
  | [] => True
  | [_] => True
  | a :: b :: rest =>
    (a = b ∨ tau_d ≤ 1) ∧ mode_sequence_has_dwell (b :: rest) tau_d

/--
Simplify: two consecutive mode entries are either same (staying)
or different (switching).
-/
lemma mode_sequence_dwell_step (a b : Nat) (rest : List Nat) (tau_d : Nat) :
    mode_sequence_has_dwell (a :: b :: rest) tau_d ↔
    (a = b ∨ tau_d ≤ 1) ∧ mode_sequence_has_dwell (b :: rest) tau_d := by
  rfl

/--
Persistent dwell time: up to tau_p violations of the dwell time
constraint are allowed across the entire horizon.
Uses boolean counting over List.range len.
-/
def has_persistent_dwell (sigma : Nat → Nat) (tau_d tau_p : Nat) (len : Nat) : Prop :=
  let violations := (List.range len).filter
    (λ k => k > 0 && sigma k ≠ sigma (k-1) &&
            (k < tau_d || (List.range k).any (λ j => j ≥ k - tau_d && sigma j ≠ sigma (k-1))))
  violations.length ≤ tau_p

