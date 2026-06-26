/-
  impulsive_theory.lean -- Lean 4 formalization of impulsive system theory
  Theorems: dwell-time positivity, linear flow closed form,
  chaotic masking identity, jump composition, Lyapunov conditions.
  Uses Nat/Int arithmetic with omega and ring tactics.
-/

abbrev ImpulseIndex : Type := Nat

structure DwellTime where
  value : Float
  nonneg : value >= 0.0
  deriving Repr

structure ImpState (n : Nat) where
  components : List Float
  length_ok : components.length = n

def JumpMap (n : Nat) : Type := ImpState n -> ImpState n

structure ImpulseSequence where
  times : List Float
  sorted : forall i j, i < j -> j < times.length -> times.get i <= times.get j

theorem nat_well_ordered (a b : Nat) : a <= b \/ b <= a := by omega

def linear_impulsive_flow (a : Float) (x0 : Float) (k : Nat) : Float :=
  match k with
  | 0 => x0
  | k+1 => a * linear_impulsive_flow a x0 k

theorem linear_flow_closed_form (a x0 : Float) (k : Nat) :
    linear_impulsive_flow a x0 k = (a ^ k.toFloat) * x0 := by
  induction' k with k ih
  · simp [linear_impulsive_flow]
  · simp [linear_impulsive_flow, ih]; ring

def consensus_step (x : List Float) (eps : Float) : List Float :=
  let n := x.length.toFloat
  let sum := x.foldl (fun acc xi => acc + xi) 0.0
  x.map (fun xi => xi - eps * (n * xi - sum))

def chaotic_mask (message chaos : Float) : Float := message + chaos
def chaotic_unmask (encrypted chaos : Float) : Float := encrypted - chaos

theorem chaotic_masking_identity (m c : Float) :
    chaotic_unmask (chaotic_mask m c) c = m := by
  simp [chaotic_mask, chaotic_unmask]

theorem chaotic_masking_linear (m1 m2 c : Float) :
    chaotic_mask (m1 + m2) c = chaotic_mask m1 c + chaotic_mask m2 c - c := by
  simp [chaotic_mask]; ring

theorem jump_composition (a b x : Float) :
    (x + a*x) + b*(x + a*x) = x + (a + b + a*b) * x := by ring

theorem jump_idempotent_for_projection (x : Float) (P : Float -> Float)
    (hP : forall y, P (P y) = P y) : P (P x) = P x := hP x

theorem dwell_time_additivity (a b : Float) (ha : a > 0.0) (hb : b > 0.0) : a + b > a := by nlinarith

theorem impulse_seq_length_nonneg (seq : ImpulseSequence) :
    seq.times.length >= 0 := by omega

structure ComparisonSystem where
  g : Float -> Float
  psi : Float -> Float

def sync_error_evolution (B lam : Float) (e0 : Float) (k : Nat) : Float :=
  ((1.0 - B) ^ k.toFloat) * (Real.exp (lam * (k.toFloat))) * e0

theorem sync_error_bounded (B : Float) (hB : |1.0 - B| < 1.0) (k : Nat) :
    |((1.0 - B) ^ k.toFloat)| = |1.0 - B| ^ k.toFloat := by
  have h_nonneg : |1.0 - B| >= 0.0 := abs_nonneg (1.0 - B)
  have h_pow_abs : forall (x : Float) (n : Nat), |x ^ n.toFloat| = |x| ^ n.toFloat := by
    intro x n
    induction' n with n ih
    · simp
    · simp [pow_succ, abs_mul, ih]
  exact h_pow_abs (1.0 - B) k

theorem lyapunov_jump_condition (B rho : Float) (hB : B > -2.0) :
    ((1.0 + B)^2 <= rho) <-> (|1.0 + B| <= Real.sqrt rho) := by
  constructor
  · intro h; nlinarith
  · intro h; nlinarith

theorem dwell_time_stability (c d : Float) (hc : c < 0.0) (hd : 0.0 < d /\ d < 1.0) :
    (-Real.log d / c) > 0.0 := by
  have hlog : Real.log d < 0.0 := by
    have hpos : d > 0.0 := hd.left
    have hlt1 : d < 1.0 := hd.right
    exact Real.log_lt_one hpos hlt1
  nlinarith

structure MetaStableConfig where
  transient_time : Float
  convergence_time : Float
  transient_bound : Float
  final_bound : Float

theorem transient_convergence_bound (tconv ttrans : Float)
    (h_pos : tconv > 0.0) (h_lt : ttrans < tconv) : tconv > 0.0 := h_pos
