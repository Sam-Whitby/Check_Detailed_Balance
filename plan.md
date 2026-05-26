# Implementation Plan: Continuous-Limit VMMC with Gaussian Displacement Proposals

*Revised: 26 May 2026*

---

## 1. Problems with the Current `vmmc_continuous.wl`

**Problem 1 — The LJ potential is sub-lattice.**
`sigLJ = 2^(-1/6)` places the LJ zero-crossing at `d² ≈ 0.79`, below the minimum particle separation on the lattice (`d² = 1`). Every particle pair sits in the attractive tail; there is no repulsive core and no resolved particle diameter. The model is a nearest-neighbour attractive lattice gas, not an LJ fluid.

**Problem 2 — K directions is an unmotivated free parameter.**
Directions are chosen from K uniformly-spaced angles, rounded to grid coordinates. K and nGrid are two independent quantities that both need to grow to reach the continuous limit, with no principled relationship between them. Rounding collapses many angles to the same grid vector at small step sizes.

**Goal:** Replace both with a single physical length scale `physLen` (particle diameter in lattice units) and a Gaussian displacement proposal. K ceases to exist as a user parameter.

---

## 2. Physical Design

### 2.1 The `physLen` parameter

`physLen` is the particle diameter in lattice units:

```
sigLJ   = physLen                         (* LJ zero-crossing = particle diameter *)
d_min²  = 2^(1/3) · physLen²             (* LJ energy minimum in squared lattice units *)
```

For `physLen = 5`: the LJ minimum is at `d ≈ 5.6` lattice units. At adjacent sites (`d² = 1`), the LJ energy is approximately `4 · epsLJ · (5^12 − 5^6) ~ 10^9 epsLJ` — strongly repulsive, correctly enforcing a soft hard-core at one particle diameter.

For the interaction cutoff: during the symbolic check, `$maxD2 = Infinity` always (see §3). For numerical simulation, the user sets `$maxD2 = Ceiling[2 · physLen²]` to capture the full LJ well.

### 2.2 Gaussian displacement proposal

**Physical motivation.** In overdamped Brownian dynamics (Ermak & McCammon 1978; Allen & Tildesley), the displacement per timestep is drawn from `N(0, σ)` where `σ = √(2Dδt)` and `D = kT/γ`. In reduced units where `epsLJ` sets the energy scale, `physLen` sets the length scale, and `τ_LJ = physLen · √(m/epsLJ)` sets the natural timescale:

```
sigStep = physLen · √(2 / (numBeta · epsLJ))
```

This is the natural starting point (one BD timestep at the natural LJ timescale, m = 1, γ = 1). The user overrides `sigStep` directly in Section 0 to tune acceptance rates; the BD formula is the principled initial value, not a constraint.

The cluster translation vector `(dx, dy)` is drawn as:

```mathematica
dx = Round[RandomVariate[NormalDistribution[0, sigStep]]];
dy = Round[RandomVariate[NormalDistribution[0, sigStep]]];
dir = {dx, dy}
```

`Round` is necessary and unavoidable. Mathematica's precision machinery (`SetPrecision`, `WorkingPrecision`) controls significant figures in floating-point arithmetic, not quantisation to a grid. There is no way to make a normally distributed real number land on lattice sites without explicit rounding or direct use of a discrete distribution.

**Symmetry.** `p(dx, dy) = p(−dx, −dy)` holds exactly because the Gaussian is even. This is the only requirement on the proposal for superdetailed balance — the Whitelam-Geissler link mechanism is unchanged.

**Zero displacement.** When `dir = {0, 0}`, the move is a no-op. This contributes to `T(i→i)` and is correctly handled by the checker.

### 2.3 The role of K

K is no longer a user parameter. On a periodic lattice of side `nGrid`, there are at most `(2·floor(nGrid/2) + 1)²` distinct integer `(dx, dy)` pairs. The Gaussian proposal visits all of them with non-zero probability. The effective number of directions is determined by the lattice geometry and the checker's NMax setting (§4), not by any user-supplied constant.

---

## 3. The Abstraction Principle for the Checker

The checker proves detailed balance for **all** values of physical parameters simultaneously. Any parameter that appears in the Algorithm body must be an unbound symbol (abstract) during the symbolic BFS. Parameters that appear only in concrete numerical functions (`$couplingJConcrete`, `$fieldFConcrete`) are invisible to the BFS and need not be abstract — the abstract `couplingJ` atom already proves DB for all coupling functions, implicitly covering all physLen and epsLJ values.

**Parameters that must be abstract during BFS:**
- `sigStep` — appears in `NormalDistribution[0, sigStep]` inside Algorithm.
- `physLen`, `epsLJ` — appear only in `$couplingJConcrete`; abstracting them is principled but adds no proof generality. They are included for completeness.

**Parameters that are always Infinity/symbolic by convention:**
- `$maxD2 = Infinity` — set in the algorithm file unconditionally. This makes the energy sum include all site pairs, proving DB for all interaction ranges simultaneously. On the small lattices the checker uses (2×2, 3×3), `Infinity` and any physically reasonable cutoff include exactly the same bonds (all of them). For numerical simulation, the user overrides `$maxD2` to a concrete cutoff before running.
- `β` — always a free symbol in the checker (never assigned `numBeta`).

**Why Gaussian weights also cancel.** For direction `(a, b)` and its reverse `(−a, −b)`, the Gaussian weight `w(a, sigStep) · w(b, sigStep)` is identical on both sides of the DB equation (Gaussian symmetry). It factors out of each direction-grouped sub-expression:

```
w(a, s) · w(b, s) · [T_VMMC(i→j; (a,b)) · exp(−β·E_i) − T_VMMC(j→i; (−a,−b)) · exp(−β·E_j)]
```

The bracketed term is zero by VMMC superdetailed balance. FullSimplify sees `erfc_factor(sigStep) · 0` for each direction group and simplifies trivially. No Erf algebra is required. The proof is general for all `sigStep > 0`.

---

## 4. The `$checkerAbstractParams` Mechanism

Algorithm files declare which parameters the checker must clear before BFS:

```mathematica
$checkerAbstractParams = {physLen, epsLJ, sigStep}
```

**In `check.wls` (after loading the file, before BFS):**

```mathematica
If[ListQ[$checkerAbstractParams],
   Scan[ClearAll, $checkerAbstractParams]]
```

This resets all listed symbols to unbound. After clearing, check.wls verifies they are unbound:

```mathematica
If[AnyTrue[$checkerAbstractParams, ValueQ],
   Print["ERROR: abstract parameter still bound after clearing — check for implicit concrete assignments"]; Exit[1]]
```

The BFS then executes with these as free symbolic atoms. If the Algorithm relies on a concrete value of a cleared parameter (e.g., a conditional branch on `physLen > 3`), the BFS will produce unevaluated or incorrect transitions, which will surface as a DB failure or a Mathematica error. The checker does not attempt static analysis of the Algorithm body for numeric literals; this is documented as the algorithm author's responsibility.

**In `animate.wls`:** No special handling needed. `physLen`, `epsLJ`, and `sigStep` are concrete in Section 0 of the algorithm file; animate.wls uses them as-is.

**`sigStep` requires no special treatment** beyond inclusion in `$checkerAbstractParams`. The earlier plan attempted to keep `sigStep` unbound in the file and store its formula in a separate association — this was necessary only because animate.wls would then have assigned it a random value (step 7). With `sigStep` assigned concretely in Section 0, animate.wls uses the physical formula value directly, and the checker clears it before BFS. The treatment is identical to `physLen` and `epsLJ`.

**`symParams`** must declare the abstract parameters for FullSimplify assumptions:

```mathematica
symParams = <|"eps" -> {sigStep, physLen, epsLJ}, "couplings" -> {}|>
```

This adds `sigStep > 0, physLen > 0, epsLJ > 0` to the assumption set.

---

## 5. Changes to `dbc_core.wl`

### 5.1 Intercepting `RandomVariate[NormalDistribution[...]]`

Inside the checker's BFS execution context, `RandomVariate` is intercepted and converted to `seqBernoulli`:

```mathematica
RandomVariate[NormalDistribution[mu_, sigma_]] :=
  With[{
    nMax = $dbcCurrentNGrid,   (* NMax=Infinity is the global default; see §5.2 *)
    vals = Range[Round[mu] - nMax, Round[mu] + nMax]
  },
  With[{
    rawW = Table[
      CDF[NormalDistribution[mu, sigma], k + 1/2] -
      CDF[NormalDistribution[mu, sigma], k - 1/2],
      {k, vals}]
  },
  seqBernoulli[rawW / Total[rawW], vals]]]
```

When `sigma` is abstract, `rawW` contains symbolic Erfc expressions in `sigma`. `Total[rawW]` is a symbolic sum that normalises correctly; the seqBernoulli BFS works algebraically with symbolic weights.

For `n` values in `vals`, the seqBernoulli generates a BFS sub-tree with `2n − 1` nodes (linear in n, not exponential).

### 5.2 The `NMax = Infinity` default

`NMax = Infinity` is the global default for all `check.wls` invocations. Rather than truncating by sigmas, the Gaussian support per dimension is capped at `floor(nGrid / 2)` — the maximum distinct displacement on a periodic torus of side `nGrid`. This is stored in the global:

```mathematica
$dbcCurrentNGrid = Round[Sqrt[Length[currentState]]];
```

set before each state's BFS traversal. The effective values:

| nGrid | floor(nGrid/2) | Distinct (dx,dy) pairs |
|-------|---------------|----------------------|
| 2 | 1 | 9 |
| 3 | 1 | 9 |
| 4 | 2 | 25 |
| 5 | 2 | 25 |

The BFS trees remain small regardless of `physLen` or `sigStep`. This is correct: the checker verifies algorithm structure on small systems, not large ones.

The `TruncSigmas=N` option (sigma-based truncation) is retained as an alternative for algorithms where NMax=Infinity is insufficient (e.g., testing that sigStep-dependent truncation does not violate DB). Default: `TruncSigmas = 4` when not using `NMax = Infinity`.

### 5.3 The `Round` pass-through

```mathematica
seqBernoulli /: Round[seqBernoulli[w_, v_]] :=
  seqBernoulli[w, Round /@ v]   (* no-op since v is already integer-valued *)
```

This absorbs the `Round` wrapper in `Round[RandomVariate[NormalDistribution[0, sigStep]]]` without modifying the integer-valued seqBernoulli outcomes.

### 5.4 Extension to other distributions

The same template applies to any distribution with an available CDF:

```mathematica
RandomVariate[dist_] :=
  With[{vals = ..., rawW = Table[CDF[dist, k+1/2] - CDF[dist, k-1/2], {k, vals}]},
  seqBernoulli[rawW / Total[rawW], vals]]
```

Priority additions after Normal: `ExponentialDistribution`, `UniformDistribution` (redundant with contToken but useful for API uniformity).

---

## 6. Changes to `vmmc_continuous.wl`

Sections 1–2 and 4–10 are unchanged. Section 3 (K-direction set) is **deleted**. `numDirections`, `stepSize`, and `$dirVectors` are removed.

**Section 0 — User parameters:**

```mathematica
(* ---- Physical parameters — concrete for numerical runs ---- *)
physLen = 5      (* particle diameter in lattice units; positive integer *)
epsLJ   = 1      (* LJ well depth in kT units *)
sigLJ   = physLen

(* Displacement proposal std dev: one BD timestep at the natural LJ timescale.
   Override sigStep directly to tune acceptance rate. *)
sigStep = physLen * Sqrt[2.0 / (numBeta * epsLJ)]

(* Interaction cutoff.
   $maxD2 = Infinity: include all pairs; correct for symbolic check; override for large sims.
   For numerical runs with physLen=5: set $maxD2 = Ceiling[2 * physLen^2] = 50 *)
$maxD2 = Infinity

(* ---- Checker interface ---- *)
(* Parameters cleared to unbound symbols before BFS; concrete values above are for numerical runs *)
$checkerAbstractParams = {physLen, epsLJ, sigStep}

(* Abstract parameters declared as positive reals for FullSimplify assumptions *)
symParams = <|"eps" -> {sigStep, physLen, epsLJ}, "couplings" -> {}|>

(* ---- Concrete coupling (numerical MCMC / animation only) ---- *)
$abstractFunctions = True

$couplingJConcrete[a_Integer, b_Integer, d2_Integer] :=
  If[d2 == 0 || d2 > $maxD2, 0,
     4 * epsLJ * ((sigLJ^2 / d2)^6 - (sigLJ^2 / d2)^3)]

$concreteParams   = <||>
$couplingFormulaStr = "4*epsLJ*((sigLJ^2/d2)^6-(sigLJ^2/d2)^3); min at d2=2^(1/3)*physLen^2"
```

**Section 8 — Algorithm:**

```mathematica
Algorithm[state_List] :=
  Module[{nGrid, occupied, seed, dx, dy, dir, cluster, newState, dest},
    nGrid    = Round[Sqrt[Length[state]]];
    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];

    seed = RandomChoice[occupied];

    dx = Round[RandomVariate[NormalDistribution[0, sigStep]]];
    dy = Round[RandomVariate[NormalDistribution[0, sigStep]]];
    dir = {dx, dy};
    If[dir === {0, 0}, Return[state]];

    cluster = $vmmcBuildCluster[state, nGrid, seed, dir];
    If[cluster === None, Return[state]];

    newState = state;
    Do[newState[[cluster[[i]]]] = 0, {i, Length[cluster]}];
    Do[
      dest = $applyDir[cluster[[i]], dir, nGrid];
      If[newState[[dest]] =!= 0, Return[state, Module]];
      newState[[dest]] = state[[cluster[[i]]]],
      {i, Length[cluster]}];
    newState]
```

The Algorithm body contains no physics parameters other than `sigStep` (cleared to abstract by the checker). The literal `0` in `NormalDistribution[0, sigStep]` is a structural constant (the mean of the symmetric distribution), not a physics parameter.

**DynamicSymParams** is updated:

```mathematica
DynamicSymParams[states_List] :=
  Module[{types, nGrid, d2Vals, couplingAtoms},
    types  = Sort[DeleteCases[Union @@ states, 0]];
    nGrid  = Round[Sqrt[Length[states[[1]]]]];
    d2Vals = Sort @ DeleteDuplicates @ Select[
      Flatten @ Table[$torusD2[s1, s2, nGrid], {s1, nGrid^2}, {s2, nGrid^2}],
      0 < # <= $maxD2 &];
    couplingAtoms = Flatten @ Table[
      If[a <= b, Table[couplingJ[a, b, d2], {d2, d2Vals}], Nothing],
      {a, types}, {b, types}];
    <|"couplings" -> couplingAtoms,
      "numericParams" -> {}|>]
```

`numericParams` is empty because all concrete parameter values are handled by Section 0 directly (and $checkerAbstractParams for the symbolic check). No random assignment via animate.wls step 7 is needed or wanted.

---

## 7. Changes to `check.wls`

**New global default:** `NMax = Infinity` for all runs (no flag required). `$dbcCurrentNGrid` is set before each state's BFS traversal.

**New option:** `TruncSigmas=N` (default 4) — only active when user explicitly overrides `NMax` to a finite value.

**Abstract parameter clearing block** (runs after file load, before BFS):

```mathematica
If[ListQ[$checkerAbstractParams],
   Scan[ClearAll, $checkerAbstractParams];
   If[AnyTrue[$checkerAbstractParams, ValueQ],
      Print["ERROR: parameter remains bound after clearing: ",
            Select[$checkerAbstractParams, ValueQ]];
      Exit[1]]]
```

No other changes to BFS logic, symbolic DB check, or ergodicity check.

---

## 8. What the Checker Now Does

1. **Loads the algorithm file.** Clears all `$checkerAbstractParams` (sigStep, physLen, epsLJ become unbound symbols). `couplingJ` has no DownValues. `$maxD2 = Infinity` (all pairs included). `β` is always symbolic.

2. **Intercepts `RandomVariate[NormalDistribution[0, sigStep]]`** during BFS. With `NMax = Infinity`, produces `seqBernoulli[w, {-1, 0, 1}]` on a 2×2 or 3×3 lattice (floor(nGrid/2) = 1), with symbolic Erfc weights in `sigStep`.

3. **Runs BFS** over all random outcomes: seed particle choice, `(dx, dy)` from the displacement seqBernoulli, and binary link decisions in the cluster builder.

4. **Builds the symbolic transition matrix** T(i→j). Each entry is a product of: symbolic Erfc weights (functions of `sigStep`), and rational link-probability expressions in `β` and abstract `couplingJ[a,b,d2]` atoms.

5. **Checks detailed balance.** FullSimplify verifies `T(i→j)·exp(−β·E_i) − T(j→i)·exp(−β·E_j) = 0` under the assumptions `β > 0`, `sigStep > 0`, `physLen > 0`, `epsLJ > 0`, and all `couplingJ` atoms real. The Erfc weights cancel direction-by-direction (§3); FullSimplify resolves the remaining VMMC algebra identically to the fixed-direction case.

6. **Checks ergodicity** by comparing discovered states to the falling factorial (S)_N.

The proof holds simultaneously for all `sigStep > 0`, all `physLen > 0`, all `epsLJ > 0`, and all coupling functions `couplingJ`.

---

## 9. Algorithm File Contract

Required definitions are unchanged:

```mathematica
energy[state_]       (* pair energy using couplingJ[a, b, d2] *)
Algorithm[state_]    (* MCMC move; may use RandomVariate[NormalDistribution[...]] *)
BitsToState[bits_]
numBeta
```

For algorithms using Gaussian proposals, additionally:

```mathematica
$checkerAbstractParams = {sigStep, ...}   (* parameters cleared to unbound before BFS *)
symParams = <|"eps" -> {sigStep, ...}, "couplings" -> {}|>  (* for FullSimplify assumptions *)
$maxD2 = Infinity   (* default; override for numerical simulation performance *)
```

**Author responsibilities:**
- Concrete literals in Algorithm bodies must be structural constants only (e.g., `0` as a distribution mean, `{0,0}` as the zero-displacement test). Physics parameters must be declared in `$checkerAbstractParams`.
- Gaussian proposals must be symmetric: mean = 0 (or an integer). Non-zero or symbolic means break `p(dx,dy) = p(−dx,−dy)`.
- The checker will detect if a cleared parameter is still bound after clearing (§7), but does not statically analyse the Algorithm body for unlisted numeric literals.

Existing algorithms (`kawasaki_1d.wl`, `vmmc_2d.wl`, etc.) are unaffected. The `RandomVariate` interception is only active inside the BFS execution context and does not run at file-load time.

---

## 10. Implementation Order

1. **`dbc_core.wl`**: Add `RandomVariate[NormalDistribution[...]]` interception and `Round` pass-through UpValue. Add `$dbcCurrentNGrid` global (set before each state's BFS). Make NMax=Infinity the default (no TruncSigmas computation unless explicitly requested). Test with a minimal synthetic algorithm: one particle, `dx = Round[RandomVariate[NormalDistribution[0, sigStep]]]`, verify seqBernoulli weights match known Gaussian CDF differences.

2. **`check.wls`**: Add `$checkerAbstractParams` clearing and verification block. Set `$dbcCurrentNGrid` before each state's BFS. Add `TruncSigmas` as an opt-in override. Test on the synthetic algorithm from step 1.

3. **`vmmc_continuous.wl`**: Rewrite Section 0 and Algorithm per §6. Delete Section 3. Verify coupling symmetry checker passes (exact integer arithmetic from `physLen = 5`). Run checker: `SeedBitStrings=11110101011110011 Mode=Both` and verify PASS.

4. **Validate symbolic cancellation**: On the 3×3 seed, inspect the raw DB expressions for a single state pair to confirm Erfc factors appear identically in numerator and denominator and are cancelled by FullSimplify. If FullSimplify struggles (unlikely given the structure), implement an explicit pre-cancellation step in the FastChecker path.

5. **`animate.wls`**: Verify that concrete `sigStep`, `physLen`, `epsLJ` from Section 0 are used correctly. Verify `$maxD2` is overridden to `Ceiling[2 * physLen^2]` before the animation run (either in the algorithm file as a comment instruction, or added as a command-line option override).

6. **`README.md`**: Update examples table, add `$maxD2` override note to animation section, document `$checkerAbstractParams` convention, document `TruncSigmas` option.
