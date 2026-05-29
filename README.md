# DetailedBalanceChecker

Symbolically proves or disproves detailed balance for Monte Carlo algorithms written in Mathematica. Uses computer algebra (Mathematica's `FullSimplify`) rather than statistical tests alone.

---

## How it works

1. **Randomness → bits.** Every random call is intercepted and reduced to reads from a binary tape:
   - `RandomReal[]` — *interval tracking*: a latent U ∈ [lo,hi]; each comparison `U < p` reads one bit and narrows the interval, weighted by conditional probability.
   - `RandomChoice[list]` / `RandomInteger[{lo,hi}]` — rejection sampling over ⌈log₂(n)⌉ bits.
   - `RandomChoice[weights → elements]` — sequential Bernoulli decomposition.
   - `RandomVariate[NormalDistribution[μ,σ]]` — truncated to `{Round[μ]−nMax,…,Round[μ]+nMax}` (nMax=⌊nGrid/2⌋) via sequential Bernoulli with CDF-difference weights.

2. **BFS over bit sequences.** Starting from a seed state, all possible bit sequences are enumerated. Each complete path gives a (start, end, weight) triple. Weights for the same (start, end) pair are summed to form the symbolic transition matrix T.

3. **Symbolic DB check.** For each state pair (i,j): `T(i→j)·exp(−β E(i)) = T(j→i)·exp(−β E(j))`. β is a free symbol so Boltzmann factors cancel algebraically.

4. **Numerical MCMC check.** The algorithm is run as a genuine Markov chain. KL divergence between sampled frequencies and the Boltzmann distribution is computed; KL < 0.02 is a pass.

**Ergodicity check.** For N labeled particles on S sites, the theoretical reachable state count is (S)_N = S·(S−1)·…·(S−N+1). The number of states discovered by BFS is compared to this value.

---

## Algorithm file structure

Algorithm files for 2D LJ colloidal systems use a three-file split:

```
vmmc_2d_grid.wl          ← shared geometry (torus, encoding, BitsToState)
lj_colloidal_system.wl   ← physical system (LJ energy, parameters)
your_algorithm.wl        ← kinetics (what the LLM writes and optimises)
```

The algorithm file loads the two shared files and defines only `$vmmcBuildCluster`, `Algorithm[]`, and `DynamicSymParams[]`:

```mathematica
$dir = DirectoryName[$InputFileName];
Get[$dir <> "vmmc_2d_grid.wl"];
Get[$dir <> "lj_colloidal_system.wl"];

(* ... define Algorithm[], $vmmcBuildCluster[], DynamicSymParams[] ... *)
```

Simpler 1D algorithms (Kawasaki, jump dynamics) are self-contained single files.

---

## Writing an algorithm file

Every `.wl` file must define:

```mathematica
energy[state_]       (* bare energy, no β factor *)
Algorithm[state_]    (* MCMC move; use RandomReal[], RandomInteger[], RandomChoice[] *)
BitsToState[bits_]   (* {0,1,...} list → seed state, or None if invalid *)
numBeta              (* numeric β for the numerical check *)
```

For Metropolis acceptance use exactly:
```mathematica
If[RandomReal[] < MetropolisProb[dE], newState, state]
```

### Optional definitions

```mathematica
DynamicSymParams[states_List] := ...      (* per-component symbolic parameters *)
DisplayState[state_] := ...               (* human-readable state string *)
ValidStateIDs[maxId_] := ...              (* restrict enumeration to valid IDs *)
$checkerAbstractParams = {"name1", ...}   (* scalar params cleared before BFS *)
$symmetryGroup = {"translation", "D4"}    (* enables G-orbit dedup + G-invariance check — see below *)
```

### Abstract scalar parameters (`$checkerAbstractParams`)

Declare parameters (e.g. `physLen`, `epsLJ`) that must be **unbound symbols** during symbolic BFS so the proof holds for all values simultaneously:

```mathematica
physLen = 1           (* concrete value for numerical runs *)
epsLJ   = 1
$checkerAbstractParams = {"physLen", "epsLJ"}   (* string names — NOT symbol refs *)
```

String names are required because `{physLen, epsLJ}` evaluates to `{1, 1}` before the checker can intercept. The checker saves concrete values, clears the symbols, runs BFS, then restores them in a `Block` for numerical MCMC.

### Abstract coupling functions (`$pairEnergy`)

When `couplingJ` should have no DownValues during BFS (each call treated as a free real atom), provide a concrete implementation as `$pairEnergy`:

```mathematica
$pairEnergy[a_Integer, b_Integer, d2_Integer] :=
  $jPairSym[a, b] * Exp[-lambdaJ * d2]    (* must satisfy J[a,b,d2] = J[b,a,d2] *)
```

`check.wls` / `animate.wls` detect `$pairEnergy` and activate it via `Block[{couplingJ}, couplingJ := $pairEnergy]` for numerical runs. No flag is needed.

---

## physLen / nGrid scaling for LJ colloidal runs

For N particles at packing fraction η on an nGrid×nGrid lattice:

```
physLen = 2 * nGrid * Sqrt[eta / (N * Pi)]
```

Example: η=0.2, N=10, nGrid=20 → physLen ≈ 2.3

Set `$maxD2` to capture the desired interaction range:
- `Ceiling[2*physLen^2]` — cutoff just past the LJ minimum at d²=2^(1/3)·physLen²≈1.26·physLen²
- `Ceiling[6.25*physLen^2]` — standard 2.5σ cutoff (full attractive tail)
- `Infinity` — include all pairs (correct for the checker on small lattices)

---

## Usage

### Checker
```bash
wolframscript -file check.wls <algorithm.wl> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `NGrid=N` | off | Only test N×N grid states; iterates IDs lazily (no large list materialised) |
| `MaxComponents=N` | unlimited | Stop after N distinct connected components are found |
| `MaxBitString=XXXX` | `11111111` | Largest bit string to enumerate (overridden by `NGrid`) |
| `SeedBitStrings=X,Y` | off | Test specific bit strings directly (BFS from each discovers the full component) |
| `Mode=Symbolic\|Numerical\|Both` | `Both` | Which checks to run |
| `NSteps=N` | `50000` | MCMC steps for numerical check |
| `MaxBitDepth=N` | `20` | BFS depth cap per state |
| `FastChecker=1` | off | Exp-polynomial fast checker (see below) |
| `SZChecker=1` | off | Schwartz-Zippel fallback instead of FullSimplify (requires `FastChecker=1`; see below) |
| `SZRepeats=N` | `30` | Number of random rational evaluations per expression for SZChecker |
| `SZOnly=1` | off | Skip FS even for suspected violations; `$dbcFS` cases still go to FS (requires `SZChecker=1`) |
| `SZPure=1` | off | Skip FastChecker entirely; apply SZ directly to all expressions (fastest probabilistic mode; implies `SZOnly`) |
| `PhysFidelity=1` | off | Compare T_MC to physical Glauber dynamics; adds M1 and M2 columns (see below) |
| `Verbose=True` | `False` | Per-state BFS progress |

`NGrid` and `MaxComponents` together give an intuitive interface: `NGrid=2 MaxComponents=6` means "check the first 6 distinct 2×2 systems, stop there." Both options use lazy ID iteration internally so no large list is held in memory, even when the total state count for that grid size is in the millions.

**Output format.** The `ID` column shows the decimal integer corresponding to the seed bit string (more compact than the raw binary). For 2D states the `State` column prints each grid row on its own line, separated by a blank line, so the layout visually matches the lattice:

```
ID       State                 #States  Ergodic       Symbolic  #Fail
41       {0,2}                 3        PASS (24)     PASS      0
         {3,1}
```

### Report (single seed state)
```bash
wolframscript -file report.wls <algorithm.wl> BitString=XXXXX [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `BitString=XXXXX` | *(required)* | Seed bit string |
| `NSteps=N` | `50000` | MCMC steps |
| `ShowTrees=True\|False` | `True` | Include decision trees |
| `ShowT=True\|False` | `True` | Include transition matrix |
| `DoNumerical=True\|False` | `True` | Run numerical MCMC |
| `PrintMatrix=True\|False` | `False` | Print T entries as LaTeX |

### Animation
```bash
wolframscript -file animate.wls <algorithm.wl> Sites=<n> N=<n> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `Sites=N` | *(required)* | Total lattice sites |
| `N=N` | *(required)* | Number of labeled particle types |
| `Steps=N` | `200` | MCMC steps |
| `Beta=f` | from `.wl` file | Inverse temperature |
| `FPS=f` | `10` | Animation frame rate |
| `Simple=1` | off | Fast 2-colour mode (holes vs particles) |
| `RecordEvery=N` | `1` | Record state every N steps |
| `NoParams=1` | off | Hide parameter panel |

Any `name=value` argument not in the table is applied as a Mathematica assignment after the algorithm file loads (e.g. `physLen=2`, `'$maxD2=8'`).

---

## Quick examples

```bash
# Check all 2×2 components of vmmc_continuous.wl (lazy iteration)
wolframscript -file check.wls examples3/vmmc_continuous.wl NGrid=2 Mode=Both

# Check the first 3×3 component only (lazily — no large list built upfront)
wolframscript -file check.wls examples3/vmmc_continuous.wl \
  NGrid=3 MaxComponents=1 Mode=Symbolic FastChecker=1

# Check a specific seed state (BFS discovers the full component)
wolframscript -file check.wls examples3/vmmc_continuous.wl \
  SeedBitStrings=11110101011110011 Mode=Both

# Check 2D Kawasaki (all 2×2 components)
wolframscript -file check.wls examples3/kawasaki_2d.wl NGrid=2 Mode=Symbolic

# Check VMMC with user-defined field (fast polynomial checker)
wolframscript -file check.wls examples3/vmmc_2d_field.wl \
  MaxBitString=1111111111 Mode=Symbolic FastChecker=1

# Check vmmc_lattice — 4 canonical single-particle-type components, symmetry-aware
# (MaxComponents=4 avoids the large 3024-state 4-particle component; wall time ~2:55)
wolframscript -file check.wls examples3/vmmc_lattice.wl \
  NGrid=3 MaxComponents=4 FastChecker=1 Mode=Symbolic

# Schwartz-Zippel fallback: faster than FullSimplify for coupling-heavy expressions
wolframscript -file check.wls examples3/vmmc_lattice.wl \
  NGrid=3 MaxComponents=4 FastChecker=1 SZChecker=1 Mode=Symbolic

# SZPure: skip FastChecker entirely — fastest probabilistic mode, no FullSimplify
wolframscript -file check.wls examples3/vmmc_lattice.wl \
  NGrid=3 MaxComponents=4 SZPure=1 Mode=Symbolic

# Animate on a 20×20 grid (physLen=2, cutoff just past LJ minimum)
wolframscript -file animate.wls examples3/vmmc_continuous.wl \
  Sites=400 N=10 Steps=2000 Beta=1 FPS=8 Simple=1 NoParams=1 physLen=2 '$maxD2=8'

# Report for a specific seed
wolframscript -file report.wls examples3/vmmc_continuous.wl \
  BitString=11110101011110011 ShowTrees=False ShowT=False

# Animate VMMC with field on a 6×6 grid
wolframscript -file animate.wls examples3/vmmc_2d_field.wl \
  Sites=36 N=10 Steps=500 Beta=1 FPS=15 Simple=1
```

---

## Example files (`examples3/`)

| File | Description | DB | Ergodic |
|------|-------------|-----|---------|
| `kawasaki_1d.wl` | 1D Kawasaki nearest-neighbour swap on a periodic ring | PASS | PASS |
| `kawasaki_2d.wl` | 2D Kawasaki on a periodic square lattice | PASS | PASS |
| `vmmc_2d_field.wl` | VMMC with user-defined field and coupling functions | PASS | PASS |
| `vmmc_continuous.wl` | VMMC with Gaussian proposal, LJ energy; three-file structure | PASS | PASS |
| `vmmc_lattice.wl` | VMMC with uniform-box displacement (no Gaussian/Erfc); checker-friendly, converges to continuous limit for large physLen | PASS | PASS |
| `jump_1d_weighted.wl` | 1D jump dynamics using `RandomChoice[weights → ...]` | PASS | PASS |
| `kawasaki_1d_fail.wl` | Sign-reversed dE in Metropolis | FAIL | PASS |
| `kawasaki_1d_nonergodic.wl` | Only type-1 particles move; type-2+ frozen | PASS | FAIL |
| `cluster_1d_fail.wl` | Cluster slides only rightward (asymmetric proposal) | FAIL | FAIL |
| `vmmc_2d_edit.wl` | VMMC with only 3 of 4 directions (asymmetric proposal) | FAIL on 3×3 | — |
| `vmmc_biased.wl` | Intentionally broken: rightward displacements proposed 2× as often (breaks proposal symmetry) | FAIL | PASS |

`kawasaki_1d_nonergodic.wl` demonstrates the key case: **DB PASS + Ergodicity FAIL**. Detailed balance cannot detect non-ergodicity; the two checks are independent.

---

## Symmetry-aware checking

For VMMC algorithms on the periodic square torus whose energy and proposal are invariant under the lattice symmetry group G, G-related state pairs produce identical DB expressions.  The checker exploits this via two features:

1. **Canonical neighbour oracle** — substituted automatically during BFS; makes G-related states produce syntactically identical seqBernoulli trees so `$dbcDedup` collapses them without any extra verification logic.
2. **G-invariance check** — after BFS, verifies T(s→s') = T(g(s)→g(s')) for all state pairs and each group generator.  Black-box: only needs the computed T matrix, never reads algorithm source.

### Architecture: oracle as standard interface

The canonical oracle is exposed through a **standard interface function** `$vmmcCandidates` (defined in `vmmc_2d_grid.wl`).  Algorithm files call `$vmmcCandidates`; the checker substitutes `$dbcCanonicalCandidates` via `Block` during BFS:

```mathematica
(* In check.wls, when $symmetryGroup is declared: *)
$bfsAlg = Function[s, Block[{$vmmcCandidates = $dbcCanonicalCandidates}, $alg[s]]]
```

This is identical in spirit to how `couplingJ` / `$pairEnergy` works: the algorithm always calls the standard function; the checker injects a specialised implementation without modifying the algorithm file.

In real (non-checker) runs, `$vmmcCandidates` returns plain unordered neighbours — correct by Whitelam-Geissler since DB holds for any fixed deterministic ordering.  Canonical ordering is only needed for the checker's dedup speedup.

### How it works (dedup speedup)

The seqBernoulli tree is syntactically determined by the *order* candidates are visited.  `$dbcCanonicalCandidates` sorts by the topological key `(d²(p,q), d²(pPost,q), d²(pRev,q), type)` — all G-invariant scalars.  G-related states therefore produce the same canonical order, the same tree, and hash-identical DB expressions.  `$dbcDedup` collapses the entire G-orbit to one representative evaluated once.

For D4 + translations on a 3×3 grid, |G|=72: the 72-state two-particle component compresses from 288 non-trivial pairs to **2 unique expressions (99% hash-collapsed)**.

### How it works (G-invariance check)

After BFS, `CheckGInvariance` applies group generators to the computed T matrix and checks `T(s→s') === T(g(s)→g(s'))` using syntactic equality (valid because canonical ordering makes G-related entries syntactically identical).

G = (ℤ/nℤ)² ⋊ D4 is generated by **4 generators**: translate(1,0), translate(0,1), rotate90 (90° CW), reflect (main diagonal).  If T is invariant under every generator, it is invariant under all group products — by induction on word length.  Checking 4 generators certifies invariance under all |G| ≤ 72 elements.

The check is O(N) precomputation + O(N²) syntactic comparisons × 4 generators.  Each generator's action on all states is precomputed once into an `Association` (`gMap`), making per-pair lookups O(1) hash reads rather than O(1) genFn calls repeated N² times.  Measured: **4.3 s** for a 504-state component on 4 kernels.

```
  G-orbit dedup: 2/288 unique  (99% hash-collapsed)
  G-invariance: PASS
```

### How to use it

Algorithm files that satisfy the structural conditions below declare:

```mathematica
$symmetryGroup = {"translation", "D4"}   (* or {"translation"} for translations only *)
```

And use `$vmmcCandidates` (not `$dbcCanonicalCandidates`) in `$vmmcBuildCluster`:

```mathematica
cands = $vmmcCandidates[p, pPost, pRev, state, nGrid, inCluster];
```

See `vmmc_lattice.wl` for a complete example.

### Conditions for a valid declaration

1. Seed selection is `RandomChoice[occupied]` — uniform over occupied sites, no position-dependent weighting.
2. Displacement list contains every `(dx,dy)` paired with `(-dx,-dy)` (proposal symmetry `P(d)=P(-d)`).
3. Energy uses only minimum-image periodic distances — no external field.

If any condition fails, the G-invariance check will report FAIL, telling you the declaration is wrong.  The DB check is unaffected (it remains correct regardless).

### Safety guarantees

| Scenario | DB check | G-inv check | Dedup |
|----------|----------|-------------|-------|
| Algorithm satisfies DB and G (correct, correct declaration) | PASS | PASS | ~\|G\|-fold speedup |
| Algorithm violates DB | FAIL | PASS or FAIL | No spurious collapse (non-G-inv expressions hash-distinct) |
| Algorithm satisfies DB but not G (rare) | PASS | FAIL | Low dedup ratio |
| Wrong `$symmetryGroup` declaration | PASS | FAIL | Warning printed |

`$dbcDedup` only collapses **syntactically identical** expressions.  A non-G-invariant algorithm (different T entries for G-orbit pairs) produces hash-distinct DB expressions → no collapse → both pairs evaluated separately → correct result.  G-invariance failure cannot produce false negatives in the DB check.

### vmmc_biased.wl: a counterexample

`vmmc_biased.wl` duplicates rightward displacements so `P(dx>0,dy) = 2·P(dx<0,dy)`.  This breaks DB (correctly caught: FAIL).  D4 is broken — a 90° rotation maps rightward moves to upward moves, which are unbiased.  If `$symmetryGroup = {"translation", "D4"}` were incorrectly declared:

- With the **old** approach (algorithm directly calls `$dbcCanonicalCandidates`): the canonical oracle makes the cluster-building part G-invariant, but the direction-choice part still has different weights for rightward vs upward moves.  T(s→s') ≠ T(rot(s)→rot(s')) syntactically → different hashes → `$dbcDedup` does NOT collapse them → both evaluated → FAIL correctly reported.
- G-invariance check would also report FAIL (rotate90 generator violated), warning that the declaration is incorrect.

The G-invariance check is the definitive diagnostic: if it fails, the declared symmetry is wrong and the speedup did not apply.

---

## Checker trust model and known failure modes

The checker is a formal verifier, but it operates on a *model* of the algorithm's randomness.  Understanding what it trusts and where it can fail is essential.

### What is verified (black-box)

- **Ergodicity**: state count from BFS vs theoretical (S)_N.
- **Detailed balance**: FullSimplify proves T(i→j)π(i) = T(j→i)π(j) for every pair, with β free.
- **G-invariance** (when `$symmetryGroup` declared): T(s→s') === T(g(s)→g(s')) for all generators and all state pairs.
- **Coupling symmetry**: couplingJ(a,b,d²) = couplingJ(b,a,d²) via Simplify.

### What is trusted (not verified)

1. **Algorithm determinism given a bit tape.** The checker feeds a fixed bit sequence and records the output. If the algorithm uses mutable global state, system time, or external calls, the recorded path may differ from actual execution. `CheckAlgorithmSafety` detects known unsafe calls (`AbsoluteTime`, `RandomWord`, etc.) but cannot detect custom C functions or I/O.

2. **Energy function correctness.** The checker uses `energy[state]` exactly as written. An incorrect energy (wrong formula, wrong units) invalidates both the symbolic and numerical checks simultaneously. `CheckEnergySafety` detects random calls inside the energy function.

3. **`$symmetryGroup` accuracy** (prior to G-invariance check). Before the G-invariance check runs (per component), the oracle substitution trusts that the algorithm genuinely uses `$vmmcCandidates`. After the check passes, this is verified.

### Known false-negative risks

| Source | Description | Mitigation |
|--------|-------------|------------|
| `MaxBitDepth` truncation | Paths requiring >20 bits are silently excluded. Default covers 8-choice displacement + typical VMMC depth; increase if algorithm has deep paths. | `MaxBitDepth=N` option; WARNING printed when triggered |
| `TimeLimit` truncation | BFS exits after 120 s per state; long paths dropped. | WARNING printed; increase via `BuildTreeAT` TimeLimit |
| Hash collisions in `$dbcDedup` | Two distinct expressions with the same hash → one evaluated for both. Probability negligible (Mathematica 64-bit hash) but non-zero. | No mitigation; practically impossible |
| `MaxComponents` limit | Only the first N components are checked; broken components beyond that are missed. | Set `MaxComponents` appropriately; run with `Mode=Symbolic` first |
| Hidden mutable global state | Algorithm caches results in DownValues across calls; BFS path for state A contaminates path for state B. | Write stateless algorithms; `CheckAlgorithmSafety` catches some patterns |
| Distance-matrix degeneracies | Two distinct geometric configurations with identical pairwise distances produce the same canonical T expression; G-invariance check sees them as equal. Only relevant for multi-particle systems with near-regular arrangements. | Test multiple components; inspect with `FailFast=1` |
| Self-symmetric states (stabilizers) | States invariant under some g (e.g. cluster at every corner of a 2×2 square) trivially pass T(s→s') = T(g(s)→g(s')). G-invariance check cannot distinguish "invariant because G-invariant" from "invariant because stabilizer". | No false negative: stabilizer states correctly contribute zero violations |
| SZChecker false-zero (probabilistic) | With k=30 random evaluations and random range ±p/q, p,q∈[1,50], a non-zero degree-d coupling polynomial evaluates to zero at all k points with probability ≤ (d/50)³⁰. For d≤20: ≤ 10⁻¹². | Increase `SZRepeats` for extremely high assurance; default k=30 is sufficient for all practical purposes |
| SZPure false-zero (probabilistic) | Same as SZChecker above, applied to every expression. The transcendental abstraction (`$dbcAbstractTrans`) handles Erf/Erfc/Bessel coefficients exactly, so there is no additional risk from those terms. Only the pure polynomial coupling part is probabilistic. | Same mitigation; `SZPure` routes `$dbcFS` (structurally anomalous) expressions to FullSimplify as final safety net |

### Known false-positive risks (spurious failures)

| Source | Description | Mitigation |
|--------|-------------|------------|
| `FullSimplify` unable to prove zero | Rare; happens when Piecewise implicit conditions are complex. `$dbcCheckOneExpr` with `deepCheck=True` is applied as a second pass to catch these. | FastChecker + FullSimplify fallback handles nearly all cases |
| Abstract params not cleared | If `$checkerAbstractParams` is incomplete, a parameter retains its numeric value during BFS, reducing the symbolic check to a single-point evaluation. | Check that all free parameters are listed; run `Verbose=True` to inspect |

## Symbolic checkers

### Standard (default)
`FullSimplify[PiecewiseExpand[expr], {β > 0, ...}]`. Correct for all algorithm types.

### FastChecker (`FastChecker=1`)
After `PiecewiseExpand`, DB expressions reduce to sums of `c·exp(−β·L)` terms. DB holds iff all coefficient groups sum to zero — verified by `Expand[...] === 0` (microseconds). Falls back to `FullSimplify` (or SZChecker if enabled) for inconclusive cases. Works directly for Metropolis acceptance; falls back gracefully for Barker/heat-bath.

All checkers share the same efficiency pipeline: trivial-zero filter → syntactic deduplication (`$dbcDedup`; collapses G-orbit pairs when canonical oracle is active) → threshold-based `ParallelMap`.

### SZPure checker (`SZPure=1`)

The fastest probabilistic mode. Bypasses FastChecker's Piecewise case enumeration entirely and applies Schwartz-Zippel directly to every deduplicated DB expression:

1. For each unique expression, substitute `k` random rationals for all coupling parameters and abstract scalar parameters (e.g. step size σ); β remains symbolic.
2. Abstract any transcendental function values (`Erf[concrete]`, `Erfc[concrete]`, Bessel functions, etc.) to fresh algebraic symbols — same argument gets the same symbol. This makes Erf/Erfc cancellations exact without FullSimplify.
3. `$dbcIsExpZero` checks whether every Exp-basis coefficient (now a rational polynomial in the Erf symbols) vanishes exactly.
4. All `k` evaluations confirm zero → PASS. Any non-zero result → violation. Unexpected structure → FullSimplify fallback.

`SZPure=1` implies `FastChecker=1` and `SZOnly=1` (automatically enabled). Use `SZRepeats=N` to control `k` (default 30).

**When to use.** SZPure works for all algorithm types in this library, including Gaussian-proposal algorithms (`vmmc_continuous`). Transcendental abstraction handles Erf/Erfc coefficients without FullSimplify. For coupling-polynomial algorithms it is trivially fast; for Gaussian-proposal algorithms it is also fast because Mathematica's evaluator auto-applies `Erf[-x]→-Erf[x]` before the check, canonicalising all arguments so same-argument cancellations appear directly.

**False-negative risk.** The probabilistic false-zero rate is ≤ (d/50)^k per expression (d = coupling polynomial degree). At d≤20, k=30: < 10⁻³⁸. Increase `SZRepeats` for mission-critical audits.

### Schwartz-Zippel checker (`FastChecker=1 SZChecker=1`)

When FastChecker cannot decide an expression, SZChecker replaces the `FullSimplify` fallback. The Schwartz-Zippel approach:

1. Draw `k` independent random rational assignments for the coupling parameters (e.g. `couplingJ[a,b,d²]`).
2. For each assignment, substitute into the DB expression and call `PiecewiseExpand[·, β>0]`. All Piecewise conditions that compare coupling values collapse to True/False; β remains symbolic.
3. Call `$dbcIsExpZero` on the result — a rational-coefficient check on the Exp-basis, microseconds per call.
4. If all `k` evaluations confirm zero: expression is certified identically zero (PASS). If any evaluation is non-zero: violation detected (FAIL). If the structure is unexpected after substitution, fall back to `FullSimplify` for that expression.

**Correctness.** The expression is a rational polynomial in the coupling atoms times Exp[−β·r] factors; for a non-zero polynomial, the Schwartz-Zippel lemma guarantees that a random rational point is non-zero with probability ≥ 1−d/N, where d is the polynomial degree and N=50 is the random range. With k=30 independent evaluations the false-zero probability is ≤ (d/50)^30 ≤ 10⁻¹² per expression (d≤20), effectively zero.

**Transcendental extension.** After coupling substitution, function arguments become concrete numerics. `$dbcAbstractTrans` replaces each unique `Erf[concrete]`, `Erfc[concrete]`, Bessel function value, etc. with a fresh algebraic symbol (same argument → same symbol). `$dbcIsExpZero` then treats the Erf symbols as free polynomial variables, making same-argument cancellations exact. Mathematica's built-in `Erf[-x]→-Erf[x]` identity fires at substitution time, so forward and backward Erf arguments are automatically canonicalised before abstraction — the checker never needs to know this identity explicitly.

**When it helps.** SZChecker is fast when many unique expressions survive FastChecker — each evaluation costs ~1 ms (substitution + abstraction + rational coefficient check) versus ~1–30 s for `FullSimplify`. For VMMC on a 3×3 grid (~60 unique post-FastChecker expressions), SZChecker adds negligible overhead; for larger grids where FullSimplify is the bottleneck, it provides the largest speedup.

**When it makes no difference.** If FastChecker resolves all expressions (common for small components), the fallback is never reached and SZChecker adds no cost.

**Compatibility.** Requires `FastChecker=1`. The `SZRepeats=N` flag sets k (default 30). Expressions where SZ returns `$dbcFS` still fall back to `FullSimplify`, so coverage is never lost.

### What still requires FullSimplify

After transcendental abstraction, `$dbcSZCheckOne` falls back to `FullSimplify` only when it returns `$dbcFS` — i.e., when `$dbcIsExpZero` cannot determine the structure. Remaining cases:

| Situation | Why $dbcFS | Likelihood |
|-----------|------------|------------|
| Term with 2+ distinct `Exp[...]` factors after `$dbcMergeExp` | `$dbcSplitTerm` sees multiple Exp subexpressions and cannot split cleanly | Rare; would require an expression like `Exp[a]·Exp[b]` that `$dbcMergeExp` failed to merge. Not seen in practice. |
| `Sin`/`Cos`/`Tan` of numeric argument in coefficient | Excluded from abstraction (algebraic relations between different argument values; sin²+cos²=1 etc.) | Not seen in any current algorithm. If needed, would require a dedicated trigonometric extension. |
| Residual `Piecewise` after `PiecewiseExpand` | `$dbcSplitTerm` encounters Piecewise at the top level, finds Exp factors in both branches, returns `$dbcFS` | Should not occur when all coupling params and abstract params are substituted; only possible if a conditional depends on β symbolically in a way `PiecewiseExpand` cannot resolve. |
| Custom transcendental not in the abstraction list | Any function not in `{Erf, Erfc, FresnelS/C, SinIntegral, CosIntegral, ExpIntegralEi, LogIntegral, BesselJ/Y/I/K, ExpIntegralE}` | Extend the list in `$dbcAbstractTrans` as needed. |

In practice, all algorithms in `examples3/` now have zero `$dbcFS` fallbacks when using `SZChecker=1` or `SZPure=1`.

### FastChecker internals: Piecewise case feasibility

`$dbcCheckOneExpr` iterates Piecewise cases and checks each one under its *full* implicit condition (raw condition AND NOT all prior cases). A case can only contribute a violation if its full condition is satisfiable. This feasibility check (`$dbcFeasible`) previously used Mathematica's `Reduce` over the reals, which performs full quantifier elimination and is slow (~0.95 s) on conditions with `!=` constraints over 5+ real-valued variables — causing almost every expression to time out and fall back to `FullSimplify`.

Two fixes applied:

**Quick path.** A conjunction of pure `Unequal` atoms between distinct expressions is always satisfiable over the reals (each `a != b` constraint has a measure-1 solution set; no finite conjunction is infeasible). Detected in O(n) with pattern matching; returns True in under 2 ms.

**`FindInstance` fallback.** For conditions that don't hit the quick path, `FindInstance` locates one example rather than characterising the full solution set. This is typically 10–100× faster than `Reduce` for polynomial inequality conditions.

### FastChecker internals: implicit equality extraction

When `FullSimplify` produces a simplified residual that the fast path still cannot decide, `$dbcCheckOneExpr` is re-run with `deepCheck=True` on the simplified expression. In this mode, `LogicalExpand` expands `Not[A ∧ B] → ¬A ∨ ¬B` before `Simplify`, exposing implicit equalities hidden inside negated conjunctions (e.g. `J₁ < J₂ ∧ J₃ = J₄` implied by elimination of prior cases). `$dbcSubstEqualities` then substitutes those equalities into the case value; algebraic cancellation confirms zero without further `FullSimplify` calls.

This was required to avoid a false positive in `vmmc_continuous.wl` at NGrid=3: `FullSimplify` treated each Piecewise case independently and missed the implicit `J₂₃₁ = J₂₃₂` constraint implied by the prior-case exclusion, reporting 1296 spurious violations.

### Performance: parallel BFS

`BuildTreeAT` uses level-synchronous parallel BFS: the entire BFS wave frontier is sent to `ParallelMap` in one call, partitioned into ≤ nK chunks.  This gives at most one round-trip overhead per wave, amortised across all states in the wave.  The parallel path is taken whenever `nK > 0` and the wave has ≥ 2 states; single-state waves use the sequential path automatically.

Measured on a 504-state component (vmmc_lattice.wl, 3×3, 4 kernels):

| Mode | BFS time |
|------|----------|
| Sequential | 96.3 s |
| Parallel (4 kernels) | 28.9 s — **3.3× speedup** |

### Performance: parallel FullSimplify

When the fast check is inconclusive for many unique expressions, the `FullSimplify` fallback previously ran sequentially on the main kernel — the dominant bottleneck for large components. `CheckDetailedBalanceFast` now runs a second `ParallelMap` over all unique expressions needing `FullSimplify` before the violation scan, which then performs only fast Association lookups. Measured speedup for `vmmc_continuous.wl` at NGrid=3, 504-state component: 29 min → 12.5 min on 4 kernels.

---

## Physical fidelity metrics (`PhysFidelity=1`)

Adds two columns to the checker output that measure how closely the algorithm's transition matrix T_MC matches the physical reference dynamics.

### Reference: T_phys (Glauber single-particle)

For a 2D periodic square lattice with N particles on an nGrid×nGrid torus, the physical reference is Glauber single-particle dynamics. For each state pair (i, j) differing by exactly one nearest-neighbor hop (torus distance d²=1, same particle type):

```
T_phys(i→j) = (1/N) · (1/z) · 1/(1 + exp(β·ΔE))
```

where z=4 (2D coordination number) and ΔE = E(j)−E(i). Diagonal entries make rows sum to 1. All other entries are zero. T_phys is built exactly from the enumerated state space — no simulation required.

This is the correct physical kinetics for a single colloid undergoing thermally activated nearest-neighbor hops. It satisfies detailed balance by construction.

### Metric 1 — Eigenvalue spectrum ratio

```
M1 = |λ₂^MC/λ₃^MC − λ₂^phys/λ₃^phys|
```

Eigenvalues are sorted by real part (descending); λ₁=1 for any ergodic chain. M1 measures whether the algorithm reproduces the correct *hierarchy* of relaxation timescales. M1=0 means the ratio of the two slowest relaxation rates matches Glauber exactly; larger M1 means the algorithm's dynamics have a distorted timescale structure.

Returns `N/A` when fewer than 3 states or when λ₃ is near zero.

### Metric 2 — Row-normalised KL divergence (fidelity score)

```
F = −Σ_i π_i Σ_j T_phys(i→j) · log(T_phys(i→j) / T_MC(i→j))
```

where π_i ∝ exp(−β E(i)) is the Boltzmann weight. F=0 is a perfect match; F<0 indicates the algorithm departs from physical dynamics (more negative = further from physical). Returns `-Inf` (hard failure) if T_MC(i→j)=0 for any transition where T_phys(i→j)>0 — the algorithm completely misses a physically required move.

### Interpretation

| Algorithm | Expected M2 | Reason |
|-----------|-------------|--------|
| Glauber single-particle | 0 | Exact match |
| Kawasaki (2D swap) | slightly negative | Swaps two particles; different proposal distribution |
| VMMC | more negative | Cluster moves; completely different kinetic pathway |

M2 is not a pass/fail criterion — it quantifies the *kinetic fidelity* of the algorithm relative to Brownian dynamics. A more negative M2 for VMMC is expected and reflects the algorithmic acceleration: VMMC decorrelates faster than physical dynamics by construction.

### Scope

`PhysFidelity=1` is only meaningful for 2D square-lattice systems (returns `N/A` for 1D or non-square lattices). It requires the energy function to evaluate numerically — the same environment used for the numerical MCMC check.

---

## Supported random calls

| Call | Behaviour |
|------|-----------|
| `RandomReal[]` / `Random[]` | Interval-tracking token |
| `RandomInteger[{lo, hi}]` | Rejection sampling over ⌈log₂(hi−lo+1)⌉ bits |
| `RandomChoice[list]` | Uniform choice via rejection sampling |
| `RandomChoice[weights → elements]` | Sequential Bernoulli decomposition |
| `RandomVariate[NormalDistribution[μ, σ]]` | Truncated Gaussian, sequential Bernoulli, nMax=⌊nGrid/2⌋ |
| `RandomPermutation[n]` / `RandomPermutation[list]` | Knuth shuffle via intercepted `RandomInteger[]` |
| `RandomSample[list]` / `RandomSample[list, k]` | Partial Knuth shuffle |

---

## 2D encoding (dbc_core.wl)

The bijective integer encoding that maps bit strings to 2D square-lattice states lives at the end of `dbc_core.wl` (tagged "2D SQUARE LATTICE SPECIFIC"). For 3D or non-square geometries, replace `$decode` and update `BitsToState` in the algorithm file. The BFS core is geometry-agnostic.
