# DetailedBalanceChecker

Symbolically proves or disproves detailed balance for Monte Carlo algorithms written in Mathematica. Uses computer algebra rather than statistical tests alone.

---

## How it works

1. **Randomness → bits.** Every random call is intercepted and reduced to reads from a binary tape:
   - `RandomInteger[{lo,hi}]` / `RandomChoice[list]` — rejection sampling over ⌈log₂(n)⌉ bits.
   - `RandomReal[]` — *interval tracking*: a latent variable U ∈ [lo,hi]; each comparison `U < p` reads one bit and narrows the interval to [lo,p] or [p,hi], with weight (p−lo)/(hi−lo) or (hi−p)/(hi−lo) respectively.
   - `RandomChoice[weights → elements]` — decomposed into n−1 sequential Bernoulli trials.

2. **BFS over bit sequences.** Starting from a seed state, all possible bit sequences are enumerated. Each complete path gives a (start, end, weight) triple. Weights for the same (start, end) pair are summed to form the symbolic transition matrix T.

3. **Symbolic detailed-balance check.** For each state pair (i, j):
   ```
   T(i→j)·exp(−β E(i)) − T(j→i)·exp(−β E(j))  =?  0
   ```
   β is kept as a free symbol so Boltzmann factors cancel algebraically.

4. **Numerical MCMC check.** The algorithm is run as a genuine Markov chain; the empirical distribution is compared to the Boltzmann distribution via KL divergence. KL < 0.02 is a pass.

---

## Ergodicity check

In addition to detailed balance, the checker tests **ergodicity** (irreducibility of the Markov chain) using state counting. For N labeled particles on S lattice sites with hard-sphere exclusion, the theoretical number of reachable states is the falling factorial:

(S)_N = S · (S−1) · … · (S−N+1)

Both S and N are inferred from any discovered state; no external parameters are needed. After BFS, the checker compares the number of discovered states to this theoretical value. A shortfall indicates the chain is confined to a subset of configuration space — a form of non-ergodicity that detailed balance cannot detect.

The ergodicity result appears as a column in the output table (`PASS (k)` or `FAIL (k/m)` where k is found and m is theoretical) and is summarised separately at the end of the run.

**Key point:** ergodicity and detailed balance are independent. `kawasaki_1d_nonergodic.wl` demonstrates this — it passes detailed balance but fails ergodicity because type-2+ particles are permanently frozen, creating disconnected sectors of configuration space.

## Continuous-limit VMMC (`vmmc_continuous.wl`)

`vmmc_continuous.wl` generalises `vmmc_2d_field.wl` to support K uniformly-spaced translation directions and an extended interaction cutoff $maxD2, so that the same code describes both checkable lattice systems and large fine-grained systems that converge to continuous off-lattice VMMC comparable to MPCD+MD.

**Two continuous-limit control parameters:**
- `numDirections` K — number of directions uniformly spaced over [0, 2π). K=4: compass (checker default); K=8: compass+diagonal; larger K approaches isotropic diffusion as nGrid grows.
- `nGrid` — inferred from the state as `Round[Sqrt[Length[state]]]`. Increasing nGrid reduces lattice spacing b=L_box/nGrid; as nGrid→∞ with fixed K, the discrete walk approaches Brownian motion.

**Energy:** Abstract `couplingJ[type1, type2, d2]` summed over all site pairs within squared grid distance ≤ `$maxD2`. During the symbolic check, `couplingJ` has no DownValues — each `(a,b,d2)` triple is a free real atom, proving detailed balance for all coupling functions simultaneously. A concrete Lennard-Jones-like implementation activates for numerical MCMC runs.

**VMMC acceptance:** Rigid cluster translation — intra-cluster distances are preserved, so ΔE_intra=0 and no post-cluster Metropolis step is needed. Superdetailed balance is satisfied by the Whitelam-Geissler link-probability mechanism.

**Physical comparison:** At large nGrid and K, dimensionless dynamical exponents (cluster diffusion D∝n^{-ν}, MSD scaling α, intermediate scattering function shape) can be compared to MPCD+MD reference data to grade physical fidelity.

```bash
# Check vmmc_continuous.wl (K=4, 3×3 seed, symbolic + numerical)
wolframscript -file check.wls examples3/vmmc_continuous.wl \
  SeedBitStrings=11110101011110011 Mode=Both
```

## Symbolic checkers

Two symbolic checking methods are available, selectable per run.

Both checkers share the same efficiency pipeline before simplification:
- **Trivial-zero filter:** pairs where both transition weights are zero are skipped immediately without any algebra.
- **Syntactic deduplication:** the remaining expressions are hashed; each unique expression is simplified only once and the result broadcast to all duplicate pairs. For periodic lattices this gives a speedup proportional to the translational symmetry (e.g. 9× for a 3×3 grid, reducing ~2500 pairs to ~280 unique ones).
- **Threshold-based parallelism:** if the number of unique expressions exceeds the kernel count, `ParallelMap` distributes work across subkernels with per-expression scheduling for natural load balancing; otherwise evaluation is sequential to avoid dispatch overhead.

### Standard checker (default)
Uses `FullSimplify[PiecewiseExpand[expr], {β > 0, ...}]`. Correct for all algorithm types.

### FastChecker (`FastChecker=1`)
An Exp-polynomial algebraic checker. After `PiecewiseExpand`, detailed balance expressions for Boltzmann algorithms reduce — within each feasible sign case — to sums of terms `c·exp(−β·L)` where c is a rational coefficient and L is a polynomial in coupling constants. Detailed balance holds iff all coefficient sums vanish, verified by `Expand[...] === 0` (microseconds).

**Algorithm:**
1. `PiecewiseExpand` splits the expression into sign cases (e.g. `Jpair12 > 0`).
2. Each case is tested for feasibility (`Reduce` on linear real arithmetic); infeasible cases are skipped.
3. Equality constraints (e.g. `Jpair12 == 0`) are substituted into the value.
4. The value is normalised: `(E^a)^n → E^(n·a)`, products of `Exp` are merged, and the result is expanded to a sum of `c·E^(poly)` terms.
5. Terms are grouped by exponent polynomial. Detailed balance holds iff all coefficient groups sum to zero.
6. **Fallback:** if any step is inconclusive, the expression falls back to `FullSimplify`. The result is always exact.

**Applicability:** Works directly for all algorithms using Metropolis acceptance. Falls back gracefully for other acceptance functions (Barker, heat-bath) or unusual expression structures.

**Timings** (4-core machine, `MaxBitString=1111111111`, `Mode=Symbolic`):

| Algorithm | Standard | FastChecker |
|-----------|----------|-------------|
| `vmmc_2d.wl` | ~58s | ~50s |
| `vmmc_2d_field.wl` | ~98s | ~59s |

FastChecker is faster for `vmmc_2d_field.wl` because its abstract field and coupling terms produce expressions that `FullSimplify` must work hard on but which the polynomial coefficient check resolves trivially.

---

## Usage

### Checker
```bash
wolframscript -file check.wls <algorithm.wl> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `MaxBitString=XXXX` | `11111111` | Largest bit string tested |
| `SeedBitStrings=X,Y` | off | Test specific bit strings only (comma-separated); bypasses `MaxBitString` enumeration. BFS from each seed discovers the full connected component. |
| `Mode=Symbolic\|Numerical\|Both` | `Both` | Which checks to run |
| `NSteps=N` | `50000` | MCMC steps for numerical check |
| `MaxBitDepth=N` | `20` | BFS depth cap per state |
| `FastChecker=1` | off | Use Exp-polynomial fast checker (see above) |
| `Verbose=True` | `False` | Print per-state BFS progress |

### Report (single seed state)
```bash
wolframscript -file report.wls <algorithm.wl> BitString=XXXXX [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `BitString=XXXXX` | *(required)* | Seed bit string, e.g. `11001` |
| `NSteps=N` | `50000` | MCMC steps |
| `MaxBitDepth=N` | `20` | BFS depth cap |
| `ShowTrees=True\|False` | `True` | Include decision trees in report |
| `ShowT=True\|False` | `True` | Include transition matrix in report |
| `DoNumerical=True\|False` | `True` | Run numerical MCMC and include scatter plot |
| `PrintMatrix=True\|False` | `False` | Print transition matrix entries as LaTeX |

### Animation
```bash
wolframscript -file animate.wls <algorithm.wl> Sites=<n> N=<n> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `Sites=N` | *(required)* | Total lattice sites |
| `N=N` | *(required)* | Number of labeled particle types |
| `Steps=N` | `200` | MCMC steps to run |
| `Beta=f` | from `.wl` file | Inverse temperature |
| `FPS=f` | `10` | Animation frame rate |
| `Simple=1` | off | Fast 2-colour mode (holes vs particles) |
| `RecordEvery=N` | `1` | Record state every N steps |
| `Jpair<a><b>=f` | random | Coupling constant, e.g. `Jpair12=-1.0` |
| `NoParams=1` | off | Hide the right-hand parameter panel (faster rendering on large systems) |

---

## Algorithm file format

Each `.wl` file must define:

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
where `dE = energy[newState] - energy[state]`.

### Optional definitions

```mathematica
symParams = <|"eps" -> {...}, "couplings" -> {...}|>
DynamicSymParams[states_List] := ...   (* per-component symbolic parameters *)
DisplayState[state_] := ...            (* human-readable state string *)
ValidStateIDs[maxId_] := ...           (* restrict enumeration to valid IDs *)
```

### Abstract field/coupling functions (`vmmc_2d_field.wl` style)

`vmmc_2d_field.wl` separates symbolic and numerical modes:

**Symbolic check:** `fieldF` and `couplingJ` have no DownValues. Each evaluation is treated as an independent free real atom, proving detailed balance for *all possible* functions simultaneously.

**Numerical MCMC / animation:** Concrete implementations (`$fieldFConcrete`, `$couplingJConcrete`) are activated inside a `Block`.

To customise, edit **Section 0 only** in `vmmc_2d_field.wl`:

```mathematica
$fieldFConcrete[x_Integer, y_Integer, L_Integer] :=
  fieldAmp * Sin[Pi * (x + 1) / L]

$couplingJConcrete[a_Integer, b_Integer, d2_Integer] :=
  $jPairSym[a, b] * Exp[-lambdaJ * d2]    (* must satisfy [a,b,d2] = [b,a,d2] *)

$concreteParams = <|fieldAmp -> 1.0, lambdaJ -> 0.5|>
$maxD2 = 1        (* max squared bond distance *)
$abstractFunctions = True
```

---

## Supported random calls

| Call | Behaviour |
|------|-----------|
| `RandomReal[]` / `Random[]` | Interval-tracking token |
| `RandomInteger[{lo, hi}]` | Rejection sampling over ⌈log₂(hi−lo+1)⌉ bits |
| `RandomInteger[n]` | Uniform over {0,…,n} |
| `RandomChoice[list]` | Uniform choice via rejection sampling |
| `RandomChoice[weights → elements]` | Sequential Bernoulli decomposition |

`RandomVariate`, `RandomSample`, `RandomPermutation` are **not** supported.

---

## Example files (`examples3/`)

| File | Description | DB | Ergodic |
|------|-------------|----|----|
| `kawasaki_1d.wl` | 1D Kawasaki (nearest-neighbour swap) on a periodic ring | PASS | PASS |
| `kawasaki_2d.wl` | 2D Kawasaki on a periodic square lattice | PASS | PASS |
| `vmmc_2d.wl` | Virtual Move Monte Carlo (Whitelam–Geissler) on a 2D torus | PASS | PASS |
| `vmmc_2d_field.wl` | VMMC with user-defined field and coupling functions | PASS | PASS |
| `jump_1d_weighted.wl` | 1D jump dynamics using `RandomChoice[weights → ...]` | PASS | PASS |
| `kawasaki_1d_fail.wl` | Sign-reversed dE in Metropolis | FAIL | PASS |
| `kawasaki_1d_nonergodic.wl` | Kawasaki where only type-1 particles move; type-2+ frozen | PASS | FAIL |
| `cluster_1d_fail.wl` | Cluster slides only rightward (asymmetric proposal) | FAIL | FAIL |
| `vmmc_2d_edit.wl` | VMMC with only 3 of 4 directions (asymmetric proposal) | FAIL on 3×3, PASS on 2×2 | — |
| `vmmc_continuous.wl` | VMMC with K uniformly-spaced directions and extended interaction cutoff $maxD2; same code is the continuous-limit reference | PASS | PASS |

---

## Quick examples

```bash
# Check 2D Kawasaki symbolically + numerically
wolframscript -file check.wls examples3/kawasaki_2d.wl

# Check VMMC with fast polynomial checker
wolframscript -file check.wls examples3/vmmc_2d_field.wl \
  MaxBitString=1111111111 Mode=Symbolic FastChecker=1

# Target a specific 3×3 state directly (bit string for {1,2,0,...} on a 3×3 grid)
# Avoids enumerating ~125k bit strings to reach 3×3 states via MaxBitString
wolframscript -file check.wls examples3/vmmc_2d.wl \
  SeedBitStrings=11110101011110011 Mode=Symbolic

# Report for a specific seed state
wolframscript -file report.wls examples3/vmmc_2d.wl BitString=101001 \
  ShowTrees=False ShowT=False

# Animate VMMC with field on a 6×6 grid
wolframscript -file animate.wls examples3/vmmc_2d_field.wl \
  Sites=36 N=10 Steps=500 Beta=1 FPS=15 Simple=1
```
