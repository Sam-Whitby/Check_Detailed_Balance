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

3. **Symbolic DB check.** For each state pair (i,j): `T(i→j)·exp(−β E(i)) − T(j→i)·exp(−β E(j)) =? 0`. β is a free symbol so Boltzmann factors cancel algebraically.

4. **Numerical MCMC check.** The algorithm runs as a genuine Markov chain; empirical distribution vs Boltzmann distribution via KL divergence. KL < 0.02 is a pass.

**Ergodicity check.** For N labeled particles on S sites, the theoretical reachable state count is (S)_N = S·(S−1)·…·(S−N+1). A BFS shortfall indicates a non-ergodic chain.

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
DynamicSymParams[states_List] := ...   (* per-component symbolic parameters *)
DisplayState[state_] := ...            (* human-readable state string *)
ValidStateIDs[maxId_] := ...           (* restrict enumeration to valid IDs *)
$checkerAbstractParams = {"name1", ...} (* scalar params cleared before BFS *)
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
| `MaxBitString=XXXX` | `11111111` | Largest bit string tested |
| `SeedBitStrings=X,Y` | off | Test specific bit strings directly (BFS from each discovers the full component) |
| `Mode=Symbolic\|Numerical\|Both` | `Both` | Which checks to run |
| `NSteps=N` | `50000` | MCMC steps for numerical check |
| `MaxBitDepth=N` | `20` | BFS depth cap per state |
| `FastChecker=1` | off | Exp-polynomial fast checker (see below) |
| `Verbose=True` | `False` | Per-state BFS progress |

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
# Check vmmc_continuous.wl (Gaussian VMMC, 3×3 seed)
wolframscript -file check.wls examples3/vmmc_continuous.wl \
  SeedBitStrings=11110101011110011 Mode=Both

# Animate on a 20×20 grid (physLen=2, cutoff just past LJ minimum)
wolframscript -file animate.wls examples3/vmmc_continuous.wl \
  Sites=400 N=10 Steps=2000 Beta=1 FPS=8 Simple=1 NoParams=1 physLen=2 '$maxD2=8'

# Check 2D Kawasaki symbolically + numerically
wolframscript -file check.wls examples3/kawasaki_2d.wl

# Check VMMC with user-defined field (fast polynomial checker)
wolframscript -file check.wls examples3/vmmc_2d_field.wl \
  MaxBitString=1111111111 Mode=Symbolic FastChecker=1

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
|------|-------------|----|----|
| `kawasaki_1d.wl` | 1D Kawasaki nearest-neighbour swap on a periodic ring | PASS | PASS |
| `kawasaki_2d.wl` | 2D Kawasaki on a periodic square lattice | PASS | PASS |
| `vmmc_2d_field.wl` | VMMC with user-defined field and coupling functions | PASS | PASS |
| `vmmc_continuous.wl` | VMMC with Gaussian proposal, LJ energy; three-file structure | PASS | PASS |
| `jump_1d_weighted.wl` | 1D jump dynamics using `RandomChoice[weights → ...]` | PASS | PASS |
| `kawasaki_1d_fail.wl` | Sign-reversed dE in Metropolis | FAIL | PASS |
| `kawasaki_1d_nonergodic.wl` | Only type-1 particles move; type-2+ frozen | PASS | FAIL |
| `cluster_1d_fail.wl` | Cluster slides only rightward (asymmetric proposal) | FAIL | FAIL |
| `vmmc_2d_edit.wl` | VMMC with only 3 of 4 directions (asymmetric proposal) | FAIL on 3×3 | — |

`kawasaki_1d_nonergodic.wl` demonstrates the key case: **DB PASS + Ergodicity FAIL**. Detailed balance cannot detect non-ergodicity; the two checks are independent.

---

## Symbolic checkers

### Standard (default)
`FullSimplify[PiecewiseExpand[expr], {β > 0, ...}]`. Correct for all algorithm types.

### FastChecker (`FastChecker=1`)
After `PiecewiseExpand`, DB expressions reduce to sums of `c·exp(−β·L)` terms. DB holds iff all coefficient groups sum to zero — verified by `Expand[...] === 0` (microseconds). Falls back to `FullSimplify` for inconclusive cases. Works directly for Metropolis acceptance; falls back gracefully for Barker/heat-bath.

Both checkers share the same efficiency pipeline: trivial-zero filter → syntactic deduplication (speedup ∝ translational symmetry) → threshold-based `ParallelMap`.

---

## Supported random calls

| Call | Behaviour |
|------|-----------|
| `RandomReal[]` / `Random[]` | Interval-tracking token |
| `RandomInteger[{lo, hi}]` | Rejection sampling over ⌈log₂(hi−lo+1)⌉ bits |
| `RandomChoice[list]` | Uniform choice via rejection sampling |
| `RandomChoice[weights → elements]` | Sequential Bernoulli decomposition |
| `RandomVariate[NormalDistribution[μ, σ]]` | Truncated Gaussian, sequential Bernoulli, nMax=⌊nGrid/2⌋ |

`RandomSample`, `RandomPermutation` are **not** supported.

---

## 2D encoding (dbc_core.wl)

The bijective integer encoding that maps bit strings to 2D square-lattice states lives at the end of `dbc_core.wl` (tagged "2D SQUARE LATTICE SPECIFIC"). For 3D or non-square geometries, replace `$decode` and update `BitsToState` in the algorithm file. The BFS core is geometry-agnostic.
