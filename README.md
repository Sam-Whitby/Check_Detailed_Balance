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

   If the algorithm declares `$symmetryGroup`, BFS operates on *canonical representatives* only — one state per symmetry orbit. States are canonicalised to the lexicographic minimum over the orbit before being added to the BFS queue. This reduces the number of states processed by a factor equal to the orbit size (up to 8·nGrid² for full D4 + translation symmetry), without affecting correctness.

3. **Symbolic DB check.** For each canonical state pair (i,j): `|orbit(i)|·T(i→j)·exp(−β E(i)) = |orbit(j)|·T(j→i)·exp(−β E(j))`. The orbit-size factors (integers) correct for the fact that canonical states may represent different numbers of physical states. β is a free symbol so Boltzmann factors cancel algebraically.

4. **Numerical MCMC check.** Canonical states are expanded back to their full orbits for MCMC; visit counts and Boltzmann weights are aggregated by orbit. KL < 0.02 is a pass.

**Ergodicity check.** For N labeled particles on S sites, the theoretical reachable state count is (S)_N = S·(S−1)·…·(S−N+1). With symmetry enabled, this full count is compared against the number of physical states represented by the canonical set.

---

## Symmetry reduction

Algorithm files can declare `$symmetryGroup` to enable BFS state canonicalization:

```mathematica
$symmetryGroup = {"translation"}        (* torus translations only *)
$symmetryGroup = {"translation", "D4"}  (* translations + 4 rotations + 4 reflections *)
```

**What counts as a valid symmetry group:**
- The energy function must be invariant under the declared symmetries.
- The proposal distribution must be invariant (i.e. the move kernel commutes with the symmetry).
- Both conditions hold for nearest-neighbour or isotropic energy functions on a square torus with uniform or Gaussian proposals.

**Field check.** If the algorithm defines `fieldF` or `$fieldFConcrete` (a spatially varying external field), symmetry reduction is automatically disabled regardless of `$symmetryGroup`. Fields break translational symmetry.

**Orbit-size correction.** When two canonical states have orbits of different sizes (e.g. a fully symmetric configuration vs. a generic one), the DB condition must account for this. The checker scales the transition matrix entry T(i→j) by |orbit(i)| before the algebraic check, so integer factors appear instead of rational functions of β — much easier for `FullSimplify`.

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
$symmetryGroup = {"translation", "D4"} (* symmetry group for BFS canonicalization *)
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
| `PhysFidelity=1` | off | Compare T_MC to physical Glauber dynamics; adds M1 and M2 columns (see below) |
| `Verbose=True` | `False` | Per-state BFS progress |

`NGrid` and `MaxComponents` together give an intuitive interface: `NGrid=2 MaxComponents=6` means "check the first 6 distinct 2×2 systems, stop there." Both options use lazy ID iteration internally so no large list is held in memory, even when the total state count for that grid size is in the millions.

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
# Check all 2×2 components of vmmc_continuous.wl (symmetry-reduced, lazy iteration)
wolframscript -file check.wls examples3/vmmc_continuous.wl NGrid=2 Mode=Both

# Check the first 3×3 component only (lazily — no large list built upfront)
wolframscript -file check.wls examples3/vmmc_continuous.wl \
  NGrid=3 MaxComponents=1 Mode=Symbolic FastChecker=1

# Check a specific seed state (BFS discovers the full component)
wolframscript -file check.wls examples3/vmmc_continuous.wl \
  SeedBitStrings=11110101011110011 Mode=Both

# Check 2D Kawasaki (all 2×2 components, symmetry-reduced)
wolframscript -file check.wls examples3/kawasaki_2d.wl NGrid=2 Mode=Symbolic

# Check VMMC with user-defined field (symmetry auto-disabled; fast polynomial checker)
wolframscript -file check.wls examples3/vmmc_2d_field.wl \
  MaxBitString=1111111111 Mode=Symbolic FastChecker=1

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

| File | Description | Symmetry | DB | Ergodic |
|------|-------------|----------|----|---------|
| `kawasaki_1d.wl` | 1D Kawasaki nearest-neighbour swap on a periodic ring | — | PASS | PASS |
| `kawasaki_2d.wl` | 2D Kawasaki on a periodic square lattice | translation, D4 | PASS | PASS |
| `vmmc_2d_field.wl` | VMMC with user-defined field and coupling functions | — (field) | PASS | PASS |
| `vmmc_continuous.wl` | VMMC with Gaussian proposal, LJ energy; three-file structure | translation, D4 | PASS | PASS |
| `jump_1d_weighted.wl` | 1D jump dynamics using `RandomChoice[weights → ...]` | — | PASS | PASS |
| `kawasaki_1d_fail.wl` | Sign-reversed dE in Metropolis | — | FAIL | PASS |
| `kawasaki_1d_nonergodic.wl` | Only type-1 particles move; type-2+ frozen | — | PASS | FAIL |
| `cluster_1d_fail.wl` | Cluster slides only rightward (asymmetric proposal) | — | FAIL | FAIL |
| `vmmc_2d_edit.wl` | VMMC with only 3 of 4 directions (asymmetric proposal) | — | FAIL on 3×3 | — |

`kawasaki_1d_nonergodic.wl` demonstrates the key case: **DB PASS + Ergodicity FAIL**. Detailed balance cannot detect non-ergodicity; the two checks are independent.

`vmmc_2d_field.wl` demonstrates automatic symmetry disabling: it declares `$symmetryGroup = {"translation", "D4"}` but the presence of `fieldF` causes the checker to ignore this declaration.

---

## Symbolic checkers

### Standard (default)
`FullSimplify[PiecewiseExpand[expr], {β > 0, ...}]`. Correct for all algorithm types.

### FastChecker (`FastChecker=1`)
After `PiecewiseExpand`, DB expressions reduce to sums of `c·exp(−β·L)` terms. DB holds iff all coefficient groups sum to zero — verified by `Expand[...] === 0` (microseconds). Falls back to `FullSimplify` for inconclusive cases. Works directly for Metropolis acceptance; falls back gracefully for Barker/heat-bath.

Both checkers share the same efficiency pipeline: trivial-zero filter → syntactic deduplication (speedup ∝ translational symmetry) → threshold-based `ParallelMap`.

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

`RandomSample`, `RandomPermutation` are **not** supported.

---

## 2D encoding (dbc_core.wl)

The bijective integer encoding that maps bit strings to 2D square-lattice states lives at the end of `dbc_core.wl` (tagged "2D SQUARE LATTICE SPECIFIC"). For 3D or non-square geometries, replace `$decode` and update `BitsToState` in the algorithm file. The BFS core is geometry-agnostic.
