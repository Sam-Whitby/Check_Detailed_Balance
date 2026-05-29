# Algebraic Verification of Translation Invariance via Symbolic τ

## The idea

Rather than testing translation invariance by running the BFS from multiple states and comparing results, the checker can verify it algebraically from a **single BFS run** by introducing a symbolic translation offset τ into the state provided to the algorithm.

The checker provides particle positions as `canonical_pos + τ` where τ is an unevaluated Mathematica symbol (a 2D offset `{τr, τc}`). The algorithm runs normally on this state. For any physically meaningful computation:

```
(canonical_i + τ) - (canonical_j + τ)  =  canonical_i - canonical_j
```

τ cancels algebraically in every pairwise quantity — distances, displacements, energy differences — before any nonlinear function (Exp, Min, Mod) is applied. The BFS tree's leaf probabilities are symbolic expressions in `{β, ε, J, ...}` exactly as they are now. If the algorithm is translation-invariant, τ never appears in any of them.

After the BFS, a single check suffices:

```mathematica
allFree = AllTrue[Values[matrix], FreeQ[#, τr] && FreeQ[#, τc] &]
```

If `allFree` is True: translation invariance is proven algebraically for every microstate reachable in the BFS. If False: the specific transitions and expressions containing τ identify exactly where and how the algorithm breaks translation invariance.

## Why this is stronger than testing translated copies

- **Coverage**: a single BFS from the seed state explores the entire ergodic component. Every reachable microstate is tested simultaneously — not just specific translated copies of the seed.
- **Algebraic certainty**: τ is a free symbol, not a sampled value. The check is not probabilistic. If τ does not appear, it cannot appear for any value of τ, including all nGrid² lattice translations.
- **Zero extra BFS runs**: compared to the current approach (which trusts the user's declaration), the overhead is one FreeQ pass over the existing leaf expressions — effectively free.

## What is required

The algorithm must receive particle positions in a form that allows τ to propagate through arithmetic. Two specific requirements:

**1. Positions must be expressible as symbolic Mathematica quantities.**
The checker constructs the state with positions `canonical_pos + {τr, τc}` instead of pure integers. The algorithm's energy and acceptance computations receive these and perform arithmetic on them. For pairwise interactions, τ cancels automatically.

**2. Periodic boundary conditions must wrap the displacement, not the absolute coordinate.**
The correct PBC formula computes the displacement between two particles and then wraps it:
```mathematica
displacement = Mod[pos_i - pos_j, nGrid]   (* τ cancels inside, then wrap — τ-free *)
```
The incorrect formula wraps absolute coordinates first:
```mathematica
Mod[pos_i, nGrid] - Mod[pos_j, nGrid]      (* τ does not cancel — mathematically wrong anyway *)
```
This is not an extra constraint on code style: wrapping the displacement is the physically correct formula for minimum-image convention. Wrapping absolute coordinates before differencing is incorrect physics (as `Mod[a+τ,n] - Mod[b+τ,n] ≠ a - b` in general — verified numerically).

## Failure modes and what they mean

**τ appears in a probability expression**: The algorithm uses absolute position in its energy or acceptance calculation (e.g., an external field). Correctly detected and reported.

**BFS fails to evaluate**: The algorithm contains a conditional branch on an absolute coordinate (`If[pos > boundary, ...]`). With symbolic `pos = canonical + τ`, Mathematica leaves the conditional unevaluated and the BFS stalls. This is also correctly interpreted as position dependence — a translation-invariant algorithm never branches on absolute position.

**Ergodicity**: As with the existing checker, the verification covers only the ergodic component reachable from the seed state. If some states are not reachable, they are not covered. This is unchanged from the current situation.

## Checking for τ in expressions: FreeQ

The check for τ in a Mathematica expression uses `FreeQ`:

```mathematica
FreeQ[expr, τr] && FreeQ[expr, τc]
```

`FreeQ` is a **pure structural pattern-match** on the expression tree. It walks the tree and returns False the moment it encounters τ, or True if it completes the walk without finding it. It performs no algebraic manipulation whatsoever — no simplification, no substitution, no rewriting.

Measured timings on a realistic sum of 80 symbolic exp/cos terms:

| Operation | Time |
|-----------|------|
| `FreeQ[expr, τr]` | ~1–2 microseconds |
| `FullSimplify[expr]` | >10 seconds (did not complete) |

The difference is approximately **7 orders of magnitude**. `FreeQ` is O(n) in the number of nodes in the expression tree — it is one of the cheapest operations in Mathematica.

One subtlety: `FreeQ` is syntactic, not semantic. If τ appears multiplied by zero (`τr * 0`), Mathematica automatically evaluates this to `0` at expression-construction time, so τ never enters the tree and `FreeQ` correctly returns True. Pathological cases where τ survives structurally with a non-trivially-zero coefficient (e.g. `τr * Sin[π * n]` for integer n) do not arise in practice from the arithmetic of VMMC energy computations.

## Implementation difficulty

**Easy parts:**
- The FreeQ check is a single line added after the BFS, essentially free computationally.
- The BFS machinery is unchanged — τ is just another free symbol alongside β, ε, J.
- Adding `{τr, τc}` to canonical positions before constructing the seed state is trivial.

**The hard part — state representation:**

The current state representation is a flat integer vector (`{0, 1, 0, 2, ...}`) where `state[[site_index]]` gives the particle type at that absolute site. This representation is incompatible with symbolic positions, for two reasons:

1. `state[[canonical_site + τ]]` — Mathematica cannot index a list at a symbolic position. Any place the algorithm looks up the occupancy of a derived site (e.g. "the neighbor of site i") would fail.
2. If the algorithm iterates `Do[f[state[[i]]], {i, L}]`, the loop variable i remains an integer, but any computed neighbor `i + 1` or `Mod[i + nGrid, L] + 1` would contain τ if i was derived from a symbolic position.

The resolution requires changing the state representation to a list of `{position, type}` pairs, where positions carry the τ offset. Spatial queries (neighbor lookup, distance computation) then operate on differences between positions, which cancel τ naturally. This is a moderate rewrite of the state interface — not a drop-in change for any algorithm currently written to accept a flat occupancy vector.

**D4 generator testing — easy by comparison:**

The D4 generator check requires no change to the state representation. It uses the existing integer-indexed states, applies precomputed permutations (already in `$dbcAllGroupPerms`), runs BFS from those permuted states, and compares the resulting transition matrices using the existing SZPure machinery. This can be added to the current codebase with minimal changes.

## Implementation in the existing checker

The change is localised to the TrustSymmetry path in `check.wls` and the state construction in `dbc_core.wl`:

1. **State construction**: when building the seed state to pass to the algorithm, add `{τr, τc}` to every particle's position. The state encoding needs to be a list of `{position, type}` pairs (rather than a flat occupancy vector indexed by absolute site) so that positions containing τ can be represented.

2. **BFS**: unchanged. The BFS intercepts random numbers exactly as now. τ is just another free symbol in the resulting expressions, handled automatically by Mathematica's symbolic engine.

3. **Post-BFS check**: after `$dbcExpandOrbitsToMatrix`, scan all matrix values for τ:
```mathematica
τFree = AllTrue[Values[matrix], FreeQ[#, τr] && FreeQ[#, τc] &];
If[τFree,
   (* Translation invariance algebraically verified — proceed with TrustSymmetry *),
   (* Report which transitions contain τ and their expressions *)
]
```

4. **If verified**: use TrustSymmetry for translations with full algebraic justification. No user declaration needed.

5. **If not verified**: fall back to BFS without translational TrustSymmetry, and report the specific violation.

---

## D4 (rotational) invariance: why the same trick does not apply

The τ approach works because translation is **additive**: τ appears linearly in every position, and linear terms cancel in differences. No further structure is needed.

D4 is a **discrete non-additive group**. Its elements (four rotations and four reflections) cannot be parameterised by a single free symbol that cancels algebraically. Specifically:

- A rotation R acts as `R(pos_i) - R(pos_j) = R(pos_i - pos_j)` (linearity holds), and `|R(d)|² = |d|²` (distance preserved). But to have Mathematica verify this for a symbolic group element ρ ∈ {1,...,8}, it would need to know that all eight D4 matrices are orthogonal — a structural fact that cannot be expressed as a single free symbol cancellation.
- D4 is finite with only 8 elements and 2 generators. There is no "generic" D4 element analogous to a continuous translation offset.

The correct approach for D4 is **explicit generator testing**:

1. Run BFS from each canonical orbit representative (already done for TrustSymmetry).
2. For each orbit representative s₀: also run BFS from r(s₀) (90° rotation) and s(s₀) (one reflection).
3. Compare: does T(r(s₀) → r(t)) = T(s₀ → t) for all transitions t? Same for s.
4. If consistent for both generators across all orbit representatives: D4 invariance is exactly verified.

**Why 2 generators are sufficient (group theory)**

D4 has 8 elements but is generated by r and s alone. Every element can be written as a composition: {e, r, r², r³, s, rs, r²s, r³s}. If T(r·i → r·j) = T(i→j) and T(s·i → s·j) = T(i→j) hold for all states i,j, then for any g = g₁g₂...gₖ expressed as a word in the generators:

```
T(g₁g₂...gₖ · i → g₁g₂...gₖ · j)
  = T(g₂...gₖ · i → g₂...gₖ · j)   [generator invariance at g₁]
  = ...
  = T(i → j)
```

So r², r³, rs, r²s, r³s never need to be tested individually — they follow by induction. This requires the generator invariance to hold **for all states i**, which is why BFS must be run from every orbit representative, not just one.

**Cost**: 2 extra BFS runs per orbit representative. For a component with 8 orbit representatives: 16 extra BFS runs total, compared to 504 for full enumeration without TrustSymmetry. The 8× D4 speedup is earned with ~3% overhead and the verification is exact.

Unlike the translation check (which is free via FreeQ), D4 requires these extra runs because there is no symbolic shortcut — the cancellation trick depends on τ being additive, which rotation is not.
