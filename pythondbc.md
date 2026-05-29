# Python Detailed Balance Checker (SZPure)

## Goal

A Python function `check_detailed_balance(algorithm, energy, seed_state, ...)` that verifies detailed balance probabilistically (Schwartz-Zippel style) for any MCMC algorithm written in standard Python. No SymPy required.

---

## User API

```python
import random

# User writes algorithm normally, calling random.random() wherever needed
def vmmc_step(state, beta, params):
    candidate = propose(state)
    delta_E = energy(candidate, params) - energy(state, params)
    if random.random() < min(1.0, math.exp(-beta * delta_E)):
        return candidate
    return state

def energy(state, params):
    return ...   # float

# Run the checker
result = check_detailed_balance(
    algorithm  = lambda state: vmmc_step(state, beta, params),
    energy     = lambda state: energy(state, params),
    seed_state = initial_state,
    rng_target = 'random.random',   # dotted path to intercept
    beta_range = (0.1, 10.0),       # SZ sampling range for beta
    k          = 100,               # SZ evaluation points
)
```

The user provides a **closed-over** `algorithm(state) -> state` that calls `random.random()` freely. No other modification required.

---

## Core Mechanism: RNG Interception

`unittest.mock.patch` replaces the named RNG function for the duration of the call:

```python
from unittest.mock import patch

def run_with_script(algorithm, state, script):
    """Run algorithm with a pre-determined sequence of RNG values."""
    it = iter(script)
    with patch(rng_target, side_effect=lambda: next(it)):
        return algorithm(state)
```

This intercepts every call regardless of call depth. The algorithm is unmodified.

---

## Execution-Path BFS (per state, per parameter point)

For a given state `s` and fixed numerical parameters, enumerate all distinct execution paths and their transition probabilities.

Each path is characterised by the sequence of RNG values drawn and the **threshold** at each draw — the value at which the algorithm branches to a different outcome. Paths are the leaves of a binary tree.

**Threshold discovery at each branch point:**

```python
def find_threshold(algorithm, state, prefix):
    """
    Returns threshold t such that:
      run(script = prefix + [u < t]) -> outcome_A
      run(script = prefix + [u >= t]) -> outcome_B
    Returns None if this RNG call is a no-op (outcome identical either side).
    """
    lo_outcome = run_with_script(algorithm, state, prefix + [0.0])
    hi_outcome = run_with_script(algorithm, state, prefix + [1.0])
    if lo_outcome == hi_outcome:
        return None, lo_outcome   # call doesn't branch; outcome determined already
    lo, hi = 0.0, 1.0
    while hi - lo > 1e-14:
        mid = (lo + hi) / 2
        if run_with_script(algorithm, state, prefix + [mid]) == lo_outcome:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2, lo_outcome, hi_outcome
```

**BFS over paths:**

```python
def build_transition_row(algorithm, state):
    """Returns {next_state: probability} for all reachable outcomes from state."""
    T = defaultdict(float)
    queue = [( [], 1.0 )]   # (script_prefix, probability_of_reaching_here)
    while queue:
        prefix, prob = queue.pop()
        t, *outcomes = find_threshold(algorithm, state, prefix)
        if t is None:
            # No further branching: this path terminates
            next_state = outcomes[0]
            T[next_state] += prob
        else:
            lo_outcome, hi_outcome = outcomes
            # Low branch: u drawn uniformly from [0, t)  -> probability = t
            # High branch: u drawn uniformly from [t, 1) -> probability = 1-t
            if isinstance(lo_outcome, list):   # path continues (more RNG calls)
                queue.append((prefix + [t * 0.5],     prob * t))
                queue.append((prefix + [(1 + t) * 0.5], prob * (1 - t)))
            else:
                T[lo_outcome] += prob * t
                T[hi_outcome] += prob * (1 - t)
    return T
```

*Note:* `find_threshold` re-runs the algorithm O(log(1/ε)) times per branch point. For algorithms with k branch points the total runs per state are O(k · 2^k · log(1/ε)). For typical VMMC steps (k ≤ 4) this is O(64) runs — negligible.

---

## State-Space BFS

Enumerate the connected component reachable from `seed_state`:

```python
def enumerate_component(algorithm, seed_state):
    rows = {}
    queue = [seed_state]
    visited = {seed_state}
    while queue:
        s = queue.pop()
        row = build_transition_row(algorithm, s)
        rows[s] = row
        for ns in row:
            if ns not in visited:
                visited.add(ns)
                queue.append(ns)
    return rows   # full transition matrix as nested dict
```

---

## SZPure Detailed Balance Check

Repeat `k` times with independently sampled parameter points:

```python
def check_detailed_balance(algorithm, energy, seed_state,
                           rng_target='random.random',
                           beta_range=(0.5, 5.0), k=100):
    violations = []
    for _ in range(k):
        beta = random.uniform(*beta_range)
        # Rebuild transition matrix at this beta (algorithm is closed over beta,
        # so re-running it with a new beta closure is sufficient if the user
        # parameterises correctly — see note below)
        T = enumerate_component(algorithm, seed_state)
        
        # Check detailed balance for every directed pair
        for s, row in T.items():
            for t, p_fwd in row.items():
                p_rev = T.get(t, {}).get(s, 0.0)
                pi_s = math.exp(-beta * energy(s))
                pi_t = math.exp(-beta * energy(t))
                lhs = p_fwd * pi_s
                rhs = p_rev * pi_t
                if abs(lhs - rhs) > tol * max(abs(lhs), abs(rhs), 1e-30):
                    violations.append((s, t, beta, lhs, rhs))
    return violations
```

**Note on beta as a free parameter:** If the algorithm closes over `beta` as a Python variable, the caller must re-create the algorithm closure for each sampled `beta`. A clean convention is to pass `algorithm_factory(beta) -> callable`, analogous to how `DynamicSymParams` works in the Mathematica code.

---

## Limitations vs. Mathematica SZPure

| Aspect | Python (this plan) | Mathematica SZPure |
|--------|--------------------|--------------------|
| Interception | `patch` — identical semantics | `Block` — identical semantics |
| Threshold discovery | Binary search, O(k·2^k·log ε⁻¹) algorithm runs | Symbolic; thresholds exact by construction |
| Floating-point thresholds | ε~1e-14 precision; degenerate thresholds possible | Exact symbolic split |
| Symbolic parameters | Must re-run per beta sample | Single BFS, k evaluations |
| Variable-length RNG usage | Needs depth/call-count bound | Handled by BFS termination |
| Ergodicity check | Free: check BFS covers all expected states | Same |

The main practical cost is re-running `enumerate_component` for each of the k beta samples, rather than building one symbolic matrix and evaluating it k times. For small state spaces (≲1000 states) this is acceptable; for large ones, caching the tree structure and re-evaluating thresholds at new beta values would recover most of the efficiency.
