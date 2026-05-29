# Conversation Summary — 26 May 2025

This document records the ideas, clarifications, and implementation decisions from a design session around the DetailedBalanceChecker project and a broader research programme using LLM-in-loop algorithm discovery for Monte Carlo methods in polymer physics.

---

## 1. Ergodicity checking — implementation

### Conceptual basis

The existing checker uses BFS over bit sequences to discover all states reachable from a seed state. If a Markov chain is non-ergodic, the BFS will find only a subset of the states accessible in principle. For N labeled particles on S lattice sites with hard-sphere exclusion, the theoretical full state count is the falling factorial:

(S)_N = S · (S−1) · … · (S−N+1)

Both S and N are readable from any single discovered state (S = length of state array, N = count of non-zero entries), so no external parameters are needed. The check is: did BFS find exactly (S)_N states?

### Implementation

Three changes were made:

**`dbc_core.wl`** — added `CheckErgodicity[allStates_List]`:
- Reads S and N from `First[allStates]`
- Computes the falling factorial as the theoretical count
- Returns `<|"ergodic" -> Bool, "found" -> k, "theoretical" -> m|>`

**`check.wls`** — added an Ergodic column to the output:
- `colErgo = 14` added to column widths; included in `totalCols`
- `ergoPassCount`, `ergoFailCount` counters added
- `CheckErgodicity[allStates]` called after BFS, result formatted as `PASS (k)` or `FAIL (k/m)`
- Ergodic column printed in the result row and summarised at the end
- `*** ERGODICITY VIOLATION DETECTED ***` warning added to summary

**`examples3/kawasaki_1d_nonergodic.wl`** — new example file:
- Standard 1D Kawasaki algorithm modified so only type-1 particles can swap with holes
- Type-2+ particles are permanently frozen
- **Passes detailed balance** (T = 0 both ways for all cross-sector pairs; Metropolis within the type-1 sector)
- **Fails ergodicity** (from {1,2,0,0}, only 3 of 12 states are reachable)
- Provides the key test case: DB PASS + ergodicity FAIL

### Test commands

```bash
# Ergodic + DB PASS (standard Kawasaki)
wolframscript -file check.wls examples3/kawasaki_1d.wl MaxBitString=11111 Mode=Symbolic

# Non-ergodic + DB PASS  (frozen type-2 particles; the interesting case)
wolframscript -file check.wls examples3/kawasaki_1d_nonergodic.wl SeedBitStrings=11101 Mode=Symbolic

# Non-ergodic + DB FAIL  (rightward-only cluster)
wolframscript -file check.wls examples3/cluster_1d_fail.wl SeedBitStrings=11101 Mode=Symbolic
```

`SeedBitStrings=11101` decodes to ID 29 = state {1,2,0,0} on a 4-site ring (2 labeled particles, 2 holes).

---

## 2. Detailed balance vs superdetailed balance

**Detailed balance** is a condition on aggregate transition probabilities between state pairs:

T(i→j) · π(i) = T(j→i) · π(j)

It is a necessary and sufficient condition for a Markov chain to have π as its stationary distribution. It says nothing about which specific random paths contribute to T(i→j).

**Superdetailed balance** is a stronger condition on individual random paths ω (specific sequences of random choices):

P(ω) · π(i) = P(ω̄) · π(j)

where ω̄ is the time-reverse of path ω. Superdetailed balance implies detailed balance but not vice versa.

VMMC satisfies superdetailed balance via the frustrated-link construction. This is the property that the code-checker's symbolic transition matrix implicitly verifies when it confirms detailed balance — each path weight is constructed so that superdetailed balance holds path-by-path.

**For kinetics:** neither detailed balance nor superdetailed balance specifies the physical transition rates. A chain satisfying DB with Metropolis rates approximates overdamped Brownian dynamics; with Glauber rates it approximates the linear-response (Fermi-function) dynamics. Superdetailed balance constrains the relationship between forward and reverse path weights but does not determine the absolute timescale or whether collective moves correctly reproduce the Stokes drag on clusters.

---

## 3. CBMC vs VMMC — comparison

**CBMC (Configurational-Bias Monte Carlo, Siepmann & Frenkel 1992):**
- Deletes a chain segment and reattempts growth bead-by-bead
- At each bead, generates k trial positions; selects one with probability proportional to the Boltzmann weight of the new bond
- Accepts the full regrowth with ratio W_new/W_old (Rosenbluth weights)
- Ergodic by construction — any chain conformation is reachable in one move
- Kinetically meaningless — the move has no physical analog in Brownian dynamics
- Used for: equilibrium phase diagrams, chain insertion/deletion (grand canonical), dense polymer melts where single-monomer moves are inefficient

**VMMC (Virtual Move Monte Carlo, Whitelam & Geissler 2007):**
- Proposes a cluster translation; recruits neighbours probabilistically via link probabilities w_fwd = max(1 − e^{β(e_init − e_fwd)}, 0)
- Satisfies superdetailed balance via frustrated links
- Approximates overdamped Brownian dynamics — collective diffusion of bound clusters emerges naturally
- Kinetically meaningful — cluster diffusion coefficient scales approximately correctly with cluster size (though with known issues for large clusters; see §6)
- Used for: nucleation, cluster dynamics, self-assembly kinetics

**Summary:** CBMC is faster at equilibrium sampling of chain conformations; VMMC is more physically realistic for dynamics. CBMC cannot be used to study kinetics; VMMC can.

---

## 4. Physical T_phys: exact calculation from a pairwise interaction matrix

For a system of N labeled particles on an S-site lattice (d-dimensional hypercubic, L^d sites, z = 2d neighbours per site), the physical single-particle transition matrix T_phys can be computed exactly as follows.

**For each ordered state pair (i, j):**

1. States must differ by exactly one particle making one nearest-neighbour hop to a vacant site. Otherwise T_phys(i→j) = 0.

2. Compute ΔE = E(state_j) − E(state_i) using the pairwise nearest-neighbour energy:
   E = Σ_{⟨s,s'⟩} J[state[s], state[s']] with J[0, k] = 0.

3. Assign the Glauber transition rate:
   T_phys(i→j) = (1/N) · (1/z) · 1/(1 + e^{βΔE})

   The (1/N)(1/z) factor accounts for: choose one of N particles uniformly, then one of z directions uniformly. The Glauber acceptance probability 1/(1 + e^{βΔE}) is the physical rate from overdamped Langevin dynamics (Fokker-Planck in the linear response regime). It satisfies detailed balance: P_Glauber(ΔE)/P_Glauber(−ΔE) = e^{−βΔE}.

4. Diagonal entry: T_phys(i→i) = 1 − Σ_{j≠i} T_phys(i→j).

T_phys is compared to T_MC (from the checker) via KL divergence or eigenvalue spectrum to grade physical fidelity.

**Note on Metropolis vs Glauber:** Metropolis acceptance min(1, e^{−βΔE}) also satisfies detailed balance and can be used as an alternative reference. Glauber is the correct physical rate for a particle coupled to a thermal bath via a Markovian process (it arises from solving the Fokker-Planck equation for overdamped motion in a symmetric double-well potential).

---

## 5. Stokesian dynamics in T_phys — why exact inclusion is impossible

The natural next question is whether the configuration-dependent Stokes drag on particles (hydrodynamic interactions, HI) can be incorporated into T_phys exactly. The answer is no, for three independent fundamental reasons:

**Reason 1 — The coarse-graining problem.** HI are governed by the Smoluchowski equation in continuous space. Mapping it to a discrete Markov chain requires integrating over transition-state surfaces in (3N−1)-dimensional configuration space — analytically intractable beyond N = 1. Kramers theory is the standard approximation, but it is exact only for high barriers, which do not exist for nearest-neighbour lattice hops.

**Reason 2 — Off-diagonal mobility requires multi-particle moves.** The mobility tensor M_{ij} couples the velocity of particle i to forces on particle j (solvent-mediated drag). On a lattice, this can only be represented as simultaneous two-particle hops, with rates proportional to off-diagonal M_{ij}. Including these grows as N² and makes it impossible to satisfy detailed balance cleanly when single-particle and two-particle moves both connect the same state pair.

**Reason 3 — No solvent degrees of freedom.** HI is mediated by the solvent. Any M_{ij} assigned to a lattice model is imported from continuous-space Stokesian dynamics theory, not derived from the model. The only exact approaches are those that model the solvent explicitly (Lattice Boltzmann, MPCD).

**Best practical approximation:** Use the Rotne-Prager-Yamakawa (RPY) tensor to compute configuration-dependent self-mobility D_i(config), and modify single-particle hop rates by D_i(config)/D_0. This is exact in the dilute limit (leading order in a/r) and captures the dominant physical effect (particles slow down when close together). Cross-mobility corrections require two-particle moves and are only approximate.

---

## 6. Simulation tools for physically accurate polymer dynamics

For a system where physically accurate dynamics are required (correct cluster diffusion, reptation, Zimm dynamics, dynamic structure factor):

| Method | HI | Reptation | Timescale | Tool |
|---|---|---|---|---|
| All-atom MD | Yes (explicit solvent) | Too slow for long chains | fs–ns | GROMACS, LAMMPS |
| Kremer-Grest CG MD | No (Langevin thermostat) | Yes, N >> N_e | µs accessible | LAMMPS, GROMACS |
| Brownian dynamics + RPY | Approximate | No (single chain) | µs | HOOMD-blue |
| DPD | Long-wavelength HI | Approximately | Fast | LAMMPS, ESPResSo |
| **MPCD + MD** | **Exact (FDT)** | **Yes** | **Mesoscale** | **ESPResSo** |
| Lattice Boltzmann + MD | Long-wavelength | Approximately | Large scale | ESPResSo, LAMMPS |

**Recommendation:** ESPResSo with MPCD for the best combination of physical accuracy and tractability. MPCD (Multi-Particle Collision Dynamics / Stochastic Rotation Dynamics) exactly satisfies the fluctuation-dissipation theorem, correctly reproduces Zimm dynamics (D ∝ N^{−0.6}) and cluster diffusion (D ∝ 1/R_H) without parameter fitting.

For reptation specifically: HI are screened in dense melts, so Kremer-Grest with a Langevin thermostat (no explicit solvent) correctly captures reptation (D ∝ N^{−2} for N >> N_e ≈ 85 beads in the standard model).

---

## 7. MPCD+MD as a ground truth for grading — assessment

**The idea:** Run MPCD+MD to get physically accurate dynamics, extract observables, use them to grade MC algorithms produced by the LLM+checker loop.

**What is correct:**
- The lattice→continuum limit is correct for statics and for long-time/large-scale dynamics (universality class is preserved)
- MPCD+MD is the appropriate physical reference for colloidal systems
- The grading framework has clear precedent (coarse-grained force field validation)

**Flaw 1 — Parameter mapping is non-unique.** The lattice coupling constants J_ab do not map uniquely to a continuous pair potential V(r) plus a timescale D_0 = k_BT/(6πηa). Absolute diffusion coefficients are not directly comparable.

**Flaw 2 — Irreducible lattice artifacts.** Rotational symmetry breaking, hard-core vs soft-core repulsion, and absence of intra-cell relaxation create systematic deviations that no MC algorithm on the lattice can eliminate. These cannot be attributed to algorithmic failure.

**Flaw 3 — "Better than VMMC" must be defined carefully.** VMMC fails cluster-size scaling (D_cluster is approximately independent of cluster size rather than D_cluster ∝ n^{−1/3}). But this failure is shared by any single-particle-based lattice algorithm without Stokes-corrected rates — it is not an algorithmic flaw unique to VMMC.

**The fix — compare dimensionless exponents, not absolute values:**
- Cluster size scaling: Does D_cluster ∝ n^{−ν} with the right ν?  
  Physical (Stokes): ν ≈ 1/3 for compact clusters. VMMC: ν ≈ 0.
- MSD scaling exponent α: ⟨Δr²⟩ ∝ t^α (Zimm: α ≈ 0.54, Rouse: α = 0.5, free: α = 1)
- Shape of the intermediate scattering function F(q,t)

These dimensionless observables characterise the dynamical universality class and are independent of the timescale mapping between MC steps and physical time. Comparing them between MPCD+MD and the MC algorithm is scientifically defensible and unambiguous.

---

## 8. Can the code-checker be adapted to off-lattice systems?

**Short answer: No, the checker cannot be directly adapted to off-lattice systems.** The BFS-over-bits approach requires:
1. A finite discrete state space (BFS terminates only if there are finitely many reachable states)
2. `RandomReal[]` appearing only in comparisons (interval tracking works only for branching, not for state generation)
3. The same state reachable by multiple paths (so transition weights can accumulate)

All three fail for continuous-space systems: the state space is uncountably infinite, displacements are generated by `RandomReal[]`, and two random paths almost surely lead to different points in R^{dN}.

The original VMMC was verified analytically (superdetailed balance proof by Whitelam & Geissler), not by a code checker. Any off-lattice algorithm would require the same analytic approach.

**The correct role of the checker in this programme:** verify the lattice version exactly as a necessary condition. If the algorithm fails on the lattice, it will fail off-lattice too. Passing on the lattice, combined with an analytic proof that the structure is preserved in the continuous generalisation, constitutes a complete verification.

---

## 9. Unified code for lattice checking and continuous-space physical simulation

**The answer is yes**, with one specific design constraint.

**The constraint:** All randomness in the move proposal must use discrete random choices (`RandomChoice[direction_list]` with a fixed step size b), not `RandomReal[]` for position generation. This is natural for VMMC-style algorithms, which already use fixed-magnitude cluster translations.

**The architecture:** Switch from a site-indexed state (index = lattice site, value = particle type) to a particle-indexed state (index = particle, value = {position, type}):

```
Current:  {0, 1, 2, 0}                          (* site-indexed *)
Unified:  {{0.0, 0.0, 1}, {1.0, 0.0, 2}}        (* particle-indexed, real positions *)
```

| Component | Change required |
|---|---|
| `BitsToState` | Map to grid-point configurations (positions = n × b); checker mode |
| `energy` | Continuous pair potential V(‖r_i − r_j‖) instead of J[a,b] at adjacent sites |
| Neighbour test | Distance cutoff instead of site adjacency |
| Move proposal | `RandomChoice[directions] × b` — same in both modes |
| Acceptance | Unchanged — `RandomReal[]` in comparison only |
| BFS core | Unchanged |
| Detailed balance algebra | Unchanged |

**Checker mode:** `BitsToState` maps to positions on a regular grid with spacing b. Since all moves displace by exactly b, positions remain on the grid. State space is finite. BFS, interval tracking, and symbolic detailed balance check all apply without modification.

**Simulation mode:** Initialise at arbitrary real positions from a physical distribution. Run with a physically motivated pair potential (WCA, Yukawa, etc.) and any b. Compare to MPCD+MD for dynamical validation.

**The same cluster-building code, link probability formulae, and acceptance criterion execute in both modes.** Only the state initialisation, the energy function, and whether `BitsToState` is called differ.

**Limitation:** Fixed step size means short-time dynamics is anisotropic and step-like. Long-time dynamical exponents (which are universal) are correctly reproduced. The comparison to MPCD+MD should focus on dimensionless exponents and scaling, not absolute short-time rates.

---

## 10. Research programme summary

The full programme as clarified in this conversation:

1. **LLM proposes** a new cluster-move algorithm in a unified particle-indexed language that works in both lattice and continuous modes.

2. **Checker verifies** detailed balance and ergodicity exactly on the grid-restricted version (lattice checker, existing infrastructure with ergodicity extension added today).

3. **Physical simulation** runs the same code with a physically motivated pair potential (WCA etc.) and arbitrary initial positions — no checker involved.

4. **MPCD+MD comparison** grades the algorithm on dimensionless dynamical observables: cluster diffusion exponent ν (target: ν ≈ 1/3 for Stokes), MSD scaling α, intermediate scattering function shape. This is the fidelity signal.

5. **LLM iterates** guided by the combined formal verdict (pass/fail) and fidelity score (continuous).

6. **Analytic generalisation** — once a lattice algorithm passes both checks, the continuous-space version is written and verified analytically following the Whitelam-Geissler superdetailed balance argument. The lattice checker provides the computationally verified reference case; the analytic proof covers the continuous generalisation.

**The novel contribution relative to existing literature** is the combination of:
- Exact formal verification (detailed balance + ergodicity) via the symbolic checker
- Physical fidelity grading against MPCD+MD using dimensionless dynamical exponents
- Automated algorithm discovery via LLM in an executable loop
- The fourth mechanism (symbolic algebraic witnesses from Mathematica) on top of the three standard execution-loop benefits (grounding, offloading, error localisation)

No existing literature (FunSearch, AlphaEvolve, Clover, LLMLOOP) has applied this combination to physics algorithm design.
