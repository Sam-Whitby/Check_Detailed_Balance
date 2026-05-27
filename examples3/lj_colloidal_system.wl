(* ================================================================
   lj_colloidal_system.wl — Lennard-Jones colloidal system physics
   ================================================================
   Defines the PHYSICAL SYSTEM: what the energy landscape looks like.
   Load from an algorithm file after vmmc_2d_grid.wl:
     Get[DirectoryName[$InputFileName] <> "lj_colloidal_system.wl"]

   This file is independent of the kinetics (Algorithm[]).  It is also
   the correct place to define the energy used by any T_phys ground-truth
   comparison (Glauber rates etc.), since both the MC algorithm and the
   physical reference share the same energy landscape.

   PARAMETERS (edit this section to change the model):
     physLen — particle diameter in lattice units.
               sigLJ = physLen: LJ zero-crossing at d² = physLen².
               LJ minimum at d² = 2^(1/3)·physLen² ≈ 1.26·physLen².
     epsLJ   — LJ well depth in kT units.
     numBeta — inverse temperature for numerical MCMC.

   INTERACTION CUTOFF:
     $maxD2 = Infinity  is correct for the symbolic checker (all
       pairs are included; small lattices have finitely many anyway).
     For production runs set $maxD2 to include the desired range:
       Ceiling[2*physLen^2]    — cutoff just past the LJ minimum
       Ceiling[6.25*physLen^2] — standard 2.5σ cutoff (full tail)

   physLen/nGrid SCALING for physical runs:
     For N particles at packing fraction η on an nGrid×nGrid grid:
       physLen = 2 * nGrid * Sqrt[eta / (N * Pi)]
     Example: η=0.2, N=10, nGrid=20 → physLen ≈ 2.3
     Set $maxD2 = Ceiling[2 * physLen^2] for that run.

   CHECKER INTERFACE:
     $checkerAbstractParams = {"physLen","epsLJ","sigStep"} — these
     are cleared to free symbols before symbolic BFS (so the proof
     holds for all parameter values simultaneously) and restored as
     concrete numeric values inside a Block for numerical MCMC.
     String names are required; symbol references evaluate too early.
   ================================================================ *)


(* ---- Physical parameters (concrete values for numerical runs) ---- *)

physLen = 1    (* particle diameter in lattice units *)
epsLJ   = 1    (* LJ well depth in kT units *)
numBeta = 1    (* inverse temperature *)

sigLJ   := physLen                                   (* delayed: always = physLen *)
sigStep := physLen * Sqrt[2.0 / (numBeta * epsLJ)]  (* BD step at natural LJ timescale *)

$maxD2 = Infinity


(* ---- Checker interface ---- *)

$checkerAbstractParams = {"physLen", "epsLJ", "sigStep"}

(* ---- Concrete pair energy (numerical MCMC and animation only) ----
   $pairEnergy is activated inside a Block by check.wls / animate.wls.
   During symbolic BFS, couplingJ has no DownValues — each call
   couplingJ[a,b,d2] is a free real atom, proving DB for all couplings. *)
$pairEnergy[a_Integer, b_Integer, d2_Integer] :=
  If[d2 == 0 || d2 === Infinity || d2 > Ceiling[2 * physLen^2], 0,
     4 * epsLJ * ((sigLJ^2/d2)^6 - (sigLJ^2/d2)^3)]

$couplingFormulaStr = "4*epsLJ*((physLen^2/d2)^6-(physLen^2/d2)^3)"
$concreteParams     = <||>


(* ================================================================
   Energy function
   ================================================================
   Total pair energy summed over occupied pairs within $maxD2.
   Uses couplingJ (abstract during BFS, concrete during numerical MCMC).
   Min/Max canonicalises argument order so the same symbolic atom is
   used regardless of bond orientation. *)

energy[state_List] :=
  Module[{nGrid = Round[Sqrt[Length[state]]], occ},
    occ = Flatten[Position[state, _?(# > 0 &)]];
    Total @ Flatten @ Table[
      With[{d2 = $torusD2[occ[[i]], occ[[j]], nGrid]},
        If[d2 <= $maxD2,
           couplingJ[Min[state[[occ[[i]]]], state[[occ[[j]]]]],
                     Max[state[[occ[[i]]]], state[[occ[[j]]]]],
                     d2],
           0]],
      {i, Length[occ]}, {j, i+1, Length[occ]}]]


(* ================================================================
   Virtual pair energy for VMMC cluster builder
   ================================================================
   Energy between particle typeI at virtual site vI and particle
   typeJ at site qSite.  Returns Infinity on hard-core overlap,
   0 outside the cutoff. *)

$virtualPairEnergy[typeI_, typeJ_, vI_, qSite_, nGrid_] :=
  Which[
    vI === qSite, Infinity,
    True,
      With[{d2 = $torusD2[vI, qSite, nGrid]},
        If[d2 > 0 && d2 <= $maxD2,
           couplingJ[Min[typeI, typeJ], Max[typeI, typeJ], d2], 0]]]
