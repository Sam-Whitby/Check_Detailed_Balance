(* ================================================================
   vmmc_biased.wl — Intentionally broken VMMC: rightward bias
   ================================================================

   PURPOSE
   -------
   This file is a deliberately incorrect VMMC algorithm, kept for
   testing purposes.  It is structurally identical to vmmc_lattice.wl
   except for one change: displacements with dx > 0 (rightward) are
   duplicated in the proposal list, making them twice as probable as
   their negatives.  This breaks proposal symmetry P(d) = P(−d) and
   therefore violates detailed balance.

   THE BREAK
   ---------
   For any pair of states (s, s') where s' = s with a cluster shifted
   right by some dx > 0:

     T(s → s')  ∝  (2 / |$displacements|) × (build prob) × (accept)
     T(s'→ s)   ∝  (1 / |$displacements|) × (build prob) × (accept)

   Because π(s) = π(s') (same particle configuration up to translation,
   same energy on a periodic lattice), detailed balance requires
   T(s → s') = T(s' → s), which fails by a factor of 2.

   WHAT THE CURRENT CHECKER SHOULD SEE
   ------------------------------------
   The current checker (no symmetry reduction) evaluates every pair
   (s_i, s_j) in the connected component.  Pairs connected by a
   purely rightward cluster translation will show a non-zero residual:
   T(s→s')π(s) − T(s'→s)π(s') = (1/|D|)·(build prob)·π ≠ 0.
   The checker will report FAIL and list the violating pairs.

   WHAT A TRANSLATION-SYMMETRY CHECKER SHOULD SEE
   -----------------------------------------------
   Translation invariance (G = Z_L × Z_L) IS preserved by the biased
   algorithm: shifting the entire lattice origin does not change which
   direction is "right."  The biased proposal is still spatially
   uniform in the sense that seed choice and displacement are drawn
   independently of absolute position.

   Consequence: within every translation orbit the DB violation is
   identical.  An orbit representative will be a pair connected by a
   rightward cluster move, and the checker will find the violation.
   A translation-based speedup would still correctly report FAIL.

   WHAT A D4 + TRANSLATION SPEEDUP CHECKER WOULD SEE  ← false negative
   -----------------------------------------------------------------------
   D4 contains 90° rotations.  Under a 90° rotation, the direction
   (+1,0) maps to (0,+1): a rightward move becomes an upward move.
   The biased algorithm is NOT D4-invariant:

     T(R₉₀·s → R₉₀·s')   [cluster moves upward]
     ≠  T(s → s')          [cluster moves rightward, biased probability]

   The D4 orbit of a "rightward move" pair contains both rightward-
   and upward-move pairs.  A D4 orbit representative may be chosen
   to be the upward-move pair.  The checker evaluates that pair,
   finds detailed balance satisfied (vertical moves are unbiased),
   and incorrectly certifies the entire orbit as passing — including
   the rightward-move pair that actually violates DB.

   This is a clean false negative: the D4-symmetry-accelerated checker
   would report PASS on a broken algorithm.  It illustrates why step 1
   of the symmetry speedup implementation (verifying G-invariance
   before using G to reduce pairs) is not optional.
   ================================================================ *)

$dir = DirectoryName[$InputFileName];
Get[$dir <> "vmmc_2d_grid.wl"];
Get[$dir <> "lj_colloidal_system.wl"];

$checkerAbstractParams = {"physLen", "epsLJ"};

(* ---- Biased displacement list ----------------------------------------
   Standard uniform box, then all rightward displacements (dx > 0) again.
   RandomChoice draws uniformly from this list, so P(dx>0, dy) = 2×P(dx<0, dy).
   Proposal symmetry P(d) = P(−d) is broken for every d with dx ≠ 0. *)
$nStep = Max[1, Round[physLen]];
$allDisps = DeleteCases[
  Flatten[Table[{dx, dy}, {dx, -$nStep, $nStep}, {dy, -$nStep, $nStep}], 1],
  {0, 0}];
$displacements = Join[$allDisps, Select[$allDisps, #[[1]] > 0 &]];


(* ================================================================
   Cluster builder  (identical to vmmc_lattice.wl / vmmc_continuous.wl)
   ================================================================ *)

$vmmcBuildCluster[state_, nGrid_, seed_, dir_] :=
  Module[{
    cluster    = {seed},
    inCluster  = <|seed -> True|>,
    queue      = {seed},
    frustrated = False,
    p, pType, pPost, pRev, nbrs, q, qType,
    eInit, eFwd, eRev, wFwd, wRev, r1, r2
  },
    While[queue =!= {} && !frustrated,
      p     = First[queue]; queue = Rest[queue];
      pType = state[[p]];
      pPost = $applyDir[p,  dir, nGrid];
      pRev  = $applyDir[p, {-dir[[1]], -dir[[2]]}, nGrid];

      nbrs = DeleteDuplicates @ Join[
               $neighborsD2[p,     nGrid],
               $neighborsD2[pPost, nGrid],
               $neighborsD2[pRev,  nGrid]];

      Do[
        q = nbrs[[k]];
        If[state[[q]] =!= 0 && !KeyExistsQ[inCluster, q],
          qType = state[[q]];

          eInit = $virtualPairEnergy[pType, qType, p,     q, nGrid];
          eFwd  = $virtualPairEnergy[pType, qType, pPost, q, nGrid];
          eRev  = $virtualPairEnergy[pType, qType, pRev,  q, nGrid];

          wFwd = Piecewise[{
              {1,                                eFwd === Infinity},
              {1 - Exp[\[Beta] (eInit - eFwd)],  eInit < eFwd}},
            0];
          wRev = Piecewise[{
              {1,                                eRev === Infinity},
              {1 - Exp[\[Beta] (eInit - eRev)],  eInit < eRev}},
            0];

          r1 = RandomReal[];
          If[r1 <= wFwd,
            r2 = RandomReal[];
            If[r2 > Piecewise[{
                  {1,
                      eFwd === Infinity && eRev === Infinity},
                  {1 - Exp[\[Beta] (eInit - eRev)],
                      eFwd === Infinity && eInit < eRev},
                  {0,
                      eFwd === Infinity},
                  {1,
                      eRev === Infinity && eInit < eFwd},
                  {Min[(1 - Exp[\[Beta] (eInit - eRev)]) /
                        (1 - Exp[\[Beta] (eInit - eFwd)]), 1],
                   eInit < eFwd && eInit < eRev},
                  {0,
                      eInit < eFwd}},
                0],
              frustrated = True; Break[],
              AppendTo[cluster, q];
              inCluster[q] = True;
              AppendTo[queue, q]
            ]
          ]
        ],
        {k, Length[nbrs]}
      ]
    ];
    If[frustrated, None, cluster]
  ]


(* ================================================================
   Algorithm
   ================================================================ *)

Algorithm[state_List] :=
  Module[{nGrid, occupied, seed, dir, cluster, newState, dest},
    nGrid    = Round[Sqrt[Length[state]]];
    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];

    seed = RandomChoice[occupied];
    dir  = RandomChoice[$displacements];   (* biased: dx>0 chosen 2× as often *)

    cluster = $vmmcBuildCluster[state, nGrid, seed, dir];
    If[cluster === None, Return[state]];

    newState = state;
    Do[newState[[cluster[[i]]]] = 0, {i, Length[cluster]}];
    Do[
      dest = $applyDir[cluster[[i]], dir, nGrid];
      If[newState[[dest]] =!= 0, Return[state, Module]];
      newState[[dest]] = state[[cluster[[i]]]],
      {i, Length[cluster]}];
    newState
  ]


(* ================================================================
   Symbolic parameters  (identical to vmmc_lattice.wl)
   ================================================================ *)

DynamicSymParams[states_List] :=
  Module[{types, nGrid, d2Vals, couplingAtoms},
    types  = Sort[DeleteCases[Union @@ states, 0]];
    nGrid  = Round[Sqrt[Length[states[[1]]]]];
    d2Vals = Module[{halfN = Floor[nGrid/2]},
      Sort @ DeleteDuplicates @ Select[
        Flatten @ Table[
          Min[dr, nGrid-dr]^2 + Min[dc, nGrid-dc]^2,
          {dr, 0, halfN}, {dc, 0, halfN}],
        0 < # <= $maxD2 &]];
    couplingAtoms = Flatten @ Table[
      If[a <= b, Table[couplingJ[a, b, d2], {d2, d2Vals}], Nothing],
      {a, types}, {b, types}];
    <|"couplings"     -> couplingAtoms,
      "numericParams" -> {}|>]
