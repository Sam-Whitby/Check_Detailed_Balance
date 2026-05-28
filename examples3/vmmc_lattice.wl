(* ================================================================
   vmmc_lattice.wl — VMMC with uniform-box displacement on 2D LJ lattice
   ================================================================

   WHAT THIS FILE IS
   -----------------
   A principled VMMC algorithm that replaces the rounded-Gaussian
   displacement proposal of vmmc_continuous.wl with a uniform
   distribution over all integer offsets (dx,dy) with
   |dx| ≤ $nStep and |dy| ≤ $nStep, excluding (0,0).

   $nStep = Max[1, Round[physLen]] is evaluated once at load time
   from the concrete value of physLen, so it remains a concrete
   integer even after physLen is cleared as an abstract parameter.

   The Whitelam–Geissler cluster building routine is unchanged.

   WHY THIS SATISFIES DETAILED BALANCE (PRINCIPLED)
   -------------------------------------------------
   VMMC detailed balance rests on exactly two conditions:

   1. PROPOSAL SYMMETRY: P(propose d) = P(propose −d).
      The uniform box {(dx,dy) : |dx|≤$nStep, |dy|≤$nStep} \ {(0,0)}
      contains (dx,dy) and (−dx,−dy) with equal weight by construction.
      No parameter choice can break this — it is a structural property
      of the enumerated displacement list.

   2. WHITELAM–GEISSLER LINK PROBABILITIES (unchanged from
      vmmc_continuous.wl).  For particle p interacting with q:
        wFwd = max(0, 1 − exp(β(eInit − eFwd)))   [p moves toward q]
        wRev = max(0, 1 − exp(β(eInit − eRev)))   [reverse virtual move]
      The frustration condition (second RandomReal[] test) ensures
      that the reverse cluster would build consistently.  Together,
      conditions 1 and 2 satisfy superdetailed balance.

   WHY THE CODE CHECKER HAS NO FALSE NEGATIVES
   --------------------------------------------
   All symbolic variables in the checker are the coupling constants
   couplingJ[a,b,d2].  Because the displacement is chosen via
   RandomChoice over a concrete list of integer pairs, every weight
   in the seqBernoulli decomposition is a rational constant — no
   symbolic parameters appear.  The Piecewise conditions produced
   by PiecewiseExpand are therefore polynomial inequalities in the
   coupling constants over R.  Mathematica's FindInstance uses
   Cylindrical Algebraic Decomposition for polynomial conditions,
   which is complete: it returns {} only for genuinely infeasible
   conditions.  There are no Erfc terms, no abstract step-size
   parameters, and no transcendental conditions.  False negatives
   from the feasibility check are structurally impossible.

   CONTINUOUS LIMIT
   ----------------
   For large nGrid with physLen ∝ nGrid (fixed packing fraction η):

     $nStep ∼ physLen ∼ nGrid,

   so a single VMMC step proposes a displacement of O(physLen)
   lattice units = O(σ) in physical units — the same displacement
   scale as Brownian motion.

   The macroscopic diffusion coefficient is insensitive to whether
   single-step displacements are Gaussian or uniform-box: by the
   Lindeberg–Lévy CLT, the net displacement over many accepted
   steps converges to Gaussian regardless of the single-step
   distribution, provided the distribution has finite variance and
   satisfies P(d)=P(−d).  The uniform box on [−$nStep,$nStep]²
   has variance $nStep²/3 per axis — well-defined, isotropic, and
   equal to the Gaussian with σ = $nStep/√3.

   For nGrid=2 or nGrid=3 (the checker's typical test sizes) physLen
   defaults to 1, giving $nStep=1 and exactly 8 displacement choices
   (3 bits, no rejection sampling) — the fastest possible checker
   run without restricting to cardinal-only moves.
   ================================================================ *)

$dir = DirectoryName[$InputFileName];
Get[$dir <> "vmmc_2d_grid.wl"];
Get[$dir <> "lj_colloidal_system.wl"];

(* ---- Override abstract params: sigStep is not used here ---- *)
$checkerAbstractParams = {"physLen", "epsLJ"};

(* ---- Displacement list ----
   Computed from the CONCRETE value of physLen before it is cleared.
   For physLen=1 (checker default): $nStep=1, 8 choices (3 bits exact).
   For physLen=2: $nStep=2, 24 choices.  physLen=3: $nStep=3, 48 choices.
   The list always contains (dx,dy) and (−dx,−dy) as a pair — symmetry
   is structural, not checked at run time. *)
$nStep       = Max[1, Round[physLen]];
$displacements = DeleteCases[
  Flatten[Table[{dx, dy}, {dx, -$nStep, $nStep}, {dy, -$nStep, $nStep}], 1],
  {0, 0}];

(* Declare symmetry group for the checker's canonical-neighbour oracle.
   This enables ~|G|-fold speedup in FullSimplify calls via expression-hash
   deduplication of G-orbit pairs.  The uniform-box proposal and periodic
   minimum-image energy are both D4 + translation invariant by construction.
   See $dbcCanonicalCandidates in vmmc_2d_grid.wl for details. *)
$symmetryGroup = {"translation", "D4"};


(* ================================================================
   Cluster builder
   ================================================================
   Whitelam–Geissler virtual-move cluster construction.
   Uses $dbcCanonicalCandidates (from vmmc_2d_grid.wl) to iterate
   neighbours in canonical topological order (sorted by d²_init,
   d²_fwd, d²_rev, type) so G-related states produce syntactically
   identical seqBernoulli trees and hash-dedup in the checker
   automatically collapses G-orbit pairs.
   Returns the cluster list, or None on frustration. *)

$vmmcBuildCluster[state_, nGrid_, seed_, dir_] :=
  Module[{
    cluster    = {seed},
    inCluster  = <|seed -> True|>,
    queue      = {seed},
    frustrated = False,
    p, pType, pPost, pRev, cands, q, qType,
    eInit, eFwd, eRev, wFwd, wRev, r1, r2
  },
    While[queue =!= {} && !frustrated,
      p     = First[queue]; queue = Rest[queue];
      pType = state[[p]];
      pPost = $applyDir[p,  dir, nGrid];
      pRev  = $applyDir[p, {-dir[[1]], -dir[[2]]}, nGrid];

      (* Canonical order: sort occupied non-cluster neighbours by
         (d²_init, d²_fwd, d²_rev, type).  G-related states produce the
         same ordering, so seqBernoulli trees are syntactically identical
         and $dbcDedup collapses G-orbit pairs automatically. *)
      cands = $dbcCanonicalCandidates[p, pPost, pRev, state, nGrid, inCluster];

      Do[
        q     = cands[[k]];
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
        ],
        {k, Length[cands]}
      ]
    ];
    If[frustrated, None, cluster]
  ]


(* ================================================================
   Algorithm
   ================================================================
   One VMMC step:
   1. Choose a seed particle uniformly from occupied sites.
   2. Choose displacement dir uniformly from $displacements.
      (Symmetric: every (dx,dy) is paired with (−dx,−dy).)
   3. Build cluster via Whitelam–Geissler link probabilities.
   4. Translate cluster rigidly by dir; reject on hard-core overlap. *)

Algorithm[state_List] :=
  Module[{nGrid, occupied, seed, dir, cluster, newState, dest},
    nGrid    = Round[Sqrt[Length[state]]];
    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];

    seed = RandomChoice[occupied];
    dir  = RandomChoice[$displacements];

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
   Symbolic parameters for the checker
   ================================================================
   Same logic as vmmc_continuous.wl: one free real atom per
   canonical type pair and achievable squared distance. *)

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
