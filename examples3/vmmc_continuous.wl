(* ================================================================
   vmmc_continuous.wl
   Continuous-limit VMMC on a periodic 2D torus — Gaussian proposal
   ================================================================

   Virtual Move Monte Carlo (Whitelam-Geissler) on a 2D periodic
   lattice with a Gaussian displacement proposal.

   Physical design:
     physLen — particle diameter in lattice units.  Sets the LJ
       length scale: sigLJ = physLen, LJ minimum at d²=2^(1/3)·physLen².
       For physLen=5 this gives strong hard-core repulsion at d²=1.

     sigStep — displacement proposal std dev.  Default value is the
       natural BD timestep at the LJ timescale:
           sigStep = physLen · √(2/(numBeta·epsLJ))
       Override sigStep directly to tune acceptance rate.

   Gaussian proposal:
     (dx, dy) drawn independently from N(0, sigStep), rounded to
     integers.  Symmetry p(dx,dy)=p(-dx,-dy) holds exactly (Gaussian
     is even).  Zero displacement (dx=dy=0) returns the state unchanged.

   Checker:
     $checkerAbstractParams = {physLen, epsLJ, sigStep} — cleared to
     unbound symbols before BFS so the proof is valid for all values.
     $maxD2 = Infinity — all site pairs included (correct for small
     lattices used by the checker).
     couplingJ has no DownValues during symbolic check — each
     (a,b,d2) triple is a free real atom.
     dbc_core.wl intercepts RandomVariate[NormalDistribution[0,sigStep]]
     and converts it to seqBernoulli over {-nMax,...,nMax} where
     nMax=floor(nGrid/2), with symbolic Erfc weights.

   Numerical / animation:
     $couplingJConcrete activates inside Block in check.wls/animate.wls.
     Override $maxD2 = Ceiling[2*physLen^2] for physical runs.

   REFERENCES
     Whitelam & Geissler, J. Chem. Phys. 127, 154101 (2007) — VMMC
     Ermak & McCammon, J. Chem. Phys. 69, 1352 (1978) — BD
   ================================================================ *)


(* ================================================================
   SECTION 0 — USER PARAMETERS
   ================================================================

   Edit this section to change the model.  Everything else is derived
   automatically.
   ================================================================ *)

(* ---- Physical parameters — concrete for numerical runs ---- *)

(* Particle diameter in lattice units.  Sets the LJ length scale:
     sigLJ = physLen   (LJ zero-crossing = one particle diameter)
     LJ minimum at d² = 2^(1/3)·physLen²
   For physLen=5: nearest-neighbour energy ≈ 4·epsLJ·(5^12−5^6)·10^9,
   strongly repulsive, enforcing a soft hard-core at one diameter. *)
physLen = 5

(* LJ well depth in kT units. *)
epsLJ = 1

sigLJ = physLen      (* LJ zero-crossing = particle diameter *)

(* Displacement proposal std dev: one BD timestep at the natural LJ timescale
   (m=1, γ=1).  Override sigStep directly to tune acceptance rate. *)
sigStep = physLen * Sqrt[2.0 / (numBeta * epsLJ)]

(* Interaction cutoff.
   $maxD2 = Infinity: include all pairs; correct for symbolic check (small
   lattices have finitely many pairs anyway).
   For numerical runs with physLen=5: set $maxD2 = Ceiling[2*physLen^2] = 50. *)
$maxD2 = Infinity

(* ---- Checker interface ---- *)

(* Parameters cleared to unbound symbols before BFS; concrete values above
   are for numerical runs only. *)
$checkerAbstractParams = {physLen, epsLJ, sigStep}

(* fixedParams: declared as Reals for FullSimplify; assigned concrete values
   (not random) during numerical MCMC.  Values are captured here before
   $checkerAbstractParams clears these symbols for BFS.
   sigStep formula uses numBeta and epsLJ; evaluated after numBeta=1 (line 405). *)
symParams = <|"fixedParams" -> <|physLen -> 5, epsLJ -> 1,
                                  sigStep -> physLen * Sqrt[2.0 / (numBeta * epsLJ)]|>|>

(* ---- Abstract-functions flag ---- *)
(* couplingJ has NO DownValues during the symbolic check; each call
   couplingJ[a,b,d2] is a free real atom.  Do NOT set this to False. *)
$abstractFunctions = True

(* ---- Concrete coupling (numerical MCMC / animation only) ---- *)
$couplingJConcrete[a_Integer, b_Integer, d2_Integer] :=
  If[d2 == 0 || d2 === Infinity || d2 > Ceiling[2 * physLen^2], 0,
     4 * epsLJ * ((sigLJ^2/d2)^6 - (sigLJ^2/d2)^3)]

$concreteParams   = <||>
$couplingFormulaStr = "4*epsLJ*((sigLJ^2/d2)^6-(sigLJ^2/d2)^3); min at d2=2^(1/3)*physLen^2"


(* ================================================================
   SECTION 1 — Bijective integer encoding   (identical to vmmc_2d.wl)
   ================================================================ *)

$cL[L_]        := $cL[L]       = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]     := $cLPre[L]    = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_] := $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

$unrankCombo[rank_, L_, N_] :=
  Module[{pos = ConstantArray[0, N], x = L - 1, r = rank},
    Do[While[Binomial[x, i] > r, x--]; pos[[i]] = x; r -= Binomial[x, i]; x--,
       {i, N, 1, -1}]; pos]

$rankPerm[perm_List] :=
  Module[{n = Length[perm], elems = Range[Length[perm]], rank = 0, idx},
    Do[idx = FirstPosition[elems, perm[[i]]][[1]] - 1;
       rank += idx * Factorial[n - i]; elems = Delete[elems, idx + 1],
       {i, n}]; rank]

$unrankPerm[k_, n_] :=
  Module[{elems = Range[n], perm = {}, r = k, idx},
    Do[idx = Quotient[r, Factorial[i - 1]]; r = Mod[r, Factorial[i - 1]];
       AppendTo[perm, elems[[idx + 1]]]; elems = Delete[elems, idx + 1],
       {i, n, 1, -1}]; perm]

$decode[id_Integer] :=
  Module[{L = 0, N = 0, r, rpos, rperm, pos, perm, arr},
    While[$cLPre[L + 1] <= id, L++];
    r = id - $cLPre[L];
    While[$cLNPre[L, N + 1] <= r, N++];
    r -= $cLNPre[L, N];
    rpos = Quotient[r, Factorial[N]]; rperm = Mod[r, Factorial[N]];
    pos  = $unrankCombo[rpos, L, N]; perm  = $unrankPerm[rperm, N];
    arr  = ConstantArray[0, L];
    Do[arr[[pos[[i]] + 1]] = perm[[i]], {i, N}]; arr]


(* ================================================================
   SECTION 2 — Periodic torus helpers
   ================================================================ *)

(* Row and column of site s on an nGrid×nGrid torus (1-indexed) *)
$row[s_, nGrid_] := Ceiling[s / nGrid]
$col[s_, nGrid_] := Mod[s - 1, nGrid] + 1

(* Translate site s by integer grid displacement {dr, dc} on the torus *)
$applyDir[s_, {dr_, dc_}, nGrid_] :=
  Mod[$row[s, nGrid] - 1 + dr, nGrid] * nGrid +
  Mod[$col[s, nGrid] - 1 + dc, nGrid] + 1

(* Minimum-image squared distance between sites s1 and s2 in grid units.
   Uses the torus minimum image convention (shortest path around the ring). *)
$torusD2[s1_, s2_, nGrid_] :=
  With[{dr0 = Abs[$row[s1, nGrid] - $row[s2, nGrid]],
        dc0 = Abs[$col[s1, nGrid] - $col[s2, nGrid]]},
    With[{dr = Min[dr0, nGrid - dr0], dc = Min[dc0, nGrid - dc0]},
      dr^2 + dc^2]]



(* ================================================================
   SECTION 4 — Interaction shells   (memoised)
   ================================================================ *)

(* All sites within squared distance $maxD2 of site s on nGrid torus *)
$neighborsD2[s_, nGrid_] := $neighborsD2[s, nGrid] =
  Select[Range[nGrid^2],
    Function[q, With[{d2 = $torusD2[s, q, nGrid]},
      d2 > 0 && d2 <= $maxD2]]]

(* All unique undirected bonds {s1, s2, d2} with s1<s2, d2≤$maxD2 *)
$uniqueBondsExt[nGrid_] := $uniqueBondsExt[nGrid] =
  Flatten[
    Table[
      With[{d2 = $torusD2[s1, s2, nGrid]},
        If[d2 > 0 && d2 <= $maxD2, {{s1, s2, d2}}, {}]],
      {s1, nGrid^2}, {s2, s1 + 1, nGrid^2}],
    2]


(* ================================================================
   SECTION 5 — Energy
   ================================================================ *)

(* Total pair energy: sum of couplingJ over all bonds within $maxD2.
   couplingJ[a,b,d2] is always called with a≤b (canonical order) so
   that the same symbolic atom is used regardless of bond orientation.
   During symbolic check: couplingJ has no DownValues → free atoms.
   During numerical MCMC: activated via Block in check.wls. *)
energy[state_List] :=
  With[{nGrid = Round[Sqrt[Length[state]]]},
    Total @ Map[
      Function[bond,
        With[{t1 = state[[bond[[1]]]], t2 = state[[bond[[2]]]], d2 = bond[[3]]},
          If[t1 != 0 && t2 != 0,
             couplingJ[Min[t1, t2], Max[t1, t2], d2], 0]]],
      $uniqueBondsExt[nGrid]]]


(* ================================================================
   SECTION 6 — Virtual pair energy for VMMC link weights
   ================================================================

   Energy between a particle of type typeI at virtual site vI and a
   particle of type typeJ at site qSite.  Returns Infinity for
   hard-core overlap (vI === qSite) and 0 outside the cutoff.
   Uses $torusD2 so it is correct for any step direction, including
   diagonals and multi-unit steps. *)

$virtualPairEnergy[typeI_, typeJ_, vI_, qSite_, nGrid_] :=
  Which[
    vI === qSite, Infinity,
    True,
      With[{d2 = $torusD2[vI, qSite, nGrid]},
        If[d2 > 0 && d2 <= $maxD2,
           couplingJ[Min[typeI, typeJ], Max[typeI, typeJ], d2], 0]]]


(* ================================================================
   SECTION 7 — VMMC cluster builder
   ================================================================

   Whitelam-Geissler link-weight logic, generalised to K directions
   and arbitrary $maxD2.  Identical in structure to vmmc_2d_field.wl:
   the neighbour shell of each cluster particle p includes all sites
   within $maxD2 of p, pPost=p+dir, and pRev=p-dir, ensuring no
   interaction is missed regardless of the direction's magnitude.

   Returns the cluster (list of site indices) or None on frustration. *)

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

      (* Union of neighbour shells of p, pPost, and pRev *)
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
   SECTION 8 — Algorithm
   ================================================================

   One VMMC step with Gaussian displacement proposal:
   1. Choose seed particle uniformly from all occupied sites.
   2. Draw (dx, dy) independently from N(0, sigStep), rounded to integers.
      If dir = {0,0}, return state (no-op move).
   3. Build cluster via Whitelam-Geissler link probabilities.
   4. Apply rigid cluster translation; reject on hard-sphere collision.

   No post-cluster Metropolis step is needed: the link-probability
   mechanism exactly accounts for all cluster–noncluster pair energy
   changes via superdetailed balance.  Intra-cluster distances are
   preserved by rigid translation, so ΔE_intra = 0 identically.

   Symmetry: p(dx,dy) = p(−dx,−dy) holds exactly because the Gaussian
   is even.  This is the only requirement for superdetailed balance. *)

Algorithm[state_List] :=
  Module[{nGrid, occupied, seed, dx, dy, dir, cluster, newState, dest},
    nGrid    = Round[Sqrt[Length[state]]];
    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];

    seed = RandomChoice[occupied];

    dx  = Round[RandomVariate[NormalDistribution[0, sigStep]]];
    dy  = Round[RandomVariate[NormalDistribution[0, sigStep]]];
    dir = {dx, dy};
    If[dir === {0, 0}, Return[state]];

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
   SECTION 9 — Dynamic symbolic parameters
   ================================================================

   Returns the symbolic coupling atoms that appear in this component's
   transition matrix, so check.wls can include them in the
   Element[#, Reals] assumptions passed to FullSimplify.

   "couplings"     — couplingJ[a,b,d2] for every canonical type pair
                     (a≤b) and every achievable d2 ∈ [1,$maxD2] on
                     this component's lattice.
   "numericParams" — scalar symbols used by $couplingJConcrete
                     (epsLJ, sigLJ, etc.) for the numerical MCMC Block. *)

DynamicSymParams[states_List] :=
  Module[{types, nGrid, d2Vals, couplingAtoms},
    types  = Sort[DeleteCases[Union @@ states, 0]];
    nGrid  = Round[Sqrt[Length[states[[1]]]]];
    (* Achievable squared distances on this torus ($maxD2=Infinity → all pairs) *)
    d2Vals = Sort @ DeleteDuplicates @ Select[
      Flatten @ Table[$torusD2[s1, s2, nGrid],
                      {s1, nGrid^2}, {s2, nGrid^2}],
      # > 0 &];
    couplingAtoms = Flatten @ Table[
      If[a <= b, Table[couplingJ[a, b, d2], {d2, d2Vals}], Nothing],
      {a, types}, {b, types}];
    <|"couplings"     -> couplingAtoms,
      "numericParams" -> {}|>]


(* ================================================================
   SECTION 10 — Checker interface   (identical to vmmc_2d_field.wl)
   ================================================================ *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2], state, sqrtM},
    If[id == 0, Return[None]];
    state = $decode[id];
    sqrtM = Sqrt[Length[state]];
    If[!IntegerQ[sqrtM], Return[None]];
    state]

DisplayState[state_List] :=
  With[{nGrid = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table["{" <> StringRiffle[ToString /@ state[[(r-1)*nGrid+1 ;; r*nGrid]], ","] <> "}",
            {r, 1, nGrid}],
      "|"]]

ValidStateIDs[maxId_Integer] :=
  Module[{L = 1, ids = {}},
    While[$cLPre[L^2] <= maxId,
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2 + 1] - 1, maxId]]];
      L++];
    ids]

numBeta = 1
