(* ================================================================
   vmmc_2d_field.wl
   2D VMMC with user-defined field and coupling functions on a
   fully-periodic (torus) square lattice.
   ================================================================

   The field f(x,y) and coupling J(a,b,d²) are supplied by the user
   in Section 0 below as concrete functions used for numerical runs.

   For the SYMBOLIC detailed-balance check the concrete definitions
   are NOT used.  Instead, fieldF and couplingJ have no DownValues
   during the check, so every site-value fieldF[x,y,L] and every
   pair-value couplingJ[a,b,d²] is treated as an independent free
   real atom by FullSimplify.  This proves detailed balance for ALL
   possible field and coupling functions simultaneously — not just
   the specific sinusoidal/exponential forms defined below.

   Concrete functions are activated (via Block) only when numerical
   MCMC verification or animation runs are performed.

   To use a different model, edit Section 0 only.

   STATE / ENCODING
      Identical to vmmc_2d.wl: flat L²-array, state[[s]]=0 (empty) or
      k∈{1,...,N} (labeled particle).  BitsToState filters to perfect-
      square lengths.

   REFERENCES
      S. Whitelam & P. L. Geissler, J. Chem. Phys. 127, 154101 (2007)
   ================================================================ *)


(* ================================================================
   SECTION 0 — USER INTERFACE: field and coupling functions
   ================================================================

   HOW THIS WORKS
   ──────────────
   Two modes are used automatically:

   SYMBOLIC (detailed-balance check):
     fieldF and couplingJ have NO DownValues.  Expressions like
     fieldF[0,0,2] and couplingJ[1,2,1] remain as independent
     unevaluated atoms — like free symbols — so the checker proves
     detailed balance for ALL possible field and coupling functions
     simultaneously, without any assumption about their specific form.

   NUMERICAL (MCMC verification and animation):
     $fieldFConcrete and $couplingJConcrete provide concrete
     implementations that are activated inside a Block.
     $concreteParams supplies numeric values for their parameters.

   TO CUSTOMISE THE MODEL
   ──────────────────────
   Edit only this section:
     $fieldFConcrete    — your field function f(x, y, L)
     $couplingJConcrete — your coupling J(a, b, d²); MUST be symmetric
     $concreteParams    — numeric values for the parameters used above
     $maxD2             — maximum squared-distance for bonds
   ================================================================ *)

(* ---- CONCRETE field: sinusoidal in the x (column) direction ---- *)
(* For an L×L lattice with 0-indexed column x ∈ {0,...,L-1}:
     x=0 → fieldAmp · Sin[π/L]
     x=1 → fieldAmp · Sin[2π/L]  ...                                *)
$fieldFConcrete[x_Integer, y_Integer, L_Integer] :=
  fieldAmp * Sin[\[Pi] * (x + 1) / L]

(* ---- CONCRETE coupling: exponential decay with per-pair amplitude ---- *)
(* $jPairSym[a,b] = Jpair<lo><hi> is symmetric by construction.
   couplingJ[a,b,d²] = Jpair<lo><hi> · exp(−λ d²)  ✓ symmetric     *)
$couplingJConcrete[a_Integer, b_Integer, d2_Integer] :=
  $jPairSym[a, b] * Exp[-lambdaJ * d2]

(* ---- Numeric parameter values for MCMC and animation runs ---- *)
$concreteParams = <|fieldAmp -> 1.0, lambdaJ -> 0.5|>

(* ---- Nearest-neighbour only by default (fastest symbolic check) ---- *)
$maxD2 = 1

(* ---- Flag: tells check.wls / report.wls / animate.wls that fieldF and
   couplingJ have no DownValues and need Block-based activation for
   numerical runs. Do NOT remove this line. ---- *)
$abstractFunctions = True

(* ---- Human-readable formula strings shown in the animation panel ---- *)
$couplingFormulaStr = "Jpair_ab * Exp[-lambdaJ * d2]"
$fieldFormulaStr    = "fieldAmp * Sin[Pi * (x+1) / L]"


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
   SECTION 2 — Periodic lattice helpers   (identical to vmmc_2d.wl)
   ================================================================ *)

(* Translate site s by direction vector {dRow, dCol} on the L×L torus *)
$applyDir[s_, {dr_, dc_}, L_] :=
  Mod[Ceiling[s/L] - 1 + dr, L]*L + Mod[Mod[s-1, L] + dc, L] + 1

(* Row and column of site s (1-indexed) *)
$row[s_, L_] := Ceiling[s/L]
$col[s_, L_] := Mod[s-1, L] + 1

(* 0-indexed coordinates for fieldF (x = column index, y = row index) *)
$siteX[s_, L_] := $col[s, L] - 1   (* x ∈ {0,...,L-1} *)
$siteY[s_, L_] := $row[s, L] - 1   (* y ∈ {0,...,L-1} *)


(* ================================================================
   SECTION 3 — Minimum-image squared distance on the L×L torus
   ================================================================ *)

$torusD2[s1_, s2_, L_] :=
  With[{dr0 = Abs[$row[s1,L] - $row[s2,L]],
         dc0 = Abs[$col[s1,L] - $col[s2,L]]},
    With[{dr = Min[dr0, L - dr0],
           dc = Min[dc0, L - dc0]},
      dr^2 + dc^2]]


(* ================================================================
   SECTION 4 — Neighbour shells (memoised)
   ================================================================

   Includes all sites within squared distance $maxD2.
   Restricting this shell to active bond distances is critical for
   checker speed: a site with zero coupling would still trigger a
   RandomReal[] draw in the cluster builder (consuming a BFS bit)
   even though the link weight is always 0. *)

$neighborsD2[s_, L_] := $neighborsD2[s, L] =
  Select[Range[L^2],
    Function[q, Module[{d2 = $torusD2[s, q, L]}, d2 > 0 && d2 <= $maxD2]]]


(* ================================================================
   SECTION 5 — Unique undirected bonds (memoised)
   ================================================================

   Each bond {s1, s2, d²} with s1 < s2 at squared distance ≤ $maxD2.
   Used by energy[] to sum all pair contributions. *)

$uniqueBondsExt[L_] := $uniqueBondsExt[L] =
  Flatten[
    Table[
      With[{d2 = $torusD2[s1, s2, L]},
        If[d2 > 0 && d2 <= $maxD2, {{s1, s2, d2}}, {}]],
      {s1, L^2}, {s2, s1+1, L^2}],
    2]


(* ================================================================
   SECTION 6 — Energy function
   ================================================================

   Total energy = pair energy + field energy.
   Pair energy: sum over all undirected bonds within $maxD2, using
     couplingJ[type_a, type_b, d²].
   Field energy: each occupied site s contributes fieldF[x, y, L]
     where x = $siteX[s,L] and y = $siteY[s,L] are 0-indexed integers. *)

energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    (* Pair energy: sum over bonds, skip holes.
       Min/Max canonicalises argument order so couplingJ is always called
       with a ≤ b — ensuring the same symbolic atom is used regardless of
       which end of a bond carries which particle type. *)
    Total[Map[Function[bond,
      With[{t1 = state[[bond[[1]]]], t2 = state[[bond[[2]]]], d2 = bond[[3]]},
        If[t1 != 0 && t2 != 0,
           couplingJ[Min[t1,t2], Max[t1,t2], d2], 0]]],
      $uniqueBondsExt[L]]] +
    (* Field energy: each occupied site contributes fieldF *)
    Total @ Table[
      If[state[[s]] != 0,
         fieldF[$siteX[s, L], $siteY[s, L], L],
         0],
      {s, L^2}]]


(* ================================================================
   SECTION 7 — Virtual pair energy for VMMC link weights
   ================================================================

   Pair energy between a particle of type typeI at virtual site vI
   and a particle of type typeJ at site qSite.
   Uses couplingJ at the torus distance; Infinity at hard-core overlap. *)

$virtualPairEnergy[typeI_, typeJ_, vI_, qSite_, L_] :=
  Which[
    vI === qSite, Infinity,
    True,
      With[{d2 = $torusD2[vI, qSite, L]},
        (* Min/Max canonicalises argument order — same convention as energy[] *)
        If[d2 > 0 && d2 <= $maxD2,
           couplingJ[Min[typeI,typeJ], Max[typeI,typeJ], d2], 0]]]


(* ================================================================
   SECTION 8 — VMMC cluster builder
   ================================================================

   Whitelam-Geissler link-weight logic.  Identical structure to
   vmmc_2d.wl; uses $neighborsD2 (respects $maxD2) and
   $virtualPairEnergy (calls couplingJ). *)

$vmmcBuildCluster[state_, L_, seed_, dir_] :=
  Module[{
    cluster   = {seed},
    inCluster = <|seed -> True|>,
    queue     = {seed},
    frustrated = False,
    p, pType, pPost, pRev, nbrs, q, qType,
    eInit, eFwd, eRev, wFwd, wRev, r1, r2
  },

    While[queue =!= {} && !frustrated,
      p     = First[queue]; queue = Rest[queue];
      pType = state[[p]];
      pPost = $applyDir[p,  dir, L];
      pRev  = $applyDir[p, {-dir[[1]], -dir[[2]]}, L];

      (* Union of neighbour shells of p, pPost, and pRev *)
      nbrs = DeleteDuplicates @ Join[
               $neighborsD2[p,     L],
               $neighborsD2[pPost, L],
               $neighborsD2[pRev,  L]];

      Do[
        q = nbrs[[k]];
        If[state[[q]] =!= 0 && !KeyExistsQ[inCluster, q],
          qType = state[[q]];

          eInit = $virtualPairEnergy[pType, qType, p,     q, L];
          eFwd  = $virtualPairEnergy[pType, qType, pPost, q, L];
          eRev  = $virtualPairEnergy[pType, qType, pRev,  q, L];

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
   SECTION 9 — Algorithm
   ================================================================ *)

Algorithm[state_List] :=
  Module[{L, occupied, seed, dirs, dir, cluster, newState, dEField, dE},
    L = Round[Sqrt[Length[state]]];

    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];
    seed = RandomChoice[occupied];

    dirs = {{0,1},{0,-1},{1,0},{-1,0}};
    dir  = RandomChoice[dirs];

    cluster = $vmmcBuildCluster[state, L, seed, dir];
    If[cluster === None, Return[state]];

    newState = state;
    Do[newState[[cluster[[i]]]] = 0, {i, Length[cluster]}];
    Do[With[{dest = $applyDir[cluster[[i]], dir, L]},
         If[newState[[dest]] =!= 0, Return[state, Module]];
         newState[[dest]] = state[[cluster[[i]]]]
       ], {i, Length[cluster]}];

    (* Post-cluster Metropolis for field energy change.
       ΔE_field = Σ_{p in cluster} [fieldF(dest(p)) − fieldF(p)]
       fieldF and $siteX/$siteY always produce integer arguments. *)
    dEField = Total @ Table[
      fieldF[$siteX[$applyDir[cluster[[i]], dir, L], L],
             $siteY[$applyDir[cluster[[i]], dir, L], L], L] -
      fieldF[$siteX[cluster[[i]], L],
             $siteY[cluster[[i]], L], L],
      {i, Length[cluster]}];
    If[RandomReal[] >= MetropolisProb[dEField], Return[state]];

    dE = energy[newState] - energy[state];
    MetropolisProb[dE];

    newState
  ]


(* ================================================================
   SECTION 10 — Checker interface   (identical to vmmc_2d.wl)
   ================================================================ *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2], state, sqrtM},
    If[id == 0, Return[None]];
    state = $decode[id];
    sqrtM = Sqrt[Length[state]];
    If[!IntegerQ[sqrtM], Return[None]];
    state]

DisplayState[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table["{" <> StringRiffle[ToString /@ state[[(r-1)*L+1 ;; r*L]], ","] <> "}",
            {r, 1, L}],
      "|"]]

ValidStateIDs[maxId_Integer] :=
  Module[{L = 1, ids = {}},
    While[$cLPre[L^2] <= maxId,
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2+1]-1, maxId]]];
      L++];
    ids]

(* DynamicSymParams: returns two sets of parameters.
   "couplings" — symbolic atoms for the FullSimplify-based checker:
     fieldF[x,y,L] for every site on this component's lattice, and
     couplingJ[a,b,d2] for every active type pair (a < b) and bond
     distance d2 ∈ {1,...,$maxD2}.  Because fieldF and couplingJ have
     no DownValues these remain as independent unevaluated atoms — each
     treated as a distinct free real by FullSimplify.
   "numericParams" — scalar symbols for the numerical MCMC Block:
     the $jPairSym coupling symbols plus the parameters that appear in
     $fieldFConcrete and $couplingJConcrete ($concreteParams keys).
     These are assignable symbols that the Block mechanism can bind to
     numeric values. *)
DynamicSymParams[states_List] :=
  Module[{types, L, fieldAtoms, couplingAtoms, jPairSyms, concreteKeys},
    types = Sort[DeleteCases[Union @@ states, 0]];
    (* All states in a connected component share the same lattice size *)
    L = Round[Sqrt[Length[states[[1]]]]];
    (* Independent field atoms: one per lattice site *)
    fieldAtoms = Flatten @ Table[
      fieldF[x, y, L],
      {x, 0, L - 1}, {y, 0, L - 1}];
    (* Independent coupling atoms: one per canonical pair (a<b) per d2 *)
    couplingAtoms = Flatten @ Table[
      If[a < b, Table[couplingJ[a, b, d2], {d2, 1, $maxD2}], Nothing],
      {a, types}, {b, types}];
    (* $jPairSym symbols needed by $couplingJConcrete *)
    jPairSyms = Flatten @ Table[
      If[a < b, {$jPairSym[a, b]}, Nothing],
      {a, types}, {b, types}];
    (* Parameters used by the concrete implementations *)
    concreteKeys = If[AssociationQ[$concreteParams], Keys[$concreteParams], {}];
    <|"couplings"     -> Join[fieldAtoms, couplingAtoms],
      "numericParams" -> Join[jPairSyms, concreteKeys]|>]

numBeta = 1
