(* ================================================================
   2D VMMC (Virtual Move Monte Carlo) on a periodic square lattice
   ================================================================

   State:   flat array of length L^2 (row-major).
            State[[s]] = 0 (empty) or k ∈ {1,...,N} (labeled particle).
            L is inferred as Sqrt[Length[state]]; only IDs whose decoded
            array has perfect-square length are accepted by BitsToState.

   Encoding: same bijective integer map as kawasaki_2d.wl.
             BitsToState additionally filters to perfect-square lengths.

   Site numbering (row-major, 1-indexed):
     1  2  ... L
     L+1 ...  2L
     ...
     (L-1)L+1 ... L^2

   Move:    Choose a random occupied site as the seed particle. Propose
            a virtual displacement of one lattice unit in a randomly chosen
            direction d ∈ {right, left, down, up}. Grow a cluster using
            the Whitelam–Geissler VMMC link probabilities:

              wFwd = max(1 − exp(e_init − e_fwd),  0)   [forward  link weight]
              wRev = max(1 − exp(e_init − e_rev),  0)   [reverse  link weight]

            where e_init is the pre-move pair energy, e_fwd is the pair
            energy after the forward virtual move, and e_rev is the pair
            energy after the reverse virtual move. For each occupied
            non-cluster neighbour q of cluster particle p:
              • with probability wFwd a link is attempted;
              • of those, fraction (1 − wRev/wFwd) are frustrated links
                → entire move is rejected;
              • fraction wRev/wFwd → q joins the cluster.
            Cluster growth continues recursively (BFS) until no new
            neighbours are found. If any frustrated link was detected the
            move is rejected; otherwise all cluster particles translate
            rigidly by one lattice unit in direction d.

   Hard-sphere exclusion:
            Two particles on the same lattice site have infinite repulsive
            energy. This is hard-coded, not controlled by any coupling
            symbol. In the link-probability calculation, same-site virtual
            energy is treated as +Infinity, giving wFwd = 1 (mandatory link
            attempt). Whether that link succeeds or is frustrated depends on
            the inter-particle coupling, as described above.

   Nearest-neighbour coupling:
            Same J_ab symbols as kawasaki_2d.wl, summed over all horizontal
            and vertical bonds. energy[] and MetropolisProb[] are called to
            evaluate the total energy change of the proposed cluster move;
            for pure VMMC the link-probability mechanism already enforces
            superdetailed balance, so MetropolisProb is not applied as a
            final acceptance gate (doing so would double-count bond energy
            changes and break detailed balance).

   References:
            S. Whitelam & P. L. Geissler, J. Chem. Phys. 127, 154101 (2007)
            L. O. Hedges, libVMMC (https://github.com/lohedges/vmmc)
   ================================================================ *)


(* ---- Bijective integer encoding (identical to kawasaki_2d.wl) ----------- *)

$cL[L_]       := $cL[L]    = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]    := $cLPre[L] = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_]:= $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

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


(* ---- 2D lattice helpers -------------------------------------------------- *)

(* Right and down neighbours (as in kawasaki_2d.wl) *)
$right2D[s_, L_] := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, (r-1)*L + Mod[c,L] + 1]
$down2D[s_, L_]  := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, Mod[r,L]*L + c]

(* Left and up neighbours, needed for VMMC reverse-direction moves *)
$left2D[s_, L_]  := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, (r-1)*L + Mod[c-2,L] + 1]
$up2D[s_, L_]    := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, Mod[r-2,L]*L + c]

(* All four nearest neighbours of site s on the L×L torus *)
$allNeighbors2D[s_, L_] := {$right2D[s,L], $left2D[s,L], $down2D[s,L], $up2D[s,L]}

(* Translate site s by direction vector {dRow, dCol} on the L×L torus *)
$applyDir[s_, {dr_, dc_}, L_] :=
  Mod[Ceiling[s/L] - 1 + dr, L]*L + Mod[Mod[s-1, L] + dc, L] + 1


(* ---- Coupling function --------------------------------------------------- *)

(* Nearest-neighbour only: couplingJ = 0 for d² ≠ 1 or if either site is a hole. *)
couplingJ[a_Integer, b_Integer, d2_Integer] :=
  If[a == 0 || b == 0 || d2 != 1, 0, $jPairSym[a, b]]

DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, {$jPairSym[a, b]}, Nothing],
        {a, types}, {b, types}]|>]


(* ---- Energy -------------------------------------------------------------- *)

(* Unique undirected bonds on the L×L torus, memoised per L.
   Sorting each {s, neighbour} pair before DeleteDuplicates ensures that the
   same undirected bond (e.g. {1,2} from s=1 and {2,1} from s=2) is counted
   exactly once for every L, including L=2 where right(s)=left(s). *)
$uniqueBonds2D[L_] := $uniqueBonds2D[L] =
  DeleteDuplicates[Sort /@ Flatten[
    Table[{{s, $right2D[s,L]}, {s, $down2D[s,L]}}, {s, L^2}], 1]]

energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    Total[couplingJ[state[[#[[1]]]], state[[#[[2]]]], 1] & /@ $uniqueBonds2D[L]]]


(* ---- VMMC helpers -------------------------------------------------------- *)

(* Virtual pair energy between a particle of type typeI at virtual site vI
   and a particle of type typeJ currently at site qSite.
   Same-site occupancy (vI == qSite) is hard-coded as Infinity, encoding
   the infinite hard-sphere repulsion that lattice exclusion imposes.
   NN adjacency gives couplingJ at d²=1; beyond NN gives 0. *)
$virtualPairEnergy[typeI_, typeJ_, vI_, qSite_, L_] :=
  Which[
    vI === qSite,                           Infinity,
    MemberQ[$allNeighbors2D[vI, L], qSite], couplingJ[typeI, typeJ, 1],
    True,                                   0]


(* ---- VMMC cluster builder ------------------------------------------------ *)

(* Grow a VMMC cluster via BFS, starting from seed site with proposed
   rigid displacement dir = {dRow, dCol} on the L×L torus.
   Returns the list of cluster site indices on success, or None if any
   frustrated link is encountered (which rejects the entire move).

   For each cluster particle p:
     pPost = p + dir  (virtual forward position)
     pRev  = p - dir  (virtual reverse position)
   For each occupied non-cluster neighbour q of p:
     eInit = J(typeP, typeQ)             current NN bond energy
     eFwd  = virtualPairEnergy(pPost, q) energy if p makes forward move
     eRev  = virtualPairEnergy(pRev,  q) energy if p makes reverse move
     wFwd  = max(1 - exp(eInit - eFwd), 0)  forward  link weight
     wRev  = max(1 - exp(eInit - eRev), 0)  reverse link weight
   Link logic (Whitelam-Geissler superdetailed balance):
     draw r1 ~ U[0,1]:
       if r1 > wFwd  → no link (neighbour stays out of cluster)
       else draw r2 ~ U[0,1]:
         if r2 > wRev/wFwd → frustrated link → reject entire move
         else              → q joins the cluster                      *)
$vmmcBuildCluster[state_, L_, seed_, dir_] :=
  Module[{
    cluster    = {seed},
    inCluster  = <|seed -> True|>,   (* Association for O(1) lookup *)
    queue      = {seed},             (* BFS processing queue *)
    frustrated = False,
    p, pType, pPost, pRev, nbrs, q, qType,
    eInit, eFwd, eRev, wFwd, wRev, r1, r2
  },
    While[queue =!= {} && !frustrated,
      (* Dequeue next cluster particle to process *)
      p     = First[queue]; queue = Rest[queue];
      pType = state[[p]];
      pPost = $applyDir[p, dir, L];                     (* virtual fwd site *)
      pRev  = $applyDir[p, {-dir[[1]], -dir[[2]]}, L];  (* virtual rev site *)
      (* Union of current, forward-virtual, and reverse-virtual neighbour shells.
         This is essential: a particle q not currently bonded to p but adjacent
         to pPost (or pRev) changes its pair energy with p upon the move and
         MUST be subjected to a link test.  Without it, bond-forming moves are
         never suppressed, violating detailed balance for repulsive (J>0) or
         mixed-sign interactions.  DeleteDuplicates prevents double-testing,
         which is also critical on L=2 where right=left and up=down. *)
      nbrs  = DeleteDuplicates[Join[
                  $allNeighbors2D[p, L],
                  $allNeighbors2D[pPost, L],
                  $allNeighbors2D[pRev, L]]];

      (* Test every particle in the combined neighbour shell for link formation *)
      Do[
        q = nbrs[[k]];
        If[state[[q]] =!= 0 && !KeyExistsQ[inCluster, q],
          qType = state[[q]];
          (* eInit via $virtualPairEnergy (not $pairJ) because q may be a future
             neighbour rather than a current one, in which case eInit = 0. *)
          eInit = $virtualPairEnergy[pType, qType, p,     q, L];  (* current bond *)
          eFwd  = $virtualPairEnergy[pType, qType, pPost, q, L]; (* fwd virtual  *)
          eRev  = $virtualPairEnergy[pType, qType, pRev,  q, L]; (* rev virtual  *)
          (* Piecewise (not Max) so FullSimplify can resolve sign-cases of J.
             Leading {1, e===Infinity} guards prevent Exp[β(finite-Infinity)]
             with symbolic β, which would produce ComplexInfinity warnings.
             Ratio case guide: eFwd=∞,eRev=∞→1; eFwd=∞,eInit<eRev→wRev;
             eFwd=∞→0; eRev=∞,eInit<eFwd→1; else→Min[wRev/wFwd,1]. *)
          wFwd = Piecewise[{
              {1,                               eFwd === Infinity},
              {1 - Exp[\[Beta] (eInit - eFwd)], eInit < eFwd}},
            0];
          wRev = Piecewise[{
              {1,                               eRev === Infinity},
              {1 - Exp[\[Beta] (eInit - eRev)], eInit < eRev}},
            0];
          r1 = RandomReal[];
          If[r1 <= wFwd,
            r2 = RandomReal[];
            If[r2 > Piecewise[{
                  {1,                               eFwd === Infinity && eRev === Infinity},
                  {1 - Exp[\[Beta] (eInit - eRev)], eFwd === Infinity && eInit < eRev},
                  {0,                               eFwd === Infinity},
                  {1,                               eRev === Infinity && eInit < eFwd},
                  {Min[(1 - Exp[\[Beta] (eInit - eRev)]) / (1 - Exp[\[Beta] (eInit - eFwd)]), 1],
                   eInit < eFwd && eInit < eRev},
                  {0,                               eInit < eFwd}},
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


(* ---- Algorithm ----------------------------------------------------------- *)

Algorithm[state_List] :=
  Module[{L, occupied, seed, dirs, dir, cluster, newState, dE},
    L = Round[Sqrt[Length[state]]];

    (* Select a random occupied site as the seed of the proposed cluster move *)
    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];
    seed = RandomChoice[occupied];

    (* Choose one of the four lattice directions at random *)
    dirs = {{0,1},{0,-1},{1,0},{-1,0}};   (* right, left, down, up *)
    dir  = RandomChoice[dirs];

    (* Attempt to build the VMMC cluster; None signals a frustrated link *)
    cluster = $vmmcBuildCluster[state, L, seed, dir];
    If[cluster === None, Return[state]];  (* frustrated link: reject move *)

    (* Apply the rigid cluster translation:
       first clear all cluster sites, then place at destinations.
       Any residual non-cluster particle at a destination site signals an
       overlap that was not resolved during cluster growth → reject.        *)
    newState = state;
    Do[newState[[cluster[[i]]]] = 0,          {i, Length[cluster]}];
    Do[With[{dest = $applyDir[cluster[[i]], dir, L]},
         If[newState[[dest]] =!= 0, Return[state, Module]];
         newState[[dest]] = state[[cluster[[i]]]]
       ], {i, Length[cluster]}];

    (* MetropolisProb called for checker interface compliance only;
       VMMC superdetailed balance is enforced by the link probabilities. *)
    dE = energy[newState] - energy[state];
    MetropolisProb[dE];

    newState
  ]


(* ---- Checker interface (identical to kawasaki_2d.wl) -------------------- *)

(* Accept only IDs whose decoded array has perfect-square length (L x L grid) *)
BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2], state, sqrtM},
    If[id == 0, Return[None]];
    state = $decode[id];
    sqrtM = Sqrt[Length[state]];
    If[!IntegerQ[sqrtM], Return[None]];
    state]

(* Display state as a grid: {1,2}|{3,0} for a 2x2 state *)
DisplayState[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table["{" <> StringRiffle[ToString /@ state[[(r-1)*L+1 ;; r*L]], ","] <> "}",
            {r, 1, L}],
      "|"]]

(* Return all integer IDs that decode to L×L (perfect-square) arrays with
   ID ≤ maxId.  check.wls uses this to skip the ~99% of bit strings that
   BitsToState would immediately reject, reducing pre-3×3 overhead from
   ~500K BitsToState calls to the small set of valid grid-state IDs. *)
ValidStateIDs[maxId_Integer] :=
  Module[{L = 1, ids = {}},
    While[$cLPre[L^2] <= maxId,
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2 + 1] - 1, maxId]]];
      L++];
    ids]

numBeta = 1
