(* ================================================================
   vmmc_2d_grid.wl — Shared 2D square-lattice infrastructure
   ================================================================
   Load from an algorithm file with:
     Get[DirectoryName[$InputFileName] <> "vmmc_2d_grid.wl"]

   Provides:
     - Periodic torus geometry helpers ($row, $col, $applyDir, $torusD2)
     - Neighbour shell cache ($neighborsD2)
     - Checker interface (BitsToState, DisplayState, ValidStateIDs)

   The bijective integer encoding ($decode etc.) lives in dbc_core.wl
   and is loaded below if not already available.

   THIS FILE IS 2D SQUARE LATTICE SPECIFIC.
   For 3D or non-square geometries, replace this file.  The BFS core
   in dbc_core.wl is geometry-agnostic.
   ================================================================ *)

(* Load dbc_core.wl if not already loaded (needed when called from
   animate.wls, which does not load dbc_core.wl directly). *)
If[!ValueQ[$dbcDir],
  Get[DirectoryName[$InputFileName] <> "../dbc_core.wl"]];


(* ================================================================
   Periodic torus helpers
   ================================================================ *)

$row[s_, nGrid_] := Ceiling[s / nGrid]
$col[s_, nGrid_] := Mod[s - 1, nGrid] + 1

$applyDir[s_, {dr_, dc_}, nGrid_] :=
  Mod[$row[s, nGrid] - 1 + dr, nGrid] * nGrid +
  Mod[$col[s, nGrid] - 1 + dc, nGrid] + 1

(* Minimum-image squared distance between sites s1 and s2. *)
$torusD2[s1_, s2_, nGrid_] :=
  With[{dr0 = Abs[$row[s1, nGrid] - $row[s2, nGrid]],
        dc0 = Abs[$col[s1, nGrid] - $col[s2, nGrid]]},
    With[{dr = Min[dr0, nGrid - dr0], dc = Min[dc0, nGrid - dc0]},
      dr^2 + dc^2]]


(* ================================================================
   Neighbour shell cache  (memoised)
   ================================================================
   All sites within squared distance $maxD2 of site s on an
   nGrid×nGrid torus.  $maxD2 must be set before first call. *)

$neighborsD2[s_, nGrid_] := $neighborsD2[s, nGrid] =
  If[TrueQ[$maxD2 === Infinity],
    Select[Range[nGrid^2], # =!= s &],
    Module[{rMax = Ceiling[Sqrt[$maxD2]]},
      DeleteDuplicates @ Select[
        Flatten @ Table[$applyDir[s, {dr, dc}, nGrid],
                        {dr, -rMax, rMax}, {dc, -rMax, rMax}],
        Function[q, q =!= s && $torusD2[s, q, nGrid] <= $maxD2]]]]


(* ================================================================
   Canonical neighbour oracle for symmetry-aware checking
   ================================================================
   $dbcCanonicalCandidates returns the VMMC candidate neighbours of
   cluster particle p in a canonical (symmetry-invariant) order.

   The three distance arguments d²(p,q), d²(pPost,q), d²(pRev,q) are
   topological scalars preserved by every isometry of the torus (all
   translations, and all D4 rotations/reflections).  Sorting by the
   tuple (d²_init, d²_fwd, d²_rev, type) therefore assigns the same
   ordering to geometrically G-related states, so the seqBernoulli
   random-number trees produced by $vmmcBuildCluster are syntactically
   identical for G-orbit pairs.  The existing $dbcDedup hash step then
   collapses those pairs automatically, giving an ~|G|-fold speedup
   in FullSimplify calls with no extra verification logic.

   Arguments:
     p         — site index of the current cluster particle
     pPost     — site p would move to ($applyDir[p, dir, nGrid])
     pRev      — site p would come from ($applyDir[p, -dir, nGrid])
     state     — current state vector
     nGrid     — grid side length
     inCluster — Association of sites already in the cluster

   Returns a list of site indices: occupied, not in cluster, sorted
   by the canonical key (d²(p,q), d²(pPost,q), d²(pRev,q), type).

   This function should be used in $vmmcBuildCluster for algorithm
   files that declare $symmetryGroup.  It replaces the raw
   DeleteDuplicates@Join[$neighborsD2[...]] computation + If-guard. *)

$dbcCanonicalCandidates[p_, pPost_, pRev_, state_, nGrid_, inCluster_] :=
  Module[{rawNbrs, occupied},
    rawNbrs = DeleteDuplicates @ Join[
                $neighborsD2[p,     nGrid],
                $neighborsD2[pPost, nGrid],
                $neighborsD2[pRev,  nGrid]];
    occupied = Select[rawNbrs,
      Function[q, state[[q]] =!= 0 && !KeyExistsQ[inCluster, q]]];
    SortBy[occupied, {
      $torusD2[p,    #, nGrid] &,
      $torusD2[pPost, #, nGrid] &,
      $torusD2[pRev,  #, nGrid] &,
      state[[#]] &}]]


(* ================================================================
   Checker interface
   ================================================================ *)

(* BitsToState: interpret bit string as integer id → decode to state.
   Returns None for id=0 or non-square lattice sizes. *)
BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2], state, sqrtM},
    If[id == 0, Return[None]];
    state = $decode[id];
    sqrtM = Sqrt[Length[state]];
    If[!IntegerQ[sqrtM], Return[None]];
    state]

(* DisplayState: human-readable row-by-row grid string. *)
DisplayState[state_List] :=
  With[{nGrid = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table["{" <> StringRiffle[ToString /@ state[[(r-1)*nGrid+1 ;; r*nGrid]], ","] <> "}",
            {r, 1, nGrid}],
      "|"]]

(* ValidStateIDs: restrict bit-string enumeration to valid square-lattice IDs. *)
ValidStateIDs[maxId_Integer] :=
  Module[{L = 1, ids = {}},
    While[$cLPre[L^2] <= maxId,
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2 + 1] - 1, maxId]]];
      L++];
    ids]
