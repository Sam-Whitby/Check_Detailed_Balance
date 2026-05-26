(* ================================================================
   1D Kawasaki dynamics on a periodic ring -- PASSES detailed balance
   ================================================================

   State:   flat array of length L.
            State[[i]] = 0 (empty) or k ∈ {1,...,N} (labeled particle).
            L and N are inferred from the bit string by BitsToState.

   Encoding: bit string → integer ID → unique labeled array.
             ID encodes (L, N, particle positions, label permutation).
             See $decode below.

   Move:    Pick bond b ∈ {0,...,L-1} uniformly; swap sites b+1 and
            Mod[b+1,L]+1; accept/reject with Metropolis.
            The proposal is symmetric: every bond has equal weight.

   Energy:  Nearest-neighbour pairwise coupling J_ab between particle
            types a < b. Holes (type 0) do not contribute to energy.
   ================================================================ *)


(* ---- Bijective integer encoding ----------------------------------------- *)
(* Maps every labeled-particle array of any length to a unique integer.
   ID = countPrefix[L] + countNPrefix[L,N] + rankCombo * N! + rankPerm
   where countPrefix accumulates array counts over all smaller lengths. *)

$cL[L_]       := $cL[L]    = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]    := $cLPre[L] = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_]:= $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

(* Combinatorial number system: rank/unrank sorted 0-indexed N-combinations *)
$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

$unrankCombo[rank_, L_, N_] :=
  Module[{pos = ConstantArray[0, N], x = L - 1, r = rank},
    Do[While[Binomial[x, i] > r, x--]; pos[[i]] = x; r -= Binomial[x, i]; x--,
       {i, N, 1, -1}]; pos]

(* Factorial number system: rank/unrank permutations of {1,...,N} *)
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

energy[state_List] :=
  With[{L = Length[state]},
    Total[Table[couplingJ[state[[i]], state[[Mod[i, L] + 1]], 1], {i, L}]]]


(* ---- Algorithm ----------------------------------------------------------- *)

Algorithm[state_List] :=
  Module[{L, b, s1, s2, newState, dE},
    L  = Length[state];
    b  = RandomInteger[{0, L - 1}];          (* bond index *)
    s1 = b + 1; s2 = Mod[b + 1, L] + 1;     (* the two sites *)
    newState = ReplacePart[state, {s1 -> state[[s2]], s2 -> state[[s1]]}];
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ---- Checker interface --------------------------------------------------- *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2]},
    If[id == 0, None, $decode[id]]]

numBeta = 1
