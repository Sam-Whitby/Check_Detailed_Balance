(* ================================================================
   DetailedBalanceChecker  -  Core Library  (Approach 3: Interval Tracking)
   ================================================================
   Load with:  Get["path/to/dbc_core.wl"]

   Entry point:
     RunFullCheck[allStates, symAlg, numAlg, symEnergy, numEnergy, opts]

   Key change vs the original library:
     RandomReal[] is now modelled via interval tracking rather than a
     single independent Bernoulli trial per comparison.  Each call to
     RandomReal[] creates a fresh latent variable U ~ Uniform[0,1]
     represented by a token $dbc$irand[j, acceptTestI] that carries its
     own interval [lo, hi], initially {0,1}.

     When the token is compared to a threshold p (e.g. U < p), the
     conditional probability P(U < p | U in [lo,hi]) = (p-lo)/(hi-lo)
     is used to weight the bit, and the interval is narrowed to [lo,p]
     (accepted) or [p,hi] (rejected).  Subsequent comparisons of the
     same variable correctly condition on all previous comparisons.

     This implements inverse-CDF / sequential conditioning exactly, so
     chains that compare the same RandomReal[] token multiple times
     (or use RandomChoice[weights->elements] decomposed into sequential
     Bernoulli trials) are handled correctly and symbolically.

   RandomChoice[weights -> elements] is now supported:
     It is decomposed into sequential Bernoulli trials, each using a
     fresh RandomReal[] interval-tracking variable.  The probability of
     selecting element k is exactly weights[[k]] / Total[weights].

   See README.md for the full interface description.
   ================================================================ *)


$dbcDir = DirectoryName[$InputFileName];

(* ----------------------------------------------------------------
   $jPairSym
   Returns the canonical per-type-pair coupling symbol Jpair<lo><hi>.
   Symmetric by construction: $jPairSym[a,b] = $jPairSym[b,a].
   Available to all algorithm files; defined here so it does not need
   to be redefined in each .wl file.
   ---------------------------------------------------------------- *)
$jPairSym[a_Integer, b_Integer] :=
  ToExpression["Jpair" <> ToString[Min[a, b]] <> ToString[Max[a, b]]]


(* ================================================================
   SECTION 0 – SYSTEM PARAMETER UTILITIES
   ================================================================ *)

(* ----------------------------------------------------------------
   RingDist
   Minimum image convention distance between sites a and b on an
   L-site periodic ring: min(|a-b|, L-|a-b|).
   ---------------------------------------------------------------- *)
RingDist[a_Integer, b_Integer, L_Integer] :=
  Min[Abs[a - b], L - Abs[a - b]]

(* ----------------------------------------------------------------
   BuildRingEnergy
   Constructs an energy function for particles on an L-site periodic
   ring with optional pairwise interactions.

   params = Association with keys:
     "L"         -> lattice size (integer)
     "eps"       -> list of length L (site energies, symbolic or numeric)
     "couplings" -> list of coupling strengths indexed by ring distance d
                    couplings[[d]] = coupling at distance d  (d=1,2,...)
                    Hard-sphere exclusion (d=0) is handled by the
                    algorithm itself, not the energy function.

   For single-particle state (Integer):
     E = eps[[s]]

   For multi-particle state (List of integers):
     E = sum_i eps[[p_i]]
       + sum_{i<j} couplings[[RingDist[p_i, p_j, L]]]
       (distances beyond Length[couplings] contribute 0)
   ---------------------------------------------------------------- *)
BuildRingEnergy[params_Association] :=
  With[{L         = params["L"],
        eps       = params["eps"],
        couplings = Lookup[params, "couplings", {}]},
    Function[state,
      Which[
        IntegerQ[state],
          eps[[state]],
        ListQ[state],
          Total[eps[[#]] & /@ state] +
          Total[
            Table[
              With[{d = RingDist[state[[i]], state[[j]], L]},
                If[d >= 1 && d <= Length[couplings], couplings[[d]], 0]],
              {i, Length[state]}, {j, i + 1, Length[state]}],
            2]
      ]
    ]
  ]

(* ----------------------------------------------------------------
   MakeRingParams
   Creates a symbolic parameter Association for an L-site ring.
   prefix    = short string making symbol names unique (e.g. "rk")
   nCoupling = number of coupling distances to model (default 0)
   Returns <|"L" -> L, "eps" -> {...}, "couplings" -> {...}|>
   where eps and couplings contain unassigned globally-unique symbols.
   Symbol names: \[Epsilon]<prefix><i> for site i,
                 J<prefix><d> for coupling at distance d.
   ---------------------------------------------------------------- *)
MakeRingParams[L_Integer, prefix_String, nCoupling_Integer : 0] :=
  <|"L"         -> L,
    "eps"       -> Table[
                     ToExpression["\[Epsilon]" <> prefix <> ToString[i]],
                     {i, L}],
    "couplings" -> Table[
                     ToExpression["J" <> prefix <> ToString[d]],
                     {d, 1, nCoupling}]|>

(* ----------------------------------------------------------------
   MakeNumericSubs
   Generates reproducible random numerical substitution rules for the
   symbolic parameters in params (as returned by MakeRingParams).
   Site energies: drawn uniformly from (-2, 2).
   Couplings:     drawn uniformly from (-1, 1).
   seed = integer for SeedRandom reproducibility.
   Returns a list of rules {sym -> numVal, ...}.
   ---------------------------------------------------------------- *)
MakeNumericSubs[params_Association, seed_Integer : 1] := Module[
  {eps       = Lookup[params, "eps",       {}],
   couplings = Lookup[params, "couplings", {}]},
  SeedRandom[seed];
  Join[
    Thread[eps       -> RandomReal[{-2, 2}, Length[eps]]],
    Thread[couplings -> RandomReal[{-1, 1}, Length[couplings]]]
  ]
]


(* ================================================================
   SECTION 1 – PRIMITIVES
   ================================================================ *)

(* ----------------------------------------------------------------
   MetropolisProb
   Standard Metropolis acceptance as a symbolic Piecewise.
   deltaE is the BARE energy difference (no beta); the global
   symbol \[Beta] appears in the result and stays unassigned during
   the symbolic check.
   ---------------------------------------------------------------- *)
MetropolisProb[deltaE_] :=
  Piecewise[{{1, deltaE <= 0}, {Exp[-\[Beta] deltaE], deltaE > 0}}]

(* ----------------------------------------------------------------
   RunWithBits
   Run alg[state, readBit] against a fixed bit list.
   readBit[] returns bits[[1]], bits[[2]], ... and throws $OutOfBits
   if the list is exhausted.
   Returns  {outcomes, nBitsConsumed}  or  $OutOfBits.
   outcomes = {{p1,s1},{p2,s2},...}
   ---------------------------------------------------------------- *)
RunWithBits[alg_, state_, bits_List] := Module[
  {pos = 0, readBit, raw},
  readBit[] := (
    pos++;
    If[pos > Length[bits], Throw[$OutOfBits, $dbc$tag], bits[[pos]]]
  );
  raw = Catch[alg[state, readBit], $dbc$tag, ($OutOfBits &)];
  If[raw === $OutOfBits, $OutOfBits,
    {If[ListQ[raw] && Length[raw] > 0 && ListQ[raw[[1]]],
       raw, {{1, raw}}], pos}
  ]
]


(* ================================================================
   SECTION 2 – TREE BUILDING
   ================================================================ *)

(* ----------------------------------------------------------------
   BuildTreeData
   BFS over bit sequences for every starting state.
   Returns  Association[ state -> { {bits, outcomes}, ... } ]
   where outcomes = {{p1,s1},...}  (raw leaves of the decision tree).
   ---------------------------------------------------------------- *)
Options[BuildTreeData] = {
  "MaxBitDepth" -> 20,
  "TimeLimit"   -> 60.,
  "Verbose"     -> True
}

BuildTreeData[allStates_List, alg_, OptionsPattern[]] := Module[
  {maxDepth = OptionValue["MaxBitDepth"],
   tlim     = N @ OptionValue["TimeLimit"],
   verbose  = OptionValue["Verbose"],
   result   = <||>,
   queue, bits, res, outcomes, k, t0, timedOut, leaves},

  Do[
    If[verbose, Print["  Tree for state: ", s]];
    queue    = {{}};
    leaves   = {};
    t0       = AbsoluteTime[];
    timedOut = False;

    While[queue =!= {} && !timedOut,
      If[AbsoluteTime[] - t0 > tlim, timedOut = True; Break[]];
      bits  = First[queue]; queue = Rest[queue];
      res   = RunWithBits[alg, s, bits];
      Which[
        res === $OutOfBits && Length[bits] < maxDepth,
          queue = Join[queue, {Append[bits, 0], Append[bits, 1]}],
        res === $OutOfBits,
          Print["  WARNING: MaxBitDepth=", maxDepth,
                " reached for state ", s, " at prefix ", bits,
                " -- path excluded (algorithm may not halt on this input)."],
        True,
          {outcomes, k} = res;
          AppendTo[leaves, {bits, outcomes}]
      ]
    ];

    If[timedOut, Print["  WARNING: Time limit reached for state ", s]];
    result[s] = leaves,
    {s, allStates}
  ];
  result
]

(* ----------------------------------------------------------------
   TreeDataToMatrix
   Derive the transition matrix from raw tree leaves.
   Returns  Association[ {from,to} -> symbolicProbability ]
   ---------------------------------------------------------------- *)
TreeDataToMatrix[allStates_List, treeData_Association] := Module[
  {stateSet = Association[# -> True & /@ allStates], matrix = <||>},
  Do[
    Do[
      With[{bits = leaf[[1]], outcomes = leaf[[2]]},
        Do[
          With[{p = out[[1]], ns = out[[2]]},
            If[KeyExistsQ[stateSet, ns],
              matrix[{s, ns}] =
                Lookup[matrix, Key[{s, ns}], 0] + p * (1/2)^Length[bits],
              Print["  WARNING: Invalid state ", ns, " returned from ", s]
            ]
          ],
          {out, outcomes}
        ]
      ],
      {leaf, treeData[s]}
    ],
    {s, allStates}
  ];
  matrix
]

(* BuildTransitionMatrix kept for API compatibility *)
Options[BuildTransitionMatrix] = Options[BuildTreeData]
BuildTransitionMatrix[allStates_List, alg_, opts : OptionsPattern[]] :=
  TreeDataToMatrix[allStates,
    BuildTreeData[allStates, alg,
      "MaxBitDepth" -> OptionValue["MaxBitDepth"],
      "TimeLimit"   -> OptionValue["TimeLimit"],
      "Verbose"     -> OptionValue["Verbose"]]]


(* ================================================================
   SECTION 1b – GENERALISED PRIMITIVES  (interval-tracking interface)
   ================================================================ *)

(* ----------------------------------------------------------------
   RunWithBitsAT
   Runs the algorithm on a fixed bit tape, intercepting all native
   Mathematica random calls (RandomReal, RandomInteger, RandomChoice).

   The algorithm takes ONE argument:  alg[state]
   It uses native random calls, which are shadowed via Block to be
   deterministic given the bit tape.

   RandomReal[] / Random[] use INTERVAL TRACKING (Approach 3):
     Each call returns a fresh token $dbc$irand[j, acceptTestI] that
     represents a latent variable U_j ~ Uniform[0,1] with a current
     interval [lo_j, hi_j] (initially {0,1}).  When U_j is compared
     to a threshold p via an UpValue-dispatched acceptTestI call:
       - conditional probability condP = (p - lo)/(hi - lo) is computed
       - one bit is read from the tape
       - bit=1 (accept): weight *= condP, interval -> [lo, p]
       - bit=0 (reject): weight *= 1-condP, interval -> [p, hi]
     Multiple comparisons on the same variable correctly condition on
     all prior narrowings, implementing inverse-CDF sampling exactly.

   RandomChoice[weights -> elements] decomposes into sequential
   Bernoulli trials using independent fresh RandomReal[] variables,
   so that P(element k) = weights[[k]] / Total[weights] exactly.

   Returns  {nextState, pathWeight}  or  $OutOfBits.
   pathWeight = prod_{RandomInteger bits}(1/2) *
                prod_{RandomReal comparisons}(condP or 1-condP)
   ---------------------------------------------------------------- *)
RunWithBitsAT[alg_, state_, bits_List] := Module[
  {pos = 0, weight = 1,
   nReals = 0, intervals = {},
   readBit, acceptTestI, makeRealVar, seqBernoulli, makeContToken, result},

  (* --- Fair bit read: each bit contributes factor 1/2 --- *)
  readBit[] := (
    pos++;
    If[pos > Length[bits],
      Throw[$OutOfBits, $dbc$tag],
      weight *= (1/2); bits[[pos]]
    ]
  );

  (* --- Interval-tracking acceptance test for variable j ---
     Computes the conditional probability that U_j < p given that
     U_j lies in the current interval [lo, hi].
     Returns True (accept, U < p) or False (reject, U >= p).
     Updates the interval and path weight.

     Deterministic short-circuit (no bit consumed):
       p <= lo  =>  U >= p with certainty  =>  Return[False]
       p >= hi  =>  U <  p with certainty  =>  Return[True]

     For a fresh variable lo=0, hi=1 these fire whenever the threshold
     p is exactly 0 (reject certainly) or exactly 1 (accept certainly).
     This covers e.g. wFwd=0 for zero-coupling neighbour pairs, and
     wFwd=1 for hard-sphere exclusion, without consuming any BFS bit.

     To handle thresholds that are symbolically 0/1 but not yet reduced
     to a number — e.g. Max[0,...] rather than Piecewise, or a Piecewise
     whose conditions have not yet been evaluated — we first attempt a
     lightweight PiecewiseExpand, and then use TrueQ so that an
     unevaluated symbolic comparison never triggers the short-circuit
     incorrectly. *)
  acceptTestI[j_, p_] := Module[{lo, hi, condP, pR, pVal},
    lo = intervals[[j, 1]];
    hi = intervals[[j, 2]];
    (* Attempt to reduce symbolic threshold without full Simplify overhead *)
    pVal = If[NumericQ[p], p, PiecewiseExpand[p]];
    (* Deterministic cases: no bit consumed, interval unchanged *)
    If[TrueQ[pVal <= lo], Return[False, Module]];
    If[TrueQ[pVal >= hi], Return[True,  Module]];
    condP = (pVal - lo) / (hi - lo);
    pos++;
    If[pos > Length[bits], Throw[$OutOfBits, $dbc$tag]];
    pR = condP /. {r_Real :> Rationalize[r]};
    If[bits[[pos]] == 1,
      weight *= pR;
      intervals[[j]] = {lo, pVal};
      True,
      weight *= (1 - pR);
      intervals[[j]] = {pVal, hi};
      False
    ]
  ];

  (* --- Create a fresh interval-tracking token for a new RandomReal[] --- *)
  makeRealVar[] := (
    nReals++;
    AppendTo[intervals, {0, 1}];
    $dbc$irand[nReals, acceptTestI]
  );

  (* --- Sequential Bernoulli decomposition for RandomChoice[w->e] ---
     Element k is chosen with probability weights[[k]] / Total[weights].
     Uses n-1 independent fresh RandomReal[] variables. *)
  seqBernoulli[ws_List, elems_List] :=
    Module[{n = Length[elems], remainW = Total[ws], chosen, var, p},
      chosen = n;  (* default: last element *)
      Do[
        var = makeRealVar[];
        p   = ws[[i]] / remainW;
        If[var < p,
          chosen = i;
          Break[],
          remainW -= ws[[i]]
        ],
        {i, 1, n - 1}
      ];
      elems[[chosen]]
    ];

  (* --- Create a continuous uniform token for RandomReal[{lo,hi}] ---
     Carries seqBernoulli so UpValues can branch into discrete outcomes
     when Floor / Round / Ceiling is applied to the token. *)
  makeContToken[lo_, hi_] := $dbc$contToken["Uniform", lo, hi, seqBernoulli];

  (* --- Run the algorithm with random-call interception --- *)
  result = Catch[
    Block[{
      (* RandomReal[] -> interval-tracking comparison token.
         RandomReal[{lo,hi}] -> continuous uniform token for use with
         Floor/Round/Ceiling (UpValues defined below discretise it). *)
      RandomReal = Function[
        Module[{args = {##}},
          Which[
            args === {},
              makeRealVar[],
            MatchQ[args, {{_, _}}],
              makeContToken[args[[1,1]], args[[1,2]]],
            True,
              Throw[$dbc$cantHandle[
                "RandomReal[" <> ToString[args] <>
                "]: only RandomReal[] and RandomReal[{lo,hi}] are supported"],
                $dbc$tag]]]],
      Random = Function[{}, makeRealVar[]],

      (* RandomInteger: rejection sampling for all ranges.
         Read k = IntegerLength[n-1, 2] bits; if value falls outside
         [0, n-1] throw $dbc$outOfRange so BuildTreeAT silently discards
         the path.  The missing probability fraction is uniform across all
         starting states, so the unnormalised T still satisfies DB exactly. *)
      RandomInteger = Function[
        Module[{args = {##}},
          Which[
            (* RandomInteger[] or RandomInteger[1] -> uniform {0,1} *)
            args === {} || args === {1},
              readBit[],
            (* RandomInteger[{lo,hi}] *)
            MatchQ[args, {{_Integer, _Integer}}],
              Module[{lo = args[[1,1]], hi = args[[1,2]], n, k, val},
                n = hi - lo + 1;
                Which[
                  n == 1, lo,
                  n > 1,
                    k   = IntegerLength[n - 1, 2];
                    val = $dbc$readBitsAsInt[k, readBit];
                    If[val >= n,
                      Throw[$dbc$outOfRange, $dbc$tag],
                      lo + val]
                ]],
            (* RandomInteger[n] -> uniform {0,...,n} *)
            MatchQ[args, {_Integer?NonNegative}],
              Module[{n = args[[1]] + 1, k, val},
                Which[
                  n == 1, 0,
                  n > 1,
                    k   = IntegerLength[n - 1, 2];
                    val = $dbc$readBitsAsInt[k, readBit];
                    If[val >= n,
                      Throw[$dbc$outOfRange, $dbc$tag],
                      val]
                ]],
            True,
              Throw[$dbc$cantHandle[
                "RandomInteger[" <> ToString[args] <> "]: unsupported form"],
                $dbc$tag]
          ]]],

      (* RandomChoice[list]: unweighted uniform choice via bits.
         RandomChoice[weights -> elements]: sequential Bernoulli. *)
      RandomChoice = Function[
        Module[{args = {##}},
          Which[
            (* Unweighted list *)
            Length[args] == 1 && !MatchQ[args[[1]], _Rule],
              Module[{list = args[[1]], n, k, idx},
                n = Length[list];
                Which[
                  n == 0, Throw[$dbc$cantHandle["RandomChoice[]: empty list"], $dbc$tag],
                  n == 1, list[[1]],
                  True,
                    k   = IntegerLength[n - 1, 2];
                    idx = $dbc$readBitsAsInt[k, readBit];
                    If[idx >= n,
                      Throw[$dbc$outOfRange, $dbc$tag],
                      list[[1 + idx]]]
                ]],
            (* Weighted form: RandomChoice[weights -> elements] *)
            Length[args] == 1 && MatchQ[args[[1]], Rule[_List, _List]],
              Module[{ws = args[[1, 1]], elems = args[[1, 2]]},
                If[Length[ws] =!= Length[elems],
                  Throw[$dbc$cantHandle[
                    "RandomChoice[w->e]: weights and elements must have equal length"],
                    $dbc$tag]];
                If[Length[elems] == 0,
                  Throw[$dbc$cantHandle["RandomChoice[w->e]: empty list"], $dbc$tag]];
                If[Length[elems] == 1,
                  elems[[1]],
                  seqBernoulli[ws, elems]]],
            True,
              Throw[$dbc$cantHandle[
                "RandomChoice[...]: unsupported form (only plain list or weights->elements)"],
                $dbc$tag]
          ]]],

      (* RandomVariate: UniformDistribution[{lo,hi}] -> continuous uniform token.
         All other distributions throw $dbc$cantHandle. *)
      RandomVariate = Function[
        Module[{args = {##}},
          Which[
            MatchQ[args, {HoldPattern[UniformDistribution[{_, _}]]}],
              makeContToken[args[[1,1,1]], args[[1,1,2]]],
            True,
              Throw[$dbc$cantHandle[
                "RandomVariate[" <> ToString[args] <>
                "]: only UniformDistribution[{lo,hi}] is supported; " <>
                "use RandomReal[{lo,hi}] for uniform draws"],
                $dbc$tag]]]],

      (* RandomPermutation: Knuth (Fisher-Yates) shuffle via the intercepted
         RandomInteger[], which reads fair bits from the tape. *)
      RandomPermutation = Function[
        Module[{args = {##}, list, n, perm},
          Which[
            MatchQ[args, {_Integer?Positive}],
              n    = args[[1]];
              perm = Range[n];
              Do[With[{j = RandomInteger[{1, i}]},
                   perm[[{i, j}]] = perm[[{j, i}]]],
                 {i, n, 2, -1}];
              perm,
            MatchQ[args, {_List}],
              list = args[[1]];
              n    = Length[list];
              If[n == 0, Return[{}, Module]];
              perm = Range[n];
              Do[With[{j = RandomInteger[{1, i}]},
                   perm[[{i, j}]] = perm[[{j, i}]]],
                 {i, n, 2, -1}];
              list[[perm]],
            True,
              Throw[$dbc$cantHandle[
                "RandomPermutation[" <> ToString[args] <>
                "]: only RandomPermutation[n] and RandomPermutation[list] are supported"],
                $dbc$tag]]]],

      (* RandomSample: partial Knuth shuffle via intercepted RandomInteger[].
         RandomSample[list]    = full random permutation of list.
         RandomSample[list, k] = uniformly random ordered k-subset of list. *)
      RandomSample = Function[
        Module[{args = {##}, list, k, n, perm, result},
          Which[
            MatchQ[args, {_List}],
              list = args[[1]]; n = Length[list];
              If[n == 0, Return[{}, Module]];
              perm = Range[n];
              Do[With[{j = RandomInteger[{1, i}]},
                   perm[[{i, j}]] = perm[[{j, i}]]],
                 {i, n, 2, -1}];
              list[[perm]],
            MatchQ[args, {_List, _Integer?NonNegative}],
              list = args[[1]]; k = args[[2]]; n = Length[list];
              If[k > n,
                Throw[$dbc$cantHandle[
                  "RandomSample: k=" <> ToString[k] <>
                  " exceeds list length " <> ToString[n]], $dbc$tag]];
              perm   = Range[n];
              result = Table[
                With[{j = RandomInteger[{1, i}]},
                  perm[[{i, j}]] = perm[[{j, i}]];
                  list[[perm[[i]]]]],
                {i, n, n - k + 1, -1}];
              result,
            True,
              Throw[$dbc$cantHandle[
                "RandomSample[" <> ToString[args] <>
                "]: only RandomSample[list] and RandomSample[list,k] are supported"],
                $dbc$tag]]]],

      (* Unsupported random functions: throw $dbc$cantHandle immediately. *)
      RandomWord  = Function[Throw[$dbc$cantHandle["RandomWord"],  $dbc$tag]],
      RandomPrime = Function[Throw[$dbc$cantHandle["RandomPrime"], $dbc$tag]]
    },
    alg[state]
    ],
    $dbc$tag,
    Function[{ex}, ex]   (* return thrown value as-is *)
  ];
  Which[
    result === $OutOfBits,  $OutOfBits,
    result === $dbc$outOfRange, $dbc$outOfRange,
    MatchQ[result, $dbc$cantHandle[_]], result,
    (* Detect unconsumed comparison token: RandomReal[] never compared *)
    !FreeQ[{result, weight}, $dbc$irand],
      $dbc$cantHandle[
        "RandomReal[]/Random[] result was used in an unsupported way " <>
        "(not directly in a comparison like RandomReal[] < p)"],
    (* Detect unconsumed continuous token: RandomReal[{lo,hi}] never discretised *)
    !FreeQ[{result, weight}, $dbc$contToken],
      $dbc$cantHandle[
        "RandomReal[{lo,hi}] or RandomVariate[UniformDistribution[...]] result " <>
        "was never discretised with Floor / Round / Ceiling"],
    True,
      {result, weight}
  ]
]


(* ================================================================
   SECTION 1c – RANDOM CALL INTERCEPTION SUPPORT
   ================================================================ *)

(* Helper: read k fair bits and return the integer in {0,...,2^k-1}
   they represent in big-endian binary order.
   Each readBit[] call contributes factor 1/2 to the path weight. *)
$dbc$readBitsAsInt[k_Integer, readBit_] :=
  Fold[#1 * 2 + readBit[] &, 0, Range[k]]

(* $dbc$irand[j, at] is the interval-tracking random token returned
   when RandomReal[]/Random[] is called inside RunWithBitsAT.
   j  = index of this variable in the intervals list
   at = the local acceptTestI function

   UpValues implement comparisons as calls to acceptTestI:
     U < p   ->  at[j, p]           (True if U lands in [lo, p])
     U >= p  ->  !at[j, p]          (True if U lands in [p, hi])
   Symmetric forms handle p on the left side of comparisons.

   Each comparison may narrow the interval [lo, hi] for variable j,
   correctly conditioning on all previous comparisons of the same
   variable (inverse-CDF / sequential Bernoulli semantics). *)
$dbc$irand /: Less[$dbc$irand[j_, at_], p_]         :=  at[j, p]
$dbc$irand /: LessEqual[$dbc$irand[j_, at_], p_]    :=  at[j, p]
$dbc$irand /: Greater[$dbc$irand[j_, at_], p_]      := !at[j, p]
$dbc$irand /: GreaterEqual[$dbc$irand[j_, at_], p_] := !at[j, p]
$dbc$irand /: Less[p_, $dbc$irand[j_, at_]]         := !at[j, p]
$dbc$irand /: LessEqual[p_, $dbc$irand[j_, at_]]    := !at[j, p]
$dbc$irand /: Greater[p_, $dbc$irand[j_, at_]]      :=  at[j, p]
$dbc$irand /: GreaterEqual[p_, $dbc$irand[j_, at_]] :=  at[j, p]

(* Legacy token kept for backward compatibility with old-API callers.
   In the new checker this is never created; any stray occurrence means
   an old algorithm file was loaded. *)
$dbc$rand /: Less[$dbc$rand[at_], p_]         := (at[p] == 1)
$dbc$rand /: LessEqual[$dbc$rand[at_], p_]    := (at[p] == 1)
$dbc$rand /: Greater[$dbc$rand[at_], p_]      := (at[1 - p] == 1)
$dbc$rand /: GreaterEqual[$dbc$rand[at_], p_] := (at[1 - p] == 1)
$dbc$rand /: Less[p_, $dbc$rand[at_]]         := (at[1 - p] == 1)
$dbc$rand /: LessEqual[p_, $dbc$rand[at_]]    := (at[1 - p] == 1)
$dbc$rand /: Greater[p_, $dbc$rand[at_]]      := (at[p] == 1)
$dbc$rand /: GreaterEqual[p_, $dbc$rand[at_]] := (at[p] == 1)


(* ================================================================
   SECTION 1d – CONTINUOUS UNIFORM TOKEN  ($dbc$contToken)
   ================================================================ *)

(* $dbc$contToken["Uniform", lo, hi, sb] represents a continuous latent
   variable X ~ Uniform[lo, hi).  sb is the Module-local seqBernoulli
   function captured at token creation.

   Arithmetic UpValues propagate the distribution through addition;
   Mod collapses a full-period shift back to Uniform[0, L);
   Floor / Round / Ceiling discretise by calling seqBernoulli with the
   exact rational bin probabilities 1/L each -- no approximation.

   Addition:  Uniform[lo, hi) + c  =  Uniform[lo+c, hi+c)
   Mod:       Mod[Uniform[lo, lo+L), L]  =  Uniform[0, L)
              (only when span hi-lo equals the modulus exactly)
   Floor:     Uniform[0, L) -> {0,...,L-1} each with prob 1/L
   Round:     same (for integer L the rounding-boundary effects cancel)
   Ceiling:   Uniform[0, L) -> {1,...,L} each with prob 1/L
              (Mod[Ceiling[X], L] gives {0,...,L-1} for periodic wrap) *)

$dbc$contToken /: Plus[before___, $dbc$contToken["Uniform", lo_, hi_, sb_], after___] :=
  With[{offset = Plus @@ {before, after}},
    If[!FreeQ[offset, $dbc$contToken],
      Throw[$dbc$cantHandle[
        "Adding two continuous random tokens together is not supported"], $dbc$tag],
      $dbc$contToken["Uniform", lo + offset, hi + offset, sb]]]

$dbc$contToken /: Mod[$dbc$contToken["Uniform", lo_, hi_, sb_], L_] :=
  With[{span = hi - lo},
    If[TrueQ[span === L],
      $dbc$contToken["Uniform", 0, L, sb],
      Throw[$dbc$cantHandle[
        "Mod[$dbc$contToken, " <> ToString[L] <>
        "]: interval span " <> ToString[span] <>
        " does not equal the modulus. Only Mod[X,L] where X~Uniform[lo,lo+L) is supported."],
        $dbc$tag]]]

$dbc$contToken /: Floor[$dbc$contToken["Uniform", 0, L_Integer, sb_]] :=
  (* Each integer k in {0,...,L-1} has probability 1/L exactly *)
  sb[Table[1/L, {L}], Range[0, L - 1]]

$dbc$contToken /: Round[$dbc$contToken["Uniform", 0, L_Integer, sb_]] :=
  (* For integer L, Round[Uniform[0,L)] gives each of {0,...,L-1} with prob 1/L:
     boundary bins [0,1/2) and [L-1/2,L) both map to 0 and L-1 respectively
     with combined width 1/L, matching all interior bins exactly. *)
  sb[Table[1/L, {L}], Range[0, L - 1]]

$dbc$contToken /: Ceiling[$dbc$contToken["Uniform", 0, L_Integer, sb_]] :=
  (* Ceiling[Uniform[0,L)] gives {1,...,L} each with prob 1/L.
     Use Mod[Ceiling[X],L]+1 for 1-indexed periodic sites. *)
  sb[Table[1/L, {L}], Range[1, L]]


(* ================================================================
   SECTION 2b – GENERALISED TREE BUILDING  (state discovery)
   ================================================================ *)

(* ----------------------------------------------------------------
   BuildTreeAT
   BFS over bit sequences, starting from seedState.
   New states are discovered automatically as algorithm outputs.
   The algorithm takes ONE argument: alg[state].

   Returns  Association[ state -> { {bits, nextState, pathWeight}, ... } ]
   ---------------------------------------------------------------- *)
Options[BuildTreeAT] = {
  "MaxBitDepth" -> 20,
  "TimeLimit"   -> 60.,
  "Verbose"     -> True
}

BuildTreeAT[seedState_, alg_, OptionsPattern[]] := Module[
  {maxDepth = OptionValue["MaxBitDepth"],
   tlim     = N @ OptionValue["TimeLimit"],
   verbose  = OptionValue["Verbose"],
   discovered, toProcess, result,
   s, queue, bits, res, ns, w, leaves, t0, timedOut},

  discovered = {seedState};
  toProcess  = {seedState};
  result     = <||>;

  While[toProcess =!= {},
    s         = First[toProcess];
    toProcess = Rest[toProcess];
    If[verbose, Print["  Tree for state: ", s]];
    queue    = {{}};
    leaves   = {};
    t0       = AbsoluteTime[];
    timedOut = False;

    While[queue =!= {} && !timedOut,
      If[AbsoluteTime[] - t0 > tlim, timedOut = True; Break[]];
      bits = First[queue]; queue = Rest[queue];
      res  = RunWithBitsAT[alg, s, bits];
      Which[
        res === $OutOfBits && Length[bits] < maxDepth,
          queue = Join[queue, {Append[bits, 0], Append[bits, 1]}],
        res === $OutOfBits,
          Print["  WARNING: MaxBitDepth=", maxDepth,
                " reached at prefix ", bits, " for state ", s,
                " -- path excluded."],
        (* Rejection-sampling dead-end: bit string mapped to out-of-range value.
           Silently discard this path; the missing probability is state-independent
           so the unnormalised transition matrix still satisfies detailed balance. *)
        res === $dbc$outOfRange,
          Null,
        (* Unanalysable call (e.g. RandomVariate, AbsoluteTime) *)
        MatchQ[res, $dbc$cantHandle[_]],
          Print["  ANALYSIS FAILED: algorithm contains a call that cannot be",
                " converted to readBit/acceptTest:"];
          Print["    ", res[[1]]];
          Print["  See the README for supported random-call forms."];
          Return[res, Module],
        True,
          {ns, w} = res;
          (* Guard against unevaluated algorithm calls appearing as states *)
          If[!FreeQ[ns, alg],
            Print["  ANALYSIS FAILED: algorithm returned an unevaluated call as next state."];
            Print["  Check that the algorithm's pattern matches the seed state type."];
            Return[$dbc$cantHandle[
              "Algorithm returned unevaluated call -- pattern mismatch or argument error"],
              Module]
          ];
          AppendTo[leaves, {bits, ns, w}];
          (* Discover new states reached by this path *)
          If[!MemberQ[discovered, ns],
            AppendTo[discovered, ns];
            AppendTo[toProcess, ns]
          ]
      ]
    ];

    If[timedOut, Print["  WARNING: Time limit reached for state ", s]];
    result[s] = leaves
  ];
  result
]

(* ----------------------------------------------------------------
   TreeATToMatrix
   Derive transition matrix from BuildTreeAT output.
   Each leaf contributes its full pathWeight to T[from, to].
   Returns  Association[ {from,to} -> totalProbability ]
   ---------------------------------------------------------------- *)
TreeATToMatrix[treeData_Association] := Module[
  {matrix = <||>},
  Do[
    Do[
      With[{ns = leaf[[2]], w = leaf[[3]]},
        matrix[{s, ns}] = Lookup[matrix, Key[{s, ns}], 0] + w
      ],
      {leaf, treeData[s]}
    ],
    {s, Keys[treeData]}
  ];
  matrix
]


(* ================================================================
   SECTION 3c – GENERALISED CHECKERS  (acceptTest interface)
   ================================================================ *)

(* ----------------------------------------------------------------
   RunNumericalMCMCAT
   Runs a 1-argument algorithm as a genuine Markov chain.
   The global symbol \[Beta] is temporarily assigned numBeta via Block,
   so the same algorithm code works for both symbolic tree-building
   (where \[Beta] is unassigned) and numeric MCMC.
   ---------------------------------------------------------------- *)
Options[RunNumericalMCMCAT] = {
  "NSteps"     -> 100000,
  "WarmupFrac" -> 0.1
}

RunNumericalMCMCAT[allStates_List, alg_, numBeta_, OptionsPattern[]] := Module[
  {nSteps  = OptionValue["NSteps"],
   nWarmup = Round[OptionValue["NSteps"] * OptionValue["WarmupFrac"]],
   state, counts},

  state  = RandomChoice[allStates];
  counts = AssociationThread[allStates -> 0];

  (* \[Beta] is set via Block so MetropolisProb evaluates numerically.
     The algorithm uses native random calls (RandomReal[], RandomInteger[],
     etc.) which are NOT intercepted here -- they run as genuine random calls. *)
  Do[Block[{\[Beta] = numBeta}, state = alg[state]], {nWarmup}];
  Do[
    Block[{\[Beta] = numBeta}, state = alg[state]];
    If[KeyExistsQ[counts, state], counts[state]++],
    {nSteps - nWarmup}];
  counts
]

(* ----------------------------------------------------------------
   BoltzmannWeightsAT
   Compute Boltzmann weights using the unified bare energy function.
   energy[s] must return a real number (no beta).
   ---------------------------------------------------------------- *)
BoltzmannWeightsAT[allStates_List, energy_, numBeta_] := Module[
  {ws, Z},
  ws = N[Exp[-numBeta * energy[#]] & /@ allStates];
  Z  = Total[ws];
  AssociationThread[allStates -> ws / Z]
]

(* ----------------------------------------------------------------
   CheckAlgorithmSafety
   Scans the DownValues of alg for calls that CANNOT be automatically
   converted to readBit/acceptTest during BFS.

   Automatically intercepted (safe to use freely):
     RandomReal[], Random[], RandomInteger[], RandomChoice[]

   Also intercepted (safe to use):
     RandomVariate[UniformDistribution[{lo,hi}]], RandomPermutation[n/list],
     RandomSample[list], RandomSample[list,k]

   These are NOT interceptable and will cause analysis failure:
     RandomWord, RandomPrime, RandomColor, AbsoluteTime, etc.

   Returns True if no unhandleable calls found, False with warnings otherwise.
   ---------------------------------------------------------------- *)
$unanalyzableFunctions = {
  RandomWord, RandomPrime,
  RandomColor, AbsoluteTime, SessionTime, TimeObject, DateObject, Now
};

CheckAlgorithmSafety[alg_Symbol] := Module[
  {defs, found},
  defs  = DownValues[alg];
  If[defs === {},
    Print["  SAFETY WARNING: '", alg, "' has no DownValues. ",
          "Is it defined before CheckAlgorithmSafety is called?"];
    Return[False]
  ];
  found = Select[$unanalyzableFunctions, !FreeQ[defs, #] &];
  If[found === {},
    True,
    Print["  SAFETY FAIL: algorithm '", alg,
          "' contains calls that cannot be converted to readBit/acceptTest: ",
          found];
    Print["  Note: RandomReal[], RandomInteger[], RandomChoice[], ",
          "RandomVariate[UniformDistribution[...]], RandomPermutation[], ",
          "and RandomSample[] ARE supported and are intercepted automatically."];
    False
  ]
]

CheckAlgorithmSafety[alg_] := (
  Print["  SAFETY WARNING: argument is not a named Symbol -- ",
        "cannot inspect DownValues. Pass the function name, not a value."];
  False
)

(* ----------------------------------------------------------------
   CheckEnergySafety
   Scans the DownValues of an energy function for calls that would
   make the symbolic check meaningless: random-number generators,
   time-dependent functions, etc.

   These are safe:   exact arithmetic, lookup tables, Cos, Exp, ...
   These are unsafe: any random call, AbsoluteTime, SessionTime, ...

   Returns True if the energy function looks deterministic.
   ---------------------------------------------------------------- *)
$energyUnsafeFunctions = {
  RandomReal, Random, RandomInteger, RandomChoice, RandomVariate,
  RandomSample, RandomPermutation, RandomWord, RandomPrime, RandomColor,
  AbsoluteTime, SessionTime, TimeObject, DateObject, Now
};

CheckEnergySafety[energy_Symbol] := Module[
  {defs, found},
  defs = DownValues[energy];
  If[defs === {},
    Print["  ENERGY WARNING: '", energy, "' has no DownValues. ",
          "Is it defined before RunFullCheck is called?"];
    Return[False]
  ];
  found = Select[$energyUnsafeFunctions, !FreeQ[defs, #] &];
  If[found === {},
    True,
    Print["  ENERGY SAFETY FAIL: energy function '", energy,
          "' contains calls that make the symbolic check meaningless: ",
          found];
    Print["  Energy must be a deterministic pure function of state. ",
          "Remove all random/time-dependent calls from the energy function."];
    False
  ]
]

CheckEnergySafety[energy_] := True  (* anonymous functions pass through *)


(* ================================================================
   SECTION 3 – CHECKERS
   ================================================================ *)

(* ----------------------------------------------------------------
   CheckCouplingSymmetry
   Verifies that couplingJ[a,b,d2] = couplingJ[b,a,d2] for all pairs
   in typesToCheck at all squared distances in d2Values.
   Uses Simplify (not FullSimplify) — negligible overhead.
   Returns a list of violation records; empty list means symmetric.
   Call this after loading the algorithm file whenever couplingJ is
   defined.  An asymmetric coupling immediately implies a detailed-
   balance violation regardless of the algorithm structure.
   ---------------------------------------------------------------- *)
CheckCouplingSymmetry[typesToCheck_List, d2Values_List] :=
  Module[{violations = {}, diff},
    Do[
      If[a < b,
        Do[
          diff = Simplify[couplingJ[a, b, d2] - couplingJ[b, a, d2]];
          If[diff =!= 0,
            AppendTo[violations,
              <|"a" -> a, "b" -> b, "d2" -> d2, "asymmetry" -> diff|>]],
          {d2, d2Values}]],
      {a, typesToCheck}, {b, typesToCheck}];
    violations]


(* ================================================================
   ERGODICITY CHECK (state-counting)
   ================================================================
   For N labeled particles on S lattice sites (hard-sphere, no
   double-occupation) the theoretical number of states is the
   falling factorial  (S)_N = S*(S-1)*...*(S-N+1).
   S and N are read from any one discovered state; no external
   parameters needed.
   Returns <|"ergodic"->True/False, "found"->k, "theoretical"->m|>.
   ================================================================ *)

CheckErgodicity[allStates_List] :=
  Module[{s0, S, N, theoretical, nFound},
    s0          = First[allStates];
    S           = Length[s0];
    N           = Count[s0, k_ /; k > 0];
    theoretical = Product[S - k, {k, 0, N - 1}];
    nFound      = Length[allStates];
    <|"ergodic"      -> (nFound == theoretical),
      "found"        -> nFound,
      "theoretical"  -> theoretical|>]


(* ================================================================
   FAST EXP-POLYNOMIAL DETAILED BALANCE CHECKER
   ================================================================
   Alternative to FullSimplify.  For Boltzmann algorithms the DB
   expression, after PiecewiseExpand, reduces within each feasible
   case to a sum of terms  c·Exp[-β·L(params)]  where L is a
   polynomial in coupling constants and c is a rational constant.
   Detailed balance holds iff all coefficient sums vanish — a pure
   algebraic check requiring no FullSimplify.

   Falls back to FullSimplify per-expression whenever the structure
   deviates from this pattern.  Never approximates or guesses.

   Toggle in check.wls with FastChecker=1.
   ================================================================ *)

(* Sentinel: signals that fast check cannot determine zero-ness *)
$dbcFS = Symbol["$dbcFallbackSentinel"];

(* ---- Extract {value, condition} pairs from PiecewiseExpand output ---- *)
(* Match on head _Piecewise to avoid triggering argument-validation warnings *)
$dbcPWCases[pw_Piecewise] := Module[{cases, default},
  cases   = pw[[1]];
  default = If[Length[pw] >= 2, pw[[2]], 0];
  Append[cases, {default, True}]]
$dbcPWCases[other_] := {{other, True}}

(* ---- Feasibility: can the case condition hold given β>0, free params?
        Returns True, False, or $dbcFS (timed out / undecidable).       ---- *)
$dbcFeasible[True,  _] := True
$dbcFeasible[False, _] := False
$dbcFeasible[cond_,  symParams_List] := Module[{res},
  res = Quiet @ TimeConstrained[
    Reduce[cond && \[Beta] > 0, Append[symParams, \[Beta]], Reals],
    0.5, $dbcFS];
  Which[res === False, False,
        res === $dbcFS, $dbcFS,
        True, True]]

(* ---- Normalise and merge Exp factors ---- *)
(* Converts (E^a)^n → E^(n·a), then merges products E^a·E^b → E^(a+b),
   distributing over Plus via Expand so that e.g.
   (1 - E^a)/E^a  → E^(-a) - 1  before grouping. *)
$dbcMergeExp[expr_] := FixedPoint[Function[e,
  e
  (* (E^a)^n → E^(n·a) *)
  /. {Power[E^x_, n_]   :> E^Expand[n*x],
      Power[Exp[x_], n_] :> Exp[Expand[n*x]]}
  (* expand products/quotients over sums so Exp factors can be merged *)
  // Expand
  (* merge adjacent Exp factors *)
  /. {Times[a___, E^x_,   E^y_,   b___] :> Times[a, E^Expand[x+y],   b],
      Times[a___, Exp[x_], Exp[y_], b___] :> Times[a, Exp[Expand[x+y]], b]}],
  expr, 30]

(* ---- Split one expanded term into {coefficient, exponentPoly}
        where term = coefficient · Exp[exponentPoly].
        Returns $dbcFS if structure is unexpected.                       ---- *)
$dbcSplitTerm[0]   := {0, 0}
$dbcSplitTerm[0.]  := {0, 0}
$dbcSplitTerm[t_]  := Module[{ep},
  ep = Cases[{t}, Exp[x_] :> x, {0, Infinity}];
  Which[
    Length[ep] == 0, {t, 0},
    Length[ep] == 1, {Expand[Cancel[t / Exp[ep[[1]]]]], Expand[ep[[1]]]},
    True, $dbcFS]]

(* ---- Group split-term list by exponent polynomial (exact comparison).
        Returns list of {exponentPoly, {coeff1, coeff2, ...}}.          ---- *)
$dbcGroupByExp[splits_List] := Module[
  {polys = {}, groups = {}, matched},
  Scan[Function[s,
    matched = SelectFirst[Range @ Length[polys],
      Function[i, Expand[s[[2]] - polys[[i]]] === 0], 0];
    If[matched === 0,
      AppendTo[polys,  s[[2]]];
      AppendTo[groups, {s[[1]]}],
      groups[[matched]] = Append[groups[[matched]], s[[1]]]]],
    splits];
  MapThread[{#1, #2} &, {polys, groups}]]

(* ---- Check whether an expression is identically zero by treating
        each Exp[poly] as an independent basis element.
        Returns True, False, or $dbcFS.                                  ---- *)
$dbcIsExpZero[0]   := True
$dbcIsExpZero[0.]  := True
$dbcIsExpZero[expr_] := Module[
  {merged, expanded, terms, splits, grouped},
  merged   = $dbcMergeExp[expr];
  expanded = Expand[merged];
  If[expanded === 0, Return[True]];
  terms  = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  splits = Map[$dbcSplitTerm, terms];
  If[MemberQ[splits, $dbcFS], Return[$dbcFS]];
  grouped = $dbcGroupByExp[splits];
  If[AllTrue[grouped, Function[g, Expand[Total[g[[2]]]] === 0]],
    True, False]]

(* ---- Substitute equality constraints from a condition into an expression.
        E.g. cond = (Jpair12 == 0) && (x > 0)  →  val /. {Jpair12 -> 0}  ---- *)
$dbcSubstEqualities[val_, True]  := val
$dbcSubstEqualities[val_, False] := val
$dbcSubstEqualities[val_, cond_] := Module[
  {eqs, rules, solns},
  eqs = Cases[{cond}, HoldPattern[lhs_ == rhs_], Infinity];
  If[eqs === {}, Return[val]];
  (* Solve each equality for a symbolic param and substitute *)
  rules = Flatten @ Map[
    Function[eq,
      Quiet @ Check[Solve[eq, {}], {}]],  (* empty var list: Solve returns rules *)
    eqs];
  If[rules === {}, val /. Map[(#[[1]] -> #[[2]]) &, eqs],
     val /. rules]]

(* ---- Fast check for one DB expression.
        Returns True (zero), False (non-zero), or $dbcFS (fall back).   ---- *)
$dbcCheckOneExpr[expr_, assm_List, symParams_List] := Module[
  {pw, cases},
  pw    = PiecewiseExpand[expr, assm];
  cases = $dbcPWCases[pw];
  Catch[
    Scan[Function[vc,
      With[{val = vc[[1]], cond = vc[[2]]},
        With[{feas = $dbcFeasible[cond, symParams]},
          If[feas === False,   Return[]];          (* infeasible case: skip *)
          If[feas === $dbcFS,  Throw[$dbcFS]];     (* can't decide: fall back *)
          (* Substitute any equality constraints before the exponent check.
             E.g. cond = (Jpair12 == 0) → substitute Jpair12→0 into val. *)
          With[{val2 = $dbcSubstEqualities[val, cond]},
            With[{z = $dbcIsExpZero[val2]},
              Which[z === True,    Return[],          (* zero: continue *)
                    z === $dbcFS,  Throw[$dbcFS],     (* can't decide: fall back *)
                    True,          Throw[False]]]]]]],  (* non-zero found *)
      cases];
    True]]  (* every feasible case was zero *)

(* ================================================================
   CheckDetailedBalanceFast
   Drop-in replacement for CheckDetailedBalance.
   Uses Exp-coefficient matching; falls back to FullSimplify
   per-expression only when needed.  Never approximates.
   Same argument signature and return format as CheckDetailedBalance.

   Parallelises using ParallelMap when subkernels are available:
   expressions are built on the main kernel, then each pair's
   fast check runs on a subkernel.  Requires the helper functions
   ($dbcCheckOneExpr etc.) to be distributed first via
   $dbcDistributeFastChecker[].
   ================================================================ *)

(* ---- Deduplicate a list by Hash[].
        Returns {uniqueIdxs, canonIdx} where:
          uniqueIdxs : indices (into lst) of the first occurrence of each
                       distinct hash — the canonical representatives.
          canonIdx   : for each k, the canonical index canonIdx[[k]] gives
                       the representative for lst[[k]].
        Safe: hash collisions map to the same representative; the checker
        still sees the expression from that representative (syntactically
        identical for the colliding pair), so no false negatives.        ---- *)
$dbcDedup[lst_List] := Module[{seen = <||>, uIdxs = {}, cIdx = {}},
  Do[
    With[{h = Hash[lst[[k]]]},
      If[KeyExistsQ[seen, h],
        AppendTo[cIdx, seen[h]],
        seen[h] = k;
        AppendTo[uIdxs, k];
        AppendTo[cIdx, k]]],
    {k, Length[lst]}];
  {uIdxs, cIdx}]

(* Distribute all fast-checker helpers to subkernels.
   Called once from check.wls after LaunchKernels[].           *)
$dbcDistributeFastChecker[] := (
  DistributeDefinitions[
    $dbcFS,
    $dbcPWCases,
    $dbcFeasible,
    $dbcMergeExp,
    $dbcSplitTerm,
    $dbcGroupByExp,
    $dbcIsExpZero,
    $dbcSubstEqualities,
    $dbcCheckOneExpr])

(* Worker: evaluate one (expr, assm, symParams) triple.
   Returns {fastResult, expr, assm} so the caller can fall back if needed. *)
$dbcFastWorker[expr_, assm_, symParams_] :=
  {$dbcCheckOneExpr[expr, assm, symParams], expr, assm}

CheckDetailedBalanceFast[matrix_Association, allStates_List, symEnergy_,
                         extraAssumptions_List : {}] := Module[
  {n, pairs, assm, symParams, exprs,
   nonTrivIdx, ntExprs, ntPairs,
   uniqueIdxs, canonIdx, uniqueExprs, uniqueResults, idxToResult,
   results, violations, fsCache, nK,
   si, sj, tij, tji, ei, ej, fRes},

  n         = Length[allStates];
  pairs     = Flatten[Table[{i, j}, {i, 1, n}, {j, i+1, n}], 1];
  If[Length[pairs] == 0, Return[{}]];

  assm      = Join[{\[Beta] > 0}, extraAssumptions];
  symParams = Cases[extraAssumptions, Element[x_, Reals] :> x, Infinity];

  (* Build all expressions on main kernel — symEnergy may not be on subkernels *)
  exprs = Map[
    Function[pair,
      si  = allStates[[pair[[1]]]]; sj = allStates[[pair[[2]]]];
      tij = Lookup[matrix, Key[{si, sj}], 0];
      tji = Lookup[matrix, Key[{sj, si}], 0];
      ei  = symEnergy[si] /. r_Real :> Rationalize[r];
      ej  = symEnergy[sj] /. r_Real :> Rationalize[r];
      tij * Exp[-\[Beta] * ei] - tji * Exp[-\[Beta] * ej]],
    pairs];

  (* A: Drop pairs whose DB expression is syntactically 0 (both tij=tji=0). *)
  nonTrivIdx = Select[Range[Length[exprs]], exprs[[#]] =!= 0 &];
  If[Length[nonTrivIdx] == 0, Return[{}]];
  ntExprs = exprs[[nonTrivIdx]];
  ntPairs = pairs[[nonTrivIdx]];

  (* Dedup: run each unique expression once; broadcast results back. *)
  {uniqueIdxs, canonIdx} = $dbcDedup[ntExprs];
  uniqueExprs = ntExprs[[uniqueIdxs]];

  (* B+C: ParallelMap for natural load balancing; sequential when the workload
          is too small to justify the dispatch overhead. *)
  nK = Length[Kernels[]];
  With[{a = assm, sp = symParams},
    uniqueResults = If[nK > 0 && Length[uniqueExprs] > nK,
      ParallelMap[$dbcFastWorker[#, a, sp] &, uniqueExprs],
      Map[       $dbcFastWorker[#, a, sp] &, uniqueExprs]]];

  idxToResult = AssociationThread[uniqueIdxs -> uniqueResults];
  results = Map[Function[k, idxToResult[canonIdx[[k]]]], Range[Length[ntExprs]]];

  (* Collect violations; fall back to FullSimplify once per unique expression. *)
  fsCache    = <||>;
  violations = {};
  Do[
    With[{fRes = results[[k, 1]],
          pair = ntPairs[[k]],
          cidx = canonIdx[[k]]},
      Which[
        fRes === True,
          Null,
        True,
          If[!KeyExistsQ[fsCache, cidx],
            fsCache[cidx] = FullSimplify[
              PiecewiseExpand[ntExprs[[cidx]]], Assumptions -> assm]];
          If[fsCache[cidx] =!= 0,
            AppendTo[violations,
              <|"pair"     -> {allStates[[pair[[1]]]], allStates[[pair[[2]]]]},
                "residual" -> fsCache[cidx]|>]]]],
    {k, Length[ntExprs]}];

  violations
]

(* ----------------------------------------------------------------
   CheckDetailedBalance
   Verifies T(i->j)*pi(i) = T(j->i)*pi(j) for all i<j pairs,
   where pi(s) = Exp[-beta * symEnergy[s]].
   Uses FullSimplify with beta>0 to resolve Piecewise expressions.
   Returns list of violation records; empty list = PASS.
   ---------------------------------------------------------------- *)
CheckDetailedBalance[matrix_Association, allStates_List, symEnergy_,
                     extraAssumptions_List : {}] := Module[
  {n = Length[allStates], pairs, exprs,
   nonTrivIdx, ntExprs, ntPairs,
   uniqueIdxs, canonIdx, uniqueExprs, uniqueResults, idxToResult,
   results, violations, nK},

  pairs = Flatten[Table[{i, j}, {i, 1, n}, {j, i + 1, n}], 1];
  If[Length[pairs] == 0, Return[{}]];

  (* Build all expressions on the MAIN kernel — symEnergy is only called here.
     The results are pure symbolic data safe to send to remote kernels. *)
  exprs = Map[
    Function[ij,
      With[{si = allStates[[ij[[1]]]], sj = allStates[[ij[[2]]]]},
        With[{tij = Lookup[matrix, Key[{si, sj}], 0],
              tji = Lookup[matrix, Key[{sj, si}], 0],
              ei  = symEnergy[si] /. {r_Real :> Rationalize[r]},
              ej  = symEnergy[sj] /. {r_Real :> Rationalize[r]}},
          tij * Exp[-\[Beta] * ei] - tji * Exp[-\[Beta] * ej]]]],
    pairs];

  (* A: Drop pairs whose DB expression is syntactically 0 (both tij=tji=0). *)
  nonTrivIdx = Select[Range[Length[exprs]], exprs[[#]] =!= 0 &];
  If[Length[nonTrivIdx] == 0, Return[{}]];
  ntExprs = exprs[[nonTrivIdx]];
  ntPairs = pairs[[nonTrivIdx]];

  (* Dedup: run FullSimplify only on unique expressions. *)
  {uniqueIdxs, canonIdx} = $dbcDedup[ntExprs];
  uniqueExprs = ntExprs[[uniqueIdxs]];

  (* B+C: ParallelMap for natural per-expression load balancing; sequential
          when the workload is too small to justify the dispatch overhead. *)
  nK = Length[Kernels[]];
  With[{assm = Join[{\[Beta] > 0}, extraAssumptions]},
    uniqueResults = If[nK > 0 && Length[uniqueExprs] > nK,
      ParallelMap[FullSimplify[PiecewiseExpand[#], Assumptions -> assm] &, uniqueExprs],
      Map[        FullSimplify[PiecewiseExpand[#], Assumptions -> assm] &, uniqueExprs]]];

  idxToResult = AssociationThread[uniqueIdxs -> uniqueResults];
  results = Map[Function[k, idxToResult[canonIdx[[k]]]], Range[Length[ntExprs]]];

  violations = {};
  Do[
    If[results[[k]] =!= 0,
      AppendTo[violations,
        <|"pair"     -> {allStates[[ntPairs[[k, 1]]]], allStates[[ntPairs[[k, 2]]]]},
          "residual" -> results[[k]]|>]],
    {k, Length[ntExprs]}];
  violations
]


(* ----------------------------------------------------------------
   RunNumericalMCMC
   Run numAlg with true random bits; sample outcomes weighted by
   their returned probabilities.
   Returns  Association[ state -> visitCount ]
   ---------------------------------------------------------------- *)
Options[RunNumericalMCMC] = {
  "NSteps"     -> 100000,
  "WarmupFrac" -> 0.1
}

RunNumericalMCMC[allStates_List, numAlg_, OptionsPattern[]] := Module[
  {nSteps  = OptionValue["NSteps"],
   nWarmup = Round[OptionValue["NSteps"] * OptionValue["WarmupFrac"]],
   state, counts},

  state  = RandomChoice[allStates];
  counts = AssociationThread[allStates -> 0];

  mcmcStep[] := Module[{liveRb, raw, outs, u, cumP, ns},
    liveRb[] := RandomInteger[1];
    raw  = numAlg[state, liveRb];
    outs = If[ListQ[raw] && Length[raw] > 0 && ListQ[raw[[1]]],
               raw, {{1, raw}}];
    u = RandomReal[]; cumP = 0.; ns = outs[[-1, 2]];
    Do[cumP += N[out[[1]]]; If[u < cumP, ns = out[[2]]; Break[]], {out, outs}];
    state = ns
  ];

  Do[mcmcStep[], {nWarmup}];
  Do[mcmcStep[];
     If[KeyExistsQ[counts, state], counts[state]++,
        Print["  WARNING: unexpected state from numAlg: ", state]],
     {nSteps - nWarmup}];
  counts
]

(* ----------------------------------------------------------------
   BoltzmannWeights
   numEnergy[s] must be fully numeric WITH beta already included.
   Returns  Association[ state -> normalisedWeight ]
   ---------------------------------------------------------------- *)
BoltzmannWeights[allStates_List, numEnergy_] := Module[
  {ws, Z},
  ws = N[Exp[-numEnergy[#]] & /@ allStates];
  Z  = Total[ws];
  AssociationThread[allStates -> ws / Z]
]


(* ================================================================
   SECTION 3b – JSON EXPORT + PYTHON VISUALISATION
   ================================================================ *)

(* Compact string for a probability (symbolic or numeric) *)
$probStr[p_] := Which[
  p === 0 || p === 0.,   "0",
  p === 1 || p === 1.,   "1",
  NumericQ[p],           ToString[NumberForm[N[p], {4, 3}]],
  True,                  StringTake[ToString[p, InputForm], UpTo[60]]
]

(* Convert a tree leaf to the canonical JSON format.
   Handles both old-API leaves {bits, {{p1,s1},...}} and
   new-API leaves {bits, nextState, pathWeight}. *)
$leafToJSON[{bits_List, outcomes_List}] :=
  <|"bits"     -> bits,
    "outcomes" -> Table[<|"probStr" -> $probStr[out[[1]]], "state" -> out[[2]]|>,
                        {out, outcomes}]|>
$leafToJSON[{bits_List, ns_, w_}] :=
  <|"bits"     -> bits,
    "outcomes" -> {<|"probStr" -> $probStr[w], "state" -> ns|>}|>

(* Normalise a leaf to old-style {bits, {{p,s},...}} for Mathematica renderers *)
$normLeaf[{bits_List, ns_, w_}] := {bits, {{w, ns}}}
$normLeaf[leaf_List]            := leaf

ExportReportJSON[args_Association, outPath_String] := Module[
  {name, allStates, treeData, matrix, violations, simFreq, bw, kl, algCode,
   pass, n, dbJSON, treeJSON, matJSON},

  name       = args["name"];
  allStates  = args["allStates"];
  treeData   = args["treeData"];
  matrix     = args["matrix"];
  violations = args["violations"];
  simFreq    = args["simFreq"];
  bw         = args["bw"];
  kl         = N @ args["kl"];
  algCode    = args["algCode"];
  (* Normalise violations: must be a plain list (guards against Null / Missing) *)
  If[!ListQ[violations], violations = {}];
  n = Length[allStates];

  (* Tree data as an ordered list (same order as allStates) so Python
     can look up by index rather than by state key string, which avoids
     mismatches between Mathematica's ToString and Python's str() for
     rational and list-valued states. *)
  treeJSON = Table[($leafToJSON /@ treeData[s]), {s, allStates}];

  matJSON = Flatten[Table[
    <|"from" -> allStates[[i]],
      "to"   -> allStates[[j]],
      "str"  -> $probStr[Lookup[matrix, Key[{allStates[[i]], allStates[[j]]}], 0]]|>,
    {i, n}, {j, n}], 1];

  (* Use AnyTrue with an explicit predicate so pair matching works correctly
     regardless of how violations is structured internally. *)
  dbJSON = Flatten[Table[
    With[{si = allStates[[i]], sj = allStates[[j]]},
      <|"i" -> si, "j" -> sj,
        "pass" -> !AnyTrue[violations,
                    AssociationQ[#] && #["pair"] === {si, sj} &]|>],
    {i, n}, {j, i+1, n}], 1];

  (* Derive the overall pass verdict from the db pair results so that the
     banner and the table always agree. *)
  pass = AllTrue[dbJSON, TrueQ @ #["pass"] &];

  Export[outPath,
    <|"name"      -> name,
      "pass"      -> pass,
      "kl"        -> kl,
      "algCode"   -> algCode,
      "allStates" -> allStates,
      "simFreq"   -> (N[simFreq[#]] & /@ allStates),
      "boltzmann" -> (N[bw[#]]      & /@ allStates),
      "dbPairs"   -> dbJSON,
      "matrix"    -> matJSON,
      "treeData"  -> treeJSON|>,
    "JSON"]
]

ExportAndShowPython[args_Association] := Module[
  {safeName, jsonPath, pyScript, result, pngPath},
  safeName = StringReplace[args["name"],
               {" " -> "_", Except[WordCharacter | "_"] -> ""}];
  jsonPath = $dbcDir <> safeName <> "_report.json";
  pyScript = $dbcDir <> "show_report.py";

  If[!FileExistsQ[pyScript],
    Print["  show_report.py not found at: ", pyScript]; Return[None]];

  Print["  Exporting JSON ..."];
  ExportReportJSON[args, jsonPath];

  Print["  Running Python visualisation ..."];
  result = RunProcess[{"python3", pyScript, jsonPath}];

  If[result["ExitCode"] == 0,
    pngPath = StringTrim[result["StandardOutput"]];
    Print["  Report saved: ", pngPath];
    If[FileExistsQ[pngPath], RunProcess[{"open", pngPath}]],
    Print["  Python error:\n", result["StandardError"]]
  ]
]

(* ----------------------------------------------------------------
   ExportReportNotebook
   Exports the full check report as a native Mathematica notebook
   (.nb file) and opens it in Mathematica.app.

   All mathematical expressions are rendered as Mathematica's own
   typeset output -- inherently vector, perfectly sharp at any zoom.
   Decision-tree outcome nodes show only the destination state; hover
   over any green/orange node to see the exact path weight as a Tooltip.
   Trees are arranged in rows of up to 3, so they never overflow
   horizontally and cannot overlap each other or subsequent sections.
   ---------------------------------------------------------------- *)
ExportReportNotebook[args_Association] := Module[
  {name, allStates, treeData, matrix, violations, simFreq, bw, kl, algCode,
   showTrees, showT, doNumerical,
   pass, badge, treeGrid, matGrid, dbTab, freqPanel, cells, nb,
   safeName, nbPath},

  name        = args["name"];
  allStates   = args["allStates"];
  treeData    = args["treeData"];
  matrix      = args["matrix"];
  violations  = If[ListQ[args["violations"]], args["violations"], {}];
  simFreq     = args["simFreq"];
  bw          = args["bw"];
  kl          = args["kl"];
  algCode     = args["algCode"];
  showTrees   = Lookup[args, "showTrees",   True];
  showT       = Lookup[args, "showT",       True];
  doNumerical = Lookup[args, "doNumerical", True];
  pass        = violations === {};

  badge = Panel[
    Row[{
      Style[name, 18, Bold, GrayLevel[0.1]],
      Spacer[25],
      Style[If[pass,
               "\[FilledCircle]  DETAILED BALANCE:  PASS",
               "\[FilledCircle]  DETAILED BALANCE:  FAIL"],
            17, Bold,
            If[pass, RGBColor[0.04,0.54,0.04], RGBColor[0.74,0.04,0.04]]]
    }],
    Background -> GrayLevel[0.94],
    FrameMargins -> {{16,16},{12,12}}];

  (* Only build the sections that are requested *)
  If[showTrees,
    treeGrid = Grid[
      Partition[
        Table[
          Framed[
            Column[{DrawStateTree[s, $normLeaf /@ treeData[s]]}, Alignment -> Center],
            FrameStyle -> GrayLevel[0.82], RoundingRadius -> 5,
            Background -> GrayLevel[0.975], FrameMargins -> 8],
          {s, allStates}],
        UpTo[3]],
      Spacings -> {2, 2}, Alignment -> {Left, Top}]];

  If[showT,    matGrid   = MakeTransitionGrid[allStates, matrix]];
               dbTab     = MakeDBTable[allStates, violations];
  If[doNumerical, freqPanel = MakeFrequencyPanel[allStates, simFreq, bw, kl]];

  cells = Flatten @ {
    Cell[BoxData @ ToBoxes @ badge,
         "Output", CellMargins -> {{8,8},{4,14}}],

    Cell["Algorithm Under Test", "Section"],
    Cell[algCode, "Code"],

    If[showTrees, {
      Cell["Decision Trees  (Symbolic Execution Paths)", "Section"],
      Cell[TextData[{
        "Each tree shows all bit sequences the algorithm can consume from a \
given starting state.  ",
        StyleBox["Blue", FontWeight->"Bold",
                 FontColor->RGBColor[0.22,0.42,0.70]],
        " = root / internal bit-choice node.  ",
        StyleBox["Green", FontWeight->"Bold", FontColor->Darker[Green,0.2]],
        " outcome = moved to a new state.  ",
        StyleBox["Orange", FontWeight->"Bold", FontColor->Darker[Orange,0.1]],
        " outcome = stayed.  Edge labels are bit values (0/1).  \
Hover over an outcome node to see the exact path weight."
      }], "Text"],
      Cell[BoxData @ ToBoxes @ treeGrid,
           "Output", CellMargins -> {{8,8},{4,4}}]
    }, Nothing],

    If[showT, {
      Cell["Symbolic Transition Matrix  T[ i \[Rule] j ]", "Section"],
      Cell["Entry (i, j) = total probability of transitioning FROM state i \
TO state j, accumulated over all execution paths.  \[Beta] is kept \
symbolic throughout.", "Text"],
      Cell[BoxData @ ToBoxes @ matGrid,
           "Output", CellMargins -> {{8,8},{4,4}}]
    }, Nothing],

    {Cell["Detailed Balance Check:  \
T(i\[Rule]j)\[CenterDot]\[Pi](i) = T(j\[Rule]i)\[CenterDot]\[Pi](j)",
           "Section"],
    Cell["Every pair of distinct states is tested using FullSimplify with \
\[Beta] > 0.  \[Pi](s) = Exp[\[Minus]\[Beta] E(s)].", "Text"],
    Cell[BoxData @ ToBoxes @ dbTab,
         "Output", CellMargins -> {{8,8},{4,4}}]},

    If[doNumerical, {
      Cell["Numerical MCMC Validation", "Section"],
      Cell["The algorithm is run as a live Markov chain with genuinely random \
bits.  Simulated state frequencies are compared to the analytical Boltzmann \
distribution Exp[\[Minus]\[Beta]E(s)] / Z.", "Text"],
      Cell[BoxData @ ToBoxes @ freqPanel,
           "Output", CellMargins -> {{8,8},{4,18}}]
    }, Nothing]
  };

  nb = Notebook[cells,
    WindowTitle  -> "DetailedBalanceChecker \[LongDash] " <> name,
    StyleDefinitions -> "Default.nb",
    WindowSize   -> {1200, 900},
    WindowMargins -> {{Automatic, Automatic}, {Automatic, 0}},
    Editable     -> True,
    Background   -> White
  ];

  safeName = StringReplace[name,
               {" " -> "_", Except[WordCharacter | "_"] -> ""}];
  nbPath   = $dbcDir <> safeName <> "_report.nb";

  Print["  Exporting notebook: ", nbPath];
  Export[nbPath, nb];
  RunProcess[{"open", nbPath}];
  nbPath
]


(* ================================================================
   SECTION 4 – VISUALISATION PRIMITIVES
   ================================================================ *)

(* ----------------------------------------------------------------
   DrawStateTree
   Render the BFS decision tree for one starting state as a Graph.
   leaves = { {bits_List, outcomes}, ... }
   outcomes = {{p1,s1},{p2,s2},...}
   Green outcome nodes  = algorithm moved to a NEW state.
   Orange outcome nodes = algorithm stayed in the SAME state.

   Each bit edge in the tree represents one bit consumed from the tape.
   For RandomReal[] comparisons this bit encodes whether the latent
   U variable landed in the lower (bit=1, accept) or upper (bit=0,
   reject) sub-interval after conditioning; the edge weight shown is
   the conditional probability for that sub-interval.
   ---------------------------------------------------------------- *)

(* String vertex IDs -- safer than lists as Graph vertex names *)
$bKey[b_List] := If[b === {}, "root", "b" <> StringJoin[ToString /@ b]]
$oKey[b_List, i_Integer] := $bKey[b] <> "o" <> ToString[i]

(* Compact probability label *)
$probLabel[p_] := Which[
  p === 1 || p === 1.,   Style["p=1", 7, Darker[Green,0.3]],
  p === 0 || p === 0.,   Style["p=0", 7, Red],
  NumericQ[p],           Style["p=" <> ToString[NumberForm[N[p],{3,2}]], 7, GrayLevel[0.3]],
  True,                  Style[TraditionalForm[p], 7, GrayLevel[0.2]]
]

DrawStateTree[startState_, leaves_List] := Module[
  {leafBits, allPfx, intPfx, leafPfx,
   bitEdges, outEdges, allEdges,
   vLabels, eLabels, vStyle, vSize, nLeaves, maxDepth, g},

  If[leaves === {},
    Return @ Framed[Style["(no paths)", 9, Gray],
                   FrameStyle -> LightGray, ImageSize -> 140]
  ];

  leafBits = #[[1]] & /@ leaves;
  nLeaves  = Length[leafBits];
  maxDepth = Max[Length /@ leafBits];

  (* All bit-sequence prefixes appearing as nodes *)
  allPfx = DeleteDuplicates @ Flatten[
    Table[Take[b, k], {b, leafBits}, {k, 0, Length[b]}], 1];
  intPfx = DeleteDuplicates @ Flatten[
    Table[Take[b, k], {b, leafBits}, {k, 0, Length[b]-1}], 1];
  leafPfx = Complement[allPfx, intPfx];

  (* Edges between bit nodes *)
  bitEdges = DeleteDuplicates @ Flatten[
    Table[DirectedEdge[$bKey @ Take[b,k-1], $bKey @ Take[b,k]],
          {b, leafBits}, {k, 1, Length[b]}], 1];

  (* Edges to outcome nodes *)
  outEdges = Flatten @ Table[
    Table[DirectedEdge[$bKey[leaf[[1]]], $oKey[leaf[[1]],i]],
          {i, Length[leaf[[2]]]}],
    {leaf, leaves}];

  allEdges = Join[bitEdges, outEdges];

  (* Vertex labels *)
  vLabels = Flatten @ {
    $bKey[{}] -> Placed[Style["S=" <> ToString[startState], 9, Bold, White], Center],
    Table[$bKey[b] -> Placed[Style["?", 8, White], Center],
          {b, Complement[intPfx, {{}}]}],
    Table[$bKey[b] -> Placed[Style["\[DownArrow]", 9, GrayLevel[0.35]], Center],
          {b, leafPfx}],
    Flatten @ Table[
      With[{bits = leaf[[1]], outs = leaf[[2]]},
        Table[$oKey[bits,i] -> Placed[
          (* Show only the destination state inside the node.
             Hovering reveals the exact path weight as a Tooltip. *)
          Tooltip[
            Style["\[RightArrow]" <> ToString[outs[[i,2]]], 8,
                  If[outs[[i,2]] =!= startState, Darker[Green,0.2], Darker[Orange,0.15]]],
            Column[{
              Style["Next state:", 9, Bold, GrayLevel[0.2]],
              Style[ToString[outs[[i,2]]], 9, GrayLevel[0.1]],
              Spacer[4],
              Style["Path weight:", 9, Bold, GrayLevel[0.2]],
              Style[TraditionalForm[outs[[i,1]]], 10]
            }, Spacings -> 0.3]
          ],
          Center],
        {i, Length[outs]}]
      ],
      {leaf, leaves}]
  };

  (* Edge labels: 0 / 1 on bit edges *)
  eLabels = DeleteDuplicates @ Flatten[
    Table[DirectedEdge[$bKey @ Take[b,k-1], $bKey @ Take[b,k]] ->
            Placed[Style[ToString[b[[k]]], 9, Bold, RGBColor[0.2,0.3,0.7]], Automatic],
          {b, leafBits}, {k, 1, Length[b]}], 1];

  (* Vertex colours *)
  vStyle = Flatten @ {
    $bKey[{}] -> RGBColor[0.22,0.42,0.70],
    Table[$bKey[b] -> GrayLevel[0.52], {b, Complement[intPfx, {{}}]}],
    Table[$bKey[b] -> GrayLevel[0.70], {b, leafPfx}],
    Flatten @ Table[
      With[{bits = leaf[[1]], outs = leaf[[2]]},
        Table[$oKey[bits,i] ->
              If[outs[[i,2]] =!= startState,
                 RGBColor[0.25,0.68,0.38],
                 RGBColor[0.88,0.58,0.20]],
              {i, Length[outs]}]],
      {leaf, leaves}]
  };

  (* Vertex sizes *)
  vSize = Flatten @ {
    $bKey[{}] -> 0.65,
    Table[$bKey[b] -> 0.40, {b, Complement[intPfx, {{}}]}],
    Table[$bKey[b] -> 0.35, {b, leafPfx}],
    Flatten @ Table[Table[$oKey[leaf[[1]],i] -> 0.60, {i, Length[leaf[[2]]]}],
                    {leaf, leaves}]
  };

  Graph[
    allEdges,
    VertexLabels -> vLabels,
    EdgeLabels   -> eLabels,
    VertexStyle  -> vStyle,
    VertexSize   -> vSize,
    GraphLayout  -> {"LayeredDigraphEmbedding",
                     "RootVertex"   -> $bKey[{}],
                     "Orientation"  -> Top},
    ImageSize    -> {Max[300, 200*nLeaves], Max[240, 120*(maxDepth+2)]},
    Background   -> GrayLevel[0.97],
    PlotLabel    -> Style["State " <> ToString[startState], 10, Bold,
                          GrayLevel[0.3]]
  ]
]

(* ----------------------------------------------------------------
   MakeTransitionGrid
   Styled Grid showing the symbolic transition matrix.
   ---------------------------------------------------------------- *)
MakeTransitionGrid[allStates_List, matrix_Association] := Module[
  {n = Length[allStates], hdr, rows},
  hdr = Prepend[
    Style[ToString[#], Bold, 10, GrayLevel[0.2]] & /@ allStates,
    Style["T[i\[Rule]j]", 10, Italic, GrayLevel[0.4]]];
  rows = Table[
    Prepend[
      Table[
        With[{p = Lookup[matrix, Key[{allStates[[i]], allStates[[j]]}], 0]},
          If[p === 0,
             Style["0", 9, GrayLevel[0.7]],
             (* Pane with ShrinkToFit scales long expressions to fit the cell
                while preserving full vector resolution for zooming. *)
             Pane[
               Style[TraditionalForm @ FullSimplify[p], 11],
               {200, Automatic},
               ImageSizeAction -> "ShrinkToFit"
             ]]],
        {j, n}],
      Style[ToString[allStates[[i]]], Bold, 10, GrayLevel[0.2]]],
    {i, n}];
  Grid[
    Prepend[rows, hdr],
    Frame      -> All,
    FrameStyle -> GrayLevel[0.82],
    Background -> {None, None, Flatten @ {
      Table[{1,j} -> RGBColor[0.84,0.90,1.00], {j,n+1}],
      Table[{i,1} -> RGBColor[0.84,0.90,1.00], {i,n+1}],
      Table[{i+1,i+1} -> RGBColor[1.00,0.98,0.84], {i,n}]}},
    Spacings   -> {2, 1},
    Alignment  -> Center
  ]
]

(* ----------------------------------------------------------------
   MakeDBTable
   Colour-coded table of detailed-balance pair results.
   ---------------------------------------------------------------- *)
MakeDBTable[allStates_List, violations_] := Module[
  {viols, n, pairs, hdr, rows},
  (* Normalise: guard against Null / Missing from the new API path *)
  viols  = If[ListQ[violations], violations, {}];
  n      = Length[allStates];
  pairs  = Flatten[Table[{allStates[[i]],allStates[[j]]},
                         {i,n},{j,i+1,n}], 1];
  hdr = Style[#, Bold, 10] & /@
        {"State i", "State j",
         "T(i\[Rule]j)\[CenterDot]\[Pi](i)  \[Minus]  T(j\[Rule]i)\[CenterDot]\[Pi](j)",
         "Result"};
  rows = Table[
    With[{pass = !AnyTrue[viols, AssociationQ[#] && #["pair"] === pair &]},
      {Style[ToString[pair[[1]]], 10],
       Style[ToString[pair[[2]]], 10],
       Style[If[pass, "= 0   (FullSimplify, \[Beta] > 0)",
                      ToString[TraditionalForm @
                        First[Select[viols, #["pair"]===pair&],
                              <|"residual"->"?"|>]["residual"]]],
             9, If[pass, Darker[Green,0.2], Darker[Red,0.1]]],
       Style[If[pass, "\[Checkmark] PASS", "\[Times] FAIL"],
             10, Bold, If[pass, Darker[Green], Red]]}],
    {pair, pairs}];
  Grid[
    Prepend[rows, hdr],
    Frame      -> All,
    FrameStyle -> GrayLevel[0.82],
    Background -> {None, None,
      Table[{i+1,4} -> If[!AnyTrue[viols, AssociationQ[#] && #["pair"]===pairs[[i]] &],
                           RGBColor[0.88,1.00,0.88],
                           RGBColor[1.00,0.88,0.88]],
            {i, Length[pairs]}]},
    Spacings   -> {2, 0.9},
    Alignment  -> {Left, Center}
  ]
]

(* ----------------------------------------------------------------
   MakeFrequencyPanel
   Scatter plot of simulated vs Boltzmann frequencies with a y=x
   reference line and Pearson correlation coefficient.
   x-axis: theoretical (Boltzmann) frequency.
   y-axis: simulated (MCMC) frequency.
   ---------------------------------------------------------------- *)
MakeFrequencyPanel[allStates_List, simFreq_Association,
                   bw_Association, kl_Real] := Module[
  {simD, bwD, pass, pearson, maxVal, plot, verdict},
  simD    = N[simFreq[#]] & /@ allStates;
  bwD     = N[bw[#]]      & /@ allStates;
  pass    = kl < 0.02;
  pearson = If[Length[allStates] > 1 && StandardDeviation[bwD] > 0,
               Correlation[simD, bwD], 1.];
  maxVal  = Max[Max[simD], Max[bwD]] * 1.15;
  plot = Show[
    ListPlot[
      Transpose[{bwD, simD}],
      PlotStyle  -> Directive[RGBColor[0.20, 0.45, 0.80], PointSize[0.018]],
      AxesLabel  -> {Style["Theoretical frequency", 10],
                     Style["Simulated frequency",   10]},
      PlotLabel  -> Style["Simulated vs Boltzmann state frequencies", 12, Bold],
      PlotRange  -> {{0, maxVal}, {0, maxVal}},
      AspectRatio -> 1,
      ImageSize  -> {380, 380},
      Background -> White],
    Graphics[{Dashed, GrayLevel[0.55],
              Line[{{0, 0}, {maxVal, maxVal}}]}]];
  verdict = Style[
    "Pearson r = " <> ToString[NumberForm[pearson, {5, 4}]] <>
    "     KL = " <> ToString[NumberForm[kl, {5, 4}]] <> "     " <>
    If[pass, "\[Checkmark] Consistent with Boltzmann",
             "\[Times] Significant deviation from Boltzmann"],
    11, Bold, If[pass, Darker[Green], Red]];
  Column[{plot, Spacer[6], verdict}, Alignment -> Left, Spacings -> 0.3]
]


(* ================================================================
   SECTION 5 – REPORT WINDOW
   ================================================================ *)

(* ----------------------------------------------------------------
   MakeReportWindow
   Assembles all results into a Mathematica notebook and opens it.
   args = Association with keys:
     name, allStates, treeData, matrix, violations,
     simFreq, bw, kl, algCode
   ---------------------------------------------------------------- *)
MakeReportWindow[args_Association] := Module[
  {name, allStates, treeData, matrix, violations,
   simFreq, bw, kl, algCode, pass,
   badge, trees, matGrid, dbTab, freqPanel, cells, nb},

  name       = args["name"];
  allStates  = args["allStates"];
  treeData   = args["treeData"];
  matrix     = args["matrix"];
  violations = args["violations"];
  simFreq    = args["simFreq"];
  bw         = args["bw"];
  kl         = args["kl"];
  algCode    = args["algCode"];
  pass       = violations === {};

  (* ---- Pass/fail banner ---- *)
  badge = Panel[
    Row[{
      Style[name, 18, Bold, GrayLevel[0.1]],
      Spacer[25],
      Style[If[pass,
               "\[FilledCircle]  DETAILED BALANCE:  PASS",
               "\[FilledCircle]  DETAILED BALANCE:  FAIL"],
            17, Bold,
            If[pass, RGBColor[0.04,0.54,0.04], RGBColor[0.74,0.04,0.04]]]
    }],
    Background -> GrayLevel[0.94],
    FrameMargins -> {{16,16},{12,12}}];

  (* ---- Trees: one per state, side by side ---- *)
  trees = Table[
    Framed[
      Column[{DrawStateTree[s, $normLeaf /@ treeData[s]]}, Alignment -> Center],
      FrameStyle -> GrayLevel[0.82], RoundingRadius -> 5,
      Background -> GrayLevel[0.975], FrameMargins -> 8],
    {s, allStates}];

  matGrid  = MakeTransitionGrid[allStates, matrix];
  dbTab    = MakeDBTable[allStates, violations];
  freqPanel = MakeFrequencyPanel[allStates, simFreq, bw, kl];

  (* ---- Notebook cells ---- *)
  cells = {
    Cell[BoxData @ ToBoxes @ badge,
         "Output", CellMargins -> {{8,8},{4,14}}],

    Cell["Algorithm Under Test", "Section"],
    Cell[algCode, "Code"],

    Cell["Decision Trees  (Symbolic Execution Paths)", "Section"],
    Cell[TextData[{
      "Each tree shows all bit sequences the algorithm can consume from a \
given starting state.  ",
      StyleBox["Blue", FontWeight->"Bold",
               FontColor->RGBColor[0.22,0.42,0.70]],
      " = root / internal bit-choice node.  ",
      StyleBox["Green", FontWeight->"Bold", FontColor->Darker[Green,0.2]],
      " outcome = algorithm moved to a new state.  ",
      StyleBox["Orange", FontWeight->"Bold", FontColor->Darker[Orange,0.1]],
      " outcome = algorithm stayed.  Edge labels are bit values (0 / 1).  \
Each bit either selects an integer (for RandomInteger / RandomChoice) or \
narrows the interval of a latent uniform variable U (for RandomReal[]).  \
Hover over an outcome node to see the exact path weight."
    }], "Text"],
    Cell[BoxData @ ToBoxes @
           Grid[Partition[trees, UpTo[3]], Spacings -> {2, 2}, Alignment -> {Left, Top}],
         "Output", CellMargins -> {{8,8},{4,4}}],

    Cell["Symbolic Transition Matrix  T[ i \[Rule] j ]", "Section"],
    Cell["Entry (i, j) = total probability of transitioning FROM state i \
TO state j, accumulated over all execution paths.  \[Beta] is kept symbolic \
throughout.", "Text"],
    Cell[BoxData @ ToBoxes @ matGrid,
         "Output", CellMargins -> {{8,8},{4,4}}],

    Cell["Detailed Balance Check:  \
T(i\[Rule]j)\[CenterDot]\[Pi](i) = T(j\[Rule]i)\[CenterDot]\[Pi](j)", "Section"],
    Cell["Every pair of distinct states is tested using FullSimplify with \
\[Beta] > 0 so that Piecewise Metropolis expressions are resolved exactly.  \
\[Pi](s) = Exp[\[Minus]\[Beta] E(s)].", "Text"],
    Cell[BoxData @ ToBoxes @ dbTab,
         "Output", CellMargins -> {{8,8},{4,4}}],

    Cell["Numerical MCMC Validation", "Section"],
    Cell["The algorithm is run as a live Markov chain with genuinely random \
bits.  Simulated state frequencies are compared to the analytical Boltzmann \
distribution Exp[\[Minus]\[Beta]E(s)] / Z.", "Text"],
    Cell[BoxData @ ToBoxes @ freqPanel,
         "Output", CellMargins -> {{8,8},{4,18}}]
  };

  (* ---- Open notebook window ---- *)
  nb = Quiet @ Check[
    CreateDocument[
      cells,
      WindowTitle  -> "DetailedBalanceChecker \[LongDash] " <> name,
      WindowSize   -> {1100, 860},
      WindowMargins -> {{Automatic, Automatic}, {Automatic, 0}},
      Editable     -> False,
      Background   -> White,
      StyleDefinitions -> "Default.nb"
    ],
    (Print["  (No Mathematica frontend available; window not opened.)"]; None)
  ];
  nb
]


(* ================================================================
   SECTION 6 – TOP-LEVEL ENTRY POINT
   ================================================================ *)

(* ----------------------------------------------------------------
   RunFullCheck
   Orchestrates BFS, symbolic check, numerical MCMC, and the
   graphical report window.

   Arguments:
     allStates   list of all valid system states
     symAlg      algorithm for symbolic check  (uses MetropolisProb)
     numAlg      algorithm for numerical MCMC  (fully numeric)
     symEnergy   bare energy, symbolic in couplings, no beta
     numEnergy   numeric energy WITH beta folded in

   Options:
     "SystemName"       display name
     "AlgorithmCode"    string shown in the Algorithm section
                        (Automatic = extracted from symAlg DownValues)
     "MaxBitDepth"      BFS depth cap per state (default 20)
     "TimeLimit"        seconds per state BFS (default 60)
     "Verbose"          print BFS progress (default True)
     "NSteps"           MCMC steps (default 100 000)
     "WarmupFrac"       warm-up fraction (default 0.1)
     "OpenWindow"       open graphical report window (default True)
   ---------------------------------------------------------------- *)
Options[RunFullCheck] = Join[
  Options[BuildTreeData],
  Options[RunNumericalMCMC],
  {"SystemName"    -> "Unnamed system",
   "AlgorithmCode" -> Automatic,
   "OpenWindow"    -> True}
]

RunFullCheck[allStates_List, symAlg_, numAlg_,
             symEnergy_, numEnergy_?(Head[#] =!= Rule && Head[#] =!= RuleDelayed &),
             OptionsPattern[]] := Module[
  {name     = OptionValue["SystemName"],
   n        = Length[allStates],
   algCode, treeData, matrix, violations,
   counts, bw, simFreq, kl, pass},

  algCode = OptionValue["AlgorithmCode"];
  If[algCode === Automatic,
    algCode = StringTrim @
              ToString[InputForm[DownValues[symAlg]], OutputForm]];

  (* ---- Terminal progress ---- *)
  Print[StringRepeat["=", 62]];
  Print["DETAILED BALANCE CHECKER"];
  Print[StringRepeat["=", 62]];
  Print["System : ", name];
  Print["States : ", n, "  --  ", allStates];
  Print[StringRepeat["=", 62]];

  (* ---- 1. Symbolic tree building ---- *)
  Print["\n[1/3] Building decision trees (symAlg) ..."];
  treeData = BuildTreeData[allStates, symAlg,
               "MaxBitDepth" -> OptionValue["MaxBitDepth"],
               "TimeLimit"   -> OptionValue["TimeLimit"],
               "Verbose"     -> OptionValue["Verbose"]];
  matrix   = TreeDataToMatrix[allStates, treeData];
  Print["      Transition matrix: ", Length[matrix], " non-zero entries."];

  (* ---- 2. Detailed balance check ---- *)
  Print["\n[2/3] Checking detailed balance for ",
        Binomial[n,2], " pairs ..."];
  violations = CheckDetailedBalance[matrix, allStates, symEnergy];
  pass       = violations === {};
  If[pass,
    Print["      RESULT: PASS -- all pairs satisfy detailed balance exactly."],
    Print["      RESULT: FAIL -- ", Length[violations], " violation(s) found."]
  ];

  (* ---- 3. Numerical MCMC ---- *)
  Print["\n[3/3] Running ", OptionValue["NSteps"],
        " MCMC steps (numAlg) ..."];
  counts  = RunNumericalMCMC[allStates, numAlg,
              "NSteps"     -> OptionValue["NSteps"],
              "WarmupFrac" -> OptionValue["WarmupFrac"]];
  bw      = BoltzmannWeights[allStates, numEnergy];
  simFreq = N[# / Total[counts]] & /@ counts;
  kl      = Total @ Table[
    With[{p = simFreq[s], q = N@bw[s]},
      If[p > 0 && q > 0, p * Log[p/q], 0.]], {s, allStates}];

  Print["      KL divergence (sim || Boltzmann) = ",
        NumberForm[kl,{5,4}]];
  Print["      Numerical: ",
    If[kl < 0.02, "CONSISTENT with Boltzmann.",
                  "WARNING -- significant deviation from Boltzmann."]];

  Print["\n", StringRepeat["=", 62]];
  Print["OVERALL: ", If[pass, "PASS", "FAIL"]];
  Print[StringRepeat["=", 62]];

  (* ---- Open graphical window / Python fallback ---- *)
  If[TrueQ @ OptionValue["OpenWindow"],
    With[{reportArgs = <|
        "name"       -> name,
        "allStates"  -> allStates,
        "treeData"   -> treeData,
        "matrix"     -> matrix,
        "violations" -> violations,
        "simFreq"    -> simFreq,
        "bw"         -> bw,
        "kl"         -> kl,
        "algCode"    -> algCode|>},
      Print["\nOpening report window ..."];
      nb = MakeReportWindow[reportArgs];
      If[nb === None,
        Print["  (No Mathematica frontend -- falling back to Python.)"];
        ExportAndShowPython[reportArgs]
      ]
    ]
  ];

  <|"pass"       -> pass,
    "violations" -> violations,
    "treeData"   -> treeData,
    "matrix"     -> matrix,
    "kl"         -> kl|>
]
