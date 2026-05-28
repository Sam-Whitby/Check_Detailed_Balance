(* ================================================================
   template.wl  --  Template for the detailed-balance checker
   ================================================================
   Copy this file, fill in the required definitions, and run:

     wolframscript -file check.wls your_algorithm.wl MaxBitString=<bits>
     wolframscript -file report.wls your_algorithm.wl BitString=<bits>

   REQUIRED:  energy, Algorithm, BitsToState, numBeta
   OPTIONAL:  DynamicSymParams or symParams  (symbolic energy parameters)
              DisplayState                   (custom state display)
              ValidStateIDs                  (fast enumeration of valid seeds)
   ================================================================ *)


(* ---- Energy function ---------------------------------------------------- *)
(*
   Return the bare energy of state (a real number, no β factor).
   Must be deterministic -- no random calls.
   May reference symbolic variables declared in symParams / DynamicSymParams.
*)
energy[state_] :=
  0   (* replace with your energy *)


(* ---- MCMC algorithm ------------------------------------------------------- *)
(*
   Take a state, return the next state (same type).
   Use native Mathematica random calls: RandomInteger[], RandomReal[],
   RandomChoice[].

   For Metropolis acceptance use exactly:

     If[RandomReal[] < MetropolisProb[dE], newState, state]

   where dE = energy[newState] - energy[state].
   MetropolisProb[dE] is intercepted by the checker; symbolically it acts as:

     Piecewise[{{1, dE <= 0}, {Exp[-\[Beta]*dE], dE > 0}}]

   with \[Beta] kept as a free symbol during the symbolic check so that the
   Boltzmann factor cancels algebraically in the detailed-balance ratio.
*)
Algorithm[state_] :=
  state   (* replace with your algorithm *)


(* ---- Symbolic energy parameters (optional) ------------------------------ *)
(*
   Declare any symbolic coupling constants or site energies here so the
   checker can treat them as free symbols during FullSimplify and assign
   random numeric values for the numerical MCMC check.

   Option A -- static parameters (the same for every connected component):

     symParams = <|"eps" -> {eps1, eps2}, "couplings" -> {J1, J2}|>

   Option B -- dynamic parameters (generated per component from particle types
   actually present, e.g. for systems where the relevant J symbols change
   depending on which species are in the component):

     DynamicSymParams[states_List] :=
       Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
         <|"couplings" ->
           Flatten @ Table[
             If[a < b, ToExpression["J" <> ToString[a] <> ToString[b]], Nothing],
             {a, types}, {b, types}]|>]

   Both may be used together.  If neither is defined the energy is treated as
   fully numeric (symbolic check is skipped for any component whose energy
   evaluates to a number without extra symbols).
*)


(* ---- BitsToState --------------------------------------------------------- *)
(*
   Convert a bit string (list of 0s and 1s) to a seed state, or return None
   to skip this bit string.  The checker calls BitsToState for every bit
   string up to MaxBitString, then uses BFS from the returned seed to discover
   the full connected component automatically.

   Guidelines:
   - Accept lists of any length; return None for patterns that don't
     correspond to a valid state (wrong length, impossible configuration, etc.).
   - Different bit strings that return non-None should map to states in
     DIFFERENT connected components, to avoid redundant BFS work.
   - The returned state must be accepted by both energy and Algorithm.

   Simple example: fixed 4-site ring, state = list of occupancies {0,1}
     BitsToState[bits_List] :=
       If[Length[bits] =!= 4, None, bits]

   Labeled-particle example (kawasaki_1d.wl encoding):
     BitsToState[bits_List] :=
       Module[{id = FromDigits[bits, 2]},
         If[id == 0, None, $decode[id]]]
*)
BitsToState[bits_List] :=
  None   (* replace with your BitsToState *)


(* ---- Inverse temperature ------------------------------------------------ *)
(*
   Numeric value of β = 1/(k_B T) used for:
     - the numerical MCMC check (RunNumericalMCMCAT)
     - the Boltzmann weights the MCMC distribution is compared against
   The symbolic check uses \[Beta] as a free symbol regardless of numBeta.
   Set to 1 unless you have a specific reason to use a different value.
*)
numBeta = 1


(* ---- DisplayState (optional) -------------------------------------------- *)
(*
   If defined, the checker calls DisplayState[state] instead of
   ToString[state] for the State column in the output table.
   Useful for 2D states, e.g.:
     DisplayState[state_List] :=
       With[{L = Round[Sqrt[Length[state]]]},
         StringJoin @ Riffle[
           Table["{" <> StringRiffle[ToString /@ state[[(r-1)*L+1 ;; r*L]], ","] <> "}",
                 {r, 1, L}], "|"]]
*)


(* ---- ValidStateIDs (optional optimisation) ------------------------------ *)
(*
   If defined, check.wls calls ValidStateIDs[maxId] instead of exhaustively
   iterating every bit string from 1 to maxId.  Return the list of integer IDs
   (in the same bijective encoding used by BitsToState) that actually decode to
   valid seed states.  This can dramatically reduce pre-check overhead when only
   a small fraction of bit strings map to valid states.

   Example: 2D algorithms where valid states have perfect-square length (L×L):
     ValidStateIDs[maxId_Integer] :=
       Module[{L = 1, ids = {}},
         While[$cLPre[L^2] <= maxId,
           ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2 + 1] - 1, maxId]]];
           L++];
         ids]
*)


(* ---- $symmetryGroup (optional) ----------------------------------------- *)
(*
   Declare the symmetry group of your algorithm so the checker can exploit
   G-orbit deduplication.  When declared, the canonical-neighbour oracle in
   $dbcCanonicalCandidates (vmmc_2d_grid.wl) sorts VMMC candidate neighbours
   in an order that is invariant under the declared symmetries.  This makes
   seqBernoulli trees syntactically identical for G-related state pairs, so
   the existing $dbcDedup hash step automatically collapses them — giving an
   ~|G|-fold speedup in FullSimplify calls with no extra verification logic.

   Supported values:
     $symmetryGroup = {"translation"}           (* translational invariance only *)
     $symmetryGroup = {"translation", "D4"}     (* translations + 4-fold rotation + reflection *)

   Requirements for the declaration to be valid:
     1. Seed selection is uniform over occupied sites (no position weighting).
     2. Displacement list contains every (dx,dy) paired with (-dx,-dy).
     3. Energy uses only minimum-image periodic distances (no external field).

   If any of these fail, the declaration is incorrect.  An incorrect declaration
   does NOT corrupt the DB check — each pair is still evaluated — but the dedup
   ratio will be lower than claimed.  $dbcDedup simply collapses fewer pairs.

   The check.wls script prints the declared group and dedup statistics per
   component so you can verify the speedup empirically.

   For algorithms with an external field (e.g. vmmc_2d_field.wl) DO NOT declare
   $symmetryGroup: the field breaks translational invariance, so G-related pairs
   produce hash-distinct DB expressions and no orbit dedup occurs.

   This feature is only meaningful for 2D algorithms using $dbcCanonicalCandidates.
   For 1D or other geometries, omit $symmetryGroup entirely.
*)
(* $symmetryGroup = {"translation", "D4"} *)   (* uncomment if applicable *)
