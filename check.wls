#!/usr/bin/env wolframscript
(* ================================================================
   check.wls  --  Exhaustive detailed-balance checker
   ================================================================
   Usage:
     wolframscript -file check.wls <alg_file> [options]

   The .wl file must define:
     energy[state_]       bare energy function (no beta)
     Algorithm[state_]    MCMC move function using native random calls
     BitsToState[bits_]   {0,1,...} list -> state, or None
     numBeta              numeric inverse temperature

   Optional (for symbolic check with energy parameters):
     symParams = <|"eps" -> {...}, "couplings" -> {...}|>
       -- static symbolic parameters (always the same across components)

     DynamicSymParams[states_List] -> <|"couplings" -> {...}|>
       -- called once per component; return the J symbols needed for
          that component (e.g. auto-generated from particle types).
          Takes precedence over / merges with static symParams.

   Options (key=value, no spaces around =):
     MaxBitString=XXXX    Largest bit string to test (default 11111111)
                          All strings from length 1 through len(XXXX),
                          stopping at XXXX within the final length.
     Mode=Symbolic        Symbolic detailed-balance check only
     Mode=Numerical       Numerical MCMC check only
     Mode=Both            Both checks (default)
     NSteps=N             MCMC steps per seed state (default 50000)
     MaxBitDepth=N        BFS tree depth cap per state (default 20)
     Verbose=True/False   Print per-state BFS progress (default False)
     SeedBitStrings=X,Y  Comma-separated list of bit strings to test
                          directly, bypassing MaxBitString enumeration.
                          BFS from each seed discovers the full component.
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ---- Launch parallel kernels ---- *)
$nKernels = 0;
Quiet[Check[LaunchKernels[]; $nKernels = Length[Kernels[]], $nKernels = 0]];
If[$nKernels > 0,
  Print["Parallel: ", $nKernels, " kernels launched."],
  Print["Note: parallel kernel launch failed — running single-threaded."]];

(* Distribute fast-checker helpers to subkernels (no-op if no kernels) *)
If[$nKernels > 0, $dbcDistributeFastChecker[]];

(* Clean-exit helper: always close sub-kernels before exiting so no
   WolframKernel processes are left consuming CPU after the script ends. *)
$dbcExit[code_Integer] := (If[$nKernels > 0, CloseKernels[]]; Exit[code]);

(* ---- Parse command-line arguments ---- *)
args = $ScriptCommandLine;
If[Length[args] >= 1 && StringMatchQ[args[[1]], ___ ~~ "check.wls"],
  args = Rest[args]];

posArgs = Select[args, !StringMatchQ[#, __ ~~ "=" ~~ __] &];
kvArgs  = Select[args,  StringMatchQ[#, __ ~~ "=" ~~ __] &];

If[Length[posArgs] < 1,
  Print["Usage: wolframscript -file check.wls <alg_file> [options]"];
  Print[""];
  Print["Required definitions in the .wl file:"];
  Print["  energy[state_]       bare energy (no beta)"];
  Print["  Algorithm[state_]    MCMC move function"];
  Print["  BitsToState[bits_]   {0,1,...} -> state, or None if invalid"];
  Print["  numBeta              numeric inverse temperature"];
  Print[""];
  Print["Options (key=value, no spaces around =):"];
  Print["  NGrid=N              Only test N×N grid states (recommended for 2D)"];
  Print["  MaxComponents=N      Stop after N distinct components (default: all)"];
  Print["  MaxBitString=XXXX    Largest bit string to enumerate (default 11111111)"];
  Print["  SeedBitStrings=X,Y  Comma-separated bit strings to test directly"];
  Print["  Mode=Symbolic|Numerical|Both  (default Both)"];
  Print["  NSteps=N             MCMC steps per system (default 50000)"];
  Print["  MaxBitDepth=N        BFS depth cap (default 20)"];
  Print["  Verbose=True/False   Per-state BFS output (default False)"];
  Print["  FailFast=1           Stop on first violation; print T, E, residual"];
  $dbcExit[1]];

(* ---- Option defaults ---- *)
maxBitString    = "11111111";
mode            = "Both";
nSteps          = 50000;
maxBitDepth     = 20;
verbose         = False;
fastChecker     = False;   (* FastChecker=1: use Exp-polynomial fast path *)
failFast        = False;   (* FailFast=1: stop+print debug info on first violation *)
seedBitStrings  = None;    (* SeedBitStrings=X,Y: test specific bit strings *)
physFidelity    = False;   (* PhysFidelity=1: compare T_MC to T_phys (Glauber) *)
maxComponents   = Infinity; (* MaxComponents=N: stop after N distinct components *)
nGrid           = None;    (* NGrid=N: only test N×N grid states *)

Do[
  If[StringMatchQ[opt, "MaxBitString="   ~~ __], maxBitString   = StringDrop[opt, 13]];
  If[StringMatchQ[opt, "SeedBitStrings=" ~~ __], seedBitStrings = StringDrop[opt, 15]];
  If[StringMatchQ[opt, "Mode="           ~~ __], mode           =              StringDrop[opt,  5]];
  If[StringMatchQ[opt, "NSteps="         ~~ __], nSteps         = ToExpression[StringDrop[opt,  7]]];
  If[StringMatchQ[opt, "MaxBitDepth="    ~~ __], maxBitDepth    = ToExpression[StringDrop[opt, 12]]];
  If[StringMatchQ[opt, "Verbose="        ~~ __], verbose        = ToExpression[StringDrop[opt,  8]]];
  If[StringMatchQ[opt, "FastChecker="    ~~ __], fastChecker    = (StringDrop[opt, 12] =!= "0")];
  If[StringMatchQ[opt, "FailFast="       ~~ __], failFast       = (StringDrop[opt,  9] =!= "0")];
  If[StringMatchQ[opt, "PhysFidelity="   ~~ __], physFidelity   = (StringDrop[opt, 12] =!= "0")];
  If[StringMatchQ[opt, "MaxComponents="  ~~ __], maxComponents  = ToExpression[StringDrop[opt, 14]]];
  If[StringMatchQ[opt, "NGrid="          ~~ __], nGrid          = ToExpression[StringDrop[opt,  6]]],
  {opt, kvArgs}];

(* ---- Validate MaxBitString ---- *)
If[!StringMatchQ[maxBitString, ("0" | "1")..],
  Print["ERROR: MaxBitString must contain only 0s and 1s, got: ", maxBitString];
  $dbcExit[1]];

(* ---- Resolve file path ---- *)
algFile = posArgs[[1]];
If[!FileExistsQ[algFile], algFile = scriptDir <> algFile];
If[!FileExistsQ[algFile],
  Print["ERROR: File not found: ", posArgs[[1]]]; $dbcExit[1]];

(* ---- Load algorithm file ---- *)
ClearAll[Algorithm, energy, BitsToState, numBeta, symParams, DynamicSymParams, DisplayState,
         fieldF, couplingJ,
         $fieldFConcrete, $pairEnergy, $concreteParams,
         $fieldFreeParams, $couplingFreeParams,
         $checkerAbstractParams, $symmetryGroup];
Print["Loading: ", algFile];
Get[algFile];

(* ---- Clear abstract parameters declared by the algorithm file ----
   Algorithm files declare $checkerAbstractParams = {"name1", "name2", ...}
   (string names, which survive evaluation unlike symbol references) to list
   parameters that must be unbound during symbolic BFS.
   Concrete values are captured BEFORE clearing so that numerical MCMC can
   restore them inside a Block. *)
$abstractParamRestore = <||>;
If[ListQ[$checkerAbstractParams] && AllTrue[$checkerAbstractParams, StringQ],
  (* Save current numeric values before clearing *)
  $abstractParamRestore = Association @ Flatten @ Map[
    Function[name,
      With[{val = Symbol[name]},
        If[NumericQ[val], {name -> val}, {}]]],
    $checkerAbstractParams];
  Scan[ClearAll, $checkerAbstractParams];
  bound = Select[$checkerAbstractParams, ToExpression["ValueQ[" <> # <> "]"] &];
  If[bound =!= {},
    Print["ERROR: parameter remains bound after clearing: ", bound];
    $dbcExit[1]];
  Print["Abstract params cleared: ", $checkerAbstractParams]];

(* Build symbolic assumptions for abstract parameters from their saved numeric values.
   A positive parameter (e.g. sigStep, physLen, epsLJ) gets {Element[sym,Reals], sym>0};
   negative gets {Element[sym,Reals], sym<0}; zero gets {Element[sym,Reals]}.
   These assumptions let FullSimplify apply real-valued identities (e.g. Erf[-x]=-Erf[x])
   that are otherwise blocked when the parameter's sign is unknown. *)
$abstractParamAssumptions = Flatten @ KeyValueMap[
  Function[{name, val},
    With[{sym = Symbol[name]},
      Which[val > 0, {Element[sym, Reals], sym > 0},
            val < 0, {Element[sym, Reals], sym < 0},
            True,    {Element[sym, Reals]}]]],
  $abstractParamRestore];

(* ---- Validate required symbols ---- *)
If[Length[DownValues[Algorithm]] == 0,
  Print["ERROR: 'Algorithm' not defined in ", posArgs[[1]]]; $dbcExit[1]];
If[Length[DownValues[energy]] == 0,
  Print["ERROR: 'energy' not defined in ", posArgs[[1]]]; $dbcExit[1]];
If[Length[DownValues[BitsToState]] == 0,
  Print["ERROR: 'BitsToState' not defined in ", posArgs[[1]]]; $dbcExit[1]];
If[!ValueQ[numBeta] || !NumericQ[numBeta],
  Print["ERROR: 'numBeta' (numeric) not defined in ", posArgs[[1]]]; $dbcExit[1]];

$alg         = Algorithm;
$energy      = energy;
$numBeta     = numBeta;
$hasDynParams = Length[DownValues[DynamicSymParams]] > 0;

(* ---- Verify couplingJ symmetry (if defined) ---- *)
(* Check type pairs {1..5} at d²∈{1,2,4,5,8,9,10} — covers all
   physically relevant cases.  Uses Simplify only: negligible overhead.
   An asymmetric coupling immediately implies a DB violation.
   When couplingJ has no DownValues (abstract during BFS), check $pairEnergy
   instead via a temporary Block. *)
$checkCouplingSymmetry[] := (
  $coupSymViolations = CheckCouplingSymmetry[{1,2,3,4,5}, {1,2,4,5,8,9,10}];
  If[$coupSymViolations =!= {},
    Print["ERROR: couplingJ is NOT symmetric — detailed balance cannot hold."];
    Print["  First violation: a=", $coupSymViolations[[1]]["a"],
          "  b=", $coupSymViolations[[1]]["b"],
          "  d2=", $coupSymViolations[[1]]["d2"],
          "  J[a,b]-J[b,a] = ", $coupSymViolations[[1]]["asymmetry"]];
    $dbcExit[1]];
  Print["couplingJ symmetry: PASS"]);
Which[
  Length[DownValues[couplingJ]] > 0,
    $checkCouplingSymmetry[],
  Length[DownValues[$pairEnergy]] > 0,
    Block[{couplingJ},
      couplingJ[a_Integer, b_Integer, d2_Integer] := $pairEnergy[a, b, d2];
      $checkCouplingSymmetry[]]];

(* ---- Symmetry group and oracle substitution ----
   When $symmetryGroup is declared AND the algorithm uses $vmmcCandidates (the
   standard candidate interface from vmmc_2d_grid.wl), the checker wraps $alg in
   Block[{$vmmcCandidates = $dbcCanonicalCandidates}, ...] during BFS.  This mirrors
   the couplingJ/$pairEnergy substitution pattern: the algorithm calls the standard
   function; the checker injects the canonical implementation without touching source.

   Canonical order makes G-related state pairs produce syntactically identical
   seqBernoulli trees, so $dbcDedup automatically collapses them — ~|G|-fold speedup.

   The checker also runs CheckGInvariance per component: for each generator g of G,
   it verifies T(s→s') === T(g(s)→g(s')) using syntactic equality.  This is a
   black-box check that does NOT read algorithm source; it only needs the T matrix.
   Checking 4 generators suffices for all |G| ≤ 72 elements (group-generation theorem).

   A G-invariance failure means the $symmetryGroup declaration is wrong.  It does NOT
   corrupt the DB check: $dbcDedup only collapses syntactically identical expressions,
   so non-G-invariant T entries are always evaluated separately.  G-invariance failure
   is reported as a diagnostic alongside the (always-correct) DB result. *)
$hasSym   = ValueQ[$symmetryGroup] && ListQ[$symmetryGroup] && Length[$symmetryGroup] > 0;
$hasField = Length[DownValues[fieldF]] > 0 || Length[DownValues[$fieldFConcrete]] > 0;

(* $algCallsOracle: True if the cluster builder actually calls $vmmcCandidates.
   Checked via DownValues so no algorithm source reading is needed. *)
$algCallsOracle = Length[DownValues[$vmmcBuildCluster]] > 0 &&
                  !FreeQ[DownValues[$vmmcBuildCluster], $vmmcCandidates];

(* $useCanonicalOracle: True when all conditions are met to substitute oracle + check G-inv *)
$useCanonicalOracle = $hasSym && !$hasField && $algCallsOracle;

If[$hasSym,
  Print["Symmetry group: ", StringRiffle[$symmetryGroup, ", "],
    Which[
      $hasField,
        "  (WARNING: field active — G-orbit dedup inapplicable)",
      $useCanonicalOracle,
        "  (oracle: $vmmcCandidates \[Rule] $dbcCanonicalCandidates during BFS; G-invariance check active)",
      True,
        "  (declared; algorithm uses raw neighbour loop — oracle substitution inactive)"]]];

(* $bfsAlg: algorithm used for BFS.  When oracle is active, wraps $alg so every
   call to $vmmcCandidates inside the algorithm is redirected to $dbcCanonicalCandidates.
   In real (non-BFS) runs $vmmcCandidates uses plain unordered neighbours — any ordering
   is correct by Whitelam-Geissler.  The substitution only affects checker BFS. *)
$bfsAlg = If[$useCanonicalOracle,
  Function[s, Block[{$vmmcCandidates = $dbcCanonicalCandidates}, $alg[s]]],
  $alg];

(* Static symbolic parameters (may be None if only DynamicSymParams is used) *)
$staticSysParams = If[ValueQ[symParams], symParams, None];
If[AssociationQ[$staticSysParams],
  $staticSymVars = Join[Lookup[$staticSysParams, "eps",       {}],
                        Lookup[$staticSysParams, "couplings", {}]],
  $staticSymVars = {}];

(* Fixed symbolic parameters: kept free during the symbolic check so that
   transition matrix entries contain exact symbolic expressions, but assigned
   their declared numeric values (not random) during numerical MCMC.
   Declared in the .wl file as:  symParams = <|"fixedParams" -> <|sym -> val, ...|>|> *)
$fixedParams  = If[AssociationQ[$staticSysParams],
                   Lookup[$staticSysParams, "fixedParams", <||>], <||>];
$fixedSymVars = Keys[$fixedParams];
$fixedNumVals = Values[$fixedParams];

(* ---- Set up iteration over candidate states ---- *)
maxBitVec = ToExpression /@ Characters[maxBitString];
maxLen     = Length[maxBitVec];
maxVal     = FromDigits[maxBitVec, 2];

(* When NGrid=N is set, ensure maxVal covers all N×N states so that the lazy
   iterator reaches them even if the user did not also set MaxBitString. *)
If[nGrid =!= None,
  maxVal = Max[maxVal, $cLPre[nGrid^2 + 1] - 1]];

(* Lazy mode: iterate state IDs by grid-size block, generating each bit string
   on the fly.  No large list is materialised in memory.
   Active when NGrid is set, or when ValidStateIDs is defined (2D algorithms).
   SeedBitStrings always uses the eager path (the explicit seed list is small).
   Non-2D algorithms (no ValidStateIDs, no NGrid) also use the eager path. *)
$useLazy = seedBitStrings === None &&
           (nGrid =!= None || Length[DownValues[ValidStateIDs]] > 0);

If[$useLazy,
  (* Grid-size blocks to iterate: just {nGrid^2} when NGrid is set, else all
     square sizes up to maxVal.  Computed once; IDs generated lazily in the loop. *)
  $lazyGridL2s =
    If[nGrid =!= None,
      {nGrid^2},
      Module[{L = 1, sizes = {}},
        While[$cLPre[L^2] <= maxVal, AppendTo[sizes, L^2]; L++];
        sizes]],
  (* Eager path: materialise a small list (SeedBitStrings or non-2D enumeration) *)
  $eagerBitStrings =
    If[seedBitStrings =!= None,
      Map[ToExpression /@ Characters[#] &, StringSplit[seedBitStrings, ","]],
      Join[
        Flatten[Table[IntegerDigits[i, 2, len],
                      {len, 1, maxLen - 1}, {i, 0, 2^len - 1}], 1],
        Table[IntegerDigits[i, 2, maxLen], {i, 0, maxVal}]]]];

(* ---- Column widths ---- *)
(* ID column: width for decimal integer representation of maxVal *)
colBits  = Max[6, IntegerLength[maxVal, 10] + 2];
colState = 22; colN = 9; colSym = 10; colFail = 8; colKL = 12; colNum = 10; colErgo = 14;
colM1 = 14; colM2 = 14;
totalCols = colBits + colState + colN + colErgo +
  If[mode =!= "Numerical", colSym + colFail, 0] +
  If[mode =!= "Symbolic",  colKL  + colNum,  0] +
  If[physFidelity, colM1 + colM2, 0];

(* ---- Print header ---- *)
Print[];
Print[StringRepeat["=", totalCols]];
Print["DETAILED BALANCE CHECKER"];
Print["File:          ", algFile];
Which[
  seedBitStrings =!= None,
    Print["SeedBitStrings: ", seedBitStrings,
          "  (", Length[$eagerBitStrings], " seed states)"],
  nGrid =!= None,
    Print["NGrid:         ", nGrid, "x", nGrid,
          "  (", $cLPre[nGrid^2+1] - $cLPre[nGrid^2], " state IDs -- all ",
          nGrid, "x", nGrid, " states)",
          If[IntegerQ[maxComponents], "  MaxComponents: " <> ToString[maxComponents], ""]],
  $useLazy,
    Print["MaxBitString:  ", maxBitString, "  (lazy: ",
          ToString[Length[$lazyGridL2s]], " grid size(s), IDs generated on-the-fly)",
          If[IntegerQ[maxComponents], "  MaxComponents: " <> ToString[maxComponents], ""]],
  True,
    Print["MaxBitString:  ", maxBitString,
          "  (", Length[$eagerBitStrings], " bit strings)",
          If[IntegerQ[maxComponents], "  MaxComponents: " <> ToString[maxComponents], ""]]];
Print["Mode:          ", mode];
If[mode =!= "Symbolic", Print["NSteps:        ", nSteps]];
Print["MaxBitDepth:   ", maxBitDepth];
If[fastChecker && mode =!= "Numerical",
  Print["SymChecker:    FastChecker (Exp-polynomial; FullSimplify fallback)"],
  If[mode =!= "Numerical",
    Print["SymChecker:    FullSimplify (standard)"]]];
If[failFast && mode =!= "Numerical",
  Print["FailFast:      Enabled (stop on first violation, print debug info)"]];
If[physFidelity, Print["PhysFidelity:  Enabled (T_phys Glauber; M1=EigRatio, M2=KLFidelity)"]];
Print[StringRepeat["=", totalCols]];
Print[];

Print[StringPadRight["ID",      colBits],
      StringPadRight["State",   colState],
      StringPadRight["#States", colN],
      StringPadRight["Ergodic", colErgo],
      If[mode =!= "Numerical",
         StringPadRight["Symbolic", colSym] <> StringPadRight["#Fail", colFail], ""],
      If[mode =!= "Symbolic",
         StringPadRight["KL", colKL] <> StringPadRight["Numerical", colNum], ""],
      If[physFidelity,
         StringPadRight["M1(EigRatio)", colM1] <> StringPadRight["M2(Fidelity)", colM2], ""]];
Print[StringRepeat["-", totalCols]];

(* ================================================================
   Progress bar helpers
   ================================================================ *)
(* In eager mode we know the total; in lazy mode we track component count. *)
$nBitStrings = If[$useLazy, Infinity, Length[$eagerBitStrings]];
$checkStart  = AbsoluteTime[];
$barWidth    = 28;

formatDuration[secs_] :=
  Which[secs < 60,   ToString[Round[secs]] <> "s",
        secs < 3600, ToString[Floor[secs/60]] <> "m " <>
                     ToString[Round[Mod[secs, 60]]] <> "s",
        True,        ToString[Floor[secs/3600]] <> "h " <>
                     ToString[Floor[Mod[secs/60, 60]]] <> "m"];

writeProgress[idx_] := If[verbose, Module[{elapsed, elStr},
  elapsed = AbsoluteTime[] - $checkStart;
  elStr   = formatDuration[elapsed];
  If[$useLazy,
    (* Lazy mode: show component count; total unknown until MaxComponents is reached *)
    WriteString[$Output,
      "\r[" <> ToString[testedCount] <> " component" <>
      If[testedCount =!= 1, "s", ""] <>
      If[IntegerQ[maxComponents],
         "/" <> ToString[maxComponents], " found"] <>
      "  elapsed: " <> elStr <> "]   "],
    (* Eager mode: traditional idx/total progress bar *)
    Module[{pct, filled, arrow, bar, eta, etaStr},
      pct    = If[$nBitStrings > 0, idx / $nBitStrings, 0.];
      filled = Round[pct * $barWidth];
      arrow  = If[filled < $barWidth, ">", ""];
      bar    = "[" <> StringRepeat["=", filled] <> arrow <>
               StringRepeat[" ", Max[0, $barWidth - filled - StringLength[arrow]]] <> "]";
      eta    = If[idx > 0 && elapsed > 0,
                  elapsed * ($nBitStrings - idx) / idx,
                  Infinity];
      etaStr = If[!NumberQ[eta] || eta === Infinity, "?", formatDuration[eta]];
      WriteString[$Output,
        "\r" <> bar <> " " <> ToString[idx] <> "/" <> ToString[$nBitStrings] <>
        "  " <> ToString[Round[pct * 100]] <> "%  elapsed: " <> elStr <>
        "  ETA: " <> etaStr <> "   "]]]]];  (* trailing spaces clear previous longer text *)

(* Before printing a permanent row: clear the progress line first *)
clearProgress[] := If[verbose, WriteString[$Output, "\r" <> StringRepeat[" ", 78] <> "\r"]];

(* ================================================================
   Main loop
   ================================================================ *)
coveredAssoc   = <||>;   (* Association for O(1) covered-state lookup *)
testedCount    = 0;
symPassCount   = 0; symFailCount   = 0;
numPassCount   = 0; numFailCount   = 0;
ergoPassCount  = 0; ergoFailCount  = 0;
gInvPassCount  = 0; gInvFailCount  = 0;
$loopIdx       = 0;

(* ---- Loop body: process one bit string.
        Uses Throw[Null, $dbcSkip] in place of Continue[].
        Returns True when MaxComponents has been reached (signal outer loop to stop). ---- *)
$loopBody[bitString_List] := Catch[

  state = Check[BitsToState[bitString], $Failed];
  If[state === $Failed || state === None,
    writeProgress[$loopIdx]; Throw[Null, $dbcSkip]];
  If[KeyExistsQ[coveredAssoc, state],
    writeProgress[$loopIdx]; Throw[Null, $dbcSkip]];

  bitsStr = ToString[FromDigits[bitString, 2]];
  stateStr = If[Length[DownValues[DisplayState]] > 0,
               DisplayState[state], ToString[state]];

  (* ---- BFS: discover connected component ---- *)
  treeData = BuildTreeAT[state, $bfsAlg,
               "MaxBitDepth" -> maxBitDepth,
               "TimeLimit"   -> 120.,
               "Verbose"     -> verbose];

  If[MatchQ[treeData, $dbc$cantHandle[_]],
    clearProgress[];
    With[{eRows = StringSplit[stateStr, "|"]},
      Print[StringPadRight[bitsStr, colBits], StringPadRight[eRows[[1]], colState],
            StringPadRight["?", colN], "ERROR: ", treeData[[1]]];
      Do[Print[StringRepeat[" ", colBits], eRows[[r]]], {r, 2, Length[eRows]}]];
    writeProgress[$loopIdx];
    testedCount++; Throw[Null, $dbcSkip]];

  allStates = Sort @ Keys[treeData];
  matrix    = TreeATToMatrix[treeData];
  nStates   = Length[allStates];
  coveredAssoc = Join[coveredAssoc, Association[# -> True & /@ allStates]];
  testedCount++;

  (* ---- Ergodicity check ---- *)
  ergoResult = CheckErgodicity[allStates];
  ergoStr    = If[ergoResult["ergodic"],
    ergoPassCount++; "PASS (" <> ToString[ergoResult["found"]] <> ")",
    ergoFailCount++; "FAIL (" <> ToString[ergoResult["found"]] <>
                     "/" <> ToString[ergoResult["theoretical"]] <> ")"];

  (* ---- Per-component symbolic parameters ---- *)
  (* If DynamicSymParams is defined, call it with this component's states
     to get the symbolic variables needed for this component.
     "couplings" / "eps" — symbolic atoms for FullSimplify assumptions.
     "numericParams"     — scalar symbols for the numerical MCMC Block. *)
  If[$hasDynParams,
    dynP            = DynamicSymParams[allStates];
    dynVars         = Join[Lookup[dynP, "eps", {}], Lookup[dynP, "couplings", {}]];
    dynNumVars      = Lookup[dynP, "numericParams", dynVars];
    compSymVars     = Union[$staticSymVars, dynVars];
    compNumericVars = Union[$staticSymVars, dynNumVars],
    compSymVars     = $staticSymVars;
    compNumericVars = $staticSymVars];
  compSymAssumptions = Join[
    Map[Element[#, Reals] &, Union[compSymVars, $fixedSymVars]],
    $abstractParamAssumptions];

  (* ---- Symbolic check ---- *)
  symResult = "-"; failStr = "-";
  If[mode =!= "Numerical",
    With[{viol = If[fastChecker,
                  CheckDetailedBalanceFast[matrix, allStates, $energy, compSymAssumptions, failFast],
                  CheckDetailedBalance[    matrix, allStates, $energy, compSymAssumptions, failFast]]},
      symResult = If[viol === {}, "PASS", "FAIL"];
      failStr   = If[viol === {}, "0",    ToString[Length[viol]]];
      If[symResult === "PASS", symPassCount++, symFailCount++];
      (* FailFast: print violation details and exit immediately *)
      If[failFast && viol =!= {},
        With[{v = viol[[1]],
              $ds = If[Length[DownValues[DisplayState]] > 0, DisplayState, ToString]},
          clearProgress[];
          Print[];
          Print["=== FAIL-FAST: first DB violation ==="];
          Print["  Component seed: ", bitsStr, "  (", nStates, " states)"];
          Print["  State 1 : ", $ds[v["pair"][[1]]]];
          Print["  State 2 : ", $ds[v["pair"][[2]]]];
          Print["  T(1->2) = ", v["tij"]];
          Print["  T(2->1) = ", v["tji"]];
          Print["  E(1)    = ", v["ei"]];
          Print["  E(2)    = ", v["ej"]];
          Print["  Residual (T(1->2)*Exp[-b*E1] - T(2->1)*Exp[-b*E2]) ="];
          Print["    ", v["residual"]]];
        Print[];
        Print[StringRepeat["=", totalCols]];
        Print["*** DETAILED BALANCE VIOLATION DETECTED (FailFast) ***"];
        $dbcExit[1]]]];

  (* ---- Numerical MCMC check ---- *)
  klStr = "-"; numResult = "-";
  If[mode =!= "Symbolic" && nStates > 1,
    allBlockVars = Union[compNumericVars, $fixedSymVars];
    $runNumericalMCMC[] := Module[{},
      mcmcCounts = RunNumericalMCMCAT[allStates, $alg, $numBeta, "NSteps" -> nSteps];
      boltzmann  = BoltzmannWeightsAT[allStates, $energy, $numBeta]];
    $activateAndRun[] :=
      If[Length[DownValues[$pairEnergy]] > 0 || Length[DownValues[$fieldFConcrete]] > 0,
        Block[{fieldF, couplingJ},
          If[Length[DownValues[$fieldFConcrete]] > 0,
            fieldF[x_Integer, y_Integer, L_Integer] := $fieldFConcrete[x, y, L]];
          If[Length[DownValues[$pairEnergy]] > 0,
            couplingJ[a_Integer, b_Integer, d2_Integer] := $pairEnergy[a, b, d2]];
          $runNumericalMCMC[]],
        $runNumericalMCMC[]];
    $restoreSyms     = Map[Symbol, Keys[$abstractParamRestore]];
    $restoreVals     = Values[$abstractParamRestore];
    allBlockVarsFull = Union[allBlockVars, $restoreSyms];
    If[allBlockVarsFull =!= {},
      SeedRandom[42];
      compNumVals = If[compNumericVars =!= {}, RandomReal[{-1, 1}, Length[compNumericVars]], {}];
      Block[Evaluate[allBlockVarsFull],
        If[compNumericVars =!= {}, MapThread[Set, {compNumericVars, compNumVals}]];
        If[$fixedSymVars =!= {}, MapThread[Set, {$fixedSymVars, $fixedNumVals}]];
        If[$restoreSyms =!= {}, MapThread[Set, {$restoreSyms, $restoreVals}]];
        $activateAndRun[]],
      $activateAndRun[]];
    simFreq = N[# / Total[mcmcCounts]] & /@ mcmcCounts;
    kl = Total @ Table[
      With[{p = simFreq[s], q = N @ boltzmann[s]},
        If[p > 0 && q > 0, p * Log[p / q], 0.]], {s, allStates}];
    klStr     = Which[Abs[kl] < 5.*^-4, "~0", True, ToString[NumberForm[kl, {4, 3}]]];
    numResult = If[kl < 0.02, "PASS", "FAIL"];
    If[numResult === "PASS", numPassCount++, numFailCount++]];

  (* ---- Physical fidelity metrics (if requested) ---- *)
  (* Builds T_phys (Glauber single-particle) and compares to T_MC.
     Uses the same Block structure as numerical MCMC so all symbolic
     atoms (couplingJ, physLen, epsLJ, \[Beta]) evaluate numerically. *)
  m1Str = "-"; m2Str = "-";
  If[physFidelity,
    $pfRestoreSyms = Map[Symbol, Keys[$abstractParamRestore]];
    $pfRestoreVals = Values[$abstractParamRestore];
    $pfAllBlockVars = Union[compNumericVars, $fixedSymVars, $pfRestoreSyms, {\[Beta]}];
    $computePhysFidelity[] := Module[{tPhys, tMCmat, bw, m1, m2},
      tPhys  = BuildTPhys[allStates, $energy, $numBeta];
      If[tPhys === None, m1Str = "N/A"; m2Str = "N/A"; Return[]];
      tMCmat = $tAssocToNumericMatrix[matrix, allStates];
      bw     = N[Exp[-$numBeta * $energy[#]] & /@ allStates];
      m1     = ComputeMetric1[tMCmat, tPhys];
      m2     = ComputeMetric2[tMCmat, tPhys, bw];
      m1Str  = If[m1 === Indeterminate, "N/A",
                  If[Abs[m1] < 1.*^-6, "~0",
                     ToString[Round[m1, 0.001]]]];
      m2Str  = Which[m2 === -Infinity, "-Inf",
                     m2 === Indeterminate, "N/A",
                     True, ToString[Round[m2, 0.001]]]];
    $pfActivateAndRun[] :=
      If[Length[DownValues[$pairEnergy]] > 0 || Length[DownValues[$fieldFConcrete]] > 0,
        Block[{fieldF, couplingJ},
          If[Length[DownValues[$fieldFConcrete]] > 0,
            fieldF[x_Integer, y_Integer, L_Integer] := $fieldFConcrete[x, y, L]];
          If[Length[DownValues[$pairEnergy]] > 0,
            couplingJ[a_Integer, b_Integer, d2_Integer] := $pairEnergy[a, b, d2]];
          $computePhysFidelity[]],
        $computePhysFidelity[]];
    If[$pfAllBlockVars =!= {},
      SeedRandom[42];
      $pfCompVals = If[compNumericVars =!= {}, RandomReal[{-1, 1}, Length[compNumericVars]], {}];
      Block[Evaluate[$pfAllBlockVars],
        If[compNumericVars =!= {}, MapThread[Set, {compNumericVars, $pfCompVals}]];
        If[$fixedSymVars =!= {}, MapThread[Set, {$fixedSymVars, $fixedNumVals}]];
        If[$pfRestoreSyms =!= {}, MapThread[Set, {$pfRestoreSyms, $pfRestoreVals}]];
        \[Beta] = $numBeta;
        $pfActivateAndRun[]],
      \[Beta] = $numBeta;
      $pfActivateAndRun[]]];

  (* ---- Print result row ---- *)
  clearProgress[];
  stateRows = StringSplit[stateStr, "|"];
  Print[StringPadRight[bitsStr,             colBits],
        StringPadRight[stateRows[[1]],      colState],
        StringPadRight[ToString[nStates],   colN],
        StringPadRight[ergoStr,             colErgo],
        If[mode =!= "Numerical",
           StringPadRight[symResult, colSym] <> StringPadRight[failStr, colFail], ""],
        If[mode =!= "Symbolic",
           StringPadRight[klStr, colKL] <> StringPadRight[numResult, colNum], ""],
        If[physFidelity,
           StringPadRight[m1Str, colM1] <> StringPadRight[m2Str, colM2], ""]];
  Do[Print[StringRepeat[" ", colBits], stateRows[[r]]], {r, 2, Length[stateRows]}];
  If[Length[stateRows] > 1, Print[]];
  If[$hasSym && mode =!= "Numerical" && $dbcLastDedupTotal > 0,
    Print["  G-orbit dedup: ", $dbcLastDedupUnique, "/", $dbcLastDedupTotal,
          " unique  (",
          Round[100. * (1. - $dbcLastDedupUnique / N[$dbcLastDedupTotal])],
          "% hash-collapsed)"]];

  (* ---- G-invariance check ----
     Verifies T(s→s') === T(g(s)→g(s')) for each declared generator g.
     Only runs when the canonical oracle was active during BFS (i.e. the
     algorithm uses $vmmcCandidates).  Syntactic === comparison is valid
     because canonical ordering makes G-related pairs produce identical
     expressions; no FullSimplify is needed. *)
  If[$useCanonicalOracle,
    With[{gviol = CheckGInvariance[matrix, allStates, $symmetryGroup]},
      If[gviol === {},
        gInvPassCount++;
        Print["  G-invariance: PASS"],
        gInvFailCount++;
        With[{v1 = gviol[[1]]},
          Print["  G-invariance: FAIL  (", Length[gviol], " generator(s) violated)",
                "  — ", v1["gen"],
                If[v1["issue"] === "orbit-escapes-component",
                  ": g(s) not in component for s=" <> ToString[v1["s"]],
                  ": T mismatch at s=" <> ToString[v1["s"]] <> " s2=" <> ToString[v1["s2"]]]]];
        Print["  WARNING: $symmetryGroup declaration may be incorrect; ",
              "DB check is still correct (dedup only collapses hash-identical exprs)."]]]];

  writeProgress[$loopIdx];
  (* Return True when MaxComponents is reached — signals the outer loop to stop *)
  IntegerQ[maxComponents] && testedCount >= maxComponents,

  $dbcSkip];  (* Catch tag — Throw[Null, $dbcSkip] acts as Continue[] in the outer loop *)

(* ---- Dispatch: lazy (NGrid / ValidStateIDs) or eager (SeedBitStrings / 1D) ---- *)
If[$useLazy,
  (* Lazy path: iterate state IDs by grid-size block without materialising a list.
     The outer Do steps through grid sizes; the inner Do steps through IDs in that block. *)
  $lazyDone = False;
  Do[
    If[$lazyDone, Break[]];
    Module[{L2 = $lazyGridL2s[[gIdx]], idStart, idStop},
      idStart = $cLPre[L2];
      idStop  = Min[maxVal, $cLPre[L2 + 1] - 1];
      If[idStart > maxVal, Continue[]];   (* this grid size lies beyond maxVal *)
      Do[
        If[$lazyDone, Break[]];
        $loopIdx++;
        If[TrueQ[$loopBody[IntegerDigits[id, 2]]],
          $lazyDone = True; Break[]],
        {id, idStart, idStop}]],
    {gIdx, Length[$lazyGridL2s]}],
  (* Eager path: iterate the pre-materialised list (SeedBitStrings or non-2D enumeration).
     TrueQ guards against the $dbcSkip symbol returned when Throw[] is caught. *)
  Do[
    $loopIdx++;
    If[TrueQ[$loopBody[bitString]], Break[]],
    {bitString, $eagerBitStrings}]];

(* Final newline to clear progress bar before the summary *)
If[verbose, WriteString[$Output, "\n"]];

(* ================================================================
   Summary
   ================================================================ *)
Print[];
Print[StringRepeat["-", totalCols]];
Print["Systems tested (distinct connected components): ", testedCount];
Print["Ergodicity -- PASS: ", ergoPassCount, "  FAIL: ", ergoFailCount];
If[mode =!= "Numerical",
  Print["Symbolic   -- PASS: ", symPassCount, "  FAIL: ", symFailCount]];
If[mode =!= "Symbolic" && (numPassCount + numFailCount) > 0,
  Print["Numerical  -- PASS: ", numPassCount, "  FAIL: ", numFailCount]];
If[$useCanonicalOracle && (gInvPassCount + gInvFailCount) > 0,
  Print["G-invariance -- PASS: ", gInvPassCount, "  FAIL: ", gInvFailCount]];
If[ergoFailCount > 0,
  Print["*** ERGODICITY VIOLATION DETECTED ***"]];
If[$useCanonicalOracle && gInvFailCount > 0,
  Print["*** G-INVARIANCE VIOLATION DETECTED — $symmetryGroup declaration is incorrect ***"]];
If[symFailCount == 0 && numFailCount == 0 && testedCount > 0,
  Print["All tested systems satisfy detailed balance."],
  If[testedCount == 0,
    Print["No valid states found for the given parameters."],
    Print["*** DETAILED BALANCE VIOLATION DETECTED ***"]]];
If[IntegerQ[maxComponents] && testedCount >= maxComponents,
  Print["(Stopped after MaxComponents=", maxComponents, " -- further components not checked.)"]];
Print[StringRepeat["=", totalCols]];

(* Shut down sub-kernels cleanly so no WolframKernel processes linger. *)
If[$nKernels > 0, CloseKernels[]];
