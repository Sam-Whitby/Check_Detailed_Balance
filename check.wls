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
  Print["  MaxBitString=XXXX    Largest bit string (default 11111111)"];
  Print["  SeedBitStrings=X,Y  Comma-separated bit strings to test directly"];
  Print["  Mode=Symbolic|Numerical|Both  (default Both)"];
  Print["  NSteps=N             MCMC steps per system (default 50000)"];
  Print["  MaxBitDepth=N        BFS depth cap (default 20)"];
  Print["  Verbose=True/False   Per-state BFS output (default False)"];
  Exit[1]];

(* ---- Option defaults ---- *)
maxBitString    = "11111111";
mode            = "Both";
nSteps          = 50000;
maxBitDepth     = 20;
verbose         = False;
fastChecker     = False;   (* FastChecker=1: use Exp-polynomial fast path *)
seedBitStrings  = None;    (* SeedBitStrings=X,Y: test specific bit strings *)

Do[
  If[StringMatchQ[opt, "MaxBitString="   ~~ __], maxBitString   = StringDrop[opt, 13]];
  If[StringMatchQ[opt, "SeedBitStrings=" ~~ __], seedBitStrings = StringDrop[opt, 15]];
  If[StringMatchQ[opt, "Mode="           ~~ __], mode           =              StringDrop[opt,  5]];
  If[StringMatchQ[opt, "NSteps="         ~~ __], nSteps         = ToExpression[StringDrop[opt,  7]]];
  If[StringMatchQ[opt, "MaxBitDepth="    ~~ __], maxBitDepth    = ToExpression[StringDrop[opt, 12]]];
  If[StringMatchQ[opt, "Verbose="        ~~ __], verbose        = ToExpression[StringDrop[opt,  8]]];
  If[StringMatchQ[opt, "FastChecker="    ~~ __], fastChecker    = (StringDrop[opt, 12] =!= "0")],
  {opt, kvArgs}];

(* ---- Validate MaxBitString ---- *)
If[!StringMatchQ[maxBitString, ("0" | "1")..],
  Print["ERROR: MaxBitString must contain only 0s and 1s, got: ", maxBitString];
  Exit[1]];

(* ---- Resolve file path ---- *)
algFile = posArgs[[1]];
If[!FileExistsQ[algFile], algFile = scriptDir <> algFile];
If[!FileExistsQ[algFile],
  Print["ERROR: File not found: ", posArgs[[1]]]; Exit[1]];

(* ---- Load algorithm file ---- *)
ClearAll[Algorithm, energy, BitsToState, numBeta, symParams, DynamicSymParams, DisplayState,
         fieldF, couplingJ,
         $fieldFConcrete, $couplingJConcrete, $abstractFunctions, $concreteParams,
         $fieldFreeParams, $couplingFreeParams];
Print["Loading: ", algFile];
Get[algFile];

(* ---- Validate required symbols ---- *)
If[Length[DownValues[Algorithm]] == 0,
  Print["ERROR: 'Algorithm' not defined in ", posArgs[[1]]]; Exit[1]];
If[Length[DownValues[energy]] == 0,
  Print["ERROR: 'energy' not defined in ", posArgs[[1]]]; Exit[1]];
If[Length[DownValues[BitsToState]] == 0,
  Print["ERROR: 'BitsToState' not defined in ", posArgs[[1]]]; Exit[1]];
If[!ValueQ[numBeta] || !NumericQ[numBeta],
  Print["ERROR: 'numBeta' (numeric) not defined in ", posArgs[[1]]]; Exit[1]];

$alg         = Algorithm;
$energy      = energy;
$numBeta     = numBeta;
$hasDynParams = Length[DownValues[DynamicSymParams]] > 0;

(* ---- Verify couplingJ symmetry (if defined) ---- *)
(* Check type pairs {1..5} at d²∈{1,2,4,5,8,9,10} — covers all
   physically relevant cases.  Uses Simplify only: negligible overhead.
   An asymmetric coupling immediately implies a DB violation.
   In abstract-function mode ($abstractFunctions=True), couplingJ has no
   DownValues; we check $couplingJConcrete instead via a temporary Block. *)
$checkCouplingSymmetry[] := (
  $coupSymViolations = CheckCouplingSymmetry[{1,2,3,4,5}, {1,2,4,5,8,9,10}];
  If[$coupSymViolations =!= {},
    Print["ERROR: couplingJ is NOT symmetric — detailed balance cannot hold."];
    Print["  First violation: a=", $coupSymViolations[[1]]["a"],
          "  b=", $coupSymViolations[[1]]["b"],
          "  d2=", $coupSymViolations[[1]]["d2"],
          "  J[a,b]-J[b,a] = ", $coupSymViolations[[1]]["asymmetry"]];
    Exit[1]];
  Print["couplingJ symmetry: PASS"]);
Which[
  Length[DownValues[couplingJ]] > 0,
    $checkCouplingSymmetry[],
  Length[DownValues[$couplingJConcrete]] > 0,
    Block[{couplingJ},
      couplingJ[a_Integer, b_Integer, d2_Integer] := $couplingJConcrete[a, b, d2];
      $checkCouplingSymmetry[]]];

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

(* ---- Generate candidate bit strings ---- *)
maxBitVec = ToExpression /@ Characters[maxBitString];
maxLen     = Length[maxBitVec];
maxVal     = FromDigits[maxBitVec, 2];

(* SeedBitStrings=X,Y overrides MaxBitString enumeration entirely.
   Each comma-separated token is treated as a literal bit string.
   BFS from each seed will discover the full component as usual. *)
allBitStrings =
  If[seedBitStrings =!= None,
    (* Parse comma-separated bit strings and convert chars to integers *)
    Map[ToExpression /@ Characters[#] &, StringSplit[seedBitStrings, ","]],
    (* Otherwise: enumerate from MaxBitString (with ValidStateIDs if available) *)
    If[Length[DownValues[ValidStateIDs]] > 0,
      IntegerDigits[#, 2] & /@ ValidStateIDs[maxVal],
      Join[
        Flatten[Table[IntegerDigits[i, 2, len],
                      {len, 1, maxLen - 1}, {i, 0, 2^len - 1}], 1],
        Table[IntegerDigits[i, 2, maxLen], {i, 0, maxVal}]]]];

(* ---- Column widths ---- *)
(* When SeedBitStrings is used, base column width on the longest seed string *)
colBits  = If[seedBitStrings =!= None,
  Max[14, Max[StringLength /@ StringSplit[seedBitStrings, ","]] + 2],
  Max[14, StringLength[maxBitString] + 2]];
colState = 26; colN = 9; colSym = 10; colFail = 8; colKL = 12; colNum = 10; colErgo = 14;
totalCols = colBits + colState + colN + colErgo +
  If[mode =!= "Numerical", colSym + colFail, 0] +
  If[mode =!= "Symbolic",  colKL  + colNum,  0];

(* ---- Print header ---- *)
Print[];
Print[StringRepeat["=", totalCols]];
Print["DETAILED BALANCE CHECKER"];
Print["File:          ", algFile];
If[seedBitStrings =!= None,
  Print["SeedBitStrings: ", seedBitStrings,
        "  (", Length[allBitStrings], " seed states)"],
  Print["MaxBitString:  ", maxBitString,
        "  (", Length[allBitStrings], " bit strings generated)"]];
Print["Mode:          ", mode];
If[mode =!= "Symbolic", Print["NSteps:        ", nSteps]];
Print["MaxBitDepth:   ", maxBitDepth];
If[fastChecker && mode =!= "Numerical",
  Print["SymChecker:    FastChecker (Exp-polynomial; FullSimplify fallback)"],
  If[mode =!= "Numerical",
    Print["SymChecker:    FullSimplify (standard)"]]];
Print[StringRepeat["=", totalCols]];
Print[];

Print[StringPadRight["Bits",    colBits],
      StringPadRight["State",   colState],
      StringPadRight["#States", colN],
      StringPadRight["Ergodic", colErgo],
      If[mode =!= "Numerical",
         StringPadRight["Symbolic", colSym] <> StringPadRight["#Fail", colFail], ""],
      If[mode =!= "Symbolic",
         StringPadRight["KL", colKL] <> StringPadRight["Numerical", colNum], ""]];
Print[StringRepeat["-", totalCols]];

(* ================================================================
   Progress bar helpers
   ================================================================ *)
$nBitStrings = Length[allBitStrings];
$checkStart  = AbsoluteTime[];
$barWidth    = 28;

formatDuration[secs_] :=
  Which[secs < 60,   ToString[Round[secs]] <> "s",
        secs < 3600, ToString[Floor[secs/60]] <> "m " <>
                     ToString[Round[Mod[secs, 60]]] <> "s",
        True,        ToString[Floor[secs/3600]] <> "h " <>
                     ToString[Floor[Mod[secs/60, 60]]] <> "m"];

writeProgress[idx_] := If[verbose, Module[{elapsed, pct, filled, arrow, bar, eta, etaStr, elStr},
  elapsed = AbsoluteTime[] - $checkStart;
  pct     = If[$nBitStrings > 0, idx / $nBitStrings, 0.];
  filled  = Round[pct * $barWidth];
  arrow   = If[filled < $barWidth, ">", ""];
  bar     = "[" <> StringRepeat["=", filled] <> arrow <>
            StringRepeat[" ", Max[0, $barWidth - filled - StringLength[arrow]]] <> "]";
  eta     = If[idx > 0 && elapsed > 0,
               elapsed * ($nBitStrings - idx) / idx,
               Infinity];
  etaStr  = If[!NumberQ[eta] || eta === Infinity, "?", formatDuration[eta]];
  elStr   = formatDuration[elapsed];
  WriteString[$Output,
    "\r" <> bar <> " " <> ToString[idx] <> "/" <> ToString[$nBitStrings] <>
    "  " <> ToString[Round[pct * 100]] <> "%  elapsed: " <> elStr <>
    "  ETA: " <> etaStr <> "   "]]];  (* trailing spaces clear previous longer text *)

(* Before printing a permanent row: clear the progress line first *)
clearProgress[] := If[verbose, WriteString[$Output, "\r" <> StringRepeat[" ", 78] <> "\r"]];

(* ================================================================
   Main loop
   ================================================================ *)
coveredAssoc = <||>;   (* Association for O(1) covered-state lookup *)
testedCount  = 0;
symPassCount  = 0; symFailCount  = 0;
numPassCount  = 0; numFailCount  = 0;
ergoPassCount = 0; ergoFailCount = 0;
$loopIdx     = 0;

Do[
  $loopIdx++;
  state = Check[BitsToState[bitString], $Failed];
  If[state === $Failed || state === None,
    writeProgress[$loopIdx]; Continue[]];
  If[KeyExistsQ[coveredAssoc, state],
    writeProgress[$loopIdx]; Continue[]];

  bitsStr  = StringJoin[ToString /@ bitString];
  stateStr = If[Length[DownValues[DisplayState]] > 0,
               DisplayState[state], ToString[state]];

  (* ---- BFS: discover connected component ---- *)
  treeData = BuildTreeAT[state, $alg,
               "MaxBitDepth" -> maxBitDepth,
               "TimeLimit"   -> 120.,
               "Verbose"     -> verbose];

  If[MatchQ[treeData, $dbc$cantHandle[_]],
    clearProgress[];
    Print[StringPadRight[bitsStr,  colBits], StringPadRight[stateStr, colState],
          StringPadRight["?", colN], "ERROR: ", treeData[[1]]];
    writeProgress[$loopIdx];
    testedCount++; Continue[]];

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
  compSymAssumptions = Map[Element[#, Reals] &, Union[compSymVars, $fixedSymVars]];

  (* ---- Symbolic check ---- *)
  symResult = "-"; failStr = "-";
  If[mode =!= "Numerical",
    With[{viol = If[fastChecker,
                  CheckDetailedBalanceFast[matrix, allStates, $energy, compSymAssumptions],
                  CheckDetailedBalance[    matrix, allStates, $energy, compSymAssumptions]]},
      symResult = If[viol === {}, "PASS", "FAIL"];
      failStr   = If[viol === {}, "0",    ToString[Length[viol]]];
      If[symResult === "PASS", symPassCount++, symFailCount++]]];

  (* ---- Numerical MCMC check ---- *)
  klStr = "-"; numResult = "-";
  If[mode =!= "Symbolic" && nStates > 1,
    allBlockVars = Union[compNumericVars, $fixedSymVars];
    $runNumericalMCMC[] := (
      mcmcCounts = RunNumericalMCMCAT[allStates, $alg, $numBeta, "NSteps" -> nSteps];
      boltzmann  = BoltzmannWeightsAT[allStates, $energy, $numBeta]);
    $activateAndRun[] :=
      If[ValueQ[$abstractFunctions] && TrueQ[$abstractFunctions],
        Block[{fieldF, couplingJ},
          If[Length[DownValues[$fieldFConcrete]] > 0,
            fieldF[x_Integer, y_Integer, L_Integer] := $fieldFConcrete[x, y, L]];
          If[Length[DownValues[$couplingJConcrete]] > 0,
            couplingJ[a_Integer, b_Integer, d2_Integer] := $couplingJConcrete[a, b, d2]];
          $runNumericalMCMC[]],
        $runNumericalMCMC[]];
    If[allBlockVars =!= {},
      SeedRandom[42];
      compNumVals = If[compNumericVars =!= {}, RandomReal[{-1, 1}, Length[compNumericVars]], {}];
      Block[Evaluate[allBlockVars],
        If[compNumericVars =!= {}, MapThread[Set, {compNumericVars, compNumVals}]];
        If[$fixedSymVars =!= {}, MapThread[Set, {$fixedSymVars, $fixedNumVals}]];
        $activateAndRun[]],
      $activateAndRun[]];
    simFreq = N[# / Total[mcmcCounts]] & /@ mcmcCounts;
    kl = Total @ Table[
      With[{p = simFreq[s], q = N @ boltzmann[s]},
        If[p > 0 && q > 0, p * Log[p / q], 0.]], {s, allStates}];
    klStr     = Which[Abs[kl] < 5.*^-4, "~0", True, ToString[NumberForm[kl, {4, 3}]]];
    numResult = If[kl < 0.02, "PASS", "FAIL"];
    If[numResult === "PASS", numPassCount++, numFailCount++]];

  (* ---- Print result row ---- *)
  clearProgress[];
  Print[StringPadRight[bitsStr,           colBits],
        StringPadRight[stateStr,          colState],
        StringPadRight[ToString[nStates], colN],
        StringPadRight[ergoStr,           colErgo],
        If[mode =!= "Numerical",
           StringPadRight[symResult, colSym] <> StringPadRight[failStr, colFail], ""],
        If[mode =!= "Symbolic",
           StringPadRight[klStr, colKL] <> StringPadRight[numResult, colNum], ""]];

  writeProgress[$loopIdx],

  {bitString, allBitStrings}];

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
If[ergoFailCount > 0,
  Print["*** ERGODICITY VIOLATION DETECTED ***"]];
If[symFailCount == 0 && numFailCount == 0 && testedCount > 0,
  Print["All tested systems satisfy detailed balance."],
  If[testedCount == 0,
    Print["No valid states found for the given MaxBitString."],
    Print["*** DETAILED BALANCE VIOLATION DETECTED ***"]]];
Print[StringRepeat["=", totalCols]];
