#!/usr/bin/env wolframscript
(* ================================================================
   report.wls  --  Generate a PNG visualisation report for a
                   single bit-string seed of any .wl algorithm.
   ================================================================
   Usage:
     wolframscript -file report.wls <alg_file> BitString=XXXXX [options]

   Options (key=value, no spaces around =):
     BitString=XXXXX        Bit string seed (required, e.g. BitString=11101)
     NSteps=N               MCMC steps for numerical check (default 50000)
     MaxBitDepth=N          BFS tree depth cap (default 20)
     Verbose=True/False     Print per-state BFS progress (default False)
     ShowTrees=True/False   Include decision trees in report (default True)
     ShowT=True/False       Include transition matrix in report (default True)
     DoNumerical=True/False Run numerical MCMC and include frequency plot (default True)
     PrintMatrix=True/False Print transition matrix as LaTeX to terminal (default False)
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* Launch parallel kernels for FullSimplify parallelisation *)
$nKernels = 0;
Quiet[Check[LaunchKernels[]; $nKernels = Length[Kernels[]], $nKernels = 0]];
If[$nKernels > 0,
  Print["Parallel: ", $nKernels, " kernels launched."],
  Print["Note: parallel kernel launch failed — running single-threaded."]];

(* ---- Parse command-line arguments ---- *)
args = $ScriptCommandLine;
If[Length[args] >= 1 && StringMatchQ[args[[1]], ___ ~~ "report.wls"],
  args = Rest[args]];

posArgs = Select[args, !StringMatchQ[#, __ ~~ "=" ~~ __] &];
kvArgs  = Select[args,  StringMatchQ[#, __ ~~ "=" ~~ __] &];

If[Length[posArgs] < 1,
  Print["Usage: wolframscript -file report.wls <alg_file> BitString=XXXXX [options]"];
  Print[""];
  Print["Required definitions in the .wl file:"];
  Print["  energy[state_]       bare energy (no beta)"];
  Print["  Algorithm[state_]    MCMC move function"];
  Print["  BitsToState[bits_]   {0,1,...} -> state, or None if invalid"];
  Print["  numBeta              numeric inverse temperature"];
  Print[""];
  Print["Options (key=value, no spaces around =):"];
  Print["  BitString=XXXXX        Bit string seed (required)"];
  Print["  NSteps=N               MCMC steps (default 50000)"];
  Print["  MaxBitDepth=N          BFS depth cap (default 20)"];
  Print["  Verbose=True/False     Per-state BFS output (default False)"];
  Print["  ShowTrees=True/False   Include decision trees (default True)"];
  Print["  ShowT=True/False       Include transition matrix (default True)"];
  Print["  DoNumerical=True/False Run numerical MCMC + frequency plot (default True)"];
  Print["  PrintMatrix=True/False Print transition matrix as LaTeX to terminal (default False)"];
  Exit[1]];

(* ---- Option defaults ---- *)
bitStringArg = None;
nSteps       = 50000;
maxBitDepth  = 20;
verbose      = False;
showTrees    = True;
showT        = True;
doNumerical  = True;
printMatrix  = False;

Do[
  If[StringMatchQ[opt, "BitString="   ~~ __], bitStringArg = StringDrop[opt, 10]];
  If[StringMatchQ[opt, "NSteps="      ~~ __], nSteps       = ToExpression[StringDrop[opt,  7]]];
  If[StringMatchQ[opt, "MaxBitDepth=" ~~ __], maxBitDepth  = ToExpression[StringDrop[opt, 12]]];
  If[StringMatchQ[opt, "Verbose="     ~~ __], verbose      = ToExpression[StringDrop[opt,  8]]];
  If[StringMatchQ[opt, "ShowTrees="   ~~ __], showTrees    = ToExpression[StringDrop[opt, 10]]];
  If[StringMatchQ[opt, "ShowT="       ~~ __], showT        = ToExpression[StringDrop[opt,  6]]];
  If[StringMatchQ[opt, "DoNumerical=" ~~ __], doNumerical  = ToExpression[StringDrop[opt, 12]]];
  If[StringMatchQ[opt, "PrintMatrix=" ~~ __], printMatrix  = ToExpression[StringDrop[opt, 12]]],
  {opt, kvArgs}];

If[bitStringArg === None,
  Print["ERROR: BitString argument is required (e.g. BitString=11101)"];
  Exit[1]];

If[!StringMatchQ[bitStringArg, ("0" | "1")..],
  Print["ERROR: BitString must contain only 0s and 1s, got: ", bitStringArg];
  Exit[1]];

bitString = ToExpression /@ Characters[bitStringArg];

(* ---- Resolve and load algorithm file ---- *)
algFile = posArgs[[1]];
If[!FileExistsQ[algFile], algFile = scriptDir <> algFile];
If[!FileExistsQ[algFile],
  Print["ERROR: File not found: ", posArgs[[1]]]; Exit[1]];

ClearAll[Algorithm, energy, BitsToState, numBeta, symParams, DynamicSymParams, DisplayState,
         fieldF, couplingJ,
         $fieldFConcrete, $couplingJConcrete, $abstractFunctions, $concreteParams,
         $fieldFreeParams, $couplingFreeParams];
Print["Loading: ", algFile];
Get[algFile];

If[Length[DownValues[Algorithm]] == 0,
  Print["ERROR: 'Algorithm' not defined in ", posArgs[[1]]]; Exit[1]];
If[Length[DownValues[energy]] == 0,
  Print["ERROR: 'energy' not defined in ", posArgs[[1]]]; Exit[1]];
If[Length[DownValues[BitsToState]] == 0,
  Print["ERROR: 'BitsToState' not defined in ", posArgs[[1]]]; Exit[1]];
If[!ValueQ[numBeta] || !NumericQ[numBeta],
  Print["ERROR: 'numBeta' (numeric) not defined in ", posArgs[[1]]]; Exit[1]];

$alg     = Algorithm;
$energy  = energy;
$numBeta = numBeta;

(* ---- Verify couplingJ symmetry (if defined) ---- *)
(* In abstract-function mode ($abstractFunctions=True), couplingJ has no
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

(* ---- Map bit string to seed state ---- *)
state = Check[BitsToState[bitString], $Failed];
If[state === $Failed || state === None,
  Print["ERROR: BitString ", bitStringArg,
        " does not map to a valid state (BitsToState returned None or failed)."];
  Exit[1]];

Print["Seed state: ", state];

(* ---- BFS: discover connected component ---- *)
Print["Building decision trees (BFS) ..."];
treeData = BuildTreeAT[state, $alg,
             "MaxBitDepth" -> maxBitDepth,
             "TimeLimit"   -> 300.,
             "Verbose"     -> verbose];

If[MatchQ[treeData, $dbc$cantHandle[_]],
  Print["ERROR: BFS failed -- ", treeData[[1]]]; Exit[1]];

allStates = Sort @ Keys[treeData];
matrix    = TreeATToMatrix[treeData];
nStates   = Length[allStates];
Print["States in connected component: ", nStates, "  --  ", allStates];

(* ---- Per-component symbolic parameters ---- *)
$hasDynParams    = Length[DownValues[DynamicSymParams]] > 0;
$staticSysParams = If[ValueQ[symParams], symParams, None];
If[AssociationQ[$staticSysParams],
  $staticSymVars = Join[Lookup[$staticSysParams, "eps",       {}],
                        Lookup[$staticSysParams, "couplings", {}]],
  $staticSymVars = {}];

(* Fixed symbolic parameters: kept free during symbolic check, assigned their
   declared numeric values (not random) during numerical MCMC. *)
$fixedParams  = If[AssociationQ[$staticSysParams],
                   Lookup[$staticSysParams, "fixedParams", <||>], <||>];
$fixedSymVars = Keys[$fixedParams];
$fixedNumVals = Values[$fixedParams];

If[$hasDynParams,
  dynP          = DynamicSymParams[allStates];
  dynVars       = Join[Lookup[dynP, "eps", {}], Lookup[dynP, "couplings", {}]];
  dynNumVars    = Lookup[dynP, "numericParams", dynVars];
  compSymVars    = Union[$staticSymVars, dynVars];
  compNumericVars = Union[$staticSymVars, dynNumVars],
  compSymVars    = $staticSymVars;
  compNumericVars = $staticSymVars];
compSymAssumptions = Map[Element[#, Reals] &, Union[compSymVars, $fixedSymVars]];

(* ---- Symbolic detailed-balance check ---- *)
Print["Running symbolic detailed balance check ..."];
violations = CheckDetailedBalance[matrix, allStates, $energy, compSymAssumptions];
If[violations === {},
  Print["Symbolic: PASS"],
  Print["Symbolic: FAIL -- ", Length[violations], " violation(s) found."]];

(* ---- Print transition matrix as LaTeX (if PrintMatrix=True) ---- *)
If[printMatrix,
  Print["\n--- Transition matrix T[i -> j]  (LaTeX) ---"];
  Print["States: ", allStates];
  Print[];
  Do[
    Do[
      With[{entry = Lookup[matrix, Key[{si, sj}], 0]},
        If[entry =!= 0,
          Print["T[", si, " -> ", sj, "] ="];
          Print[ToString[TeXForm[entry]]];
          Print[]]],
      {sj, allStates}],
    {si, allStates}];
  Print["--- end of transition matrix ---\n"]];

(* ---- Numerical MCMC check (skipped if DoNumerical=False) ---- *)
If[doNumerical,
  Print["Running numerical MCMC (", nSteps, " steps) ..."];
  allBlockVars = Union[compNumericVars, $fixedSymVars];
  $runNumericalMCMC[] := (
    mcmcCounts = RunNumericalMCMCAT[allStates, $alg, $numBeta, "NSteps" -> nSteps];
    bw = BoltzmannWeightsAT[allStates, $energy, $numBeta]);
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
    With[{p = simFreq[s], q = N @ bw[s]},
      If[p > 0 && q > 0, p * Log[p / q], 0.]], {s, allStates}];
  Print["KL divergence: ", NumberForm[kl, {5, 4}]];
  Print["Numerical: ", If[kl < 0.02, "PASS", "FAIL"]],
  (* DoNumerical=False: set sentinels so reportArgs is well-formed *)
  simFreq = None; bw = None; kl = None];

(* ---- Generate notebook report ---- *)
algCode = ReadString[algFile];
sysName = FileBaseName[algFile] <> "  bits=" <> bitStringArg;

reportArgs = <|
  "name"        -> sysName,
  "allStates"   -> allStates,
  "treeData"    -> treeData,
  "matrix"      -> matrix,
  "violations"  -> violations,
  "simFreq"     -> simFreq,
  "bw"          -> bw,
  "kl"          -> kl,
  "algCode"     -> algCode,
  "showTrees"   -> showTrees,
  "showT"       -> showT,
  "doNumerical" -> doNumerical|>;

Print["\nGenerating notebook report ..."];
ExportReportNotebook[reportArgs];
