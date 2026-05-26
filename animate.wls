#!/usr/bin/env wolframscript
(* ================================================================
   animate.wls  --  Run an algorithm and animate it via Python/matplotlib
   ================================================================
   Usage:
     wolframscript -file animate.wls <alg_file> \
       Sites=<n> N=<n> [Steps=200] [Beta=1] [FPS=10] \
       [Simple=1] [RecordEvery=1] [J12=0.5 ...]

   Runs the algorithm for the specified number of steps, collects
   the state and energy at each recorded step, writes the data to
   a temporary JSON file, then launches animate_plot.py to display
   an animated matplotlib figure with:
     - Left panel:  imshow of the lattice (colour = particle type)
     - Right panel: energy vs step (growing line)

   For 2D algorithms (those defining DisplayState), the lattice is
   shown as an L x L grid. For 1D, it is shown as a 1 x Sites row.

   Arguments:
     Sites=<n>          Total lattice sites (e.g. 6 for a 1D ring,
                        9 for a 3x3 grid).
     N=<n>              Number of labeled particles (types 1..N).
     Steps=<n>          MCMC steps to run (default 200).
     Beta=<f>           Inverse temperature (overrides numBeta in file).
     FPS=<f>            Animation frame rate in frames per second (default 10).
     Simple=1           Use fast 2-colour mode: holes vs particles.
                        Renders much faster, enabling high frame rates.
     RecordEvery=<n>    Record state every n steps (default 1).
     J<a><b>=<f>        Set coupling constant J_ab (random if omitted).

   Examples:
     wolframscript -file animate.wls examples3/kawasaki_1d.wl \
       Sites=6 N=3 Steps=500 FPS=30 Simple=1

     wolframscript -file animate.wls examples3/kawasaki_2d.wl \
       Sites=9 N=4 Steps=300 FPS=10 J12=1.5

     wolframscript -file animate.wls examples3/cluster_2d.wl \
       Sites=9 N=3 Steps=400 FPS=15 RecordEvery=2
   ================================================================ *)

(* Real MetropolisProb for actual execution (not BFS interception).
   Clamp large positive exponents to avoid underflow warnings. *)
MetropolisProb[dE_] := If[dE <= 0, 1., If[$animBeta * dE >= 700., 0., Exp[-$animBeta * dE]]]

(* $jPairSym: symmetric coupling symbol Jpair<lo><hi>.
   Defined here (without loading dbc_core.wl) so that energy[] evaluates
   numerically during animation runs. *)
$jPairSym[a_Integer, b_Integer] :=
  ToExpression["Jpair" <> ToString[Min[a, b]] <> ToString[Max[a, b]]]

Off[General::munfl]

scriptDir = DirectoryName[$InputFileName];

(* ---- Parse arguments ---- *)
args = $ScriptCommandLine;
If[Length[args] >= 1 && StringMatchQ[args[[1]], ___ ~~ "animate.wls"],
  args = Rest[args]];

posArgs = Select[args, !StringMatchQ[#, __ ~~ "=" ~~ __] &];
kvArgs  = Select[args,  StringMatchQ[#, __ ~~ "=" ~~ __] &];

If[Length[posArgs] < 1 ||
   !MemberQ[kvArgs, _?(StringMatchQ[#, "Sites=" ~~ __] &)] ||
   !MemberQ[kvArgs, _?(StringMatchQ[#, "N=" ~~ __] &)],
  Print["Usage: wolframscript -file animate.wls <alg_file> Sites=<n> N=<n> [options]"];
  Print["  Sites=<n>        total lattice sites"];
  Print["  N=<n>            number of labeled particles (types 1..N)"];
  Print["  Steps=<n>        MCMC steps (default 200)"];
  Print["  Beta=<f>         inverse temperature (default from .wl file)"];
  Print["  FPS=<f>          animation frame rate in frames per second (default 10)"];
  Print["  Simple=1         fast 2-colour mode: holes vs particles"];
  Print["  RecordEvery=<n>  record every nth step (default 1)"];
  Print["  J<a><b>=<f>      coupling constant (random if omitted)"];
  Exit[1]];

(* ---- Defaults ---- *)
nSites      = 4;
nPart       = 2;
nSteps      = 200;
animBeta    = None;
fps         = 10.;
simpleMode  = False;
recordEvery = 1;
noParams    = False;

Do[
  If[StringMatchQ[opt, "Sites="       ~~ __], nSites      = ToExpression[StringDrop[opt,  6]]];
  If[StringMatchQ[opt, "N="           ~~ __], nPart       = ToExpression[StringDrop[opt,  2]]];
  If[StringMatchQ[opt, "Steps="       ~~ __], nSteps      = ToExpression[StringDrop[opt,  6]]];
  If[StringMatchQ[opt, "Beta="        ~~ __], animBeta    = ToExpression[StringDrop[opt,  5]]];
  If[StringMatchQ[opt, "FPS="         ~~ __], fps         = N[ToExpression[StringDrop[opt,  4]]]];
  (* Legacy Delay= support: convert seconds -> fps *)
  If[StringMatchQ[opt, "Delay="       ~~ __], fps         = N[1 / ToExpression[StringDrop[opt, 6]]]];
  If[StringMatchQ[opt, "Simple="      ~~ __], simpleMode  = (StringDrop[opt, 7] =!= "0")];
  If[StringMatchQ[opt, "RecordEvery=" ~~ __], recordEvery = ToExpression[StringDrop[opt, 12]]];
  If[StringMatchQ[opt, "NoParams="    ~~ __], noParams    = (StringDrop[opt, 9] =!= "0")];
  If[StringMatchQ[opt, "NoParams"],           noParams    = True],
  {opt, kvArgs}];

(* ---- Load algorithm file ---- *)
algFile = posArgs[[1]];
If[!FileExistsQ[algFile], algFile = scriptDir <> algFile];
If[!FileExistsQ[algFile],
  Print["ERROR: File not found: ", posArgs[[1]]]; Exit[1]];

ClearAll[Algorithm, energy, BitsToState, numBeta, symParams, DynamicSymParams, DisplayState,
         AnimationSetupCouplings,
         fieldF, couplingJ,
         $fieldFConcrete, $couplingJConcrete, $abstractFunctions, $concreteParams,
         $fieldFreeParams, $couplingFreeParams];
Get[algFile];

(* ---- Activate concrete function implementations (abstract-function mode) ----
   In files where fieldF and couplingJ have no DownValues ($abstractFunctions=True),
   permanently install the concrete implementations for this animation run. *)
If[ValueQ[$abstractFunctions] && TrueQ[$abstractFunctions],
  If[Length[DownValues[$fieldFConcrete]] > 0,
    fieldF[x_Integer, y_Integer, L_Integer] := $fieldFConcrete[x, y, L]];
  If[Length[DownValues[$couplingJConcrete]] > 0,
    couplingJ[a_Integer, b_Integer, d2_Integer] := $couplingJConcrete[a, b, d2]]];

If[Length[DownValues[Algorithm]] == 0,
  Print["ERROR: 'Algorithm' not defined in ", algFile]; Exit[1]];

$animBeta = If[animBeta =!= None, animBeta, numBeta];
(* Assign \[Beta] (Unicode β) so that VMMC link weights — which use \[Beta]
   directly in Piecewise expressions — receive the correct numeric value. *)
\[Beta] = $animBeta;

(* ---- Build initial state ---- *)
If[nPart > nSites,
  Print["ERROR: N=", nPart, " exceeds Sites=", nSites]; Exit[1]];
SeedRandom[1234];
initialState = RandomSample[Join[Range[nPart], ConstantArray[0, nSites - nPart]], nSites];

(* ---- Determine display grid shape ---- *)
(* 2D if DisplayState is defined (convention for 2D algorithm files) *)
is2D   = Length[DownValues[DisplayState]] > 0;
gridL  = If[is2D, Round[Sqrt[nSites]], nSites];
nRows  = If[is2D, gridL, 1];
nCols  = If[is2D, gridL, nSites];

(* ---- Assign coupling constants ---- *)
(* 1. Parse ALL J-prefixed key=val args into a string-keyed Association.
      Accepts both simple J12=-1.0 and extended Jd112=-1.0, Jdsq212=-0.5, etc. *)
jAllSupplied = Association @ Map[
  Function[opt, With[{p = StringSplit[opt, "=", 2]}, p[[1]] -> ToExpression[p[[2]]]]],
  Select[kvArgs, StringMatchQ[#,
    "J" ~~ (LetterCharacter | DigitCharacter).. ~~ "=" ~~ __] &]];

(* 2. Build symbol-keyed jSupplied for backward-compat Lookup (keys are Mathematica symbols) *)
jSupplied = Association @ KeyValueMap[Rule[ToExpression[#1], #2] &, jAllSupplied];

types = Sort[DeleteCases[initialState, 0]];

(* 3. Pre-collect DynamicSymParams symbol names as strings BEFORE any assignments.
      At this point all symbols are unbound, so ToString gives their names correctly.
      Must happen before steps 4-6 which set symbol values.
      For abstract-function files ($abstractFunctions=True), DynamicSymParams returns
      function-application atoms (like fieldF[0,0,2]) under "couplings" — these cannot
      be assigned random numeric values directly.  Instead we use the "numericParams"
      key, which holds the scalar symbols that need values (Jpair12, fieldAmp, etc.). *)
$dynSymNames = {};
If[Length[DownValues[DynamicSymParams]] > 0,
  With[{dynP = DynamicSymParams[{initialState}]},
    $dynSymNames = Map[ToString,
      If[ValueQ[$abstractFunctions] && TrueQ[$abstractFunctions],
         Join[Lookup[dynP, "eps", {}], Lookup[dynP, "numericParams", {}]],
         Join[Lookup[dynP, "eps", {}], Lookup[dynP, "couplings",     {}]]]]]];

(* 4. Set all explicitly supplied J symbols immediately *)
KeyValueMap[Function[{name, val},
  With[{key = name}, Set[Evaluate[ToExpression[key]], val]]], jAllSupplied];

(* 5. Auto-assign any missing J<a><b> with random values (backward compat for simple files).
      Work with string keys throughout; only call ToExpression[key] for Set. *)
SeedRandom[42];
jAssigned = jAllSupplied;   (* start from explicitly supplied; will add auto-assigned ones *)
Do[If[a < b,
  With[{key = "J" <> ToString[a] <> ToString[b]},
    If[!KeyExistsQ[jAllSupplied, key],
      With[{val = RandomVariate[NormalDistribution[0, 1]]},
        Set[Evaluate[ToExpression[key]], val];
        AssociateTo[jAssigned, key -> val]]]]],
  {a, types}, {b, types}];

(* 6. Algorithm-specific coupling hook: fills in any missing extended-family symbols.
      Define AnimationSetupCouplings[jArgs, types] in the .wl file to handle e.g.
      Jd1<a><b>, Jdsq2<a><b>, Jd2<a><b>, Jdsq5<a><b>. *)
If[Length[DownValues[AnimationSetupCouplings]] > 0,
  AnimationSetupCouplings[jAllSupplied, types]];

(* 7. Assign any remaining symbols from DynamicSymParams not yet set.
      Uses $dynSymNames captured in step 3 (before any assignments), so names are
      always correct strings regardless of whether symbols now have numeric values.
      Rebuilds jAssigned so the params table shows all declared symbols. *)
If[$dynSymNames =!= {},
  SeedRandom[99];
  Scan[Function[name,
    With[{sym = ToExpression[name]},
      If[!NumericQ[sym],
        Set[Evaluate[sym], RandomVariate[NormalDistribution[0, 1]]]]]],
    $dynSymNames];
  jAssigned = Association @ Map[
    Function[name, name -> N[ToExpression[name]]],
    $dynSymNames]];
(* For files without DynamicSymParams, jAssigned has correct entries from step 5. *)

(* ---- Run algorithm and collect data ---- *)
Print["Running ", nSteps, " steps of ", algFile, " ..."];
state = initialState;
stepsRec    = {0};
statesRec   = {Normal[state]};      (* ensure plain list, not packed *)
energiesRec = {N[energy[state]]};

Do[
  state = Algorithm[state];
  If[Mod[i, recordEvery] == 0,
    AppendTo[stepsRec,    i];
    AppendTo[statesRec,   Normal[state]];
    AppendTo[energiesRec, N[energy[state]]]],
  {i, 1, nSteps}];

Print["Recorded ", Length[stepsRec], " frames. Launching plot..."];

(* ---- Write JSON data file ---- *)
dataFile = FileNameJoin[{$TemporaryDirectory, "dbc_animate_data.json"}];
Export[dataFile,
  <|"algo_file"         -> algFile,
    "grid_rows"         -> nRows,
    "grid_cols"         -> nCols,
    "n_types"           -> If[types === {}, 0, Max[types]],
    "fps"               -> fps,
    "simple"            -> simpleMode,
    "steps"             -> stepsRec,
    "states"            -> statesRec,
    "energies"          -> energiesRec,
    "params"            -> Join[<|"beta" -> $animBeta|>, jAssigned],
    "coupling_formula"  -> If[ValueQ[$couplingFormulaStr], $couplingFormulaStr, ""],
    "field_formula"     -> If[ValueQ[$fieldFormulaStr],    $fieldFormulaStr,    ""]|>];

(* ---- Launch Python animation ---- *)
pyScript = FileNameJoin[{scriptDir, "animate_plot.py"}];
If[!FileExistsQ[pyScript],
  Print["ERROR: animate_plot.py not found at ", pyScript]; Exit[1]];

(* Find a python3 that has matplotlib; prefer python3.11, fall back to python3 *)
pythonCmd = SelectFirst[{"python3.11", "python3", "python"},
  Function[cmd, RunProcess[{cmd, "-c", "import matplotlib"}, "ExitCode"] === 0],
  "python3"];

result = RunProcess[Join[{pythonCmd, pyScript, dataFile},
                        If[TrueQ[noParams], {"--no-params"}, {}]]];
If[result["ExitCode"] != 0,
  (* Only print non-empty streams; a non-zero exit with empty output usually
     means the user simply closed the animation window. *)
  If[StringLength[result["StandardError"]] > 0,
    Print["Python error:\n", result["StandardError"]]];
  If[StringLength[result["StandardOutput"]] > 0,
    Print["Standard output:\n", result["StandardOutput"]]]];
