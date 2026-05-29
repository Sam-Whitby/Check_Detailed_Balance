#!/usr/bin/env wolframscript
(* show_expressions.wls
   Print raw DB expressions (before FullSimplify) for vmmc_lattice at NGrid=3.
   Shows 72-state (2-particle) and 504-state (3-particle) components. *)

$dir = DirectoryName[$InputFileName];
If[$dir === "", $dir = NotebookDirectory[]];

Get[$dir <> "dbc_core.wl"];
Get[$dir <> "examples3/vmmc_lattice.wl"];

LaunchKernels[4];

(* ---- Abstract-param setup ---- *)
$abstractParamRestore = Association @ Flatten @ Map[
  Function[name,
    With[{val = Symbol[name]},
      If[NumericQ[val], {name -> val}, {}]]],
  $checkerAbstractParams];
Scan[ClearAll, $checkerAbstractParams];
$abstractParamAssumptions = Flatten @ KeyValueMap[
  Function[{name, val},
    With[{sym = Symbol[name]},
      Which[val > 0, {Element[sym, Reals], sym > 0},
            val < 0, {Element[sym, Reals], sym < 0},
            True, {Element[sym, Reals]}]]],
  $abstractParamRestore];

$dbcDistributeBFS[];

tex[e_] := ToString[TeXForm[e]];   (* render to LaTeX string *)

(* ---- Analyse one component ---- *)
analyseComponent[seedState_, label_] := Module[
  {treeData, allStates, tMatrix, dynP, dynVars, assm,
   n, pairs, nonTrivPairs, allExprs, caseHist,
   simple, medianE, complex, sortedByCases},

  Print["\n", StringRepeat["=", 70]];
  Print["COMPONENT: ", label];
  Print["Seed: ", seedState];

  $bfsAlg = Function[s,
    Block[{$vmmcCandidates = $dbcCanonicalCandidates}, Algorithm[s]]];

  {bfsTime, treeData} = AbsoluteTiming[
    BuildTreeAT[seedState, $bfsAlg,
      "MaxBitDepth" -> 20, "TimeLimit" -> 120, "Verbose" -> False]];
  allStates = Sort @ Keys[treeData];
  tMatrix   = TreeATToMatrix[treeData];
  Print["BFS: ", Length[allStates], " states, ", Length[tMatrix],
        " T-entries  (", bfsTime, " s)"];

  dynP    = DynamicSymParams[allStates];
  dynVars = Join[Lookup[dynP, "eps", {}], Lookup[dynP, "couplings", {}]];
  assm    = Join[{\[Beta] > 0}, Map[Element[#, Reals] &, dynVars],
                  $abstractParamAssumptions];
  Print["Coupling atoms: ", dynVars];

  n = Length[allStates];
  pairs = Flatten[Table[{i,j}, {i,1,n}, {j,i+1,n}], 1];
  nonTrivPairs = Select[pairs, Function[p,
    Lookup[tMatrix, Key[{allStates[[p[[1]]]], allStates[[p[[2]]]]}], 0] =!= 0 ||
    Lookup[tMatrix, Key[{allStates[[p[[2]]]], allStates[[p[[1]]]]}], 0] =!= 0]];
  Print["Non-trivial pairs: ", Length[nonTrivPairs]];

  (* Check: does FastChecker (Expand===0) succeed on the raw expression? *)
  Print["Computing expressions and fast-checker hit rates..."];
  allExprs = Map[Function[p,
    Module[{si, sj, tij, tji, ei, ej, raw, pw, nC, fastHit},
      si  = allStates[[p[[1]]]]; sj = allStates[[p[[2]]]];
      tij = Lookup[tMatrix, Key[{si, sj}], 0];
      tji = Lookup[tMatrix, Key[{sj, si}], 0];
      ei  = energy[si] /. r_Real :> Rationalize[r];
      ej  = energy[sj] /. r_Real :> Rationalize[r];
      raw = tij * Exp[-\[Beta] * ei] - tji * Exp[-\[Beta] * ej];
      pw  = PiecewiseExpand[raw];
      nC  = If[MatchQ[pw, Piecewise[__]], Length[pw[[1]]], 1];
      (* Quick check: can Expand alone prove it zero in every case? *)
      fastHit = If[MatchQ[pw, Piecewise[__]],
        AllTrue[pw[[1]], Expand[#[[1]]] === 0 &] && Expand[pw[[2]]] === 0,
        Expand[pw] === 0];
      <|"p"->p, "si"->si, "sj"->sj, "tij"->tij, "tji"->tji,
        "ei"->ei, "ej"->ej, "raw"->raw, "pw"->pw,
        "nCases"->nC, "expandZero"->fastHit|>]],
    nonTrivPairs];

  caseHist = Sort @ Tally[allExprs[[All, "nCases"]]];
  Print["Piecewise case-count distribution: ", caseHist];
  Print["Expand===0 (trivially zero): ",
        Count[allExprs, x_ /; x["expandZero"]], "/", Length[allExprs]];

  (* ---- Print one representative multi-case expression ---- *)
  sortedByCases = SortBy[allExprs, #["nCases"] &];
  complex = Last[sortedByCases];

  Print["\n--- Most complex expression (", complex["nCases"], " Piecewise cases) ---"];
  printOne[complex];
]

printOne[rec_] := Module[{pw = rec["pw"]},
  Print["  si = ", rec["si"]];
  Print["  sj = ", rec["sj"]];
  Print["  E(si) = $", tex[rec["ei"]], "$"];
  Print["  E(sj) = $", tex[rec["ej"]], "$"];
  Print["  T(si->sj) = $", tex[rec["tij"]], "$"];
  Print["  T(sj->si) = $", tex[rec["tji"]], "$"];
  Print[""];
  If[MatchQ[pw, Piecewise[__]],
    Do[
      Print["  Case ", c, " condition: $", tex[pw[[1,c,2]]], "$"];
      Print["  Case ", c, " value:     $", tex[pw[[1,c,1]]], "$"];
      Print[],
      {c, Length[pw[[1]]]}];
    Print["  Default: $", tex[pw[[2]]], "$"],
    Print["  (scalar) $", tex[pw], "$"]]]

(* ---- Run for both components ---- *)
analyseComponent[{1, 2, 0, 0, 0, 0, 0, 0, 0}, "72-state (2-particle)"];
analyseComponent[{1, 2, 3, 0, 0, 0, 0, 0, 0}, "504-state (3-particle)"];

Print["\nDone."];
