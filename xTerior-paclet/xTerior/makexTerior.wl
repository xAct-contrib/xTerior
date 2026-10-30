(* ::Package:: *)

(* Build the authored documentation before installing the paclet. *)
Needs["PacletTools`"];
With[{result = PacletBuild[DirectoryName[$InputFileName]]},
  If[MatchQ[result, _Success],
    xTerior = result["PacletArchive"];
    PacletInstall[xTerior, ForceVersionInstall -> True],
    result
  ]
]
