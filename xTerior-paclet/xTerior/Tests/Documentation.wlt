(* Static documentation checks; no xAct session is needed. *)
docPacletRoot = DirectoryName[DirectoryName[$TestFileName]];
docSourceFiles = FileNames["*.wl", FileNameJoin[{docPacletRoot, "Kernel"}], Infinity];
docExpectedNames = Sort@DeleteDuplicates@Flatten[
	StringCases[Import[#, "Text"],
		RegularExpression["([A-Za-z$][A-Za-z0-9$]*)::usage\\s*="] -> "$1"] & /@ docSourceFiles];
docSymbolFiles = FileNames["*.nb", FileNameJoin[{docPacletRoot, "Documentation", "English", "ReferencePages", "Symbols"}]];
docGuideFile = FileNameJoin[{docPacletRoot, "Documentation", "English", "Guides", "xTerior.nb"}];
docNotebooks = Get /@ Append[docSymbolFiles, docGuideFile];
docKnownURIs = Join[("xTerior/ref/" <> #) & /@ docExpectedNames, {"xTerior/guide/xTerior"}];
docLinks = DeleteDuplicates@Cases[docNotebooks,
	HoldPattern[ButtonData -> uri_String] /; StringStartsQ[uri, "paclet:xTerior/"] :> StringDrop[uri, 7], Infinity];

VerificationTest[Sort[FileBaseName /@ docSymbolFiles], docExpectedNames,
	TestID -> "Every explicit usage message has exactly one symbol page"]

VerificationTest[And @@ (MatchQ[#, _Notebook] & /@ docNotebooks), True,
	TestID -> "All documentation files are notebook expressions"]

VerificationTest[Complement[docLinks, docKnownURIs], {},
	TestID -> "All local paclet links resolve to authored pages"]

VerificationTest[And @@ Table[
	With[{nb = Get[file], name = FileBaseName[file]},
		!FreeQ[nb, Cell[name, "ObjectName", ___]] &&
		!FreeQ[nb, Cell[_, "Usage", ___]] &&
		!FreeQ[nb, Cell[_, "ExampleInitialization", ___]] &&
		!FreeQ[nb, Cell[_, "Input", ___]] &&
		!FreeQ[nb, Cell["xTerior/ref/" <> name, "Categorization", ___]]],
	{file, docSymbolFiles}], True,
	TestID -> "Each symbol has title usage initialization examples and URI"]

VerificationTest[And @@ (DuplicateFreeQ[Cases[#, HoldPattern[CellID -> id_] :> id, Infinity]] & /@ docNotebooks), True,
	TestID -> "Cell identifiers are unique within every notebook"]

VerificationTest[FreeQ[docNotebooks, "XXXX" | "TODO" | $Failed | _MakeBoxes], True,
	TestID -> "No template placeholders failed expressions or unresolved box generation"]

VerificationTest[And @@ (MatchQ[
	Quiet[Check[ToExpression[#, StandardForm, HoldComplete], $Failed]], _HoldComplete] & /@
	Cases[docNotebooks, Cell[BoxData[boxes_], "Input" | "ExampleInitialization", ___] :> boxes, Infinity]), True,
	TestID -> "Every example and initialization cell contains parseable input boxes"]

VerificationTest[FreeQ[docNotebooks, BoxData[_HoldComplete]], True,
	TestID -> "No unevaluated parser wrappers remain inside notebook boxes"]

VerificationTest[DeleteDuplicates@Cases[
	Cases[docNotebooks, cell:Cell[_, "Output", ___] :> cell, Infinity],
	s_Symbol /; MemberQ[docExpectedNames, SymbolName[s]] &&
		!MemberQ[{"xAct`xTerior`", "xAct`FormsToTensors`"}, Context[s]] :> HoldForm[s], Infinity], {},
	TestID -> "Cached output preserves the contexts of documented package symbols"]

VerificationTest[Sort@Cases[docNotebooks, Cell[uri_String, "Categorization", CellLabel -> "URI", ___] :> uri, Infinity], Sort[docKnownURIs],
	TestID -> "Documentation metadata supplies each expected URI exactly once"]
