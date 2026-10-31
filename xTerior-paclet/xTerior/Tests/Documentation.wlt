(* Static documentation checks; no xAct session is needed. *)
docPacletRoot = DirectoryName[DirectoryName[$TestFileName]];
docSourceFiles = FileNames["*.wl", FileNameJoin[{docPacletRoot, "Kernel"}], Infinity];
docExpectedNames = Sort@DeleteDuplicates@Flatten[
	StringCases[Import[#, "Text"],
		RegularExpression["([A-Za-z$][A-Za-z0-9$]*)::usage\\s*="] -> "$1"] & /@ docSourceFiles];
docSymbolFiles = FileNames["*.nb", FileNameJoin[{docPacletRoot, "Documentation", "English", "ReferencePages", "Symbols"}]];
docGuideFile = FileNameJoin[{docPacletRoot, "Documentation", "English", "Guides", "xTerior.nb"}];
docTutorialDirectory = FileNameJoin[{docPacletRoot, "Documentation", "English", "Tutorials"}];
docTutorialFiles = FileNames["*.nb", docTutorialDirectory];
docNotebooks = Get /@ Join[docSymbolFiles, {docGuideFile}, docTutorialFiles];
docKnownURIs = Join[("xTerior/ref/" <> #) & /@ docExpectedNames, {"xTerior/guide/xTerior"},
	("xTerior/tutorial/" <> FileBaseName[#]) & /@ docTutorialFiles];
docLinks = DeleteDuplicates@Cases[docNotebooks,
	HoldPattern[ButtonData -> uri_String] /; StringStartsQ[uri, "paclet:xTerior/"] :> First[StringSplit[StringDrop[uri, 7], "#"]], Infinity];
docHeldInputQ[_HoldComplete] := True;
docHeldInputQ[Null] := True;
docHeldInputQ[list_List] := Length[list] > 0 && And @@ (docHeldInputQ /@ list);
docHeldInputQ[_] := False;
docParseInput[BoxData[boxes_]] := Quiet[Check[ToExpression[boxes, StandardForm, HoldComplete], $Failed]];
docParseInput[text_String] := Quiet[Check[ToExpression[text, InputForm, HoldComplete], $Failed]];
docParseInput[_] := $Failed;

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

VerificationTest[FreeQ[docNotebooks, "XXXX" | "TODO" | $Failed | _MakeBoxes | _StringReplace], True,
	TestID -> "No template placeholders failed expressions or unresolved box generation"]

VerificationTest[And @@ (docHeldInputQ[docParseInput[#]] & /@
	Cases[docNotebooks, Cell[content_, "Input" | "ExampleInitialization", ___] :> content, Infinity]), True,
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

VerificationTest[MemberQ[docTutorialFiles, FileNameJoin[{docTutorialDirectory, "ExteriorCalculus.nb"}]], True,
	TestID -> "The exterior calculus tutorial is included"]

VerificationTest[And @@ Table[
	With[{nb = Get[file], uri = "xTerior/tutorial/" <> FileBaseName[file]},
		!FreeQ[nb, Cell[_, "Title", ___]] &&
		!FreeQ[nb, Cell[_, "Section", ___]] &&
		!FreeQ[nb, Cell[_, "Input", ___]] &&
		!FreeQ[nb, Cell["Tech Note", "Categorization", ___]] &&
		!FreeQ[nb, Cell[uri, "Categorization", ___]]],
	{file, docTutorialFiles}], True,
	TestID -> "Tutorials have a title sections examples and tech-note metadata"]

VerificationTest[!FreeQ[Get[docGuideFile],
	Cell[TextData[{ButtonBox["Exterior Calculus with xTerior", ___,
		ButtonData -> "paclet:xTerior/tutorial/ExteriorCalculus", ___]}], "GuideTutorial", ___]], True,
	TestID -> "The main guide links to the exterior calculus tutorial"]
