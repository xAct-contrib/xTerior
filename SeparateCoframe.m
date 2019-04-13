(* ::Package:: *)

(* ::Input::Initialization:: *)
SeparateCoframeAux[form_, list_IndexList:IndexList[]] := 
	Module[{auxind, formind, lengthformind, compind, deg, tensorhead, mfd, i, j, tensor, tensorind, result},
		formind = FindIndices[form];
		lengthformind = Length@formind;
		If[list =!= IndexList[],
			auxind = list,
			auxind = formind
		];

		(*If there is a lower index, pick the symbol*)
		For[i = 1, i<= Length@ auxind, i++,
			If[Length@auxind[[i]] != 0,
				auxind[[i]] = auxind[[i,2]]
			]
		];

		(*Delete Duplicates*)
		auxind = DeleteDuplicates@auxind;
		deg = Deg[form];
		mfd = First@ManifoldsOf[form];
		compind = IndexList@@ GetIndicesOfVBundle[Tangent[mfd], deg, auxind];
		(* Unique names for each tensor*)
		tensorhead = GiveSymbol[Head@form, "Tens"];

		(*Use this to avoid dummy indices in DefTensor*)
		tensorind = IndexList@@ GetIndicesOfVBundle[Tangent[mfd], lengthformind, Join[auxind, compind]];
		tensor = tensorhead@@ Join[tensorind, -#&/@ compind];
		If[Not@xTensorQ@tensorhead,
			DefTensor[tensor, mfd,
				PrintAs->ToString[
					StringForm["\!\(\*SuperscriptBox[\(\[InvisiblePrefixScriptBase]\), \(\[CircleTimes]\)]\)`1`",
						PrintAs[Evaluate@Head@ form]
					],
					StandardForm
				]
			];

		(*Inherit symmetry from form and respect the antisymmetry in the "component indices"*)
			xUpSet[SymmetryGroupOfTensor[tensorhead], 
					StrongGenSet[Range[lengthformind + deg],
						Join[SymmetryGroupOfTensor[Head@ form][[2]],
							Replace[
								Table[
									Times[-1, Cycles[List[j, j + 1]]], {j, lengthformind + 1, lengthformind + deg - 1}
								] // OutputForm, head_[arg__]:> GenSet[arg], {0,1}
							]
						] //.{}->Sequence[]
					]
				]
			];

		(*Redefine tensor for return*)
		tensor = tensorhead@@ Join[formind, -#&/@ compind];
		result = (1/Factorial[deg]) tensor Apply[Wedge, Coframe[mfd]/@ compind];
		Return@Validate@ {result, FindIndices[Evaluate@ result]}
	]

