(*Threads*)

FormsToTensors[expr_Plus, frame_] := FormsToTensors[#, frame] & /@ expr
FormsToTensors[expr_List, frame_] := FormsToTensors[#, frame] & /@ expr

(*Better idea to implement thread over Wedge? Simple map seems to fail*)
FormsToTensors[expr_Wedge, frame_] := Wedge @@ FormsToTensors[List @@ expr, frame]

FormsToTensors[expr_Times, frame_] := FormsToTensors[#, frame] & /@ expr
FormsToTensors[expr_Diff, frame_] := FormsToTensors[#, frame] & /@ expr


(*Main function*)

FormsToTensors[form_, frame : (Coframe | dx)] := Module[
  {head, tensorhead, deg, i, dummies, indices, ind, basis, out, vb, list1, cycles},
  head = Head@form;
  deg = Deg@form;

  (* Check that form has tensor Head and non-trivial degree, otherwise return form *)
  If[xTensorQ@head && deg !=  0,
   	tensorhead = GiveSymbol[head, Tensor];
   	dummies = {};
   	indices = FindIndices@form;
   	basis = 1;
   	vb = (HostsOf@head // Last);

   	(*Check if tensor exists; user may have already defined a tensor with this name. We can resolve this with a special mark as suggested*)
   	If[xTensorQ[tensorhead],
   		(*For the moment I am redefining it*)
     	ClearAll[GiveSymbol[head, Tensor] // Evaluate]
      ];

   	(*Make tensor entry*)
	xTensorQ[tensorhead] ^= True;
	DependenciesOfTensor[tensorhead] ^= DependenciesOfTensor@head; 
	PrintAs[tensorhead // Evaluate] ^=  PrintAs[head // Evaluate];
	HostsOf[tensorhead] ^= {DependenciesOfTensor[tensorhead], vb};
	SlotsOfTensor[tensorhead] ^= SlotsOfTensor[head];
	DefInfo[TestTensor] ^= {"tensor", ""};
	TensorID[TestTensor] ^= {};
	xAct`xTensor`Dagger[TestTensor] ^= TestTensor;

   	(*Initialize SG variables*)
   	list1 = {};
   	cycles = {};

   	(*Make slots and SG*)
  	For[i = 1, i <= deg, i++,
    	SlotsOfTensor[tensorhead] ^= Append[SlotsOfTensor[tensorhead], -vb];
    	list1 = Prepend[list1, Length@SlotsOfTensor@head + deg + 1 - i];
    	cycles = Prepend[cycles, 
      					If[i != deg, -Cycles[{Length@SlotsOfTensor@head + deg - i, Length@SlotsOfTensor@head + deg + 1 - i}],
      						Nothing
      	  				  ]
       					];

    	(* Initialize index stuff *)
    	ind = DummyIn@vb;
    	dummies = Append[dummies, ind];
    	basis = Wedge[basis, frame[DependenciesOfTensor@head // First][ind]];
       ];

   	(*Join SG*)
   	SymmetryGroupOfTensor[tensorhead] ^= JoinSGS[SymmetryGroupOfTensor[head], StrongGenSet[list1, GenSet @@ cycles]];
   		
   	(* Output *)
   	out = Wedge[(1/Factorial[deg]) tensorhead @@ Join[indices, ChangeIndex[dummies] /. List -> IndexList // Evaluate], basis];
   	Return@out,
   	Return@form
   ]
]