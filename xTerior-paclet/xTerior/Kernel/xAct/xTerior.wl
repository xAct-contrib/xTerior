(* ::Package:: *)

(*
	xTerior
	Exterior Calculus for xAct
	Alfonso Garc\[IAcute]a-Parrado
	agparrado@uco.es
	Universidad de C\[OAcute]rdoba, Spain

	Leo C. Stein
	leo.stein@gmail.com
	U of MS, MS, United States

	(c) 2013-2019, under GPL

	http://www.xAct.es/
	http://groups.google.com/group/xAct
	https://github.com/xAct-contrib
	xTerior is a package for exterior calculus in xAct.

	xTerior is distributed under the GNU General Public License, and runs on top of xTensor and xCoba which are free packages for fast
	manipulation of abstract and component tensor expressions. All packages can be downloaded from http://www.xact.es/
*)


(* ::Input::Initialization:: *)
xAct`xTerior`$xTensorVersionExpected = {"1.1.2", {2015, 8, 23}};
xAct`xTerior`$Version = {"0.9.1", {2019, 5, 17}};

(******************************************************************************)
(********************* 1. Initialization **************************************)
(******************************************************************************)

(************************ 1.1 GPL *********************************************)
(* xTerior: exterior calculus in Differential Geometry *)

(* Copyright (C) 2013-2019 Alfonso Garcia-Parrado Gomez-Lobo and Leo C. Stein *)

(* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License,or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place-Suite 330, Boston, MA 02111-1307, USA. 
*)


(*********************** 1.2 Info Package ************************************)

(* :Title: xTerior *)

(* :Author: Alfonso Garcia-Parrado Gomez-Lobo and Leo C. Stein *)

(* :Summary: exterior calculus in Differential Geometry *)

(* :Brief Discussion:
   - xTerior extends xAct to work with differentiable forms in general manifolds.
   - Introduces the exterior algebra, the exterior derivative, the Hodge dual, the connection and curvature forms for\.07\.14\.07 an arbitrary connection, the exterior covariant derivative.
   
*)
  
(* :Context: xAct`xTerior` *)

(* :Package Version: 0.9.1 *)

(* :Copyright: Alfonso Garcia-Parrado Gomez-Lobo and Leo C. Stein (2013-2019) *)

(* :History: See xTerior.History *)

(* :Keywords: *)

(* :Source: xTerior.nb *)

(* :Warning: *)

(* :Mathematica Version: 9.0 and later *)

(* :Limitations:
	- ?? *)


(*********************** 1.3 BeginPackage *********************************)

With[{
		xAct`xTerior`Private`xTeriorSymbols = 
		DeleteCases[
				Join[
					Names[
						"xAct`xTerior`*"
					], 
					Names[
						"xAct`xTerior`Private`*"
					]
				], "$Version" | "xAct`xTerior`$Version" | "$xTensorVersionExpected" | "xAct`xTerior`$xTensorVersionExpected"
		]
	},
	Unprotect /@ xAct`xTerior`Private`xTeriorSymbols;
	Clear /@ xAct`xTerior`Private`xTeriorSymbols;
]

If[Unevaluated[xAct`xCore`Private`$LastPackage] === xAct`xCore`Private`$LastPackage, xAct`xCore`Private`$LastPackage = "xAct`xTerior`"];

(* Temporary fix for the conflicts with Wolfram 14.1 *)

xAct`xTerior`Diff;

(* Explicit (not hidden) import of xTensor, xPerm and xCore: Alfonso: do we need xCoba ? *)
BeginPackage["xAct`xTerior`", {"xAct`xCoba`", "xAct`xTensor`", "xAct`xPerm`", "xAct`xCore`"}]



If[Not @ OrderedQ @ Map[
		Last, {
			xAct`xTerior`$xTensorVersionExpected, xAct`xTensor`$Version
		}
	],
	Throw @ Message[
		General::versions, "xTensor", xAct`xTensor`$Version, xAct`xTerior`$xTensorVersionExpected
	]
]



Print[xAct`xCore`Private`bars]
Print["Package xAct`xTerior`  version ", xAct`xTerior`$Version[[1]], ", ", xAct`xTerior`$Version[[2]]];
Print["Copyright (C) 2013-2019, Alfonso Garcia-Parrado Gomez-Lobo and Leo C. Stein, under the General Public License."];



Off[General::shdw]
xAct`xTerior`Disclaimer[] := Print["These are points 11 and 12 of the General Public License:\n\n
BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM `AS IS\.b4 
WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. 
SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n\n
IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO 
MAY MODIFY AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, 
INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF 
DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM 
TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES."]
On[General::shdw]


(* If xTerior is not being called from other package then write this GPL short disclaimer: *)
If[xAct`xCore`Private`$LastPackage === "xAct`xTerior`",
	Unset[
		xAct`xCore`Private`$LastPackage
	];
	
	Print[
		xAct`xCore`Private`bars
	];
	
	Print[
		"These packages come with ABSOLUTELY NO WARRANTY; for details type Disclaimer[]. This is free software, and you are welcome to redistribute it under certain conditions. See the General Public License for details."
	];
	
	Print[
		xAct`xCore`Private`bars
	]
]


(************************* 1.4. Non-standard setup ***********************************)

(* Screen all dollar indices: *)
$PrePrint = ScreenDollarIndices;


(* Switch off messages issued by ManifoldOfCovD acting on PD. *)
Unprotect @ PD;
	PD /: ManifoldOfCovD @ PD = .
Protect @ PD;


(***************************** 1.5. Usage messages  ********************************************)

(* Definition and undefinition of a differential form (just a wrapper for DefTensor with the option GradeOfTensor\[Rule]{Wedge})*)
DefDiffForm::usage = "DefDiffForm[form[inds], mani, Deg] defines a tensor valued differential form of degree deg on the manifold mani";
UndefDiffForm::usage = "UndefDiffForm[form] undefines the differential form form";

(* Grade of a differential form *)
Deg::usage = "Deg[form] returns the grade of a differential form";

DiffFormQ::usage = "DiffFormQ is an option for LieToCovD which specifies whether the expression which is acted 
upon should be regarded as a differential form or not. This is an option added by xTerior and it is not present any other package using LieToCovD.";

(* Exterior derivative and exterior covariant derivative *)
Diff::usage="Diff[form] computes the exterior derivative of form. Diff[form,covd] computes the exterior 
covariant derivative of form with respect to the covariant derivative covd.";

FindPotential::usage = "FindPotential[form, point, chart] uses the Poincar\[EAcute] 
lemma to compute a potential for a closed form form (no checks are carried out to ensure that the form is actually closed). 
The form must be written in some explicit coordinates belonging to chart 
and the argument point is the point, assumed to be in the same coordinate chart as form, which defines a star-shaped region 
where the potential is differentiable. A change in the point will give 
in general a different potential";

(* Computation of the exterior covariant derivative *)
ChangeExtD::usage = "ChangeExtD[expr, cd1, cd2] expresses the exterior covariant derivative taken 
with respect to the connection defined by the covariant derivative cd1 in terms of the exterior 
covariant derivative taken with respect to the connection defined by the covariant derivative cd2";

(* Computation of the generalized (index-free) covariant derivative *)
ChangeGenCovD::usage = "ChangeGenCovD[expr,cd1,cd2] expresses the generalized covariant derivative taken with respect 
to the connection defined by the covariant derivative cd1 in terms of the 
generalized covariant derivative taken with respect to the connection defined by the covariant derivative cd2";

(* Hodge dual *)
Hodge::usage = "Hodge[metric][expr] is the Hodge dual of expr defined with respect to metric";
ExpandHodgeDual::usage="ExpandHodgeDual[expr,Coframe[mani],g] expands out all the Hodge duals of the exterior powers of Coframe[mani], 
defined with respect to the metric g. If the manifold tag 
mani is dropped, then all the instances of Coframe are expanded. The Coframe label can be replaced by dx if we are using the holonomic coframe.";

(* Co-differential *)
Codiff::usage = "Codiff[metric][form] is the co-differential of form computed with respect to metric";

(* Expansion of the co-differential *)
CodiffToDiff::usage = "CodiffToDiff[expr] replaces all the instances of the co-differential in expr by their expansion in terms of the exterior derivative.";

(* Interior contraction *)
Int::usage = "Int[v][form] is the interior contraction of form with the vector (rank 1-tensor) v";

(* Lie derivative on forms *)
CartanD::usage = "CartanD[v][form] is the Cartan derivative of form with respect to the vector (rank 1-tensor) v. CartanD[v][form,covd] is the Cartan derivative with respect to the covariant derivative covd.";

(* Cartan formula for Lie derivatives *)
CartanDToDiff::usage = "CartanDToDiff[expr] replaces the Cartan derivative of all the differential forms in expr by their expansion obtained by means of the Cartan formula";

(* Put derivations into canonical order *)
SortDerivations::usage = "SortDerivations[expr] brings expr into a new expression where all the derivations 
(exterior derivative, Lie derivative and interior contraction) are written in a canonical order. The default left-to-right order is defined by the variable $DerivationSortOrder";
$DerivationSortOrder::usage="$DerivationSortOrder is a global variable which encodes the default ordering of the three derivatives Int, LieD and Diff. The default is {LieD, Int, Diff}";

(* Variational derivative on forms *)
FormVarD::usage = "FormVarD[form,metric][expr] computes the variational derivative of expr, which must be a n-form with n the manifold dimension, with respect to form. 
In the computation, exterior derivatives are transformed into co-differentials taken with respect of metric.";

(* Canonical forms on the frame bundle *)
Coframe::usage = "Coframe[mani] is the set of canonical 1-forms defined in the frame bundle arising from the manifold mani";
dx::usage = "dx[mani] represents a holonomic co-frame in the manifold mani.";

(* Holonomic & anholonomic frames *)
PDFrame::usage = "PDFrame[mani] represents a holonomic frame in the manifold mani.";
eFrame::usage = "eFrame[mani] represents a frame in the manifold mani.";

(* Koszul connection *)
Koszul::usage = "Koszul is the head used to introduce the Koszul operator. The latter is constructed by prepending this head to the xTensor symbol used to represent a covariant derivative";

ExpandKoszul::usage="ExpandKoszul[expr,cd1,cd2] expands Koszul derivations of frame or co-frame elements with respect to cd1 into Christoffel tensors relating the covariant 
derivatives cd1 and cd2";

ExpandCovD::usage="ExpandCovD[expr,cd,chart] expands all the generalized covariant derivatives of rank 1 tensors in terms of the Koszul derivatives with respect to elements 
of a holonomic frame associated to a coordinate chart";

(* The connection 1-form *)
ConnectionForm::usage = "ConnectionForm[cd,vbundle] represents the connection 1-form associated 
to the covariant derivatives cd defined in the bundle vbundle. If vbundle is the tangent bundle of a differentiable manifold then ConnectionForm is automatically 
replaced by ChristoffelForm (see 
the on-line help of ChristoffelForm for further details).";

(* The curvature 2-form *)
CurvatureForm::usage = "CurvatureForm[cd,vbundle] represents the curvature 2-form associated to the covariant derivative cd. If vbundle is the tangent bundle of a differentiable manifold then 
CurvatureForm is replaced by RiemannForm";

(* Connection 2-form for a connection in a frame bundle*)
ChristoffelForm::usage="ChristoffelForm[cd] is the connection 1-form associated to the covariant derivative cd which is a covariant derivative in the tangent bundle of a manifold";

(* Curvature 2-form for a connection in a frame bundle *)
RiemannForm::usage = "RiemannForm[cd] is the curvature 2-form associated to the covariant derivative cd which is a covariant derivative in the frame bundle of a manifold";
(* Transformation of the connection form to the connection tensor *)

ConnectionFormToTensor::usage = "ConnectionFormToTensor[expr,covd,frame] transforms all instances of the connection form into the (A)Christoffel tensor which relates the covariant derivative 
defining the connection form and covd. The variable frame can take the value of either Coframe or dx. If the (A)Christoffel tensor does not exist it is created automatically.";
CurvatureFormToTensor::usage="CurvatureFormToTensor[expr,frame] transforms all the instances of the curvature form into the related Riemann or FRiemann tensor, inserting the corresponding frame 
(either Coframe or dx).";

ChangeCurvatureForm::usage = "ChangeCurvatureForm[curvature, cd1, cd2] writes the curvature 2-form curvature[cd1, vbundle] in terms of the curvature 2-form curvature[cd2, vbundle]";

(* The torsion 2-form *)
TorsionForm::usage = "TorsionForm[cd] represents the torsion 2-form arising from the covariant derivative cd (cd must be defined on the tangent bundle of a manifold)";

(* Cartan structure equations *)
UseCartan::usage = "UseCartan[expr, covd] expands all the instances of the Diff using the Cartan structure equations for the connection arising from covd. 
In this way it is possible to expand the exterior derivative of a co-frame, a torsion 2-form and the curvature 2-form. If covd is the Levi-Civita connection of a metric, then the exterior 
derivatives of that metric and its volume element are expanded too. UseCartan[expr,PD] expands all instances of the exterior derivative in terms of partial derivatives defined in the list of 
manifolds returned by ManifoldsOf[expr]. It is possible to specify a custom list of manifolds as a third argument in the form UseCartan[expr,PD,{M1,M2,..}]";

(* Zero forms *)
ZeroDegQ::usage = "ZeroDegQ[expr] returns True if the degree of expr is zero";
UseDimensionStart::usage = "UseDimensionStart[] is an instruction that, when issued, 
makes any expression with degree greater than the manifold dimension equal to zero.";

UseDimensionStop::usage = "UseDimensionStop[] cancels the action of UseDimensionStart[]";

(* Integration of differential forms over manifolds with boundary *)
FormIntegrate::usage = "FormIntegrate[form, M] represents the integration of the given differential form, of degree n, over the n-dimensional manifold M with boundary.";

UseStokes::usage="UseStokes[expr, M] converts subexpressions FormIntegrate[Diff[form], M] in expr into FormIntegrate[form, ManifoldBoundary[M]], 
reducing n-dimensional integration to (n-1)-
integration. UseStokes[expr, form] converts subexpressions FormIntegrate[form, ManifoldBoundary[M]] in expr into FormIntegrate[Diff[form], M]. UseStokes[expr] performs UseStokes[expr, M] wherever 
possible in expr, to reduce dimensionality of integration.";


(************************************** 1.6 BeginPrivate *********************************)
Begin["`Private`"]


(******************************************************************************)
(********************* 2. Basic structures ************************************)
(******************************************************************************)


(***************** 2.1 Definition of the wedge product ************************)
(*
The wedge product is an associative, anticommutative (actually supercommutative) graded product with identity 1. 
The scalars are those objects of grade 0, including the identity 1, so that the scalars are actually 
in this case true elements of the algebra. The product by scalar and the product of scalars are both Times, 
 so we do not need to specify them.
*)

DefProduct[Wedge,
	AssociativeProductQ -> True,
	CommutativityOfProduct -> "SuperCommutative",
	GradedProductQ -> True,
	IdentityElementOfProduct -> 1,
	ScalarsOfProduct -> (SameQ[Grade[#, Wedge], 0]&),
	DefInfo -> Null
];


(* Relation between Wedge and Times. *)
Wedge /: GradeOfProduct[Times, Wedge] = 0;


(* When working with the exterior algebra the grade is typically called degree: *)
Deg[expr_] := Grade[expr, Wedge];


(* Behavior of the wedge product with respect to Dagger *)
Unprotect @ Dagger;
Dagger @ expr_Wedge := Dagger /@ expr;
Protect @ Dagger;


(* Vanishing of forms whose degree is higher than the manifold(s) dimension. *)
(*Code added by Jos\[EAcute]*)
$UseDimensionsQ=False;

$DimensionsZeroForms = {};

SetZeroForm[form_] := 
	If[GradeOfTensor[form, Wedge] > Plus@@ (DimOfManifold/@DependenciesOfTensor[form]),
		form[___] = 0; 
		AppendTo[$DimensionsZeroForms, form]
	];
	
UnsetZeroForm[form_]:=Unset[form[___]];


UseDimensionStart[] :=
	Module[{},
		If[$UseDimensionsQ,
			Return[]
		];
		
		$UseDimensionsQ = True;
		
		(*Forms whose degree is greater than the dimension*)
		
		SetZeroForm /@ $Tensors;
		
		(*Expressions which are wedge products*)
		
		HoldPattern @ Wedge[expr__] :=
			0 /; (Plus @@ (Grade[#, Wedge]& /@ {expr}) > (Plus @@ (DimOfManifold /@ Union @ Flatten[DependenciesOf /@ {expr}])));
		
		(*Expressions which are exterior derivatives*)
		
		HoldPattern @ Diff[expr_, PD] :=
			0 /; (1 + Plus @@ (Grade[#, Wedge]& /@ {expr}) > (Plus @@ (DimOfManifold /@ DependenciesOf[expr])));
		
		HoldPattern @ Diff[expr_, covd_] :=
			0 /; (1 + Plus @@ (Grade[#, Wedge]& /@ {expr}) > (Plus @@ (DimOfManifold /@ DependenciesOf[expr])));
	]


UseDimensionStop[] :=
	Module[{},
		If[!$UseDimensionsQ,
			Return[]
		];
		$UseDimensionsQ = False;
		UnsetZeroForm /@ Union @ $DimensionsZeroForms;
		HoldPattern @ Wedge[expr__] = .;
		HoldPattern @ Diff[expr_, PD] = .;
		HoldPattern @ Diff[expr_, covd_] = .;
	]

(* xTensions *)

xTension["xTerior`", DefTensor, "End"] :=
	DefTensorUseDimensions;
DefTensorUseDimensions[tensor_[inds___], __] :=
	If[$UseDimensionsQ,
		SetZeroForm[tensor]
	];

xTension["xTerior`", UndefTensor, "Beginning"] :=
	UndefTensorUseDimensions;
UndefTensorUseDimensions[tensor_] :=
	If[$UseDimensionsQ,
		UnsetZeroForm[tensor]
	];


(********************* 2.2 Integration of CTensor in Wedge **************************************)
(*

In this section we implement the Wedge product of CTensor objects. This code has been supplied by Jos\[EAcute] Mart\[IAcute]n-Garc\[IAcute]a.

*)

(* Grade of a CTensor with respect to Wedge *)
CTensor /: Grade[expr_CTensor, Wedge] :=
	Module[{list, grades, comps = First @ expr},
		list = Flatten[{comps}];
		grades = Grade[#, Wedge]& /@ list;
		If[Length @ Tally[grades] == 1,
			grades[[1]],
			$Failed
		]
	];


(* Contracted wedge product of CTensor objects *)

Wedge[ctensor1_CTensor[left1___, a_, right1___], ctensor2_CTensor[left2___, -a_, right2___]] :=
	Module[{n1 = Length[{left1, a}], n2 = Length[{left2, -a}], res},
		res = xAct`xCoba`Private`CTensorContract[ctensor1, ctensor2, {n1, n2}, Wedge];
		res[left1, right1, left2, right2] /; FreeQ[res, $Failed]
	]; 

Wedge[ctensor1_CTensor[left1___, -a_, right1___], ctensor2_CTensor[left2___, a_, right2___]] :=
	Module[{n1 = Length[{left1, a}], n2 = Length[{left2, -a}], res},
		res = xAct`xCoba`Private`CTensorContract[ctensor1, ctensor2, {n1, n2}, Wedge];
		res[left1, right1, left2, right2] /; FreeQ[res, $Failed]
	]; 


(* ::Input::Initialization:: *)
(*signatureOrZero[indices_]:=If[DuplicateFreeQ[indices],Signature[Ordering[indices]],0];*)

simplifyBasisWedge[expr_] :=
	expr /. wed_Wedge :> simplifyBasisWedge1[wed]; 

(*simplifyBasisWedge1[HoldPattern[Wedge[factors:((head:((xAct`xTerior`Coframe|xAct`xTerior`dx)[_]))[_]..)]]]:=With[{indices=First/@{factors}},signatureOrZero[indices] Wedge@@(head/@Sort[indices])];*)

simplifyBasisWedge1[HoldPattern[expr:Wedge[factors:Blank[]..]]] :=
	Module[{indices = FindFreeIndices[expr]},
		If[indices === IndexList[],
			expr,
			$Failed
		]
	]; 

(* Wedge product of general CTensor objects *)
CTensorWedge[] :=1; 

CTensorWedge[ctensors__CTensor] :=
	CTensor[simplifyBasisWedge[xAct`xCoba`Private`tensorproduct[Wedge] @@ #1], Join @@ #2, Plus @@ #3]& @@ Transpose[List @@@ {ctensors}]; 

CTensorWedge[___, Zero, ___] :=Zero; 

Wedge[ctensor1_CTensor[inds1___], ctensor2_CTensor[inds2___]] :=
	CTensorWedge[ctensor1, ctensor2][inds1, inds2] /; xAct`xTensor`Private`TakePairs[
		{inds1, inds2}] === {}

Wedge[ctensor1_CTensor, ctensor2_CTensor] := CTensorWedge[ctensor1, ctensor2]


(* Wedge product of frame and co-frame with CTensor objects. We first have a definition for Coframe with c - indices, avoiding the use of ToCTensor : *)
Wedge[basis: (Coframe | dx)[_][_?CIndexQ], CTensor[array_, bases_List, addweight_][b__]] :=
	CTensor[Wedge[basis, array], bases, addweight][b]

Wedge[CTensor[array_, bases_List, addweight_][b__], basis:(Coframe | dx)[_][_?CIndexQ]] :=
	CTensor[Wedge[array, basis], bases, addweight][b]


(* Then we have another definition for the general case of Wedge, in which the Coframes inside the CTensor uniquely identify the frame to use in ToCTensor : *)
FormBases[CTensor[array_, __]] :=
	DeleteDuplicates[xAct`xPerm`Private`nosign /@ Cases[array, (Coframe | dx)[_][{_, frame_}] :> frame, {0, Infinity}]]

Wedge[(basis:(Coframe | dx)[_])[ind_], ctensor_CTensor[inds___]] :=
	With[{frames = FormBases[ctensor]},
		Wedge[ToCTensor[basis, frames][ind], ctensor[inds]] /; Length[frames] === 1
	]

Wedge[ctensor_CTensor[inds___], (basis:(Coframe | dx)[_])[ind_]] :=
	With[{frames = FormBases[ctensor]},
		Wedge[ctensor[inds], ToCTensor[basis, frames][ind]] /; Length[frames] === 1
	]


(* Finally we need another definition for the case in which the CTensor is a 0 - form : *)
CTensor /: (basis:(Coframe | dx)[_])[a_] (ctensor:CTensor[_, bases_List, _][l___, -a_, ___]) :=
	ToCTensor[basis, {-bases[[Length[{l, -a}]]]}][a] ctensor

CTensor /: HoldPattern[Wedge[lw___, (basis:(Coframe | dx)[_])[a_], rw___] (ctensor:CTensor[_, bases_List, _][l___, -a_, ___])] :=
	Wedge[lw, ToCTensor[basis, {-bases[[Length[{l, -a}]]]}][a], rw] ctensor


(* Wedge product of index free quantities with CTensor objects: *)
Wedge[ten_, CTensor[array_, bases_List, addweight_][b__]] :=
	CTensor[Wedge[ten, array], bases, addweight][b] /; FindFreeIndices @ ten === IndexList[]

Wedge[CTensor[array_, bases_List, addweight_][b__], ten_] :=
	CTensor[Wedge[array, ten], bases, addweight][b] /; FindFreeIndices @ ten === IndexList[]



(************************** 2.3 Definition of the CircleTimes product ****************************)

(* By default we define the CircleTimes "\[CircleTimes]" product. *)
$DefInfoQ = False;

DefProduct[
	CircleTimes,
	AssociativeProductQ -> True,
	IdentityElementOfProduct -> 1,
	GradedProductQ -> True,
	ScalarsOfProduct -> (
			SameQ[Grade[#,CircleTimes],0]&
		),
	PrintAs -> "\[CircleTimes]"
]

$DefInfoQ = True;


(* Grade of CircleTimes with respect to Times *)
CircleTimes /: GradeOfProduct[Times, CircleTimes] = 0;


(* Exterior algebra grade and tensor grade should coincide on elementary objects: *)
CircleTimes /: GradeOfTensor[head_, CircleTimes] := GradeOfTensor[head, Wedge];


(* We also need to specify the Grade of the inert - head expressions : *)
CircleTimes /: Grade[Diff[expr_, _], CircleTimes] := Grade[expr, CircleTimes] + 1;

CircleTimes /: Grade[Hodge[metricg_][expr_], CircleTimes] := DimOfVBundle[VBundleOfMetric[metricg]] - Grade[expr, CircleTimes];


(* Behavior with respect to Dagger *)
Unprotect @ Dagger;
	Dagger @ expr_CircleTimes := Dagger /@ expr;
Protect @ Dagger;

(*
The CircleTimes - grade of a form coincides with the Wedge - grade of that form. 
However, the Wedge - grade is not well defined for generic expressions having positive CircleTimes - grade unless they are forms. 
That is why we have given CircleTimes - grade in terms of Wedge - grade and not viceversa.There is no general way in which we can know 
the relations between the grades of a given expression in several products.The user must define carefully the relations. 
Usually it is safer to declare independently the grades in each product. 
Recall that objects for which a grade has not been given are taken to have degree 0.
*)



(********************** 2.4 Integration of CTensor in CircleTimes **************************************)


(* Grade of a CTensor with respect to CircleTimes *)
CTensor /: Grade[expr_CTensor, CircleTimes] :=
	Module[{list, grades, comps = First @ expr},
		list = Flatten[{comps}];
		grades = Grade[#, CircleTimes]& /@ list;
		If[Length @ Tally[grades] == 1,
			grades[[1]],
			$Failed
		]
	]; 


(* Algebra with unit element *)
xAct`xCoba`Private`tensorproduct[CircleTimes][] = 1;


(* Contracted CircleTimes product of CTensor objects *)

CircleTimes[ctensor1_CTensor[left1___, a_, right1___], ctensor2_CTensor[left2___, -a_, right2___]] :=
	Module[{n1 = Length[{left1, a}], n2 = Length[{left2, -a}], res},
		res = xAct`xCoba`Private`CTensorContract[ctensor1, ctensor2, {n1, n2}, CircleTimes];
		res[left1, right1, left2, right2] /; FreeQ[res, $Failed]
	]; 

CircleTimes[ctensor1_CTensor[left1___, -a_, right1___], ctensor2_CTensor[left2___, a_, right2___]] :=
	Module[{n1 = Length[{left1, a}], n2 = Length[{left2, -a}], res},
		res = xAct`xCoba`Private`CTensorContract[ctensor1, ctensor2, {n1, n2}, CircleTimes];
		res[left1, right1, left2, right2] /; FreeQ[res, $Failed]
	]; 


(* CircleTimes product of General CTensor objects  *)
simplifyBasisCircleTimes[expr_] :=
	expr /. wed_CircleTimes :> simplifyBasisCircleTimes1[wed]; 

simplifyBasisCircleTimes1[HoldPattern[expr:CircleTimes[factors:Blank[]..]]] :=
	Module[{indices = FindFreeIndices[expr]},
		If[indices === IndexList[],
			expr,
			$Failed
		]
	]; 

(* CircelTimes product of general CTensor objects *)
CTensorCircleTimes[] :=1; 

CTensorCircleTimes[ctensors__CTensor] :=
	CTensor[simplifyBasisCircleTimes[xAct`xCoba`Private`tensorproduct[CircleTimes] @@ #1], Join @@ #2, Plus @@ #3]& @@ Transpose[List @@@ {ctensors}]; 

CTensorCircleTimes[___, Zero, ___] :=Zero; 

CircleTimes[ctensor1_CTensor[inds1___], ctensor2_CTensor[inds2___]] :=
	CTensorCircleTimes[ctensor1, ctensor2][inds1, inds2] /; xAct`xTensor`Private`TakePairs[{inds1, inds2}] === {}

CircleTimes[ctensor1_CTensor, ctensor2_CTensor] :=
	CTensorCircleTimes[ctensor1, ctensor2]


(* CircleTimes product of index free quantities with CTensor objects: *)
CircleTimes[ten_, CTensor[array_, bases_List, addweight_][b__]] := CTensor[CircleTimes[ten, array], bases, addweight][b] /; FindFreeIndices @ ten === IndexList[]

CircleTimes[CTensor[array_, bases_List, addweight_][b__], ten_] := CTensor[CircleTimes[array, ten], bases, addweight][b] /; FindFreeIndices @ ten === IndexList[]


(**************** 2.5 Definition and undefinition of differential forms *****************************)

(* Definition and undefinition of differential forms. This is simply DefTensor with the appropriate options *)

DefDiffForm[form_, mani_, deg_, options___?OptionQ] :=
	(
		DefTensor[form, mani, GradeOfTensor -> {Wedge -> deg}, options];
		
	)

DefDiffForm[form_, mani_, deg_, sym_, options___?OptionQ] :=
	(
		DefTensor[form, mani, sym, GradeOfTensor -> {Wedge -> deg}, options
			];
		
	)

Options @ DefDiffForm := Options @ DefTensor;

UndefDiffForm := UndefTensor;

Protect[DefDiffForm, UndefDiffForm];


(****************** 2.6 Graded derivations *************************)


(* This function will be used to declare the three derivations we need, namely diff, Int[v] and lie[v].  *)

Options[DefGradedDerivation] = {PrintAs -> Identity, Master -> Null}; 

GradeOfDerivation[head_[v_, rest___], prod_] := GradeOfDerivation[head, prod] + Grade[v, prod];

DefGradedDerivation[der_, prod_?ProductQ, dergrade_:0, options:OptionsPattern[]] :=
	With[{head = SubHead[der]},
		Module[{pa},
			{pa} = OptionValue[{PrintAs}]; 
			(* DefInertHead will take care of scalar-homogeneity and linearity *)
			DefInertHead[der, LinearQ -> True, ContractThrough -> {delta}, PrintAs-> pa, DefInfo -> Null];
			(* Other properties of a derivation *)
			MakeDerivation[head, der, NoPattern[der], prod, dergrade];(* Nonatomic derivation *)
			If[der =!= head,
				(* additivity in the vector slot (but not homogeneity!) *)
				head[0][__] := 0;
				head[v_Plus][args__] := head[#][args]& /@ v; 
				(* Subscript vector argument for formatting *)
				If[pa === Identity,
					pa = PrintAs[head]
				];
				head /: MakeBoxes[head[v_][form_], StandardForm] := 
					xAct`xTensor`Private`interpretbox[
						head[v][form], 
						RowBox[{SubscriptBox[pa, MakeBoxes[v, StandardForm]], "[", MakeBoxes[form, StandardForm], "]"}]
					];
				
			]
		]
	]; 


(* This part is separated in order to avoid renaming confusion between derL and derR: *)
MakeDerivation[head_, derL_, derR_, prod_, dergrade_] :=
	With[
		{grade = GradeOfDerivation[derR, prod]},
		(* Addition of grades in algebra *)
		head /: GradeOfDerivation[head, prod] := dergrade;
		head /: Grade[derL[expr_, ___], prod] := Grade[expr, prod] + grade;
			(* The (graded) Leibniz rule *)
		derL[expr_prod, rest___] :=
			With[{sumgrades = FoldList[Plus, 0, Grade[#, Wedge]& /@ List @@ expr]},
				Sum[(-1)^(grade * sumgrades[[i]]) MapAt[derR[#, rest]&, expr, i],
					 {i, 1, Length[expr]}]
			];
	(* QUESTION: Agreement with a regular derivative when acting on scalar functions?? *)
		derL[func_?ScalarFunctionQ[args__], rest___] := 
			xAct`xTensor`Private`multiD[derR[#, rest]&, func[args]]; 
		(* Dependencies *)
		If[!AtomQ[derR],
			head /: DependenciesOfInertHead[derL] := DependenciesOf[First[derR]]
		]
	]; 


(*************** 2.7 Exterior differentiation ********************)

(* The second key ingredient for exterior algebra is the differential operator. 
This a concept defined per manifold, or equivalently per tangent-bundle, 
though in this notebook we create only one differential operator. *)

DefGradedDerivation[Diff, Wedge, +1, PrintAs -> "d"];


(* We always have PD as the covariant derivative. *)
Diff[expr_] := Diff[expr, PD]


(* Superscript formatting for covariant exterior derivatives *)
Diff /: MakeBoxes[Diff[form_, PD?CovDQ], StandardForm] :=
	xAct`xTensor`Private`interpretbox[Diff[form, PD], RowBox[{PrintAs[Diff], "[", MakeBoxes[form, StandardForm], "]"}]]; 

Diff /: MakeBoxes[Diff[form_, cd_?CovDQ], StandardForm] :=
	xAct`xTensor`Private`interpretbox[
		Diff[form, cd], RowBox[{SuperscriptBox[PrintAs[Diff], Last @ SymbolOfCovD[cd]], "[", MakeBoxes[form, StandardForm], "]"}]
	]; 


(* Thread over lists and equal *)
Diff[expr_?ArrayQ, cd_] := Diff[#, cd]& /@ expr;
Diff[expr_Equal, cd_] := Diff[#, cd]& /@ expr;


(* The exterior derivative applied twice is zero: *)
Diff[expr_Diff, PD] := 0


(* Exterior derivative of basis objects is zero. TODO: This is not correct. *)
Diff[_Basis, PD] := 0;


(* We still need definition when acting on Times *)
(* This produces expanded expressions and is much faster when there are many scalars *)

notZero = 0 =!= #&

Diff[expr_Times, cd_] :=
	Module[{grades = Grade[#, Wedge]& /@ List @@ expr, pos, scalar, form
		},
		pos = Position[grades, _?notZero, 1, Heads -> False];
		Which[
			Length[pos] > 1,
				Throw[
					Message[Diff::error1, "Found Times product of nonscalar forms: ", expr]
				],
			Length[pos] === 1,
				pos = pos[[1, 1]];
				scalar = Delete[expr, {pos}];
				form = expr[[pos]];
				scalar Diff[form, cd] + diff0[scalar, cd, form],
			Length[pos] === 0,
				diff0[expr, cd]
		]
	]; 

(* Only scalars *)

diff0[scalar_Times, cd_] :=
	Sum[MapAt[Diff[#, cd]&, scalar, i], {i, 1, Length[scalar]}]; 

diff0[scalar_, cd_] :=
	Diff[scalar, cd]; 

(* Scalars and a form *)

diff0[scalar_Times, cd_, form_] :=
	Sum[MapAt[diff0[#, cd, form]&, scalar, i], {i, 1, Length[scalar]}]; 

diff0[scalar_, cd_, form_] :=
	Wedge[Diff[scalar, cd], form]; 


(* Exterior derivative of constants vanish *)
Diff[x_?ConstantQ, cd_] := 0;


(* Behavior with respect to Dagger. *)
Diff /: HoldPattern[Dagger[Diff[expr_, cd_]]] := Diff[Dagger[expr], cd];


(* Behaviour with respect to CTensor *)
Diff[CTensor[array_, bases_List, weight_][inds__], cd_] := CTensor[Diff[array, cd], bases, weight][inds];


(************ Introduction of co-frame, eFrame and PDFrame *************************)

(*

We create the non-atomic tensor \[Theta][M] representing a co-frame. 
Note that the formatting does not contain the manifold as this information is already 
visible in the abstract index. The abstract index may belong to the tangent space of the manifold "M" or 
to any other vector bundle with base M. In the case of the abstract index belonging to TangentM we can think of 
\[Theta][M] as the set of "canonical 1-forms". 

*)


xTensorQ @ Coframe[mani_?ManifoldQ] ^= True; 

SlotsOfTensor[Coframe[mani_?ManifoldQ]] ^:={Tangent @ mani}; 

Coframe /: GradeOfTensor[Coframe[mani_?ManifoldQ], Wedge] = 1; 

SymmetryGroupOfTensor[Coframe[mani_?ManifoldQ]] ^= StrongGenSet[{}, GenSet[]]; 

DefInfo[Coframe[mani_?ManifoldQ]] ^= {"General co-frame", ""}; 

DependenciesOfTensor[Coframe[mani_?ManifoldQ]] ^:={mani}; 

HostsOf[Coframe[mani_?ManifoldQ]] ^= {}; 

TensorID[Coframe[mani_?ManifoldQ]] ^= {}; 

PrintAs[Coframe[mani_?ManifoldQ]] ^= "\[Theta]"; 


(* Holonomic Co-frame *)
xTensorQ @ dx[mani_?ManifoldQ] ^= True; 

SlotsOfTensor[dx[mani_?ManifoldQ]] ^:={Tangent @ mani}; 

dx /: GradeOfTensor[dx[mani_?ManifoldQ], Wedge] = 1; 

SymmetryGroupOfTensor[dx[mani_?ManifoldQ]] ^= StrongGenSet[{}, GenSet[]]; 

DefInfo[dx[mani_?ManifoldQ]] ^= {"General co-frame", ""}; 

DependenciesOfTensor[dx[mani_?ManifoldQ]] ^:={mani}; 

HostsOf[dx[mani_?ManifoldQ]] ^= {}; 

TensorID[dx[mani_?ManifoldQ]] ^= {}; 

PrintAs[dx[mani_?ManifoldQ]] ^= "dx"; 


(* Condition of the co-frame being holonomic. *)
Diff[dx[mani_?ManifoldQ][ind_], PD] := 0;


(* General frame: *)
xTensorQ @ eFrame[mani_?ManifoldQ] ^= True;

SlotsOfTensor[eFrame[mani_?ManifoldQ]] ^:= {-Tangent @ mani};

eFrame /: GradeOfTensor[eFrame[mani_?ManifoldQ], CircleTimes] = 1;

SymmetryGroupOfTensor[eFrame[mani_?ManifoldQ]] ^= StrongGenSet[{}, GenSet[]];

DefInfo[eFrame[mani_?ManifoldQ]] ^= {"General frame", ""};

DependenciesOfTensor[eFrame[mani_?ManifoldQ]] ^:= {mani};

HostsOf[eFrame[mani_?ManifoldQ]] ^= {};

TensorID[eFrame[mani_?ManifoldQ]] ^= {};

PrintAs[eFrame[mani_?ManifoldQ]] ^= "e";


(* General holonomic frame: *)
xTensorQ @ PDFrame[mani_?ManifoldQ] ^= True;

SlotsOfTensor[PDFrame[mani_?ManifoldQ]] ^:= {-Tangent @ mani};

PDFrame /: GradeOfTensor[PDFrame[mani_?ManifoldQ], CircleTimes] = 1;

SymmetryGroupOfTensor[PDFrame[mani_?ManifoldQ]] ^= StrongGenSet[{}, GenSet[]];

DefInfo[PDFrame[mani_?ManifoldQ]] ^= {"General frame", ""};

DependenciesOfTensor[PDFrame[mani_?ManifoldQ]] ^:= {mani};

HostsOf[PDFrame[mani_?ManifoldQ]] ^= {};

TensorID[PDFrame[mani_?ManifoldQ]] ^= {};

PrintAs[PDFrame[mani_?ManifoldQ]] ^= "\!\(\*FractionBox[\(\[PartialD]\),\(\[PartialD]x\)]\)";


(* xTensions *)
xTension["xTerior`",DefChart,"End"]:=setdiffs;

setdiffs[chartname_,__]:=Thread[ComponentValue[dx[ManifoldOfChart@chartname][{#,chartname}]&/@CNumbersOf@chartname,Diff/@ScalarsOfChart@chartname]];

xTension["xTerior`", DefChart, "End"] :=
	defFreeAlgebraChart;

defFreeAlgebraChart[basis_, __] :=
	With[{scalars = ScalarsOfChart[basis], numbers = CNumbersOf @ basis, mani = ManifoldOfChart @ basis},
		(* Define the algebra elements forming the basis of the free algebra *)
		MapThread[
			DefTensor[
				GiveSymbol[basis, Abs @ #1][], 
				mani, 
				PrintAs -> ColorString["\!\(\*FractionBox[\(\[PartialD]\), \(\[PartialD]" <> PrintAs[Evaluate @ (Head @ #2)] <> "\)]\)", BasisColor @ basis], 
				GradeOfTensor -> {CircleTimes -> 1}, 
				Master -> basis
			]&, {numbers, scalars}
		]; 
		(* Assign values *)
		Thread[ComponentValue[dx[ManifoldOfChart @ basis][{#, basis}]& /@ CNumbersOf @ basis, Diff /@ ScalarsOfChart @ basis]]; 
		Thread[ComponentValue[PDFrame[mani][{#, -basis}]& /@ numbers, GiveSymbol[basis, #][]& /@ numbers]]
	]


xTension["xTerior`", DefBasis, "End"] :=
	defFreeAlgebraBasis; 

defFreeAlgebraBasis[basis_, __] :=
	With[{numbers = CNumbersOf @ basis, mani = BaseOfVBundle @ VBundleOfBasis
		 @ basis},
		If[Not @ ChartQ @ basis,
			(* Define the algebra elements forming the basis of the free algebra*)
			DefTensor[GiveSymbol[basis, Abs @ #][], mani, GradeOfTensor -> {CircleTimes-> 1}, Master -> basis]& /@ numbers;
			(* Assign values *)
			Thread[ComponentValue[eFrame[mani][{#, -basis}]& /@ numbers, GiveSymbol[basis, #][]& /@ numbers]]
		]
	]

(**************************** 2.8 The Hodge dual *********************************)


(* The third basic ingredient of exterior algebra is Hodge duality, which requires the presence of a metric. *)

DefInertHead[Hodge[metric_], LinearQ -> True, ContractThrough -> {delta}, DefInfo -> Null]

Hodge /: PrintAs @ Hodge[metric_] :=
	If[Head[metric] === CTensor,
		"*",
		"\!\(\*SubscriptBox[\(*\), \(" <> PrintAs[metric] <> "\)]\)"
	]; 


(* Thread over lists and equal *)
Hodge[metric_][expr_List] :=
	Hodge[metric][#]& /@ expr; 

Hodge[metric_][expr_Equal] :=
	Hodge[metric][#]& /@ expr; 


(* Hodge of a CTensor object *)
Hodge[metric_][CTensor[array_, bases_List, weight_][inds__]] :=
	CTensor[Hodge[metric][array], bases, weight][inds]; 


(* Hodge of the product of two objects *)
Hodge[metric_][x_ y_] :=
	x Hodge[metric][y] /; Grade[x, Wedge] === 0


(* [Jose: This previous definition overlaps with LinearQ->True, so we might want to rethink that option in relation with the products and ScalarsOfProduct.] *)
DimOfMetric[metric_] :=
	DimOfVBundle[VBundleOfMetric[metric]]


Hodge /: Grade[Hodge[metric_][expr_], Wedge] :=
	DimOfMetric[metric] - Grade[expr, Wedge]


Hodge[metric_] @ Hodge[metric_][expr_] :=
	(-1)^(Grade[expr, Wedge] (DimOfMetric[metric] - 1) + SignatureOfMetric[metric][[2]]) expr


(* This function converts Hodge duals of product of the coframe basis (either holonomic or non-holonomic): *)
(* Expand dual of differentials of coordinate elements *)

ExpandHodgeDual[expr_, dx[mani_?ManifoldQ], met_] :=
	ExpandHodgeDual1[(expr /. Reverse /@ Flatten[List @@ TensorValues @ dx[mani]]), dx[mani], met]; 

(* Expand of the wedge product of canonical 1-forms *)

ExpandHodgeDual[expr_, (coframe:(Coframe | dx))[mani_?ManifoldQ], met_] :=
	ExpandHodgeDual1[expr, coframe[mani], met]; 

ExpandHodgeDual1[expr_, (coframe:(Coframe | dx))[mani_?ManifoldQ], met_] :=
	expr /.
		{
			HoldPattern[Hodge[met][form:Wedge[coframe[mani][_]..]] | form:Hodge[met][coframe[mani][_]]] :>
				With[{dim = DimOfMetric[met], n = Length[form], inds = Sequence @@@List @@ form},
					With[{dummies = DummyIn /@ ConstantArray[VBundleOfMetric[met], dim- n]},
						1 / (dim - n)! epsilon[met] @@ Join[inds, ChangeIndex /@ dummies] Wedge @@ (coframe[mani] /@ dummies)
					]
				],
			HoldPattern[Hodge[met][form_]] :>
				With[{dim = DimOfMetric[met]},
					With[{dummies = DummyIn /@ ConstantArray[VBundleOfMetric[met], dim]},
						form / (dim)! epsilon[met] @@ (ChangeIndex /@ dummies) Wedge @@
							 (coframe[mani] /@ dummies)
					]
				] /; (Deg[form] === 0)
		}; 


ExpandHodgeDual1[expr_, Coframe, met_] :=
	Fold[ExpandHodgeDual1[#1, Coframe[#2], met]&, expr, $Manifolds]; 

ExpandHodgeDual1[expr_, dx, met_] :=
	Fold[ExpandHodgeDual1[#1, dx[#2], met]&, expr, $Manifolds]; 



(******************** 2.10 The co-differential ***************************)



(* We introduce the co-differential. *)
DefInertHead[Codiff[metric_], LinearQ -> True, ContractThrough -> delta, DefInfo -> Null]; 


(* Formatting and basic properties *)
Codiff /: PrintAs @ Codiff[metric_] :=
	If[Head @ metric === CTensor,
		"\[Delta]",
		"\!\(\*SubscriptBox[\(\[Delta]\), \(" <> PrintAs[metric] <> "\)]\)"
	]; 

Codiff /: Grade[Codiff[metric_][expr_, ___], Wedge] := -1 + Grade[expr, Wedge]


Codiff /: MakeBoxes[Codiff[metric_][form_, PD?CovDQ], StandardForm] :=
	xAct`xTensor`Private`interpretbox[
		Codiff[metric][form, PD], 
		RowBox[{
			PrintAs[Codiff[metric]], "[", MakeBoxes[form, StandardForm], "]"
		}]
	]; 

Codiff /: MakeBoxes[Codiff[metric_][form_, cd_?CovDQ], StandardForm] :=
	xAct`xTensor`Private`interpretbox[
		Codiff[metric][form, cd], 
		RowBox[{
			SuperscriptBox[PrintAs[Codiff[metric]], Last @ SymbolOfCovD[cd]], 
			"[",
			MakeBoxes[form, StandardForm], "]"}
		]
	]; 


(* We always have PD as the covariant derivative. *)

HoldPattern[Codiff[met_][expr_]] :=
	Codiff[met][expr, PD]; 


(* Expansion of the co-differential in terms of the Hodge dual *)
CodiffToDiff[expr_] :=
	expr //. Codiff[met_][expr1_, covd_?CovDQ] :> 
		(-1)^(DimOfMetric[met] Grade[expr1, Wedge] + DimOfMetric[met] + SignatureOfMetric[met][[2]]) Hodge[met] @ Diff[Hodge[met] @ expr1, covd]


(* For convenience we program this identity: *)
Codiff[metric_][Codiff[metric_][expr_, PD]] :=0


(* Co-differential of basis objects is zero. *)
Codiff[_Basis, covd_?CovDQ] :=0; 


(* Thread over lists and equal *)
Codiff @ expr_List :=
	Codiff /@ expr; 

Codiff @ expr_Equal :=
	Codiff /@ expr; 


(******************** 2.11 Poincar\[EAcute] lemma ************************)



(* We implement the computation of the potential for an exact form (Poincar\[EAcute] lemma) *)
FindPotential[expr_Plus, point_List, chart_?ChartQ, options:OptionsPattern[Integrate]] :=
	FindPotential[#, point, chart, options]& /@ expr; 

FindPotential[expr_Times, point_List, chart_?ChartQ, options:OptionsPattern[Integrate]] :=
	FindPotential[Expand @ expr, point, chart, options]; 


(*Simplest cases for grade 1 forms *)

FindPotential[expr_Diff, point_List, chart_?ChartQ, options:OptionsPattern[Integrate]] :=
	Part[expr, 1]; 

FindPotential[factor_ expr_Diff, point_List, chart_?ChartQ, options:OptionsPattern[Integrate]] :=
	Integrate1[(factor /. Thread[Rule[ScalarsOfChart @ chart, Times[#, t]& /@ (ScalarsOfChart @ chart - point) + point]]) 
		(Part[expr, 1] - Part[point, First @ Flatten @ Position[ScalarsOfChart @ chart, Part[expr, 1]]]), 
		{t, 0, 1}, 
		options
	]; 

(* Do the actual integration *)

Integrate1 /: HoldPattern @ Plus[var__Integrate1] :=
	Integrate[Plus @@ First /@ {var}, {t, 0, 1}]; 


(* Poincare Lemma for higher degree forms *)

FindPotential[factor_. expr_Wedge, point_List, chart_?ChartQ, options:OptionsPattern[Integrate]] :=
	Integrate[(factor /. Thread[Rule[ScalarsOfChart @ chart, Times[#, t]& /@ (ScalarsOfChart @ chart - point) + point]]) 
		Sum[(-1)^(i - 1) t^(Deg @ expr - 1) *
			(Part[expr, i, 1] - 
			Part[
				point, First @ Flatten @ Position[ScalarsOfChart @ chart, 
				Part[expr, i, 1]]
			]) *
			Delete[expr, {i}], {i, 1, Length[expr]}
		], 
		{t, 0, 1}, 
		options
	]; 

(********************* 2.12. Integration of differential forms *******************)

(* Formatting of integration of differential forms: *)
FormIntegrate::dim = "Degree of form `1` is incompatible with dimension of manifold `2`."; 


InertHeadQ[FormIntegrate] ^:=
	True; 

LinearQ[FormIntegrate] ^:=
	True; 

FormIntegrate[expr_Plus, man_?ManifoldQ] :=
	FormIntegrate[#, man]& /@ expr; 

FormIntegrate[form_, EmptyManifold] :=0; 

FormIntegrate[c_?ConstantQ, man_?ManifoldQ] :=
	c /; DimOfManifold[man] == 0; 

FormIntegrate[form_, man_?ManifoldQ] :=
	If[form === 0,
		0,
		Throw[Message[FormIntegrate::dim, form, man]];
		$Failed
	] /; Deg[form] =!= DimOfManifold[man]; 

ToCanonical[FormIntegrate[form_, man_], opts___] ^:=
	FormIntegrate[ToCanonical[form, opts], man]; 

FormIntegrate /: Grade[FormIntegrate[form_, man_], Wedge] :=
	0 /; Deg[form] === DimOfManifold[man]; 


(* Formatting. Do not remove the ugly ?InertHeadQ check here. It would break typesetting *)

FormIntegrate /: MakeBoxes[FormIntegrate?InertHeadQ[form_, man_], StandardForm] :=
	RowBox[{SubscriptBox["\[Integral]", MakeBoxes[man, StandardForm]], 
	MakeBoxes[form,
		 StandardForm]}
]; 

(* This doesn't work *)

MakeExpression[RowBox[{SubscriptBox["\[Integral]", man_], form_}], StandardForm] :=
	FormIntegrate[MakeExpression[form, StandardForm], MakeExpression[man,StandardForm]]; 


(* Use the Stokes theorem. By default it converts n-dimensional integration into (n-1)-dimensional integration, but using the two-argument form, UseStokes can work in both directions: *)

UseStokes[expr_] :=
	expr /. HoldPattern[FormIntegrate[Diff[form_, PD], man_]] :> FormIntegrate[form, ManifoldBoundary[man]]; 

UseStokes[expr_, man_?ManifoldQ] :=
	expr /. HoldPattern[FormIntegrate[Diff[form_, PD], man]] :> FormIntegrate[form, ManifoldBoundary[man]]; 

UseStokes[expr_, form_] :=
	expr /. HoldPattern[FormIntegrate[form, ManifoldBoundary[man_]]] :> FormIntegrate[Diff[form], man]; 


(************************+ 2.13 Koszul operator and its properties ***************************)

(* Koszul operator is a xTension of DefCovD:*)
(* xTensions *)

(* Actual definition of the Koszul operator *)

xTension["xTerior`", DefCovD, "End"] :=
	defKoszulCovD; 

(* Properties of the Koszul operator when acting on the elements of a chart *)

xTension["xCoba`", DefChart, "End"] :=
	setKoszulValue; 


(* Actual definition of the Koszul operator. In addition we also define the corresponding "index free covariant derivative": *)
defKoszulCovD[covd_?CovDQ[_], __] :=
	With[{covdsymbol = GiveSymbol[Koszul, covd], covdsymbolf = GiveSymbol[CovD, covd]},

		(* Definition of Koszul operator *)
		xAct`xTerior`Private`DefGradedDerivation[covdsymbol[v_], CircleTimes, 0, PrintAs -> Last @ SymbolOfCovD[covd], Master -> covd];

		(* Definition of the "index free covariant derivative"*)
		DefInertHead[covdsymbolf, PrintAs -> Last @ SymbolOfCovD[covd], LinearQ-> True, Master -> covd];

		(* Properties of the Koszul operator *)
		HoldPattern[covdsymbol[v_][expr_?ConstantQ]] := 0;
		HoldPattern[covdsymbol[factor_ v_][expr_]] := factor covdsymbol[v][expr] /; Grade[factor, CircleTimes] === 0;
		HoldPattern[covdsymbol[v_][expr1_ expr2_]] := expr2 covdsymbol[v][expr1] + expr1 covdsymbol[v][expr2];

		(* Properties of the "index free covariant derivative" *)
		HoldPattern[covdsymbolf[expr_?ConstantQ]] := 0;
		HoldPattern[covdsymbolf[expr_]] := Diff[expr] /; Grade[expr, CircleTimes] === 0;
		HoldPattern[covdsymbolf[expr1_ expr2_]] := CircleTimes[covdsymbolf[expr1], expr2] + CircleTimes[covdsymbolf[expr2], expr1];
		HoldPattern[covdsymbolf[CircleTimes[expr1_, expr2_]]] := CircleTimes[covdsymbolf[expr1], expr2] + CircleTimes[expr1, covdsymbolf[expr2]]; 
		
		CircleTimes /: HoldPattern @ Grade[covdsymbolf[expr_], CircleTimes]:= Grade[expr, CircleTimes] + 1;
		Wedge /: HoldPattern[Grade[covdsymbolf[expr_], Wedge]] := Grade[expr,Wedge] + 1;
		CircleTimes /: HoldPattern[Grade[covdsymbolf[expr_], CircleTimes]] := Grade[expr, CircleTimes] + 1;

		HoldPattern[covd[{a_?NumberQ, basis_}] @ expr_] :=
			Which[
				ChartQ[basis] && Grade[expr, CircleTimes] >= 1,
					covdsymbol[PDFrame[ManifoldOfCovD @ covd][{a, basis}]][expr],
				BasisQ[basis] && Grade[expr, CircleTimes] >= 1,
					covdsymbol[eFrame[ManifoldOfCovD @ covd][{a, basis}]][expr],
				True,
					$Failed
			];
		HoldPattern[covd[{a_?NumberQ, -basis_}] @ expr_] :=
			Which[
				ChartQ[basis] && Grade[expr, CircleTimes] >= 1,
					covdsymbol[PDFrame[ManifoldOfCovD @ covd][{a, -basis}]][expr],
				BasisQ[basis] && Grade[expr, CircleTimes] >= 1,
					covdsymbol[eFrame[ManifoldOfCovD @ covd][{a, -basis}]][expr],
				True,
					$Failed
			]
	]


(* Action of the Koszul operator on the holonomic elements of a coordinate chart: *)
setKoszulValue[chart_?ChartQ, __] :=
	With[{covd = PDOfBasis[chart]},
		With[{sym0 = GiveSymbol[Koszul, covd], sym1 = GiveSymbol[PD, chart], mani = ManifoldOfChart[chart]},
			Outer[
				Set[sym0[PDFrame[mani][{#1, -chart}]] @ Part[ScalarsOfChart[chart], #2], 
				sym1[{#1, -chart}][Part[ScalarsOfChart[chart], #2]]]&, CNumbersOf@ chart, 
				Range[1, Length @ ScalarsOfChart @ chart]
			];
			Outer[Set[sym0[GiveSymbol[chart, #1][]] @ Part[ScalarsOfChart[chart], #2], 
				sym1[{#1, -chart}][Part[ScalarsOfChart[chart], #2]]]&, 
				CNumbersOf@ chart, 
				Range[1, Length @ ScalarsOfChart @ chart]
			]
		]
	]


(* Expansion of the action of Koszul operator on frame elements: *)

ExpandKoszul[expr_, covd1_?CovDQ, covd2_?CovDQ] :=
	expr /.
		{
			GiveSymbol[Koszul, covd1][eFrame[mani_?ManifoldQ][a_]][eFrame[mani_
				?ManifoldQ][b_]] :>
				Module[{c = DummyIn[Tangent[mani]]},
					Christoffel[covd1, covd2][c, a, b] eFrame[mani][-c]
				]
			,
			GiveSymbol[Koszul, covd1][eFrame[mani_?ManifoldQ][a_]][Coframe[mani_
				?ManifoldQ][b_]] :>
				Module[{c = DummyIn[Tangent[mani]]},
					-Christoffel[covd1, covd2][b, a, -c] Coframe[mani][c]
				]
			,
			GiveSymbol[Koszul, covd1][PDFrame[mani_?ManifoldQ][a_]][PDFrame[mani_
				?ManifoldQ][b_]] :>
				Module[{c = DummyIn[Tangent[mani]]},
					Christoffel[covd1, covd2][c, a, b] PDFrame[mani][-c]
				]
			,
			GiveSymbol[Koszul, covd1][PDFrame[mani_?ManifoldQ][a_]][dx[mani_?ManifoldQ
				][b_]] :>
				Module[{c = DummyIn[Tangent[mani]]},
					-Christoffel[covd1, covd2][b, a, -c] dx[mani][c]
				]
		}; 

(* Expansion with respect to a coordinate chart of the generalized covariant derivative of rank 1 tensors *)

ExpandCovD[expr_, covd_?CovDQ, chart_?ChartQ] :=
	With[{covdf = GiveSymbol[CovD, covd], M = ManifoldOfChart[chart], covdk = GiveSymbol[Koszul, covd]},
		(expr /. HoldPattern[covdf[expr1_]] :> Plus @@ (dx[M][{#, chart}] \[CircleTimes]
			 covdk[PDFrame[M][{#, -chart}]][expr1]& /@ CNumbersOf @ chart) /; Grade[
			expr1, CircleTimes] === 1)
	]

(******************************************************************************)

(********************* 3. Graded derivations **********************************)

(******************************************************************************)


(****************** 3.1 Connection forms **********************************)

(*

There are two objects we need for covariant exterior derivatives. First, 1-forms taking values in VB\[CircleTimes]-VB which represent the tensorial difference between two connections. 
Second, 2-forms taking values in VB\[CircleTimes]-VB which represent the curvature of an individual connection. 
We'll use nonatomic heads for both of these. The notation will be ConnectionForm[CD1,CD2,VB][A,-B] CurvatureForm[CD,VB][A,-B].

From a mathematical point of view the connection form should not be regarded as the tensorial difference between two connections. 
The reason for this is that the connection form is defined as a 1-form in the frame bundle regarded as a differentiable manifold 
and in this sense it is truly tensorial. One loses the tensorial character when one does a splitting of the frame bundle into the 
base manifold and the fibres. For this reason we use the following notation for the connection form:

ConnectionForm[CD,VB][A,-B]

This notation has been implemented in this notebook.

First, making these nonatomic-head tensors. We start by the general connection form.  What happens with the MastersOf a non-atomic symbol? 
Given that a non-atomic tensor cannot be undefined I guess that it makes no sense to ask about its Masters,  right ?  
Also I guess that it cannot be Servant of anything either.

*)


xTensorQ[ConnectionForm[cd_?CovDQ, _]] ^= True; 

SlotsOfTensor[ConnectionForm[_, vb_?VBundleQ]] ^:=
	{vb, -vb}; 

ConnectionForm /: GradeOfTensor[ConnectionForm[_, _], Wedge] = 1; 

SymmetryGroupOfTensor[ConnectionForm[_, _]] ^= StrongGenSet[{}, GenSet[]]; 

Dagger[ConnectionForm[cd1_, vb_]] ^:=
	ConnectionForm[cd1, Dagger @ vb]; 

DefInfo[ConnectionForm[_, _]] ^= {"nonsymmetric Connection 1-form", ""}; 

DependenciesOfTensor[ConnectionForm[cd1_, _]] ^:=
	Union @@ DependenciesOfCovD /@ {cd1}; 

HostsOf[ConnectionForm[cd1_, vb_]] ^:=
	Join[{cd1}, Union @@ HostsOf /@ {cd1, vb}];  (* Should we put Union@@HostsOf/@{cd1,vb} here? Yes but we need to add cd1 itself *)

TensorID[ConnectionForm[_, _]] ^= {}; 

PrintAs[ConnectionForm] ^= "A"; 

PrintAs[ConnectionForm[cd1_, _]] ^:=
	PrintAs[ConnectionForm] <> "[" <> Last @ SymbolOfCovD[cd1] <> "]"; 

PrintAs[ConnectionForm[PD, _]] ^:=
	PrintAs[ConnectionForm]; 


(* We introduce now the ChristoffelForm (connection form for a connection in the frame bundle).  *)
xTensorQ[ChristoffelForm[cd1_?CovDQ]] ^= True; 

SlotsOfTensor[ChristoffelForm[cd1_?CovDQ]] ^:=
	{Tangent @ ManifoldOfCovD @ cd1, -Tangent @ ManifoldOfCovD @ cd1}; 

ChristoffelForm /: GradeOfTensor[ChristoffelForm[_], Wedge] = 1; 

SymmetryGroupOfTensor[ChristoffelForm[_]] ^= StrongGenSet[{}, GenSet[]]; 

Dagger[ChristoffelForm[cd1_]] ^:=
	ChristoffelForm[cd1]; 

DefInfo[ChristoffelForm[_]] ^= {"nonsymmetric frame bundle Connection 1-form",""}; 

DependenciesOfTensor[ChristoffelForm[cd1_]] ^:=
	Union @@ DependenciesOfCovD /@ {cd1}; 

HostsOf[ChristoffelForm[cd1_]] ^:=
	Join[{cd1}, Union @@ HostsOf /@ {cd1}];  (* Should we put Union@@HostsOf/@{cd1,cd2,vb} here? Yes, but addint also cd1 *)

TensorID[ChristoffelForm[_]] ^= {}; 

PrintAs[ChristoffelForm] ^= "\[CapitalGamma]"; 

(PrintAs[ChristoffelForm[cd1_]] /; Head @ cd1 =!= CCovD) ^:=
	PrintAs[ChristoffelForm] <> "[" <> Last @ SymbolOfCovD[cd1] <> "]"; 

PrintAs[ChristoffelForm[PD]] ^:=
	PrintAs[ChristoffelForm]; 


(* Connection form when the metric is given as a CTensor *)
ChristoffelForm[exr_CCovD] :=
	Head @
		Module[{ind = DummyIn @ VBundleOfBasis[-First @ Part[Last @ exr, 2]], a1, a2},
			{a1, a2} = GetIndicesOfVBundle[VBundleOfBasis[-First @ Part[Last @exr, 2]], 2];
			Part[exr, 2][a1, -ind, -a2] ToCTensor[Coframe[BaseOfVBundle @ VBundleOfBasis[-First @ Part[Part[exr, 3], 2]]], {-First @ Part[Part[exr, 3], 2]}][ind]
		]; 


(* Automatic replacement of the connection form by the Christoffel form when the connection corresponds to a connection in a frame bundle. *)
ConnectionForm[cd1_, vb_] :=
	ChristoffelForm[cd1] /; (Tangent @ ManifoldOfCovD @ cd1 === vb); 

ChristoffelForm[cd_, tangentbundle_] :=
	ChristoffelForm[cd]; 


(* PD is regarded as a trivial connection in a principal bundle. Hence: *)
ConnectionForm[PD, vb_] :=Zero; 

ChristoffelForm[PD] :=Zero; 


(* Transformation of the connection form into the connection tensor plus co-frame. *)
ConnectionFormToTensor[expr_, covd_, frame:(Coframe | dx)] :=
	expr /.
		{
			ChristoffelForm[cd1_][ind1_, ind2_] :>
				Module[{a = DummyIn @ First @ VBundlesOfCovD @ covd, i1, i2},
					i1 = NewIndexIn[VBundleOfIndex[ind1]];
					i2 = NewIndexIn[VBundleOfIndex[ind2]];
					If[xTensorQ @ GiveSymbol[Christoffel, cd1, covd] === False,
						DefTensor[GiveSymbol[Christoffel, cd1, covd][i1, -a, i2], ManifoldOfCovD
							 @ covd]
					];
					GiveSymbol[Christoffel, cd1, covd][ind1, -a, ind2] frame[ManifoldOfCovD
						 @ covd][a]
				] /; covd =!= PD
			,
			ConnectionForm[cd1_, vbundle_][ind1_, ind2_] :>
				Module[{a = DummyIn @ Tangent @ BaseOfVBundle @ vbundle, i1, i2},
					i1 = NewIndexIn[VBundleOfIndex[ind1]];
					i2 = NewIndexIn[VBundleOfIndex[ind2]];
					If[xTensorQ @ GiveSymbol[AChristoffel, covd, cd1] === False,
						DefTensor[GiveSymbol[AChristoffel, covd, cd1][i1, -a, i2], BaseOfVBundle
							 @ vbundle]
					];
					GiveSymbol[AChristoffel, covd, cd1][ind1, -a, ind2] frame[BaseOfVBundle
						 @ vbundle][a]
				] /; covd =!= PD
		}; 


(* If the covariant derivative is PD then the code takes the following simpler form: *)
ConnectionFormToTensor[expr_, PD, frame:(Coframe | dx)] :=
	expr /.
		{
			ChristoffelForm[cd1_][ind1_, ind2_] :>
				Module[{a = DummyIn @ First @ VBundlesOfCovD @ cd1},
					GiveSymbol[Christoffel, cd1][ind1, -a, ind2] frame[ManifoldOfCovD
						 @ cd1][a]
				]
			,
			ConnectionForm[cd1_, _][ind1_, ind2_] :>
				Module[{a = DummyIn @ First @ VBundlesOfCovD @ cd1},
					GiveSymbol[AChristoffel, cd1][ind1, -a, ind2] frame[ManifoldOfCovD
						 @ cd1][a]
				]
		}

(************ 3.2 Curvature forms *************************)

(* Now for the curvature form: *)
xTensorQ[CurvatureForm[_, _]] ^= True; 

SlotsOfTensor[CurvatureForm[_, vb_?VBundleQ]] ^:=
	{vb, -vb}; 

CurvatureForm /: GradeOfTensor[CurvatureForm[_, _], Wedge] = 2; 

SymmetryGroupOfTensor[CurvatureForm[_, _]] ^= StrongGenSet[{}, GenSet[
	]]; 

Dagger[CurvatureForm[cd_, vb_]] ^:=
	CurvatureForm[cd, Dagger @ vb]; 

DefInfo[CurvatureForm[_, _]] ^= {"Curvature 2-form", ""}; 

DependenciesOfTensor[CurvatureForm[cd_, _]] ^:=
	DependenciesOfCovD[cd]; 

HostsOf[CurvatureForm[cd_, vb_]] ^:=
	Join[{cd}, Union @@ HostsOf /@ {cd, vb}];  (* Should we put Union@@HostsOf/@{cd,vb} here? Yes but we need to add cd itself *)

TensorID[CurvatureForm[_, _]] ^= {}; 

PrintAs[CurvatureForm] ^= "F"; 

PrintAs[CurvatureForm[cd_, _]] ^:=
	PrintAs[CurvatureForm] <> "[" <> Last @ SymbolOfCovD[cd] <> "]"; 


(* Case of a connection in the frame bundle (RiemannForm) *)
xTensorQ[RiemannForm[_]] ^= True; 

SlotsOfTensor[RiemannForm[cd_?CovDQ]] ^:=
	{Tangent @ ManifoldOfCovD @ cd, -Tangent @ ManifoldOfCovD @ cd}; 

RiemannForm /: GradeOfTensor[RiemannForm[_], Wedge] = 2; 

SymmetryGroupOfTensor[RiemannForm[cd_?CovDQ]] ^:=
	If[MetricOfCovD @ cd =!= Null,
		Antisymmetric[{1, 2}, Cycles],
		StrongGenSet[{}, GenSet[]]
	]; 

Dagger[RiemannForm[cd_]] ^:=
	RiemannForm[cd]; 

DefInfo[RiemannForm[_]] ^= {"Curvature 2-form in the frame bundle", ""}; 

DependenciesOfTensor[RiemannForm[cd_]] ^:=
	DependenciesOfCovD[cd]; 

HostsOf[RiemannForm[cd_]] ^:=
	Join[{cd}, Union @@ HostsOf /@ {cd}];  (* Should we put Union@@HostsOf/@{cd,vb} here? Yes but we need to add cd itself *)

TensorID[RiemannForm[_]] ^= {}; 

PrintAs[RiemannForm] ^= "R"; 

(PrintAs[RiemannForm[cd_]] /; Head @ cd =!= CCovD) ^:=PrintAs[RiemannForm] <> "[" <> Last @ SymbolOfCovD[cd] <> "]"; 


(* Riemann form of a CTensor expression *)
RiemannForm[exr_CCovD] :=
	If[Riemann[exr] =!= Zero,
		Head @
			Module[{ind1 = DummyIn @ VBundleOfBasis[-First @ Part[Last @ exr, 2]], ind2 = DummyIn @ VBundleOfBasis[-First @ Part[Last @ exr, 2]], a1,a2},
				{a1, a2} = GetIndicesOfVBundle[VBundleOfBasis[-First @ Part[Last @ exr, 2]], 2];
				Riemann[exr][-ind1, -ind2, -a1, a2] 
				Wedge[
					ToCTensor[Coframe[BaseOfVBundle@ VBundleOfBasis[-First @ Part[Part[exr, 3], 2]]], {-First @ Part[Part[exr, 3], 2]}][ind1], 
					ToCTensor[Coframe[BaseOfVBundle @ VBundleOfBasis[-First @ Part[Part[exr, 3], 2]]], {-First @ Part[Part[exr, 3], 2]}][ind2]
				]
			],
		Zero
	]; 


(* Automatic replacement of the curvature form if the covariant derivative corresponds to a covariant derivative defined in the tangent bundle. *)

CurvatureForm[cd_?CovDQ, vbundle_] :=
	RiemannForm[cd] /; vbundle === Tangent @ ManifoldOfCovD @ cd; 

RiemannForm[cd_, vbundle_] := RiemannForm[cd]; 


(* Transformation of the Curvature form into the curvature tensor plus co-frame. *)

CurvatureFormToTensor[expr_, frame:(Coframe | dx)] :=
	expr /.
		{
			HoldPattern @ CurvatureForm[cd1_, vbundle1_?VBundleQ][inds__] :>
				Module[{a = DummyIn @ First @ VBundlesOfCovD @ cd1, b = DummyIn @
					 First @ VBundlesOfCovD @ cd1},
					1 / 2 GiveSymbol[FRiemann, cd1][-a, -b, Sequence @@ Reverse @ List@ inds] Wedge[frame[ManifoldOfCovD @ cd1][a], frame[ManifoldOfCovD @ cd1][b]]
				],
			RiemannForm[cd1_][inds__] :>
				Module[{a = DummyIn @ First @ VBundlesOfCovD @ cd1, b = DummyIn @
					 First @ VBundlesOfCovD @ cd1},
					-$RiemannSign 1 / 2 GiveSymbol[Riemann, cd1][-a, -b, Sequence @@ Reverse @ List @ inds] Wedge[frame[ManifoldOfCovD @ cd1][a], frame[ManifoldOfCovD @ cd1][b]]
				]
		}


(********************* 3.3 The torsion 2-form ********************************)



(* The torsion 2-form (only for covariant derivatives on the tangent bundle) *)

xTensorQ[TorsionForm[_]] ^= True; 

SlotsOfTensor[TorsionForm[cd_]] ^:= {Tangent @ ManifoldOfCovD @ cd}; 

TorsionForm /: GradeOfTensor[TorsionForm[_], Wedge] = 2; 

SymmetryGroupOfTensor[TorsionForm[_]] ^= StrongGenSet[{}, GenSet[]]; 

Dagger[TorsionForm[cd_]] ^:= TorsionForm[cd]

DefInfo[TorsionForm[_]] ^= {"Torsion 2-form", ""}; 

DependenciesOfTensor[TorsionForm[cd_]] ^:= DependenciesOfCovD[cd]; 

HostsOf[TorsionForm[cd_]] ^:= HostsOf @ cd;  (* Should we put HostsOf@cd here? OK*)

TensorID[TorsionForm[_]] ^= {}; 


(* ::Input::Initialization:: *)
PrintAs[TorsionForm] ^= "\[GothicCapitalT]"; 

(PrintAs[TorsionForm[cd_]] /; Head @ cd =!= CCovD) ^:= PrintAs[TorsionForm] <> "[" <> Last @ SymbolOfCovD[cd] <> "]"; 


(* Torsion form of a covariant derivative of a metric given as a CTensor *)

TorsionForm[exr_CCovD] :=
	Head @
		Module[{ind1 = DummyIn @ VBundleOfBasis[-First @ Part[Last @ exr, 2]], ind2 = DummyIn @ VBundleOfBasis[-First @ Part[Last @ exr, 2]], a1},
			{a1} = GetIndicesOfVBundle[VBundleOfBasis[-First @ Part[Last @ exr, 2]], 1];
			Torsion[exr][a1, -ind1, -ind2]
			Wedge[
				ToCTensor[
					Coframe[BaseOfVBundle @ VBundleOfBasis[-First @ Part[Part[exr, 3], 2]]], 
					{-First @ Part[Part[exr, 3], 2]}
				][ind1], 
				ToCTensor[
					Coframe[BaseOfVBundle @ VBundleOfBasis[-First @ Part[Part[exr, 3], 2]]], 
					{-First @ Part[Part[exr, 3], 2]}
				][ind2]
			]
		]; 

(********************* 3.4 Change of exterior covariant derivative and generalized covariant derivative **************************)

(*
We will use a convention that there should be a correspondence 
  AChr1[CD1,CD2,VB][A,-B] -> AChristoffel[CD1,CD2][A,-\[Mu],-B] \[Theta][\[Mu]]
  FRiem2[CD,VB][A,-B] -> 1/2 FRiemann[CD][-\[Mu],-\[Nu],-B,A] \[Theta][\[Mu]]\[Wedge]\[Theta][\[Nu]]          
  (notice different ordering of A,B on two sides)
  Torsion2[CD][\[Alpha]] -> 1/2 Torsion[CD][\[Alpha],-\[Beta],-\[Gamma]] \[Theta][\[Beta]]\[Wedge]\[Theta][\[Gamma]]
  
We'll want to be able to
  1) Change from diff[form,CD1] to diff[form,CD2], introducing AChr1[CD1,CD2,...] as needed
  2) Compute diff[diff[form,CD],CD], introducing FRiem2[CD,...] as needed
  3) Change from FRiem2[CD1,VB] to FRiem2[CD2,VB], introducing diff[AChr1[CD...]] and Wedge[AChr1,AChr1]
*)

(*
1. Change from diff[form,CD1] to diff[form,CD2]. The formula is:
  d^CD1 v^A-d^CD2 v^A  = \!\(
\*SubsuperscriptBox[\(A[CD1, CD2, VB]\), \(\(\ \ \ \)\(B\)\), \(A\)]\[Wedge]
\*SuperscriptBox[\(v\), \(B\)]\)     (where VB is the VBundle of A; plus sign for +VB, minus sign for -VB)
TODO: What about densities? Alfonso: I'm not sure whether the distinction "true tensor" / "tensor density" holds in this context.
*)

(*
TODO: describe the action of the "Generalized Covariant derivative".
*)

(* ConnectionForm with two covariant derivatives as arguments represents 
the difference between two connection forms which carry only a covariant derivative in the argument. *)
ConnectionForm[cd1_, cd2_, vbundle_][inds__] :=
	ConnectionForm[cd1, vbundle][inds] - ConnectionForm[cd2, vbundle][inds]; 


(* ::Input::Initialization:: *)
ChangeExtD[expr_, cd_?CovDQ, cd_] :=
	expr; 

ChangeExtD[expr_, cd1_?CovDQ, cd2_:PD] :=
	expr /. HoldPattern[Diff[expr1_, cd1]] :> makeChangeExtD[ChangeExtD[expr1, cd1, cd2], cd1, cd2]; 

ChangeExtD[expr_, list_List, covd2_:PD] :=
	Fold[ChangeExtD[#1, #2, covd2]&, expr, list]; 

ChangeExtD[expr_, x_, _:PD] :=
	Throw @ Message[ChangeExtD::unknown, "covariant derivative", x]; 

ChangeExtD[expr_] :=
	ChangeExtD[expr, $CovDs]; 


(* This only gives a correct result for an expression which is a sum of simple tensor monomials, 
 each simple term being the product of rank 1 basis elements   *)

ChangeGenCovD[expr_, cd_?CovDQ, cd_] :=
	expr; 

ChangeGenCovD[expr_, cd1_?CovDQ, cd2_:PD] :=
	With[{covdf = GiveSymbol[CovD, cd1]},
		Module[{epk},
			epk = SeparateBasis[AIndex][expr];
			If[FindIndices[expr] == IndexList[],
				Return[expr]
			];
			epk /. HoldPattern[covdf[expr1_]] :> makeChangeCovD[expr1, cd1, cd2
				]
		]
	]

ChangeGenCovD[expr_, list_List, covd2_:PD] :=
	Fold[ChangeGenCovD[#1, #2, covd2]&, expr, list]; 

ChangeGenCovD[expr_, x_, _:PD] :=
	Throw @ Message[ChangeGenCovD::unknown, "covariant derivative", x]; 

ChangeGenCovD[expr_] :=
	ChangeGenCovD[expr, $CovDs]; 


(* ::Input::Initialization:: *)
makeChangeExtD[expr_, cd1_, cd2_] :=
	With[{vbs = Apply[Union, VBundlesOfCovD /@ DeleteCases[{cd1, cd2}, PD]]},
		Diff[expr, cd2] + Plus @@ Map[addAChr1[expr, cd1, cd2], xAct`xTensor`Private`selecton[
			Select[FindFreeIndices @ expr, GIndexQ], vbs]] // ReduceAChr1
	]; 


(* ::Input::Initialization:: *)
makeChangeCovD[expr_, cd1_, cd2_] :=
	With[{vbs = Apply[Union, VBundlesOfCovD /@ DeleteCases[{cd1, cd2}, PD]], covd2 = GiveSymbol[CovD, cd2]},
		covd2[expr] + 
		Plus @@ Map[
			addAChr2[expr, cd1, cd2], 
			xAct`xTensor`Private`selecton[
				Select[
					IndexList @@ DeleteCases[List @@ FindIndices @ expr, 
						Alternatives@@ List @ FindDummyIndices @ expr
					], 
					GIndexQ
				], 
				vbs
			]
		] // ReduceAChr1
	]; 


(* ::Input::Initialization:: *)
addAChr1[expr_, cd1_, cd2_][{oldind_?xAct`xTensor`Private`upQ, dummy_}] :=
	Wedge[
		ConnectionForm[cd1, cd2, VBundleOfIndex @ oldind][oldind, -dummy], 
		ReplaceIndex[expr, oldind -> dummy]
	]

addAChr1[expr_, cd1_, cd2_][{oldind_?xAct`xTensor`Private`downQ, dummy_}] :=
	- Wedge[
		ConnectionForm[cd1, cd2, VBundleOfIndex @ oldind][dummy, oldind], 
		ReplaceIndex[expr, oldind -> -dummy]
	]


(* ::Input::Initialization:: *)
addAChr2[expr_, cd1_, cd2_][{oldind_?xAct`xTensor`Private`upQ, dummy_}] :=
	-CircleTimes[
		ConnectionForm[cd1, cd2, VBundleOfIndex @ oldind][oldind, -dummy], 
		ReplaceIndex[expr, oldind -> dummy]
	]

addAChr2[expr_, cd1_, cd2_][{oldind_?xAct`xTensor`Private`downQ, dummy_}] :=
	CircleTimes[
		ConnectionForm[cd1, cd2, VBundleOfIndex @ oldind][dummy, oldind], 
		ReplaceIndex[expr, oldind -> -dummy]
	]


(* 

Note that
  d^CD1 v^A-d^CD2 v^A  = d^CD1 v^A-d v^A+d v^A- d^CD2 v^A=\!\(
\*SubsuperscriptBox[\(A[CD1, PD, VB]\), \(\(\ \ \ \)\(B\)\), \(A\)]\[Wedge]
\*SuperscriptBox[\(v\), \(B\)]\)+\!\(
\*SubsuperscriptBox[\(A[PD, CD2, VB]\), \(\(\ \ \ \)\(B\)\), \(A\)]\[Wedge]
\*SuperscriptBox[\(v\), \(B\)]\)
and so
  \!\(
\*SubsuperscriptBox[\(A[CD1, CD2, VB]\), \(\(\ \ \ \)\(B\)\), \(A\)] = \(
\*SubsuperscriptBox[\(A[CD1, PD, VB]\), \(\(\ \ \ \)\(B\)\), \(A\)] + 
\*SubsuperscriptBox[\(A[PD, CD2, VB]\), \(\(\ \ \ \)\(B\)\), \(A\)]\)\).
Now Subsuperscript[A[CD,PD,VB],    B, A]=0  iff FreeQ[VBundlesOfCovD[CD],VB]. This allows the simplification of A[CD1,CD2,VB] when only one of {CD1,CD2} has VB in VBundlesOfCovD.

*)
ReduceAChr1[expr_] :=
	expr //. 
	{ConnectionForm[cd1_, cd2_, vb_][a_, b_] :> 
		ConnectionForm[cd1, PD, vb][a, b] /; (cd2 =!= PD && FreeQ[VBundlesOfCovD @ cd2, vb]), 
	ConnectionForm[cd1_, cd2_, vb_][a_, b_] :> 
		ConnectionForm[PD, cd2, vb][a, b] /; (cd1 =!= PD && FreeQ[VBundlesOfCovD @ cd1, vb])
	}


(* 

2. Compute diff[diff[form,CD],CD], introducing FRiem2[CD,...] as needed. The formula:
  d^\[Del] d^\[Del] v^A = \!\(
\*SubsuperscriptBox[\(F[\[Del]\(\(,\)\(VB\)\)]\), \(\(\ \ \ \)\(B\)\), \(A\)]\[Wedge]
\*SuperscriptBox[\(v\), \(B\)]\)             (where VB is the VBundle of A; plus sign for +VB, minus sign for -VB) 

*)
Diff[Diff[expr_, PD], PD] := 0

HoldPattern @ Diff[HoldPattern[Diff[expr_, cd_]], cd_] :=
	Plus @@ Map[
		addFRiem2[expr, cd], 
		xAct`xTensor`Private`selecton[Select[FindFreeIndices @ expr, AIndexQ], VBundlesOfCovD @ cd]
	]; 

addFRiem2[expr_, cd_][{oldind_?xAct`xTensor`Private`upQ, dummy_}] :=
	Wedge[
		CurvatureForm[cd, VBundleOfIndex[oldind]][oldind, -dummy], 
		ReplaceIndex[expr, oldind -> dummy]
	]; 

addFRiem2[expr_, cd_][{oldind_?xAct`xTensor`Private`downQ, dummy_}] :=
	-Wedge[
		CurvatureForm[cd, VBundleOfIndex[oldind]][dummy, oldind], 
		ReplaceIndex[expr, oldind -> -dummy]
	]


(* Just in case, add the following definitions: *)
CurvatureForm[PD, _] := Zero; 


RiemannForm[PD] :=
	Zero; 


(* 
3. Change from FRiem2[CD1,VB] to FRiem2[CD2,VB],  introducing diff[AChr1[CD...]] and Wedge[AChr1,AChr1]. The formula:
  \!\(
\*SubsuperscriptBox[\(F[CD2, VB]\), \(\(\ \ \ \)\(B\)\), \(A\)] = \(
\*SubsuperscriptBox[\(F[CD1, VB]\), \(\(\ \ \ \)\(B\)\), \(A\)]\  + \ 
\*SuperscriptBox[\(d\), \(CD1\)]
\*SubsuperscriptBox[\(A[CD2, CD1, VB]\), \(\(\ \ \ \ \)\(B\)\), \(A\)] + 
\*SubsuperscriptBox[\(A[CD2, CD1, VB]\), \(\(\ \ \ \ \)\(C\)\), \(A\)]\[Wedge]
\*SubsuperscriptBox[\(A[CD2, CD1, VB]\), \(\(\ \ \ \ \)\(B\)\), \(C\)]\)\)
*)

ChangeCurvatureForm[expr_, cd_, cd_] :=
	expr; 

ChangeCurvatureForm[expr_, cd1_?CurvatureQ, cd2_:PD] :=
	ReduceAChr1[expr /. changeCurvatureFormRules[cd1, cd2]]; 

ChangeCurvatureForm[expr_, list_List, cd2_:PD] :=
	Fold[ChangeCurvatureForm[#1, #2, cd2]&, expr, list]; 

ChangeCurvatureForm[expr_, _, _:PD] :=
	expr; 

ChangeCurvatureForm[expr_] :=
	ChangeCurvatureForm[expr, $CovDs]; 


changeCurvatureFormRules[cd2_, cd1_] :=
	{
		CurvatureForm[cd2, vb_?VBundleQ][a_, b_] :>
			With[{c = DummyIn @ vb, A1 = ConnectionForm[cd2, cd1, vb]},
				CurvatureForm[cd1, vb][a, b] + Diff[A1[a, b], cd1] + Wedge[A1[a, -c], A1[c, b]]
			],
		RiemannForm[cd2][a_, b_] :>
			With[{c = DummyIn @ Tangent @ ManifoldOfCovD @ cd2, A1 = ChristoffelForm[cd2, cd1]},
				RiemannForm[cd1][a, b] + Diff[A1[a, b], cd1] + Wedge[A1[a, -c], A1[c, b]]
			]
	}; 

(************ 3.5 UseCartan ***************************)



(* Thread over equations and lists *)

UseCartan[expr_List, covd_] := UseCartan[#, covd]& /@ expr; 

UseCartan[expr_List, PD, {mani__?ManifoldQ}] := UseCartan[#, PD, {mani}]& /@ expr; 

UseCartan[expr_Equal, covd_] := UseCartan[#, covd]& /@ expr; 

UseCartan[expr_Equal, PD, {mani__?ManifoldQ}] := UseCartan[#, PD, {mani}]& /@ expr; 


(* Exterior derivative when covd is PD *)

UseCartan[expr_, PD, {mani__?ManifoldQ}] := (
	expr /. Diff @ expr1_ :> Module[{a = DummyIn /@ (Tangent /@ {mani})},
		Inner[dx[#1][#2] PD[-#2] @ expr1&, {mani}, a, Plus]
	] /; Deg @ expr1 === 0
); 

UseCartan[expr_, PD] := (
	expr /. Diff @ expr1_ :> Module[{a = DummyIn /@ (Tangent /@ ManifoldsOf @ expr)},
		Inner[dx[#1][#2] PD[-#2] @ expr1&, ManifoldsOf @ expr, a, Plus]
	] /; Deg @ expr1 === 0
); 


(* We encode the Cartan equations in a set of rules. We also encode the action of the exterior derivative on a function whose arguments are scalars of a coordinate chart. *)

UseCartan[expr_, covd_] := 
With[{metric = MetricOfCovD[covd], basis = BasisOfCovD[covd]},
	expr /. Flatten @ 
	{
		(* Exterior derivative of the coframe *)
		HoldPattern[Diff[Coframe[mani_][ind_?UpIndexQ], PD]] :> 
			Module[{a = DummyIn @ VBundleOfIndex @ ind},
				If[TorsionQ @ covd,
					-Wedge[ConnectionForm[covd, VBundleOfIndex @ ind][ind, -a], Coframe[mani][a]] + TorsionForm[covd][ind],
					-Wedge[ConnectionForm[covd, VBundleOfIndex @ ind][ind, -a], Coframe[mani][a]]
				]
			],
			
		(* Exterior derivative of the connection *)
		HoldPattern[Diff[(connection:(ConnectionForm | ChristoffelForm))[covd, vbundle_:Tangent @ ManifoldOfCovD @ covd][a1_?UpIndexQ, -a2_], PD]] :> 
			Module[{a = DummyIn @ VBundleOfIndex @ a1},
				CurvatureForm[covd, vbundle][a1, -a2] - Wedge[connection[covd, vbundle][a1, -a], connection[covd, vbundle][a, -a2]]
			],
			
		(* Exterior derivative of the torsion *)
		HoldPattern[Diff[TorsionForm[covd][ind_?UpIndexQ], PD]] :> 
			Module[{a = DummyIn @ VBundleOfIndex @ ind},
				Wedge[Coframe[ManifoldOfCovD @ covd][a], RiemannForm[covd][ind, -a]] - Wedge[ChristoffelForm[covd][ind, -a], TorsionForm[covd][a]]
			],
			
		(* Exterior derivative of the curvature *)
		HoldPattern[Diff[(curvature:(CurvatureForm | RiemannForm))[covd, vbundle_:Tangent @ ManifoldOfCovD @ covd][a1_?UpIndexQ, -a2_], PD]] :> 
			Module[{a = DummyIn @ VBundleOfIndex @ a1},
				Wedge[ConnectionForm[covd, vbundle][a, -a2], curvature[covd, vbundle][a1, -a]] - Wedge[curvature[covd, vbundle][a, -a2], ConnectionForm[covd, vbundle][a1, -a]]
			],
			
		(* Exterior derivative of the metric (indices downstairs) *)
		HoldPattern[Diff[metr_?MetricQ[-a1_, -a2_], PD]] :> 
			Module[{a = DummyIn @ VBundleOfIndex @ a1},
				ChristoffelForm[covd][a, -a1] metric[-a, -a2] + ChristoffelForm[covd][a, -a2] metric[-a, -a1]
			] /; MetricOfCovD[covd] === metr,
			
		(* Exterior derivative of the metric (indices upstairs) *)
		HoldPattern[Diff[metric[a1_Symbol, a2_Symbol], PD]] :> -ChristoffelForm[covd][a1, a2] - ChristoffelForm[covd][a2, a1],
		
		(* Exterior derivative when covd is the parallel derivative of a coordinate chart *)
		HoldPattern @ Diff[expr1_?ScalarQ, PD] :> 
			Inner[
				covd[{#1, -BasisOfCovD@ covd}] @ expr1 Diff @ #2&, 
				CNumbersOf[basis, VBundleOfBasis[basis]], 
				ScalarsOfChart[basis], Plus
			] /; (Deg @ expr1 === 0 && basis =!= Null),
			
		(* Exterior derivative of the volume element associated to a metric *)
		If[metric =!= Null,
			With[{eps = epsilon[metric]},
				HoldPattern[Diff[eps[inds__?DownIndexQ], PD]] :> 
					Module[{a = DummyIn[VBundleOfMetric[metric]]},
						ChristoffelForm[covd][a, -a] eps[inds]
					]
			],
			{}
		]
	}
]; 



(******************* 3.6 Interior product and Lie derivative *****************************)

(* We define the interior product of form with a vector as a derivation of degree -1. *)


DefGradedDerivation[Int[v_], Wedge, -1, PrintAs->"\[Iota]"];


(* Main properties of the interior product: *)
Int[v_][f_ form_] := f Int[v][form] /; Deg @ f === 0; 

Int[v_][f_] := 0 /; Deg @ f === 0; 

Int[f_?ScalarQ v_][form_] := f Int[v][form]; 

Int[v_[ind1_Symbol]][Coframe[mani_][ind2_]] := v[ind2] /; v =!= Coframe[mani]; 

Int[v_[ind1_Symbol]][dx[mani_][ind2_]] := v[ind2] /; v =!= dx[mani]; 

Int[v_[-ind1_Symbol]][Coframe[mani_][ind2_]] := v[-ind1] (First @ MetricsOfVBundle@ Tangent @ mani)[ind1, ind2]; 

Int[v_[-ind1_Symbol]][dx[mani_][ind2_]] := v[-ind1] (First @ MetricsOfVBundle@ Tangent @ mani)[ind1, ind2]; 

Int[Basis[{cnumber1_?IntegerQ, -basisname_?BasisQ}, ind_Symbol]][dx[mani_][{cnumber2_?IntegerQ, basisname_?BasisQ}]] := 0 /; cnumber1 =!= cnumber2; 

Int[Basis[{cnumber1_?IntegerQ, -basisname_?BasisQ}, ind_Symbol]][dx[mani_][{cnumber2_?IntegerQ, basisname_?BasisQ}]] := 1 /; cnumber1 === cnumber2; 

Int[Basis[{cnumber1_?IntegerQ, -basisname_?BasisQ}, ind_Symbol]][Coframe[mani_][{cnumber2_?IntegerQ, basisname_?BasisQ}]] := 0 /; cnumber1 =!= cnumber2; 

Int[Basis[{cnumber1_?IntegerQ, -basisname_?BasisQ}, ind_Symbol]][Coframe[mani_][{cnumber2_?IntegerQ, basisname_?BasisQ}]] := 1 /; cnumber1 === cnumber2; 


(* Interior contraction of basis objects is zero. *)
Int[v_][_Basis] := 0; 


(* Cartan operator acting on forms *)
DefGradedDerivation[CartanD[v_], Wedge, 0, PrintAs -> "L"]; 


(* e acting on basis objects is zero. *)
CartanD[v_][_Basis] := 0; 


CartanD[v_][_Basis, _] := 0; 


(* Superscript formatting for Lie derivatives arising from the exterior derivatives *)
CartanD /: MakeBoxes[CartanD[v_][form_, PD?CovDQ], StandardForm] := 
	xAct`xTensor`Private`interpretbox[
		CartanD[v][form, PD], 
		RowBox[
			{
				SubscriptBox["L", MakeBoxes[v, StandardForm]], 
				"[", 
				MakeBoxes[form, StandardForm], "]"
			}
		]
	]; 

CartanD /: MakeBoxes[CartanD[v_][form_, cd_?CovDQ], StandardForm] := 
	xAct`xTensor`Private`interpretbox[
		CartanD[v][form, cd], 
		RowBox[
			{
				SubsuperscriptBox["L", MakeBoxes[v, StandardForm], Last @ SymbolOfCovD[cd]], 
				"[", 
				MakeBoxes[form, StandardForm], 
				"]"
			}
		]
	]; 


CartanD[f_?ScalarQ v_][form_] := f CartanD[v] @ form + Wedge[Diff @ f, Int[v] @ form]; 

(* We still need definition when acting on Times *)
(* This produces expanded expressions and is much faster when there are many scalars *)

CartanD[v_][expr_Times] := 
	Module[{grades = Grade[#, Wedge]& /@ List @@ expr, pos, scalar, form},
		pos = Position[grades, _?(# =!= 0&), 1, Heads -> False]; 
		Which[
			Length[pos] > 1,
				Throw[Message[Diff::error1, "Found Times product of nonscalar forms: ", expr]],
			Length[pos] === 1,
				pos = pos[[1, 1]]; 
				scalar = Delete[expr, {pos}]; 
				form = expr[[pos]]; 
				scalar CartanD[v][form] + lie0[v][scalar, form],
			Length[pos] === 0,
				lie0[v][expr]
		]
	]; 

(* Only scalars *)

lie0[v_][expr_Times] := Sum[MapAt[CartanD[v], expr, i], {i, 1, Length[expr]}]; 

lie0[v_][expr_] := CartanD[v][expr]; 

(* Scalars and a form *)

lie0[v_][expr_Times, form_] := Sum[MapAt[lie0[v][#, form]&, expr, i], {i, 1, Length[expr]}]; 

lie0[v_][expr_, form_] := Wedge[CartanD[v][expr], form]; 


CartanD[v_][expr_?ConstantQ] := 0; 


(* Leibinitz rule for xTensor LieD when acting on wedge product expressions: *)
Unprotect @ LieD; 

LieD[v_] @ expr_Wedge := (CartanD[v] @ expr /. CartanD -> LieD); 

Protect @ LieD; 


(* Lie and diff commute (probably not true for covariant exterior derivatives). We need to pick a canonical order. Should this be done by a SortDs type function? *)
(* Cartan identity: *)
CartanDToDiff[expr_] := expr /. {
	CartanD[v_][form_] :> Diff @ Int[v] @ form + Int[v] @ Diff @ form, 
	CartanD[v_][form_, covd_?CovDQ] :> Diff[Int[v] @ form, covd] + Int[v] @ Diff[form, covd]
	}; 


(* Lie derivative of a tensor valued differential form. We overload the xTensor command LieDToCovD. Note: xTensor LieD doesn't know how to work with wedge products. *)
LieDForm[v_[ind_], covd_?CovDQ] @ ten_ := 
	Module[{a = DummyIn @ VBundleOfIndex@ ind},
		ToCanonical[LieD[v[ind], covd] @ ten + CartanD[v[ind]][ten, covd] - v[a] covd[-a] @ ten, UseMetricOnVBundle -> None]
	]; 

LieDFormToCovD[expr_, covd_] := expr //. LieD[vector_] :> LieDForm[vector, covd]; 


Unprotect @ LieDToCovD; 
	LieDToCovD[expr_, arg_:PD] := LieDFormToCovD[expr, arg] /; Deg @ expr >= 1
Protect @ LieDToCovD; 

(***************** 3.7 Commuting and sorting of the derivations ********************)

(*

There are 6 equations in all for the (super-)commutators of the three derivations (also the negations of these three equations). 
It is easy to write them down. The LHS's are just super-commutators of two derivations, in whatever order you want. 
The sign is + just for [odd,odd] so that's for [d,d], [\[Iota],\[Iota]], and [\[Iota],d]. 
Then the RHS is itself a single derivation, with the degree being the sum of degrees from the LHS. 
If there is only one vector argument, it's clear what the RHS should be. 
When there are two vector arguments, 
then the vector on the RHS is the Lie bracket of the two from the LHS, in the order that they came on the LHS.
 
The equations are:
1. d^2=0
2. Subscript[\[Iota], v] Subscript[\[Iota], w]+Subscript[\[Iota], w] Subscript[\[Iota], v]=0
3. Subscript[\[ScriptCapitalL], v] Subscript[\[ScriptCapitalL], w]-Subscript[\[ScriptCapitalL], w] Subscript[\[ScriptCapitalL], v]=Subscript[\[ScriptCapitalL], [v,w]]
4. Subscript[\[Iota], v] Subscript[\[ScriptCapitalL], w]-Subscript[\[ScriptCapitalL], w] Subscript[\[Iota], v] = Subscript[\[Iota], [v,w]]
5. Subscript[\[Iota], v]d+d Subscript[\[Iota], v]= Subscript[\[ScriptCapitalL], v]
6. Subscript[\[ScriptCapitalL], v]d-d Subscript[\[ScriptCapitalL], v]=0 ("Cartan's magic formula")

*)

(*

1. requires no action (it will just be automatically killed).
2, 3. Maybe have this sorted automatically? Any canonical order for v,w imposes a canonical order on the derivations.
4, 5, 6\[LongDash]write rules to go in either direction.

*)

SortDerivationsRule[Diff, Diff] = {}; 


SortDerivationsRule[Int, Int] = {HoldPattern[Int[v_] @ Int[w_] @ form_] :> -Int[w] @ Int[v] @ form /; !OrderedQ[{v, w}]}; 


SortDerivationsRule[CartanD, CartanD] = {
	HoldPattern[CartanD[v_] @ CartanD[w_] @ form_] :> 
		Module[{a = First @ FindFreeIndices[v]},
			CartanD[w] @ CartanD[v] @ form + CartanD[Bracket[v, w][a]] @ form
		] /; !OrderedQ[{v, w}],
		
	HoldPattern[CartanD[v_][CartanD[w_][form_, covd_?CovDQ], covd_?CovDQ]] :> 
		Module[{a = First @ FindFreeIndices[v]},
			CartanD[w][CartanD[v][form, covd], covd] + CartanD[Bracket[v, w][a]][form, covd]
		] /; !OrderedQ[{v, w}]
}; 


(* Below, SortDerivationsRule[d1,d2][expr] changes instances of d2\[CenterDot]d1 into d1\[CenterDot]d2. There is no need to check orderings. *)
SortDerivationsRule[Int, CartanD] = {
	HoldPattern[CartanD[w_] @ Int[v_] @ form_] :> 
		Module[{a = First @ FindFreeIndices[w]},
			Int[v] @ CartanD[w] @ form + Int[Bracket[w, v][a]] @ form
		],
	HoldPattern[CartanD[w_][Int[v_] @ form_, covd_?CovDQ]] :> 
		Module[{a = First @ FindFreeIndices[w]},
			Int[v] @ CartanD[w][form, covd] + Int[Bracket[w, v][a]] @ form
		]
}; 

SortDerivationsRule[CartanD, Int] = {
	HoldPattern[Int[v_] @ CartanD[w_] @ form_] :> 
		Module[{a = First @ FindFreeIndices[w]},
			CartanD[w] @ Int[v] @ form + Int[Bracket[v, w][a]] @ form
		],
	HoldPattern[Int[v_] @ CartanD[w_][form_, covd_?CovDQ]] :> 
		Module[{a = First @ FindFreeIndices[w]},
			CartanD[w][Int[v] @ form, covd] + Int[Bracket[v, w][a]] @ form
		]
}; 


(* ::Input::Initialization:: *)
SortDerivationsRule[Int, Diff] = {HoldPattern[Diff[Int[v_] @ form_, covd_?CovDQ]] :> -Int[v] @ Diff[form, covd] + CartanD[v][form, covd]}; 

SortDerivationsRule[Diff, Int] = {HoldPattern[Int[v_] @ Diff[form_, covd_?CovDQ]] :> -Diff[Int[v] @ form, covd] + CartanD[v][form, covd]}; 


(* ::Input::Initialization:: *)
SortDerivationsRule[CartanD, Diff] = {
	HoldPattern[Diff[CartanD[v_] @ form_, PD]] :> CartanD[v] @ Diff @ form, 
	HoldPattern[Diff[CartanD[v_][form_, covd_?CovDQ], covd_?CovDQ]] :> CartanD[v][Diff[form, covd], covd]
}; 

SortDerivationsRule[Diff, CartanD] = {
	HoldPattern[CartanD[v_] @ Diff[form_, PD]] :> Diff @ CartanD[v] @ form, 
	HoldPattern[CartanD[v_][Diff[form_, covd_?CovDQ], covd_?CovDQ]] :> Diff[CartanD[v][form, covd], covd]
}; 


(* Now we choose a default left-to-right order that we want derivations to have. I don't know if there is a best order (i.e. what encourages the largest number of terms to vanish). *)

$Derivations = {CartanD, Int, Diff}; 

$DerivationSortOrder = $Derivations; 


(* Code of SortDerivations *)
SortDerivations[expr_] := SortDerivations[expr, $DerivationSortOrder]

SortDerivations[expr_, order_List] := 
	Module[
		{},
	(* Make sure that order is some permutation of $Derivations *)
		If[Sort @ order =!= Sort @ $Derivations,
			Throw @ Message[SortDerivations::invalid, "order", order]; 
		]; 
		expr //. Join @@ Table[Join @@ (SortDerivationsRule[order[[i]], #]& /@ Drop[order, i]), {i, 1, Length @ order - 1}] 
		//. Join @@ (SortDerivationsRule[#, #]& /@ order)
	]; 

(*********************************************************************************************)
(************************** 4. Variational derivatives of top rank forms *********************)
(*********************************************************************************************)

(************************** 4.1 FormVarD *************************************)



(* Check that a form has top rank *)

TopRankQ[form_] := With[{manifolds = Select[DependenciesOf @ form, ManifoldQ]},
	If[Length @ manifolds != 1,
		Throw @ Message[TopRankQ::error1, "Forms must have exactly 1 manifold in dependencies."],
		TopRankQ[form, First @ manifolds]
	]
]; 

TopRankQ[form_, mani_?ManifoldQ] := Grade[form, Wedge] === DimOfManifold @ mani; 


(* 

FormVarD is supposed to be like VarD, but for top-rank forms. A derivative is not specified, because 
the exterior derivative is natural. However, the adjoint to the exterior derivative, the codifferential, 
requires a Hodge dual, which we define via a (non-degenerate!) metric. Therefore FormVarD has the notation FormVarD[form, metric][expr, rest]. 
The convention is that expr and rest are combined as \[Delta](expr)\[Wedge]rest (here I write the variational derivative as \[Delta], and always subscript the codifferential 
as Subscript[\[Delta], g]). This means that when writing the Leibniz rule for the Wedge product, we must re-order the factors, i.e.

  \[Delta](a\[Wedge]b\[Wedge]c) = \[Delta]a\[Wedge]b\[Wedge]c + a\[Wedge]\[Delta]b\[Wedge]c + a\[Wedge]b\[Wedge]\[Delta]c=
  \[Delta]a\[Wedge]b\[Wedge]c + (-1)^(|a||b|)\[Delta]b\[Wedge]a\[Wedge]c + (-1)^((|a|+|b|)|c|)\[Delta]c\[Wedge]a\[Wedge]b 
  
  etc.
The other identities we need follow. For a Hodge star,
 
  a\[Wedge](\!\(
\*SubscriptBox[\(\[Star]\), \(g\)]\(b\)\)) 
  = Subscript[\[LeftAngleBracket]a,b\[RightAngleBracket], g] Subscript[dVol, g]     

where |a|+|*b|=n, and where

Subscript[\[LeftAngleBracket]a,b\[RightAngleBracket], g] = g^(Subscript[i, 1] Subscript[j, 1])... g^(Subscript[i, k] Subscript[j, k]) 
Subscript[a, Subscript[i, 1]... Subscript[i, k]] Subscript[b, Subscript[j, 1]... Subscript[j, k]]

is the inner product on k-forms given by the metric g, and Subscript[dVol, g]=\[Star]1 is the metric volume element, given in a coordinate basis by

Subscript[dVol, g]=Sqrt[|det g|](dx^1)\[Wedge]...\[Wedge]dx^n.

With a real metric and real-valued forms, it is clear that \[LeftAngleBracket]a,b\[RightAngleBracket] = \[LeftAngleBracket]b,a\[RightAngleBracket] and so
  a\[Wedge](\!\(
\*SubscriptBox[\(\[Star]\), \(g\)]\(b\)\)) = 
  b\[Wedge](\!\(
\*SubscriptBox[\(\[Star]\), \(g\)]\(a\)\))        (for |a|+|*b|=n) 

or in the way we'll use it,
  (\!\(
\*SubscriptBox[\(\[Star]\), \(g\)]\(a\)\))\[Wedge]b = (-1)^(|a|(n-|a|)) a\[Wedge](\!\(
\*SubscriptBox[\(\[Star]\), \(g\)]\(b\)\))      (for |a|+|*b|=n) 
QUESTION: What about \[DoubleStruckCapitalC]-valued forms?
We also have that diff and codiff (of any CovD) are adjoints of each other, in the sense that

  \[Integral]  \!\(\(\(da\)\(\[Wedge]\)\)
\*SubscriptBox[\(\[Star]\), \(g\)]b\) = \[Integral]  \!\(\(\(b\)\(\[Wedge]\)\)
\*SubscriptBox[\(\[Star]\), \(g\)]da\) = \[Integral]  \!\(\(\(a\)\(\[Wedge]\)\)
\*SubscriptBox[\(\[Star]\), \(g\)]\(
\*SubscriptBox[\(\[Delta]\), \(g\)]b\)\) = \[Integral]  \!\(\(
\*SubscriptBox[\(\[Delta]\), \(g\)]\(\(b\)\(\[Wedge]\)\)\)
\*SubscriptBox[\(\[Star]\), \(g\)]\(a\n\n
The\ above\ rule\ is\ convenient\ with\ the\ inverse\ of\ Hodge\)\), \!\(\(i . e . \)\ 
\*SubsuperscriptBox[\(\[Star]\), \(g\), \(-1\)]
\*SubscriptBox[\(\[Star]\), \(g\)]\)=\!\(
\*SubscriptBox[\(\[Star]\), \(g\)]
\*SubsuperscriptBox[\(\[Star]\), \(g\), \(-1\)]\)=id. This inverse is
  Subsuperscript[\[Star], g, -1]a=\!\(\(
\*SuperscriptBox[\((\(-1\))\), \(\(|\)\(a\)\(|\)\((\(\(n\)\(-\)\) | a | )\)\)]s\)\ 
\*SubscriptBox[\(\[Star]\), \(g\)]\(a\n
where\ s\)\)=SignDetOfMetric[g]

TODO:
  1. Act on Subscript[\[Star], g2] and Subscript[\[Delta], g2] where g2 is a second metric on the same manifold, and is also non-degenerate. Subscript[\[Star], g2] can be converted to Subscript[\[Star], g] by including a ratio of volume elements (the ratio is coordinate-free).
  2. What to do with Lie, Int?
  3. What about coframe, connection, torsion, curvature, etc.?


*)

InvHodge[metric_][expr_] := With[{k = Grade[expr, Wedge], n = DimOfMetric@ metric, s = SignDetOfMetric @ metric},
	(-1)^(k (n - k)) s Hodge[metric] @ expr
]; 


(* 

We will use it as follows:

\[Integral] da \[Wedge] c = \[Integral] a \[Wedge] 
\!\(
\*SubscriptBox[\(\[Star]\), \(g\)]
\*SubscriptBox[\(\[Delta]\), \(g\)]
\*SubsuperscriptBox[\(\[Star]\), \(g\), \(-1\)]\(c\)\)

and

  \[Integral] Subscript[\[Delta], g]b \[Wedge] c = \[Integral] \!\(\(\(b\)\(\ \)\(\[Wedge]\)\)\ 
\*SubscriptBox[\(\[Star]\), \(g\)]d
\*SubsuperscriptBox[\(\[Star]\), \(g\), \(-1\)]c\)

*)
(* TODO:More checks that form is actually on same manifold as metric, etc. *)

(* Generate rest. Replace dummies in expr. This does not act on scalar arguments of functions *)

FormVarD[form_, met_][expr_] := 
If[TopRankQ[expr] && ScalarQ[expr],
	FormVarD[form, met][ReplaceDummies @ expr, 1],
	Throw @ Message[General::error1, "Can only take variational derivative of top-rank form."]
]; 

(* Thread over Plus *)

FormVarD[form_, met_][expr_Plus, rest_] := FormVarD[form, met][#, rest]& /@ expr; 

FormVarD[form_, met_][expr_SeriesData, rest_] := xAct`xTensor`Private`SeriesDataMap[FormVarD[form, met][#, rest]&, expr]; 

(* FormVarD on products: sum of VarDtake's of elements *)

FormVarD[form_, met_][expr_Times, rest_] := With[{grades = Grade[#, Wedge]& /@ List @@ expr},
	If[Length @ Position[grades, _?(# =!= 0&), 1, Heads -> False] > 1,
		Throw[Message[FormVarD::error1, "Found Times product of nonscalar forms: ", expr]]
	]; 
Sum[FormVarDtake[form, met, rest, List @@ expr, count], {count, Length @ expr}]
]; 

(* FormVarD on wedges: sum of FormVarDtake's of elements.
Note the use of sumgrades for reordering the Wedge, as described above. *)

FormVarD[form_, met_][expr_Wedge, rest_] := With[{grades = Grade[#, Wedge]& /@ List @@ expr},
	With[{sumgrades = FoldList[Plus, 0, grades]},
		Sum[(-1)^(grades[[count]] * sumgrades[[count]]) FormVarDtake[form, met, rest, List @@ expr, count], {count, Length @ expr}]
	]
]; 

(* FormVarD element n of a list of Wedge factors (no sign--it was included above) *)

FormVarDtake[form_, met_, rest_, list_List, n_Integer] := FormVarD[form, met][list[[n]], Wedge[Wedge @@ Drop[list, {n}], rest]]; 

           (* Scalar functions. Multiargument generalization contributed by Leo. multiD is not enough here.
 Since operating on a scalar function, don't need extra Wedge's *)

FormVarD[form_, met_][func_?ScalarFunctionQ[args__], rest_] := 
With[{repargs = ReplaceDummies /@ {args}},
	Plus @@ MapThread[FormVarD[form, met][#1, rest (Derivative @@ #2) @ func @@ repargs]&, {repargs, IdentityMatrix @ Length @ repargs}]
]; 

(* Remove Scalar head because in general the result is not a scalar *)

FormVarD[form_, met_][Scalar[expr_], rest_] := FormVarD[form, met][ReplaceDummies[expr], rest]; 

(* Constants *)

FormVarD[_, _][x_?ConstantQ, _] := 0; 

                  (* Same tensor: metric. Do not use ContractMetric, which hides the metric.
Note: This part is identical to the code in VarD, since it's only metric being Wedged with rest. *)

FormVarD[metric_[a_, b_], met_][metric_Symbol?MetricQ[c_, d_], rest_] := 
	xAct`xTensor`Private`metricsign[a, b, c, d] ToCanonical[
		rest (metric[ChangeIndex @ a, c] 
		metric[ChangeIndex@ b, d] + metric[ChangeIndex @ a, d] metric[ChangeIndex @ b, c]) / 2,
		UseMetricOnVBundle -> None
	]; 

                  (* Same tensor. Place indices in proper delta positions. QUESTION: could this be problematic for spinors?
Note: This part is identical to the code in VarD, since it's only deltas being Wedged with rest. *)

FormVarD[form_[inds1___], met_][form_?xTensorQ[inds2___], rest_] := 
	With[{clist = ChangeIndex /@ IndexList[inds1]},
		ToCanonical[
			ImposeSymmetry[
				Inner[xAct`xTensor`Private`varddelta, clist, IndexList[inds2], Times], 
				clist, 
				SymmetryGroupOfTensor[form[inds1]]
			]rest, 
			UseMetricOnVBundle -> None
		]
]; 

(* A different tensor *)

FormVarD[form1_[inds1___], met_][form2_?xTensorQ[inds2___], rest_] := 0 /;
	 !ImplicitTensorDepQ[form2, form1]; 

(* Hodge identity *)

FormVarD[form_, met_][Hodge[met_][expr_], rest_] := 
	With[{k = Grade[expr, Wedge], n = DimOfMetric @ met},
		(-1)^(k (n - k)) FormVarD[form, met][expr, Hodge[met] @ rest]
	]; 

(* diff \[Rule] Replaced by Diff to adjust to the new notation. Dropped cd. Added back PD. *)

FormVarD[form_, met_][Diff[expr_, PD], rest_] := 
	-FormVarD[form, met][expr, Hodge[met] @ Codiff[met][InvHodge[met] @ rest]]; 

(* codiff \[Rule] Replaced by Codiff to adjust to the new notation. Dropped cd and replaced ExtCovDiff by Diff . Added back covd *)

FormVarD[form_, met_][Codiff[met_][expr_, covd_?CovDQ], rest_] := 
	-FormVarD[form, met][expr, Hodge[met] @ Diff[InvHodge[met] @ rest, covd]]; 


(****************************************************************)
(****************** 5. End private and package ******************)
(****************************************************************)

End[]; 

EndPackage[]; 
