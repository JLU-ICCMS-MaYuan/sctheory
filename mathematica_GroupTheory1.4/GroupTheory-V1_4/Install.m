(****m* /Install.m
!
! NAME
!  Install.m
! AUTHOR
!  W. Hergert, M. Geilhufe, S. Schenk
! MODIFICATION HISTORY
!  2/23/12 : initial documentation 
! USAGE
!  see single modules
!
! ERROR MESSAGES
!  see single modules
! GTPack MODULES
!
! --- Installation of symmetry elements ---
!
!  GTInstallAxis installs symmetry elements according to a certain axis denoted by a label and defined by cartesian coordinates
!  GTReinstallAxes reinstalls all symmetry elements. According to the chosen convention, symmetry elements are either defined via active or passive rotations
!  GTWhichAxes informs about the actual setting for symmmetry transformations (active or passive)
!
! --- Installation of point groups ---
!
!  GTInstallGroup gives a faithful representation of a crystallographic point group
!  GTTableToGroup gives a faithful representation of an arbitrary group from its multiplication table
!  GTGroupFromGenerators installs a finite group from a list of generators
! 
! --- Standard representation ---
!
!  GTChangeRepresentation changes the used standard representation
!  GTWhichRepresentation prints the used standard representation
!  
!--------------------------------------------------------------------------------
!
***)
BeginPackage["GroupTheory`Install`",{"GroupTheory`Symbols`","GroupTheory`Basic`","GroupTheory`Auxiliary`","GroupTheory`CrystalStructure`"}]

(*--------------------------- Installation of symmetry elements ----------------------------*)
 GTInstallAxis          ::usage = "GTInstallAxis[\*StyleBox[\"{label, coordinates}\", \"TI\"]] installs symmetry elements with respect to a given rotation axis denoted by \*StyleBox[\"label\", \"TI\"] and defined by a \*StyleBox[\"direction\", \"TI\"] in Cartesian coordinates."
 GTReinstallAxes 		::usage = "GTReinstallAxes[\*StyleBox[\"convention\", \"TI\"]] reinstalls all symmetry elements. According to the chosen \*StyleBox[\"convention\", \"TI\"], symmetry elements are either defined via active or passive rotations."
 GTWhichAxes            ::usage = "GTWhichAxes[] informs about the actual setting for symmmetry transformations (active or passive)"
 GTIcosahedronAxes      ::usage=  "GTIcosahedronAxes[] installs axes for icosahedron group."

(*--------------------------- Installation of point groups ----------------------------*)
 GTGroupFromGenerators  ::usage = "GTGroupFromGenerators[\*StyleBox[\"generators\", \"TI\"]] installs a finite group from a list of \*StyleBox[\"generators\", \"TI\"]."
 GTInstallGroup         ::usage = "GTInstallGroup[\*StyleBox[\"group\",\"TI\"]] gives a faithful representation of a crystallographic group."
 GTTableToGroup         ::usage = "GTTableToGroup[\*StyleBox[\"list of elements,multiplication table\", \"TI\"]] gives a faithful representation of an arbitrary group from a given \*StyleBox[\"list of elements\", \"TI\"] and a \*StyleBox[\"multiplication table\", \"TI\"], using permutation matrices."
 
(*--------------------------- Standard representation ---------------------------------*)
 GTChangeRepresentation ::usage = "GTChangeRepresentation[\*StyleBox[\"representation\", \"TI\"]] changes the currently used standard \*StyleBox[\"representation\", \"TI\"]."
 GTWhichRepresentation  ::usage = "GTWhichRepresentation gives the currently used standard representation."
  
(*
 GTSetIndex             ::usage = "Switching between matrix representation of group elements. Set 1 for O(3), 2 for SU(2), 3 for SO(2) and 4 for permutation matrices." 
 GTWhichIndex           ::usage = "Which matrix representation is used? Gives 1 for O(3), 2 for SU(2), 3 for SO(2) and 4 for permutation matrices."
*)

(*--------------------------- Options ----------------------------*)
 GTInstallGroup = GroupTheory`Auxiliary`Private`GTInstallGroup
 
Options[GTChangeRepresentation]    ={GOVerbose->True}
Options[GTInstallAxis]             ={GOVerbose ->True}
Options[GTInstallGroup]            ={GORepresentation->"O(3)",GOVerbose->True}
Options[GTTableToGroup]            ={GOVerbose->True}
Options[GTWhichRepresentation]     ={GOVerbose->True}
Options[GTIcosahedronAxes]         = {GOVerbose -> True}

GTInstallGroup = GroupTheory`Auxiliary`Private`GTInstallGroup
GTSymbolInfoToSU2Matrix = GroupTheory`Basic`Private`GTSymbolInfoToSU2Matrix
GTSymbolInfoToSymbol = GroupTheory`Basic`Private`GTSymbolInfoToSymbol
GTSymbolInfoToO3Matrix = GroupTheory`Basic`Private`GTSymbolInfoToO3Matrix
GTSymbolInfoToEulerAngles = GroupTheory`Basic`Private`GTSymbolInfoToEulerAngles
 
 

Begin["`Private`"] 

(* Begin Private Context *)
(*--------------------------- external Modules -------------------*)
(* only non public functions must have this *)

(*--------------------------- Modules ----------------------------*)



(****g* /GTSIcosahedronAxes
! NAME
!  GTIcosahedronAxes
! AUTHOR
!  M. Geilhufe, W. Hergert
! PACKAGE
!  Install.m
! MODIFICATION HISTORY
!   * 06.08.2016 : first version
!   * 10.10.2018 : implementation in Install.m and check of header
! USAGE
!  GTIcosahedronAxes[] installs axes for icosahedron group.
! INPUT
!  no input argument
! OUTPUT
!  no direct output, but all the axes are installed an can be used symbollically.
! GTPack OPTIONS
!  o GOVerbose 
!     - False : Messages are switched of
!     - True  : Complains due to precision, but this does not matter
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTInstallAxis
! GTPack NOTEBOOKS 
!  GTIcosaeder_new.nb in Wolfram_Devel/4_Discrete_Symmetry_Groups/Nanotubes_Buckyballs.
! DESCRIPTION
! 
! LITERATURE
!  Altman, Herzig
! TODO
! 
! PROBLEMS
!  ! Two axes are only defined and not installed !! (what does this mean? 10.10.18 W.)
!--------------------------------------------------------------------------------
! SOURCE
*)


GTIcosahedronAxes[OptionsPattern[]] := Module[
	{verb, cn5, sn10, \[Sigma], \[Rho], \[Rho]s, ne1, ne2, ne3, ne4, ne5, ne6, ne7, ne8, ne9, ne10, ne11, ne12, neg,
     neh, nei, nej, nek, nel, nem, nen, neo, nep, neq, ner, ne\[Epsilon], ne\[Phi], ne\[Eta], ne\[CurlyPhi], list1, list2, 
     list3, i},
     verb = OptionValue[GOVerbose];
     (*--- parameter definition ---*)
     cn5      = Cos[\[Pi]/5];
     sn10     = Sin[\[Pi]/10];
     \[Sigma] = ArcTan[2, 1 + Sqrt[5]];
     \[Rho]   = ArcTan[2, 3 - Sqrt[5]];
     \[Rho]s  = ArcTan[2, 3 + Sqrt[5]];
    (*--- first set of axes ---*)
     ne1  = {"1", {Cos[\[Sigma]], 0, Sin[\[Sigma]]}};
     ne2  = {"2", {-Cos[\[Sigma]], 0, Sin[\[Sigma]]}};
     ne3  = {"3", {0, -Sin[\[Sigma]], Cos[\[Sigma]]}};
     ne4  = {"4", {0, Sin[\[Sigma]], Cos[\[Sigma]]}};
     ne5  = {"5", {Sin[\[Rho]], Cos[\[Rho]], 0}};
     ne6  = {"6", {-Sin[\[Rho]], Cos[\[Rho]], 0}};
     ne7  = {"7", {0, Sin[\[Rho]], Cos[\[Rho]]}};
     ne8  = {"8", {0, -Sin[\[Rho]], Cos[\[Rho]]}};
     ne9  = {"9", {-Sin[\[Sigma]], Cos[\[Sigma]], 0}};
     ne10 = {"\[Kappa]", {Sin[\[Sigma]], Cos[\[Sigma]], 0}};
     ne11 = {"\[Lambda]", {-Cos[\[Rho]], 0, Sin[\[Rho]]}};
     ne12 = {"\[Mu]", {Cos[\[Rho]], 0, Sin[\[Rho]]}};
    (*--- second set of axes ---*)  
     ne\[Epsilon] = {"\[Epsilon]", {1, 1, 1}/Sqrt[3]};
     ne\[Phi] = {"\[Phi]", {-1, -1, 1}/Sqrt[3]};
     ne\[CurlyPhi] = {"\[CurlyPhi]", -{1, -1, -1}/Sqrt[3]};
     ne\[Eta] = {"\[Eta]", -{-1, 1, -1}/Sqrt[3]};
    (*--- third set of axes ---*)
     neg = {"g", {sn10, -1/2, cn5}};
     neh = {"s", {sn10, 1/2, cn5}};
     nei = {"t", {-sn10, -1/2, cn5}};
     nej = {"j", {-sn10, 1/2, cn5}};
     nek = {"k", {-cn5, -sn10, 1/2}};
     nel = {"l", {-cn5, sn10, 1/2}};
     nem = {"m", {cn5, -sn10, 1/2}};
     nen = {"n", {cn5, sn10, 1/2}};
     neo = {"o", {-1/2, -cn5, sn10}};
     nep = {"p", {-1/2, cn5, sn10}};
     neq = {"q", {1/2, -cn5, sn10}};
     ner = {"r", {1/2, cn5, sn10}};
   (*--- lists with all definitions ---*)
     list1 = {ne1, ne2, ne3, ne4, ne5, ne6, ne7, ne8, ne9, ne10, ne11, ne12} // FullSimplify // ToRadicals;
     list2 = {neg, neh, nei, nej, nek, nel, nem, nen, neo, nep, neq, ner} // FullSimplify // ToRadicals;
     list3 = {ne\[Epsilon], ne\[Phi], ne\[CurlyPhi], ne\[Eta]} // FullSimplify // ToRadicals;
     If[verb,
       Print["26 new axes are defined and ready for installation. Needs some time!"],
       None
     ];
   (*--- Installation of all axes ---*)
     If[verb,
        Print["Warnings with respect to precision do not matter!"],
        Off[N::meprec]
     ];
     Do[
       GTInstallAxis[list1[[i, 1]], list1[[i, 2]], GOVerbose -> verb]
     , {i, 1, Length[list1]}];
     Do[
     	GTInstallAxis[list2[[i, 1]], list2[[i, 2]], GOVerbose -> verb]
     , {i, 1, Length[list2]}];
     Do[
     	GTInstallAxis[list3[[i, 1]], list3[[i, 2]], GOVerbose -> verb]
     ,{i, 3, Length[list3]}];
     If[verb,
        None,
        On[N::meprec]
     ];
  ]




(*
***)




(****g* /GTWhichRepresentation
!
! NAME
!  GTWhichRepresentation.m
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Install.m
! MODIFICATION HISTORY
!  2/23/12 : initial documentation 
! USAGE
!  GTWhichRepresentation prints the used standard representation.
! EXAMPLE
! 
! INPUT
!   
! OUTPUT
!  
! ERROR MESSAGES
!  none
! GTPack MODULES
!  none
! GTPack NOTEBOOKS 
!  none
! DESCRIPTION
! 
! LITERATURE
!
! TODO
! 
! PROBLEMS
!  
! Release
!
!--------------------------------------------------------------------------------
! SOURCE

*)

GTWhichRepresentation[OptionsPattern[]] := Module[{rep}, 
	rep = Which[
		grpdgrp==1,	If[OptionValue[GOVerbose],Print["The standard representation is O(3)"]];
					"O(3)",
		grpdgrp==2,	If[OptionValue[GOVerbose],Print["The standard representation is SU(2)"]];
					"SU(2)",
		grpdgrp==3,	If[OptionValue[GOVerbose],Print["The standard representation is O(2)"]];
					"O(2)",
		grpdgrp==4,	If[OptionValue[GOVerbose],Print["The standard representation is given by permutation matrices."]];
					"own",
		grpdgrp==5,	If[OptionValue[GOVerbose],Print["The standard representation is SU(2)xS."]];
					"SU(2)xS"];
		Return[rep];
]

(*
***)


(****g* /GTChangeRepresentation
!
! NAME
!  GTChangeRepresentation.m
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Install.m
! MODIFICATION HISTORY
!  2/23/12 : initial documentation 
!  September 2014, SU(2)xS
! USAGE
!  GTChangeRepresentation[representation] changes the used standard representation.
! EXAMPLE
!  GTChangeRepresentation[Rep_]
! INPUT
!  Rep="O(3)" -> O(3) , Rep="SU(2)" -> SU(2) , Rep="O(2)" -> O(2), "Own" -> own definition, Rep="SU(2)xS" -> SU(2)xS 
! OUTPUT
!  
! ERROR MESSAGES
!  none
! GTPack MODULES
!  none
! GTPack NOTEBOOKS 
!  none
! DESCRIPTION
! 
! LITERATURE
!
! TODO
! 
! PROBLEMS
! 
! Release
! 
!--------------------------------------------------------------------------------
! SOURCE

*)

GTChangeRepresentation[Rep_,OptionsPattern[]] := Module[{sonst}, 
	Unprotect[grpdgrp];
    Clear[grpdgrp];
    sonst=True;
    Which[	Rep=="O(3)",
    			grpdgrp=1;
    			sonst=False;
    			If[OptionValue[GOVerbose],Print["The standard representation has changed to O(3)"],None],
	 	  	Rep=="SU(2)",
	 	  		grpdgrp=2;
	 	  		sonst=False;
	 	  		If[OptionValue[GOVerbose],Print["The standard representation has changed to SU(2)"],None],
	 		Rep=="O(2)",
	 			grpdgrp=3;
	 			sonst=False;
	 			If[OptionValue[GOVerbose],Print["The standard representation has changed to O(2)"],None],
			Rep=="Permutation",
				grpdgrp=4;
				sonst=False;
				If[OptionValue[GOVerbose],Print["The standard representation has changed to permutation matrices."],None],
			Rep=="SU(2)xS",
				grpdgrp=5;
				sonst=False;
				If[OptionValue[GOVerbose],Print["The standard representation has changed to SU(2)xS."],None]
    ];
    If[sonst,Print[$ErrorChangeRepresentation];Abort[],None];
    Protect[grpdgrp];
]
(*
***)

(****g* /GTChangeRepInst
!
! NAME
!  GTChangeRepIns
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Install.m
! MODIFICATION HISTORY
!  5/22/14 : initial documentation 
! USAGE
!  Switching between matrix representations. It is possible to use O(3), SU(2), O(2). Special case of GTChangeRepresentation, only used by GTInstallGroup
! EXAMPLE
!  GTChangeRepInst[Rep_]
! INPUT
!  Rep=1 -> O(3) , Rep="U(2)" -> U(2) , Rep="O(2)" -> O(2)
! OUTPUT
!  
! ERROR MESSAGES
!  none
! GTPack MODULES
!  none
! GTPack NOTEBOOKS 
!  none
! DESCRIPTION
! 
! LITERATURE
!
! TODO
! 
! PROBLEMS
!  
! Release
!
!--------------------------------------------------------------------------------
! SOURCE

*)
GTChangeRepInst[Rep_,verb_] := Module[{sonst}, 
	Unprotect[grpdgrp];
    Clear[grpdgrp];
    sonst=True;
    Which[	Rep=="O(3)",
    			grpdgrp=1;
    			sonst=False;
    			If[verb,Print["The standard representation has changed to O(3)"],None],
	 	  	Rep=="SU(2)",
	 	  		grpdgrp=2;
	 	  		sonst=False;
	 	  		If[verb,Print["The standard representation has changed to SU(2)"],None],
	 		Rep=="O(2)",
	 			grpdgrp=3;
	 			sonst=False;
	 			If[verb,Print["The standard representation has changed to O(2)"],None],
	 		Rep=="SU(2)xS",
	 			grpdgrp=5;
	 			sonst=False;
	 			If[verb,Print["The standard representation has changed to SU(2)xS"],None]
    ];
    If[sonst,Print["Error: Representation not known. Choose one of the following: O(3), SU(2) O(2) or SU(2)xS."];Abort[],None];
    Protect[grpdgrp];
]

(*
***)

(****g* /GTGroupFromGenerators
! NAME
!  GTGroupFromGenerators
! AUTHOR
!  W. Hergert
! PACKAGE
!  test.m -> basic.m
! MODIFICATION HISTORY
!  10/26/2012  : first version
!  24/02/2014  : use matrix multiplication instead of small circle
! USAGE
!  GTGroupFromGenerators[generators] gives a finite group from a list of generators.
! INPUT
!  set of generators
! OUTPUT
!  group
! ERROR MESSAGES
!  
! GTPack MODULES
!  GTGetSymbol, GTOrderOfElement
! GTPack NOTEBOOKS 
!  Wolfram : GT_Basic.nb
! DESCRIPTION
!  besides the methods to install groups from a list or from a multiplication table
!  this command allows to start from a set of generators. This might be helpfull if,
!  for example, you want to have C_3v with a different direction of the 3fold axis.
! LITERATURE
!
! TODO
!  not fully tested
! PROBLEMS
!  It should be tested inside if the result is a group and then the multiplication table 
!  should be printed, but the output should be a group or an empty list. GTGroupQ has 
!  to be changed a little bit.
! 
! Release
!
!--------------------------------------------------------------------------------
! SOURCE
*)


GTGroupFromGenerators[gen_] := Module[{grp, grp1, ln, ord, t, i, j, ordo, ordn,whinp,grpout},
 whinp=GTWhichInput[gen[[1]]];
 grp = GTGetMatrix[gen];
 ln = Length[grp];
 grp1 = {};
 ord = Map[GTOrderOfElement[#] &, grp];
 Do[t = grp[[i]]; grp1 = Append[grp1, GTSimplify[t]];
  Do[t = GTSimplify[t.t];
   grp1 = Append[grp1, t], {j, 1, ord[[i]] - 1}], {i, 1, ln}];
 grp = Union[grp1, SameTest -> (Re[Chop[N[Norm[#1 - #2]]]] < 10^-6 &)];
 ordn = Length[grp]; ordo = 0;
 While[ordn > ordo, ordo = Length[grp];
  grp1 = {};
  Do[Do[grp1 = Append[grp1, GTSimplify[grp[[i]].grp[[j]]]], {i, 1, ordo}], {j, 1, 
    ordo}]; grp1 = Union[grp1, SameTest -> (Re[Chop[N[Norm[#1 - #2]]]] < 10^-6 &)]; ordo = ordn;
  ordn = Length[grp1];
  grp = grp1];
  grpout=SortBy[grp1, GTOrderOfElement[#]*10^(1 - Det[GTGetMatrix[#]]) &];
  (*grpout=Sort[grp1, GTOrderOfElement[#1] < GTOrderOfElement[#2] &];*)
  Return[GTWhichOutput[GTSimplify[grpout],whinp]]];
 
(*
***) 

(****g* /GTInstallAxis
! NAME
!  GTInstallAxis
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  test.m -> Basic.m
! MODIFICATION HISTORY
!   02/24/2014: first version
! USAGE
!   GTInstallAxis[{label,coordinates}] installs symmetry elements according to a certain axis denoted by a label and defined by cartesian coordinates.
! INPUT
!   label, axis
! OPTIONS
!  GOVerbose
! OUTPUT
!  
! ERROR MESSAGES
!  
! GTPack MODULES
!  GTSimplify, GTSymbolInfoToSymbol, GTSymbolInfoToO3Matrix, GTSymbolInfoToSU2Matrix, GTSymbolInfoToEulerAngles
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! PROBLEMS
!
! Release
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTConvertToSU2Info[syminf_] := {syminf[[1]], Abs[syminf[[2]]], If[syminf[[2]] > 0, syminf[[3]], 2*Abs[syminf[[2]]] - syminf[[3]]], syminf[[4]]}

GTInstallAxis[label_, axis0_,OptionsPattern[]] := Module[{dge,inv,pos, newlist, axis,list, O3name, O2matrix, O3matrix, SU2name, SU2matrix,SU2xSmatrix,elmnew,elmSU2new,elmSU2xSnew,quaternion, ea, elmO2new},
  (*--- add the new axis to the axeslist ---*)
  axis = ToRadicals@FullSimplify@Normalize[axis0];
  pos = Position[axeslist, label];
  If[Length[pos] > 0, Print["Error: The axis label " <> label <> " already exists!"]; Abort[]];
  pos = Position[axeslist, axis];
  If[Length[pos] > 0, Print["Error: The axis already exists and is labeled with " <> axeslist[[pos[[1, 1]], 1]] <> "!"]; Abort[]];
  newlist = Append[axeslist, {ToString[label], axis}];
  Unprotect[axeslist];
  axeslist = newlist;
  Protect[axeslist];
  If[OptionValue[GOVerbose],Print["The axis " <> label <>" was added to the variable axeslist!"]];
  
  (*--- add the new axis to the axeslist ---*)
  Unprotect[elm,elmSU2,elmSU2xS,elmo2,Nelm,Nelmo2,NelmSU2,NelmSU2xS];
  list = Flatten[{list0 /. d -> 1, list0 /. d -> -1, list0 /. d -> I, list0 /. d -> -I} /. axistmp -> axis, 1];
  (*--- O(3)-name, O(3)-matrix, SU(2)-name, SU(2)-matrix ---*)
  Do[
  	If[Element[list[[i,4]], Reals],dge=1;inv=list[[i,4]],dge=-1;inv=list[[i,4]]/I];
  	O3name = ToString[GTSymbolInfoToSymbol[list[[i]]]];
  	O3matrix = GTSimplify@GTSymbolInfoToO3Matrix[list[[i]]];
  	O2matrix = {O3matrix[[1, 1 ;; 2]], O3matrix[[2, 1 ;; 2]]};
  	SU2name = ToString[GTSymbolInfoToSymbol[GTConvertToSU2Info@list[[i]]]];
  	SU2matrix = GTSimplify@GTSymbolInfoToSU2Matrix[GTConvertToSU2Info@list[[i]]];
  	SU2xSmatrix = {{SU2matrix[[1,1]],SU2matrix[[1,2]],0},{SU2matrix[[2,1]],SU2matrix[[2,2]],0},{0,0,inv}};
  	quaternion = GTSimplify@{-Re[SU2matrix[[1,1]]],{Im[SU2matrix[[2,1]]],-Re[SU2matrix[[2,1]]],Im[SU2matrix[[1,1]]]}};
  	ea= GTSimplify@GTSymbolInfoToEulerAngles[list[[i]]];
  	elmnew = {O3name,O3matrix,quaternion,ea,dge};
  	elmO2new = {O3name,O2matrix,quaternion,ea,dge};
  	elmSU2new = {O3name,SU2matrix,quaternion,ea,dge};
  	elmSU2xSnew = {O3name,SU2xSmatrix,quaternion,ea,dge};
  	elm = Append[elm,elmnew];
  	If[{0,0,1}.(O3matrix.{0,0,1})==1,elmo2 = Append[elmo2,elmO2new]];
  	elmSU2 = Append[elmSU2,elmSU2new];
  	elmSU2xS = Append[elmSU2xS,elmSU2xSnew];
  	
    Nelm = Append[Nelm,Rationalize[N[N[elmnew, 15], 8]]];
  	If[{0,0,1}.(O3matrix.{0,0,1})==1,Nelmo2 = Append[Nelmo2,Rationalize[N[N[elmO2new, 15], 8]]]];
  	NelmSU2 = Append[NelmSU2,Rationalize[N[N[elmSU2new, 15], 8]]];
  	NelmSU2xS = Append[NelmSU2xS,Rationalize[N[N[elmSU2xSnew, 15], 8]]];
  	
  	, {i, 1, Length[list]}];
  	  	
  	
  	  	
  Protect[elm,elmSU2,elmSU2xS,elmo2,Nelm,Nelmo2,NelmSU2,NelmSU2xS];  
]
(*
***)

(****g* /GTInstallGroup
!
! NAME
!  GTInstallGroup.m
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Install.m
! MODIFICATION HISTORY
!  2/23/12    : initial documentation 
!  29/05/2015 : extension to input of point group name in HM notation
!  24/08/2016 : extension to input of space groups  
! USAGE
!  GTInstallGroup[group] gives a faithful representation of a crystallographic point group.
! INPUT
!  group
!  The identifier "key" of the point group can be taken from the Palette to GTPack.
!  OptionsPattern: GOVerbose
! OUTPUT
!  list of matrices corresponding to the representation 
! ERROR MESSAGES
!  none
! GTPack MODULES
!  GTGroupFromGenerators, GTChangeRepInst, GTSpaceGroups
! GTPack NOTEBOOKS 
!  none
! DESCRIPTION
! 
! LITERATURE
!
! TODO
! 
! PROBLEMS
!  
! Release
!
!--------------------------------------------------------------------------------
! 
! SOURCE
*)
Unprotect[generators]
Clear[generators]

generators ={{{"C1"}, {Ee}}, {{"Ci", "S2"}, {IEe}}, {{"C2"}, {C2z}}, {{"Cs","C1h"}, {IC2z}}, {{"C3"}, {C3z}}, {{"S4"}, {IC4y}}, {{"C4"}, 
	{C4z}}, {{"D2", "V"}, {C2x, C2y}}, {{"C2h"}, {C2z,IEe}}, {{"C2v"}, {C2z, IC2y}}, {{"C6"}, {C6z}}, {{"C3i", "S6"}, {C3z, IEe}},
	{{"D3"}, {C3z, C2x}}, {{"C3v"}, {C3z,IC2y}}, {{"C3h"}, {IC6z}}, {{"C4h"}, {C4z, IEe}}, {{"D2d","Vd"}, {IC4y, C2x}}, {{"D4"}, {C4y,C2z}},
	{{"C4v"}, {C4z, IC2y}}, {{"D2h", "Vh"}, {C2x, C2e,IEe}}, {{"T"}, {C3\[Alpha], C2z}}, {{"D6"}, {C6z,C2x}}, {{"C6h"}, {C6z, IEe}},
	{{"C6v"}, {C6z,IC2y}}, {{"D3d"}, {C3z, IEe, C2x}}, {{"D3h"}, {IC6z,C2x}}, {{"D4h"}, {C4z, C2a, IEe}}, {{"Th"}, {C2z,C3\[Alpha], IEe}},
	{{"O"}, {C2z, C3\[Alpha],C2a}}, {{"Td"}, {C2z, C3\[Alpha], IC4x}}, {{"D6h"}, {C6z,C2x, IEe}}, {{"Oh"}, {C4z, C3\[Alpha], IEe}}}

Protect[generators]





GTInstallGroup::badsymbol =  "The specified group `1` is not known.";


GTInstallGroup[key0_,OptionsPattern[]] := Module[{pos,key,postmp,rep,inp},
	inp=GTSpaceGroups[key0, GOVerbose -> False];
	If[Length[inp]==0,
		(*--- convert input to string ---*)
		key=GTNotationToSFL[key0];
		(*--- find the position of the group within the list of generators ---*) 
		postmp = Flatten[Position[generators, key]];
		If[Length[postmp]==0,
(*			Print[$ErrorGTInstallGroup] *)
             Message[GTInstallGroup::badsymbol,key0];Abort[],pos=First[postmp]];
		(*--- check if the user wants to install a double group containing the inversion ---*)
		rep = OptionValue[GORepresentation];
		If[OptionValue[GORepresentation]=="SU(2)",
			If[Length[Position[generators[[pos,2]],IEe]]==0,
				rep="SU(2)",
			
				If[OptionValue[GOVerbose],Print["Warning: the installed group contains the inversion IEe. The representation is switched to SU(2)xS!"]];
				rep="SU(2)xS"]
		],
	
	If[Length[gtspacegroupelements[[inp[[1]],2]]]==0,Print["In future, this will install the "<>inp[[2]]<>" space group "<>ToString[inp[[3]]]]];
	Return[gtspacegroupelements[[inp[[1]],2]]]
	];
	
	GTChangeRepInst[rep,OptionValue[GOVerbose]]; 
 	GTGroupFromGenerators[Flatten[{generators[[pos, 2]],DEe}]]
 ]



(*
***)

(****g* /GTReinstallAxes
!
! NAME
!  GTReinstallAxes
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Install.m
! MODIFICATION HISTORY
!  6/26/15 : initial documentation 
! USAGE
!  Reinstalls element definitions for passive or active rotation convention. 
! EXAMPLE
!  
! INPUT
!  string
! OUTPUT
!  none 
! ERROR MESSAGES
!  none
! GTPack MODULES
!  GTInstallAxis
! GTPack NOTEBOOKS 
!  none
! DESCRIPTION
! 
! LITERATURE
!
! TODO
! 
! PROBLEMS
!  
! Release
!
!--------------------------------------------------------------------------------
! SOURCE 

*)

GTReinstallAxes[act_] := Module[
  {axesold, elmold, elmo2old, elmSU2old, elmSU2xSold,Nelmold, Nelmo2old, NelmSU2old, NelmSU2xSold,tmp1,tmp2,tmp3,tmp4,tmp1b},
  Unprotect[gtactvpsv];
  error=True;
  Which[
  	act == "active", 
  	    If[gtactvpsv==-1,Print["The active definition of rotation matrices is already installed."];Return[]];
  		gtactvpsv = -1; 
  		error=False,
  	act == "passive", 
  	    If[gtactvpsv==1,Print["The passive definition of rotation matrices is already installed."];Return[]];
  		gtactvpsv = 1; 
  		error=False
  ];
  If[error,Print[$ErrorReinstallAxes];Abort[]];
  Protect[gtactvpsv]; 
  axesold = axeslist;

  elmold = elm[[1 ;; 4]];
  elmo2old = elmo2[[1 ;; 4]];
  elmSU2old = elmSU2[[1 ;; 4]];
  elmSU2xSold = elmSU2xS[[1 ;; 4]];
  
  Nelmold = Nelm[[1 ;; 4]];
  Nelmo2old = Nelmo2[[1 ;; 4]];
  NelmSU2old = NelmSU2[[1 ;; 4]];
  NelmSU2xSold = NelmSU2xS[[1 ;; 4]];
  
  Unprotect[elm,elmo2,elmSU2,elmSU2xS,Nelm,Nelmo2,NelmSU2,NelmSU2xS,axeslist];

  elm = elmold;
  elmo2 = elmo2old;
  elmSU2 = elmSU2old;
  elmSU2xS = elmSU2xSold;

  Nelm = Nelmold;
  Nelmo2 = Nelmo2old;
  NelmSU2 = NelmSU2old;
  NelmSU2xS = NelmSU2xSold;

  axeslist={};
  Protect[elm,elmo2,elmSU2,elmSU2xS,Nelm,Nelmo2,NelmSU2,NelmSU2xS,axeslist];
    
  tmp1 = PrintTemporary["Installation of the Mathematica Group Theory Package GTPack"];
  tmp1b = PrintTemporary["for "<>act<>" rotations."];
  tmp2 = PrintTemporary["--------------------------------------------------------------"];
  tmp3 = PrintTemporary["install symmetry elements"];
  Do[
  	tmp4 = PrintTemporary[ToString[axesold[[i,1]]]<>"-Axis"];
  	GTInstallAxis[ToString[axesold[[i,1]]],axesold[[i,2]],GOVerbose->None];
  	NotebookDelete[tmp4];
  ,{i,1,Length[axesold]}];
  NotebookDelete[tmp1];
  NotebookDelete[tmp2];
  NotebookDelete[tmp3];
  ]
  
(*
***)


(****g* /GTWhichAxes
!
! NAME
!  GTWhichAxes
! AUTHOR
!  W. Hergert
! PACKAGE
!  Install.m
! MODIFICATION HISTORY
!  7/15/15 : initial documentation 
! USAGE
!  With GTReinstallAxes it is possible to switch between an active or passive defintion 
!  of symmetry transformations. The actual setting is stored in the protected variable
!  gtactvpsv. Sometime it might be good for the user to remember which definition he is 
!  actually using. Because the name of the variable is hard to remember, this small command 
!  will be perhaps a good alternative.  
! EXAMPLE
!  
! INPUT
!  no input
! OUTPUT
!  actual setting of the dfinition fpor symmetry transformations
! ERROR MESSAGES
!  none
! GTPack MODULES
!  none
! GTPack NOTEBOOKS 
!  none
! DESCRIPTION
! 
! LITERATURE
!
! TODO
! 
! PROBLEMS
!  
! Release
!
!--------------------------------------------------------------------------------
! SOURCE 

*)


GTWhichAxes[] := Module[{},
 (* If[gtactvpsv == -1, Print["active definition"], *)
  If[gtactvpsv == -1, Print["active definition"], 
   Print["passive definition"]]
  ]
  
 (*
***)
 
  
(****g* /GTSetIndex
!
! NAME
!  GTSetIndex.m
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Install.m
! MODIFICATION HISTORY
!  2/23/12 : initial documentation 
! USAGE
!  Switching between matrix representation of group elements. 
! EXAMPLE
!  GTSetIndex[ind]
! INPUT
!  ind=1 for O(3), ind=2 for SU(2), ind=3 for O(2), ind=4 for permutation matrices
! OUTPUT
!  none 
! ERROR MESSAGES
!  none
! GTPack MODULES
!  none
! GTPack NOTEBOOKS 
!  none
! DESCRIPTION
! 
! LITERATURE
!
! TODO
! 
! PROBLEMS
!  
! Release
!
!--------------------------------------------------------------------------------
! SOURCE 

*)

GTSetIndex[ind_] := Module[{},
	Unprotect[grpdgrp];
	Clear[grpdgrp];
	grpdgrp=ind;
    Protect[grpdgrp];
]

(*
***)

(****g* /GTTableToGroup
!
! NAME
!  GTTableToGroup.m
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Install.m
! MODIFICATION HISTORY
!  * 23.02.2012 : initial documentation
!  * 31.12.2012 : renamed from GTTableToRepresentation to GTTableToGroup 
!  * 10.4.2015  : The Unprotect command was at a wrong place. Conflicts if Verbose->False
!  * 10.02.2016 : output of table changed
!  * 07.09.2018 : srious bug detected and removed (W)
!
! USAGE
!  GTTableToGroup[list of elements,multiplication table] gives a faithful representation of an arbitrary group from a given list of elements and a multiplication table, using permutation matrices.
!
! INPUT
!  * list of elements
!  * multiplication table
! OUTPUT
!  * group in a matrix representation(permutation matrices)
! GTPack OPTIONS
!  * GOVerbose
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTSetIndex, GTChangeRepresentation
! GTPack NOTEBOOKS 
!  none
! DESCRIPTION
! 
! LITERATURE
!
! TODO
! 
! PROBLEMS
!  07.09.2018: The module did'nt really work. It worked accedentialy for the Klein's 4-group because the
!  indentity element was in all diagonal elements of the multiplication table.If this is not the case  no real
!  identity element will be generated. It is assumed, that the (1,1) element of the multiplication table contains 
!  the identity element. All generated matrices will be multiplied by the inverse of the matrix corresponding to this
!  element. Then all is fine.
! Release
!
!--------------------------------------------------------------------------------
! SOURCE
!
*)

GTTableToGroup[el_, mt_,OptionsPattern[]] := 
 Module[{sym2num, mtidx, mats, mat2sym,i,j,k,mtprint,mt1,grps1,corr}, 
  (*--- Check for consistency ---*)
  (*--- Check for group axioms ---*)
  GTChangeRepresentation["Permutation",GOVerbose->OptionValue[GOVerbose]];
  sym2num[a_] := Position[el, a][[1, 1]];
  mtidx = Map[sym2num, mt, {2}];
  mats = Table[
    Table[0, {i, Length[el]}, {j, Length[el]}], {k, Length[el]}];
  Do[
   Do[
   	(*--- Iteration over columns of the MT ---*)
    mats[[i]][[j]][[mtidx[[j]][[i]]]] = 1, {j, 1, Length[el]} 
    ],
    (*--- Iteration over matrices ---*)
    {i, 1, Length[el]} 
   ];
   (*--- 7.9.18 correction to get really a group---*)
   corr = Inverse[mats[[1]]]; mats = Inverse[corr].# & /@ mats;
  (*--- Map matrix to symbol ---*)
  mat2sym[mat_] := el[[Position[mats, mat][[1, 1]]]]; 
  If[OptionValue[GOVerbose],
     Do[
     	Print[mat2sym[mats[[i]]], "=", MatrixForm [mats[[i]]]]
     , {i, 1, Length[el]}];
     mtprint = Table[mat2sym[mats[[j]].mats[[i]]], {j, Length[el]}, {i, Length[el]}];
     mt1 = Transpose[Join[{el}, mtprint]];
     grps1 = Prepend[el, " "]; mt1 = Join[{grps1}, mt1];
     Print[Grid[mt1, Frame -> All,Dividers -> {{2 -> Red}, {2 -> Red}}, Background -> {{1 -> Pink}, {1 -> Pink}, {1, 1} -> Yellow}]];  
     (* Print[TableForm[Table[mat2sym[mats[[j]].mats[[i]]]
     	   , {j, Length[el]}, {i, Length[el]}] , TableHeadings -> {el, el}]]*)
     ];
     Unprotect[elmown];
     Do[
  	    If[Length[Position[elmown, el[[i]]]]==0,
  		   elmown=Append[elmown,{el[[i]],mats[[i]]}]];
  	 ,{i,1,Length[mats]}];
     Protect[elmown];
  mats
  ]

(*
***)	
	
(****g* /GTWhichIndex
!
! NAME
!  GTWhichIndex.m
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Install.m
! MODIFICATION HISTORY
!  2/23/12 : initial documentation 
! USAGE
!  GTWhichIndex gives the current value of grpdgrp
! EXAMPLE
!  GTWhichIndex
! INPUT
!  none
! OUTPUT
!  ind=1 for SO(3), ind=2 for SU(2), ind=3 for O(2), ind=4 for permutation matrices
! ERROR MESSAGES
!  none
! GTPack MODULES
!  none
! GTPack NOTEBOOKS 
!  none
! DESCRIPTION
! 
! LITERATURE
!
! TODO
! 
! PROBLEMS
!  
! Release
!
!--------------------------------------------------------------------------------
! SOURCE 

*)

GTWhichIndex := grpdgrp

(*
***)





(*-------------------------- Attributes ------------------------------*)

Attributes[GTInstallGroup]={Protected, ReadProtected}
Attributes[GTInstallAxis]={Protected, ReadProtected}
Attributes[GTTableToGroup]={Protected, ReadProtected}
Attributes[GTChangeRepresentation]={Protected, ReadProtected}
(* Attributes[GTGroupFromGenerators]={Protected, ReadProtected} *)
Attributes[GTWhichRepresentation]={Protected, ReadProtected}
 
End[] (* End Private Context *)

EndPackage[]
