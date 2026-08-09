(* ::Package:: *)

(****m* /Basic.m
!
! NAME
!  Basic.m
! AUTHOR
!  W. Hergert, M. Geilhufe, S. Schenk, M. Daene
! MODIFICATION HISTORY
!  2/24/12 : initial documentation 
! USAGE
!  Contains basic modules of group theory.
! 
! ERROR MESSAGES
!  
! GTPack MODULES
!
! --- classes ---
! GTClasses gives the conjugacy classes of a group
! GTClassMult calculates the class constants
! GTClassMultTable calculates a table of class constants
! 
! --- characterization of groups and elements ---
! GTCenter determines the center of a group
! GTConjugacyClass constructs the set of all subgroups conjugated to the subgroup or group
! GTConjugateElement gives the conjugate element of an element with respect to a third element
! GTGenerators gives the generators of a certain group
! GTGetInvSubgroup gives an invariant subgroup of index n.
! GTGetSubGroups finds all subgroups of a group
! GTGroupOrder gives the order of a group
! GTInverseElement gives the inverse element of a symmetry element
! GTInvSubGroups gives the invariant subgroups of a group 
! GTLeftCosets gives the left coset of a group according to a subgroup
! GTMultTable gives the multiplication table for a group
! GTNormalizer constructs the normalizer of group with respect to its supergroup	
! GTOrderOfElement gives the order of a group element
! GTProductGroup forms the product of group1 and group2 and checks, if the product of the two groups forms a correct direct or semidirect product
! GTQuotientGroup gives the multiplication table of the quotient group of group and normal divisor
! GTRightCosets gives the right coset of a group according to a subgroup
!   
! --- logical group characterization ---  
! GTAbelianQ gives True if a group is an abelian group, and gives False otherwise
! GTCyclicQ gives True if group is a cyclic group, and gives False otherwise
! GTGroupQ gives True if a set of symmetry elements forms a group, and gives False otherwise
! GTInvSubGroupQ gives True if the group with smaller order is an invariant subgroup of the group with larger order, and gives False otherwise
! GTProductGroupQ gives True if group is a product group of the groups group1 and group2, and False otherwise
! GTQuotientGroupQ gives True if group1 and group2 can form a quotient group, and False otherwise
! GTSelfAdjointQ gives True if an element of a group is self-adjoint, and gives False otherwise
! GTSubGroupQ gives True if the group with smaller order is a subgroup of the group with larger order, and gives False otherwise
!
! --- Symbols, matrices, quaternions and Euler angles ---
! GTAllSymbols gives a list of all implemented symbols for point group elements 
! GTGetEulerAngles gives the Euler angles of a symmetry element
! GTGetMatrix gives a representation matrix of a symmetry element
! GTGetQuaternion gives the quaternion of a symmetry element
! GTGetRotationMatrix gives the 3-dimensional rotation matrix of a symmetry element
! GTGetSU2Matrix gives the SU(2) matrix of a symmetry element
! GTGetSymbol gives the symbol of a symmetry element
! GTgmat performs a multiplication
! GTSetMultiplication multplies each symmetry element of set1 with each symmetry element of set2 and gives the set of distinct elements
! GTType transformes a point group element to a space group element
! GTWhichInput determines the type of a symmetry element (Symbol, Matrix, Quaternion, Euler Angles or space group elements)
! GTWhichOutput transforms a symmetry element to either symbol, matrix, quaternion, a set of Euler Angles or a space group element, according to the given type index
! 
! --- Coordinate transformations ---  
! GTTransformation applies a coordinate transformation to a given vector
! GTTransformationOperator applies a coordinate transformation to a certain function with given arguments
!
! --- Magnetic groups --- 
!
! GTMagnetic sets the variable GTmaggrpq and determines if magnetic point group have to be taken into account
! GTMagneticQ  gives True if if magnetic groups are taken into account, and gives False otherwise
!
! --- Internal --- 
!
! GTSymbolToSymbolInfo 
!
! GTPack NOTEBOOKS 
!  none in the moment
!
! DESCRIPTION
!  Basic.m contains all functions for basic group theoretical calculations that are supported by the GroupTheory package. There is a group multiplication 
!  defined by the operator \[SmallCircle].
!
! LITERATURE
! 
! TODO
!  
! PROBLEMS
!
!--------------------------------------------------------------------------------
!
***)

BeginPackage["GroupTheory`Basic`",{"GroupTheory`Symbols`","GroupTheory`Install`","GroupTheory`Auxiliary`"}]

(*-------------------------- classes -----------------------------*)
 GTClasses                ::usage = "GTClasses[\*StyleBox[\"group\", \"TI\"]] gives the classes of conjugate elements of a \*StyleBox[\"group\", \"TI\"]."
 GTClassMult              ::usage = "GTClassMult[\*StyleBox[\"classes,i,j\", \"TI\"]] calculates the class constants \*Cell[BoxData[FormBox[SubscriptBox[\"C\", \"ijk\"],TraditionalForm]],\"InlineMath\"]."
 GTClassMultTable         ::usage = "GTClassMultTable[\*StyleBox[\"classes\", \"TI\"]] calculates a table of class constants."
 
(*-------------------------- characterization of groups and elements -------------*)
 GTCenter                 ::usage = "GTCenter[\*StyleBox[\"group\", \"TI\"]] gives the center of a \*StyleBox[\"group\", \"TI\"]."
 GTConjugacyClass         ::usage=  "GTConjugacyClass[\*StyleBox[\"group,subgroub\", \"TI\"]] constructs the set of all subgroups conjugated to a \*StyleBox[\"subgroub\", \"TI\"] of a \*StyleBox[\"group\", \"TI\"]."  
 GTConjugateElement       ::usage = "GTConjugateElement[\*StyleBox[\"T,X\", \"TI\"]] gives the conjugate element \*Cell[TextData[{Cell[BoxData[FormBox[\"X\", TraditionalForm]],\"InlineMath\"],Cell[BoxData[FormBox[\"T\", TraditionalForm]],\"InlineMath\"],Cell[BoxData[FormBox[SuperscriptBox[\"X\",RowBox[{\"-\",\"1\"}]], TraditionalForm]],\"InlineMath\"]}],\"Usage\"]."
 GTGenerators             ::usage = "GTGenerators[\*StyleBox[\"group\", \"TI\"]] gives the generators of a certain \*StyleBox[\"group\", \"TI\"]."
 GTGetInvSubgroup         ::usage = "GTGetInvariantSubgroup[\*StyleBox[\"group, classes, index n\", \"TI\"]]  gives an invariant subgroup  of \*StyleBox[\"index n\", \"TI\"]."
 GTGetSubGroups           ::usage = "GTGetSubGroups[\*StyleBox[\"group\", \"TI\"]] finds all subgroups of a \*StyleBox[\"group\", \"TI\"]."
 GTGroupOrder             ::usage = "GTGroupOrder[\*StyleBox[\"group\", \"TI\"]] gives the order of a \*StyleBox[\"group\", \"TI\"]."
 GTInverseElement         ::usage = "GTInverseElement[\*StyleBox[\"symmetry element\", \"TI\"]] gives the inverse element of a certain \*StyleBox[\"symmetry element\", \"TI\"]."
 GTInvSubGroups           ::usage = "GTInvSubGroups[\*StyleBox[\"group\", \"TI\"]] gives the invariant subgroups of a \*StyleBox[\"group\", \"TI\"]." 
 GTLeftCosets             ::usage = "GTLeftCosets[\*StyleBox[\"group,subgroub\", \"TI\"]] gives the left coset of a \*StyleBox[\"group\", \"TI\"] according to a \*StyleBox[\"subgroup\", \"TI\"]."
 GTMultTable              ::usage = "GTMultTable[\*StyleBox[\"group\", \"TI\"]] gives the multiplication table for a \*StyleBox[\"group\", \"TI\"]."
 GTNormalizer             ::usage=  "GTNormalizer[\*StyleBox[\"supergroup,group\", \"TI\"]] constructs the normalizer of a \*StyleBox[\"group\", \"TI\"] with respect to its \*StyleBox[\"supergroup\", \"TI\"]."	
 GTOrderOfElement         ::usage = "GTOrderOfElement[\*StyleBox[\"element\", \"TI\"]] gives the order of a group \*StyleBox[\"element\", \"TI\"]."
 GTProductGroup           ::usage=  "GTProductGroup[\*StyleBox[\"group1,groub2\", \"TI\"]] forms the product of \*StyleBox[\"group1\", \"TI\"] and \*StyleBox[\"groub2\", \"TI\"] and checks, if the product of the two groups forms a correct direct or semidirect product."
 GTQuotientGroup          ::usage=  "GTQuotientGroup[\*StyleBox[\"group,normal divisor\", \"TI\"]] gives the multiplication table of the quotient group of \*StyleBox[\"group\", \"TI\"] and \*StyleBox[\"normal divisor\", \"TI\"]."
 GTRightCosets            ::usage = "GTRightCosets[\*StyleBox[\"group,subgroub\", \"TI\"]] gives the right coset of a \*StyleBox[\"group\", \"TI\"] according to a \*StyleBox[\"subgroup\", \"TI\"]."
   
(*-------------------------- logical group characterization -----------------------------------*)  
 GTAbelianQ               ::usage = "GTAbelianQ[\*StyleBox[\"group\", \"TI\"]] gives True if \*StyleBox[\"group\", \"TI\"] is an abelian group, and gives False otherwise."
 GTCyclicQ                ::usage = "GTCyclicQ[\*StyleBox[\"group\", \"TI\"]] gives True if \*StyleBox[\"group\", \"TI\"] is a cyclic group, and gives False otherwise."
 GTGroupQ                 ::usage = "GTGroupQ[\*StyleBox[\"group\", \"TI\"]] gives True if \*StyleBox[\"set of symmetry elements\", \"TI\"] forms a group, and gives False otherwise."
 GTInvSubGroupQ           ::usage = "GTInvSubGroupQ[\*StyleBox[\"group1,group2\", \"TI\"]] gives True if the \*StyleBox[\"group\", \"TI\"] with smaller order is an invariant subgroup of the \*StyleBox[\"group\", \"TI\"] with larger order, and gives False otherwise."
 GTProductGroupQ          ::usage=  "GTProductGroupQ[\*StyleBox[\"group,group1,group2\", \"TI\"]] gives True if \*StyleBox[\"group\", \"TI\"] is a product group of the groups \*StyleBox[\"group1\", \"TI\"] and \*StyleBox[\"group2\", \"TI\"], and False otherwise."
 GTQuotientGroupQ         ::usage=  "GTQuotientGroupQ[\*StyleBox[\"group1,group2\", \"TI\"]] gives True if \*StyleBox[\"group1\", \"TI\"] and \*StyleBox[\"group2\", \"TI\"] can form a quotient group, and False otherwise."
 GTSelfAdjointQ           ::usage = "GTSelfAdjointQ[\*StyleBox[\"group,element\", \"TI\"]] gives True if an \*StyleBox[\"element\", \"TI\"] of a \*StyleBox[\"group\", \"TI\"] is self-adjoint, and gives False otherwise."
 GTSubGroupQ              ::usage = "GTSubGroupQ[\*StyleBox[\"group1,group2\", \"TI\"]] gives True if the \*StyleBox[\"group\", \"TI\"] with smaller order is a subgroup of the \*StyleBox[\"group\", \"TI\"] with larger order, and gives False otherwise."

 (*------------------- Symbols, matrices, quaternions and Euler angles -------*)
 GTAllSymbols             ::usage = "GTAllSymbols gives a list of all implemented symbols for point group elements." 
 GTGetEulerAngles         ::usage = "GTGetEulerAngles[\*StyleBox[\"symmetry element\", \"TI\"]] gives the Euler angles corresponding to a \*StyleBox[\"symmetry element\", \"TI\"]."
 GTGetMatrix              ::usage = "GTGetMatrix[\*StyleBox[\"symmetry element\", \"TI\"]] gives a matrix corresponding to a \*StyleBox[\"symmetry element\", \"TI\"]."
 GTGetQuaternion          ::usage = "GTGetQuaternion[\*StyleBox[\"symmetry element\", \"TI\"]] gives the quaternion corresponding to a \*StyleBox[\"symmetry element\", \"TI\"]."
 GTGetRotationMatrix      ::usage = "GTGetRotationMatrix[\*StyleBox[\"symmetry element\", \"TI\"]] gives the 3-dimensional rotation matrix of a \*StyleBox[\"symmetry element\", \"TI\"]."
 GTGetSU2Matrix           ::usage = "GTGetSU2Matrix[\*StyleBox[\"symmetry element\", \"TI\"]]  gives the SU(2) matrix of a certain \*StyleBox[\"symmetry element\", \"TI\"]."
 GTGetSymbol              ::usage = "GTGetSymbol[\*StyleBox[\"symmetry element\", \"TI\"]] gives the symbol of a \*StyleBox[\"symmetry element\", \"TI\"]."
 GTgmat                   ::usage = "GTgmat[\*StyleBox[\"element1,element2\", \"TI\"]] performs a multiplication of \*StyleBox[\"element1\", \"TI\"] and \*StyleBox[\"element2\", \"TI\"]."
 GTSetMultiplication      ::usage = "GTSetMultiplication[set1,set2] multplies each symmetry element of set1 with each symmetry element of set2 and gives the set of distinct elements."
 GTType                   ::usage = "GTType[\*StyleBox[\"point group element\", \"TI\"]] transformes a \*StyleBox[\"point group element\", \"TI\"] to a space group element."
 GTWhichInput             ::usage = "GTWhichInput[\*StyleBox[\"symmetry element\", \"TI\"]] determines the type of a \*StyleBox[\"symmetry element\", \"TI\"] (Symbol, Matrix, Quaternion, Euler Angles or space group elements)."
 GTWhichOutput            ::usage = "GTWhichOutput[\*StyleBox[\"symmetry element,type index\", \"TI\"]] transforms an \*StyleBox[\"symmetry element\", \"TI\"] to either a symbol, a matrix, a quaternion, a set of Euler Angles or a space group element, according to the given \*StyleBox[\"type index\", \"TI\"]."
 
(*-------------------------- Coordinate transformations --------------------------*)  
 GTTransformation         ::usage = "GTTransformation[\*StyleBox[\"transformation,vector\", \"TI\"]] applies a coordinate \*StyleBox[\"transformation\", \"TI\"] to a given \*StyleBox[\"vector\", \"TI\"]."
 GTTransformationOperator ::usage = "GTTransformationOperator[\*StyleBox[\"transformation,function,arguments\", \"TI\"]] applies a coordinate \*StyleBox[\"transformation\", \"TI\"] to a certain \*StyleBox[\"function\", \"TI\"] with given \*StyleBox[\"arguments\", \"TI\"]."

 GTMagnetic               ::usage = "GTMagnetic[\*StyleBox[\"logical\", \"TI\"]] sets the variable GTmaggrpq and determines if magnetic point group are taken into account."
 GTMagneticQ              ::usage = "GTMagneticQ[] gives True if if magnetic groups are taken into account, and gives False otherwise."
 (*--------------------------------------------------*)
 (*-------------------------- Internal --------------*)
 (*--------------------------------------------------*)

 GTSymbolToSymbolInfo     ::usage = ""
 
 (*
***)
 GTGenerators = GroupTheory`Auxiliary`Private`GTGenerators
 
 GTGroupFromGenerators = GroupTheory`Install`Private`GTGroupFromGenerators
 
 GTWhichRepresentation = GroupTheory`Install`Private`GTWhichRepresentation
(*--------------------------- Options --------------------------------*)
 Options[GTCenter]              ={GOFast->GOFastValue} 
 Options[GTClasses]             ={GOFast->GOFastValue}
 Options[GTClassMult]           ={GOFast->GOFastValue}
 Options[GTClassMultTable]      ={GOVerbose->True,GOFast->GOFastValue}
 Options[GTConjugacyClass]      ={GOFast->GOFastValue}
 Options[GTGenerators]          ={GOFast->GOFastValue}
 Options[GTGetSubGroups]        ={GOGroupOrder->1,GOFast->GOFastValue}
 Options[GTgmat]                ={GOVerbose->False}
 Options[GTGroupOrder]          ={GOVerbose->True,GOFast->GOFastValue}
 Options[GTGroupQ]              ={GOMethod->"Numeric"}
 Options[GTGetInvSubGroup]      ={MaxIterations -> 10000, GOVerbose -> True};
 Options[GTInvSubGroupQ]        ={GOVerbose->True}
 Options[GTInverseElement]      ={GOLattice->{{1,0,0},{0,1,0},{0,0,1}}}
 Options[GTInvSubGroups]        ={GOFast->GOFastValue}
 Options[GTLeftCosets]          ={GOFast->GOFastValue}
 Options[GTMultTable]           ={GOVerbose->True,GOFast->GOFastValue}
 Options[GTNormalizer]          ={GOFast->GOFastValue}
 Options[GTProductGroup]        ={GOFast->GOFastValue}
 Options[GTProductGroupQ]       ={GOVerbose -> False,GOFast->GOFastValue}
 Options[GTQuotientGroup]       ={GOVerbose -> True,GOQuotientGroup->False,GOFast->GOFastValue}
 Options[GTQuotientGroupQ]      ={GOVerbose -> True,GOFast->GOFastValue}
 Options[GTRightCosets]         ={GOFast->GOFastValue}
 Options[GTSelfAdjointQ]        ={GOFast->GOFastValue}
 Options[GTSubGroupQ]        ={GOFast->GOFastValue}
 
 
 
 
Begin["`Private`"] (* Begin Private Context *) 


(****c* /GTAbelianQ
! NAME
!  GTAbelianQ
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  Wolfram : 10/25/2012
!  12/31/2012: Wolfram: old version removed, output of multiplication table suppressed
! USAGE
!  GTAbelianQ[group] gives True if group is an abelian group, and gives False otherwise.
! INPUT
!  group
! OUTPUT
!  logical
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetMatrix, GTMultTable
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!
! PROBLEMS
! 
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTAbelianQ[group_]:=Module[{group1,mtables},
	group1 = GTGetMatrix[group];
	mtables=GTMultTable[group1,GOVerbose->False];
	If[mtables===Transpose[mtables],True,False]]
			 
(*
***)


(****c* /GTAllSymbols
! NAME
!  GTAllSymbols
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  15.10.2015: Rewriting the command using ToExpression instead of GTGetMatrix -> much faster!
! USAGE
!  GTAllSymbols gives a list of all implemented symbols for point group elements.
! INPUT
!  
! OUTPUT
!  list of all implemented symbols
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!   
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTAllSymbols[] := Table[ToExpression[elm[[i, 1]]], {i, 1, Length[elm]}]

(*
***)

(****c* /GTCenter
! NAME
!  GTCenter
! AUTHOR
!  W. Hergert
! PACKAGE
!  basic.m
! MODIFICATION HISTORY
!   13/11/2014 : moved to basic.m
!   13/11/2014 : variables are defined properly
!   11/30/2012  : first version
! USAGE
!  GTCenter[group] determines the center of a group.
! INPUT
!  group
! OUTPUT
!  center
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTSelfAdjointQ
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  not fully tested
!  documentation  
! RELEASE
!
! PROBLEMS
! 
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTCenter[group_,OptionsPattern[]] := Module[{ce, l, i}, 
  If[OptionValue[GOFast],None,If[GTGroupQ[group],None,Print[$ErrorNoGroup];Abort[]];];
  ce = {}; 
  l = Length[group]; 
  Do[If[GTSelfAdjointQ[group, group[[i]]],ce = Append[ce, group[[i]]], Null] , {i, 1, l}]; Return[ce]
 ]

(*
***) 

(****r* /GTClasses
! NAME
!  GTClasses
! AUTHOR
!  W. Hergert, M. Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  15/05/2013 Change of Structure, Error Messages
!  11/05/2020 fixed issue with doubly occuring elements for group C5
! USAGE
!  GTClasses[group] gives the conjugacy classes of a groups.
! INPUT
!  group - list of group elements
! OUTPUT
!  list of classes
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetMatrix, GTGroupQ
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!
! LITERATURE
!
! TODO
!  
! RELEASE
! 
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTClasses[group_,OptionsPattern[]] :=Module[{g, tb1, grpma,grpmai, i, j,run,dimmat,classes,clout,whinp,identityposition},
	If[OptionValue[GOFast],run=True,run=GTGroupQ[group]];
	If[run,
			whinp = GTWhichInput[group[[1]]];
			g=Length[group];
			grpma=GTGetMatrix[group];
			grpmai = Table[Inverse[grpma[[i]]], {i, 1, g}];
			dimmat=Length[grpma[[1]]];
			tb1=GTSimplify[Table[grpma[[i]].grpma[[j]].grpmai[[i]],{j,1,g},{i,1,g}]];
			tb1 = Table[Union[tb1[[i]],SameTest -> (Norm[N[#1 - #2], "Frobenius"] < 0.0001 &)], {i, 1, Length[tb1]}];
            classes = Union[tb1];
            identityposition = Position[classes,  n_ /; (Norm[N[n - IdentityMatrix[dimmat]], "Frobenius"] < 0.0001), {2}][[1]];
            clout=RotateLeft[classes,identityposition-1];
			Return[GTWhichOutput[clout,whinp]]
			,
	   		Print[$ErrorNoGroup];
	   		Return[]]    	
    	]

(*
***)

(****r* /GTClassMult
! NAME
!  GTClassMult
! AUTHOR
!  W. Hergert, M. Daene
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  
! USAGE
!  GTClassMult[classes,i,j] calculates the class constant C_{ijk}.
! INPUT
!  classes - list of classes
!  i       - index of class i
!  j       - index of class j
! OUTPUT
!  list
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetMatrix, GTGroupQ, GTWhichInput
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!
! LITERATURE
!
! TODO
!  Seems to work, but check the module!  
! RELEASE
!
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTClassMult[cla0_ ,ic1_ ,ic2_,OptionsPattern[]] := Module[{run,cla,grp, cl1, cl2, tb, i, j},
	
    If[Length[Quiet[Flatten[Position[{GTWhichInput[cla0[[ic1]]]},0]]]]>0,
		    		Print[$ErrorInput];
		    		Return[],None];
    If[Length[Quiet[Flatten[Position[{GTWhichInput[cla0[[ic2]]]},0]]]]>0,
	    		Print[$ErrorInput];
	    		Return[],None];

  	grp = Quiet[Flatten[cla0,1]];
  	If[OptionValue[GOFast],run=True,run=GTGroupQ[grp]];
 	If[run, 		
 				cla=GTGetMatrix[cla0];
 				cl1=cla[[ic1]];
				cl2=cla[[ic2]];
				tb=Flatten[Table[cl1[[i]].cl2[[j]],{i,1,Length[cl1]},{j,1,Length[cl2]}],1];
				Table[Length[Position[tb,cla[[i,1]]]],{i,1,Length[cla]}],
 				Print[$ErrorNoGroup];
	 			Return[]]
]

(*
***)

(****r* /GTClassMultTable
! NAME
!  GTClassMultTable
! AUTHOR
!  W. Hergert, M. Daene
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  10.2.2016 - new form of print of the table
! USAGE
!  GTClassMultTable[classes] calculates a table of class constants.
! INPUT
!  classes - list of classes
! OUTPUT
!  list
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTClassMult, GTGetSymbol, GTGroupQ
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!
! LITERATURE
!
! TODO
!  Seems to work, but check the module!  
! RELEASE
!
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

(* GTClassMultTable[cl_] :=Table[GTClassMult[cl, i, j], {i, 1, Length[cl]}, {j, 1, Length[cl]}] *)

GTClassMultTable[cl_, OptionsPattern[]] := Module[{l,i,j,k,grps,v,mt,mt1,grps1,t,a,b, group,run},
	group = Quiet[Flatten[cl,1]];
	If[OptionValue[GOFast],run=True,run=GTGroupQ[group]];
	If[run, 
 		l = Length[cl];
		If[OptionValue[GOVerbose],Do[Print["C", i, " = ", GTGetSymbol[cl[[i]]]], {i, 1, l}];,None]; 

		grps = Table[0, {l}]; 
		Do[grps[[k]] = StringJoin["C", ToString[k]], {k, 1, l}];
		t = Table[0, {l}, {l}]; 
		mt = t; 
		Do[Do[ v = GTClassMult[cl, i, j,GOFast->True]; t[[i, j]] = v;
		a = {};
		Do[b = v[[k]] StringJoin["C", ToString[k]]; 
		a = Append[a, b], {k, 1, l}]; 
		mt[[i, j]] = Plus @@ a;, {i, 1, l}], {j, 1, l}];
  		If[OptionValue[GOVerbose], 
  		  (* Print[TableForm[mt, TableHeadings -> {grps, grps}]]*)
  		  mt1 = Transpose[Join[{grps}, mt]];
          grps1 = Prepend[grps, " "]; mt1 = Join[{grps1}, mt1];
          Print[Grid[mt1, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}},
          	    Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} -> GTCornerColor}]]
          	              	                     
  		  ,None
  		];  		 
  		Return[mt],
  		Print[$ErrorNoGroup];
  		Return[]]]
	
(*
***)

(****c* /GTConjugateElement
! NAME
!  GTConjugateElement
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  
! USAGE
!  GTConjugateElement[T,X] gives the conjugate element XTX^{-1}
! INPUT
!  element1, element2
! OUTPUT
!  element
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTgmat, GTInverseElement
! GTPack NOTEBOOKS 
!   
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!    
! RELEASE
!
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTConjugateElement[elem1_,elem2_] := GTgmat[GTgmat[elem2,elem1],GTInverseElement[elem2]]

(*
***)


(****c* /GTConjugacyClass
! NAME
!  GTConjugacyClass
! AUTHOR
!  W. Hergert
! PACKAGE
!  Basic.m
! MODIFICATION HISTORY
!   August2015  : first version
!   April 2016  : transferred to Basic.m
! USAGE
!  Constructs the conjugacy class of a subgroup of a group
! INPUT
!  gr1 - group
!  gr2 - a subgroup of gr1
! OUTPUT
!  the conjugacy class
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGroupQ, GTConjugateElement
! GTPack NOTEBOOKS 
!  Reference_Symbols : GTConjugacyClass.nb , Wolfram_Book_Boxes : ConjugacyClass.nb
! DESCRIPTION
!  If GOFastValue==False the input is checked, if the two sets are really groups. The first 
!  group should be the covering group. If the order of the second group is larger, the groups
!  will be interchanged.
! LITERATURE
!  El-Batanouny p.182
! TODO
!  not fully tested
!  
! RELEASE
!  
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTConjugacyClass[gr1_, gr2_, OptionsPattern[]] :=  Module[{tt,g,cgclass,test,i},
  (*--- Check if input consists of proper groups ---*)
  If[OptionValue[GOFast]==False,
     If[GTGroupQ[gr1],
       None,
       Print[gr1, " not a group!"]; Abort[]
     ];
     If[GTGroupQ[gr2],
       None,
       Print[gr2, " not a group!"]; Abort[]
     ],
     None
  ];  
  (*--- The the second group should be the subgroup ---*) 
  If[Length[gr1] >= Length[gr2],
     None,
     tt = gr1; gr1 = gr2; gr2 = tt
   ];
  (*--- Construct of the conjugacy class ---*)
  cgclass = {};
  Do[
  	 g = gr1[[i]];
     test = Map[GTConjugateElement[#, g] &, gr2];
     cgclass = Append[cgclass, Sort[test]]
  , {i, 1, Length[gr1]}];
  cgclass = Union[cgclass];
  Return[cgclass]
]

 
(*
***)

(****c* /GTConvert
! NAME
!  GTConvert
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  08/06/2014 New Version
! USAGE
!  GTConvert[inp,i] gives the ith element of the list elm belonging to inp. 
! INPUT
!  symbol etc, index i
! OUTPUT
!  different
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack OPTIONS
!  
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTConvertMag, GTConvertNoMag, GTWhichInput
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTConvert[inp_,i_]:=If[GTmaggrpq,GTConvertMag[inp,i,grpdgrp],GTConvertNoMag[inp,i,grpdgrp]]; 

GTConvertMag[inp_,i_,grpdgrpin_] := Module[{rep0,whinp,in,strlist,mat,finalmat,mag,m,n},
  rep0 = grpdgrpin;
  Unprotect[grpdgrp];
  If[rep0 == 2,grpdgrp=5];
  If[rep0 == 3,grpdgrp=1];
  Protect[grpdgrp];
  whinp = GTWhichInput[inp];
  in = If[whinp==1,in = ToString[inp], Rationalize[N[N[inp, 15], 8]]];
  If[whinp>2,Print["Error: Quaternions, Euler angles and space group elements are not implemented for magnetic groups yet."];Abort[]];
  
  If[i==whinp,Return[inp]];
  Which[
  i == 1,Return[If[inp[[4, 4]] < 0, ToExpression[ToString[GTConvertNoMag[inp[[1 ;; 3, 1 ;; 3]],1,grpdgrp]] <> ToString["'"]],GTConvertNoMag[inp[[1 ;; 3, 1 ;; 3]],1,grpdgrp]]],
  i == 2,strlist = Characters[ToString[inp]];
         mat = GTConvertNoMag[StringJoin[If[Last@strlist == "'", mag = -1; Delete[strlist, Length[strlist]], mag = 1; strlist]],2,grpdgrp];
         finalmat = Table[0, {m, 1, 4}, {n, 1, 4}];
         finalmat[[1 ;; 3, 1 ;; 3]] = mat[[1 ;; 3, 1 ;; 3]];
         finalmat[[4, 4]] = mag;
         Return[finalmat],
  i>2,Print["Error: Quaternions, Euler angles and space group elements are not implemented for magnetic groups yet."];Abort[]
  ];
  
  
  Unprotect[grpdgrp];
  Clear[grpdgrp];
  grpdgrp = rep0;
  Protect[grpdgrp];
];

GTConvertNoMag[inp_,i_,grpdgrp_] := Module[{pos,in,index,whinp, mat},
	
 whinp = GTWhichInput[inp];
(* in = If[whinp==1,in = ToString[inp], GTSimplify[inp]];*)
 in = If[whinp==1,in = ToString[inp], Rationalize[N[N[inp, 15], 8]]];
 If[i>1,
 	If[i==whinp,Return[inp]];
 
 	Which[	grpdgrp==1,	 	pos = Position[Nelm,in];
 							If[Length[pos] > 0, index = First[pos][[1]],Print[$ErrorInputb];Abort[]]; 
 							Return[elm[[index,i]]],
 			grpdgrp==2,	 	pos = Position[NelmSU2,in];
 							If[Length[pos] > 0, index = First[pos][[1]],Print[$ErrorInputb];Abort[]]; 
 							Return[elmSU2[[index,i]]],
 		    grpdgrp==3,	 	pos = Position[Nelmo2,in];
 							If[Length[pos] > 0, index = First[pos][[1]],Print[$ErrorInputb];Abort[]]; 
 							Return[elmo2[[index,i]]],
 			grpdgrp==4,     mat = {};
 							If[i>2,Print[$ErrorOwn];Abort[]];
 							Table[If[ToString[in] == ToString[elmown[[j, 1]]],mat = elmown[[j, 2]]], {j, 1, Length[elmown]}];
 							Table[If[in == elmown[[j, 2]],mat = elmown[[j, 2]]], {j, 1, Length[elmown]}];
 							If[Length[mat]==0,Print[$ErrorInputb];Abort[],Return[mat]],
 			grpdgrp==5,	 	pos = Position[NelmSU2xS,in];
 							If[Length[pos] > 0, index = First[pos][[1]],Print[$ErrorInputb];Abort[]]; 
 							Return[elmSU2xS[[index,i]]]
 		] 	
 	,
 	Which[	grpdgrp==1,	 	pos = Position[Nelm,in];
 							If[Length[pos] > 0, index = First[pos][[1]],Print[$ErrorInputb];Abort[]]; 
 							Return[ToExpression@elm[[index,i]]],
 			grpdgrp==2,	 	pos = Position[NelmSU2,in];
 							If[Length[pos] > 0, index = First[pos][[1]],Print[$ErrorInputb];Abort[]]; 
 							Return[ToExpression@elmSU2[[index,i]]],
 			grpdgrp==3,	 	pos = Position[Nelmo2,in];
 							If[Length[pos] > 0, index = First[pos][[1]],Print[$ErrorInputb];Abort[]]; 
 							Return[ToExpression@elmo2[[index,i]]],
 		 	grpdgrp==4,	 	If[i == whinp, Return[inp]];
 		 					mat = "GTTemporaryNameXYZ67";
 							Table[If[ToString[in] == ToString[elmown[[j, 1]]],mat = elmown[[j, 1]]], {j, 1, Length[elmown]}];
 							Table[If[Norm[in-elmown[[j, 2]]]==0,mat = elmown[[j, 1]]], {j, 1, Length[elmown]}];
 							If[ToString[mat]=="GTTemporaryNameXYZ67",Print[$ErrorInputb];Abort[],Return[mat]],
 			grpdgrp==5,	 	pos = Position[NelmSU2xS,in];
 							If[Length[pos] > 0, index = First[pos][[1]],Print[$ErrorInputb];Abort[]]; 
 							Return[ToExpression@elmSU2xS[[index,i]]]
 		] 


];
 ]
 

(*
***)

(****c* /GTCyclicQ
! NAME
!  GTCyclicQ
! AUTHOR
!  W. Hergert
! PACKAGE
!  test.m -> basic.m
! MODIFICATION HISTORY
!  10/26/2012  : first version
! USAGE
!  GTCylic[group] gives True if group is a cyclic group, and gives False otherwise.
! INPUT
!  group
! OUTPUT
!  logical
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetMatrix, GTAbelianQ, GTGenerators
! GTPack NOTEBOOKS 
!  Wolfram : GTBasic.nb
! DESCRIPTION
!  Easy to construct, because a cylic group is Abelian and should have only one generator.
! LITERATURE
!
! TODO
!  not fully tested
!  documentation  
! RELEASE
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTCyclicQ[gen_] := Module[{grp}, 
	
	
(* geaendert am 10.4.2018 	grp = GTGetSymbol[gen]; *)
grp=GTGetMatrix[gen];
  If[GTAbelianQ[grp], If[Length[GTGenerators[grp]] == 1, True, False],
    False]]
(*
***) 

(****c* /GTGetInvSubGroup
! NAME
!  GTGetInvSubGroup
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   04/04/2020: first version
! USAGE
!   GTGetInvSubgroup gives an invariant subgroup  of index n
! INPUT
!   group, classes, index
! OUTPUT
!   an invariant subgroup of index n 
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTSGgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTGetIreps for the construction of a subgroup of index 2 or 3
! LITERATURE
!  -
! TODO
! -
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)


GTGetInvSubGroup[grp_, classes_, order_, OptionsPattern[]] := Module[{lg, pos, reducedclasses, reducedclassmult, lcl, numb, ln, 
   orderpos, lorderpos, nn, llfin, eq, eqint, const, opt, erg, poscl, C, it, cc, case, sets, choice, compare, grpout, groupfound},
  lg = Length[grp];
  (*remove the class containing the identity element*)
  (*classes = GTClasses[grp,GOFast->True];*)
  
  pos = If[grpdgrp ==4, {Position[classes,GTGetSymbol[IdentityMatrix[Length[GTGetMatrix[classes[[1, 1]]]]]]], 
     Position[classes, IdentityMatrix[Length[GTGetMatrix[classes[[1, 1]]]]]]}, {Position[classes, Ee], 
     Position[classes, IdentityMatrix[Length[GTGetMatrix[classes[[1, 1]]]]]], Position[classes, GTGetEulerAngles[Ee]], 
     Position[classes, GTGetQuaternion[Ee]]}];
  
  reducedclasses = Delete[classes, pos[[First@First@Position[Map[Length, pos], 1], 1, 1]]];
  (*estimate the class lengths, positions of classes with a specified length, and numbers of classes with specified length*)
  
  lcl = Table[Length[reducedclasses[[i]]],{i,1,Length[reducedclasses]}];(*Map[Length, reducedclasses];*)
  numb = Union[lcl];
  ln = Length[numb];
  orderpos = Table[Position[lcl, numb[[i]]], {i, 1, ln}];
  lorderpos = Map[Length, orderpos];
  nn = Table[Length[Position[lcl, numb[[i]]]], {i, 1, ln}];
  llfin = lg/order;
  (*the invariant subgroup is a union of classes. Set up an equation to estimate possible combinations*)
  
  eq = {1 + Table[C[i], {i, 1, ln}].numb == llfin};
  eqint = Table[Element[C[i], Integers], {i, 1, ln}];
  const = Table[0 <= C[i] <= lorderpos[[i]], {i, 1, ln}];
  opt = Solve[Flatten[{eq, eqint, const}], Table[C[i], {i, 1, ln}]];
  If[Length[opt] == 0, If[OptionValue[GOVerbose], Print["Warning: GTGetInvSubgroupOrdern could not find a normal subgroup of specified order."]]; Return[{}]];
  (*calculate the class multiplication table for the reduced classes (without identity)*)
  reducedclassmult = Table[
    erg = {};
    Do[
     poscl = Position[reducedclasses, GTgmat[reducedclasses[[i, m]], reducedclasses[[j, n]]]];
     If[Length[poscl] > 0, erg = Append[erg, First@First@poscl]];
     , {m, 1, Length[reducedclasses[[i]]]}, {n, 1, Length[reducedclasses[[j]]]}];
    Union[erg]
    , {i, 1, Length[reducedclasses]}
    , {j, 1, Length[reducedclasses]}];
  (*Try to find an invariant subgroup using random choice of possible combinations.*)
  it = 0;
  groupfound = False;
  While[Not[groupfound] && it < OptionValue[MaxIterations],
   it = it + 1;
   cc = RandomInteger[{1, Length[opt]}];
   case = Table[C[j], {j, 1, ln}] /. opt[[cc]];
   sets = Table[Subsets[orderpos[[j]], {case[[j]]}], {j, 1, ln}];
   choice = {};
   Do[If[case[[i]] > 0, choice = Append[choice, sets[[i, RandomInteger[{1, Length[sets[[i]]]}]]]]], {i, 1, ln}];
   choice = Flatten[choice];
   compare = {};
   Do[compare = Append[compare, reducedclassmult[[choice[[i]], choice[[j]]]]], {i, 1, Length[choice]}, {j, 1, Length[choice]}];
   compare = Union[Flatten[compare]];
   groupfound = Length[Union[compare, choice]] == Length[choice]
   ];
  If[Not[groupfound], If[OptionValue[GOVerbose], Print["Warning: GTGetInvSubgroupOrdern could not find a normal subgroup of specified order. In case you know a normal subgroup of specified order exists, try to increase MaxIterations."]]; Return[{}]];
  (*Construct the final output*)
  
  grpout = classes[[pos[[First@First@Position[Map[Length, pos], 1], 1, 1]]]];
  Do[
   Do[
    grpout = Append[grpout, reducedclasses[[choice[[i]], j]]]
    , {j, 1, Length[reducedclasses[[choice[[i]] ]] ]}]
   , {i, 1, Length[choice]}];
  Return[grpout]]
(*
***)

(****c* /GTGenerators
! NAME
!  GTGenerators
! AUTHOR
!  W. Hergert
! PACKAGE
!  test.m    -> basic.m
! MODIFICATION HISTORY
!  10/26/2012  : first version
!  01/01/2013: Wolfram: input as matrices symbols or quaternions. internaly \SmallCircel instead \CirclePlus
! USAGE
!  GTGenerators[group] gives the generators of a certain group.
! INPUT
!  group
! OUTPUT
!  set of generators
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGroupQ, GTWhichInput, GTOrderOfElement, GTMultTable, GTWhichOutput
! GTPack NOTEBOOKS 
!  Wolfram : GT_Basic.nb
! GTPack OPTIONS
!	GOFast
!	
! DESCRIPTION
!  To know the generators of a group makes it simpler to reinstall the group.
!
! LITERATURE
!
! TODO
!  not fully tested
!  documentation  
! RELEASE
!
! PROBLEMS
!  for Oh it gives wrong result, we get generators for O
! 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTGenerators::badarg = 
  "Error: Argument is not a group.";

GTGenerators[g0_,OptionsPattern[]] := Module[
    {ord, gs, gen, oe, t, co, ig, o1, i,g,whinp,subs,sl,log,final},
    If[OptionValue[GOFast]==False,
       If[GTGroupQ[g0],
  	 	  None, 
  	 	  Message[GTGenerators::badarg];Abort[]
  	   ],  							
  	   None
    ];
  (*---check type of input and convert to matrix---*)
    whinp = GTWhichInput[g0[[1]]];
    g     = GTGetMatrix[g0];
    ord   = Length[g];
  (*--- Search for generators ---*)
    gs    = SortBy[g, GTOrderOfElement[#] &]; 
    gen   = {gs[[ord]]};
    oe    = GTOrderOfElement[gen[[1]]]; 
    t     = Table[0, {oe}]; 
    t[[1]]= gen[[1]]; 
    Do[
       t[[i + 1]] = t[[i]]\[SmallCircle]t[[1]]
    , {i, 1, oe - 1}];
    co = Complement[g, t]; 
    ig = 1; 
    While[Length[co] > 0, 
  	      ig  = ig + 1; 
  	      o1  = Length[co]; 
          gen = Append[gen,co[[o1]]];
          t   = Union[t, Map[gen[[ig]]\[SmallCircle]# &, t]]; 
          co  = Complement[g, t];
    ]; 
  (*---minimize number of generators---*)  
    ord  = Length[gen];
    subs = Subsets[gen,{1,ord}];
    sl   = Length[subs];
    log  = True;
    i    = 1;
    While[log,
          t =GTGroupFromGenerators[subs[[i]]];
  	      If[Complement[g,t]=={},
  	         final=subs[[i]];log=False,
  	         If[i<sl,
  	         	i=i+1,
  	         	log=False
  	         ]
  	      ] 
    ];
    Return[GTWhichOutput[final,whinp]]   	
  ]
  
(*
GTGenerators[g0_,OptionsPattern[]] := Module[{ord, gs, gen, oe, t, co, ig, o1, tb, i,g,whinp},
  If[OptionValue[GOFast]==False,
  							If[GTGroupQ[g0],None,Print[$ErrorNoGroup];Abort[]],  							
  							None];
  whinp = GTWhichInput[g0[[1]]];
  g=GTGetMatrix[g0];
  ord = Length[g];
  (*--- Search Generators ---*)
  gs = SortBy[g, GTOrderOfElement[#] &]; gen = {gs[[ord]]};
  oe = GTOrderOfElement[gen[[1]]]; t = Table[0, {oe}]; 
  t[[1]] = gen[[1]]; 
  Do[t[[i + 1]] = t[[i]]\[SmallCircle]t[[1]], {i, 1, oe - 1}];
  co = Complement[g, t]; 
  ig = 1; 
  While[Length[co] > 0, 
  	ig = ig + 1; o1 = Length[co]; 
    gen = Append[gen, SortBy[g, GTOrderOfElement[#] &][[o1]]];
    t = Union[t, Map[gen[[ig]]\[SmallCircle]# &, t]]; 
    co = Complement[g, t]]; 
    tb = GTMultTable[gen,GOVerbose-> False,GOFast->True]; 
    gen = Complement[gen, Intersection[Flatten[tb], gen]
  ];
  Return[GTWhichOutput[gen,whinp]]    	
  ]

(*
***) 
*)


(****c* /GTGetEulerAngles
! NAME
!  GTGetEulerAngles
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  15/05/2013 Change of Structure, Error Messages
! USAGE
!  GTGetEulerAngles[symmetry element] gives the Euler angles of a certain symmetry element.
! INPUT
!  symmetry element
! OUTPUT
!  Euler Angles
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack OPTIONS
!  GOAbort
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTWhichInput, GTConvert
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  error message, if input is spacegroup element  
! RELEASE
!
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)

GTGetEulerAngles[inp_] := 
  Which[Length[GTWhichInput[inp]]==0, GTConvert[inp,4], 
  	    Length[GTWhichInput[inp]]>0, Map[GTGetEulerAngles, inp]]


 
(*
***)


(****c* /GTGetMatrix
! NAME
!  GTGetMatrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  15/05/2013 Change of Structure, Error Messages
!  21/02/2014 new version
! USAGE
!  GTGetmatrix[symmetry element] gives the matrix representation for a given symmetry element.
! INPUT
!  symmetry element
! OUTPUT
!  Matrix
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack OPTIONS
!  
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTWhichInput, GTConvert
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)
GTGetMatrix[inp_] := Which[Length[GTWhichInput[inp]]==0, GTConvert[inp,2], 
  	                       Length[GTWhichInput[inp]]>0, Map[GTGetMatrix, inp]] 


(*
***)

(****c* /GTGetQuaternion
! NAME
!  GTGetQuaternion
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  15/05/2013 Change of Structure, Error Messages
!  24/02/2014 New design
! USAGE
!  GTGetQuaternion[symmetry element] gives the quaternion for a given symmetry element.
! INPUT
!  symmetry element
! OUTPUT
!  Quaternion
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTWhichInput, GTConvert
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  error message if input is space group element  
! RELEASE
!
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)


GTGetQuaternion[inp_] := Which[Length[GTWhichInput[inp]]==0, GTConvert[inp,3], 
  	                       Length[GTWhichInput[inp]]>0, Map[GTGetQuaternion, inp]] 

(*
***)

(****c* /GTGetRotationMatrix
! NAME
!  GTGetRotationMatrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  08/06/2014 New Version
!  23/08/2016 Rotation matrices for space group elements
! USAGE
!  GTGetRotationMatrix[symmetry element] gives the 3-dimensional rotation matrix of a certain symmetry element.
! INPUT
!  symmetry element
! OUTPUT
!  rotation matrix
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTWhichInput, GTGetSymbol, GTGetMatrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTGetRotationMatrix[inp_]:= Module[{mat,rep,sym},
	sym=If[GTWhichInput[inp]==5,GTGetSymbol[inp[[1]]],GTGetSymbol[inp]];
	(*matdim=If[grpdgrp==3,3,1];*)
	rep = grpdgrp;
	Unprotect[grpdgrp];
	grpdgrp=1;
	mat = GTGetMatrix[sym][[1;;3,1;;3]];
	grpdgrp=rep;
	Protect[grpdgrp];
	Return[mat]
 ]

(*
***)

(****c* /GTGetSubGroups
! NAME
!  GTGetSubGroups
! AUTHOR
!  W. Hergert
! PACKAGE
!  test.m -> basic.m
! MODIFICATION HISTORY
!  10/26/2012  : first version
!  05/03/2013  : options introduced 
! USAGE
!  GTGetSubGroups[group] gives all subgroups of a group.
! INPUT
!  group     -  group to determine subgroups
! OUTPUT
!  list of subgroups
! GTPack OPTIONS
!  GOGroupOrder, GOFast 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetSymbol, GTGroupQ
! GTPack NOTEBOOKS 
!  Wolfram : GTBasic.nb
! DESCRIPTION
!  The algorithm is based on the Euler-Lagrange theorem. This gives the possible orders of subgroups. 
!  All subsets of the according lengths are checked if it is a group.
! LITERATURE
!
! TODO
!  not fully tested
!  documentation  
! RELEASE
!
! PROBLEMS
!  The algorithm becomes of course slow for groups of higher order. Up to ord(G)=12 it works well.
!  It helps to search only subgroups of a certain order.
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTGetSubGroups[g0_,OptionsPattern[]] := Module[{grp, n, div, sub, n1, sg, test, lt, j, g}, 
  If[!OptionValue[GOFast],
  		   	If[!GTGroupQ[g0],Print[$ErrorNoGroup];Abort[]];
  		   	If[OptionValue[GOGroupOrder]>Length[g0],Print[$ErrorGTGetSubgroupsOrder];Abort[]];
    ];	
  g=GTGetSymbol[g0];
  grp=GTGetSymbol[g0];
  n = Length[g]; 
  div = Divisors[n]; 
  n1 = Length[div];nsu=OptionValue[GOGroupOrder];
  sg = {}; 
  If[nsu === 1, 
   Do[sub = Subsets[Complement[grp, {Ee}], {div[[i]] - 1}] ; 
    lt = Length[sub]; 
    Do[test = Union[sub[[j]], {Ee}]; 
     If[GTGroupQ[test], sg = Append[sg, test] , Null], {j, 1, 
      lt}], {i, 2, n1 - 1}], Null];
  If[nsu > 1, 
   If[Intersection[div, {nsu}] === {}, 
    Print["According to Euler-Lagrange theorem no subgroup"], 
    sub = Subsets[Complement[grp, {Ee}], {nsu - 1}] ; 
    lt = Length[sub]; 
    Do[test = Union[sub[[j]], {Ee}]; 
     If[GTGroupQ[test], sg = Append[sg, test] , Null], {j, 1, lt}]]]; 
     Return[sg]  
  ]

(*
***)

(****c* /GTGetSU2Matrix
! NAME
!  GTGetSU2Matrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  08/06/2014 New Version
! USAGE
!  GTGetSU2Matrix[symmetry element] gives the SU(2) matrix of a certain symmetry element.
! INPUT
!  symmetry element
! OUTPUT
!  SU2Matrix
! ERROR MESSAGES
!  $ErrorInputb
! GTPack OPTIONS
!  
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetMatrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTGetSU2Matrix[inp_]:= Module[{mat,rep},
	rep = grpdgrp;
	Unprotect[grpdgrp];
	grpdgrp=5;
	mat = GTGetMatrix[inp][[1;;2,1;;2]];
	Unprotect[grpdgrp];
	grpdgrp=rep;
	Protect[grpdgrp];
	Return[mat]
 ]

(*
***)

(****c* /GTGetSymbol
! NAME
!  GTGetSymbol
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  15/05/2013 Change of Structure, Error Messages
!  02/21/2014 New Version
! USAGE
!  GTGetSymbol[symmetry element] gives the symbol of a certain symmetry element.
! INPUT
!  symmetry element
! OUTPUT
!  Symbol
! ERROR MESSAGES
!  $ErrorInput
!  $ErrorNoSymbol
!  Error: Symbols can not be achieved from space group elements!
! GTPack OPTIONS
!  
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTWhichInput, GTGetSymbol, GTConvert
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTGetSymbol[inp_] :=	Which[	Length[GTWhichInput[inp]]==0,
									GTConvert[inp,1], 
									
  								Length[GTWhichInput[inp]]>0, 	
  									Map[GTGetSymbol, inp]]
 

(*
***)

(****)(****c* /GTgmat
! NAME
!  GTgmat
! AUTHOR
!  W. Hergert, M.Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  12/29/2012 : new definition of group operation, new symbol for the operation.
!               We use now \SmallCircle instead of \CirclePlus
!  05/03/2013 : OptionPattern added
!  07/06/2013 : Rewritten using GTWhichInput
!  02/24/2014 : Rename to GTgmat
! USAGE
!  GTgmat[element1,element2] performs a multiplication of element1 and element2.
! INPUT
!  element1,element2 are pointgroup or spacegroup elements. Every input form is allowed!
! OUTPUT
!  If cg1 and cg2 are both pointgroup elements the result is a pointgroup element. If one of 
!  the elements is a spacegroup element the result is a spacegroup element.
! ERROR MESSAGES
!  $ErrorInput
!  Warning: Arguments are not point group operations.
! GTPack OPTIONS
!  GOVerbose
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetMatrix, GTWhichInput, GTWhichOutput
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  The formulation is general. Both point group elements as well as space group elements can be used. 
! 
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTgmat[cg1_,cg2_,OptionsPattern[]]:= Module[{h1,h2,c1,c2,c1mat,c2mat,result,c1t,c2t,materg},
	(*--- check the input ---*)
	h1 = GTWhichInput[cg1];
	h2 = GTWhichInput[cg2];
	
	(*--- error messages for invalid input ---*)
	If[h1==0,Print[$ErrorInput];Abort[],None];
	If[h2==0,Print[$ErrorInput];Abort[],None];
	
	(*--- Multiplication ---*)
	Which[	h1<5,
				Which[	h2<5,
							result=GTGetMatrix[cg1].GTGetMatrix[cg2];
							Return[GTWhichOutput[result,h1]],
						h2==5,
							If[OptionValue[GOVerbose],Print["Warning: Arguments are not point group operations."],None];
							c1=GTType[cg1];
							c2=cg2;
							c1mat=GTGetMatrix[c1[[1]]];
							c2mat=GTGetMatrix[c2[[1]]];
							materg=GTWhichOutput[c1mat.c2mat,GTWhichInput[c1[[1]]]];
							c1t=c1[[2]];
							c2t=c2[[2]];
							result=AngleBracket[materg,c1mat.c2t+c1t]];
							Return[result],
			h1==5,
				Which[	h2<5,
							If[OptionValue[GOVerbose],Print["Warning: Arguments are not point group operations."],None];
							c1=cg1;
							c2=GTType[cg2];
							c1mat=GTGetMatrix[c1[[1]]];
							c2mat=GTGetMatrix[c2[[1]]];
							materg=GTWhichOutput[c1mat.c2mat,GTWhichInput[c1[[1]]]];
							c1t=c1[[2]];
							c2t=c2[[2]];
							result=AngleBracket[materg,c1mat.c2t+c1t];
							Return[result],
						h2==5,
							If[OptionValue[GOVerbose],Print["Warning: Arguments are not point group operations."],None];
							c1=cg1;
							c2=cg2;
							c1mat=GTGetMatrix[c1[[1]]];
							c2mat=GTGetMatrix[c2[[1]]];
							materg=GTWhichOutput[c1mat.c2mat,GTWhichInput[c1[[1]]]];
							c1t=c1[[2]];
							c2t=c2[[2]];
							result=AngleBracket[materg,c1mat.c2t+c1t]];
							Return[result]];
]

HoldPattern[cg1_ \[SmallCircle] cg2_] := GTgmat[cg1,cg2]

(*
***)




(****c* /GTGroupOrder
! NAME
!  GTGroupOrder
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  12/31/2012 : Wolfram: Checks first if the list is a group
!  01/01/2013 : Wolfram: argument mode to suppress printout
!  05/03/2013 :Options introduced
! USAGE
!  GTGroupOrder[group] gives the order of a group.
! INPUT
!  group
! OUTPUT
!  integer
! ERROR MESSAGES
!  $ErrorNoGroup
! GTPack OPTONS
!  GOVerbose : supression of printout
! STANDARD OTIONS
!  
! GTPack MODULES
!   GTGroupQ
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTGroupOrder[group_,OptionsPattern[]]:=Module[{ord, run},
	If[OptionValue[GOFast],run=True,run=GTGroupQ[group]];
	If[run,
		ord=Length[group];
		If[OptionValue[GOVerbose],Print["Order of the group: ",ord],Null];
		Return[ord],
	
		Print[$ErrorNoGroup];
		Return[]]
]
(*
***)

(****c* /GTGroupQ
! NAME
!  GTGroupQ
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  10/25/2012 : Wolfram : print of multiplication table suppressed
!  12/30/2012 : Wolfram : excepts matrices, symbols and quaternions now
! USAGE
!  GTGroupQ[set of symmetry elements] gives True if set of symmetry elements forms a group, and gives False otherwise.
! INPUT
!  set of symmetry elements
! OUTPUT
!  logical
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTWhichInput, GTGetMatrix, GTSimplify
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  We have to work with the matrix representation. Using Symbols, Errors occur for angular momentum representations, using GTAngularMomentumRep.
!
!--------------------------------------------------------------------------------
! SOURCE
*)

 GTGroupQ[group_,OptionsPattern[]]:=Module[{grpt,lg,whinp,tb},
 	whinp = GTWhichInput[group];
 	If[Length[Position[whinp,0]]==0,
 		(*Numerical*)
        If[Characters[OptionValue[GOMethod]][[1]] == "N",
        	grpt = Union[N[GTGetMatrix[group]],SameTest -> (Norm[#1 - #2, "Frobenius"] < 0.00001 &)];
            lg = Length[grpt];
            tb = Flatten[Table[grpt[[i]].grpt[[j]], {i, 1, lg}, {j, 1, lg}], 1];
            Length[Complement[tb, grpt, SameTest -> (Norm[#1 - #2, "Frobenius"] < 0.00001 &)]] == 0,
        (*Analytical*)
            grpt = Union@GTSimplify@GTGetMatrix[group];
            lg = Length[grpt];
            tb = Union[GTSimplify@Simplify@Flatten[Table[grpt[[i]].grpt[[j]], {i, 1, lg}, {j, 1, lg}], 1]];
            Length[Complement[tb, grpt]] == 0],			
	Print[$ErrorInput];Return[False]]
]
(*
***)


(****c* /GTInverseElement
! NAME
!  GTInverseElement
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  01/01/2013 : Wolfram : extended to space group elements
!  April      : Wolfram : error for inverse space group element corrected 
!  March 2019 : Matthias : lattice for space group elements added 
! USAGE
!  GTInverseElement[symmetry element] gives the inverse element of a certain symmetry element.
! INPUT
!  symmetry element
! OUTPUT
!  symmetry element
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTWhichInput, GTGetMatrix, GTWhichOutput
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTInverseElement[inp_,OptionsPattern[]]:=Module[{whinp,inv,c1,c2,c1ind,invvec,latt,lt,list,sol,dc,invel},
  latt=OptionValue[GOLattice];
  whinp = GTWhichInput[inp];
  If[whinp==0,Print[$ErrorInput];Abort[],None];
  Which[	whinp<5, 
				inv = Inverse[GTGetMatrix[inp]];
				Return[GTWhichOutput[inv,whinp]],
  			whinp==5,
  				c1 = inp[[1]];
  				c2 = inp[[2]].latt;
  				c1ind = GTWhichInput[c1];
  			
  				inv = Inverse[GTGetRotationMatrix[c1]];
  				invel = Inverse[GTGetMatrix[c1]];
  				invvec = -inv.c2;
  				
  				lt = Length[invvec];
                  list = Sum[c[i] latt[[i]], {i, 1, lt}];
                  sol = Solve[Table[list[[i]] == invvec[[i]], {i, 1, lt}],Table[c[i], {i, 1, lt}]];
                  dc = Table[c[i], {i, 1, lt}] /. sol[[1]];
  				
  				Return[AngleBracket[GTWhichOutput[invel,c1ind],dc]]]
]

(*
***)

(****c* /GTInvSubGroups
! NAME
!  GTInvSubGroups
! AUTHOR
!  W. Hergert
! PACKAGE
!  basic.m
! MODIFICATION HISTORY
!  13/11/2014 : moved to basic.m
!  12/07/2012 : first version
! USAGE
!  GTInvSubGroups[group] gives the invariant subgroups of a group.
! INPUT
!  group
! OUTPUT
!  list of all invariant subgroups
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTWhichInput, GTGetMatrix, GTGroupQ, GTClasses
! GTPack NOTEBOOKS 
!  GTClasses
!  GTGroupQ
! DESCRIPTION
!  The algorithm is based on the fact, that an invariant subgroup consists
!  of complete classes. All possible combinations of classes are calculated.
!  Then it is tested if it is a group. If this is true, the group has to 
!  be an invariant subgroup
! LITERATURE
!   
! TODO
!  Write the documentation notebooks.
!  Seems to work, but check the module!  
! RELEASE
!
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTInvSubGroups[gr0_,OptionsPattern[]] := 
 Module[{lc, g, lg, ltg, n, lt, gi, gr, dim, whinp},
  (*--- Check the input ---*)
  If[OptionValue[GOFast],None,If[GTGroupQ[gr0],None,Print[$ErrorNoGroup];Abort[]]];
  whinp = GTWhichInput[gr0[[1]]];
  gr = GTGetMatrix[gr0];
  dim = Length[gr[[1]]];
  lc = GTClasses[gr];
  g = Length[gr];
  lc = Complement[lc, {{IdentityMatrix[dim]}}];
  lg = Subsets[lc];
  ltg = {};
  n = Length[lg];
  Do[lt = Flatten[Union[lg[[i]], {{IdentityMatrix[dim]}}], 1];
   gi = Length[lt];
   If[GTGroupQ[lt], 
    If[gi > 1 && gi < g, ltg = Append[ltg, GTWhichOutput[lt, whinp]], 
     Null], Null], {i, 1, n}]; Return[ltg]]

(*
***)

(****c* /GTInvSubGroupQ
! NAME
!  GTInvSubGroupQ
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  12/31/2012 : Wolfram : arbitrary order of elements, works with matrices, quaternions and symbols
!  07/06/2013 : Matthias : New module
! USAGE
!  GTInvSubSubGroupQ[group1,group2] gives True if the group with smaller order is an invariant subgroup of the group with larger order, and gives False otherwise.
! INPUT
!  group1
!  group2
! OUTPUT
!  logical
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGroupQ, GTSubGroupQ, GTGetMatrix, GTConjugateElement
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTInvSubGroupQ[group1_,group2_]:=Module[{l1,l2,group,subgroup,i1,i2,invsbgrp,tb,utb},
	
	invsbgrp = True;
	
	If[GTGroupQ[group1],None,Print[$ErrorNoGroup," Input 1"];Abort[]];
	If[GTGroupQ[group2],None,Print[$ErrorNoGroup," Input 2"];Abort[]];

	l1 = Length[group1];
	l2 = Length[group2];
	If[l1>l2,
			group=GTGetMatrix[group1];
			subgroup=GTGetMatrix[group2];
			i1=1;
			i2=2,	
			
			group=GTGetMatrix[group2];
			subgroup=GTGetMatrix[group1];
			i1=2;
			i2=1];

	If[GTSubGroupQ[group,subgroup],None,Print["Error: Argument ",i2," is not a subgroup of argument ",i1];Return[False]];
	tb=Table[GTConjugateElement[subgroup[[j]],group[[i]]],{i,1,Length[group]},{j,1,Length[subgroup]}];
	utb=Flatten[Union[tb],1];
	If[Length[Union[utb,subgroup]]==Length[subgroup],True,False]
]

(*
***)

(****c* /GTLeftCosets
! NAME
!  GTLeftCosets
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  06/11/2020 : Matthias: module simplified
!  11/25/2014 : Matthias: GOFast, optimization of the module
!  01/01/2013 : Wolfram : order of group and subgroup arbitrary, form of input symbol, matrix, quaternion
! USAGE
!  LeftCoset[group,subgroup] gives the left coset of a group according to a subgroup.
! INPUT
!  group, subgroup
! OUTPUT
!  list of element
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTGroupQ, GTWhichInput, GTGetMatrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTLeftCosets[group_, sgroup_, OptionsPattern[]] := Module[{el, pos, lc,grp1,grp2,whinp},
  If[Length[group]>Length[sgroup],
    grp1 = GTGetMatrix[group];
    grp2 = GTGetMatrix[sgroup];,
    grp2 = GTGetMatrix[group];
    grp1 = GTGetMatrix[sgroup];];
  whinp = GTWhichInput[group[[1]]];

  If[OptionValue[GOFast],None,
		If[GTGroupQ[grp1],None,Print[$ErrorNoGroupFirst];Abort[]];
		If[GTGroupQ[grp2],None,Print[$ErrorNoGroupSecond];Abort[]];
		If[GTSubGroupQ[grp1,grp2],None,Print["Error: input should contain a group and a subgroup."]];
	];	  
 
  lc = Union[Map[Sort[#] &, 
  	      Table[
  	      	el = GTgmat[grp1[[i]], grp2[[j]]];
  	      	pos = Position[grp1, el];
  	      	grp1[[pos[[1, 1]]]], {i, 1, Length[grp1]}, {j, 1, Length[grp2]}]]];
  GTWhichOutput[Table[SortBy[lc[[i]], GTOrderOfElement[#] &], {i, 1, Length[lc]}], whinp]]

(*
***)


(****c* /GTMagnetic
! NAME
!  GTMagnetic
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  06/05/2016 new command 
! USAGE
!  GTMagnetic[logical] sets the variable GTmaggrpq and determines if magnetic point group have to be taken into account.
! INPUT
!  logical
! OUTPUT
! 
! ERROR MESSAGES
!  
! GTPack OPTIONS
!  
! STANDARD OPTIONS
!
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)
  GTMagnetic[log_] := Module[{},
    Unprotect[GTmaggrpq];
    GTmaggrpq = log;
    Protect[GTmaggrpq];
  ]
(*
***)

(****c* /GTMagneticQ
! NAME
!  GTMagneticQ
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  06/05/2016 new command 
! USAGE
!  GTMagneticQ[] gives True if if magnetic groups are taken into account, and gives False otherwise.
! INPUT
!  logical
! OUTPUT
! 
! ERROR MESSAGES
!  
! GTPack OPTIONS
!  
! STANDARD OPTIONS
!
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)
  GTMagneticQ[] := GTmaggrpq
  
(*
***)


(****c* /GTMultTable
! NAME
!  GTMultTable
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  12/30/2012 : multiplication changed to \SmallCircle
!  05/03/2013 : options introduced
   02/10/2016 : form of table changed (W)
! USAGE
!  GTMultTable[group] gives the multiplication table of a group.
! INPUT
!  group
! OUTPUT
!  multiplication table
! ERROR MESSAGES
!  $ErrorNoGroup
! GTPack OPTIONS
!  GOVerbose : printout of multiplication table if TRUE
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetSymbol, GTGroupQ, GTWhichInput, GTMultTableInternal, GTGetMatrix, GTWhichOutput, GTgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)


GTMultTable[group0_,OptionsPattern[]] :=Module[{mt,group,mtprint,grps,whinp,mtout,run,mt1,grps1},
	If[OptionValue[GOFast],run=True,run=GTGroupQ[group0]];
	If[run,
		whinp=GTWhichInput[group0[[1]]];
		If[whinp<5,group=GTGetMatrix[group0],None];
		mt=GTMultTableInternal[group,whinp];
		If[OptionValue[GOVerbose],
			mtprint=Quiet[GTGetSymbol[mt]];
			grps=Quiet[GTGetSymbol[group]];
			(* Print[TableForm[mtprint,TableHeadings->{grps,grps}]]*)
			mt1 = Transpose[Join[{grps}, mtprint]];
            grps1 = Prepend[grps, " "]; mt1 = Join[{grps1}, mt1];
            Print[Grid[mt1, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}},
            	  Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} -> GTCornerColor}
            	  ]]
		    ,None
		 ],
		 Print[$ErrorNoGroup];
		 Abort[]];
	mtout = GTWhichOutput[mt,whinp];	
	Return[mtout]]

GTMultTableInternal[group_,whinp_] := Module[{lg},
		lg=Length[group];
		If[whinp<5,
			Table[group[[i]].group[[j]],{i,1,lg},{j,1,lg}],		
			
			Simplify[Table[GTgmat[group[[i]],group[[j]]],{i,1,Length[group]},{j,1,Length[group]}]]
		]
]

(*
***)


(****c* /GTNormalizer
! NAME
!  GTNormalizer
! AUTHOR
!  W. Hergert
! PACKAGE
!  Basic.m
! MODIFICATION HISTORY
!   August2015  : first version
!   April 2016  : transferred to Basic.m
! USAGE
!  Constructs the nomalizer from a group and its supergroup
! INPUT
!  gr1 - supergroup, covering group
!  gr2 - a subgroup of gr1
! OUTPUT
!  the normalizer
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! ERROR MESSAGES
!   Messages if GOFastValue==False and input not a group
! GTPack MODULES
!  GTGroupQ, GTConjugateElement
! GTPack NOTEBOOKS 
!  Reference_Symbols : GTNormalizer.nb , Wolfram_Book_Boxes : Normalizer.nb
! DESCRIPTION
!  If GOFastValue==False the input is checked, if the two sets are really groups. The first 
!  group should be the covering group. If the order of the second group is larger, the groups
!  will be interchanged.
! LITERATURE
!  El-Batanouny p.185
! TODO
!  not fully tested
!  
! RELEASE
!  
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)


GTNormalizer[gr1_, gr2_, OptionsPattern[]] := Module[{tt,norm,g,test,i},
  (*--- Check if input consists of proper groups ---*)
  If[OptionValue[GOFast]==False,
     If[GTGroupQ[gr1],
        None,
        Print[gr1, " not a group!"]; Abort[]
     ];
     If[GTGroupQ[gr2],
        None,
        Print[gr2, " not a group!"]; Abort[]
     ],
     None
  ];   
  (*--- The first group should be the supergroup ---*)
  If[Length[gr1] >= Length[gr2],
     None,
     tt = gr1; gr1 = gr2; gr2 = tt
  ];
  (*--- Construct the normalizer ---*)
  norm = {};
  Do[
  	 g = gr1[[i]];
     test = Map[GTConjugateElement[#, g] &, gr2];
     If[Sort[test] === Sort[gr2],
        norm = Append[norm, g],
        None
     ];
  ,{i, 1, Length[gr1]}];
  Return[norm]
]
  
  
(*
***)


(****c* /GTOrderOfElement
! NAME
!  GTOrderOfElement
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  
! USAGE
!  GTOrderOfElement[element] gives the order of a group element.
! INPUT
!  element
! OUTPUT
!  integer
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetMatrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTOrderOfElement[symb_] := Module[{ord, mat, idm, mat1},
 ord = 1;
 mat = N[GTGetMatrix[symb]];
 idm=N[IdentityMatrix[Length[mat[[1]]]]];
 mat1 = mat;
 While[!idm == mat1, mat1 = Chop[mat.mat1]; ord = ord + 1];
 Return[ord]]

(*
***)

(****c* /GTO2ToSymbolInfo
! NAME
!  GTO3ToSymbolInfo
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2013 first commit
! USAGE
!  GTO2ToSymbolInfo[matrix] gives the following information from a O(3) matrix: {axis,angle,multiplicity,determinant}
! INPUT
!  matrix
! OUTPUT
!  {axis,angle,multiplicity,determinant}
! ERROR MESSAGES
!
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTO3ToSymbolInfo
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)

GTO2ToSymbolInfo[mat_] := GTO3ToSymbolInfo[{{mat[[1,1]],mat[[1,2]],0},{mat[[2,1]],mat[[2,2]],0},{0,0,1}}]

(*
***)


(****c* /GTO3ToSymbolInfo
! NAME
!  GTO3ToSymbolInfo
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2013 first commit
! USAGE
!  GTO3ToSymbolInfo[matrix] gives the following information from a O(3) matrix: {axis,angle,multiplicity,determinant}
! INPUT
!  matrix
! OUTPUT
!  {axis,angle,multiplicity,determinant}
! ERROR MESSAGES
!
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
! 
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)

GTO3ToSymbolInfo[matinp_] := 
 Module[{det, es, eval, evpos, ev, angtmp, ang, mult, mat},
  (*--- calculate the determinant (proper/improper?) ---*)
  
  det = Det[matinp // N] // Rationalize;
  mat = det*matinp;
  (*--- All information can be achieved from the eigensystem of the \
rotation matrix ---*)
  es = Eigensystem[mat];
  eval = N[{Log[es[[1]]]/I}/(2 \[Pi])] // Chop // Rationalize;
  (*--- determine the eigenvector ---*)
  
  evpos = Flatten[Position[eval, 0]][[2]];
  ev = Normalize[es[[2, evpos]]] // Simplify;
  pos = Position[axeslist, ev];
  fac = 1;
  If[Length[pos] == 0, ev = -ev; fac = -1];
  (*--- determine the multiplicity of the rotation ---*)
  
  angtmp = fac*Max[eval];
  {ang, mult} = {Denominator[angtmp], Abs[Numerator[angtmp]]};
  (*--- Output ---*)
  
  If[RotationMatrix[(ang - 1)*2*\[Pi]*angtmp, ev] === mat, {ev, 
    fac*ang, mult, det}, {ev, -fac*ang, mult, det}]
  ]
(*
***)


(****c* /GTProductGroup
! NAME
!  GTProductGroup
! AUTHOR
!  W. Hergert
! PACKAGE
!  Basic.m
! MODIFICATION HISTORY
!   July 2015   : first version
!   April 2016  : transferred to Basic.m
!   July 2017   : removed mistake with identifying semi-direct product groups; shortened estimation of commutation relation
! USAGE
!  Constructs the product group from two groups
! INPUT
!  grp1 - a group
!  grp2 - a group
! OPTIONS
!  GOSymbols - True : Output of the result in symbolic form
!              False: Output in matrix form
! OUTPUT
!  if construction according tot he definitions the product group is given.
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTWhichInput,GTGroupQ, GTGetMatrix, GTSimplify, GTInvSubGroupQ, GTSetMultiplication
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  It is not checked, if a group can be written as a product of two groups.
!  The command works the other way arround. Two groups are given and it is tested,
!  if a product group can be formed. First it is tested, if the groups have only
!  Ee in common. If this is the case its tested if the two sets commute. If not, 
!  the first group has to be a normal divisor of the product group, i.e. the 
!  product is semi-direct, otherwise it is a direct product. 
! LITERATURE
!  -
! TODO
!  not fully tested
!  documentation  
! RELEASE
!
! PROBLEMS
! If arbitrary groups are used it is not guaranteed, that at the end the product 
! can be written in a symbolic form. Example : C_3z. C_4z = ??
!--------------------------------------------------------------------------------
! SOURCE
*)

GTProductGroup[gr1_, gr2_, OptionsPattern[]] :=
    Module[{whinp,is,l1,l2,gr1m,gr2m,grpm,(*test,*)i,j,lcomm}, 
    whinp = GTWhichInput[gr1[[1]]];
    If[OptionValue[GOFast] == False,
       (*--- check if input is correct ---*)
       If[GTGroupQ[gr1], 
       	  None, 
       	  Print[gr1, " not a group!"]; Abort[]
       ];
       If[GTGroupQ[gr2], 
       	  None, 
       	  Print[gr2, " not a group!"]; Abort[]
       ], 
        None
    ];
    is = Intersection[gr1, gr2];
    If[is == {Ee}, 
       None, 
       Print["Error: more elements than Ee in common : ", is]; Abort[]
    ];
    l1 = Length[gr1]; l2 = Length[gr2];
    gr1m = GTGetMatrix[gr1]; gr2m = GTGetMatrix[gr2]; 
    grpm = GTSimplify[GTSetMultiplication[gr1m, gr2m]];
    (*--- Do the elements commute? ---*)
    (*test = {};*)
    lcomm = True;
    Do[
       Do[
       	  If[gr1m[[i]].gr2m[[j]] == gr2m[[j]].gr1m[[i]],None,lcomm=False] 
             (*test = Append[test, True], 
             test = Append[test, False]*)
          (*]*)
       , {i, 1, l1}]
    , {j, 1, l2}];
    (*If[Intersection[test, {False}] == {}, 
       lcomm = True, 
       lcomm = False
    ];*)
    If[lcomm,
       If[GTInvSubGroupQ[grpm, gr1m], 
          Print["semi-direct product"],Print["direct product group"]], 
       Print["Neither direct product nor semi-direct product"]    
    ];
    Return[GTWhichOutput[grpm, whinp]]
]
(*
***)


(****c* /GTProductGroupQ
! NAME
!  GTProductGroupQ
! AUTHOR
!  W. Hergert
! PACKAGE
!  Basic.m
! MODIFICATION HISTORY
!   July 2015  : first version
!   April 2016 : implemented in basic.m from Test.m
! USAGE
!  Constructs the product group from two groups
! INPUT
!  group - the product group  
!  grp1  - a subgroup
!  grp2  - a subgroup
! OPTIONS
!  GOVerbose - True : output of information
!              False: no output 
! OUTPUT
!  logical value True or False if it is a dirct or semi-direct product group
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGroupQ, GTSubGroupQ, GTGetMatrix, GTSetMultiplication, GTInvSubGroupQ
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
! The input is tested, if we have groups. The it is checked that gr1 anfd gr2 are disjunct.
! Furthermore the subgroup and invariant subgroup relationships are tested.
! The product is constructed and compared with group.
! LITERATURE
!  -
! TODO
!  not fully tested
!  documentation  
! RELEASE
!
! PROBLEMS
! If arbitrary groups are used it is not guaranteed, that at the end the product 
! can be written in a symbolic form. Example : C_3z. C_4z = ??
!--------------------------------------------------------------------------------
! SOURCE
*)

GTProductGroupQ[group_, gr1_, gr2_,OptionsPattern[]] := 
   Module[{verb,ltest,is,lsub1,lsub2,l1,l2,grt,j,lt,lfinal,gr1m,gr2m,grpm,lcomm,test},
   verb = OptionValue[GOVerbose];lfinal = True;ltest=True;
   (*--- check whether input is correct ---*)
   If[OptionValue[GOFast]==False,
      If[GTGroupQ[group],
         None,
         Print[group, " not a group!"]; Abort[]
      ];
      If[GTGroupQ[gr1],
         None,
         Print[gr1, " not a group!"]; Abort[]
      ];
      If[GTGroupQ[gr2],
         None,
         Print[gr2, " not a group!"]; Abort[]
       ],
       None
   ];    
   lsub1 = GTSubGroupQ[group, gr1];
   lsub2 = GTSubGroupQ[group, gr2];
   If[lsub1&&lsub2,
   	  (*--- both are subgroups, continue ---*)
      is = Intersection[gr1, gr2];
      If[is == {Ee},
         None,
         lfinal = False;ltest=False;
         If[verb, 
     	    Print["Error: more elements than Ee in common : ", is], 
     	    None
         ]
      ];
      If[ltest, 
         (*--- gives the product of gr1 and gr2 really group? ---*)
         grt = GTSetMultiplication[gr1,gr2];
         If[Complement[group, grt] == {},
            lt = True,
            lt = False;lfinal=False;
            If[verb,
       	       Print["Not full product, elements left : ",Complement[group,grt]],
               None
            ]
          ];
          (*--- check for type of type of product ---*)
          If[lt,
          	l1 = Length[gr1]; l2 = Length[gr2];
            gr1m = GTGetMatrix[gr1]; gr2m = GTGetMatrix[gr2]; 
            grpm = GTSetMultiplication[gr1m, gr2m];
            (*--- Do the elements commute? ---*)
            test = {};
            Do[
               Do[
       	          If[gr1m[[i]].gr2m[[j]] == gr2m[[j]].gr1m[[i]], 
                     test = Append[test, True], 
                     test = Append[test, False]
                  ]
                , {i, 1, l1}]
            , {j, 1, l2}];
            If[Intersection[test, {False}] == {}, 
               lcomm = True, 
               lcomm = False
            ];
            If[lcomm, 
               If[verb,
               	  Print["product group"],
               	  None
               ];lfinal=True,	   
               If[GTInvSubGroupQ[grpm, gr1m], 
                  If[verb,	
                     Print["semi-direct product"],
                     None
                  ];lfinal=True,   
                  If[verb,
               	     Print["Neither direct product nor semi-direct product"],
               	     None
                  ];lfinal=False
               ]   
            ]  
          ],
          lfinal=False
        ],
        lfinal=False
     ];
     Return[lfinal]
   ]
      
      


(*
***)


(****c* /GTQuaternionToSymbolInfo
! NAME
!  GTQuaternionToSymbolInfo
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2013 first commit
! USAGE
!  GTQuaternionToSymbolInfo[quaternion] gives the following information from a quaternion: {axis,angle,multiplicity,determinant}
! INPUT
!  symbol
! OUTPUT
!  {axis,angle,multiplicity,determinant}
! ERROR MESSAGES
!
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTSymbolInfoToQuaternion 
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)


GTQuaternionToSymbolInfo[quatin_] := Module[{iquat, det, angle, axis, multtmp, mult, quat, sign,rem, dbl},
  iquat = {1, {0, 0, 0}};
  If[quatin == iquat, Return[{{1, 0, 0}, 0, 1, 1}]];
  If[quatin == -I*iquat, Return[{{1, 0, 0}, 0, 1, -1}]];
  If[quatin == -iquat, Return[{{1, 0, 0}, 0, 1, I}]];
  If[quatin == I*iquat, Return[{{1, 0, 0}, 0, 1, -I}]];
  det = If[Element[quatin, Reals], 1, -1];
  quat = If[det < 0, quatin/(-I), quatin];
  rem = Re[quat[[1]]];
  Which[
   rem > 0, 
   		angle = -2*ArcCos[rem] // FullSimplify; dbl = 1; 
   		sign = 1,
   rem == 0,
   		angle = -\[Pi];,
   rem < 0, 
   		angle = -2*ArcCos[-rem] // FullSimplify; dbl = I; 
   		sign = -1];
  axis = quat[[2]]/Sin[angle/2];
  If[Length[Position[axeslist, axis]] == 0, axis = -axis; angle = -angle
  	, If[Length[Position[axeslist, axis]] == 0, Abort[]]];
  If[Abs[angle] == \[Pi],
   Which[
   	GTSymbolInfoToQuaternion[{axis, -\[Pi], 1, I}] == quat, 
   		Return[{axis, -\[Pi], 1, I*det}],
    GTSymbolInfoToQuaternion[{axis, -\[Pi], 1, 1}] == quat, 
    	Return[{axis, -\[Pi], 1, 1*det}]]];
  multtmp = 2 \[Pi]/Abs[angle];
  mult = If[IntegerQ[multtmp], 1, Denominator[multtmp]];
  {axis, angle*sign/mult, mult, det*dbl}]

(*
***)



(****c* /GTQuotientGroup
! NAME
!  GTQuotientGroup
! AUTHOR
!  W. Hergert
! PACKAGE
!  Basic.m
! MODIFICATION HISTORY
!   July 2015  : first version
!   April 2016  : transferred to Basic.m
! USAGE
!  Constructs the multiplication table of a quotient group
! INPUT
!  grp1 - usually the group
!  grp2 - usually the normal divisor
! OPTIONS
!  -
! OUTPUT
!  multiplication table
! GTPack OPTIONS
!  GOVerbose, GOQuotientGroup, GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGroupQ, GTGroupOrder, GTQuotientGroupQ, GTLeftCosets, GTSetMultiplication
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  
! LITERATURE
!  -
! TODO
!  not fully tested
!  documentation  
! RELEASE
!
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTQuotientGroup[grp1_, grp2_, OptionsPattern[]] := 
     Module[{gr1,gr2,verb,lcs,ncs,names,tab,i,j,c1,c2,mult,group,m,printtab},
     If[OptionValue[GOFast] == False,
     (*--- check if input is correct ---*)
        If[GTGroupQ[grp1], 
           None, 
       	   Print[grp1, " not a group!"]; Abort[]
        ];
        If[GTGroupQ[grp2], 
       	   None, 
       	   Print[grp2, " not a group!"]; Abort[]
        ], 
        None
     ];
     (*--- Rename groups with respect to order ---*)	
     If[GTGroupOrder[grp1, GOVerbose -> False] >= GTGroupOrder[grp2, GOVerbose -> False],
	  	 gr1=grp1;gr2=grp2,
	  	 gr1=grp2;gr2=grp1
	  ];
      verb = OptionValue[GOVerbose];
      If[GTQuotientGroupQ[gr1, gr2, GOVerbose -> False],
         None,
         Print[$ErrorQuotientGroup]; Abort[]
      ];
      If[verb,
         Print["Invariant subgoup \[ScriptCapitalN] = ", gr2],
         None
      ];
      lcs = GTLeftCosets[gr1, gr2]; ncs = Length[lcs];
      If[verb,
      	 Print["Left cosets = ", lcs],
      	 None
      ];
      names = Table[0, {ncs}];
      Do[
         If[lcs[[i]] === gr2, 
            names[[i]] = "\[ScriptCapitalN]", 
            names[[i]] = ToString[lcs[[i, 1]]] <> ".\[ScriptCapitalN]"
         ]
      , {i, 1, ncs}];
      tab = Table[0, {ncs}, {ncs}];group={};
      Do[
         Do[
            c1 = lcs[[i]]; c2 = lcs[[j]];
            mult=GTSetMultiplication[c1,c2];
            group=Append[group,mult];
            Do[
               If[Union[mult, lcs[[m]]] == Sort[lcs[[m]]], 
                  tab[[i, j]] = names[[m]],
                  None
               ]
            , {m, 1, ncs}]
        , {i, 1, ncs}]
      , {j, 1, ncs}];
      printtab = Prepend[Table[Prepend[tab[[i]],names[[i]]],{i,1,Length[tab]}],Flatten[{"",names}]];
      Print[Grid[printtab, Frame->All, Alignment -> Center,Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}},Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} -> GTCornerColor}, ItemSize -> Full]];
     (*TableHeadings -> {names, names},*)
      If[OptionValue[GOQuotientGroup],
      	Return[tab],
      	Return[Union[group]]
      ]	
  ]

(*
***)



(****c* /GTQuotientGroupQ
! NAME
!  GTQuotientGroupQ
! AUTHOR
!  W. Hergert
! PACKAGE
!  Basic.m
! MODIFICATION HISTORY
!   July 2015  : first version
!   April 2016 : implemented in Basic.m from Test.m
! USAGE
!  GTQuotientGroupQ checks if a quotient group can be formed from two groups
!  according to the definition
! INPUT
!  grp1 - usually the group
!  grp2 - usually a subgroup of grp1
! OPTIONS
!  GOVerbose -  controls the output (standard ->True)
! OUTPUT
!  True,  False
! GTPack OPTIONS
!  GOVerbose, GOQuotientGroup, GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGroupQ, GTGroupOrder, GTSubGroupQ, GTInvSubGroupQ 
! GTPack NOTEBOOKS 
!  GTSimPack.nb
! DESCRIPTION
!  grp1 and grp2 will be renamed such, that the order of grp1 is larger than that of
!  grp2. Than it will be checked if grp2 is a subgroup and normal divisor of group1.
! LITERATURE
!  -
! TODO
!  not fully tested
!  documentation  
! RELEASE
!
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTQuotientGroupQ[grp1_, grp2_, OptionsPattern[]] := Module[
	  {verb, test,gr1,gr2},
	  If[OptionValue[GOFast] == False,
      (*--- check if input is correct ---*)
         If[GTGroupQ[grp1], 
            None, 
       	    Print[grp1, " not a group!"]; Abort[]
         ];
         If[GTGroupQ[grp2], 
       	    None, 
       	    Print[grp2, " not a group!"]; Abort[]
         ], 
         None
      ];
	  If[GTGroupOrder[grp1, GOVerbose -> False] >= GTGroupOrder[grp2, GOVerbose -> False],
	  	 gr1=grp1;gr2=grp2,
	  	 gr1=grp2;gr2=grp1
	  ];
      verb = OptionValue[GOVerbose]; test = True;
      If[GTSubGroupQ[gr1, gr2],
         None,
         test = False;
         If[verb, 
         	Print["Error: ",grp2," is not a subgroup of ",grp1], 
         	None
         ]
       ];
       If[GTInvSubGroupQ[gr1, gr2],
          test,
          test = False;
         If[verb, 
            Print["Error: ",grp2," is not an invariant subgroup of ",grp1], 
            None
         ]         
       ];
       Return[test]
 ]  
(*
***)


(****c* /GTRightCosets
! NAME
!  GTRightCosets
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  01/01/2013 : Wolfram : order of group and subgroup arbitrary, form of input symbol, matrix, quaternion
! USAGE
!  RigthCoset[group,subgroup] gives the right coset of a group according to a subgroup.
! INPUT
!  group, subgroup
! OUTPUT
!  list of element
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGroupQ, GTWhichInput, GTGetMatrix, GTWhichOutput
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!    
! RELEASE
!
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTRightCosets[group_, sgroup_, OptionsPattern[]] := Module[{el, pos, lc,grp1,grp2,whinp},
  If[Length[group]>Length[sgroup],
    grp1 = GTGetMatrix[group];
    grp2 = GTGetMatrix[sgroup];,
    grp2 = GTGetMatrix[group];
    grp1 = GTGetMatrix[sgroup];];
  whinp = GTWhichInput[group[[1]]];

  If[OptionValue[GOFast],None,
		If[GTGroupQ[grp1],None,Print[$ErrorNoGroupFirst];Abort[]];
		If[GTGroupQ[grp2],None,Print[$ErrorNoGroupSecond];Abort[]];
		If[GTSubGroupQ[grp1,grp2],None,Print["Error: input should contain a group and a subgroup."]];
	];	  
 
  lc = Union[Map[Sort[#] &, 
  	      Table[
  	      	el = GTgmat[grp2[[j]], grp1[[i]]];
  	      	pos = Position[grp1, el];
  	      	grp1[[pos[[1, 1]]]], {i, 1, Length[grp1]}, {j, 1, Length[grp2]}]]];
  GTWhichOutput[Table[SortBy[lc[[i]], GTOrderOfElement[#] &], {i, 1, Length[lc]}], whinp]]


(*
***)



(****g* /GTSelfAdjointQ
! NAME
!  GTSelfAdjointQ
! AUTHOR
!  W. Hergert
! PACKAGE
!  basic.m
! MODIFICATION HISTORY
!   13/11/2014 : moved to basic.m
!   11/30/2012  : first version
! USAGE
!  GTSelfAdjointQ[group,element] gives True if an element of a group is self-adjoint, and gives False otherwise.
! INPUT
!  group   - a group
!  element - a special element of this group
! OUTPUT
!  True/False
! GTPack OPTIONS
!  GTGroupQ, GTWhichInput, GTGetMatrix, GTgmat, GTInverseElement
! STANDARD OPTIONS
!
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  not fully tested
!  documentation  
! RELEASE
!
! PROBLEMS
!
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSelfAdjointQ[group0_, elem0_,OptionsPattern[]] := Module[{k, l, i, elem1, whinp, group, tmp, h1,elem}, 
	If[OptionValue[GOFast],None,
		If[GTGroupQ[group0],None,Print[$ErrorNoGroup];Abort[]];
	];
	
	(*--- check if there is a WhichInput=5 within the list grp ---*)
	whinp = GTWhichInput[elem0];
	elem=If[whinp<5,GTGetMatrix[elem0],elem0];	
	
	tmp = # < 5 & /@ GTWhichInput[group0];
	h1 = And @@ tmp;
	group = If[h1, GTGetMatrix[group0], group0];
	
	l = Length[group];
	k = 0;
	Do[elem1 = GTgmat[GTgmat[group[[i]], elem], GTInverseElement[group[[i]]]];
 	If[elem == elem1, k = k + 1, Null], {i, 1, l}];
	If[k == l, True, False]
]

(*
***) 


(****c* /GTSetMultiplication
! NAME
!  GTSetMultiplication
! AUTHOR
!  W. Hergert
! PACKAGE
!  Basic.m
! MODIFICATION HISTORY
!   July 2015  : first version
!   April 2016  : transferred to Basic.m
! USAGE
!  Given are two sets of symmetry elements. The result is a set of all
!  distinct elements if each element of the one set is multiplied with
!  all elements of the other set.
! INPUT
!  set1 - set of symmetry elements
!  set2 - set of symmetry elements
! OPTIONS
!  -
! OUTPUT
!  set of symmetry elements
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTWhichInput, GTGetMatrix, GTWhichOutput
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  
! LITERATURE
!  -
! TODO
!  not fully tested
!  documentation  
! RELEASE
!
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTSetMultiplication[set1_, set2_] := Module[{nc1,nc2,mult,k,l,s1,s2,whinp},
  nc1 = Length[set1]; nc2 = Length[set2];mult = {};	 
  whinp=GTWhichInput[set1[[1]]];	  
  s1=GTGetMatrix[set1];s2=GTGetMatrix[set2];	
  Do[
     Do[
        mult = Append[mult, s1[[k]].s2[[l]]];
     , {k, 1, nc1}]
  , {l, 1, nc2}]; 
  mult = Union[mult]; 
  Return[GTWhichOutput[mult,whinp]]
]

(*
***)


(****c* /GTSU2ToSymbolInfo
! NAME
!  GTSU2ToSymbolInfo
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2013 first commit
! USAGE
!  GTSU2ToSymbolInfo[symbol] gives the following information: {axis,angle,multiplicity,determinant}
! INPUT
!  symbol
! OUTPUT
!  {axis,angle,multiplicity,determinant}
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTSymbolInfoToSU2Matrix 
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSU2ToSymbolInfo[matin_] := Module[{imat, det, angle, axis, multtmp, mult, mat, sign, rem, dbl},
  imat = IdentityMatrix[2];
  If[matin == imat, Return[{{1, 0, 0}, 0, 1, 1}]];
  If[matin == -I*imat, Return[{{1, 0, 0}, 0, 1, -1}]];
  If[matin == -imat, Return[{{1, 0, 0}, 0, 1, I}]];
  If[matin == I*imat, Return[{{1, 0, 0}, 0, 1, -I}]];
  det = Det[matin];
  mat = If[det < 0, matin/(-I), matin];
  rem = Re[mat[[1, 1]]];
  Which[
   rem > 0, angle = -2*ArcCos[rem] // FullSimplify; dbl = 1; 
   sign = 1,
   rem == 0, angle = -\[Pi];,
   rem < 0, angle = -2*ArcCos[-rem] // FullSimplify; dbl = I; 
   sign = -1];
  axis = {-Im[mat[[1, 2]]], -Re[mat[[1, 2]]], -Im[mat[[1, 1]]]}/
    Sin[angle/2];
  If[Length[Position[axeslist, axis]] == 0, axis = -axis; 
   angle = -angle, If[Length[Position[axeslist, axis]] == 0, Abort[]]];
  If[Abs[angle] == \[Pi],
   Which[GTSymbolInfoToSU2Matrix[{axis, -\[Pi], 1, I}] == mat, 
    Return[{axis, -\[Pi], 1, I*det}],
    	GTSymbolInfoToSU2Matrix[{axis, -\[Pi], 1, 1}] == mat, 
    Return[{axis, -\[Pi], 1, 1*det}]]];
  multtmp = 2 \[Pi]/Abs[angle];
  mult = If[IntegerQ[multtmp], 1, Denominator[multtmp]];
  {axis, angle*sign/mult, mult, det*dbl}]

(*
***)



(****c* /GTSubGroupQ
! NAME
!  GTSubGroupQ
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  12/30/2012,Wolfram : arbitrary order in argument list allowed, works with matrices, symbols and quaternions
!  07/06/2013 Matthias -> Module is in a better shape now
!  11/06/2020 Matthias. Module simplified
! USAGE
!  GTSubGroupQ[group0_,subgroup0_] checks if subgroup0 is a subgroup of group0
! INPUT
!  2 Groups
! OUTPUT
!  logical
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetMatrix, GTGroupQ
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  We have to work with the matrix representation. Using Symbols, Errors occur for angular momentum representations, using GTAngularMomentumRep.
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSubGroupQ[group1_,group2_,OptionsPattern[]]:=Module[{test0,gm1,gm2}, (*{l1,l2,group,subgroup,intersec},*)
	If[OptionValue[GOFast],
	   If[GTGroupQ[group1],None,Print[$ErrorNoGroup," Input 1"];Abort[]];
	   If[GTGroupQ[group1],None,Print[$ErrorNoGroup," Input 1"];Abort[]];
	];
	
	test0=Length[Intersection[group1,group2]] == Min[Length[group1],Length[group2]];
	If[test0,Return[test0],
		gm1=GTGetMatrix[group1];
		gm2=GTGetMatrix[group2];
		Return[Length[Intersection[gm1,gm2]] == Min[Length[gm1],Length[gm2]]]]
]

(*
***)

(****c* /GTSymbolInfoToEulerAngles
! NAME
!  GTSymbolInfoToEulerAngles
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2013 first commit
!  02/17/2015 add gtactvpsv
! USAGE
!  GTSymbolInfoToEulerAngles[symbol] gives Euler angles from the information: {axis,angle,multiplicity,determinant}.
! INPUT
!  symbol info
! OUTPUT
!  Euler Angles
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTnumberlist, axeslist 
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSymbolInfoToEulerAngles[syminf_] := Module[{angle, n, cb, sb, a, b, c,sasq,sa},
  If[syminf[[2]] == 0, Return[{{0, 0, 0}, syminf[[4]]}]];
 (* angle = If[Element[syminf[[4]], Reals], 2*Pi*syminf[[3]]/syminf[[2]],2*Pi*syminf[[3]]/syminf[[2]] + 2*Pi];*)
  angle = -gtactvpsv*2*Pi*syminf[[3]]/syminf[[2]];
  n = syminf[[1]];
  cb = Simplify[1 - 2 (n[[1]]^2 + n[[2]]^2)*Sin[angle/2]^2];
  b = Simplify[ArcCos[cb]];
  sb = Simplify[Sin[b]];
  sasq = Simplify[Sin[angle/2]^2];
  sa=Simplify[Sin[angle]];
  If[b == 0 && Abs[n[[3]]] == 1, a = 0; c = n[[3]]*angle,
   If[b == \[Pi] && Abs[angle] == \[Pi], a = 0; 
    c = 2*ArcTan[n[[2]], n[[1]]],
    a = ArcTan[
      n[[2]]*sa + 
       2*n[[3]]*n[[1]]*sasq, -n[[1]]*sa + 
       2*n[[3]]*n[[2]]*sasq];
    c = ArcTan[n[[2]]*sa - 2*n[[3]]*n[[1]]*sasq, 
      n[[1]]*sa + 2*n[[3]]*n[[2]]*sasq];
    ]];
  {Simplify[{a, b, c}], syminf[[4]]}]

(*
***)

(****c* /GTSymbolToSymbolInfo
! NAME
!  GTSymbolToSymbolInfo
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2014 first commit
!  06/24/2014 major changes
! USAGE
!  GTSymbolToSymbolInfo[symbol] gives the following information: {axis,angle,multiplicity,determinant}
! INPUT
!  symbol
! OUTPUT
!  {axis,angle,multiplicity,determinant}
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTnumberlist, axeslist 
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSymbolToSymbolInfo[insym_] := 
 Module[{charstmp, charstmp2, chars, axis, angle, det, detf, anglefac,
    i}, charstmp2 = Characters[ToString[insym]];
  If[charstmp2[[1]] == "D", 
   charstmp = Table[charstmp2[[i]], {i, 2, Length[charstmp2]}]; 
   detf = I, detf = 1; charstmp = charstmp2];
  Which[charstmp[[1]] == "I", 
   chars = Table[charstmp[[i]], {i, 2, Length[charstmp]}];
   det = -1;
   If[chars[[1]] == "E", Return[{{1, 0, 0}, 0, 1, -1*detf}]], 
   charstmp[[1]] == "C", chars = charstmp;
   det = 1, charstmp[[1]] == "E", Return[{{1, 0, 0}, 0, 1, 1*detf}]];
  axis = axeslist[[Flatten[Position[axeslist, chars[[3]]]][[1]], 2]];
  anglefac = 1;
  Which[Length[chars] == 3, angle = ToExpression[chars[[2]]], 
   Length[chars] == 4, 
   If[Length[Position[GTnumberlist, chars[[4]]]] == 1, 
    anglefac = ToExpression@chars[[4]];
    angle = ToExpression[chars[[2]]], 
    angle = -ToExpression[chars[[2]]]], Length[chars] == 5, 
   anglefac = ToExpression@chars[[4]];
   angle = -ToExpression[chars[[2]]]];
  Return[{axis, angle*gtactvpsv, anglefac, det*detf}]]

(*
***)

(****c* /GTSymbolInfoToO3Matrix
! NAME
!  GTSymbolInfoToO3Matrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2014 first commit
!  06/24/2014 major changes
!  02/17/2015 add gtactvpsv
! USAGE
!  GTSymbolInfoToO3Matrix[symbol info] gives a rotation matrix from the information {axis,angle,multiplicity,determinant}
! INPUT
!  symbol info
! OUTPUT
!  3 dimensional rotation matrix
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!   
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
!
!	to rotation matrices for ordinary elements and duble group elements are identical
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)
GTSymbolInfoToO3Matrix[syminf_] := Module[{det, angle},
  If[Element[syminf[[4]], Reals], det = syminf[[4]], 
   det = syminf[[4]]/I];
  If[Abs[syminf[[2]]] > 0, angle = -gtactvpsv*syminf[[3]]*2*\[Pi]/syminf[[2]], 
   angle = 0];
  Return[Simplify@det*RotationMatrix[angle, syminf[[1]]]]
  ]
 
(*
***)

(****c* /GTSymbolInfoToSU2Matrix
! NAME
!  GTSymbolInfoToSU2Matrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2013 first commit
!  07/04/2014 new module
!  02/17/2015 add gtactvpsv
! USAGE
!  GTSymbolInfoToSU2Matrix[symbol info] gives a spin rotation matrix from the information {axis,angle,multiplicity,determinant}
! INPUT
!  symbol info
! OUTPUT
!  2 dimensional spin rotation matrix
! ERROR MESSAGES
!
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTSU2Matrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
!
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)
  
 GTSymbolInfoToSU2Matrix[syminf_] := Module[{angle, dble},
  dble = -syminf[[4]]^2;
  angle = -gtactvpsv*syminf[[3]]*2*\[Pi]/syminf[[2]];
  Return[dble*GTSU2Matrix[angle, syminf[[1]]]]
  ]

(*
***)

(****c* /GTSymbolInfoToSymbol
! NAME
!  GTSymbolInfoToSymbol
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2014 first commit
!  06/24/2014 major changes
!  07/21/2015 active rotation were introduced
!			  change of the line If[angle == "2", invprint = 1]; to If[angle == "2", invprint = 1;inv=""];
! USAGE
!  GTSymbolInfoToSymbol[symbol info] gives a symbol from the information {axis,angle,multiplicity,determinant}
! INPUT
!  symbol info
! OUTPUT
!  symbol
! ERROR MESSAGES
!
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!   
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
!	GTSymbolInfoToSymbol automatically formats the output to standard form!
!
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSymbolInfoToSymbol[syminf_] := 
 Module[{Axislabel, angle, inv, invprint, mult, multprint, Inver, 
   Inverprint, dble, sym, symprint, symlabel, pos}, 
  If[syminf[[3]] == 0, symlabel = "Ee"; angle = ""; Axislabel = "", 
   angle = ToString[Abs[syminf[[2]]]];
   pos = Position[axeslist, syminf[[1]]][[1, 1]];
   Axislabel = ToString[axeslist[[pos, 1]]];
   symlabel = "C"];
 (* {inv, invprint} = If[gtactvpsv*Sign[syminf[[2]]] > 0, {"", 1}, {"i", -1}];*)
  {inv, invprint} = If[Sign[syminf[[2]]] > 0, {"", 1}, {"i", -1}];
  {mult, multprint} = 
   If[syminf[[3]] > 1, {ToString[syminf[[3]]], syminf[[3]]}, {"", 
     syminf[[3]]}];
  {Inver, Inverprint, dble} = 
   If[Element[syminf[[4]], Reals], 
    If[syminf[[4]] > 0, {"", "", 0}, {"I", "I", 0}], 
    If[syminf[[4]]/I > 0, {"D", "", 1}, {"DI", "I", 1}]];
  (*If[dble == 1, invprint = -invprint];*)
  If[angle == "2", invprint = 1;inv=""];
  sym = ToExpression[
    Inver <> symlabel <> angle <> Axislabel <> mult <> inv];
  symprint = 
   Power[Subscript[
     Overscript[Inverprint <> symlabel, 
      If[dble == 1, 
       Which[StringLength[Inverprint] == 0, "_", 
        StringLength[Inverprint] == 1, "__", 
        StringLength[Inverprint] == 2, "___"], ""]], 
     angle <> Axislabel], 
    ToString[
     If[invprint*multprint == 1 || invprint*multprint == 0, "", 
      invprint*multprint]]];
  Format[sym, StandardForm] := symprint;
  sym]
(*
***)

(****c* /GTSymbolInfoToO2Matrix
! NAME
!  GTSymbolInfoToO2Matrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2013 first commit
! USAGE
!  GTSymbolInfoToO2Matrix[symbol info] gives a rotation matrix from the information {axis,angle,multiplicity,determinant}
! INPUT
!  symbol info
! OUTPUT
!  2 dimensional rotation matrix
! ERROR MESSAGES
!
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTSymbolInfoToO3Matrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
!
!	to rotation matrices for ordinary elements and duble group elements are identical
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)
  
 GTSymbolInfoToO2Matrix[syminf_] := {{1, 0, 0}, {0, 1, 0}}.GTSymbolInfoToO3Matrix[syminf].{{1, 0}, {0, 1}, {0, 0}}

(*
***)

(****c* /GTSymbolInfoToQuaternion
! NAME
!  GTSymbolInfoToQuaternion
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  02/21/2013 first commit
!  02/17/2015 add gtactvpsv
! USAGE
!  GTSymbolInfoToQuaternion[symbol info] gives a quaternion from the information {axis,angle,multiplicity,determinant}
! INPUT
!  symbol info
! OUTPUT
!  quaternion
! ERROR MESSAGES
!
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!   
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  multiplicity it needed for 5, 7 and 9 fold rotations. E.g. for C5z2 the multiplicity is 2.
!  determinant can have the following values:
!	 1: proper roations
!	-1: improper rotations
!	 I: proper rotations, double group element
!	-I: improper rotations, double group element
!
!	to rotation matrices for ordinary elements and duble group elements are identical
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSymbolInfoToQuaternion[syminf_] := Module[{angsum, angle, a, b, inf},
  If[Element[syminf[[4]], Reals],angsum = 0;
   		If[syminf[[4]] > 0, inf = 1, inf = -I];
   ,angsum = 2 \[Pi];
   	If[syminf[[4]]/I > 0, inf = 1, inf = -I];];
  angle = angsum + gtactvpsv*syminf[[2]]*syminf[[3]];
  a = Cos[angle/2];
  b = Sin[angle/2];
  inf*{a, b*syminf[[1]]}]
  
(*
***)

(****c* /GTTransformation
! NAME
!  GTTransformation
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  12/30/2012 : improvemet of description
! USAGE
!  GTTransformation[transformation,vector] applies a coordinate transformation to a given vector.
! INPUT
!  transformation - in form of <R(T),t(t)>
!  vector         - space vector
! OUTPUT
!  r'  - the transformed vector
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetMatrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTransformation[sym_, r_] :=Module[{he},
	he=Head[sym];
	If[he===AngleBracket , GTGetMatrix[sym[[1]]].r+sym[[2]] ,
	Print["Error: Input does not describe a transformation!"] ]   ]
	
(*
***)

(****c* /GTTransformationOperator
! NAME
!  GTTransformationOperator
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  07/22/2013 : initial implementation
!  05/09/2021 : fixed bug with argument (f(R^-1 {x,y,z}) instead of f(R{x,y,z}) )
! USAGE
!  GTTransformationOperator[transformation,function,arguments] applies a coordinate transformation to a certain function with given arguments.
! INPUT
!  transformation - group element
!  function
!  arguments
! OUTPUT
!  f' - transformed function
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTWhichRepresentation, GTGetRotationMatrix, GTGetMatrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTransformationOperator[sym_, fu_, xin_] := Module[{rotmat, mat, arg, dimlist, fdimlist, futmp, fudim,testfu},
  dimlist = {3, 3, 2, None, 3};
  fdimlist = {3, 2, 2, None, 2};
  
  If[grpdgrp == 4, 
  	Print[$ErrorTransformationNoPermutation]; 
    Abort[]];
    
  If[Length[xin] == dimlist[[grpdgrp]], 
  	None, GTWhichRepresentation[];
    Print[$ErrorTransformationArgumentDimension]; 
    Abort[]];
    
  rotmat = GTGetRotationMatrix[sym];
  mat = GTGetMatrix[sym];
  arg = Inverse[rotmat].xin;
  futmp = Apply[fu, arg];
  testfu=Apply[fu, Table[1, {i, 1, dimlist[[grpdgrp]]}]];
  fudim = If[VectorQ[testfu], Length[testfu],0];
  Which[
  	(*--- O(3) representation ---*)
  	grpdgrp == 1,
   		If[fudim == 0, Return[futmp],
     	 If[fudim == 3, Return[mat.futmp],
       		Print[$ErrorTransformationO3];
            Abort[];
         ];
        ],
    (*--- SU(2) representation ---*)
  	grpdgrp == 2,
   		If[fudim == 2, Return[mat.futmp],
       		Print[$ErrorTransformationSU2];
            Abort[];
        ],
    (*--- O(2) representation ---*)
  	grpdgrp == 3,
   		If[fudim == 0, Return[futmp],
     	 If[fudim == 2, Return[mat.futmp],
       		Print[$ErrorTransformationO2];
            Abort[];
         ];
        ],
    (*--- SU(2)xS representation ---*)
  	grpdgrp == 5,
     	If[fudim == 2, Return[mat[[1;;2,1;;2]].futmp],
       		Print[$ErrorTransformationSU2];
            Abort[];
        ];
   ]
  ]
  
(*
***)


(****c* /GTType
! NAME
!  GTType
! AUTHOR
!  W. Hergert, M.Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  12/29/2012 : installation in package
!  07/06/2013 : Rewritten using GTWhichInput
! USAGE
!  GTType[point group element] transforms a point group elment in a space group element
! INPUT
!  point group element
! OUTPUT
!  A spacegroup element in the standard notation. The rotational part is given in symbolic form.
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTWhichInput
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  GTType is used in the definition of the general group operation GTgmat. If one of the elements in the multiplication
!  is a spacegroup element, both elements are transformed to spacegroup elements. This transformation is done by GTType.
!  input                   result
!  spacegroup element      not changed, spacegroup element
!  pointgroup (symbol)     <symbol, {0,0,0}> 
!  pointgroup (matrix)     matrix -> symbol -> <symbol, {0,0,0}>
! LITERATURE
!  
! TODO
!  
! RELEASE
!  
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTType[symb_] :=Module[{hinp}, 
	hinp=GTWhichInput[symb];
	If[hinp==0,Print[$ErrorInput];Return[],None];
	If[hinp==5,Return[symb],
		Return[\[LeftAngleBracket] symb,{0,0,0}\[RightAngleBracket]]]]		

(*
***)
  
(****c* /GTWhichInput
! NAME
!  GTWhichInput
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  16/05/2013 Including Angle Bracket
! USAGE
!  GTWhichInput[symmetry element] determines the type of a symmetry element (Symbol, Matrix, Quaternion, Euler Angles or space group elements).
! INPUT
!  symmetry element
! OUTPUT
!  index - 0 for unidentified
!          1 for symbol
!          2 for matrix
!          3 for quaternion
!          4 for euler angles
!          5 for space group elements
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTSymbolQ, GTQuaternionQ, GTEulerAnglesQ
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  An error message should appear if in_ is a list containing elements of differnt type.
!  Extension to space group elements case 5. This will allow to react with error message 
!  in other moduls     
! RELEASE
!
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTWhichInputSingle[in_] := Module[{},
  If[GTSymbolQ[in], Return[1]];
  If[Last@Characters@ToString[in] == "'", Return[1]];
  If[MatrixQ[in], Return[2]];
  If[GTQuaternionQ[in], Return[3]];
  If[GTEulerAnglesQ[in], Return[4]];
  If[Head[in] === AngleBracket, Return[5]];
  Return[0]]

GTWhichInput[in_] := Module[{whinp},
	whinp = GTWhichInputSingle[in];
	If[whinp == 0,If[Length[in] > 0, Map[GTWhichInputSingle, in], 0],whinp]]
(*
***)

(****c* /GTWhichOutput
! NAME
!  GTWhichOutput
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  07/06/2013 initial module
! USAGE
!  GTWhichOutput[symmetry element,type index] transforms an symmetry element to either a symbol, a matrix, a quaternion, a set of Euler Angles or a space group element, according to the given type index.
! INPUT
!  symmetry element
!  type index - 1 for symbol
!               2 for matrix
!               3 for quaternion
!               4 for euler angles
!               5 for space group elements
! OUTPUT
!  symmetry element
! GTPack OPTIONS
! 
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetSymbol, GTGetMatrix, GTGetQuaternion, GTGetEulerAngles, GTType
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!
! TODO
!  
! RELEASE
! 
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTWhichOutput[element_,index_] := Which[	index==1,	GTGetSymbol[element],
											index==2,	GTGetMatrix[element],
											index==3,	GTGetQuaternion[element],
											index==4,	GTGetEulerAngles[element],
											index==5,	GTType[element]]
(*
***)

Attributes[GTAbelianQ]={Protected, ReadProtected}
Attributes[GTAllSymbols]={Protected, ReadProtected}
Attributes[GTCenter]={Protected, ReadProtected}
Attributes[GTClasses]={Protected, ReadProtected}
Attributes[GTClassMult]={Protected, ReadProtected}
Attributes[GTClassMultTable]={Protected, ReadProtected}
Attributes[GTCyclicQ]={Protected, ReadProtected}
Attributes[GTConjugateElement]={Protected, ReadProtected}
Attributes[GTGenerators]={Protected, ReadProtected}
Attributes[GTGetEulerAngles]={Protected, ReadProtected}
Attributes[GTGetMatrix]={Protected, ReadProtected}
Attributes[GTGetQuaternion]={Protected, ReadProtected}
Attributes[GTGetSubGroups]={Protected, ReadProtected}
Attributes[GTGetSymbol]={Protected, ReadProtected}
Attributes[GTSU2Matrix]={Protected, ReadProtected}
Attributes[GTGetRotationMatrix]={Protected, ReadProtected}
Attributes[GTgmat]={Protected, ReadProtected}
Attributes[GTGroupOrder]={Protected, ReadProtected}
Attributes[GTGroupQ]={Protected, ReadProtected}
Attributes[GTInverseElement]={Protected, ReadProtected}
Attributes[GTInvSubGroups]={Protected, ReadProtected}
Attributes[GTInvSubGroupQ]={Protected, ReadProtected}
Attributes[GTLeftCosets]={Protected, ReadProtected}
Attributes[GTMultTable]={Protected, ReadProtected}
Attributes[GTOrderOfElement]={Protected, ReadProtected}
Attributes[GTSubGroupQ]={Protected, ReadProtected}
Attributes[GTSelfAdjointQ]={Protected, ReadProtected}
Attributes[GTRightCosets]={Protected, ReadProtected}
Attributes[GTTransformation]={Protected, ReadProtected}
Attributes[GTTransformationOperator]={Protected, ReadProtected}
Attributes[GTType]={Protected, ReadProtected}
Attributes[GTWhichInput]={Protected, ReadProtected}
Attributes[GTWhichOutput]={Protected, ReadProtected}

End[] (* End Private Context *)


(*/*-------------------------------------------- Commands for GTPack2.0 ---------------------------------------------*)


(*---------------- Usage ---------------*) 
 

(*---------------- Options -------------*)

(*---------------- Error messages ------*)



Begin["`Private`"] (* Begin Private Context  GTPack2.0*) 

(****q* 
! NAME
! 
! AUTHOR
! 
! PACKAGE
!  photonics.m 
! MODIFICATION HISTORY
!        
! USAGE
!  
! INPUT
! 
! OPTIONS 
!   
! OUTPUT
!  
! ERROR MESSAGES
! 
! GTPack MODULES
!  - 
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
! SOURCE
!--------------------------------------------------------------------------------
!
*)


(*
***)
 
End[] (* End Private Context GTPack2.0*)

(**/---------------------- End GTPack2.0 --------------------------------------*)



EndPackage[]
