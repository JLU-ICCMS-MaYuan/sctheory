(****m* /RepresentationTheory.m
!
! NAME
!  RepresentationTheory.m
! AUTHOR
!  M. Geilhufe, W. Hergert, S. Schenk
! MODIFICATION HISTORY
!  2/24/12 : initial documentation 
! USAGE
!  Contains modules belonging to representation theory.
! 
! ERROR MESSAGES
!  
! GTPack MODULES
!
! --- Characters ---
! GTCharacterTable           - gives the character table of a group
! GTCharacters               - Characters of a matrix representation
! GTExtraRepresentations     - extracts the extra representations from the character table of a double group
! GTReorderCharacterTable    - reorders a charcter table
! GTReality                  - estimates if a representation is potentially real, essentially complex or pseudo-real from the theorem of Frobenius and Schur
! GTSpinCharacters           - gives the character of the spinor representation for each class
! GTSOCSplitting             - calculates the splitting of states due to spin-orbit coupling
! 
! --- Projection operators ---
! GTCharProjectionOperator   - applies the character projection operator of a irreducible representation to a certain function with given arguments
! GTProjectionOperator       - gives the part of a given function with arguments which transforms like the m-th row and the n-th column of an irreducible representation 
! GTWignerProjectionOperator - applies the projection operator on spherical harmonics
!
! --- Auxiliary representations ---
! GTAngularMomentumChars     - calculates the character system using (2j+1)-dimensional representation matrices according to the angular momentum operator
! GTAngularMomentumRep       - gives a (2j+1)-dimensional matrix representation using the WignerD-function
! GTRegularRepresentation    - gives the regular representation
! GTVectorRep                - gives the vector representation
! 
! --- Irreducible representations ---  
! GTClebschGordanSum         - gives the direct sum of two representations 
! GTGetIrep                  - gives the representation matrices of an irreducible representation within the caracter table
! GTIrep                     - gives the number of times an irreduzible representation appears in a reducible representation with given characters
! GTIrepDimension            - gives the dimension of an irreducible representations
! GTIrepMatrixView           - gives a neat print of the matrices of the irreducible representation n of a group 
! GTNumberOfIreps            - gives the number of non-equivalent irreducible representations
!
! --- direct product representations ---  
! GTClebschGordanCoefficients - gives an array containing Clebsch-Gordan-Coefficients for the calculation of basis functions of Irep3 within a direct product representation of Irep1 and Irep2
! GTClebschGordanTable        - illustrates a list of Clebsch-Gordan coefficients calculated for a direct product representation
! GTDirectProductChars        - gives the character system of a direct product representation from the character systems of two given representations
! GTDirectProductRep          - gives the direct product representation of two representations 
 !GTSymmetrizedProductChars   - gives charactes of symmetrized products
! GTPack NOTEBOOKS 
!  none in the moment
! DESCRIPTION
!  RepresentationTheory.m contains all modules which are connected with representation theory of groups. 
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

BeginPackage["GroupTheory`RepresentationTheory`",{"GroupTheory`Basic`","GroupTheory`Symbols`","GroupTheory`Install`","GroupTheory`AngularMomentum`","GroupTheory`Auxiliary`"}]

(*-------------------------- characters --------------------------*)
 GTCharacterTable				::usage = "GTCharacterTable[\*StyleBox[\"group\",\"TI\"]] gives the character table of a \*StyleBox[\"group\", \"TI\"]."
 GTCharacters                   ::usage = "GTCharacters[\*StyleBox[\"matrices\",\"TI\"]] gives the characters of a representation matrix or a list of representation matrices."
 GTExtraRepresentations         ::usage = "GTExtraRepresentations[\*StyleBox[\"character table\",\"TI\"]] extracts the extra representations from the character table of a double group."
 GTReorderCharacterTable        ::usage = "GTReorderCharacterTable[group symbol, character table, order vector] reorders a charcter table."
 GTReality                      ::usage = "GTReality[\*StyleBox[\"classes,characters\",\"TI\"]] estimates if a representation is potentially real, essentially complex or pseudo-real from the theorem of Frobenius and Schur."
 GTSpinCharacters				::usage = "GTSpinCharacters[\*StyleBox[\"classes\",\"TI\"]] gives the character of the spinor representation for each \*StyleBox[\"class\"]."
 GTSOCSplitting                 ::usage = "GTSOCSplitting[\*StyleBox[\"character table\",\"TI\"]] calculates the splitting of states due to spin-orbit coupling."
 
(*-------------------------- Projection operators --------------------------*)
 GTCharProjectionOperator		::usage = "GTCharProjectionOperator[\*StyleBox[\"classes,characters,function,arguments\",\"TI\"]] applies the character projection operator corresponding to a certain irreducible representation to a given \*StyleBox[\"function\", \"TI\"] with given \*StyleBox[\"arguments\", \"TI\"]."
 GTProjectionOperator			::usage = "GTProjectionOperator[\*StyleBox[\"group,ireducible representation,m,n,function,arguments\",\"TI\"]] gives the part of a given \*StyleBox[\"function\",\"TI\"] with \*StyleBox[\"arguments\",\"TI\"] which transforms like the \*StyleBox[\"m\",\"TI\"]-th row and the \*StyleBox[\"n\",\"TI\"]-th column of an \*StyleBox[\"irreducible representation\",\"TI\"]."
 GTWignerProjectionOperator		::usage = "GTWignerProjectionOperator[\*StyleBox[\"group,irep,l,m,i,j\",\"TI\"]] applies the projection operator on spherical harmonics."

(*-------------------------- Auxiliary representations --------------------------*)
 GTRegularRepresentation        ::usage = "GTRegularRepresentation[\*StyleBox[\"group\",\"TI\"]] gives the regular representation of a \*StyleBox[\"group\",\"TI\"]."
 GTVectorRep                    ::usage=  "GTVectorRep[\*StyleBox[\"character table\",\"TI\"]]  gives the vector representation for a group with \*StyleBox[\"character table\",\"TI\"]."
(*-------------------------- Irreducible representations -------------------------------*)  
 GTClebschGordanSum				::usage = "GTClebschGordanSum[\*StyleBox[\"rep1,rep2\",\"TI\"]] gives the direct sum of two representations \*StyleBox[\"rep1\",\"TI\"] and \*StyleBox[\"rep2\",\"TI\"]."
 GTGetIreps                     ::usage = "GTGetIreps[\*StyleBox[\"group\",\"TI\"]] gives the character table and the representation matrices of a point \*StyleBox[\"group\",\"TI\"]."
 GTGetIrep    					::usage = "GTGetIrep[\*StyleBox[\"group, index, (character table)\",\"TI\"]] gives the representation matrices of an irreducible representation (denoted by its \*StyleBox[\"index\",\"TI\"] within the \*StyleBox[\"character table\",\"TI\"])."
 GTIrep        					::usage = "GTIrep[\*StyleBox[\"characters,character table\",\"TI\"]] gives the number of times \*SubscriptBox[StyleBox[\"n\", \"TI\"], StyleBox[\"p\", \"TI\"]] an irreduzible representation \*SuperscriptBox[\[CapitalGamma], \"p\"] appears in a reducible representation with given \*StyleBox[\"characters\", \"TI\"]."
 GTIrepDimension				::usage = "GTIrepDimension[\*StyleBox[\"irep\",\"TI\"]] gives the dimension of an irreducible representations."
 GTIrepMatrixView               ::usage = "GTRepMatrixView[\*StyleBox[\"group,character table,number\",\"TI\"]] illustrates the matrices of the irreducible representation \*StyleBox[\"number\",\"TI\"] of a \*StyleBox[\"group\",\"TI\"]." 
 GTNumberOfIreps				::usage = "GTNumberOfIreps[\*StyleBox[\"group\",\"TI\"]] gives the number of non-equivalent irreducible representations."

(*-------------------------- direct product representations ---------------------*)  
 GTClebschGordanCoefficients	::usage = "GTClebschGordanCoefficients[\*StyleBox[\"Irep1,Irep2,Irep3\",\"TI\"]] gives an array containing Clebsch-Gordan-Coefficients for the calculation of basis functions of \*StyleBox[\"Irep3\",\"TI\"] within a direct product representation of \*StyleBox[\"Irep1\",\"TI\"] and \*StyleBox[\"Irep2\",\"TI\"]."
 GTClebschGordanTable        	::usage = "GTClebschGordanTable[\*StyleBox[\"SuperscriptBox[\"\[CapitalGamma]\",\"1\"],SuperscriptBox[\"\[CapitalGamma]\",\"2\"],Clebsch-Gordan coefficients,names\",\"TI\"]] illustrates a list of \*StyleBox[\"Clebsch-Gordan coefficients\",\"TI\"] calculated for the direct product representation \*Cell[BoxData[FormBox[RowBox[{StyleBox[SuperscriptBox[\"\[CapitalGamma]\",\"1\"],\"TI\"],\"\[CircleTimes]\",StyleBox[SuperscriptBox[\"\[CapitalGamma]\",\"2\"],\"TI\"]}],TraditionalForm]],\"InlineMath\"]."
 GTDirectProductChars			::usage = "GTDirectProductChars[\*StyleBox[\"characters1,characters2\",\"TI\"]] gives the character system of a direct product representation from the character systems of two given representations."
 GTDirectProductRep				::usage = "GTDirectProductRep[\*StyleBox[\"rep1,rep2\",\"TI\"]] gives the direct product representation of two representations \*StyleBox[\"rep1\",\"TI\"] and \*StyleBox[\"rep2\",\"TI\"]."
 GTSymmetrizedProductChars      ::usage=  "GTSymmetrizedProductChars[character table, irepn, n , symmetry] gives the characters of the nth power of the irepn-th IREP of a group with character table."
(**************************************************************************************************)

(*--------------------------- Options ----------------------------*)
Options[GTCharacterTable]            = {GOIrepNotation->"Bethe",GOVerbose->True,GOFast->GOFastValue,GOReality->False,GOMethod->"NumericalApproximant"}
Options[GTCharProjectionOperator]    = {GOFast->GOFastValue,GOLattice->{{0,0,0},{0,0,0},{0,0,0}}}
Options[GTClebschGordanCoefficients] = {GOFast->GOFastValue}
Options[GTClebschGordanSum]          = {GOFast->GOFastValue}
Options[GTDirectProductRep]          = {GOFast->GOFastValue,GODiagonal->True}
Options[GTExtraRepresentations]      = {GOFast -> GOFastValue, GOVerbose -> True,GOIrepNotation->{}}
Options[GTGetIreps]                  = {GOIrepNotation->"Bethe",GOVerbose->True,GOFast->GOFastValue};
Options[GTGetIrep]                   = {GOFast->GOFastValue,GOMethod->"Cornwell",GOlmax->15}
Options[GTIrep]                      = {GOVerbose->True,GOFast->GOFastValue}
Options[GTIrepDimension]             = {GOFast->GOFastValue}
Options[GTNumberOfIreps]             = {GOFast->GOFastValue,GOVerbose->True}
Options[GTProjectionOperator]        = {GOFast->GOFastValue}
Options[GTReality]                   = {GOFast->GOFastValue}
Options[GTRegularRepresentation]     = {GOFast->GOFastValue}
Options[GTReorderCharacterTable]     = {GOIrepNotation -> 0}
Options[GTSOCSplitting]              = {GOFast->GOFastValue}
Options[GTSpinCharacters]            = {GOFast->GOFastValue}
Options[GTWignerProjectionOperator]  = {GOFast->GOFastValue, GOHarmonics->"Complex"}
Options[GTSymmetrizedProductChars]   = {GOVerbose->True}
Options[GTVectorRep]                 = {GOVerbose -> False}
Options[GTCharacters]                = {GOClasses -> False}


GTIrepNotation = GroupTheory`Auxiliary`Private`GTIrepNotation

Begin["`Private`"] (* Begin Private Context *) 
(*--------------------------- external Modules -------------------*)
(* only non public functions must have this *)

(*--------------------------- Symolic Operators ----------------------------*)
(*
HoldPattern[cg1_ \[CirclePlus] cg2_] := GTClebschGordanSum[cg1,cg2]
HoldPattern[rep1_ \[CircleTimes] rep2_] := GTDirectProductRep[rep1,rep2]
*)

(*--------------------------- Used Errors ----------------------------*)
(*  ---- names may include .m name 
$ErrorInput
$ErrorInvCt
$ErrorNoGroup
$ErrorProjectionNo2D
$ErrorProjectionNoPermutation
$ErrorProjectionNoSpinor
$ErrorProjectionOperator
*)

(*--------------------------- Modules ----------------------------*)

(****n* /GTAxisOrder
! NAME
!  GTAxisOrder
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m
! MODIFICATION HISTORY
!   October 2014: first version
! USAGE
!   gives the index of the rotation axis of a given element within axeslist.
! INPUT
!   rotation matrix
! OUTPUT
!  index
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetRotationMatrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTClassInfo for the denotation of elements in the character table
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
 GTAxisOrder[in_] := Module[{es, apos, ev, pos},
  es = Expand@Eigensystem[GTGetRotationMatrix[in]];
  apos = First@Flatten@Position[es[[1]], 1];
  ev = Normalize[es[[2, apos]]];
  pos = Position[axeslist, ev];
  If[Length[pos] == 0, pos = Position[axeslist, -ev]];
  First@First@pos]
(*
***)



(****n* /GTCharacterTable
! NAME
!  GTCharacterTable
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  December 2014 alignment of the characters within print out -> right
!  October 2014 Mulliken and Bouckaert notation  
!  June 2014 new algorithm
!  03/05/2013 Robodoc documentation updated
!  15/05/2013 check grpdgrp; Set Bethe notation for SU(2) representation
!		  check if input is a group
!  October 2014, Irep notation
!  23/08/2016 compatibility for space group elements
! USAGE
!  GTCharacterTable[group] gives the character table of a group.
! INPUT
!  group
! OUTPUT
!  list with classes and characters.
! GTPack OPTIONS
!  GOIrepNotation, GOVerbose, GOFast, GOReality
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTClasses, GTClassMult, GTGroupQ, GTIrepNotation, GTSimplify
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!
! LITERATURE
!
! TODO
!  see problems
! Release
!
! PROBLEMS
!  Errors occur for GTCharacterTable if representation 4 is used 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTCharacterTable[grp0_,OptionsPattern[]] := Module[{grp,k,j,i,m,n,inp,run,clm, h, g, nc, h0mat, h1mat, hmat, hh,Vt,V, Vi, cc, hma, dim, ct, vec, cts, notation,label,Liste, cla,Pr,chi,cl,grpmat,el,pos,clmprint,irepnotation,notation0,randint,posidentity},
  (*Print a warning in case magnetic groups are used.*)
  If[GTMagneticQ[], "Warning: You are using magnetic groups. These groups contain anti-unitary elements and need a special treatment concerning the irreducible representations. This feature is not implemented yet and the character table you obtain might not be the one you were aiming for. We are currently working on that!", None];  	
  (*It is possible to specify a list of irep names. A random number will be used to generate an internal string to distinguish this case from "Bethe",...*)
  randint = "list"<>ToString[RandomInteger[50000]];	
  irepnotation = If[Length[OptionValue[GOIrepNotation]]==0,OptionValue[GOIrepNotation],randint];
	
  inp = GTWhichInput[grp0[[1]]];
  
  (*Take rotational part from space group elements*)
  grp= If[inp==5,
  	    grp0[[;;,1]]
 		(*Union[grp0[[;;,1]],grp0[[;;,1]]]*)	
  		,grp0];
  
  (*--- Check the input ---*)
  If[OptionValue[GOFast],run=True,run=GTGroupQ[grp]];
  If[Not[run],Print[$ErrorNoGroup];Abort[]];	

  Which[	grpdgrp==1, (*--- O(3) ---*)
				notation0 = irepnotation,
			grpdgrp==2, (*--- SU(2) ---*)
			    notation0 = If[irepnotation==randint,randint,"Bethe"];
				If[OptionValue[GOVerbose],
					If[irepnotation=="Mulliken",
						Print["Warning: The standard representation is given by SU(2) matrices. Mulliken symbols can not be evaluated for double groups. The notation is switched to Bethe notation."]
						,None];
					If[irepnotation=="Bouckaert",
						Print["Warning: The standard representation is given by SU(2) matrices. Bouckaert notation can not be evaluated for double groups. The notation is switched to Bethe notation."]
						,None]
					,None],
			grpdgrp==3, (*--- O(2) ---*)
				notation0 = irepnotation,
			grpdgrp==4, (* ---Own ---*)
				notation0 = irepnotation,
			grpdgrp==5, (*--- SU(2)xS ---*)
			    notation0 = If[irepnotation==randint,randint,"Bethe"];
				If[OptionValue[GOVerbose],
					If[irepnotation=="Mulliken",
						Print["Warning: The standard representation is given by SU(2)xS matrices. Mulliken symbols can not be evaluated for double groups. The notation is switched to Bethe notation."]
						,None];
					If[irepnotation=="Bouckaert",
						Print["Warning: The standard representation is given by SU(2)xS matrices. Bouckaert notation can not be evaluated for double groups. The notation is switched to Bethe notation."]
						,None]
					,None]];

  notation = If[notation0=="Bethe"||notation0=="Mulliken"||notation0=="Bouckaert"||notation0==randint,notation0,Print["Warning: Notation unknown. The notation is switched to Bethe notation."];"Bethe"];

  (*--- some preliminaries ---*)
  grpmat=GTGetMatrix[grp];
  clm = GTClasses[grpmat, GOFast -> True];
  g = Length[grp];
  h = Length[clm];
  nc = Map[Length, clm];					
  If[g > 1,
	  (*--- calculate the class constant matrices ---*)
	  h0mat = Table[Table[Flatten[Table[GTSimplify[clm[[j, m]].clm[[i, n]]], {m, 1, nc[[j]]}, {n, 1,nc[[i]]}], 1], {i, 1, h}], {j, 1, h}];
	  h1mat = Table[Table[Length[Position[h0mat[[j, i]], GTSimplify[clm[[k, 1]]]]], {j, 1,h}, {i, 1, h}], {k, 1, h}];
	  hmat = Table[h1mat[[k, i, j]], {i, 1, h}, {j, 1, h}, {k, 1, h}];
	  (*--- find the matrix V that diagonalizes the class constant matrices ---*)
	  (*If[h>8,*)
	  Which[
	  	OptionValue[GOMethod]=="NumericalApproximant",
	    hh = Transpose[Sum[SetAccuracy[RandomReal[], 50]*hmat[[i]], {i, 2, h}]];
	    Vt = Chop[Eigenvectors[hh]];
	    V = ToRadicals@RootApproximant[Vt];
	    Vi = ToRadicals@RootApproximant[Inverse[Vt]];,
	    OptionValue[GOMethod]=="Analytic",
	    V = Eigenvectors[Transpose[Sum[(RandomInteger[10]+1)/(RandomInteger[100]+1)*hmat[[i]], {i, 2, h}]]];
	    Vi = Inverse[V];,
	    OptionValue[GOMethod]=="Numeric",
	    V = Eigenvectors[Transpose[Sum[RandomReal[]*hmat[[i]], {i, 2, h}]]];
	    Vi = Inverse[V];
	  ];
	  
	  
	  (*--- calculate the class constants ---*)
	  cc = Transpose@Table[	hma = hmat[[k]];
	     				   	vec = Diagonal[V.hma.Vi]
	     		  	,{k, 1, h}];
	  (*--- calculate the dimensions of the irreducible representations ---*)
	  dim = Sqrt[Table[g/(cc[[k]]/nc).Conjugate[cc[[k]]], {k, 1, h}]];
	  (*--- calculate the character table ---*)
	  ct = FullSimplify@Table[Table[cc[[i, j]]*dim[[i]]/nc[[j]], {j, 1, h}], {i, 1, h}];
	  (*--- sortby ---*)
	  cts = SortBy[Sort[ct], First];
	  posidentity = Position[cts, n_ /; (Norm[N[n - Table[1, {i, 1, Length[ct]}]], "Frobenius"] < 0.0001), {1}];	  
	  chi=Prepend[Delete[cts,posidentity],Table[1, {i, 1, Length[ct]}]];
	  
	  ,
	  (*--- For the group C1 ---*)
	  chi={{1}};
  ];

	
    (*Check if an own name of notations was specified*)
    If[notation == randint,
    	  (*Specified*)
    	  label = Table["",{i,1,Length[clm]}];
    	  If[Length[irepnotation]<Length[clm], label[[1;;Length[OptionValue[GOIrepNotation]]]] = OptionValue[GOIrepNotation], label = OptionValue[GOIrepNotation]];
    	  (*Not specified *)
         ,label = GTIrepNotation[{clm,chi},notation]];
	If[OptionValue[GOVerbose],
		clmprint = If[	inp<5,	
			Table[GTGetSymbol[clm[[i]]],{i,1,h}], 	
						(*Space group elements*)
						Table[Table[
							el=clm[[i,j]];
							pos = First@Flatten@Position[grpmat,el];
							 \[LeftAngleBracket] GTGetSymbol[grpmat[[pos]]],grp0[[pos,2]]\[RightAngleBracket]
							,{j,1,Length[clm[[i]]]}],{i,1,h}]];

		Liste = chi;
		Do[Liste[[i]] = Prepend[Liste[[i]], label[[i]]], {i, 1, h}];
		cla = Table[(*Subscript["C", i]*)Length[clm[[i]]] clmprint[[i,1]], {i, 1, h}];
		cla = Prepend[cla, ""];
		If[OptionValue[GOReality],
			Do[Liste[[i]] = Append[Liste[[i]], GTReality[clm, chi[[i]]]], {i, 1, h}];
			cla = Append[cla, "Reality"];
			
			,None];
		Pr = Chop[Prepend[Liste, cla]];
		If[OptionValue[GOReality],
		(*	Print[DisplayForm[GridBox[Pr, RowLines -> Flatten[{1, Table[0,{i,1,h-1}],1}], ColumnLines -> Flatten[{1, Table[0,{i,1,h-1}],1}], ColumnAlignments -> Center]]],
			Print[DisplayForm[GridBox[Pr, RowLines -> {1, 0}, ColumnLines -> {1, 0}, ColumnAlignments -> Right]]];*)
		   Print[Grid[Pr, Frame -> All, Dividers -> {{2 ->GTDividerColor1}, {2 -> GTDividerColor1}}, 
                 Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} ->GTCornerColor}]], 
           Print[Grid[Pr, Frame -> All, ColumnAlignments -> Center, 
                 Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} -> GTCornerColor}]]
	];
	Table[Print[Subscript["C", k], " = ", clmprint[[k]]], {k, 1, h}];,None];
	cl = If[	inp<5,	Table[GTWhichOutput[clm[[i]],inp],{i,1,h}],
						
						(*Space group elements*)
						Table[Table[
							el=clm[[i,j]];
							pos = First@Flatten@Position[grpmat,el];
							grp0[[pos]]
							,{j,1,Length[clm[[i]]]}],{i,1,h}]];
	{cl, chi,label}
 
  ]
(*
***)


(****n* /GTCharProjectionOperator
! NAME
!  GTCharProjectionOperator
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  07. Juni: New Version of the Character Projection operator
! USAGE
!  GTCharProjectionOperator[classes,characters,function,arguments] applies the character projection operator of a irreducible representation to a certain function with given arguments. 
! INPUT
!  list of classes, list of characters, function, arguments
! OUTPUT
!  function
! GTPack OPTIONS
!  GOFast, GOLattice
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTWhichInput, GTWhichOutput, GTGetMatrix, GTGroupQ, GTGetRotationMatrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
! 	The Character Projection operator works for multiple cases, illustrated below:
!
!	O(3)
!	----
!	scalar function			+
!	2d function				+
!	3d function				+
!	higher dimensional		-
!
!	SU(2)
!	----
!	scalar function			-
!	2d function	(spinor)	+
!	3d function				-
!	higher dimensional		-
!
!	O(2)
!	----
!	scalar function			+
!	2d function				+
!	3d function				-
!	higher dimensional		-
!
!	permutation matrices
!	--------------------
!	scalar function			-
!	2d function				-
!	3d function				-
!	higher dimensional		-
!
! LITERATURE
!
! TODO
!  Needs to be tested!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTCharProjectionOperator[class0_,char_,fu_,arg_,OptionsPattern[]] := Module[
	{grp,classes,classestmp,vecs,vecstmp,dimirep,run,h,tr,i,j,el,lattice},
	classestmp = Table[
		        Table[
		          el = class0[[i,j]];
		          If[GTWhichInput[el]<5, \[LeftAngleBracket] GTGetMatrix[el],  {0, 0, 0} \[RightAngleBracket],
		          	
		          	                     \[LeftAngleBracket] GTGetMatrix[el[[1]]], el[[2]] \[RightAngleBracket] ]
			    ,{j,1,Length[class0[[i]]]}]
			  ,{i,1,Length[class0]}];
	classes = classestmp[[;; , ;; , 1]];
	vecs = classestmp[[;; , ;; , 2]];
	If[OptionValue[GOFast],	run=True, run=GTGroupQ[Flatten[classes,1]]];
	If[run,
		grp = Flatten[classes,1];
		vecstmp = Flatten[vecs,1];
		lattice = OptionValue[GOLattice];
		vecs = Table[vecstmp[[i]].lattice,{i,1,Length[vecstmp]}];
		If[Length[char]>0&&Length[char[[1]]]>0,
			h = Length[char];
			dimirep = Table[Max[Abs[char[[i]]]],{i,1,h}];
			tr=Table[GTInternalCharProjectionOperator[classes,char[[i]],fu,grp,dimirep[[i]],arg,vecs],{i,1,h}], 
			
			dimirep = Max[Abs[char]];
			tr=GTInternalCharProjectionOperator[classes,char,fu,grp,dimirep,arg,vecs]],
			
		Print[$ErrorNoGroup];
		Abort[]
	 ];
	 Return[tr]
]


GTInternalCharProjectionOperator[classes_,char_,fu_,grp_,dimirep_,xyz_,vecs0_:{}] := Module[{g,pos,a,grpchar,dimfu,dimerror,tr,futmp,trtmp, rotgrp, arg, umat,i,vecs},
	    (*--- initial definitions ---*)
		g=Length[grp];
		vecs=If[Length[vecs0]==0,Table[{0,0,0},{i,1,g}],vecs0];
		pos=Table[Flatten[Position[classes,grp[[i]]]],{i,1,g}];
		grpchar=Table[char[[pos[[i,1]]]],{i,1,g}];
		
		(*--- Check the dimension of the function ---*)
		If[VectorQ[fu[a,a,a]],dimfu = Length[fu[a,a,a]],dimfu=0];
		dimerror = True;
		
		Which[
			
		(*---------------------------*)
		(*--- O(3) representation ---*)
		(*---------------------------*)
				grpdgrp==1, 
				Which[	dimfu==0,
								dimerror=False;
								(*ft[{x_,y_,z_}] = fu[x,y,z];*)
								(*tr=(dimirep/g)*Sum[Conjugate[grpchar[[i]]]*Apply[fu,Inverse[grp[[i]]].xyz],{i,1,g}]//Simplify,*)
								tr=(dimirep/g)*Sum[arg=Inverse[grp[[i]]].xyz+vecs[[i]];Conjugate[grpchar[[i]]]*fu[arg[[1]],arg[[2]],arg[[3]]],{i,1,g}]//Simplify,
						  	dimfu==2,
								dimerror=False;
								futmp[x_,y_,z_]={fu[x,y,z][[1]],fu[x,y,z][[2]],0};
								trtmp=(dimirep/g)*Sum[Conjugate[grpchar[[i]]]*grp[[i]].Apply[futmp,Inverse[grp[[i]]].xyz+vecs[[i]]],{i,1,g}];
								If[Simplify[trtmp[[3]]]!=0,Print["Warning: Projection is going out of plane!"],None];
								tr={trtmp[[1]],trtmp[[2]]}//Simplify,
							dimfu==3,
							    dimerror=False;
								tr=(dimirep/g)*Sum[Conjugate[grpchar[[i]]]*grp[[i]].Apply[fu,Inverse[grp[[i]]].xyz+vecs[[i]]],{i,1,g}]//Simplify
						  ],
		(*----------------------------*)
		(*--- SU(2) representation--- *)
		(*----------------------------*)
				grpdgrp==2,
					rotgrp = Table[GTGetRotationMatrix[grp[[i]]],{i,1,g}];
					Which[	dimfu==0,
								dimerror=False;
								Print[$ErrorProjectionNoSpinor];
								Return[0],
						  	dimfu==2,
								dimerror=False;
								tr=(dimirep/g)*Sum[Conjugate[grpchar[[i]]] grp[[i]].Apply[fu,Inverse[rotgrp[[i]]].xyz+vecs[[i]]],{i,1,g}]//Simplify,
							dimfu==3,
							    dimerror=False;
							    Print[$ErrorProjectionNoSpinor];
								Return[0]
						  ],
		(*---------------------------*)
		(*--- O(2) representation ---*)
		(*---------------------------*)
				grpdgrp==3,
					Which[	dimfu==0,
								dimerror=False;
								tr=(dimirep/g)*Sum[Conjugate[grpchar[[i]]] Apply[fu,Inverse[grp[[i]]].{x,y}],{i,1,g}]//Simplify,
						  	dimfu==2,
								dimerror=False;
								tr=(dimirep/g)*Sum[Conjugate[grpchar[[i]]] grp[[i]].Apply[fu,Inverse[grp[[i]]].{x,y}],{i,1,g}]//Simplify,
							dimfu==3,
							    dimerror=False;
							    Print[$ErrorProjectionNo2D];
								Abort[]
						  ],
		(*----------------------------*)
		(*--- Premutation matrices ---*)
		(*----------------------------*)
				grpdgrp==4,
					dimerror=False;
					Print[$ErrorProjectionNoPermutation];
					Abort[],
					
		(*------------------------------*)
		(*--- SU(2)xS representation ---*)
		(*------------------------------*)
				grpdgrp==5,
					rotgrp = Table[GTGetRotationMatrix[grp[[i]]],{i,1,g}];
					Which[	dimfu==0,
								dimerror=False;
								Print[$ErrorProjectionNoSpinor];
								Return[0],
						  	dimfu==2,
								dimerror=False;
								umat = Table[grp[[i,1;;2,1;;2]],{i,1,Length[grp]}];
								tr=(dimirep/g)*Sum[Conjugate[grpchar[[i]]] umat[[i]].Apply[fu,Inverse[rotgrp[[i]]].xyz+vecs[[i]]],{i,1,g}]//Simplify,
							dimfu==3,
							    dimerror=False;
							    Print[$ErrorProjectionNoSpinor];
								Return[0]
						  ]					
					];				

	If[dimerror,Print[$ErrorProjectionOperator],None];
	Return[Simplify[tr]]
]

(*
***)


(****g* /GTClassInfo
! NAME
!  GTClassInfo
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  test.m -> Auxiliary.m
! MODIFICATION HISTORY
!   October 2014: first version
! USAGE
!   gives a list of numbers to distinguish the elements of a class
! INPUT
!   character table
! OUTPUT
!  list
! GTPack OPTIONS
!
! STANDARD OPTIONS
!  
! GTPack MODULES
!   GTAxisOrder, GTGetMatrix, GTOrderOfElement, GTInverseElement
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed for the denotation of elements in the character table
! LITERATURE
!  -
! TODO
! -
! Release
!
! PROBLEMS
! - 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTClassInfo[ct_] := 
 Module[{ncl, orders, nelm, axespos, axespostmp, nap, ap, iee,
    ieepos, ergtmp, ergtmp2, pos, m},
  ncl = Length[ct[[1]]];
  orders = 
   Table[m = GTGetMatrix[ct[[1, i, 1]]]; 
    Det[m]*GTOrderOfElement[m], {i, 1, ncl}];
  nelm = Table[Length[ct[[1, i]]], {i, 1, ncl}];
  axespos = Table[0, {i, 1, ncl}];
  axespostmp = 
   Table[If[orders[[i]] == 2, GTAxisOrder[ct[[1, i, 1]]], 0], {i, 1, 
     ncl}];
  nap = 0;
  iee = Table[0, {i, 1, ncl}];
  If[Length[Position[ct[[1]], IEe]] > 0, 
   ieepos = First@Flatten[Position[ct[[1]], IEe]]; 
   iee[[ieepos]] = 1/1000];
  Do[If[axespostmp[[i]] > 0, nap = nap + 1], {i, 1, ncl}];
  While[Norm[axespostmp] > 0,
   ap = First@Flatten@Position[axespostmp, Max[axespostmp]];
   axespostmp[[ap]] = 0;
   axespos[[ap]] = nap;
   nap = nap - 1;
   ];
  ergtmp = orders + nelm/100 + I*axespos + iee;
  ergtmp2 = 
   If[Length[Position[ergtmp, -199/100]] > 1, 
    ergtmp[[First@Flatten[Position[ergtmp, -199/100]]]] = -199/100 + 
      1/10000; ergtmp, ergtmp];
  ergtmp2 = 
   If[Length[Position[ergtmp, -99/50]] > 1, 
    ergtmp[[First@Flatten[Position[ergtmp, -99/50]]]] = -99/50 + 
      1/10000; ergtmp, ergtmp];
  (*--- Ausnahme T ---*)
  
  ergtmp2 = 
   If[Length[Position[ergtmp, 76/25]] > 1, 
    ergtmp[[First@Flatten[Position[ergtmp, 76/25]]]] = 
     76/25 + 1/10000; ergtmp, ergtmp];
  (*--- Ausnahme C6v ---*)
  
  ergtmp2 = 
   If[Length[Position[ergtmp, -197/100]] > 1, 
    ergtmp[[First@Flatten[Position[ergtmp, -197/100]]]] = -197/100 + 
      1/10000; ergtmp, ergtmp];
  (*--- Ausnahme Th ---*)
  
  ergtmp2 = 
   If[Length[Position[ergtmp, -149/25]] > 1, 
    ergtmp[[First@Flatten[Position[ergtmp, -149/25]]]] = -149/25 + 
      1/10000; ergtmp, ergtmp];
  (*--- Ausnahme D6h ---*)
  
  ergtmp2 = 
   If[Length[Position[ergtmp, -299/50]] > 1, 
    ergtmp[[First@Flatten[Position[ergtmp, -299/50]]]] = -299/50 + 
      1/10000; ergtmp, ergtmp];
  Do[
   pos = Position[ergtmp2, ergtmp2[[i]]];
   If[Length[pos] > 1 && Length[ct[[1, i]]] == 1 && 
     GTInverseElement[ct[[1, i, 1]]] == ct[[1, pos[[2, 1]], 1]], 
    ergtmp2[[i]] = ergtmp2[[i]] + I/100], {i, 1, ncl}];
  Return[ergtmp2]
  ]
(*
***)



(****n* /GTClebschGordanCoefficients
! NAME
!  GTClebschGordanCoefficients
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  
! USAGE
!  GTClebschGordanCoefficients[Irep1, Irep2, Irep3] gives an array containing Clebsch-Gordan-Coefficients for the calculation of basis functions of Irep3 within a direct product representation of Irep1 and Irep2. 
! INPUT
!  matrices of representation 1, matrices of representation 2, matrices of representation 3
! OUTPUT
!  Clebsch-Gordan-Coefficients in matrix-form
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTGroupQ, GTDirectProductRep
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!
! LITERATURE
!
! TODO
!
! Release
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTClebschGordanCoefficients[rep10_, rep20_, rep3_,OptionsPattern[]] := Module[{rep1, rep2, d1, d2, d3, dp, g, Es, npq, dprep,A, CGK}, 
  If[OptionValue[GOFast]==False,
  							If[GTGroupQ[rep10],None,Print["Error: first representation is not a group"];Abort[]];
  							If[GTGroupQ[rep20],None,Print["Error: second representation is not a group"];Abort[]];
  							If[GTGroupQ[rep3],None,Print["Error: third representation is not a group"];Abort[]];

  							If[Length[rep10]==Length[rep20],None,Print["Error: group order of first and second representation does not match!"];Abort[]];
  							If[Length[rep10]==Length[rep3],None,Print["Error: group order of first and third representation does not match!"];Abort[]],
  							
  							None];
(*--- estimation of the length of the representations ---*)
  d1 = Length[rep10[[1]]]; 
  d2 = Length[rep20[[1]]]; 
  d3 = Length[rep3[[1]]]; 
  dp = d1*d2;
  If[d1 == 1, rep1 = Flatten[rep10, 1], rep1 = rep10]; 
  If[d2 == 1, rep2 = Flatten[rep20, 1], rep2 = rep20]; 
  g = Length[rep3]; 
  dprep = Simplify@GTDirectProductRep[rep1,rep2,GOFast->True];
  (*--- calculation of the projection operator ---*)
  A[l_, k_] := Sum[dprep[[i]]*Conjugate[rep3[[i, l, k]]], {i, 1, g}]*d3/g; 
  (*--- looking for an eigenvector with eigenvalue 1 ---*)
  Es = Eigensystem[Simplify@A[1, 1]]; 
  npq = Sum[Es[[1, i]], {i, 1, dp}]; 
  CGK = Simplify@Table[Table[A[i, 1].Es[[2, a]], {i, 1, d3}], {a, 1, npq}]; 
  Simplify@Map[Orthogonalize, CGK]]
(*
***)

(****n* /GTClebschGordanSum
! NAME
!  GTClebschGordanSum
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  December 2014: New Algorithm
! USAGE
!  GTClebschGordanSum[rep1,rep2] gives the direct sum of two representations rep1 and rep2.
! INPUT
!  list with several matrix-representations
! OUTPUT
!  representation matrices according to a Clebsch-Gordan-Sum
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! ERROR MESSAGES
!  
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  GTGroupQ
! DESCRIPTION
!
! LITERATURE
!
! TODO
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTClebschGordanSum[rep1_,rep2_, OptionsPattern[]] := Module[{g,ls,l1,l2, i,j, zero, Trzero,T,cgsum},
	If[OptionValue[GOFast]==False,
  							If[GTGroupQ[rep1],None,Print["Error: first representation is not a group"];Abort[]];
  							If[GTGroupQ[rep2],None,Print["Error: second representation is not a group"];Abort[]],
  							
  							None];
	
	If[Length[rep1]==Length[rep2],
		g=Length[rep1];	
		l1=Length[rep1[[1]]];
		l2=Length[rep2[[1]]];
		ls=l1+l2;
		
		zero = Table[0, {i, 1, l1}, {j, 1, l2}];
		Trzero = Transpose[zero];
		
		cgsum = Table[
			Flatten[{Table[Flatten@Append[rep1[[T, i]], zero[[i]]], {i, 1, l1}],Table[Flatten@Prepend[rep2[[T, i]], Trzero[[i]]], {i, 1, l2}]}, 1]
			, {T, 1, g}];
	 	Return[cgsum];
		,
	
		Print["Error: group order of representation 1 is not equal to group order of representation 2"];Abort[]]
 ]

HoldPattern[cg1_ \[CirclePlus] cg2_] := GTClebschGordanSum[cg1,cg2]
(*
***)

(****n* /GTClebschGordanTable
! NAME
!  GTClebschGordanTable
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  02/03/2014 development
! USAGE
!  GTClebschGordanTable[\[CapitalGamma]^1,\[CapitalGamma]^2,Clebsch-Gordan coefficients,names] illustrates a list of Clebsch-Gordan coefficients calculated for the direct product representation \[CapitalGamma]^1\[CircleTimes]\[CapitalGamma]^2.
! INPUT
!  matrices of representation 1, matrices of representation 2, list of Clebsch-Gordan coefficients, names of ireps
! OUTPUT
!  
! ERROR MESSAGES
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
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTClebschGordanTable[G1_, G2_, Inpcg_, Irep_] := Module[{l1, l2, listcg, IrepList, AlphaList, Numberlist, list2,   empty, mlist, nlist, list3},
  l1 = Length[G1[[1]]];
  l2 = Length[G2[[1]]];
  listcg = Transpose[Flatten[Simplify@Inpcg, 2]];
  IrepList = Flatten@Table[Table[Table[If[k == 1, Irep[[i]], ""], {k, 1, Length[Inpcg[[i, j]]]}], {j,1, Length[Inpcg[[i]]]}], {i, 1, Length[Inpcg]}];
  AlphaList = Flatten@Table[Table[Table[If[k == 1, j, ""], {k, 1, Length[Inpcg[[i, j]]]}], {j, 1,Length[Inpcg[[i]]]}], {i, 1, Length[Inpcg]}];
  Numberlist = Flatten@Table[Table[Table[k, {k, 1, Length[Inpcg[[i, j]]]}], {j, 1,Length[Inpcg[[i]]]}], {i, 1, Length[Inpcg]}];
  list2 = Prepend[Prepend[Prepend[listcg, Numberlist], AlphaList],IrepList];
  empty = {"", ""};
  mlist = Flatten@Prepend[Prepend[Table[Table[i, {j, 1, l2}], {i, 1, l1}], "j"], empty];
  nlist = Flatten@Prepend[Prepend[Table[Table[j, {j, 1, l2}], {i, 1, l1}], "k"], empty];
  list3 = Table[Prepend[Prepend[list2[[i]], nlist[[i]]], mlist[[i]]], {i, 1,Length[list2]}];
  (*Print@Grid[list3,Dividers -> {{None, None, Black}, {None, None, None, Black}}]*)
  Print@Grid[list3,Frame -> All,Dividers -> {{{2 ->GTDividerColor1}, {2 ->GTDividerColor1}, {2 ->GTDividerColor1}}, {{2 ->GTDividerColor1}, {2 ->GTDividerColor1}, {2 ->GTDividerColor1}, {2 ->GTDividerColor1}}},
  	Background -> {{1 -> GTBackGroundColor1,2 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1,2 -> GTBackGroundColor1,3 -> GTBackGroundColor1},{ 
  		{1, 1} ->GTCornerColor, {1, 2} ->GTCornerColor, {2, 1} ->GTCornerColor, {2, 2} ->GTCornerColor}}]
  ]
(*
***)

(****n* /GTDirectProductChars
! NAME
!  GTDirectProductChars
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  
! USAGE
!  GTDirectProductChars[characters1,characters2] gives the character system of a direct product representation from the character systems of two given representations.
! INPUT
!  characters of representation 1, characters of representation 2
! OUTPUT
!  list with characters of the direct product group
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! ERROR MESSAGES
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
! Release
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)
			   
GTDirectProductChars[ch1_,ch2_] := Module[{ch},
	If[Length[ch1]==Length[ch2],ch=ch1*ch2
		,Print["Error: characters don't fit"]]]
(*
***)

(****n* /GTCosetRepresentative
! NAME
!  GTCosetRepresentative
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m
! MODIFICATION HISTORY
!   6/04/2020: first version
! USAGE
!   GTCosetRepresentative determines a coset representative for a left coset decomposition a group into a sub group of index 2 or 3
! INPUT
!   group
!   a subgroup of index 2 or 3
! OUTPUT
!  a group element q which serves as coset representative according to G = H + q*H (+ q^2*H)
!  permutation for the map between G and G = H + q H (+ q^2 H)
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTOrderOfElement, GTgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTSGGetIreps to establish a subgroup chain with sub groups of index 2 or 3
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
GTCosetRepresentative[grp_, grpn_] := 
 Module[{n, gln, comp, qel, reordgroup, perm, qelsq, Tm, groupfound},
  (*Calculate Index and Length of group*)
  n = Length[grp]/Length[grpn];
  gln = Length[grpn];
  (*Estimate the complement of group G and subgroup H to write G=H+qH (+q^2 H);*)
  comp = Complement[grp, grpn];
  reordgroup = {};
  groupfound = False;
  Tm = 0;
  While[Not[groupfound] && Tm < gln, Tm = Tm + 1;
   qel = comp[[Tm]];
   qelsq = GTgmat[qel, qel];
   (*write G=H+qH+(q^2H);*)
   reordgroup = Which[
     (*Index 2*)
     n == 2, 
     Flatten[{grpn, Table[GTgmat[qel, grpn[[T]]], {T, 1, gln}]},1],
     (*Index 3*)
     n == 3, 
     Flatten[{grpn, Table[GTgmat[qel, grpn[[T]]], {T, 1, gln}], 
       Table[GTgmat[qelsq, grpn[[T]]], {T, 1, gln}]},1]];
   groupfound = Length[Complement[reordgroup, grp]] === 0;];
  If[Tm > gln, 
   Print["Error: GTCosetRepresentative: Could not identify coset representative!"]; Abort[]];
  (*Find the corresponding Permutation of the reordered group*)
  
  perm = FindPermutation[reordgroup, grp];
  Return[{qel, perm}];]
(*
***)

(****n* /GTDirectProductRep
! NAME
!  GTDirectProductRep
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  
! USAGE
!  GTDirectProductRep[rep1,rep2] calculates the direct product representation of two representations rep1 and rep2.
! INPUT
!  representation 1, representation2
! OUTPUT
!  direct product representation
! GTPack OPTIONS
!  GOFast, GODiagonal
! STANDARD OPTIONS
!
! ERROR MESSAGES
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
!--------------------------------------------------------------------------------
! SOURCE
*)
			   
GTDirectProductRep[rep1_,rep2_,OptionsPattern[]] := Module[{},
	If[OptionValue[GOFast]==False,
  						If[GTGroupQ[rep1],None,Print["Error: first representation is not a group"];Abort[]];
						If[GTGroupQ[rep2],None,Print["Error: second representation is not a group"];Abort[]];
						If[OptionValue[GODiagonal]==True,					
  						If[Length[rep1]==Length[rep2],None,Print["Error: the orders of the representations do not match."];Abort[]]],
  						None];
	
 
	If[OptionValue[GODiagonal]==True,	
		Return[Table[KroneckerProduct[rep1[[i]],rep2[[i]]],{i,1,Length[rep1]}]],
        Return[Flatten[Table[KroneckerProduct[rep1[[i]], rep2[[j]]], {i, 1, Length[rep1]}, {j, 1,Length[rep2]}], 1]]]
]

HoldPattern[rep1_ \[CircleTimes] rep2_] := GTDirectProductRep[rep1,rep2]
(*
***)

(****n* /GTExtraRepresentations
! NAME
!  GTExtraRepresentations
! AUTHOR
!  M.Geilhufe, W. Hergert
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  03/21/2016 : first implementation
!  24.02.2017 : listes of names of extrarepresentations can be provided.
! USAGE
!  GTExtraRepresentations[double group] extracts the extra representations from the 
!  character table of a double group
! INPUT
!  character table of the double group
! OUTPUT
!  extra representation of the double group
! GTPack OPTIONS
!  GOFast, GOIrepNotation, GOVerbose
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTGroupQ
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

  GTExtraRepresentations[ct_, OptionsPattern[]] := Module[
    {names,fast,verb,cl,grp,deepos,extrachars,extranames,i},
  (*--- options ---*)
     names = OptionValue[GOIrepNotation];
     fast  = OptionValue[GOFast];
     verb  = OptionValue[GOVerbose];
     cl    = ct[[1]];
     grp   = Flatten[cl, 1];
  (*--- test input ---*)   
     If[fast,
        None,
        If[Not[grpdgrp == 2 || grpdgrp == 5],
           Print[$ErrorGTSOCSplitting]; Abort[];
        ];
        If[Not[GTGroupQ[grp]],
           Print[$ErrorNoGroup]; Abort[]
        ];
     ];
     deepos     = First@First@Position[cl, DEe];
     extrachars = {};
     extranames = {};
     Do[
        If[ct[[2, i, deepos]] > 0,
           None,
           extrachars = Append[extrachars, ct[[2, i]]];
           extranames = Append[extranames, ct[[3, i]]]]
     , {i, 1, Length[cl]}];
  (*--- provide list of names ---*)
     If[names == {},
        None ,
        extranames = names
     ];      
     If[verb, 
        Print[Grid[
              Prepend[Table[Prepend[extrachars[[i]], extranames[[i]]], {i, 1, Length[extrachars]}], 
              Flatten@{"extra-representations", cl[[;; , 1]]}], 
              Frame -> All,Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
              Background -> {{1 -> GTBackGroundColor1}, {1 -> 
              GTBackGroundColor1}, {1, 1} -> GTCornerColor}]
        ];, 
        None
     ];
     Return[{cl, extrachars, extranames}]
  ]
(*
***)

(****n* /GTGetIreps
! NAME
!  GTSGGetIreps
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m
! MODIFICATION HISTORY
!   12/06/2018: first version
! USAGE
!   GTGetIreps gives the character table and the representation matrices of a point group.
! INPUT
!   group
! OUTPUT
!   {{classes, characters, names of ireps},{G,representation matrices }}
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTPrintGeneral, GTPrintChars, GTPrintIreps, GTgmat, GTClasses, GTGetInvSubGroup, GTCosetRepresentative, GTGetIreps, GTInverseElement, GTInduceOrbitLengthOne, GTInduceOrbitLengthTwoThree
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
!  Aroyo, Mois I., et al. "Bilbao Crystallographic Server. II. Representations of crystallographic point groups and space groups." Acta Crystallographica Section A 62.2 (2006): 115-128.
! TODO
! -
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTGetIreps[grp_, OptionsPattern[]] := 
 Module[{ireps, charsfinal, classes, grpn, irepsn, qel, perm, qelsq, 
   qelqu, qelinv, qelsqinv, gl, gln, orbitq, pos, orbitqsq, irep, 
   irepq, irepqsq, lengthone, chars, charscompinitial, charscomp, irepsfinal, clt, el,
   randint,irepnotation,run,notation0,notation, label},
  
  (**********************************)
  (* Preliminary Consideration      *)
  (**********************************) 
  (*Print a warning in case magnetic groups are used.*)
  If[GTMagneticQ[], "Warning: You are using magnetic groups. These groups contain anti-unitary elements and need a special treatment concerning the irreducible representations. This feature is not implemented yet and the character table you obtain might not be the one you were aiming for. We are currently working on that!", None];
  (*It is possible to specify a list of irep names. A random number will be used to generate an internal string to distinguish this case from "Bethe",...*)
  randint = "list"<>ToString[RandomInteger[50000]];	
  irepnotation = If[Length[OptionValue[GOIrepNotation]]==0,OptionValue[GOIrepNotation],randint];
  
  (*--- Check the input ---*)
  If[OptionValue[GOFast],run=True,run=GTGroupQ[grp]];
  If[Not[run],Print[$ErrorNoGroup];Abort[]];	

  Which[	grpdgrp==1, (*--- O(3) ---*)
				notation0 = irepnotation,
			grpdgrp==2, (*--- SU(2) ---*)
			    notation0 = If[irepnotation==randint,randint,"Bethe"];
				If[OptionValue[GOVerbose],
					If[irepnotation=="Mulliken",
						Print["Warning: The standard representation is given by SU(2) matrices. Mulliken symbols can not be evaluated for double groups. The notation is switched to Bethe notation."]
						,None];
					If[irepnotation=="Bouckaert",
						Print["Warning: The standard representation is given by SU(2) matrices. Bouckaert notation can not be evaluated for double groups. The notation is switched to Bethe notation."]
						,None]
					,None],
			grpdgrp==3, (*--- O(2) ---*)
				notation0 = irepnotation,
			grpdgrp==4, (* ---Own ---*)
				notation0 = irepnotation,
			grpdgrp==5, (*--- SU(2)xS ---*)
			    notation0 = If[irepnotation==randint,randint,"Bethe"];
				If[OptionValue[GOVerbose],
					If[irepnotation=="Mulliken",
						Print["Warning: The standard representation is given by SU(2)xS matrices. Mulliken symbols can not be evaluated for double groups. The notation is switched to Bethe notation."]
						,None];
					If[irepnotation=="Bouckaert",
						Print["Warning: The standard representation is given by SU(2)xS matrices. Bouckaert notation can not be evaluated for double groups. The notation is switched to Bethe notation."]
						,None]
					,None]];

  notation = If[notation0=="Bethe"||notation0=="Mulliken"||notation0=="Bouckaert"||notation0==randint,notation0,Print["Warning: Notation unknown. The notation is switched to Bethe notation."];"Bethe"];

   
  (**********************************)
  (*Irreducible representation of C1*)
  (**********************************)
  If[Length[grp] == 1,
   ireps = {{{{1}}}};
   charsfinal = {{1}};
   classes = {grp};
   
   If[notation == randint,
    	  (*Specified*)
    	  label = Table["",{i,1,Length[classes]}];
    	  If[Length[irepnotation]<Length[classes], label[[1;;Length[OptionValue[GOIrepNotation]]]] = OptionValue[GOIrepNotation], label = OptionValue[GOIrepNotation]];
    	  (*Not specified *)
         ,label = GTIrepNotation[{classes,charsfinal},notation]];
  
   If[OptionValue[GOVerbose], GTPrintIreps[grp, ireps,label]];
   If[OptionValue[GOVerbose], GTPrintChars[grp, classes, charsfinal,label]];
  
   Return[{{classes, charsfinal, label}, {grp, ireps}}]
  ];
   
   (*If[OptionValue[GOVerbose], GTPrintIreps[grp, irepsfinal, grp]];
   If[OptionValue[GOVerbose], GTPrintChars[grp, classes, charsfinal]];
  
  
   Return[{{classes, charsfinal, {Superscript["\[CapitalGamma]", ToString[1]]}}, {grp, ireps}}] ];*)
  (*****************************************)
  (*Irreducible representation of other SGs*)
  (*****************************************)
  (*check if divisible by 3 or 2;  3 is preferred to keep the group order small*)
  
  n = If[Divisible[Length[grp], 2], 2, 3];
  classes = GTClasses[grp, GOFast -> True];
  
  grpn = GTGetInvSubGroup[grp, classes, n, GOVerbose -> False];
  If[Length[grpn] == 0 && n == 2, n = 3; grpn = GTGetInvSubGroup[grp, classes, n, GOVerbose -> False]];
  (*Calculate ireps of the subgroup*)
  
  irepsn = Last@Last@GTGetIreps[grpn, GOVerbose -> False];
  (*Calculate Coset decomposition*)
  {qel, perm} = GTCosetRepresentative[grp, grpn];
  qelsq = GTgmat[qel, qel];
  qelqu = GTgmat[qel, qelsq];
  qelinv = GTInverseElement[qel];
  qelsqinv = GTInverseElement[qelsq];
  
  (*----------------------*)
  (*Induce representations*)
  (*----------------------*)

  (*Calculate orbits*)
  gl = Length[grp];
  gln = Length[grpn];
  orbitq = Table[el = GTgmat[qel, GTgmat[grpn[[T]], qelinv]];
                 First@First@Position[grpn, el], {T, 1, gln}];
  orbitqsq = Table[el = GTgmat[qelsq, GTgmat[grpn[[T]], qelsqinv]];
                   First@First@Position[grpn, el], {T, 1, gln}];
  (*Scan through all ireps, check length of orbits and induce representations*)
  (*Scan through all ireps,check length of orbits and induce representations*)

    ireps = FullSimplify[Flatten[Table[
      irep = irepsn[[h]];
      irepq = Table[irep[[orbitq[[T]]]], {T, 1, gln}];
      irepqsq = Table[irep[[orbitqsq[[T]]]], {T, 1, gln}];
      lengthone = 
       If[n == 2, Length[Union[{Map[Tr, irep], Map[Tr, irepq]}]] == 1,
         Length[Union[{Map[Tr, irep], Map[Tr, irepq], 
            Map[Tr, irepqsq]}]] == 1];
      If[lengthone,
       (*Length of orbit is 1*)
       GTInduceOrbitLengthOne[grp, grpn, qel, qelsq, qelqu, perm, irep, irepq],
       (*Length of orbit is 2 or 3*)
       {GTInduceOrbitLengthTwoThree[grp, grpn, qel, qelsq, qelqu, perm, irep,irepq, irepqsq]}]
      , {h, 1, Length[irepsn]}], 1]];
  
  chars = Table[Table[Tr[ireps[[h, T]]], {T, 1, gl}], {h, 1,Length[ireps]}]//FullSimplify;
  (*--- only keep different representations and order similarly to GTCharacterTabe ---*)
  charscompinitial = SortBy[Sort[Union[chars]],First];
  charscomp = Prepend[Delete[charscompinitial,First@Flatten@Position[charscompinitial, Table[1, {i, 1, Length[charscompinitial[[1]]]}]]],Table[1, {i, 1, Length[charscompinitial[[1]]]}]];
  
  
  (*Final Ireps*)
  irepsfinal = Table[pos = First@First@Position[chars, charscomp[[i]]];
                     ireps[[pos]], {i, 1, Length[charscomp]}];
  (*General position ireps*)
  (*Final chars little group*)
  
  clt = Transpose[charscomp];
  charsfinal = 
   FullSimplify[
    Transpose[
     Table[pos = First@First@Position[grp, classes[[i, 1]]]; 
      clt[[pos]], {i, 1, Length[classes]}]]];
  
  (*----------------------*)
  (* Print Output         *)
  (*----------------------*)
  (*Check if an own name of notations was specified*)
    If[notation == randint,
    	  (*Specified*)
    	  label = Table["",{i,1,Length[classes]}];
    	  If[Length[irepnotation]<Length[classes], label[[1;;Length[OptionValue[GOIrepNotation]]]] = OptionValue[GOIrepNotation], label = OptionValue[GOIrepNotation]];
    	  (*Not specified *)
         ,label = GTIrepNotation[{classes,charsfinal},notation]];
  
  If[OptionValue[GOVerbose], GTPrintIreps[grp, irepsfinal,label]];
  If[OptionValue[GOVerbose], GTPrintChars[grp, classes, charsfinal,label]];
  
  Return[{{classes, charsfinal, label}, {grp, irepsfinal}}]
  ]

(*
***)

(****n* /GTGetIrepEQ1
! NAME
!  GTGetIrepEQ1
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  
! USAGE
!  gives the representation matrices for one dimensional irreducible representations. The representation matrix is given by the character. The module is needed
!  by GTGetIrep.
! INPUT
!  group, character table, index
! OUTPUT
!  list of representation matrices
! GTPack OPTIONS
!
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTGetSymbol
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!
! LITERATURE
!
! TODO
!  
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTGetIrepEQ1[grp_,ct_,ind_]:=Module[{g, pos,tab,i},
	g = Length[grp];
	pos=Table[Flatten[Position[ct[[1]],grp[[i]]]],{i,1,g}];
	tab=Table[{{ct[[2,ind,pos[[i,1]]]]}},{i,1,g}];
	Return[tab]
]

(*
***)

(****n* /GTGetIrepCornwell
! NAME
!  GTGetIrepCornwell
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  23/09/2015 development of command
!  17/09/2016 Transpose the matrices!
! USAGE
!  gives the representation matrices for irreducible representations with dimensions greater than one. The module is needed
!  by GTGetIrep.
! INPUT
!  group, character table, index
! OUTPUT
!  list of representation matrices
! ERROR MESSAGES
!  
! GTPack MODULES
!  GTIrep, GTAngularMomentumChars, GTCharProjectionOperator, GTCartesianTesseralHarmonicsY, GTTransformationOperator, GTOrderOfElement, GTMakeEigenvectors
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  The module applies the character projection operator to cartesian tesseral harmonics to find an sufficient number of linear independent basis functions. 
!  By transforming the basis functions according to a certain symmetry element, the representation matrices can be found.
! LITERATURE
!	See Cornwell, Academic Press, 1984
! TODO
!  
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTGetIrepCornwell[grp_, ct_, ind_, lmax_] := Module[{idim, classes, chars, names, functions, l, m, i,T,B, ireps, f, g,coeff},
  {classes, chars, names} = ct;
  idim = chars[[ind, 1]];
  (*--- Try to find Basis functions by applying the character projection operator to cartesian tesseral harmonics ---*)
  functions = {};
  l = 0;
  ireps = Table[0, {i, 1, Length[classes]}];
  While[ireps[[ind]] == 0,
          ireps = GTIrep[GTAngularMomentumChars[classes, l, GOFast -> True], ct,GOFast -> True, GOVerbose -> False];
          l = l + 1;
          If[l > lmax, Print[$ErrorGTGetIrepCornwell];Abort[]]
   ];
  l = l - 1;
  Table[
  	f[x_, y_, z_] = Simplify[GTCartesianTesseralHarmonicY[l, m, x, y, z]];
    g = Simplify@GTCharProjectionOperator[classes, chars[[ind]], f, {x, y, z}];
    If[Not[g === 0], functions = Append[functions, g]], {m, -l, l}];
  If[Length[functions]>idim,Print[$ErrorGTGetIrepCornwellMethod];Abort[]];
  (*--- Transform each basis function with each group element to obtain a linear equation system for the representation matrices. ---*)
  coeff = Table[B[i], {i, 1, idim}];
  Table[
   Transpose@Table[
    f[x_, y_, z_] = functions[[i]];
    g = GTTransformationOperator[grp[[T]], f, {x, y, z}];
    coeff /. Flatten[Solve[Table[g == coeff.functions /.x -> RandomInteger[100000]/(RandomInteger[99999] + 1) /.y -> RandomInteger[100000]/(RandomInteger[99999] + 1) /.z -> RandomInteger[100000]/(RandomInteger[99999] + 1),{j, 1,idim}], coeff]]
    , {i, 1, idim}],
  {T, 1, Length[grp]}]
 ]

(*
***)
(****n* /GTGetIrepGT1
! NAME
!  GTGetIrepGT1
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  06/06/2013 A new version of GTGetIrepGT1 was implemented, to improve the clarity of the module. A development notebook GTGetIrepGT1 exists
!  02/02/2014 new version
!  26/08/2015 conjugate of character in character projection operator, remove errors for double groups
!  08/03/2021 Command generalized. It now works for any reducible representation to start with (Flodmark -> use regular representation; Cornwell use real-valued anguar momentum representation)
! USAGE
!  gives the representation matrices for irreducible representations with dimensions greater than one. The module is needed
!  by GTGetIrep.
! INPUT
!  group, character table, index, matrix representation for reduction (regorig)
! OUTPUT
!  list of representation matrices
! GTPACK OPTIONS
!
! STANDARD OPTIONS
!
! ERROR MESSAGES
!  
! GTPack MODULES
!  GTGetMatrix, GTOrderOfElement, GTMakeEigenvectors
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!
! LITERATURE
!	See E. Blokker, A Theory for the Construction of the Irreducible Representations of Finite Groups .
! TODO
!  
! RELEASE
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTGetIrepGT1[grp_, ct_, ind_,regorig_] := Module[{lg, dp, Smat, reg, es, ev, vv,classes,charstmp,names,chars,h,grppos},
  {classes, charstmp, names} = ct;
  chars = charstmp[[ind]];
  (*--- Group order and dimension of irreducible representation ---*)
  lg = Length[grp];
  dp = Max[chars];
  h = Length[classes];
  (*--- regular representation and sort by order of element ---*)
  reg = SortBy[regorig, GTOrderOfElement];
  (*--- Char Projection operator ---*)
  Smat = Flatten[
  	Sum[
  		Conjugate[chars[[i]]]*Sum[
       					grppos = Flatten[Position[grp, classes[[i, j]]]];
       					regorig[[grppos]]
       				, {j, 1, Length[classes[[i]]]}]
     , {i, 1, h}]*dp/lg, 1];  
  (*--- Find an eigenvector with non degenerate eigenvalue and apply other elements until dp linear independent vectors were found ---*)
  Do[
   es = Eigensystem[reg[[i]]];
   Do[
    ev = Smat.es[[2, j]];
    If[Norm[ev] == 0, Continue[]];
    vv = GTMakeEigenvectors[ev, reg];
    If[Length[vv] == dp, Break[], None],
    {j, 1, Length[es[[2]]]}],
   {i, lg, lg - 1, -1}];
  vv = Orthogonalize[RowReduce[vv][[1 ;; dp]]];
  (*--- calculate representation matrices ---*)
  Simplify@Table[vv.regorig[[i]].ConjugateTranspose[vv], {i, 1, lg}]]

(*
***)

(****n* /GTGetIrep
! NAME
!  GTGetIrep
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  15. August 2014 (WH)
! USAGE
!  GTGetIrep[group,character table,index] gives the representation matrices of an irreducible representation (denoted by its index within the character table).
! INPUT
!  group, character table, (index)
! OUTPUT
!  list representation matrices
! GTPack OPTIONS
!  GOFast, GOMethod, GOlmax
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTChangeRepresentation, GTGetIrepEQ1, GTGetIrepGT1, GTGetIrepCornwell, GTGroupQ, GTGetMatrix
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!
! LITERATURE
!
! TODO
!   There was a problem with the representation Matrix <--> Symbol .
!   Now it works but it has to be checked  --> Matthias
! RELEASE
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTGetIrep[grp_,ind_,ct_:{}, OptionsPattern[]] := Module[{eo,grpt,Irep, cti,grptmp,ctind,grpind,irepsind,posind,regorig, ll, cond, irs, classes,dl,rep},
   (*Currently the Cornwell method is only implemented for ordinary groups*)
   (*in the case of an O(2) representation we change to O(3) for convenience*)
   If[OptionValue[GOMethod]=="Cornwell",
	  		   	grptmp=grpdgrp;
   				Which[grpdgrp==1, None,
         			  grpdgrp==3, GTChangeRepresentation["O(3)",GOVerbose->False]
         			  (*,
         			  grpdgrp==2||grpdgrp>3, Print[$ErrorGTGetIrepRepr];Abort[]*)
         			  ];
         	    grpt=GTGetSymbol[grp]];
         			 
	
   If[OptionValue[GOFast],
   	None,
   	If[GTGroupQ[grp],
   		None,
   		Print[$ErrorNoGroup];
   		Abort[]];
   ];

   (*--- character table in matrix form ---*)
   grpt = GTGetMatrix[grp]; 
   If[Length[ct]==0,
   	cti=GTCharacterTable[grpt,GOFast->True,GOVerbose->False],
   	
   	cti=ct;
   	cti[[1]]=Map[GTGetMatrix[#]&,cti[[1]]];
   	];
   If[ind>Length[cti[[1]]],Print[$ErrorNumberIreps];Abort[]]; 
   If[cti[[2,ind,1]]==1,
	  Irep=GTGetIrepEQ1[grpt,cti,ind],
	  
	  Which[OptionValue[GOMethod]=="Flodmark",
	  	               (*Use regular representation*)
	  	               regorig = GTRegularRepresentation[grp, GOFast -> True];
	  		           Irep=GTGetIrepGT1[grpt,cti,ind,regorig],
	  		OptionValue[GOMethod]=="Cornwell",
	  		           (*Use real-valued angular momentum representation*)
	  	              
	  		           (*Determine optimal l for angular momentum representation*)
	  		           ll = 0;
                       dl=If[grpdgrp==1||grpdgrp==3,1,1/2];
                       cond = False;
                       classes = Map[GTGetSymbol,cti[[1]]];
                       While[ll < OptionValue[GOlmax] && Not[cond],
                             ll += dl;
                             eo=-1;
                             While[eo<1 && Not[cond],
                             	eo += 1;
                             	irs = GTIrep[GTAngularMomentumChars[classes, ll, GOEvenOdd->eo,GOFast -> True], cti, GOFast -> True, GOVerbose -> False];
                             	cond = irs[[ind]] > 0;
                             ];
                            ];
                       If[Not[cond],Print[$ErrorGTGetIrepCornwell];Abort[]];
                       rep=If[IntegerQ[ll],"Real","Complex"];
                       eo=If[IntegerQ[ll],eo,-1000];
                       regorig = GTAngularMomentumRep[GTGetSymbol[grpt], ll, GOHarmonics -> rep, GOEvenOdd->eo, GOFast -> True];
                       Irep=GTGetIrepGT1[grpt,cti,ind,regorig],
	  		           (*Irep=GTGetIrepCornwell[grpt,cti,ind,OptionValue[GOlmax]],*)
	  		OptionValue[GOMethod]=="Induction",
	  		            {ctind, {grpind, irepsind}} = GTGetIreps[grpt,GOFast->True,GOVerbose->False];
	  		            posind=Position[ctind[[2]],cti[[2,ind]]];
	  		            If[Length[posind]==0,Print[$ErrorFindIrepInCT];Abort[]]; 
	  		            Irep=irepsind[[First@First@posind]]
	  		            ,
	  		Not[OptionValue[GOMethod]=="Cornwell"]&&Not[OptionValue[GOMethod]=="Flodmark"]&&Not[OptionValue[GOMethod]=="Induction"],
	  		           Print[$ErrorGTGetIrep];
	  		           Abort[];]
	  		
   ];
   
   If[grptmp==3, GTChangeRepresentation["O(2)",GOVerbose->False]];
   Return[Irep]
]

(*
***)

(****n* /GTInduceOrbitLengthOne
! NAME
!  GTInduceOrbitLengthOne
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m
! MODIFICATION HISTORY
!   04/04/2020: first version
! USAGE
!   GTInduceOrbitLengthOne induces two (three) irreducible representations from an irreducible representation of orbit length 1
! INPUT
!   group, subgroup of index n, coset representative (G/T = H + q*H (+q^2*H)), permutation between initial G/T and G/T given by coset decomposition,
!   the irreducible representation
! OUTPUT
!   induced representations
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTGetIreps for the induction of 2 (3) irreducible representations from an irep with an orbit of length 1
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

GTInduceOrbitLengthOne[grp_, grpn_, qel_, qelsq_, qelqu_, perm_,irepsn_, irepsq_] := Module[{n, gln, posbel, Amat, Bmat, Cmat, Um, Umsq, erg, ind, dim, elfound},
  (*Calculate Index and Length of group*)
  
  n = Length[grp]/Length[grpn];
  gln = Length[grpn];
  (*--------------------*)
  (*Generate the matrix U*)
  (*--------------------*)
  (*Pick an element of the original group,here, the last one and calculate the conjugate element with q*)
  dim = Length[irepsn[[1]]];
  posbel = 0;
  elfound = False;
  While[Not[elfound] && posbel < gln, 
   posbel = posbel + 1;
   elfound = Not[irepsn[[posbel]] === irepsq[[posbel]]] && Not[irepsn[[posbel]] === IdentityMatrix[dim]] && Not[irepsq[[posbel]] === IdentityMatrix[dim]]];
  (*Determine U from D(q b q^-1)=U D(b) U^-1*)
  
  Amat = irepsn[[posbel]];
  Bmat = irepsq[[posbel]];
  (*Determine U*)
  
  
  
  
  Which[
   (*Index 2*)
   n == 2,
   Cmat = irepsn[[First@First@Position[grpn, qelsq]]];
   Um = GTGetUnitaryTrafoMatrix[Amat, Bmat, Cmat, 2];,
   (*Index 3*)
   n == 3,
   Cmat = irepsn[[First@First@Position[grpn, qelqu]]];
   Um = GTGetUnitaryTrafoMatrix[Amat, Bmat, Cmat, 3];];
  Umsq = Um.Um;
  Which[
   (*Index 2*)
   n == 2, 
   ind = Table[Flatten[{irepsn, Table[Exp[i I \[Pi]]*Um.irepsn[[T]], {T, 1, gln}]}, 1], {i, 0,1}];,
   (*Index 3*)
   n == 3, 
   ind = Table[Flatten[{irepsn, Table[Exp[-2 \[Pi] I i/3]*Um.irepsn[[T]], {T, 1, gln}], Table[Exp[-4 \[Pi] I i/3]*Umsq.irepsn[[T]], {T, 1, gln}]}, 1], {i, 0, 2}];];
  erg = Table[Permute[ind[[i]], perm], {i, 1, n}];
  Return[erg]]

(*
***)

(****z* /GTInduceOrbitLengthTwoThree
! NAME
!  GTInduceOrbitLengthOne
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m
! MODIFICATION HISTORY
!   04/04/2020: first version
! USAGE
!   GTInduceOrbitLengthTwoThree induces an irreducible representations from an irreducible representation of orbit length 2 or 3
! INPUT
!   group, subgroup of index n, coset representative q (G/T = H + q*H (+q^2*H)), q^2, q^3, permutation between initial G/T and G/T given by coset decomposition,
!   the irreducible representation, the irreducible representation conjugated with q, the irreducible representation conjugated with q^2, real space basis
!   kvector
! OUTPUT
!   induced representation
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTgmat
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
GTInduceOrbitLengthTwoThree[grp_, grpn_, qel_, qelsq_, qelqu_, perm_, irep_, irepq_, irepqsq_] := Module[{n, gln, grpnrep, pos, qmat, ind, rep, qmatsq},
  n = Length[grp]/Length[grpn];
  gln = Length[grpn];
  Which[
   (*index 2*)
   n == 2,
   grpnrep = Table[KroneckerProduct[{{1, 0}, {0, 0}}, irep[[T]]] + KroneckerProduct[{{0, 0}, {0, 1}}, irepq[[T]]], {T, 1, gln}];
   pos = First@First@Position[grpn, qelsq];
   qmat = KroneckerProduct[{{0, 1}, {0, 0}}, IdentityMatrix[Length[irep[[pos]]]]] + KroneckerProduct[{{0, 0}, {1, 0}}, irep[[pos]]];
   ind = Flatten[{grpnrep, Table[qmat.grpnrep[[T]], {T, 1, gln}]}, 1],
   (*index 3*)
   n == 3, 
   grpnrep = Table[KroneckerProduct[{{1, 0, 0}, {0, 0, 0}, {0, 0, 0}}, irep[[T]]] + KroneckerProduct[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}, irepq[[T]]] + KroneckerProduct[{{0, 0, 0}, {0, 0, 0}, {0, 0, 1}}, irepqsq[[T]]], {T, 1, gln}];
   pos = First@First@Position[grpn, qelqu];
   qmat = KroneckerProduct[{{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},IdentityMatrix[Length[irep[[pos]]]]] + KroneckerProduct[{{0, 0, 0}, {0, 0, 0}, {1, 0, 0}}, irep[[pos]]];
   qmatsq = qmat.qmat // Simplify;
   ind = Flatten[{grpnrep, Table[qmat.grpnrep[[T]], {T, 1, gln}], Table[qmatsq.grpnrep[[T]], {T, 1, gln}]}, 1]];
  rep = Permute[ind, perm];
  Return[rep];]

(*
***)


(****n* /GTIrep
! NAME
!  GTIrep
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  9/3/2017: The was a conjugate missing in the equation!!!!!!
!            Modification for "non-square" character tables -> needed for space groups
! USAGE
!  GTIrep[characters,character table] gives the number of times Subscript[n, p] an irreduzible representation \[CapitalGamma]^p appears in a reducible representation with given characters. 
! INPUT
!  characters of red. rep., character table
! OUTPUT
!  list with number of occurence for each Irep
! GTPack OPTIONS
!  GOVerbose, GOFast
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTGroupQ
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

GTIrep[ch_,ct_,OptionsPattern[]] := Module[{g,np,ct1,ct2,h,tb,tb2,tb3,nrop,l, grp, run},
	If[Length[ct]==3,
		ct1=ct[[1]];
		h=Length[ch];
		grp = Flatten[ct1,1];
		g=Length[grp];
		
		If[OptionValue[GOFast],run=True,run=GTGroupQ[grp]];
		If[run,
			ct2=ct[[2]];
			np=Simplify[Table[Sum[Conjugate[ch[[i]]]*ct2[[j,i]]*Length[ct1[[i]]],{i,1,h}]/g,{j,1,Length[ct2]}]];
			(*--- print the output ---*)

    		If[OptionValue[GOVerbose],
    			tb=Flatten[Table[Which[	np[[i]]==0,	"a",
    									np[[i]]==1,	{ct[[3,i]],"\[CirclePlus]"},
    									np[[i]]>1,	{np[[i]] ct[[3,i]],"\[CirclePlus]"}],{i,1,Length[ct2]}]];
    			tb2=Delete[tb,Position[tb,"a"]];
    			nrop=Length[Position[tb2,"\[CirclePlus]"]];
				l=Length[tb2]-Length[Position[tb2,"\[CirclePlus]"]];
				tb3=If[l==nrop,Delete[tb2,Position[tb2,"\[CirclePlus]"][[nrop]]],None];
				Print[Row[tb3]];
    		,None];    		
    	Return[np],
    	Print[$ErrorNoGroup];
        Return[]],
    
    Print[$ErrorInvCt];
    Return[]
    ]
	]

(*
***)

(****n* /GTIrepDimension
! NAME
!  GTIrepDimension
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  09/07/2015 initial commit
! USAGE
!  gives the dimension of an irreducible representations.
! INPUT
!  group
! OUTPUT
!  integer
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTGroupQ
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

GTIrepDimension[irep_,OptionsPattern[]] := Module[{ll}, 
   If[OptionValue[GOFast],None,If[GTGroupQ[irep],None,
   	Print[$ErrorNoGroup];Abort[]]];
   ll=Length[irep[[1]]];
   Return[ll]
   ];

(*
***)


(****n* /GTIrepMatrixView
! NAME
!  GTIrepMatrixView
! AUTHOR
!  W. Hergert
! PACKAGE
!  RepresentationTheory.m
! MODIFICATION HISTORY
!   January 2016  : first version
!   April 2016 : revision and transfer from Test.m to Representationtheory.m
! USAGE
!  GTIrepMatrixView[group, character table, number] gives a neat print of the matrices of the 
!  irreducible representation number of group.
! INPUT
!  group            -  group as installed with GTInstallGroup
!  character table  -  the character table of the group
!  number           -  number of irreducible representation
! OPTIONS
!  -
! OUTPUT
!  prints a table of the matrices and returns all matrices of the Irep.
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTGetIrep
! GTPack NOTEBOOKS 
! GTIrepMatrixView.nb in Wolfram_Devel/5_Representation_Theory
! DESCRIPTION
!  its only a a version of GTGetIrep for better output
! LITERATURE
!  -
! TODO
!  
! RELEASE
!
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)
 
GTIrepMatrixView[grp_, ct_, repnr_] := Module[{ir,iru,lir,liru,collect,i,j,mats,prt,lc }, 
	ir = GTGetIrep[grp, repnr, ct]; iru = Union[ir];
	lir = Length[ir];liru = Length[iru];
    Print[lir, " group elements , ", liru, " diffrent rep. matrices"];
    collect = Table[{}, {liru}];
    Do[
       Do[
          If[ir[[i]] === iru[[j]],
             collect[[j]] = Append[collect[[j]], grp[[i]]],
             None
          ]
      , {j, 1, liru}]
   , {i, 1, lir}];
   mats = Map[MatrixForm[#] &, iru]; 
   lc=Length[collect];
   If[lc > 12, 
      prt = {Take[collect, {1, 12}], Take[mats, {1, 12}]};
      Print[Grid[prt, Frame -> All,Dividers -> {2 -> Black, {2 -> GTDividerColor1}}, 
            Background -> {None, {1 -> GTBackGroundColor1}}]];
      prt = {Take[collect, {13, lc}], Take[mats, {13, lc}]};
      Print[Grid[prt, Frame -> All,Dividers -> {2 -> Black, {2 -> GTDividerColor1}}, 
            Background -> {None, {1 -> GTBackGroundColor1}}]],
      prt = {collect, mats};
      Print[Grid[prt, Frame -> All, Dividers -> {2 -> Black, {2 -> GTDividerColor1}}, 
            Background -> {None, {1 -> GTBackGroundColor1}}]]
    ];Return[ir]
 ]
 
 (*
 ***)


(****n* /GTIrepNotation
! NAME
!  GTIrepNotation
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m
! MODIFICATION HISTORY
!   October 2014: first version
! USAGE
!   gives the notation of the irreducible representations for the character table
! INPUT
!   rotation matrix
! OUTPUT
!  index
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetRotationMatrix, GTWhichRepresentation, GTChangeRepresentation, GTGetSymbol
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  interprets the character table and reads the needed notation from the list GTcttable (Symbols.m).
! LITERATURE
!  -
! TODO
! -
! RELEASE
!
! PROBLEMS
! - 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTIrepNotation[ct0_, notation_] := 
 Module[{detct, ncl, classinfo, tablepos1, tablepos2, ngroups, 
   booltable, trueq, clin, classpos, permmatrix, cttab, 
   namepos, names, nameindex,ct, clrot,postable, ch0, cl1,ch1, pos1, pos2,i, srep, lp},
  detct = Det[ct0[[2]]];
  ncl = Length[ct0[[1]]];
  nameindex = Which[
  	notation == "Bethe", Return[Table[Superscript["\[CapitalGamma]", ToString[y]],{y,1,ncl}]],
    notation == "Mulliken", 3,
    notation == "Bouckaert", 4];
  (*--- Shorten the character table and deleting extra classes (double groups) and switch to O(3) representation ---*)  
  clrot = Table[Table[GTGetRotationMatrix[ct0[[1, i, j]]], {j, 1,Length[ct0[[1, i]]]}], {i, 1, ncl}];
  postable = Table[Position[clrot, clrot[[i, 1]]], {i, 1, ncl}];
  ch0 = Transpose[ct0[[2]]];
  cl1 = {};
  ch1 = {};
  Do[
  lp = Length[postable[[i]]];
  If[lp > 0,
   Which[
    lp == 1,
    pos1 = postable[[i, 1, 1]];
    cl1 = Append[cl1, clrot[[pos1]]];
    ch1 = Append[ch1, ch0[[pos1]]],
    lp == 2,
    pos1 = postable[[i, 1, 1]];
    pos2 = postable[[i, 2, 1]];
    postable[[pos2]] = {};
    cl1 = Append[cl1, clrot[[pos1]]];
    ch1 = Append[ch1, ch0[[pos1]]];
    ]
   ]
  , {i, 1, ncl}];
 srep = GTWhichRepresentation[GOVerbose -> False];
 GTChangeRepresentation["O(3)", GOVerbose -> False];
 {GTGetSymbol[cl1],Transpose[ch1]};
     
  ct={GTGetSymbol[cl1],Transpose[ch1]};(*{GTGetSymbol[ct0[[1]]],ct0[[2]]};*)
  classinfo = GTClassInfo[ct] // Simplify;
  (*--- Where to find the character table ---*)
  
  If[Length[Position[GTcttable, {ncl, detct}]] > 0,
  	tablepos1 = First@Flatten@Position[GTcttable, {ncl, detct}],
   If[Length[Position[GTcttable, {ncl, -detct}]] > 0, 
    tablepos1 = First@Flatten@Position[GTcttable, {ncl, -detct}], 
    Print["Error: Notation can not be found for the given group. The notation is switched to Bethe notation."];  
    Return[Table[Superscript["\[CapitalGamma]", ToString[y]],{y,1,ncl}]];]];
  (*--- How many groups are possible ---*)
  
  ngroups = Length[GTcttable[[tablepos1, 2]]];
  (*--- Which group is reasonable ---*)
  booltable = Table[
    clin = GTcttable[[tablepos1, 2, i, 1]];
    Length[Union[clin, classinfo]] == Length[classinfo], {i, 1, 
     ngroups}];
  trueq = Length[Position[booltable, True]];
  tablepos2 = Which[
    trueq == 0, 
    Print["Error: Notation can not be found for the given group. The notation is switched to Bethe notation."]; 
    Return[Table[Superscript["\[CapitalGamma]", ToString[y]],{y,1,ncl}]];
    ,
    trueq == 1, First[Flatten[Position[booltable, True]]],
    trueq == 0, 
    Print["Error: Multiple notations could be possible. The algorithm is stopped at this point and the notation is switched to Bethe notation."]; 
    Return[Table[Superscript["\[CapitalGamma]", ToString[y]],{y,1,ncl}]];
    ];
  (*--- Is the order of the classes correct? ---*)  
  clin = GTcttable[[tablepos1, 2, tablepos2, 1]];
  classpos = Table[First@Flatten@Position[clin, classinfo[[j]]], {j, 1, ncl}];
  permmatrix = Table[0, {i, 1, ncl}, {j, 1, ncl}];
  Do[permmatrix[[classpos[[i]], i]] = 1, {i, 1, ncl}];
  (*--- change the standard representation back to the input value ---*)
  GTChangeRepresentation[srep, GOVerbose -> False];
  (*--- allocate names ---*)
  cttab = GTcttable[[tablepos1, 2, tablepos2, 2]];
  names = GTcttable[[tablepos1, 2, tablepos2, nameindex]];
  Table[namepos = First@Flatten@Position[cttab, permmatrix.ct[[2, i]]];names[[namepos]], {i, 1, ncl}]
  ]
(*
***)

(****n* /GTMakeEigenvectors
! NAME
!  GTMakeEigenvectors
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  02/02/2014 initial commit
! USAGE
!  gives a set of linear independent eigenvectors from a given eigenvector. The Module is only needed within GTGetIrepGT1!
! INPUT
!  eigenvector, regular representation
! OUTPUT
!  list of representation matrices
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
!	See E. Blokker, A Theory for the Construction of the Irreducible Representations of Finite Groups .
!
! TODO
!  
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTMakeEigenvectors[ev_, reg_] := Module[{vv, proj,lg}, 
   lg = Length[reg];
   vv = Table[reg[[i]].ev, {i, 1, lg}];
   proj = Complement[Orthogonalize[vv], {Table[0, {i, 1, lg}]}]
   ];

(*
***)

(****n* /GTNumberOfIreps
! NAME
!  GTNumberOfIreps
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  09/07/2015 initial commit
! USAGE
!  gives the number of non-equivalent irreducible representations.
! INPUT
!  group
! OUTPUT
!  number of irreducible representations
! GTPack OPTIONS
!  GOFast, GOVerbose
! STANDARD OPTIONS
!  
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!
! LITERATURE
!	!
! TODO
!  
! RELEASE
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTNumberOfIreps[grp_,OptionsPattern[]] := Module[{ll}, 
   If[OptionValue[GOFast],None,If[GTGroupQ[grp],None,
   	Print[$ErrorNoGroup];Abort[]]];
   ll=Length[GTClasses[grp]];
   If[OptionValue[GOVerbose],
   Print["The given group has "<>ToString[ll]<>" non-equivalent irreducible representations."],None];
   Return[ll]
   ];

(*
***)

(****n* /GTPrintMullikenCharacterTable
! NAME
!  GTPrintMullikenCharacterTable
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  
! USAGE
!  Prints the character table with Mulliken symbols
! INPUT
!  classes, characters
! OUTPUT
!  text
! GTPack OPTIONS
!
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTGetSymbol, GTGetRotationMatrix, GTGetMatrix, GTOrderOfElement  
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

GTPrintMullikenCharacterTable[cl_,chi_] := 
 Module[{nir, ctmat, properpos, improperpos, proper, nproper, 
   nimproper, order, maxorder, primes, principalaxis, 
   principalaxispos, principalclass, pcl, axespositions, 
   axespositions2f, axestmp, planargroup, zpos, rotaxes, twofpos, 
   twofpostmp, twofrot, IEepos, invnames, reflectionpos, nreflections,
    horizontal, horizontaltmp, irnamelist, irnames, irnamestmp, 
   subind, subindgu, mn, mnpos, Bchar, indtmp, indtmp2, gerade, 
   indmax, ind, names, nocon, Conj, conind, conjnamelist, d,ct},
   ct={GTGetSymbol[cl],chi};
  (*--- check for Ci, Cs, C2 ---*)
  If[ct[[2]] == {{1, 1}, {1, -1}},
   (*--- C2 ---*)
   If[Det[GTGetRotationMatrix[ct[[1, 2, 1]]]] > 0,
    Return[{"A", "B"}],
    (*--- Ci,Cs ---*)    
    If[GTGetRotationMatrix[ct[[1, 2, 1]]] == GTGetRotationMatrix[IEe], 
     Return[{Subscript["A", "g"], Subscript["A", "u"]}], 
     Return[{Superscript["A", "'"], Superscript["A", "''"]}]]
    ]
   ];
  (*--- initial variables ---*)
  nir = Length[ct[[1]]];
  ctmat = Table[GTGetMatrix[ct[[1, i]]], {i, 1, nir}];
  (*--- Properties of classes ---*)
  properpos = {};
  improperpos = {};
  proper = 
   Table[d = Det[GTGetRotationMatrix[ct[[1, T, 1]]]]; 
    If[d == 1, properpos = Append[properpos, T], 
     improperpos = Append[improperpos, T]]; d, {T, 1, nir}];
  properpos = Delete[properpos, 1];
  nproper = Length[properpos];
  nimproper = Length[improperpos];
  order = Table[GTOrderOfElement[ct[[1, T, 1]]], {T, 1, nir}];
  (*--- Conjugate Representations ---*)
  nocon = True;
  Conj = Table[0, {T, 1, nir}];
  conind = 
   Table[First@First@Position[ct[[2]], Conjugate[ct[[2, T]]]], {T, 1, 
     nir}];
  Do[If[Conj[[T]] == 0 && conind[[T]] == T,
    nocon = False;
    Conj[[T]] = "",
    Conj[[T]] = "2";
    Conj[[conind[[T]]]] = "1";
    ], {T, 1, nir}];
  (*--- principal axis ---*)
  maxorder = Max[proper*order];
  principalclass = Flatten@Position[proper*order, maxorder];
  pcl = Length[principalclass];
  axespositions = Table[
    Table[
     First@
      Flatten@Position[
        Table[(ctmat[[principalclass[[i]], j]].axeslist[[k, 
              2]]).axeslist[[k, 2]], {k, 1, Length[axeslist]}], 1],
     {j, 1, Length[ct[[1, principalclass[[i]]]]]}]
    , {i, 1, pcl}];
  axestmp = 
   Union[Flatten[
     Table[axeslist[[axespositions[[i]], 2]], {i, 1, 
       Length[axespositions]}], 1]];
  zpos = Flatten[Position[axestmp, {0, 0, 1}]];
  principalaxis = If[Length[zpos] > 0, {0, 0, 1}, First@axestmp];
  planargroup = If[Length[axestmp] > 1 && maxorder > 2, False, True];
  principalaxispos = 
   principalclass[[First@Flatten@Position[axestmp, principalaxis]]];
  rotaxes = Length@Union[Flatten[Table[
       Table[
        First@
         Flatten@Position[
           Table[(ctmat[[properpos[[i]], j]].axeslist[[k, 
                 2]]).axeslist[[k, 2]], {k, 1, Length[axeslist]}], 
           1],
        {j, 1, Length[ct[[1, properpos[[i]]]]]}]
       , {i, 1, Length[properpos]}]]];
  (*--- 2fold-rotations perpendicular to the principal axis ---*)
  
  twofrot = Flatten[Position[proper*order, 2]];
  axespositions2f = 
   Table[First@
     Flatten@Position[
       Table[(GTGetMatrix[
            ct[[1, twofrot[[i]], 1]]].axeslist[[j, 2]]).axeslist[[j, 
           2]], {j, 1, Length[axeslist]}], 1], {i, 1, 
     Length[twofrot]}];
  axestmp = 
   Table[axeslist[[axespositions2f[[i]], 2]], {i, 1, 
     Length[axespositions2f]}];
  twofpostmp = 
   If[maxorder > 2, 
    Flatten@Position[
      Table[axestmp[[i]].principalaxis, {i, 1, Length[axestmp]}], 0], 
    axespositions2f];
  twofpos = 
   If[maxorder > 2, 
    Table[twofrot[[twofpostmp[[i]]]], {i, 1, Length[twofpostmp]}]];
  (*--- Inversion ---*)
  IEepos = Flatten[Position[ct[[1]], IEe]];
  invnames = If[Length[IEepos] > 0,
    Table[If[ct[[2, T, IEepos[[1]]]] < 0, "u", "g"], {T, 1, nir}]
    , Table["", {T, 1, nir}]];
  (*--- Horizontal reflections ---*)
  
  reflectionpos = 
   If[Length[IEepos] > 0, 
    Delete[improperpos, 
     First@First@Position[improperpos, IEepos[[1]]]], improperpos];
  nreflections = Length[reflectionpos];
  horizontaltmp = {};
  If[planargroup,
   Do[If[GTGetRotationMatrix[
        ct[[1, reflectionpos[[i]], 
          1]]].principalaxis.principalaxis == -1, 
     horizontaltmp = Append[horizontaltmp, reflectionpos[[i]]]], {i, 
     1, nreflections}]];
  horizontal = {};
  Do[If[Length[ct[[1, horizontaltmp[[i]]]]] == 1, 
    horizontal = Append[horizontal, horizontaltmp[[i]]]], {i, 1, 
    Length[horizontaltmp]}];
  (*--- Dimension of representations ---*)
  
  irnamelist = {"A", "E", "T", "F", "H", "I"};
  irnamestmp = Table[irnamelist[[ct[[2, T, 1]]]], {T, 1, nir}];
  Bchar = If[maxorder == 2 && Length[twofpos] > 1,
    Table[Min[
      Table[
       ct[[2, T, twofpos[[j]]]], {j, 1, Length[twofpos]}]],
     {T, 1, nir}],
    Table[ct[[2, T, principalaxispos]], {T, 1, nir}]
    ];
  irnames = 
   Table[If[
     irnamestmp[[T]] == "A" && Bchar[[T]] == -1 && planargroup && 
      rotaxes > 1, "B", irnamestmp[[T]]], {T, 1, nir}];
  conjnamelist = {"E", "F"};
  Do[If[Conj[[T]] == "", None, 
    irnames[[T]] = conjnamelist[[ct[[2, T, 1]]]]], {T, 1, nir}];
  (*--- Names ---*)
  
  mn = Table[Length[Position[irnames, irnames[[T]]]], {T, 1, nir}];
  subindgu = Table[If[mn[[T]] > 1, invnames[[T]], ""], {T, 1, nir}];
  (*--- Primes ---*)
  
  primes = If[Length[IEepos] == 0 && Length[horizontal] > 0 && nocon,
    Table[
     If[ct[[2, T, horizontal[[1]]]] > 0, "'", "''"]
     , {T, 1, nir}],
    Table["", {T, 1, nir}]
    ];
  (*--- Indices ---*)
  
  irnamestmp = 
   Table[irnames[[T]] <> subindgu[[T]] <> primes[[T]] <> 
     Conj[[T]], {T, 1, nir}];
  mnpos = 
   Union[Table[
     Flatten@Position[irnamestmp, irnamestmp[[T]]], {T, 1, nir}]];
  indtmp = Table[indmax = 1;
    If[Length[mnpos[[i]]] > 1,
     Table[
      gerade = True;
      Do[
       If[ct[[2, mnpos[[i, j]], twofpos[[k]]]] < 0, gerade = False]
       , {k, 1, Length[twofpos]}];
      If[gerade, 1, indmax = indmax + 1; indmax]
      ,
      {j, 1, Length[mnpos[[i]]]}],
     ""]
    , {i, 1, Length[mnpos]}];
  Do[
   k = 1;
   While[Length[indtmp[[i]]] > 0 && 
     Length[Position[indtmp[[i]], indtmp[[i, 1]]]] > 1,
    k = k + 1;
    Do[
     If[Re[ct[[2, k, mnpos[[i, j]]]]] < 0, indtmp[[i, j]]++]
     , {j, 1, Length[indtmp[[i]]]}]]
   , {i, 1, Length[indtmp]}];
  indtmp2 = Table[If[Length[indtmp[[i]]] > 0,
     While[Length[Position[indtmp[[i]], 1]] == 0, 
      indtmp[[i]] = indtmp[[i]] - 1];
     indtmp[[i]]
     ,
     ""], {i, 1, Length[indtmp]}];
  ind = Table["", {T, 1, nir}];
  Do[
   Do[
    If[Length[mnpos[[i]]] > 1,
      ind[[mnpos[[i, j]]]] = indtmp2[[i, j]], 
      ind[[mnpos[[i, j]]]] = ""];
    , {j, 1, Length[mnpos[[i]]]}]
   , {i, 1, Length[mnpos]}];
  subind = Table[ToString[ind[[T]]] <> subindgu[[T]], {T, 1, nir}];
  names = 
   Table[Subsuperscript[irnames[[T]], subind[[T]], 
     primes[[T]] <> Conj[[T]]], {T, 1, nir}];
  Return[names]]

(*
***) 

(****n* /GTPrintChars
! NAME
!  GTPrintChars
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m
! MODIFICATION HISTORY
!   04/04/2020: first version
! USAGE
!   Prints the character table of a group
! INPUT
!   group, classes, characters
! OUTPUT
!   none
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTGetIreps for the final output messages
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
GTPrintChars[grpk_, classes_, chars_,label_] := Module[{Pr, labels},
  Print["Classes:"];
  Table[Print[Subscript["C", i], " = ", classes[[i]]], {i, 1, 
    Length[classes]}];
  Print["Character Table:"];
  Pr = chars;
  labels = 
   Table[Length[classes[[i]]] classes[[i, 1]], {i, 1, 
     Length[classes]}];
  Pr = Prepend[Pr, labels];
  Pr = Table[ Prepend[Pr[[i]], If[i>1,label[[i-1]],""]], {i, 1, Length[Pr]}];
  Pr[[1, 1]] = "";
  Print[Grid[Pr, Frame -> All, ItemSize -> Full, 
    Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
    Background -> {{1 -> GTBackGroundColor1}, {1 -> 
        GTBackGroundColor1}, {1, 1} -> GTCornerColor}]];]
(*
***) 

(****n* /GTPrintIreps
! NAME
!  GTPrintIreps
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m
! MODIFICATION HISTORY
!   04/04/2018: first version
! USAGE
!   Prints the representation matrices of a group
! INPUT
!   group, ireps
! OUTPUT
!   none
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTGetIreps for the final output messages
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
GTPrintIreps[grpk_, ireps_,label_] := Module[{Pr}, 
  Print["Representation matrices:"];
  Pr = Table[Table[MatrixForm[ireps[[h, T]]], {T, 1, Length[grpk]}], {h, 1, Length[ireps]}];
  Pr = Prepend[Pr, grpk];
  Pr = Table[Prepend[Pr[[i]], If[i>1,label[[i-1]],""]], {i, 1,Length[Pr]}];
  Pr[[1, 1]] = "";
  Print[Grid[Pr, Frame -> All, ItemSize -> Full, 
  Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
    Background -> {{1 -> GTBackGroundColor1}, {1 -> 
        GTBackGroundColor1}, {1, 1} -> GTCornerColor}]];
  ]
(*
***) 

(****n* /GTProjectionOperator
! NAME
!  GTProjectionOperator
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  07. Juni: New Version of the Projection operator
! USAGE
!  GTProjectionOperator[group,ireducible representation,m,n,function,arguments]	gives the part of a given function with arguments which transforms like the m-th row and the n-th column of an irreducible representation.
! INPUT
!  group, irreducible representation matrices, indices m and n, function, arguments
! OUTPUT
!  function
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTGetMatrix, GTCharProjectionOperator
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
! 	The projection operator works for multiple cases, illustrated below:
!
!	O(3)
!	----
!	scalar function			+
!	2d function				+
!	3d function				+
!	higher dimensional		-
!
!	SU(2)
!	----
!	scalar function			-
!	2d function	(spinor)	+
!	3d function				-
!	higher dimensional		-
!
!	O(2)
!	----
!	scalar function			+
!	2d function				+
!	3d function				-
!	higher dimensional		-
!
!	permutation matrices
!	--------------------
!	scalar function			-
!	2d function				-
!	3d function				-
!	higher dimensional		-
!
! LITERATURE
!
! TODO
!  Needs to be tested!
! RELEASE
!
! PROBLEMS
! 
!--------------------------------------------------------------------------------
! SOURCE
*)


GTProjectionOperator[grp0_, ireps_, m_, n_, fu_,arg_,OptionsPattern[]] :=  Module[{grp,pseudoclasses,pseudochars,g,dimirep,run,tr},
	grp = GTGetMatrix[grp0];
	If[OptionValue[GOFast],run=True,run=GTGroupQ[grp]];
	If[run,
		g = Length[grp];
		dimirep = Length[ireps[[1]]];
		pseudoclasses = Table[{grp[[i]]},{i,1,g}];
		pseudochars = Table[ireps[[i,m,n]],{i,1,g}];
		tr = GTInternalCharProjectionOperator[pseudoclasses,pseudochars,fu,grp,dimirep,arg];
		
		,Print[$ErrorNoGroup];
		Abort[]];
	Return[tr]
]

(*
***)

(****n* /GTReality
! NAME
!  GTReality
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  
! USAGE
!  GTReality[classes,characters] estimates if a representation is potentially real, essentially complex or pseudo-real from the theorem of Frobenius and Schur.
! INPUT
!  classes,characters
! OUTPUT
!  information
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!  
! ERROR MESSAGES
!  
! GTPack MODULES
!  GTGetMatrix, GTGroupQ, GTSimplify
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

GTReality[classes_, chars_,OptionsPattern[]] := Module[{h, grpm, grp2, g, su, sum},
  If[OptionValue[GOFast]==False,If[GTGroupQ[Flatten[classes,1]],None,Print[$ErrorNoGroup];Abort[]]];
  (*--- How many classes? ---*)
  h = Length[classes];
  (*--- Use the matrices of the elements ---*)
  
  grpm = GTGetMatrix[#] & /@ classes// GTSimplify;
  (*--- Calculate T*T for each element T ---*)
  
  grp2 = Flatten[
    Table[MatrixPower[#, 2] & /@ grpm[[i]], {i, 1, h}] // GTSimplify, 1];
  (*--- Calculate the group order ---*)
  g = Length[grp2];
  (*--- Estimate the character of the element T*T from the character system ---*)
  
  su = chars[[#]] & /@ 
    Map[First, 
     Flatten[Table[Position[grpm, grp2[[i]]], {i, 1, g}], 1]];
  (*--- Calculate the sum according to the therem of Frobenius-Schur ---*)
  
  sum = Sum[su[[i]], {i, 1, g}]/g//FullSimplify;
  Which[sum == 1, "potentially real", sum == 0, "essentially complex",
    sum == -1, "pseudo-real"]
  ]

(*
***)


(****n* /GTRegularRepresentation
! NAME
!  GTRegularRepresentation
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  
! USAGE
!  GTRegularRepresentation[group] gives the regular representation of a group.
! INPUT
!  group
! OUTPUT
!  regular representation matrices
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
!--------------------------------------------------------------------------------
! SOURCE
*)

GTRegularRepresentation[grp_,OptionsPattern[]] := Module[{grpt,gt,el0,el,g,pos},
	If[OptionValue[GOFast],None,If[GTGroupQ[grp],None,Print[$ErrorNoGroup];Abort[]]];
	grpt = Expand@GTGetMatrix[grp];
	g = Length[grpt];
	gt = Expand@Table[grpt[[i]].Inverse[grpt[[j]]], {i, 1, g}, {j, 1, g}];
	el0 = Table[0, {i, 1, g}, {j, 1, g}];
	Table[
 		pos = Position[gt, grpt[[t]]];
 		el = el0;
 		Do[el[[pos[[i, 1]], pos[[i, 2]]]] = 1, {i, 1, g}];
 		el, {t, 1, g}]
]

(*
***)

 
(****n* /GTReorderCharacterTable
! NAME
!  GTReorderCharacterTable
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m
! MODIFICATION HISTORY
!   January 2015  : first version
!   April 2016    : from Test.m to RepresentationTheory.m
!   09.04.2017    : Reordering of classes and Ireps possible
! USAGE
!  GTReorderCharacterTable reorders a character table with respect to the classes and/or Ireps
! INPUT
!  grpname - name of the group (String)
!  ctab    - character table
! GTPack OPTIONS
!  GOIrepNotation - alternative notation of the Ireps
! STANDARD OPTIONS
!
! OUTPUT
!  reordered character table
! ERROR MESSAGES
!  if reordering vector has the wrong length
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  Wolfram_Devel/5_Represeentation_Theory/GTReorderCharacterTable.nb
! DESCRIPTION
!  Sometimes one has to compare character tables. The claases in the character table to compare with
!  can be in a different order than the standard GTPack output. The command allows to reorder the classes.
! LITERATURE
!  -
! TODO
!  
! RELEASE
!
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTReorderCharacterTable[grpname_, ctab_, col_, OptionsPattern[]] := 
  Module[{classes, characters, names, lc, col1, go, goo, nm, n, c1, c2,i, chars, j, spw, tab,
  	      col2,go1,names1,classes1},
   go = OptionValue[GOIrepNotation];
  (*---originbal character table---*)
   classes    = ctab[[1]];
   lc         = Length[classes];
   characters = ctab[[2]];
   names      = ctab[[3]];
  (*---rule interpretation---*)
   col1 = col[[1]];
   col2 = col[[2]];
   If[col1 == {},
      col1 = Table[i, {i, 1, lc}],
      If[Length[col1] == lc,
         None,
         Print["Error : wrong length of class reordering vector"]; Abort[]
      ]
   ];
   If[col2 == {},
      col2 = Table[i, {i, 1, lc}],
      If[Length[col2] == lc,
         None,
         Print["Error : wrong length of Irep reordering vector"]; Abort[]
      ]
   ];
  (*---interprete option---*)
   goo = False;
   If[Head[go] === List,
      nm  = {names, go};
      n   = {{grpname}, {" "}};
      goo = True,
      nm  = {names};
      n   = {{grpname}}
   ];
  (*---create reordered headline of the table---*)
   c1 = Map[Length[#] &, classes];
   c2 = Map[#[[1]] &, classes];
   Do[
   	  n = Append[n, {c1[[col1[[i]]]] c2[[col1[[i]]]]}]
   , {i, 1, lc}];
  (*---setup new character Table---*)
  chars    = Table[0, {lc}, {lc}];
  go1      = names1 = Table[0, {lc}];
  classes1 = Table[0, {lc}];
  Do[
     names1[[j]]   = names[[col2[[j]]]];
     classes1[[j]] = classes[[col1[[j]]]];
     If[goo,
        go1[[j]] = go[[col2[[j]]]],
        None
     ]
  , {j, 1, lc}];
  Do[
     Do[
        chars[[j, i]] = characters[[col2[[j]], col1[[i]]]]
     , {i, 1, lc}]
  , {j, 1, lc}];
  chars = Transpose[chars];
  If[goo,
     nm = {names1, go1},
     nm = {names1}
   ];
   spw = Transpose[Join[nm, chars]];
   spw = Join[n // Transpose, spw];
   If[goo, 
      tab = Grid[spw, Dividers -> {{2 -> GTDividerColor1, 3 -> Blue}, {2 -> GTDividerColor1}}, Frame -> All, 
                      Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {{1, 1} -> GTCornerColor, 
                      {{2, lc + 1}, {2, 2}} -> GTCornerColor}}, Alignment -> Right], 
      tab = Grid[spw, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, Frame -> All, 
      	              Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} -> GTCornerColor}, 
                      Alignment -> Right]
   ];
   Print[tab]; 
   Return[{classes1, Transpose[chars], names1}]
]
   
(*
***)


(****n* /GTSpinCharacters
! NAME
!  GTSpinCharacters
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  1/2/17: GTChangeRepresentation does not work within the command. The standard representation is changed via directly specifying grpdgrp.
! USAGE
!  GTSpinCharacters[classes] gives the character of the spinor representation for each class.
! INPUT
!  classes
! OUTPUT
!  List of characters of the classes
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTGroupQ, GTGetMatrix
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

GTSpinCharacters[classes_,OptionsPattern[]] := Module[{wr, chartab, i},
  If[OptionValue[GOFast],
  	None,
  	If[GTGroupQ[Flatten[classes, 1]],None,Print[$ErrorNoGroup];Abort[]];
  ];  
  
  wr = grpdgrp;
  Unprotect[grpdgrp];
  grpdgrp=2;
  chartab = Table[Tr@GTGetMatrix@classes[[i, 1]], {i, 1, Length[classes]}];
  grpdgrp=wr;
  Protect[grpdgrp];
  Return[chartab]
  ]

(*
***)

(****n* /GTSOCSplitting
! NAME
!  GTSOCSplitting
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  
! USAGE
!  GTSOCSplitting[character table] calculates the splitting of states due to spin-orbit coupling.
! INPUT
!  classes
! OUTPUT
!  List of characters of the classes
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTGroupQ, GTSpinCharacters, GTIrep
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

GTSOCSplitting[ct_,OptionsPattern[]] := Module[{cl, deepos, ordchars, extrachars, ordnames, extranames,ireps, names, np, tb, tb2, nrop, l, tb3, i, j, spinchars,grp},
  cl = ct[[1]];
  grp=Flatten[cl,1];
  If[OptionValue[GOFast],
   	None,
   	
    If[Not[grpdgrp==2||grpdgrp==5],Print[$ErrorGTSOCSplitting];Abort[];];
   	If[Not[GTGroupQ[grp]],Print[$ErrorNoGroup];Abort[]];
  ];
  deepos = First@First@Position[cl, DEe];
  ordchars = {};
  extrachars = {};
  ordnames = {};
  extranames = {};
  Do[If[ct[[2, i, deepos]] > 0, 
    ordchars = Append[ordchars, ct[[2, i]]]; 
    ordnames = Append[ordnames, ct[[3, i]]], 
    extrachars = Append[extrachars, ct[[2, i]]]; 
    extranames = Append[extranames, ct[[3, i]]]], {i, 1, Length[cl]}];
  spinchars = GTSpinCharacters[cl,GOFast->True];
  ireps = Table[GTIrep[ordchars[[i]]*spinchars, ct, GOFast -> True, GOVerbose -> False], {i, 1, Length[ordchars]}];
  names = Table[
    np = ireps[[j]];
    tb = Flatten[
      Table[Which[np[[i]] == 0, "a", 
        np[[i]] == 1, {ct[[3, i]], "\[CirclePlus]"}, 
        np[[i]] > 1, {np[[i]] ct[[3, i]], "\[CirclePlus]"}], {i, 1, 
        Length[cl]}]];
    tb2 = Delete[tb, Position[tb, "a"]];
    nrop = Length[Position[tb2, "\[CirclePlus]"]];
    l = Length[tb2] - Length[Position[tb2, "\[CirclePlus]"]];
    tb3 = If[l == nrop, Delete[tb2, Position[tb2, "\[CirclePlus]"][[nrop]]], None], {j,1, Length[ireps]}];
  Print[Grid[Prepend[Table[Append[Prepend[ordchars[[i]], ordnames[[i]]],Row[names[[i]]]], {i, 1, Length[ordchars]}],Flatten@{"ord. irred. representations", cl[[;; , 1]], "Spin-orbit coupling"}], Frame -> All,Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
    Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} -> GTCornerColor}]];
  Print[Grid[Prepend[Table[Prepend[extrachars[[i]], extranames[[i]]], {i, 1, Length[extrachars]}],Flatten@{"extra-representations", cl[[;; , 1]]}], Frame -> All,Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}},Background -> {{1 -> GTBackGroundColor1}, {1 -> 
        GTBackGroundColor1}, {1, 1} -> GTCornerColor}]];]

(*
***)

(****n* /GTWignerProjectionOperator
! NAME
!  GTWignerProjection
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  December 2014: remove arguments theta and phi from Y[l,m]
!  April 2021: replaced by a new version using GTAngularMomentumRep
! USAGE
!  GTWignerProjectionOperator[group,irep,l,m,i,j] applies the projection operator on spherical harmonics. 
!  See Tsukerblat (2006) P. 165 and P.168
! INPUT
!  group, IREP matrices, l, m, i, j
! OUTPUT
!  sum of spherical harmonics
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTGetEulerAngles
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
(*
GTWignerProjectionOperator[grp_, ir_, l_, m_, i_, j_,OptionsPattern[]] := Module[{wig, m2, angles, erg}, 
  If[OptionValue[GOFast],None,If[GTGroupQ[grp],None,Print[$ErrorNoGroup];Abort[]]];
   Clear[\[Theta], \[Phi]];
   angles = GTGetEulerAngles[grp];
   wig = Table[If[Abs[ir[[k, i, j]]]>0,Sum[WignerD[{l, m2, m}, angles[[k,1, 1]], angles[[k, 1, 2]],angles[[k, 1, 3]]] (angles[[k, 2]]^l) Y[l,m2], {m2, -l, l}],0], {k, 1, Length[grp]}];
   erg = (Length[ir[[1, 1]]]/Length[grp])*Sum[Conjugate[ir[[k, i, j]]] wig[[k]], {k, 1, Length[grp]}]//Simplify;
   Return[Expand[erg]]
  ]
*)
GTWignerProjectionOperator[group_, irep_, l_, m_, i_, j_, OptionsPattern[]] :=
 Module[{ar, proj},
  ar = GTAngularMomentumRep[group, l, GOHarmonics -> OptionValue[GOHarmonics]];
  proj = Total[ar[[ ;; ]]*Conjugate[irep[[;; , i, j]]]] Length[irep[[1]]]/
     Length[irep];
  Total[proj.Table[Y[l, mm], {mm, -l, l, 1}]]
  ]

(*
***)


(****n* /GTSymmetrizeProductChars
! NAME
!  GTSymmetrizeProductChars
! AUTHOR
!  W. Hergert
! PACKAGE
!  RepresentationTheory.m
! MODIFICATION HISTORY
!  * 15.07.2016 : first version
!  * 10.10.2018 : implementation in package, check of headers
! USAGE
!  GTSymmetrizedProductChars[ct,irep,power,symm] gives the characters of the nth power of an irreducible 
!  represantation of a given group
! INPUT
!  o ct    - character table of the group under consideration
!  o irep  - number of the irreducibel representation  
!  o power - the power of the Irep to be investigated
!  o sym   - decide if symmetri or antisymmetric power is considered, can be only "sym" od "asymm" 
! OUTPUT
!  table
! GTPack OPTIONS
! o GOVerbose 
!
!     - False : no additional information (standard)
!     - True  : additional information about character table
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTMagnetic, GTInstallGroup, GTGroupQ, GTInvSubGroupQ
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  - 
! TODO
! -
! RELEASE
!  1.0.1
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)


 GTSymmetrizedProductChars[ct_, irep_, power_, sign_,  OptionsPattern[]] := Module[
 	{classes, chars, names, pch,  part, ns, head, pchar, g, gg, qs, rs, rsu, pos, prod, i, j, k, l, fac, su, ncl, tab, tn, no, 
     sort,verb,gpw},
    {classes, chars, names} = ct;
     ncl  = Length[classes];
     verb = OptionValue[GOVerbose];
  (*--- test of input ---*)
    If[sign == "symm" || sign == "asymm",
       None,
       Print["Error: ", sign, " not valid, has to be symm or asymm"]; Abort[]
    ];
 (*--- select the names ---*)
    pch = names[[irep]];
    If[sign == "symm",
       sort = "[T]"^power,
       sort = "{T}"^power
    ];
 (*--- calculate the partitions ---*)  
    part = Tally[#] & /@ IntegerPartitions[power];
    ns   = Length[part];
 (*--- calculate the characters --*)
    head = {}; pchar = {};
 (*--- all classes ---*)    
    Do[
 (*--- calculate the powers of elements ---*)    
       g = classes[[i, 1]]; head = Append[head, g];
       gpw = {g}; gg = g;
       Do[
          gg = g\[SmallCircle]gg; gpw = Append[gpw, gg]
       , {l, 2, power}];
 (*--- sum over all partitions ---*)
       su = 0;
       Do[
 (*--- contribution of one partition ---*)
          prod = 1; 
          rsu  = 0;
          Do[
             qs   = part[[k, j, 1]];
             rs   = part[[k, j, 2]]; 
             rsu  = rsu + rs; 
             gg   = gpw[[qs]];
             pos  = Flatten[Position[classes, gg]][[1]];
             prod = prod* chars[[irep, pos]]^rs/qs^rs/rs!;
          , {j, 1, Length[part[[k]]]}];
          If[sign == "asymm",
             fac = (-1)^(rsu - power),
             fac = 1
          ];
          su = su + prod*fac
       , {k, 1, ns}];
       pchar = Append[pchar, su]
    , {i, 1, ncl}];
 (*--- construct table ---*)
    If[verb, 
       tab = Append[Prepend[chars, head], pchar] // Transpose;
       tn  = {sort, names, pch} // Flatten; no = Length[tn];
       tab = Prepend[tab, tn] // Transpose;
       Print[Grid[tab, Frame -> All, Alignment->Right,Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1, no -> GTDividerColor1}}, 
                       Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1, no -> GTCornerColor}, 
                       {{1, 1} -> GTCornerColor, {no, 1} -> GTCornerColor}}
                 ]
           ],
      None
    ]; 
    Return[pchar]
  ]

(*
***)


(****n* /GTVectorRep
! NAME
!  GTVectorRep
! AUTHOR
!  W. Hergert
! PACKAGE
!   RepresentationTheory.m
! MODIFICATION HISTORY
!  * 07.08.2016 : first version
!  * 10.10.2018 : implementation in package, check of header
! USAGE
!  GTVectorReps[character table] calculates the vector representation to a group with character table.
! INPUT
!   character table
! OUTPUT
!  vector representation
! GTPack OPTIONS
!  o GOVerbose 
!
!     - False : no additional information (standard)
!     - True  : additional information about character table
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  o GTGetMatrix
!  o GTCharacters
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  see book
! TODO
!  -
! RELEASE
!  1.0.1
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTVectorRep[ct_, OptionsPattern[]] := Module[{selem, irs, verb,vrep},
  verb  = OptionValue[GOVerbose];
  selem = GTGetMatrix[#] & /@ Flatten[ct[[1]]];
  irs   = GTCharacters[selem, GOClasses -> True];
  vrep  = GTIrep[irs, ct, GOVerbose -> verb];
  Return[vrep]
]


(*
***)


(****n* /GTCharacters
! NAME
!  GTCharacters
! AUTHOR
!  W. Hergert
! PACKAGE
!   RepresentationTheory.m
! MODIFICATION HISTORY
!  * 01.08.2016 : first version
!  * 20.02.2019 : command dissapeared somehow. Put here again. Probably not transferred from test.m
! USAGE
!  GTCharacters[matrices] calculates the characters of a representation matrix or a list of representation matrices.
! INPUT
!  matrix or list of matrices
! OUTPUT
!  charcters
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTGetMatrix, GTClasses
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.1
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTCharacters[elements_, OptionsPattern[]] := Module[{cls,elm1,cc,lc,elm2,i,chars},
   cls = OptionValue[GOClasses];
   elm1=elements;
   If[Length[Dimensions[elements]]>2,
      If[cls,
         cc = GTClasses[elm1]; lc = Length[cc];
         elm2 = Table[0, {lc}];
         Do[
            elm2[[i]] = cc[[i, 1]]
         , {i, 1, lc}],
         elm2=elm1
      ];
      chars = Map[Tr[#] &, elm2],
      chars = Tr[GTGetMatrix[elements]]
   ];
   Return[chars]
]

(*
***)

(*-------------------------- Attributes ------------------------------*)

Attributes[GTCharacterTable]={Protected, ReadProtected}
Attributes[GTCharProjectionOperator]={Protected, ReadProtected}
Attributes[GTProjectionOperator]={Protected, ReadProtected}
Attributes[GTClebschGordanSum]={Protected, ReadProtected}
Attributes[GTGetIreps]={Protected, ReadProtected}
Attributes[GTGetIrep]={Protected, ReadProtected}
Attributes[GTIrep]={Protected, ReadProtected}
Attributes[GTClebschGordanCoefficients]={Protected, ReadProtected}
Attributes[GTDirectProductChars]={Protected, ReadProtected}
Attributes[GTDirectProductRep]={Protected, ReadProtected}
Attributes[GTReality]={Protected, ReadProtected}
Attributes[GTWignerProjectionOperator]={Protected, ReadProtected}

End[]


EndPackage[]
