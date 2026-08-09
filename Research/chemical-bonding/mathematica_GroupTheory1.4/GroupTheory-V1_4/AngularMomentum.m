(* ::Package:: *)

(****m* /AngularMomentum.m
!
! NAME
!  AnularMomentum.m
! AUTHOR
!  M. Geilhufe
! MODIFICATION HISTORY
!  06/23/20 : initially created and documented  
! USAGE
!  all modules connected to angular momentum
!
! ERROR MESSAGES
!  
! GTPack MODULES
!
!
! GTPack NOTEBOOKS 
!  none in the moment
!
! DESCRIPTION
!  AngularMomentum.m contains all modules connected to angular momentum.
! LITERATURE
! 
! TODO
!  
!
! PROBLEMS
!
!--------------------------------------------------------------------------------
!
***)
BeginPackage["GroupTheory`AngularMomentum`",{"GroupTheory`Basic`","GroupTheory`Symbols`","GroupTheory`Auxiliary`"}]

 GTAngularMomentumChars			::usage = "GTAngularMomentumChars[\*StyleBox[\"classes,j\",\"TI\"]] calculates the character system using (2\*StyleBox[\"j\",\"TI\"]+1)-dimensional representation matrices according to the angular momentum operator."
 GTAngularMomentumRep			::usage = "GTAngularMomentumRep[\*StyleBox[\"group,j\",\"TI\"]] gives a (2\*StyleBox[\"j\",\"TI\"]+1)-dimensional matrix representation using the \*ButtonBox[StyleBox[\"WignerD\", \"SR\"], BaseStyle->\"Link\", ButtonData->\"paclet:ref/WignerD\"]-function."
 GTJz  	          ::usage = "GTJz[\*StyleBox[\"J\",\"TI\"]] gives the z component of the total angular momentum operator for angular momentum \*StyleBox[\"J\",\"TI\"]."
 GTJx  	          ::usage = "GTJz[\*StyleBox[\"J\",\"TI\"]] gives the x component of the total angular momentum operator for angular momentum \*StyleBox[\"J\",\"TI\"]."
 GTJy  	          ::usage = "GTJz[\*StyleBox[\"J\",\"TI\"]] gives the y component of the total angular momentum operator for angular momentum \*StyleBox[\"J\",\"TI\"]."
 GTJplus  	      ::usage = "GTJz[\*StyleBox[\"J\",\"TI\"]] gives the raising operator for angular momentum \*StyleBox[\"J\",\"TI\"]."
 GTJminus  	      ::usage = "GTJz[\*StyleBox[\"J\",\"TI\"]] gives the lowering operator for angular momentum \*StyleBox[\"J\",\"TI\"]."
 GTJMatrix        ::usage = "GTJTransform[\*StyleBox[\"J, element\",\"TI\"]] gives the representation matrix of an \*StyleBox[\"element\",\"TI\"] for an irreducible representation \*StyleBox[\"J\",\"TI\"] of O(3)."
 GTJTransform     ::usage = "GTJTransform[\*StyleBox[\"J, m, element\",\"TI\"]] applies a symmetry transformation to the basis functions of an irreducible representation of O(3) with angular momentum \*StyleBox[\"J\",\"TI\"]."
 GTPauliMatrix    ::usage = "GTJz[\*StyleBox[\"x, J\",\"TI\"]] gives the generalized Pauli matrix for the component \*StyleBox[\"x\",\"TI\"] and angular momentum \*StyleBox[\"J\",\"TI\"]."
 
 Options[GTAngularMomentumChars]      = {GOFast->GOFastValue,GOEvenOdd->-1000}
 Options[GTAngularMomentumRep]        = {GOFast->GOFastValue,GOHarmonics->"Complex",GOEvenOdd->-1000}

Begin["`Private`"] (* Begin Private Context *) 

(****n* /GTAngularMomentumChars
! NAME
!  GTAngularMomentumChars
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  15/05/2013 Change of Structure, Error Messages
! USAGE 
!  GTAngularMomentumChars[classes,j] calculates the character system using (2j+1)-dimensional representation matrices according to the angular momentum operator.
! INPUT
!  classes, j
! OUTPUT
!  list of characters
! GTPack OPTIONS
!  GOFast
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGroupQ, GTGetEulerAngles
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

GTAngularMomentumChars[classes_,J_,OptionsPattern[]] := Module[{grp,run},
    If[OptionValue[GOFast],	run=True,
    						
    						grp = Flatten[classes,1];
    						run=GTGroupQ[grp]];
	If[run,
		 Return[Map[Tr, GTAngularMomentumRep[classes[[;; , 1]], J, GOFast -> True,GOEvenOdd->OptionValue[GOEvenOdd]]]],
    Print[$ErrorNoGroup];Return[]]]
 

(*
***)

(****n* /GTAngularMomentumRep
! NAME
!  GTAngularMomentumRep
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  RepresentationTheory.m 
! MODIFICATION HISTORY
!  15/05/2013 Change of Structure, Error Messages
!  20/01/2020 Tranformation coefficients for tesseral harmonics are revised
!  07/05/2021 Bug fix with real representation
! USAGE
!  GTAngularMomentumRep[group,j] gives a (2j+1)-dimensional matrix representation using the WignerD-function.
! INPUT
!  group, j
! OUTPUT
!  matrix representation
! GTPack OPTIONS
!  GOFast, GOHarmonics
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTGetEulerAngles, GTGroupQ
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
GTalm[l_, m_, ms_] := Which[
  m > 0,
  (KroneckerDelta[m, -ms] + (-1)^m*KroneckerDelta[m, ms])/Sqrt[2],
  m == 0,
  KroneckerDelta[m, ms],
  m < 0,
  (KroneckerDelta[m, ms] - (-1)^m*KroneckerDelta[m, -ms]) I/Sqrt[2]]
  
(*GTalm[l_, m_, ms_] := Which[m > 0, (-1)^m*(KroneckerDelta[m, ms] + KroneckerDelta[-m, ms])/
   Sqrt[2], m == 0, KroneckerDelta[m, ms], 
  m < 0, (-1)^m*(KroneckerDelta[m, -ms] - KroneckerDelta[m, ms])/(I*
     Sqrt[2])]*)

GTAngularMomentumRep[grp_, l_, OptionsPattern[]] := Module[{eang, g,run,amr, trafo, invtrafo,amrtmp,exp},
	If[OptionValue[GOFast],	run=True,
    						
    						run=GTGroupQ[grp]];
	If[run,
			g = Length[grp];
 			If[	grpdgrp == 4, (*--- permutation matrices ---*) 
 				Print["Error: The standard representation is given by permutation matrices. Matrices for an angular momentum representation can not be evaluated from permutation matrices!"];
 				Return[]
 				];
 			eang = GTGetEulerAngles[grp];		
 			If[IntegerQ[l],	
 				(*If l is an integer multiply with (-1)^l for improper rotations*)
 				exp=OptionValue[GOEvenOdd];
 				If[exp==-1000,exp=l];
 				amr=Simplify[Table[(If[Element[eang[[i,2]],Reals],eang[[i,2]],eang[[i,2]]/I]^exp)*GTJMatrix[l,grp[[i]]],{i,1,g}]],
 				amr=Simplify[Table[GTJMatrix[l,grp[[i]]],{i,1,g}]]];
 			(*--- In the case of real spherical harmonics ---*)
 			If[OptionValue[GOHarmonics]=="Real",
 				trafo = Table[GTalm[l,ms,m],{m,l,-l,-1},{ms,l,-l,-1}];
 				(*trafo = Table[GTalm[l,ms,m],{m,-l,l},{ms,-l,l}];*)
 				
 				invtrafo = Inverse[trafo];
 				amrtmp = amr;
 				amr = Table[invtrafo . amrtmp[[T]] . trafo,{T,1,g}]// GTSimplify
 				,None];
 			Return[amr],
 			
 			Print[$ErrorNoGroup];
			Return[]]]
 					
(*
***)

(****p* /GTJz
! NAME
!  GTJz
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  AngularMomentum.m
! MODIFICATION HISTORY
!   23/06/2020: first version
! USAGE
!   GTJz gives the z component of the total angular momentum operator for angular momentum J."
! INPUT
!   angular momentum J
! OUTPUT
!   corresponding matrix 
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
!  -
! TODO
! -
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTJz[J_] := DiagonalMatrix[Table[m, {m, J, -J, -1}]]

  (*
***)

(****p* /GTJplus
! NAME
!  GTJplus
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  AngularMomentum.m
! MODIFICATION HISTORY
!   23/06/2020: first version
! USAGE
!   GTJz gives the raising operator for angular momentum J."
! INPUT
!   angular momentum J
! OUTPUT
!   corresponding matrix 
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
!  -
! TODO
! -
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTJplus[J_] := DiagonalMatrix[Table[Sqrt[J (J + 1) - m (m + 1)], {m, J - 1, -J, -1}], 1]

  (*
***)

(****p* /GTJminus
! NAME
!  GTJminus
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  AngularMomentum.m
! MODIFICATION HISTORY
!   23/06/2020: first version
! USAGE
!   GTJz gives the lowering operator for angular momentum J."
! INPUT
!   angular momentum J
! OUTPUT
!   corresponding matrix 
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
!  -
! TODO
! -
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)


GTJminus[J_] := DiagonalMatrix[Table[Sqrt[J (J + 1) - m (m - 1)], {m, J, -J + 1, -1}], -1]

  (*
***)

(****p* /GTJx
! NAME
!  GTJx
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  AngularMomentum.m
! MODIFICATION HISTORY
!   23/06/2020: first version
! USAGE
!   GTJz gives the x component of the total angular momentum operator for angular momentum J."
! INPUT
!   angular momentum J
! OUTPUT
!   corresponding matrix 
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
!  -
! TODO
! -
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)


GTJx[J_] := (GTJplus[J] + GTJminus[J])/2

  (*
***)

(****p* /GTJy
! NAME
!  GTJy
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  AngularMomentum.m
! MODIFICATION HISTORY
!   23/06/2020: first version
! USAGE
!   GTJz gives the y component of the total angular momentum operator for angular momentum J."
! INPUT
!   angular momentum J
! OUTPUT
!   corresponding matrix 
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
!  -
! TODO
! -
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTJy[J_] := (GTJplus[J] - GTJminus[J])/(2 I)

  (*
***)

(****p* /GTJTransform
! NAME
!  GTJTransform
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  AngularMomentum.m
! MODIFICATION HISTORY
!   23/06/2020: first version
! USAGE
!   GTJTransform applies a symmetry transformation to the basis functions of an irreducible representation of O(3) with angular momentum J."
! INPUT
!   angular momentum J, m, element
! OUTPUT
!   sum of basis functions 
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
!  -
! TODO
! -
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTJTransform[J_, m_, el_] := Module[{su2mat, a, b,ms,k},
  su2mat = GTGetSU2Matrix[el];
  a = su2mat[[1, 1]];
  b = su2mat[[2, 1]];
  Which[
  	Simplify[Abs[b]] == 0,
        a^(J + m) Conjugate[a]^(J - m) Y[J, m],
    Simplify[Abs[a]] == 0,
        (-Conjugate[b])^(J + m) b^(J - m) Y[J, -m],
    Simplify[Abs[a]] > 0 && Simplify[Abs[b]] > 0,
        Sum[Sum[Sqrt[(J - m)! (J + m)! (J - ms)! (J + ms)!] a^(J + m - k) (-Conjugate[b])^k b^(ms - m + k) Conjugate[a]^(J - ms - k) Y[J,ms]/((J + m - k)! k! (ms - m + k)! (J - ms - k)!), {ms, J - k, m - k, -1}], {k, 0, J + m}]]
 ]

  (*
***)

(****p* /GTJMatrix
! NAME
!  GTJMatrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  AngularMomentum.m
! MODIFICATION HISTORY
!   23/06/2020: first version
! USAGE
!   GTJMatrix gives the representation matrix of an element for an irreducible representation J of O(3)."
! INPUT
!   angular momentum J
! OUTPUT
!   corresponding matrix 
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
!  -
! TODO
! -
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)


GTJMatrix[J_, el_] := Module[{m1,ms},Table[GTJTransform[J, m1, el] /. Y[Jj_, m_] :> Table[KroneckerDelta[m, ms], {ms, Jj, -Jj, -1}], {m1, J, -J, -1}]]

  (*
***)

(****p* /GTPauliMatrix
! NAME
!  GTPauliMatrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  AngularMomentum.m
! MODIFICATION HISTORY
!   23/06/2020: first version
! USAGE
!   GTPauliMatrix gives the generalized Pauli matrix for the component x and angular momentum J."
! INPUT
!   angular momentum J
! OUTPUT
!   corresponding matrix 
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
!  -
! TODO
! -
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTPauliMatrix[x_, J_] := {IdentityMatrix[2 J + 1], 2 GTJx[J], 2 GTJy[J], 2 GTJz[J]}[[x + 1]]

  (*
***)

Attributes[GTAngularMomentumChars]={Protected, ReadProtected}
Attributes[GTAngularMomentumRep]={Protected, ReadProtected}
Attributes[GTJz]={Protected, ReadProtected}
Attributes[GTJplus]={Protected, ReadProtected}
Attributes[GTJminus]={Protected, ReadProtected}
Attributes[GTJx]={Protected, ReadProtected}
Attributes[GTJy]={Protected, ReadProtected}
Attributes[GTPauliMatrix]={Protected, ReadProtected}

End[]

EndPackage[]
