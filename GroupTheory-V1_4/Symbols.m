(*
!
! 30.5. 21                            GTDivideColor2 corrected
!
*)


(****m* /Symbols.m
!
! NAME
!  Symbols.m
! AUTHOR
!  W. Hergert, M. Geilhufe
! MODIFICATION HISTORY
!   30.5.2021 :GTDividerColor2 corrected

! USAGE
!  Contains the lists of the group elements.
! 
! ERROR MESSAGES
!  
! GTPack MODULES
!	elm{} - list of elements
!	elm2d{} - list of 2d elements
!   elmown{} - empty list that can be filled by GTTableToElements 
!	elmpalette{} - Names of group elements for the palette
!
! GTPack NOTEBOOKS 
!  none in the moment
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
***)


BeginPackage["GroupTheory`Symbols`"]



   kx                       ::usage = "bla"  
   ky                       ::usage = "bla"  
   kz                       ::usage = "bla"  
   p                        ::usage = "bla"   
   P                        ::usage = "bla"   
   ci                       ::usage = "bla"  
   xswf                     ::usage = "internal variable in GTSymmetryBasisFunction"
   yswf                     ::usage = "internal variable in GTSymmetryBasisFunction"
   zswf                     ::usage = "internal variable in GTSymmetryBasisFunction" 
 (*  x                       ::usage = "bla"  
   y                       ::usage = "bla"  
   z                       ::usage = "bla" 	*)
   

(*$Context="GroupTheory`PointGroup`"*)

Off[General::"spell1"]; Off[General::"spell"];

(*--- general symbols for crystal field ---*)
Unprotect[A,Y,S,r]
Clear[A,Y,S,r]
Protect[A,Y,S,r]

(*--- Error Messages ---*)

Unprotect[$ErrorCrystalFieldA,$ErrorEAToPermMat,$ErrorEAToSU2,$ErrorChangeRepresentation,$ErrorDimensionKvec,$ErrorGTGetSubgroupsOrder,$ErrorGTGetIrep,$ErrorGTGetIrepRepr,$ErrorGTGetIrepCornwell,
$ErrorGTGetIrepCornwellMethod,$ErrorGTInstallGroup,$ErrorInput,$ErrorInputb,$ErrorInvCt,$ErrorNoGroup,$ErrorNoGroupFirst,$ErrorNoGroupSecond,
$ErrorNOQuotientGroup,$ErrorNoSymbol,$ErrorNumberIreps,$ErrorProjectionOperator,$ErrorProjectionNoSpinor,$ErrorProjectionNo2D,$ErrorProjectionNoPermutation,
$ErrorReinstallAxes,$ErrorTransformationOperator,$ErrorTransformationO3,$ErrorTransformationO2,$ErrorTransformationSU2,$ErrorTransformationNoPermutation,$ErrorOwn,$ErrorGTStevensTheta,$ErrorGTStevensThetaLval];

$ErrorChangeRepresentation="Error: Representation not known. Choose one of the following: O(3), SU(2), O(2), Permutation or SU(2)xS.";
$ErrorCrystalFieldA = "Error: number of elements in list of vectors does not match with the list of charges!"
$ErrorDimensionKvec = "Error: dimension of reciprocal basis and k-vector does not fit!"
$ErrorEAToPermMat = "Error: The standard representation is given by permutation matrices. Permutation matrices can not be evaluated from euler angles!"
$ErrorEAToSU2 = "Error: The standard representation is given by SU(2) matrices. SU(2) can not be evaluated from euler angles! This Problem will be fixed soon!"
$ErrorGTCrystalFieldSubGroup ="Error: group2 is not a subgroup of group1!"
$ErrorFindIrepInCT = "Error: GTGetIrep can't find the specified irreducible representation using the induction method. Please specify a different method."
$ErrorGTGetSubgroupsOrder ="Error: The requested order of the subgroup is larger than the order of the group!"
$ErrorGTGetIrep ="Error: The specified method is not known! Use either Flodmark or Cornwell."
$ErrorGTGetIrepRepr ="Error: The specified method only works for O(3) matrices or O(2) matrices. Use Flodmark's method by specifying the option Method."
$ErrorGTGetIrepCornwell ="Error: Unable to find matrices! Try to increase lmax by setting the option GOlmax."
$ErrorGTGetIrepCornwellMethod ="Error: Unable to find matrices! The projection is not unique. Use Flodmark's method by specifying the option Method."
(* $ErrorGTInstallGroup = "Error: the specified group is not known!" *)
$ErrorGTStevensTheta = "Error: Element not found! Specify one of the following: Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tu, Yb"
$ErrorGTStevensThetaLval = "Error: no data for specified l value. Specify one of the following: l = 0, 2, 4, 6"
$ErrorInput = "Error: Unidentified input! The input is neither a symbol, a matrix, a quaternion nor a set of Euler angles!";
$ErrorInputb = "Error: Input not known!";
$ErrorInvCt = "Error: Invalid character table!";
$ErrorNoGroup = "Error: Input is not a group!";
$ErrorNoGroupFirst = "Error: First input is not a group!";
$ErrorNoGroupSecond = "Error: Second input is not a group!";
$ErrorNOQuotientGroup = "Error: cCnnot build quotient group!"
$ErrorNoSymbol = "Error: No symbol can be allocated to the input!";
$ErrorNumberIreps = "Error: The irep index is greater than the number of irreducible representations!";
$ErrorProjectionOperator = "Error: The input function has an invalid dimension! The character projection operator is implemented for one dimensional, two dimension and three dimension functions as well as for spinor functions (if double groups are used)!";
$ErrorProjectionNoSpinor ="Error: Projection operators are only implemented for spinor functions, if double groups are used!";
$ErrorProjectionNo2D ="Error: Projection operators are not implemented for three dimensional functions, if 2D groups are used! Use O(3) representation!";
$ErrorProjectionNoPermutation ="Error: Projection operators are not implemented for permutation matrices. Use O(3), SU(2), O(2), or SU(2)xS!";
$ErrorReinstallAxes = "Error: You have to specify either active or passive.";
$ErrorGTSOCSplitting = "Error: double groups have to be used!";
$ErrorTransformationOperator = "Error: The input function has an invalid dimension! The transformation operator is implemented for one dimensional, two dimension and three dimension functions as well as for spinor functions (if double groups are used)!";
$ErrorTransformationArgumentDimension = "Error: dimension of arguments does not fit to the dimension of the standard representation!"
$ErrorTransformationO3 ="Error: The given function has to be a scalar function or a 3-dimensional vector field!";
$ErrorTransformationSU2 ="Error: The given function has to be a 2-component spinor function!";
$ErrorTransformationO2 ="Error: The given function has to be a scalar function or a 2-dimensional vector field!";
$ErrorTransformationNoPermutation ="Error: The transformation operator is not implemented for permutation matrices. Use O(3), SU(2) or O(2)!";
$ErrorOwn = "Error: only GTGetSymbol and GTGetMatrix are implemented for groups, installed using GTTableToGroup!";

Protect[$ErrorCrystalFieldA,$ErrorEAToPermMat,$ErrorEAToSU2,$ErrorChangeRepresentation,$ErrorDimensionKvec,$ErrorGTGetSubgroupsOrder,$ErrorGTGetIrep,$ErrorGTGetIrepRepr,$ErrorGTGetIrepCornwell,
$ErrorGTGetIrepCornwellMethod,$ErrorGTInstallGroup,$ErrorInput,$ErrorInputb,$ErrorInvCt,$ErrorNoGroup,$ErrorNoGroupFirst,$ErrorNoGroupSecond,
$ErrorNOQuotientGroup,$ErrorNoSymbol,$ErrorNumberIreps,$ErrorProjectionOperator,$ErrorProjectionNoSpinor,$ErrorProjectionNo2D,$ErrorProjectionNoPermutation,
$ErrorReinstallAxes,$ErrorTransformationOperator,$ErrorTransformationNoPermutation,$ErrorTransformationO3,$ErrorTransformationO2,$ErrorTransformationSU2,$ErrorOwn,$ErrorGTStevensTheta,$ErrorGTStevensThetaLval];

Unprotect[GOFastValue];
GOFastValue=False;
Protect[GOFastValue];

(*--- active or passive rotations; passive=1, active=-1 ---*)
Unprotect[gtactvpsv];
gtactvpsv = 1;
Protect[gtactvpsv];

(*--- Representation ---*)
Unprotect[grpdgrp];
Clear[grpdgrp];
grpdgrp=1;
Protect[grpdgrp];

(*--- magnetic ---*)
Unprotect[GTmaggrpq];
Clear[GTmaggrpq];
GTmaggrpq=False;
Protect[GTmaggrpq];

(*--- crystal structures ---*)
Unprotect[spcgrp];
Clear[spcgrp];
spcgrp={};
Protect[spcgrp];

(*--- Wyckoff positions ---*)
Unprotect[wyckoff];
Clear[wyckoff];
wyckoff={};
Protect[wyckoff];

Unprotect[GTnumberlist]
GTnumberlist = Table[ToString[i], {i, 1, 9}];
Protect[GTnumberlist]

(*--- Axes ---*)
Unprotect[ex,ey,ez,ea,eb,ec,ed,ee,ef,e\[Alpha],e\[Beta],e\[Gamma],eA,eB,eC,eD, axeslist,axistmp];
Clear[ex,ey,ez,ea,eb,ec,ed,ee,ef,e\[Alpha],e\[Beta],e\[Gamma],eA,eB,eC,eD, axeslist,axistmp];
ex={1,0,0};
ey={0,1,0};
ez={0,0,1};
ea={1/Sqrt[2],1/Sqrt[2],0};
eb={1/Sqrt[2],-(1/Sqrt[2]),0};
ec={1/Sqrt[2],0,1/Sqrt[2]};
ed={-(1/Sqrt[2]),0,1/Sqrt[2]};
ee={0,1/Sqrt[2],1/Sqrt[2]};
ef={0,1/Sqrt[2],-(1/Sqrt[2])};
e\[Alpha]={-(1/Sqrt[3]),-(1/Sqrt[3]),1/Sqrt[3]};
e\[Beta]={-(1/Sqrt[3]),1/Sqrt[3],-(1/Sqrt[3])};
e\[Gamma]={1/Sqrt[3],-(1/Sqrt[3]),-(1/Sqrt[3])};
e\[Delta]={1/Sqrt[3],1/Sqrt[3],1/Sqrt[3]};
eA={1/2,Sqrt[3]/2,0};
eB={-(1/2),Sqrt[3]/2,0};
eC={-(Sqrt[3]/2),1/2,0};
eD={-(Sqrt[3]/2),-(1/2),0};
axeslist ={}
Protect[ex,ey,ez,ea,eb,ec,ed,ee,ef,e\[Alpha],e\[Beta],e\[Gamma],eA,eB,eC,eD, axeslist,axistmp];

(*--- Symbols ---*)
Format[Ee,StandardForm]:="Ee";
Format[IEe,StandardForm]:="IEe";
Format[DEe,StandardForm]:=Overscript["Ee", "__"];
Format[DIEe,StandardForm]:=Overscript["IEe", "___"];

Unprotect[elm,elmo2,elmSU2,elmSU2xS,elmown,elmpalette,list0,Nelm,Nelmo2,NelmSU2,NelmSU2xS]
Clear[elm,elmo2,elmSU2,elmSU2xS,elmown,elmpalette,list0]
elm = {
{ToString[Ee],IdentityMatrix[3],{1,{0,0,0}},{{0,0,0},1},1},
{ToString[DEe],IdentityMatrix[3],{-1,{0,0,0}},{{0,0,0},I},-1},
{ToString[IEe],-IdentityMatrix[3],{1,{0,0,0}},{{0,0,0},-1},1},
{ToString[DIEe],-IdentityMatrix[3],{-1,{0,0,0}},{{0,0,0},-I},-1}
};

elmo2 = {
{ToString[Ee],IdentityMatrix[2],{1,{0,0,0}},{{0,0,0},1},1},
{ToString[DEe],IdentityMatrix[2],{-1,{0,0,0}},{{0,0,0},I},-1},
{ToString[IEe],-IdentityMatrix[2],{1,{0,0,0}},{{0,0,0},-1},1},
{ToString[DIEe],-IdentityMatrix[2],{-1,{0,0,0}},{{0,0,0},-I},-1}
};

elmSU2 = {
{ToString[Ee],IdentityMatrix[2],{1,{0,0,0}},{{0,0,0},1},1},
{ToString[DEe],-IdentityMatrix[2],{-1,{0,0,0}},{{0,0,0},I},-1},
{ToString[IEe],IdentityMatrix[2],{1,{0,0,0}},{{0,0,0},-1},1},
{ToString[DIEe],-IdentityMatrix[2],{-1,{0,0,0}},{{0,0,0},-I},-1}
};

elmSU2xS = {
{ToString[Ee],IdentityMatrix[3],{1,{0,0,0}},{{0,0,0},1},1},
{ToString[DEe],DiagonalMatrix[{-1,-1,1}],{-1,{0,0,0}},{{0,0,0},I},-1},
{ToString[IEe],DiagonalMatrix[{1,1,-1}],{1,{0,0,0}},{{0,0,0},-1},1},
{ToString[DIEe],DiagonalMatrix[{-1,-1,-1}],{-1,{0,0,0}},{{0,0,0},-I},-1}
};

elmown={};

elmpalette =
  {Ee,C2x,C2y,C2z,C2a,C2b,C2c,C2d,C2e,C2f,C2A,C2B,C2C,C2D,C3\[Alpha],C3\[Beta],C3\[Gamma],C3\[Delta],
  	C3\[Alpha]i,C3\[Beta]i,C3\[Gamma]i,C3\[Delta]i, C3z,C3zi,C4x,C4y,C4z,C4xi,C4yi,C4zi, C6z,C6zi,
   IEe,IC2x,IC2y,IC2z,IC2a,IC2b,IC2c,IC2d,IC2e,IC2f,IC2A,IC2B,IC2C,IC2D,IC3\[Alpha],IC3\[Beta],IC3\[Gamma],IC3\[Delta],
  	IC3\[Alpha]i,IC3\[Beta]i,IC3\[Gamma]i,IC3\[Delta]i, IC3z,IC3zi,IC4x,IC4y,IC4z,IC4xi,IC4yi,IC4zi, IC6z,IC6zi,
   DEe,DC2x,DC2y,DC2z,DC2a,DC2b,DC2c,DC2d,DC2e,DC2f,DC2A,DC2B,DC2C,DC2D,DC3\[Alpha],DC3\[Beta],DC3\[Gamma],DC3\[Delta],
  	DC3\[Alpha]i,DC3\[Beta]i,DC3\[Gamma]i,DC3\[Delta]i, DC3z,DC3zi,DC4x,DC4y,DC4z,DC4xi,DC4yi,DC4zi, DC6z,DC6zi, 
   DIEe,DIC2x,DIC2y,DIC2z,DIC2a,DIC2b,DIC2c,DIC2d,DIC2e,DIC2f,DIC2A,DIC2B,DIC2C,DIC2D,DIC3\[Alpha],DIC3\[Beta],DIC3\[Gamma],DIC3\[Delta],
  	DIC3\[Alpha]i,DIC3\[Beta]i,DIC3\[Gamma]i,DIC3\[Delta]i, DIC3z,DIC3zi,DIC4x,DIC4y,DIC4z,DIC4xi,DIC4yi,DIC4zi, DIC6z,DIC6zi};

list0 = Append[Union@Flatten[Table[
     mmax = (n - 1)/2;
     Table[
      Table[{axistmp, s*Numerator[n/m], Denominator[n/m], d}, {s, -1, 1, 
        2}], {m, 1, mmax}]
     , {n, 2, 9}], 2], {axistmp, 2, 1, d}]

Nelm = Rationalize[N[N[elm, 15], 8]]
Nelmo2 = Rationalize[N[N[elmo2, 15], 8]]
NelmSU2 = Rationalize[N[N[elmSU2, 15], 8]]
NelmSU2xS = Rationalize[N[N[elmSU2xS, 15], 8]]
     
Protect[elm,elmSU2,elmSU2xS,elmown,elmpalette,list0,Nelm,Nelmo2,NelmSU2,NelmSU2xS]     

(*--- Notation of point groups in Schoenflies form ---*)

Unprotect[C1h,S2,S6,Oh,Th,Td,D6h,D4h,C6h,C6v,D6,D3h,D3d,C4h,C4v,D4,D2h,D2d,Vd,Vh,S1,
	C6,C3i,C3h,C3v,D3,C4,S4,C2h,C2v,D2,Cs,Ci,C2,C1,C3]
	
Clear[C1h,S2,S6,Oh,Th,Td,D6h,D4h,C6h,C6v,D6,D3h,D3d,C4h,C4v,D4,D2h,D2d,Vd,Vh,S1,
	C6,C3i,C3h,C3v,D3,C4,S4,C2h,C2v,D2,Cs,Ci,C2,C1,C3]

Format[C1h,StandardForm]:=Subscript["C","1h"];
Format[S1,StandardForm]:=Subscript["S","1"];
Format[S2,StandardForm]:=Subscript["S","2"];
Format[S6,StandardForm]:=Subscript["S","6"];
Format[Oh,StandardForm]:=Subscript["O","h"];
Format[Th,StandardForm]:=Subscript["T","h"];
Format[Td,StandardForm]:=Subscript["T","d"];
Format[D6h,StandardForm]:=Subscript["D","6h"];
Format[D4h,StandardForm]:=Subscript["D","4h"];
Format[C6h,StandardForm]:=Subscript["C","6h"];
Format[C6v,StandardForm]:=Subscript["C","6v"];
Format[D6,StandardForm]:=Subscript["D","6"];
Format[D3h,StandardForm]:=Subscript["D","3h"];
Format[D3d,StandardForm]:=Subscript["D","3d"];
Format[C4h,StandardForm]:=Subscript["C","4h"];
Format[C4v,StandardForm]:=Subscript["C","4v"];
Format[D4,StandardForm]:=Subscript["D","4"];
Format[D2h,StandardForm]:=Subscript["D","2h"];
Format[D2d,StandardForm]:=Subscript["D","2d"];
Format[C6,StandardForm]:=Subscript["C","6"];
Format[C3i,StandardForm]:=Subscript["C","3i"];
Format[C3h,StandardForm]:=Subscript["C","3h"];
Format[C3v,StandardForm]:=Subscript["C","3v"];
Format[D3,StandardForm]:=Subscript["D","3"];
Format[C4,StandardForm]:=Subscript["C","4"];
Format[S4,StandardForm]:=Subscript["S","4"];
Format[C2h,StandardForm]:=Subscript["C","2h"];
Format[C2v,StandardForm]:=Subscript["C","2v"];
Format[D2,StandardForm]:=Subscript["D","2"];
Format[Cs,StandardForm]:=Subscript["C","s"];
Format[Ci,StandardForm]:=Subscript["C","i"];
Format[C2,StandardForm]:=Subscript["C","2"];
Format[C1,StandardForm]:=Subscript["C","1"];
Format[C3,StandardForm]:=Subscript["C","3"];
Format[Vd,StandardForm]:=Subscript["V","d"];
Format[Vh,StandardForm]:=Subscript["V","h"];

Protect[C1h,S2,S6,Oh,Th,Td,D6h,D4h,C6h,C6v,D6,D3h,D3d,C4h,C4v,D4,D2h,D2d,Vd,Vh,S1,
	C6,C3i,C3h,C3v,D3,C4,S4,C2h,C2v,D2,Cs,Ci,C2,C1,C3]
(*	
Unprotect[\[LeftBracketingBar]1\[RightBracketingBar],
	      \[LeftBracketingBar]-1\[RightBracketingBar],
	      \[LeftBracketingBar]2\[RightBracketingBar],
	      \[LeftBracketingBar]m\[RightBracketingBar],
	      \[LeftBracketingBar]3\[RightBracketingBar],
	      \[LeftBracketingBar]-4\[RightBracketingBar],
	      \[LeftBracketingBar]4\[RightBracketingBar],
	      \[LeftBracketingBar]2/m\[RightBracketingBar],
	      \[LeftBracketingBar]mm2\[RightBracketingBar],
	       \[LeftBracketingBar]6\[RightBracketingBar],
	       \[LeftBracketingBar]-3\[RightBracketingBar],
	       \[LeftBracketingBar]32\[RightBracketingBar],
	       \[LeftBracketingBar]3m\[RightBracketingBar],
	       \[LeftBracketingBar]-6\[RightBracketingBar],
	       \[LeftBracketingBar]4/m\[RightBracketingBar],
	       \[LeftBracketingBar]42m\[RightBracketingBar],
	       \[LeftBracketingBar]422\[RightBracketingBar],
	       \[LeftBracketingBar]4mm\[RightBracketingBar],
	       \[LeftBracketingBar]mmm\[RightBracketingBar],
	       \[LeftBracketingBar]23\[RightBracketingBar],
	       \[LeftBracketingBar]622\[RightBracketingBar],
	       \[LeftBracketingBar]6/m\[RightBracketingBar],
	       \[LeftBracketingBar]6mm\[RightBracketingBar],
	       \[LeftBracketingBar]-3m\[RightBracketingBar],
	       \[LeftBracketingBar]-6m2\[RightBracketingBar],
	       \[LeftBracketingBar]4/mmm\[RightBracketingBar],
	       \[LeftBracketingBar]m3\[RightBracketingBar],
	       \[LeftBracketingBar]43\[RightBracketingBar],
	       \[LeftBracketingBar]-43m\[RightBracketingBar],
	       \[LeftBracketingBar]6/mmm\[RightBracketingBar],
	       \[LeftBracketingBar]m3m\[RightBracketingBar]
]
	
	
Clear[    \[LeftBracketingBar]1\[RightBracketingBar],
	      \[LeftBracketingBar]-1\[RightBracketingBar],
	      \[LeftBracketingBar]2\[RightBracketingBar],
	      \[LeftBracketingBar]m\[RightBracketingBar],
	      \[LeftBracketingBar]3\[RightBracketingBar],
	      \[LeftBracketingBar]-4\[RightBracketingBar],
	      \[LeftBracketingBar]4\[RightBracketingBar],
	      \[LeftBracketingBar]222\[RightBracketingBar],
	       \[LeftBracketingBar]2/m\[RightBracketingBar],
	       \[LeftBracketingBar]mm2\[RightBracketingBar],
	       \[LeftBracketingBar]6\[RightBracketingBar],
	       \[LeftBracketingBar]-3\[RightBracketingBar],
	       \[LeftBracketingBar]32\[RightBracketingBar],
	       \[LeftBracketingBar]3m\[RightBracketingBar],
	       \[LeftBracketingBar]-6\[RightBracketingBar],
	       \[LeftBracketingBar]4/m\[RightBracketingBar],
	       \[LeftBracketingBar]42m\[RightBracketingBar],
	       \[LeftBracketingBar]422\[RightBracketingBar],
	       \[LeftBracketingBar]4mm\[RightBracketingBar],
	       \[LeftBracketingBar]mmm\[RightBracketingBar],
	       \[LeftBracketingBar]23\[RightBracketingBar],
	       \[LeftBracketingBar]622\[RightBracketingBar],
	       \[LeftBracketingBar]6/m\[RightBracketingBar],
	       \[LeftBracketingBar]6mm\[RightBracketingBar],
	       \[LeftBracketingBar]-3m\[RightBracketingBar],
	       \[LeftBracketingBar]-6m2\[RightBracketingBar],
	       \[LeftBracketingBar]4/mmm\[RightBracketingBar],
	       \[LeftBracketingBar]m3\[RightBracketingBar],
	       \[LeftBracketingBar]43\[RightBracketingBar],
	       \[LeftBracketingBar]-43m\[RightBracketingBar],
	       \[LeftBracketingBar]6/mmm\[RightBracketingBar],
	       \[LeftBracketingBar]m3m\[RightBracketingBar]
]
	*)	
Format[\[LeftBracketingBar]1\[RightBracketingBar],StandardForm]:= "1";
Format[\[LeftBracketingBar]-1\[RightBracketingBar],StandardForm]:= GTOverScript[1];
Format[\[LeftBracketingBar]2\[RightBracketingBar],StandardForm]:= "2";
Format[\[LeftBracketingBar]m\[RightBracketingBar],StandardForm]:= "m";
Format[\[LeftBracketingBar]2/m\[RightBracketingBar],StandardForm]:= "2/m";
Format[\[LeftBracketingBar]222\[RightBracketingBar],StandardForm]:= "222";
Format[\[LeftBracketingBar]mm2\[RightBracketingBar],StandardForm]:= "mm2";
Format[\[LeftBracketingBar]mmm\[RightBracketingBar],StandardForm]:= "mmm";
Format[\[LeftBracketingBar]4\[RightBracketingBar],StandardForm]:= "4";
Format[\[LeftBracketingBar]-4\[RightBracketingBar],StandardForm]:= GTOverScript[4];
Format[\[LeftBracketingBar]4/m\[RightBracketingBar],StandardForm]:= "4/m";
Format[\[LeftBracketingBar]422\[RightBracketingBar],StandardForm]:= "422";
Format[\[LeftBracketingBar]4mm\[RightBracketingBar],StandardForm]:= "4mm";
Format[\[LeftBracketingBar]-42m\[RightBracketingBar],StandardForm]:= GTOverScript[4]<>"2m";
Format[\[LeftBracketingBar]4/mmm\[RightBracketingBar],StandardForm]:= "4/mmm";
Format[\[LeftBracketingBar]3\[RightBracketingBar],StandardForm]:= "3";
Format[\[LeftBracketingBar]-3\[RightBracketingBar],StandardForm]:= GTOverScript[3];
Format[\[LeftBracketingBar]32\[RightBracketingBar],StandardForm]:= "32";
Format[\[LeftBracketingBar]3m\[RightBracketingBar],StandardForm]:= "3m";
Format[\[LeftBracketingBar]-3m\[RightBracketingBar],StandardForm]:= GTOverScript[3]<>"m";
Format[\[LeftBracketingBar]6\[RightBracketingBar],StandardForm]:= "6";
Format[\[LeftBracketingBar]-6\[RightBracketingBar],StandardForm]:= GTOverScript[6];
Format[\[LeftBracketingBar]6/m\[RightBracketingBar],StandardForm]:= "6/m";
Format[\[LeftBracketingBar]622\[RightBracketingBar],StandardForm]:= "622";
Format[\[LeftBracketingBar]6mm\[RightBracketingBar],StandardForm]:= "6mm";
Format[\[LeftBracketingBar]-6m2\[RightBracketingBar],StandardForm]:= GTOverScript[6] <> "m2";
Format[\[LeftBracketingBar]6/mmm\[RightBracketingBar],StandardForm]:= "6/mmm";
Format[\[LeftBracketingBar]23\[RightBracketingBar],StandardForm]:= "23";
Format[\[LeftBracketingBar]m-3\[RightBracketingBar],StandardForm]:= "m"<>GTOverScript[3];
Format[\[LeftBracketingBar]432\[RightBracketingBar],StandardForm]:= "432";
Format[\[LeftBracketingBar]-43m\[RightBracketingBar],StandardForm]:= GTOverScript[4] <> "3m";
Format[\[LeftBracketingBar]m-3m\[RightBracketingBar],StandardForm]:= "m"<>GTOverScript[3]<>"m";



(*	
Protect[\[LeftBracketingBar]1\[RightBracketingBar],
	      \[LeftBracketingBar]-1\[RightBracketingBar],
	      \[LeftBracketingBar]2\[RightBracketingBar],
	      \[LeftBracketingBar]m\[RightBracketingBar],
	      \[LeftBracketingBar]3\[RightBracketingBar],
	      \[LeftBracketingBar]-4\[RightBracketingBar],
	      \[LeftBracketingBar]4\[RightBracketingBar],
	      \[LeftBracketingBar]2/m\[RightBracketingBar],
	      \[LeftBracketingBar]mm2\[RightBracketingBar],
	       \[LeftBracketingBar]6\[RightBracketingBar],
	       \[LeftBracketingBar]-3\[RightBracketingBar],
	       \[LeftBracketingBar]32\[RightBracketingBar],
	       \[LeftBracketingBar]3m\[RightBracketingBar],
	       \[LeftBracketingBar]-6\[RightBracketingBar],
	       \[LeftBracketingBar]4/m\[RightBracketingBar],
	       \[LeftBracketingBar]42m\[RightBracketingBar],
	       \[LeftBracketingBar]422\[RightBracketingBar],
	       \[LeftBracketingBar]4mm\[RightBracketingBar],
	       \[LeftBracketingBar]mmm\[RightBracketingBar],
	       \[LeftBracketingBar]23\[RightBracketingBar],
	       \[LeftBracketingBar]622\[RightBracketingBar],
	       \[LeftBracketingBar]6/m\[RightBracketingBar],
	       \[LeftBracketingBar]6mm\[RightBracketingBar],
	       \[LeftBracketingBar]-3m\[RightBracketingBar],
	       \[LeftBracketingBar]-6m2\[RightBracketingBar],
	       \[LeftBracketingBar]4/mmm\[RightBracketingBar],
	       \[LeftBracketingBar]m3\[RightBracketingBar],
	       \[LeftBracketingBar]43\[RightBracketingBar],
	       \[LeftBracketingBar]-43m\[RightBracketingBar],
	       \[LeftBracketingBar]6/mmm\[RightBracketingBar],
	       \[LeftBracketingBar]m3m\[RightBracketingBar]
]
*)
(* *)

(*rRule=r -> Sqrt[x^2+y^2+z^2]; 

Protect[Ee, C3\[Alpha], C3\[Beta], C3\[Gamma], C3\[Delta], C3\[Alpha]i, 
    C3\[Beta]i, C3\[Gamma]i, C3\[Delta]i, C2x, C2y, C2z, C4x, C4y, C4z, C4xi, 
    C4yi, C4zi, C2a, C2b, C2c, C2d, C2e, C2f, C3z, C3zi, C6z, C6zi, C2A, C2B, 
    C2C, C2D, IEe, IC3\[Alpha], IC3\[Beta], IC3\[Gamma], IC3\[Delta], 
    IC3\[Alpha]i, IC3\[Beta]i, IC3\[Gamma]i, IC3\[Delta]i, IC2x, IC2y, IC2z, 
    IC4x, IC4y, IC4z, IC4xi, IC4yi, IC4zi, IC2a, IC2b, IC2c, IC2d, IC2e, IC2f,
     IC3z, IC3zi, IC6z, IC6zi, IC2A, IC2B, IC2C, IC2D,
    rRule];

Protect[DEe, DC3\[Alpha], DC3\[Beta], DC3\[Gamma], DC3\[Delta], DC3\[Alpha]i, 
    DC3\[Beta]i, DC3\[Gamma]i, DC3\[Delta]i, DC2x, DC2y, DC2z, DC4x, DC4y, DC4z, DC4xi, 
    DC4yi, DC4zi, DC2a, DC2b, DC2c, DC2d, DC2e, DC2f, DC3z, DC3zi, DC6z, DC6zi, DC2A, DC2B, 
    DC2C, DC2D, DIEe, DIC3\[Alpha], DIC3\[Beta], DIC3\[Gamma], DIC3\[Delta], 
    DIC3\[Alpha]i, DIC3\[Beta]i, DIC3\[Gamma]i, DIC3\[Delta]i, DIC2x, DIC2y, DIC2z, 
    DIC4x, DIC4y, DIC4z, DIC4xi, DIC4yi, DIC4zi, DIC2a, DIC2b, DIC2c, DIC2d, DIC2e, DIC2f,
    DIC3z, DIC3zi, DIC6z, DIC6zi, DIC2A, DIC2B, DIC2C, DIC2D];
*)
Unprotect[xez]
Clear[xez]
xez = {\[Xi], \[Eta], \[Zeta]}
Protect[xez]

Unprotect[GTcttable]
Clear[GTcttable]
GTcttable = {
    (*--- C1 ---*)
    {{1, 1}, {{{101/100}, {{1}}, {"A"}, {"1"}}}},
    (*--- Cs, Ci, C2 ---*)
    {{2, -2}, {
      (*--- Cs ---*)
      {{101/
        100, -(199/100)}, {{1, 1}, {1, -1}}, {Superscript["A", "'"], 
        Superscript["A", "''"]}, {"1", "2"}},
      (*--- Ci ---*)
      {{101/
        100, -(1989/1000)}, {{1, 1}, {1, -1}}, {Subscript["A", "g"], 
        Subscript["A", "u"]}, {"1", "2"}},
      (*--- C2 ---*)
      {{101/100, 201/100 + I}, {{1, 1}, {1, -1}}, {"A", 
        "B"}, {"1", "2"}}}
     },
    (*--- C3 ---*)
    {{3, 
      3 I Sqrt[
       3]}, {{{101/100, 301/100 + I/100, 301/
        100}, {{1, 1, 1}, {1, -(1/2) I (-I + Sqrt[3]), 
         1/2 I (I + Sqrt[3])}, {1, 
         1/2 I (I + Sqrt[3]), -(1/2) I (-I + Sqrt[3])}}, {"A", 
        Superscript["E", "1"], Superscript["E", "2"]}, {"1", "3", 
        "2"}}}},
    (*--- C4, S4 ---*)
    {{4, 16 I}, {
      (*--- C4 ---*)
      {{101/100, 201/100 + I, 401/100 + I/100, 401/
        100}, {{1, 1, 1, 1}, {1, -1, -I, I}, {1, -1, I, -I}, {1, 
         1, -1, -1}}, {"A", Superscript["E", "1"], 
        Superscript["E", "2"], "B"}, {"1", "3", "4", "2"}},
      (*--- S4 ---*)
      {{101/100, 
        201/100 + I, -(399/100) + I/100, -(399/100)}, {{1, 1, 1, 
         1}, {1, -1, -I, I}, {1, -1, I, -I}, {1, 1, -1, -1}}, {"A", 
        Superscript["E", "1"], Superscript["E", "2"], "B"}, {"1", "3",
         "4", "2"}}}
     },
    (*--- C2v, C2h, D2 ---*)
    {{4, 16}, {
      (*--- C2v ---*)
      {{101/100, 
        201/100 + I, -(19899/10000), -(199/100)}, {{1, 1, 1, 
         1}, {1, -1, -1, 1}, {1, -1, 1, -1}, {1, 
         1, -1, -1}}, {Subscript["A", "1"], Subscript["B", "2"], 
        Subscript["B", "1"], Subscript["A", "2"]}, {"1", "3", "4", 
        "2"}},
      (*--- C2h ---*)
      {{101/100, -(1989/1000), 
        201/100 + I, -(199/100)}, {{1, 1, 1, 1}, {1, -1, -1, 
         1}, {1, -1, 1, -1}, {1, 1, -1, -1}}, {Subscript["A", "g"], 
        Subscript["B", "u"], Subscript["A", "u"], 
        Subscript["B", "g"]}, {"1", "3", "4", "2"}},
      (*--- D2 ---*)
      {{101/100, 201/100 + 3 I, 201/100 + 2 I, 
        201/100 + I}, {{1, 1, 1, 1}, {1, -1, -1, 1}, {1, -1, 
         1, -1}, {1, 1, -1, -1}}, {Subscript["A", "1"], 
        Subscript["B", "1"], Subscript["B", "2"], 
        Subscript["B", "3"]}, {"1", "3", "4", "2"}}}
     },
    (*--- C3v, D3 ---*)
    {{3, -6}, {
      (*--- C3v ---*)
      {{101/100, 151/
        50, -(197/100)}, {{1, 1, 1}, {1, 1, -1}, {2, -1, 
         0}}, {Subscript["A", "1"], Subscript["A", "2"], "E"}, {"1", 
        "2", "3"}},
      (*--- D3 ---*)
      {{101/100, 151/50, 
        203/100 + I}, {{1, 1, 1}, {1, 1, -1}, {2, -1, 0}}, {Subscript[
         "A", "1"], Subscript["A", "2"], "E"}, {"1", "2", "3"}}}
     },
    (*--- C3h, C3i, C6 ---*)
    {{6, 216}, {
      (*--- C3h ---*)
      {{101/100, -(599/100), 
        301/100 + I/100, -(199/100), 301/
        100, -(599/100) + I/100}, {{1, 1, 1, 1, 1, 1}, {1, -1, 1, -1, 
         1, -1}, {1, 1/2 (1 - I Sqrt[3]), -(1/2) I (-I + Sqrt[3]), -1,
          1/2 I (I + Sqrt[3]), 1/2 (1 + I Sqrt[3])}, {1, 
         1/2 (1 + I Sqrt[3]), 
         1/2 I (I + Sqrt[3]), -1, -(1/2) I (-I + Sqrt[3]), 
         1/2 (1 - I Sqrt[3])}, {1, -(1/2) I (-I + Sqrt[3]), 
         1/2 I (I + Sqrt[3]), 1, -(1/2) I (-I + Sqrt[3]), 
         1/2 I (I + Sqrt[3])}, {1, 
         1/2 I (I + Sqrt[3]), -(1/2) I (-I + Sqrt[3]), 1, 
         1/2 I (I + Sqrt[3]), -(1/2) I (-I + Sqrt[3])}}, {Superscript[
         "A", "'"], Superscript["A", "''"], Superscript["E", "''1"], 
        Superscript["E", "''2"], Superscript["E", "'1"], 
        Superscript["E", "'2"]}, {"1", "2", "4", "3", "6", "5"}},
      (*--- C3i ---*)
      {{101/100, -(599/100) + I/100, 
        301/100 + I/100, -(1989/1000), 301/
        100, -(599/100)}, {{1, 1, 1, 1, 1, 1}, {1, -1, 1, -1, 
         1, -1}, {1, 1/2 (1 + I Sqrt[3]), 
         1/2 I (I + Sqrt[3]), -1, -(1/2) I (-I + Sqrt[3]), 
         1/2 (1 - I Sqrt[3])}, {1, 
         1/2 (1 - I Sqrt[3]), -(1/2) I (-I + Sqrt[3]), -1, 
         1/2 I (I + Sqrt[3]), 
         1/2 (1 + I Sqrt[3])}, {1, -(1/2) I (-I + Sqrt[3]), 
         1/2 I (I + Sqrt[3]), 1, -(1/2) I (-I + Sqrt[3]), 
         1/2 I (I + Sqrt[3])}, {1, 
         1/2 I (I + Sqrt[3]), -(1/2) I (-I + Sqrt[3]), 1, 
         1/2 I (I + Sqrt[3]), -(1/2) I (-I + Sqrt[3])}}, {Subscript[
         "A", "g"], Subscript["A", "u"], 
        Subsuperscript["E", "u", "2"], Subsuperscript["E", "u", "1"], 
        Subsuperscript["E", "g", "1"], 
        Subsuperscript["E", "g", "2"]}, {"1", "2", "4", "3", "6", 
        "5"}},
      (*--- C6 ---*)
      {{101/100, 601/100 + I/100, 301/100 + I/100, 
        201/100 + I, 301/100, 601/
        100}, {{1, 1, 1, 1, 1, 1}, {1, -1, 1, -1, 1, -1}, {1, 
         1/2 (1 + I Sqrt[3]), 
         1/2 I (I + Sqrt[3]), -1, -(1/2) I (-I + Sqrt[3]), 
         1/2 (1 - I Sqrt[3])}, {1, 
         1/2 (1 - I Sqrt[3]), -(1/2) I (-I + Sqrt[3]), -1, 
         1/2 I (I + Sqrt[3]), 
         1/2 (1 + I Sqrt[3])}, {1, -(1/2) I (-I + Sqrt[3]), 
         1/2 I (I + Sqrt[3]), 1, -(1/2) I (-I + Sqrt[3]), 
         1/2 I (I + Sqrt[3])}, {1, 
         1/2 I (I + Sqrt[3]), -(1/2) I (-I + Sqrt[3]), 1, 
         1/2 I (I + Sqrt[3]), -(1/2) I (-I + Sqrt[3])}}, {"A", "B", 
        Superscript["E", "'2"], Superscript["E", "'1"], 
        Superscript["E", "''1"], Superscript["E", "''2"]}, {"1", "2", 
        "4", "3", "6", "5"}}}
     },
    (*--- C4h ---*)
    {{8, -4096}, {
      {{101/100, -(1989/1000), 201/100 + 3 I, -(399/100), 
        401/100 + 2 I, -(399/100), 
        401/100 + I, -(199/100)}, {{1, 1, 1, 1, 1, 1, 1, 
         1}, {1, -1, -1, -I, -I, I, I, 1}, {1, -1, -1, I, I, -I, -I, 
         1}, {1, -1, 1, -1, 1, -1, 1, -1}, {1, -1, 1, 1, -1, 
         1, -1, -1}, {1, 1, -1, -I, I, I, -I, -1}, {1, 1, -1, 
         I, -I, -I, I, -1}, {1, 1, 1, -1, -1, -1, -1, 1}}, {Subscript[
         "A", "g"], Subsuperscript["E", "u", "1"], 
        Subsuperscript["E", "u", "2"], Subscript["A", "u"], 
        Subscript["B", "u"], Subsuperscript["E", "g", "2"], 
        Subsuperscript["E", "g", "1"], Subscript["B", "g"]}, {"1", 
        "7", "8", "5", "6", "4", "3", "2"}}}},
    (*--- C4v, D4, D2d ---*)
    {{5, -64}, {
      (*--- C4v ---*)
      {{101/100, -(19799/10000), -(99/50), 201/50, 
        201/100 + I}, {{1, 1, 1, 1, 1}, {1, -1, -1, 1, 1}, {1, -1, 
         1, -1, 1}, {1, 1, -1, -1, 1}, {2, 0, 0, 0, -2}}, {Subscript[
         "A", "1"], Subscript["A", "2"], Subscript["B", "2"], 
        Subscript["B", "1"], "E"}, {"1", "2", Superscript["1", "'"], 
        Superscript["2", "'"], "3"}},
      (*--- D4 ---*)
      {{101/100, 101/50 + 2 I, 101/50 + 3 I, 201/50, 
        201/100 + I}, {{1, 1, 1, 1, 1}, {1, -1, -1, 1, 1}, {1, -1, 
         1, -1, 1}, {1, 1, -1, -1, 1}, {2, 0, 0, 0, -2}}, {Subscript[
         "A", "1"], Subscript["A", "2"], Subscript["B", "2"], 
        Subscript["B", "1"], "E"}, {"1", "2", Superscript["1", "'"], 
        Superscript["2", "'"], "3"}},
      (*--- D2d ---*)
      {{101/100, 101/50 + 2 I, -(199/50), -(99/50), 
        201/100 + I}, {{1, 1, 1, 1, 1}, {1, -1, -1, 1, 1}, {1, -1, 
         1, -1, 1}, {1, 1, -1, -1, 1}, {2, 0, 0, 0, -2}}, {Subscript[
         "A", "1"], Subscript["B", "1"], Subscript["A", "2"], 
        Subscript["B", "2"], "E"}, {"1", Superscript["2", "'"], 
        "2",Superscript["1", "'"], "3"}}}
     },
    (*--- D2h ---*)
    {{8, 4096}, {
      {{101/100, -(1989/1000), 201/100 + 3 I, 
        201/100 + 2 I, -(19899/10000), 
        201/100 + I, -(199/100), -(199/100) + I/100}, {{1, 1, 1, 1, 1,
          1, 1, 1}, {1, -1, -1, -1, -1, 1, 1, 1}, {1, -1, -1, 1, 
         1, -1, -1, 1}, {1, -1, 1, -1, 1, -1, 1, -1}, {1, -1, 1, 
         1, -1, 1, -1, -1}, {1, 1, -1, -1, 1, 1, -1, -1}, {1, 1, -1, 
         1, -1, -1, 1, -1}, {1, 1, 1, -1, -1, -1, -1, 1}}, {Subscript[
         "A", "1g"], Subscript["B", "3u"], Subscript["B", "1u"], 
        Subscript["B", "2u"], Subscript["A", "1u"], 
        Subscript["B", "3g"], Subscript["B", "1g"], 
        Subscript["B", "2g"]}, {"1", Superscript["4", "'"], 
        Superscript["1", "'"], Superscript["3", "'"], 
        Superscript["2", "'"], "3", "2", "4"}}}
     },
    (*--- T ---*)
    {{4, -12 I Sqrt[3]}, {
      {{101/100, 203/100 + I, 30401/10000, 76/
        25}, {{1, 1, 1, 1}, {1, 1, -(1/2) I (-I + Sqrt[3]), 
         1/2 I (I + Sqrt[3])}, {1, 1, 
         1/2 I (I + Sqrt[3]), -(1/2) I (-I + Sqrt[3])}, {3, -1, 0, 
         0}}, {"A", Superscript["E", "1"], Superscript["E", "2"], 
        "T"}, {"1", "3", "2", "4"}}
      }},
    (*--- D3h, D3d, C6v, D6 ---*)
    {{6, 288}, {
      (*--- D3h ---*)
      {{101/100, -(299/50), 151/50, 
        203/100 + I, -(197/100), -(199/100)}, {{1, 1, 1, 1, 1, 
         1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, 1, -1, -1}, {1, 1, 
         1, -1, -1, 1}, {2, -1, -1, 0, 0, 2}, {2, 1, -1, 0, 
         0, -2}}, {Subsuperscript["A", "1", "'"], 
        Subsuperscript["A", "2", "''"], 
        Subsuperscript["A", "1", "''"], Subsuperscript["A", "2", "'"],
         Superscript["E", "'"], Superscript["E", "''"]}, {"1", 
        Superscript["2", "'"], Superscript["1", "'"], "2", "3", 
        Superscript["3", "'"]}},
      (*--- D3d ---*)
      {{101/100, -(299/50), 151/50, -(197/100), 
        203/100 + I, -(1989/1000)}, {{1, 1, 1, 1, 1, 1}, {1, -1, 
         1, -1, 1, -1}, {1, -1, 1, 1, -1, -1}, {1, 1, 1, -1, -1, 
         1}, {2, -1, -1, 0, 0, 2}, {2, 1, -1, 0, 0, -2}}, {Subscript[
         "A", "1g"], Subscript["A", "1u"], Subscript["A", "2u"], 
        Subscript["A", "2g"], Subscript["E", "g"], 
        Subscript["E", "u"]}, {"1", Superscript["1", "'"], 
        Superscript["2", "'"], "2", "3", Superscript["3", "'"]}},
      (*--- C6v ---*)
      {{101/100, 301/50, 151/
        50, -(19699/10000), -(197/100), 
        201/100 + I}, {{1, 1, 1, 1, 1, 1}, {1, -1, 1, -1, 
         1, -1}, {1, -1, 1, 1, -1, -1}, {1, 1, 1, -1, -1, 
         1}, {2, -1, -1, 0, 0, 2}, {2, 1, -1, 0, 0, -2}}, {Subscript[
         "A", "1"], Subscript["B", "2"], Subscript["B", "1"], 
        Subscript["A", "2"], Subscript["E", "2"], 
        Subscript["E", "1"]}, {"1", Superscript["2", "'"], 
        Superscript["1", "'"], "2", "3", Superscript["3", "'"]}},
      (*--- D6 ---*)
      {{101/100, 301/50, 151/50, 203/100 + I, 
        203/100 + 3 I, 
        201/100 + 2 I}, {{1, 1, 1, 1, 1, 1}, {1, -1, 1, -1, 
         1, -1}, {1, -1, 1, 1, -1, -1}, {1, 1, 1, -1, -1, 
         1}, {2, -1, -1, 0, 0, 2}, {2, 1, -1, 0, 0, -2}}, {Subscript[
         "A", "1"], Subscript["B", "1"], Subscript["B", "2"], 
        Subscript["A", "2"], Subscript["E", "2"], 
        Subscript["E", "1"]}, {"1", Superscript["2", "'"], 
        Superscript["1", "'"], "2", "3", Superscript["3", "'"]}}}},
    (*--- D4h ---*)
    {{10, 131072}, {{
       {101/100, 101/50 + I, -(19799/10000), 
        101/50 + 3 I, -(99/50), -(199/50), 201/50, -(1989/1000), 
        201/100 + 2 I, -(199/100)}, {{1, 1, 1, 1, 1, 1, 1, 1, 1, 
         1}, {1, -1, -1, -1, -1, 1, 1, 1, 1, 1}, {1, -1, -1, 1, 
         1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, 1, -1, 1, -1, 
         1, -1}, {1, -1, 1, 1, -1, 1, -1, -1, 1, -1}, {1, 1, -1, -1, 
         1, 1, -1, -1, 1, -1}, {1, 1, -1, 1, -1, -1, 1, -1, 
         1, -1}, {1, 1, 1, -1, -1, -1, -1, 1, 1, 1}, {2, 0, 0, 0, 0, 
         0, 0, -2, -2, 2}, {2, 0, 0, 0, 0, 0, 0, 
         2, -2, -2}}, {Subscript["A", "1g"], Subscript["A", "2g"], 
        Subscript["B", "2g"], Subscript["A", "2u"], 
        Subscript["B", "2u"], Subscript["B", "1u"], 
        Subscript["A", "1u"], Subscript["B", "1g"], 
        Subscript["E", "u"], Subscript["E", "g"]}, {"1", "4", "3", 
        Superscript["4", "'"], Superscript["3", "'"], 
        Superscript["2", "'"], Superscript["1", "'"], "2", 
        Superscript["5", "'"], "5"}}}
     },
    (*--- Th ---*)
    {{8, -6912}, {{{101/100, 
        203/100 + I, -(197/100), -(59599/10000), 30401/
        10000, -(149/25), 76/
        25, -(1989/1000)}, {{1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, -1, -1, 
         1, -1, 1, -1}, {1, 1, -1, 1/2 (1 - I Sqrt[3]), 
         1/2 I (I + Sqrt[3]), 
         1/2 (1 + I Sqrt[3]), -(1/2) I (-I + Sqrt[3]), -1}, {1, 1, -1,
          1/2 (1 + I Sqrt[3]), -(1/2) I (-I + Sqrt[3]), 
         1/2 (1 - I Sqrt[3]), 1/2 I (I + Sqrt[3]), -1}, {1, 1, 
         1, -(1/2) I (-I + Sqrt[3]), -(1/2) I (-I + Sqrt[3]), 
         1/2 I (I + Sqrt[3]), 1/2 I (I + Sqrt[3]), 1}, {1, 1, 1, 
         1/2 I (I + Sqrt[3]), 
         1/2 I (I + Sqrt[3]), -(1/2) I (-I + Sqrt[3]), -(1/2)*I*(-I + Sqrt[3]), 1}, {3, -1, -1, 0, 0, 0, 0, 3}, {3, -1, 
         1, 0, 0, 0, 0, -3}}, {Subscript["A", "g"], 
        Subscript["A", "u"], Subsuperscript["E", "u", "1"], 
        Subsuperscript["E", "u", "2"], Subsuperscript["E", "g", "2"], 
        Subsuperscript["E", "g", "1"], Subscript["T", "g"], 
        Subscript["T", "u"]}, {"1", "5", "6", "7", "3", "2", "4", 
        "8"}}}},
    (*--- Td, O ---*)
    {{5, -96}, {
      (*--- Td ---*)
      {{101/100, 203/100 + I, -(197/50), -(97/50), 77/
        25}, {{1, 1, 1, 1, 1}, {1, 1, -1, -1, 1}, {2, 2, 0, 
         0, -1}, {3, -1, -1, 1, 0}, {3, -1, 1, -1, 0}}, {Subscript[
         "A", "1"], Subscript["A", "2"], "E", Subscript["T", "2"], 
        Subscript["T", "1"]}, {"1", "2", "3", "4", "5"}},
      (*--- O ---*)
      {{101/100, 203/100 + I, 103/50 + 2 I, 203/50, 77/
        25}, {{1, 1, 1, 1, 1}, {1, 1, -1, -1, 1}, {2, 2, 0, 
         0, -1}, {3, -1, -1, 1, 0}, {3, -1, 1, -1, 0}}, {Subscript[
         "A", "1"], Subscript["A", "2"], "E", Subscript["T", "1"], 
        Subscript["T", "2"]}, {"1", "2", "3", "4", "5"}}}},
    (*--- D6h ---*)
    {{12, 
      5308416}, {{{101/100, -(59799/10000), 151/50, -(299/50), 301/50,
         203/100 + I, -(19699/10000), 
        203/100 + 3 I, -(197/100), -(1989/1000), 
        201/100 + 2 I, -(199/100)}, {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1}, {1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1}, {1, -1, 
         1, -1, 1, 1, -1, 1, -1, -1, 1, -1}, {1, -1, 1, 1, -1, -1, 1, 
         1, -1, 1, -1, -1}, {1, -1, 1, 1, -1, 1, -1, -1, 1, 
         1, -1, -1}, {1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1}, {1, 
         1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1}, {1, 1, 1, 1, 
         1, -1, -1, -1, -1, 1, 1, 1}, {2, -1, -1, -1, -1, 0, 0, 0, 0, 
         2, 2, 2}, {2, -1, -1, 1, 1, 0, 0, 0, 0, -2, -2, 2}, {2, 
         1, -1, -1, 1, 0, 0, 0, 0, 2, -2, -2}, {2, 1, -1, 1, -1, 0, 0,
          0, 0, -2, 2, -2}}, {Subscript["A", "1g"], 
        Subscript["A", "2u"], Subscript["A", "1u"], 
        Subscript["B", "1g"], Subscript["B", "2g"], 
        Subscript["B", "1u"], Subscript["B", "2u"], 
        Subscript["A", "2g"], Subscript["E", "2g"], 
        Subscript["E", "1u"], Subscript["E", "1g"], 
        Subscript["E", "2u"]}, {"1", "8", "7", "3", "4", "9", "10", 
        "2", "5", "12", "6", "11"}}}},
    (*--- Oh ---*)
{{10,294912}, 
	{
	{{101/100, 203/100 + I, -(197/100),103/50 + 2 I, -(197/50), -(97/50), 203/50, -(148/25), 77/25, -(1989/1000)},
     {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {1, 1, -1, -1,1, 1, -1, -1, 1, -1}, {1, 1, -1, 1, -1, -1, 1, -1, 1, -1}, {1, 
     1, 1, -1, -1, -1, -1, 1, 1, 1}, {2, 2, -2, 0, 0, 0, 0, 
     1, -1, -2}, {2, 2, 2, 0, 0, 0, 0, -1, -1, 2}, {3, -1, -1, -1, 
     1, -1, 1, 0, 0, 3}, {3, -1, -1, 1, -1, 1, -1, 0, 0, 3}, {3, -1, 
     1, -1, -1, 1, 1, 0, 0, -3}, {3, -1, 1, 1, 1, -1, -1, 0, 
     0, -3}}, {Subscript["A", "1g"], Subscript["A", "2u"], 
    Subscript["A", "1u"], Subscript["A", "2g"], Subscript["E", "u"], 
    Subscript["E", "g"], Subscript["T", "1g"], Subscript["T", "2g"], 
    Subscript["T", "1u"], Subscript["T", "2u"]}, {"1", 
    Superscript["2", "'"], Superscript["1", "'"], "2", 
    Superscript["12", "'"], "12", Superscript["15", "'"], 
    Superscript["25", "'"], "15", "25"}}}}
    } // Simplify;
Protect[GTcttable]


(*--- begin spacegroup notations ---*)
(*--- Triclinic ---*)
Format[C11, StandardForm] := GTSF[C, 1, 1];
Format[\[LeftBracketingBar]Verbatim[P1]\[RightBracketingBar], StandardForm] := "P1";
Format[Ci1, StandardForm] := GTSF[C, i, 1];
Format[\[LeftBracketingBar]Verbatim[P - 1]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[1];

(*--- Monoclinic ---*)
Format[C21, StandardForm] := GTSF[C, 2, 1];
Format[\[LeftBracketingBar]Verbatim[P2]\[RightBracketingBar], StandardForm] := "P2";
Format[C22, StandardForm] := GTSF[C, 2, 2];
Format[\[LeftBracketingBar]Verbatim[P2_ 1]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[2, 1];
Format[C23, StandardForm] := GTSF[C, 2, 3];
Format[\[LeftBracketingBar]Verbatim[C2]\[RightBracketingBar], StandardForm] := "C2";
Format[Cs1, StandardForm] := GTSF[C, s, 1];
Format[\[LeftBracketingBar]Verbatim[Pm]\[RightBracketingBar], StandardForm] := "Pm";
Format[Cs2, StandardForm] := GTSF[C, s, 2];
Format[\[LeftBracketingBar]Verbatim[Pc]\[RightBracketingBar], StandardForm] := "Pc";
Format[Cs3, StandardForm] := GTSF[C, s, 3];
Format[\[LeftBracketingBar]Verbatim[Cm]\[RightBracketingBar], StandardForm] := "Cm";
Format[Cs4, StandardForm] := GTSF[C, s, 4];
Format[\[LeftBracketingBar]Verbatim[Cc]\[RightBracketingBar], StandardForm] := "Cc";
Format[C2h1, StandardForm] := GTSF[C, 2 h, 1];
Format[\[LeftBracketingBar]Verbatim[P2/m]\[RightBracketingBar], StandardForm] := "P2/m";
Format[C2h2, StandardForm] := GTSF[C, 2 h, 2];
Format[\[LeftBracketingBar]Verbatim[P2_ 1/m]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[2, 1] <> "/m";
Format[C2h3, StandardForm] := GTSF[C, 2 h, 3];
Format[\[LeftBracketingBar]Verbatim[C2/m]\[RightBracketingBar], StandardForm] := "C2/m";
Format[C2h4, StandardForm] := GTSF[C, 2 h, 4];
Format[\[LeftBracketingBar]Verbatim[P2/c]\[RightBracketingBar], StandardForm] := "P2/c";
Format[C2h5, StandardForm] := GTSF[C, 2 h, 5];
Format[\[LeftBracketingBar]Verbatim[P2_ 1/c]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[2, 1] <> "/c";
Format[C2h6, StandardForm] := GTSF[C, 2 h, 6];
Format[\[LeftBracketingBar]Verbatim[C2/c]\[RightBracketingBar], StandardForm] := "C2/c";

(*--- Orthorhomic ---*)
Format[D21, StandardForm] := GTSF[D, 2, 1];
Format[\[LeftBracketingBar]Verbatim[P222]\[RightBracketingBar], StandardForm] := "P222";
Format[D22, StandardForm] := GTSF[D, 2, 2];
Format[\[LeftBracketingBar]Verbatim[P222_ 1]\[RightBracketingBar], StandardForm] := "P22" <> GTSubScript[2, 1];
Format[D23, StandardForm] := GTSF[D, 2, 3];
Format[\[LeftBracketingBar]Verbatim[P2_ 12_ 12]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[2, 1] <> GTSubScript[2, 1] <> "2";
Format[D24, StandardForm] := GTSF[D, 2, 4];
Format[\[LeftBracketingBar]Verbatim[P2_ 12_ 12_ 1]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[2, 1] <> GTSubScript[2, 1] <> GTSubScript[2, 1];
Format[D25, StandardForm] := GTSF[D, 2, 5];
Format[\[LeftBracketingBar]Verbatim[C222_ 1]\[RightBracketingBar], StandardForm] := "C22" <> GTSubScript[2, 1];
Format[D26, StandardForm] := GTSF[D, 2, 6];
Format[\[LeftBracketingBar]Verbatim[C222]\[RightBracketingBar], StandardForm] := "C222";
Format[D27, StandardForm] := GTSF[D, 2, 7];
Format[\[LeftBracketingBar]Verbatim[F222]\[RightBracketingBar], StandardForm] := "F222";
Format[D28, StandardForm] := GTSF[D, 2, 8];
Format[\[LeftBracketingBar]Verbatim[I222]\[RightBracketingBar], StandardForm] := "I222";
Format[D29, StandardForm] := GTSF[D, 2, 9];
Format[\[LeftBracketingBar]Verbatim[I2_ 12_ 12_ 1]\[RightBracketingBar], StandardForm] := "I" <> GTSubScript[2, 1] <> GTSubScript[2, 1] <> GTSubScript[2, 1];
Format[C2v1, StandardForm] := GTSF[C, 2 v, 1];
Format[\[LeftBracketingBar]Verbatim[Pmm2]\[RightBracketingBar], StandardForm] := "Pmm2";
Format[C2v2, StandardForm] := GTSF[C, 2 v, 2];
Format[\[LeftBracketingBar]Verbatim[Pmc2_ 1]\[RightBracketingBar], StandardForm] := "Pmc" <> GTSubScript[2, 1];
Format[C2v3, StandardForm] := GTSF[C, 2 v, 3];
Format[\[LeftBracketingBar]Verbatim[Pcc2]\[RightBracketingBar], StandardForm] := "Pcc2";
Format[C2v4, StandardForm] := GTSF[C, 2 v, 4];
Format[\[LeftBracketingBar]Verbatim[Pma2]\[RightBracketingBar], StandardForm] := "Pma2";
Format[C2v5, StandardForm] := GTSF[C, 2 v, 5];
Format[\[LeftBracketingBar]Verbatim[Pca2_ 1]\[RightBracketingBar], StandardForm] := "Pca" <> GTSubScript[2, 1];
Format[C2v6, StandardForm] := GTSF[C, 2 v, 6];
Format[\[LeftBracketingBar]Verbatim[Pnc2]\[RightBracketingBar], StandardForm] := "Pnc2";
Format[C2v7, StandardForm] := GTSF[C, 2 v, 7];
Format[\[LeftBracketingBar]Verbatim[Pmn2_ 1]\[RightBracketingBar], StandardForm] := "Pmn" <> GTSubScript[2, 1];
Format[C2v8, StandardForm] := GTSF[C, 2 v, 8];
Format[\[LeftBracketingBar]Verbatim[Pba2]\[RightBracketingBar], StandardForm] := "Pba2";
Format[C2v9, StandardForm] := GTSF[C, 2 v, 9];
Format[\[LeftBracketingBar]Verbatim[Pna2_ 1]\[RightBracketingBar], StandardForm] := "Pna" <> GTSubScript[2, 1];
Format[C2v10, StandardForm] := GTSF[C, 2 v, 10];
Format[\[LeftBracketingBar]Verbatim[Pnn2]\[RightBracketingBar], StandardForm] := "Pnn2";
Format[C2v11, StandardForm] := GTSF[C, 2 v, 11];
Format[\[LeftBracketingBar]Verbatim[Cmm2]\[RightBracketingBar], StandardForm] := "Cmm2";
Format[C2v12, StandardForm] := GTSF[C, 2 v, 12];
Format[\[LeftBracketingBar]Verbatim[Cmc2_ 1]\[RightBracketingBar], StandardForm] := "Cmc" <> GTSubScript[2, 1];
Format[C2v13, StandardForm] := GTSF[C, 2 v, 13];
Format[\[LeftBracketingBar]Verbatim[Ccc2]\[RightBracketingBar], StandardForm] := "Ccc2";
Format[C2v14, StandardForm] := GTSF[C, 2 v, 14];
Format[\[LeftBracketingBar]Verbatim[Amm2]\[RightBracketingBar], StandardForm] := "Amm2";
Format[C2v15, StandardForm] := GTSF[C, 2 v, 15];
Format[\[LeftBracketingBar]Verbatim[Aem2]\[RightBracketingBar], StandardForm] := "Aem2";
Format[C2v16, StandardForm] := GTSF[C, 2 v, 16];
Format[\[LeftBracketingBar]Verbatim[Ama2]\[RightBracketingBar], StandardForm] := "Ama2";
Format[C2v17, StandardForm] := GTSF[C, 2 v, 17];
Format[\[LeftBracketingBar]Verbatim[Aea2]\[RightBracketingBar], StandardForm] := "Aea2";
Format[C2v18, StandardForm] := GTSF[C, 2 v, 18];
Format[\[LeftBracketingBar]Verbatim[Fmm2]\[RightBracketingBar], StandardForm] := "Fmm2";
Format[C2v19, StandardForm] := GTSF[C, 2 v, 19];
Format[\[LeftBracketingBar]Verbatim[Fdd2]\[RightBracketingBar], StandardForm] := "Fdd2";
Format[C2v20, StandardForm] := GTSF[C, 2 v, 20];
Format[\[LeftBracketingBar]Verbatim[Imm2]\[RightBracketingBar], StandardForm] := "Imm2";
Format[C2v21, StandardForm] := GTSF[C, 2 v, 21];
Format[\[LeftBracketingBar]Verbatim[Iba2]\[RightBracketingBar], StandardForm] := "Iba2";
Format[C2v22, StandardForm] := GTSF[C, 2 v, 22];
Format[\[LeftBracketingBar]Verbatim[Ima2]\[RightBracketingBar], StandardForm] := "Ima2";
Format[D2h1, StandardForm] := GTSF[D, 2 h, 1];
Format[\[LeftBracketingBar]Verbatim[Pmmm]\[RightBracketingBar], StandardForm] := "Pmmm";
Format[D2h2, StandardForm] := GTSF[D, 2 h, 2];
Format[\[LeftBracketingBar]Verbatim[Pnnn]\[RightBracketingBar], StandardForm] := "Pnnn";
Format[D2h3, StandardForm] := GTSF[D, 2 h, 3];
Format[\[LeftBracketingBar]Verbatim[Pccm]\[RightBracketingBar], StandardForm] := "Pccm";
Format[D2h4, StandardForm] := GTSF[D, 2 h, 4];
Format[\[LeftBracketingBar]Verbatim[Pban]\[RightBracketingBar], StandardForm] := "Pban";
Format[D2h5, StandardForm] := GTSF[D, 2 h, 5];
Format[\[LeftBracketingBar]Verbatim[Pmma]\[RightBracketingBar], StandardForm] := "Pmma";
Format[D2h6, StandardForm] := GTSF[D, 2 h, 6];
Format[\[LeftBracketingBar]Verbatim[Pnna]\[RightBracketingBar], StandardForm] := "Pnna";
Format[D2h7, StandardForm] := GTSF[D, 2 h, 7];
Format[\[LeftBracketingBar]Verbatim[Pmna]\[RightBracketingBar], StandardForm] := "Pmna";
Format[D2h8, StandardForm] := GTSF[D, 2 h, 8];
Format[\[LeftBracketingBar]Verbatim[Pcca]\[RightBracketingBar], StandardForm] := "Pcca";
Format[D2h9, StandardForm] := GTSF[D, 2 h, 9];
Format[\[LeftBracketingBar]Verbatim[Pbam]\[RightBracketingBar], StandardForm] := "Pbam";
Format[D2h10, StandardForm] := GTSF[D, 2 h, 10];
Format[\[LeftBracketingBar]Verbatim[Pccn]\[RightBracketingBar], StandardForm] := "Pccn";
Format[D2h11, StandardForm] := GTSF[D, 2 h, 11];
Format[\[LeftBracketingBar]Verbatim[Pbcm]\[RightBracketingBar], StandardForm] := "Pbcm";
Format[D2h12, StandardForm] := GTSF[D, 2 h, 12];
Format[\[LeftBracketingBar]Verbatim[Pnnm]\[RightBracketingBar], StandardForm] := "Pnnm";
Format[D2h13, StandardForm] := GTSF[D, 2 h, 13];
Format[\[LeftBracketingBar]Verbatim[Pmmn]\[RightBracketingBar], StandardForm] := "Pmmn";
Format[D2h14, StandardForm] := GTSF[D, 2 h, 14];
Format[\[LeftBracketingBar]Verbatim[Pbcn]\[RightBracketingBar], StandardForm] := "Pbcn";
Format[D2h15, StandardForm] := GTSF[D, 2 h, 15];
Format[\[LeftBracketingBar]Verbatim[Pbca]\[RightBracketingBar], StandardForm] := "Pbca";
Format[D2h16, StandardForm] := GTSF[D, 2 h, 16];
Format[\[LeftBracketingBar]Verbatim[Pnma]\[RightBracketingBar], StandardForm] := "Pnma";
Format[D2h17, StandardForm] := GTSF[D, 2 h, 17];
Format[\[LeftBracketingBar]Verbatim[Cmcm]\[RightBracketingBar], StandardForm] := "Cmcm";
Format[D2h18, StandardForm] := GTSF[D, 2 h, 18];
Format[\[LeftBracketingBar]Verbatim[Cmce]\[RightBracketingBar], StandardForm] := "Cmce";
Format[D2h19, StandardForm] := GTSF[D, 2 h, 19];
Format[\[LeftBracketingBar]Verbatim[Cmmm]\[RightBracketingBar], StandardForm] := "Cmmm";
Format[D2h20, StandardForm] := GTSF[D, 2 h, 20];
Format[\[LeftBracketingBar]Verbatim[Cccm]\[RightBracketingBar], StandardForm] := "Cccm";
Format[D2h21, StandardForm] := GTSF[D, 2 h, 21];
Format[\[LeftBracketingBar]Verbatim[Cmme]\[RightBracketingBar], StandardForm] := "Cmme";
Format[D2h22, StandardForm] := GTSF[D, 2 h, 22];
Format[\[LeftBracketingBar]Verbatim[Ccce]\[RightBracketingBar], StandardForm] := "Ccce";
Format[D2h23, StandardForm] := GTSF[D, 2 h, 23];
Format[\[LeftBracketingBar]Verbatim[Fmmm]\[RightBracketingBar], StandardForm] := "Fmmm";
Format[D2h24, StandardForm] := GTSF[D, 2 h, 24];
Format[\[LeftBracketingBar]Verbatim[Fddd]\[RightBracketingBar], StandardForm] := "Fddd";
Format[D2h25, StandardForm] := GTSF[D, 2 h, 25];
Format[\[LeftBracketingBar]Verbatim[Immm]\[RightBracketingBar], StandardForm] := "Immm";
Format[D2h26, StandardForm] := GTSF[D, 2 h, 26];
Format[\[LeftBracketingBar]Verbatim[Ibam]\[RightBracketingBar], StandardForm] := "Ibam";
Format[D2h27, StandardForm] := GTSF[D, 2 h, 27];
Format[\[LeftBracketingBar]Verbatim[Ibca]\[RightBracketingBar], StandardForm] := "Ibca";
Format[D2h28, StandardForm] := GTSF[D, 2 h, 28];
Format[\[LeftBracketingBar]Verbatim[Imma]\[RightBracketingBar], StandardForm] := "Imma";

(*--- Tetragonal ---*)
Format[C41, StandardForm] := GTSF[C, 4, 1];
Format[\[LeftBracketingBar]Verbatim[P4]\[RightBracketingBar], StandardForm] := "P4";
Format[C42, StandardForm] := GTSF[C, 4, 2];
Format[\[LeftBracketingBar]Verbatim[P4_ 1]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 1];
Format[C43, StandardForm] := GTSF[C, 4, 3];
Format[\[LeftBracketingBar]Verbatim[P4_ 2]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2];
Format[C44, StandardForm] := GTSF[C, 4, 4];
Format[\[LeftBracketingBar]Verbatim[P4_ 3]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 3];
Format[C45, StandardForm] := GTSF[C, 4, 5]
Format[\[LeftBracketingBar]Verbatim[I4]\[RightBracketingBar], StandardForm] := "I4";
Format[C46, StandardForm] := GTSF[C, 4, 6];
Format[\[LeftBracketingBar]Verbatim[I4_ 1]\[RightBracketingBar], StandardForm] := "I" <> GTSubScript[4, 1];
Format[S41, StandardForm] := GTSF[S, 4, 1];
Format[\[LeftBracketingBar]Verbatim[P - 4]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4];
Format[S42, StandardForm] := GTSF[S, 4, 2];
Format[\[LeftBracketingBar]Verbatim[I - 4]\[RightBracketingBar], StandardForm] := "I" <> GTOverScript[4];
Format[C4h1, StandardForm] := GTSF[C, 4 h, 1];
Format[\[LeftBracketingBar]Verbatim[P4/m]\[RightBracketingBar], StandardForm] := "P4/m";
Format[C4h2, StandardForm] := GTSF[C, 4 h, 2];
Format[\[LeftBracketingBar]Verbatim[P4_ 2/m]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "/m";
Format[C4h3, StandardForm] := GTSF[C, 4 h, 3];
Format[\[LeftBracketingBar]Verbatim[P4/n]\[RightBracketingBar], StandardForm] := "P4/n";
Format[C4h4, StandardForm] := GTSF[C, 4 h, 4];
Format[\[LeftBracketingBar]Verbatim[P4_ 2/n]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "/n";
Format[C4h5, StandardForm] := GTSF[C, 4 h, 5];
Format[\[LeftBracketingBar]Verbatim[I4/m]\[RightBracketingBar], StandardForm] := "I4/m";
Format[C4h6, StandardForm] := GTSF[C, 4 h, 6];
Format[\[LeftBracketingBar]Verbatim[I4_ 1/a]\[RightBracketingBar], StandardForm] := "I" <> GTSubScript[4, 1] <> "/a";
Format[D41, StandardForm] := GTSF[D, 4, 1];
Format[\[LeftBracketingBar]Verbatim[P422]\[RightBracketingBar], StandardForm] := "P422";
Format[D42, StandardForm] := GTSF[D, 4, 2];
Format[\[LeftBracketingBar]Verbatim[P42_ 12]\[RightBracketingBar], StandardForm] := "P4" <> GTSubScript[2, 1] <> "2";
Format[D43, StandardForm] := GTSF[D, 4, 3];
Format[\[LeftBracketingBar]Verbatim[P4_ 122]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 1] <> "22";
Format[D44, StandardForm] := GTSF[D, 4, 4];
Format[\[LeftBracketingBar]Verbatim[P4_ 12_ 12]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 1] <> GTSubScript[2, 1] <> "2";
Format[D45, StandardForm] := GTSF[D, 4, 5];
Format[\[LeftBracketingBar]Verbatim[P4_ 222]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "22";
Format[D46, StandardForm] := GTSF[D, 4, 6];
Format[\[LeftBracketingBar]Verbatim[P4_ 22_ 12]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> GTSubScript[2, 1] <> "2";
Format[D47, StandardForm] := GTSF[D, 4, 7];
Format[\[LeftBracketingBar]Verbatim[P4_ 322]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 3] <> "22";
Format[D48, StandardForm] := GTSF[D, 4, 8];
Format[\[LeftBracketingBar]Verbatim[P4_ 32_ 12]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 3] <> GTSubScript[2, 1] <> "2";
Format[D49, StandardForm] := GTSF[D, 4, 9];
Format[\[LeftBracketingBar]Verbatim[I422]\[RightBracketingBar], StandardForm] := "I422";
Format[D410, StandardForm] := GTSF[D, 4, 10];
Format[\[LeftBracketingBar]Verbatim[I4_ 122]\[RightBracketingBar], StandardForm] := "I" <> GTSubScript[4, 1] <> "22";
Format[C4v1, StandardForm] := GTSF[C, 4 v, 1];
Format[\[LeftBracketingBar]Verbatim[P4mm]\[RightBracketingBar], StandardForm] := "P4mm";
Format[C4v2, StandardForm] := GTSF[C, 4 v, 2];
Format[\[LeftBracketingBar]Verbatim[P4bm]\[RightBracketingBar], StandardForm] := "P4bm";
Format[C4v3, StandardForm] := GTSF[C, 4 v, 3];
Format[\[LeftBracketingBar]Verbatim[P4_ 2 cm]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "cm";
Format[C4v4, StandardForm] := GTSF[C, 4 v, 4];
Format[\[LeftBracketingBar]Verbatim[P4_ 2 nm]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "nm";
Format[C4v5, StandardForm] := GTSF[C, 4 v, 5];
Format[\[LeftBracketingBar]Verbatim[P4cc]\[RightBracketingBar], StandardForm] := "P4cc";
Format[C4v6, StandardForm] := GTSF[C, 4 v, 6];
Format[\[LeftBracketingBar]Verbatim[P4nc]\[RightBracketingBar], StandardForm] := "P4nc";
Format[C4v7, StandardForm] := GTSF[C, 4 v, 7];
Format[\[LeftBracketingBar]Verbatim[P4_ 2 mc]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "mc";
Format[C4v8, StandardForm] := GTSF[C, 4 v, 8];
Format[\[LeftBracketingBar]Verbatim[P4_ 2 bc]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "bc";
Format[C4v9, StandardForm] := GTSF[C, 4 v, 9];
Format[\[LeftBracketingBar]Verbatim[I4mm]\[RightBracketingBar], StandardForm] := "I4mm";
Format[C4v10, StandardForm] := GTSF[C, 4 v, 10];
Format[\[LeftBracketingBar]Verbatim[I4cm]\[RightBracketingBar], StandardForm] := "I4cm";
Format[C4v11, StandardForm] := GTSF[C, 4 v, 11];
Format[\[LeftBracketingBar]Verbatim[I4_ 1 md]\[RightBracketingBar], StandardForm] := "I" <> GTSubScript[4, 1] <> "md";
Format[C4v12, StandardForm] := GTSF[C, 4 v, 12];
Format[\[LeftBracketingBar]Verbatim[I4_ 1 cd]\[RightBracketingBar], StandardForm] := "I" <> GTSubScript[4, 1] <> "cd";
Format[D2d1, StandardForm] := GTSF[D, 2 d, 1];
Format[\[LeftBracketingBar]Verbatim[P - 42 m]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4] <> "2m";
Format[D2d2, StandardForm] := GTSF[D, 2 d, 2];
Format[\[LeftBracketingBar]Verbatim[P - 42 c]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4] <> "2c";
Format[D2d3, StandardForm] := GTSF[D, 2 d, 3];
Format[\[LeftBracketingBar]Verbatim[P - 42_ 1 m]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4] <> GTSubScript[2, 1] <> "m";
Format[D2d4, StandardForm] := GTSF[D, 2 d, 4];
Format[\[LeftBracketingBar]Verbatim[P - 42_ 1 c]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4] <> GTSubScript[2, 1] <> "c";
Format[D2d5, StandardForm] := GTSF[D, 2 d, 5];
Format[\[LeftBracketingBar]Verbatim[P - 4 m2]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4] <> "m2";
Format[D2d6, StandardForm] := GTSF[D, 2 d, 6];
Format[\[LeftBracketingBar]Verbatim[P - 4 c2]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4] <> "c2";
Format[D2d7, StandardForm] := GTSF[D, 2 d, 7];
Format[\[LeftBracketingBar]Verbatim[P - 4 b2]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4] <> "b2";
Format[D2d8, StandardForm] := GTSF[D, 2 d, 8];
Format[\[LeftBracketingBar]Verbatim[P - 4 n2]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4] <> "n2";
Format[D2d9, StandardForm] := GTSF[D, 2 d, 9];
Format[\[LeftBracketingBar]Verbatim[I - 4 m2]\[RightBracketingBar], StandardForm] := "I" <> GTOverScript[4] <> "m2";
Format[D2d10, StandardForm] := GTSF[D, 2 d, 10];
Format[\[LeftBracketingBar]Verbatim[I - 4 c2]\[RightBracketingBar], StandardForm] := "I" <> GTOverScript[4] <> "c2";
Format[D2d11, StandardForm] := GTSF[D, 2 d, 11];
Format[\[LeftBracketingBar]Verbatim[I - 42 m]\[RightBracketingBar], StandardForm] := "I" <> GTOverScript[4] <> "2m";
Format[D2d12, StandardForm] := GTSF[D, 2 d, 12];
Format[\[LeftBracketingBar]Verbatim[I - 42 d]\[RightBracketingBar], StandardForm] := "I" <> GTOverScript[4] <> "2d";
Format[D4h1, StandardForm] := GTSF[D, 4 h, 1];
Format[\[LeftBracketingBar]Verbatim[P4/mmm]\[RightBracketingBar], StandardForm] := "P4/mmm";
Format[D4h2, StandardForm] := GTSF[D, 4 h, 2];
Format[\[LeftBracketingBar]Verbatim[P4/mcc]\[RightBracketingBar], StandardForm] := "P4/mcc";
Format[D4h3, StandardForm] := GTSF[D, 4 h, 3];
Format[\[LeftBracketingBar]Verbatim[P4/nbm]\[RightBracketingBar], StandardForm] := "P4/nbm";
Format[D4h4, StandardForm] := GTSF[D, 4 h, 4];
Format[\[LeftBracketingBar]Verbatim[P4/nnc]\[RightBracketingBar], StandardForm] := "P4/nnc";
Format[D4h5, StandardForm] := GTSF[D, 4 h, 5];
Format[\[LeftBracketingBar]Verbatim[P4/mbm]\[RightBracketingBar], StandardForm] := "P4/mbm";
Format[D4h6, StandardForm] := GTSF[D, 4 h, 6];
Format[\[LeftBracketingBar]Verbatim[P4/mnc]\[RightBracketingBar], StandardForm] := "P4/mnc";
Format[D4h7, StandardForm] := GTSF[D, 4 h, 7];
Format[\[LeftBracketingBar]Verbatim[P4/nmm]\[RightBracketingBar], StandardForm] := "P4/nmm";
Format[D4h8, StandardForm] := GTSF[D, 4 h, 8];
Format[\[LeftBracketingBar]Verbatim[P4/ncc]\[RightBracketingBar], StandardForm] := "P4/ncc";
Format[D4h9, StandardForm] := GTSF[D, 4 h, 9];
Format[\[LeftBracketingBar]Verbatim[P4_ 2/mmc]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "/mmc";
Format[D4h10, StandardForm] := GTSF[D, 4 h, 10];
Format[\[LeftBracketingBar]Verbatim[P4_ 2/mcm]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "/mcm";
Format[D4h11, StandardForm] := GTSF[D, 4 h, 11];
Format[\[LeftBracketingBar]Verbatim[P4_ 2/nbc]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "/nbc";
Format[D4h12, StandardForm] := GTSF[D, 4 h, 12];
Format[\[LeftBracketingBar]Verbatim[P4_ 2/nnm]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "/nnm";
Format[D4h13, StandardForm] := GTSF[D, 4 h, 13];
Format[\[LeftBracketingBar]Verbatim[P4_ 2/mbc]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "/mbc";
Format[D4h14, StandardForm] := GTSF[D, 4 h, 14];
Format[\[LeftBracketingBar]Verbatim[P4_ 2/mnm]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "/mnm";
Format[D4h15, StandardForm] := GTSF[D, 4 h, 15];
Format[\[LeftBracketingBar]Verbatim[P4_ 2/nmc]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "/nmc";
Format[D4h16, StandardForm] := GTSF[D, 4 h, 16];
Format[\[LeftBracketingBar]Verbatim[P4_ 2/ncm]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "/ncm";
Format[D4h17, StandardForm] := GTSF[D, 4 h, 17];
Format[\[LeftBracketingBar]Verbatim[I4/mmmm]\[RightBracketingBar], StandardForm] := "I4/mmm";
Format[D4h18, StandardForm] := GTSF[D, 4 h, 18];
Format[\[LeftBracketingBar]Verbatim[I4/mcm]\[RightBracketingBar], StandardForm] := "I4/mcm";
Format[D4h19, StandardForm] := GTSF[D, 4 h, 19];
Format[\[LeftBracketingBar]Verbatim[I4_ 1/amd]\[RightBracketingBar], StandardForm] := "I" <> GTSubScript[4, 1] <> "/amd";
Format[D4h20, StandardForm] := GTSF[D, 4 h, 20];
Format[\[LeftBracketingBar]Verbatim[I4_ 1/acd]\[RightBracketingBar], StandardForm] := "I" <> GTSubScript[4, 1] <> "/acd";

(*--- Trigonal ---*)
Format[C31, StandardForm] := GTSF[C, 3, 1];
Format[\[LeftBracketingBar]Verbatim[P3]\[RightBracketingBar], StandardForm] := "P3";
Format[C32, StandardForm] := GTSF[C, 3, 2];
Format[\[LeftBracketingBar]Verbatim[P3_ 1]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[3, 1];
Format[C33, StandardForm] := GTSF[C, 3, 3];
Format[\[LeftBracketingBar]Verbatim[P3_ 2]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[3, 2];
Format[C34, StandardForm] := GTSF[C, 3, 4];
Format[\[LeftBracketingBar]Verbatim[R3]\[RightBracketingBar], StandardForm] := "R3";
Format[C3i1, StandardForm] := GTSF[C, 3 i, 1];
Format[\[LeftBracketingBar]Verbatim[P - 3]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[3];
Format[C3i2, StandardForm] := GTSF[C, 3 i, 2];
Format[\[LeftBracketingBar]Verbatim[R - 3]\[RightBracketingBar], StandardForm] := "R" <> GTOverScript[3];
Format[D31, StandardForm] := GTSF[D, 3, 1];
Format[\[LeftBracketingBar]Verbatim[P312]\[RightBracketingBar], StandardForm] := "P312";
Format[D32, StandardForm] := GTSF[D, 3, 2];
Format[\[LeftBracketingBar]Verbatim[P321]\[RightBracketingBar], StandardForm] := "P321";
Format[D33, StandardForm] := GTSF[D, 3, 3];
Format[\[LeftBracketingBar]Verbatim[P3_ 112]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[3, 1] <> "12";
Format[D34, StandardForm] := GTSF[D, 3, 4];
Format[\[LeftBracketingBar]Verbatim[P3_ 121]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[3, 1] <> "21";
Format[D35, StandardForm] := GTSF[D, 3, 5];
Format[\[LeftBracketingBar]Verbatim[P3_ 212]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[3, 2] <> "12";
Format[D36, StandardForm] := GTSF[D, 3, 6];
Format[\[LeftBracketingBar]Verbatim[P3_ 221]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[3, 2] <> "21";
Format[D37, StandardForm] := GTSF[D, 3, 7];
Format[\[LeftBracketingBar]Verbatim[R32]\[RightBracketingBar], StandardForm] := "R32";
Format[C3v1, StandardForm] := GTSF[C, 3 v, 1];
Format[\[LeftBracketingBar]Verbatim[P3m1]\[RightBracketingBar], StandardForm] := "P3m1";
Format[C3v2, StandardForm] := GTSF[C, 3 v, 2];
Format[\[LeftBracketingBar]Verbatim[P31m]\[RightBracketingBar], StandardForm] := "P31m";
Format[C3v3, StandardForm] := GTSF[C, 3 v, 3];
Format[\[LeftBracketingBar]Verbatim[P3c1]\[RightBracketingBar], StandardForm] := "P3c1";
Format[C3v4, StandardForm] := GTSF[C, 3 v, 4];
Format[\[LeftBracketingBar]Verbatim[P31c]\[RightBracketingBar], StandardForm] := "P31c";
Format[C3v5, StandardForm] := GTSF[C, 3 v, 5];
Format[\[LeftBracketingBar]Verbatim[R3m]\[RightBracketingBar], StandardForm] := "R3m";
Format[C3v6, StandardForm] := GTSF[C, 3 v, 6];
Format[\[LeftBracketingBar]Verbatim[R3c]\[RightBracketingBar], StandardForm] := "R3c";
Format[D3d1, StandardForm] := GTSF[D, 3 d, 1];
Format[\[LeftBracketingBar]Verbatim[P - 31 m]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[3] <> "1m";
Format[D3d2, StandardForm] := GTSF[D, 3 d, 2];
Format[\[LeftBracketingBar]Verbatim[P - 31 c]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[3] <> "1c";
Format[D3d3, StandardForm] := GTSF[D, 3 d, 3];
Format[\[LeftBracketingBar]Verbatim[P - 3 m1]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[3] <> "m1";
Format[D3d4, StandardForm] := GTSF[D, 3 d, 4];
Format[\[LeftBracketingBar]Verbatim[P - 3 c1]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[3] <> "c1";
Format[D3d5, StandardForm] := GTSF[D, 3 d, 5];
Format[\[LeftBracketingBar]Verbatim[R - 3 m]\[RightBracketingBar], StandardForm] := "R" <> GTOverScript[3] <> "m";
Format[D3d6, StandardForm] := GTSF[D, 3 d, 6];
Format[\[LeftBracketingBar]Verbatim[R - 3 c]\[RightBracketingBar], StandardForm] := "R" <> GTOverScript[3] <> "c";

(*--- Hexagonal ---*)
Format[C61, StandardForm] := GTSF[C, 6, 1];
Format[\[LeftBracketingBar]Verbatim[P6]\[RightBracketingBar], StandardForm] := "P6";
Format[C62, StandardForm] := GTSF[C, 6, 2];
Format[\[LeftBracketingBar]Verbatim[P6_ 1]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 1];
Format[C63, StandardForm] := GTSF[C, 6, 3];
Format[\[LeftBracketingBar]Verbatim[P6_ 5]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 5];
Format[C64, StandardForm] := GTSF[C, 6, 4];
Format[\[LeftBracketingBar]Verbatim[P6_ 2]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 2];
Format[C65, StandardForm] := GTSF[C, 6, 5];
Format[\[LeftBracketingBar]Verbatim[P6_ 4]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 4];
Format[C66, StandardForm] := GTSF[C, 6, 6];
Format[\[LeftBracketingBar]Verbatim[P6_ 3]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 3];
Format[C3h1, StandardForm] := GTSF[C, 3 h, 1];
Format[\[LeftBracketingBar]Verbatim[P - 6]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[6];
Format[C6h1, StandardForm] := GTSF[C, 6 h, 1];
Format[\[LeftBracketingBar]Verbatim[P6/m]\[RightBracketingBar], StandardForm] := "P6/m";
Format[C6h2, StandardForm] := GTSF[C, 6 h, 2];
Format[\[LeftBracketingBar]Verbatim[P6_ 3/m]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 3] <> "/m";
Format[D61, StandardForm] := GTSF[D, 6, 1];
Format[\[LeftBracketingBar]Verbatim[P622]\[RightBracketingBar], StandardForm] := "P622";
Format[D62, StandardForm] := GTSF[D, 6, 2];
Format[\[LeftBracketingBar]Verbatim[P6_ 122]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 1] <> "22";
Format[D63, StandardForm] := GTSF[D, 6, 3];
Format[\[LeftBracketingBar]Verbatim[P6_ 522]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 5] <> "22";
Format[D64, StandardForm] := GTSF[D, 6, 4];
Format[\[LeftBracketingBar]Verbatim[P6_ 222]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 2] <> "22";
Format[D65, StandardForm] := GTSF[D, 6, 5];
Format[\[LeftBracketingBar]Verbatim[P6_ 422]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 4] <> "22";
Format[D66, StandardForm] := GTSF[D, 6, 6];
Format[\[LeftBracketingBar]Verbatim[P6_ 322]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 3] <> "22";
Format[C6v1, StandardForm] := GTSF[C, 6 V, 1];
Format[\[LeftBracketingBar]Verbatim[P6mm]\[RightBracketingBar], StandardForm] := "P6mm";
Format[C6v2, StandardForm] := GTSF[C, 6 V, 2];
Format[\[LeftBracketingBar]Verbatim[P6cc]\[RightBracketingBar], StandardForm] := "P6cc";
Format[C6v3, StandardForm] := GTSF[C, 6 v, 3];
Format[\[LeftBracketingBar]Verbatim[P6_ 3 cm]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 3] <> "cm";
Format[C6v4, StandardForm] := GTSF[C, 6 v, 4];
Format[\[LeftBracketingBar]Verbatim[P6_ 3 mc]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 3] <> "mc";
Format[D3h1, StandardForm] := GTSF[D, 3 h, 1];
Format[\[LeftBracketingBar]Verbatim[P - 6 m2]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[6] <> "m2";
Format[D3h2, StandardForm] := GTSF[D, 3 h, 2];
Format[\[LeftBracketingBar]Verbatim[P - 6 c2]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[6] <> "c2";
Format[D3h3, StandardForm] := GTSF[D, 3 h, 3];
Format[\[LeftBracketingBar]Verbatim[P - 62 m]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[6] <> "2m";
Format[D3h4, StandardForm] := GTSF[D, 3 h, 4];
Format[\[LeftBracketingBar]Verbatim[P - 62 c]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[6] <> "2c";
Format[D6h1, StandardForm] := GTSF[D, 6 h, 1];
Format[\[LeftBracketingBar]Verbatim[P6/mmm]\[RightBracketingBar], StandardForm] := "P6/mmm";
Format[D6h2, StandardForm] := GTSF[D, 6 h, 2];
Format[\[LeftBracketingBar]Verbatim[P6/mcc]\[RightBracketingBar], StandardForm] := "P6/mcc";
Format[D6h3, StandardForm] := GTSF[D, 6 h, 3];
Format[\[LeftBracketingBar]Verbatim[P6_ 3/mcm]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 3] <> "/mcm";
Format[D6h4, StandardForm] := GTSF[D, 6 h, 4];
Format[\[LeftBracketingBar]Verbatim[P6_ 3/mmc]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[6, 3] <> "/mmc";

(*--- Cubic ---*)
Format[T1, StandardForm] := GTSF[T, "", 1];
Format[\[LeftBracketingBar]Verbatim[P23]\[RightBracketingBar], StandardForm] := "P23";
Format[T2, StandardForm] := GTSF[T, "", 2];
Format[\[LeftBracketingBar]Verbatim[F23]\[RightBracketingBar], StandardForm] := "F23";
Format[T3, StandardForm] := GTSF[T, "", 3];
Format[\[LeftBracketingBar]Verbatim[I23]\[RightBracketingBar], StandardForm] := "I23";
Format[T4, StandardForm] := GTSF[T, "", 4];
Format[\[LeftBracketingBar]Verbatim[P2_ 13]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[2, 1] <> "3";
Format[T5, StandardForm] := GTSF[T, "", 5];
Format[\[LeftBracketingBar]Verbatim[I2_ 13]\[RightBracketingBar], StandardForm] := "I" <> GTSubScript[2, 1] <> "3";
Format[Th1, StandardForm] := GTSF[T, h, 1];
Format[\[LeftBracketingBar]Verbatim[Pm - 3]\[RightBracketingBar], StandardForm] := "Pm" <> GTOverScript[3];
Format[Th2, StandardForm] := GTSF[T, h, 2];
Format[\[LeftBracketingBar]Verbatim[Pn - 3]\[RightBracketingBar], StandardForm] := "Pn" <> GTOverScript[3];
Format[Th3, StandardForm] := GTSF[T, h, 1];
Format[\[LeftBracketingBar]Verbatim[Fm - 3]\[RightBracketingBar], StandardForm] := "Fm" <> GTOverScript[3];
Format[Th4, StandardForm] := GTSF[T, h, 4];
Format[\[LeftBracketingBar]Verbatim[Fd - 3]\[RightBracketingBar], StandardForm] := "Fd" <> GTOverScript[3];
Format[Th5, StandardForm] := GTSF[T, h, 5];
Format[\[LeftBracketingBar]Verbatim[Im - 3]\[RightBracketingBar], StandardForm] := "Im" <> GTOverScript[3];
Format[Th6, StandardForm] := GTSF[T, h, 6];
Format[\[LeftBracketingBar]Verbatim[Pa - 3]\[RightBracketingBar], StandardForm] := "Pa" <> GTOverScript[3];
Format[Th7, StandardForm] := GTSF[T, h, 7];
Format[\[LeftBracketingBar]Verbatim[Ia - 3]\[RightBracketingBar], StandardForm] := "Ia" <> GTOverScript[3];
Format[O1, StandardForm] := GTSF[O, "", 1];
Format[\[LeftBracketingBar]Verbatim[P432]\[RightBracketingBar], StandardForm] := "P432";
Format[O2, StandardForm] := GTSF[O, "", 2];
Format[\[LeftBracketingBar]Verbatim[P4_ 232]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 2] <> "32";
Format[O3, StandardForm] := GTSF[O, "", 3];
Format[\[LeftBracketingBar]Verbatim[F432]\[RightBracketingBar], StandardForm] := "F432";
Format[O4, StandardForm] := GTSF[O, "", 4];
Format[\[LeftBracketingBar]Verbatim[F4_ 132]\[RightBracketingBar], StandardForm] := "F" <> GTSubScript[4, 1] <> "32";
Format[O5, StandardForm] := GTSF[O, "", 5];
Format[\[LeftBracketingBar]Verbatim[I432]\[RightBracketingBar], StandardForm] := "I432";
Format[O6, StandardForm] := GTSF[O, "", 6];
Format[\[LeftBracketingBar]Verbatim[P4_ 332]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 3] <> "32";
Format[O7, StandardForm] := GTSF[O, "", 7];
Format[\[LeftBracketingBar]Verbatim[P4_ 132]\[RightBracketingBar], StandardForm] := "P" <> GTSubScript[4, 1] <> "32";
Format[O8, StandardForm] := GTSF[O, "", 8];
Format[\[LeftBracketingBar]Verbatim[I4_ 132]\[RightBracketingBar], StandardForm] := "I" <> GTSubScript[4, 1] <> "32";
Format[Td1, StandardForm] := GTSF[T, d, 1];
Format[\[LeftBracketingBar]Verbatim[P - 43 m]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4] <> "3m";
Format[Td2, StandardForm] := GTSF[T, d, 2];
Format[\[LeftBracketingBar]Verbatim[F - 43 m]\[RightBracketingBar], StandardForm] := "F" <> GTOverScript[4] <> "3m";
Format[Td3, StandardForm] := GTSF[T, d, 3];
Format[\[LeftBracketingBar]Verbatim[I - 43 m]\[RightBracketingBar], StandardForm] := "I" <> GTOverScript[4] <> "3m";
Format[Td4, StandardForm] := GTSF[T, d, 4];
Format[\[LeftBracketingBar]Verbatim[P - 43 n]\[RightBracketingBar], StandardForm] := "P" <> GTOverScript[4] <> "3n";
Format[Td5, StandardForm] := GTSF[T, d, 5];
Format[\[LeftBracketingBar]Verbatim[F - 43 c]\[RightBracketingBar], StandardForm] := "F" <> GTOverScript[4] <> "3c";
Format[Td6, StandardForm] := GTSF[T, d, 6];
Format[\[LeftBracketingBar]Verbatim[I - 43 d]\[RightBracketingBar], StandardForm] := "I" <> GTOverScript[4] <> "3d";
Format[Oh1, StandardForm] := GTSF[O, h, 1];
Format[\[LeftBracketingBar]Verbatim[Pm - 3 m]\[RightBracketingBar], StandardForm] := "Pm" <> GTOverScript[3] <> "m";
Format[Oh2, StandardForm] := GTSF[O, h, 2];
Format[\[LeftBracketingBar]Verbatim[Pn - 3 n]\[RightBracketingBar], StandardForm] := "Pn" <> GTOverScript[3] <> "n";
Format[Oh3, StandardForm] := GTSF[O, h, 3];
Format[\[LeftBracketingBar]Verbatim[Pm - 3 n]\[RightBracketingBar], StandardForm] := "Pm" <> GTOverScript[3] <> "n";
Format[Oh4, StandardForm] := GTSF[O, h, 4];
Format[\[LeftBracketingBar]Verbatim[Pn - 3 m]\[RightBracketingBar], StandardForm] := "Pn" <> GTOverScript[3] <> "m";
Format[Oh5, StandardForm] := GTSF[O, h, 5];
Format[\[LeftBracketingBar]Verbatim[Fm - 3 m]\[RightBracketingBar], StandardForm] := "Fm" <> GTOverScript[3] <> "m";
Format[Oh6, StandardForm] := GTSF[O, h, 6];
Format[\[LeftBracketingBar]Verbatim[Fm - 3 c]\[RightBracketingBar], StandardForm] := "Fm" <> GTOverScript[3] <> "c";
Format[Oh7, StandardForm] := GTSF[O, h, 7];
Format[\[LeftBracketingBar]Verbatim[Fd - 3 m]\[RightBracketingBar], StandardForm] := "Fd" <> GTOverScript[3] <> "m";
Format[Oh8, StandardForm] := GTSF[O, h, 8];
Format[\[LeftBracketingBar]Verbatim[Fd - 3 c]\[RightBracketingBar], StandardForm] := "Fd" <> GTOverScript[3] <> "c";
Format[Oh9, StandardForm] := GTSF[O, h, 9];
Format[\[LeftBracketingBar]Verbatim[Im - 3 m]\[RightBracketingBar], StandardForm] := "Im" <> GTOverScript[3] <> "m";
Format[Oh10, StandardForm] := GTSF[O, h, 10];
Format[\[LeftBracketingBar]Verbatim[Ia - 3 d]\[RightBracketingBar], StandardForm] := "Ia" <> GTOverScript[3] <> "d";
(*--- end spacegroup notations ---*)


(*--- begin space group elements *)
Unprotect[gtspacegroupelements]

gtspacegroupelements = {{1, {}}, {2, {}}, {3, {}}, {4, {}}, {5, {}}, {6, {}}, {7, {}}, {8, \
{}}, {9, {}}, {10, {}}, {11, {}}, {12, {}}, {13, {}}, 

{14, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket]}}, 

{15, {}}, {16, {}}, {17, {}}, {18, {}},

{19, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2,0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket]}}, 
  	
{20, {}}, {21, {}}, {22, \
{}}, {23, {}}, {24, {}}, {25, {}}, {26, {}}, {27, {}}, {28, {}}, {29, \
{}}, {30, {}}, {31, {}}, {32, {}}, {33, {}}, {34, {}}, {35, {}}, {36, \
{}}, {37, {}}, {38, {}}, {39, {}}, {40, {}}, {41, {}}, {42, {}}, {43, \
{}}, {44, {}}, {45, {}}, {46, {}}, {47, {}}, {48, {}}, {49, {}}, {50, \
{}}, {51, {}}, {52, {}}, {53, {}}, {54, {}}, {55, {}}, {56, {}}, {57, \
{}}, {58, {}}, {59, {}}, {60, {}}, {61, {}}, {62, {}}, {63, {}}, {64, \
{}}, {65, {}}, {66, {}}, {67, {}}, {68, {}}, {69, {}}, {70, {}}, {71, \
{}}, {72, {}}, {73, {}}, {74, {}}, {75, {}}, {76, {}}, {77, {}}, {78, \
{}}, {79, {}}, {80, {}}, {81, {}}, {82, {}}, {83, {}}, {84, {}}, {85, \
{}}, {86, {}}, {87, {}}, 

{88, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] IEe, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {3/4, 1/4, 1/2} \[RightAngleBracket], 
 \[LeftAngleBracket] IC2z, {3/4, 1/4, 1/2} \[RightAngleBracket]}}, 
	
{89, {}}, {90, {}}, {91, {}}, {92, \
{}}, {93, {}}, {94, {}}, {95, {}}, {96, {}}, {97, {}}, {98, {}}, {99, \
{}}, {100, {}}, {101, {}}, {102, {}}, {103, {}}, {104, {}}, {105, \
{}}, {106, {}}, {107, {}}, {108, {}}, {109, {}}, {110, {}}, {111, \
{}}, {112, {}}, {113, {}}, {114, {}}, {115, {}}, {116, {}}, {117, \
{}}, {118, {}}, {119, {}}, {120, {}}, {121, {}}, {122, {}}, {123, \
{}}, {124, {}}, {125, {}}, {126, {}}, {127, {}}, {128, {}}, {129, \
{}}, {130, {}}, {131, {}}, {132, {}}, {133, {}}, {134, {}}, {135, \
{}}, {136, {}}, {137, {}}, {138, {}}, {139, {}}, {140, {}}, {141, \
{}}, {142, {}}, {143, {}}, {144, {}}, {145, {}}, {146, {}}, {147, \
{}}, {148, {}}, {149, {}}, {150, {}}, {151, {}}, {152, {}}, {153, \
{}}, {154, {}}, {155, {}}, {156, {}}, {157, {}}, {158, {}}, {159, \
{}}, {160, {}}, {161, {}}, {162, {}}, {163, {}}, {164, {}}, {165, \
{}}, {166, {}}, {167, {}}, {168, {}}, {169, {}}, {170, {}}, {171, \
{}}, {172, {}}, {173, {}}, {174, {}}, {175, {}}, {176, {}}, {177, \
{}}, {178, {}}, {179, {}}, {180, {}}, {181, {}}, {182, {}}, {183, \
{}}, {184, {}}, {185, {}}, {186, {}}, {187, {}}, {188, {}}, {189, \
{}}, {190, {}}, {191, {}}, {192, {}}, {193, {}}, 

{194, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 1/2} \[RightAngleBracket]}}, 

 {195,{}}, {196, {}}, {197, {}}, {198, {}}, {199, {}}, {200, {}}, {201, \
{}}, {202, {}}, {203, {}}, {204, {}}, {205, {}}, {206, {}}, {207, \
{}}, {208, {}}, {209, {}}, {210, {}}, {211, {}}, {212, {}}, {213, \
{}}, {214, {}}, {215, {}}, {216, {}}, {217, {}}, {218, {}}, {219, \
{}}, {220, {}}, 

{221, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket]}}, 

{222, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/2, 1/2} \[RightAngleBracket], 
 \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], 
 \[LeftAngleBracket] IC3\[Beta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], 
 \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
 \[LeftAngleBracket] IC4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
	
{223, {}}, {224, {}},

{225, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] IC2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket]}}, 
 	
{226, {}}, 

{227, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], 
 \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/4, 1/4, 1/4} \[RightAngleBracket], 
 \[LeftAngleBracket] IC2z, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 1/4, 1/4} \[RightAngleBracket], 
 \[LeftAngleBracket] IC2b, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/4, 1/4, 1/4} \[RightAngleBracket], 
 \[LeftAngleBracket] IC2x, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/4, 1/4, 1/4} \[RightAngleBracket], 
 \[LeftAngleBracket] IC3\[Alpha]i, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/4, 1/4, 1/4} \[RightAngleBracket], 
 \[LeftAngleBracket] IC3\[Delta], {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/4, 1/4, 1/4} \[RightAngleBracket], 
 \[LeftAngleBracket] IC4z, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/4, 1/4, 1/4} \[RightAngleBracket]}}, 
	
{228,
 {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {0, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {3/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4yi, {3/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma]i, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta], {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4x, {1/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/4, 3/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2c, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha]i, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {3/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {3/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2e, {0, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {3/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4y, {0, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 0, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/4, 0, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/4, 1/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2f, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2d, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {3/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/2, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {3/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {3/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {1/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/4, 3/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha]i, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {3/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {3/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2e, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {3/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 1/2} \[RightAngleBracket]}}, {229, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket]}
 },

 {229, {}}, {230, {}}} 



gtspacegroupelements = {{1, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket]}}, {2, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket]}}, {3, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket]}}, 
 {4, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 0} \[RightAngleBracket]}}, {5, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {6, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket]}}, {7, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {8, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {9, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {10, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket]}}, 
 {11, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 0} \[RightAngleBracket]}}, 
 {12, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket]}}, {13, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {14, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket]}}, 
 {15, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {16, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket]}}, 
 {17, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket]}}, 
 {18, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {19, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {20, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, {21, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {22, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {23, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {24, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {25, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket]}}, 
 {26, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket]}}, 
 {27, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {28, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 0} \[RightAngleBracket]}}, 
 {29, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket]}}, 
 {30, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket]}}, 
 {31, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket]}}, 
 {32, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {33, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {34, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {35, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, {36, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {37, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {38, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket]}}, 
 {39, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket]}}, {40, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {41, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket]}}, {42, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, {43, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {3/4, 3/4, 1/4} \[RightAngleBracket]}}, {44, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {45, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket]}}, {46, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket]}}, 
 {47, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket]}}, {48, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket]}}, 
 {49, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket]}}, {50, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 0} \[RightAngleBracket]}}, 
 {51, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 0} \[RightAngleBracket]}}, {52, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket]}}, 
 {53, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket]}}, {54, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket]}}, 
 {55, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, {56, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket]}}, 
 {57, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 0} \[RightAngleBracket]}}, {58, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {59, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 0} \[RightAngleBracket]}}, {60, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {61, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, {62, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {63, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {64, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {65, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {66, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {67, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {68, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket]}}, 
 {69, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {70, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {1/2, 3/4, 1/4} \[RightAngleBracket]}}, {71, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {72, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {73, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {74, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {75, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket]}}, 
 {76, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 3/4} \[RightAngleBracket]}}, 
 {77, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {78, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/4} \[RightAngleBracket]}}, 
 {79, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {80, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 1/4} \[RightAngleBracket]}}, 
 {81, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket]}}, 
 {82, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {83, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket]}}, 
 {84, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 1/2} \[RightAngleBracket]}}, {85, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 0} \[RightAngleBracket]}}, 
 {86, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 0, 1/2} \[RightAngleBracket]}}, {87, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {88, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {3/4, 3/4, 3/4} \[RightAngleBracket]}}, {89, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket]}}, 
 {90, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket]}}, {91, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/4} \[RightAngleBracket]}}, 
 {92, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {93, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket]}}, {94, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket]}}, 
 {95, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 3/4} \[RightAngleBracket]}}, {96, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {1/2, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {97, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {98, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {99, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket]}}, {100, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {101, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket]}}, {102, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket]}}, 
 {103, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket]}}, {104, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {105, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket]}}, {106, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {107, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {108, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {109, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {0, 1/2, 1/4} \[RightAngleBracket]}}, {110, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 0, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {0, 1/2, 3/4} \[RightAngleBracket]}}, {111, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket]}}, 
 {112, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket]}}, {113, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {114, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {115, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket]}}, 
 {116, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket]}}, {117, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {118, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {119, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {120, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2b, {1/2, 1/2, 0} \[RightAngleBracket]}}, {121, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {122, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 1/2, 1/4} \[RightAngleBracket]}}, 
 {123, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket]}}, 
 {124, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {125, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {126, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {127, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {128, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {129, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket]}}, 
 {130, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {131, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {132, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket]}}, 
 {133, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {134, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {135, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {136, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket]}}, 
 {137, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {138, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket]}}, 
 {139, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {140, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2b, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {141, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {3/4, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4z, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {3/4, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/4, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {1/4, 1/4, 3/4} \[RightAngleBracket]}}, {142, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {3/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/4, 1/4, 1/4} \[RightAngleBracket]}}, {143, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket]}}, 
 {144, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 2/3} \[RightAngleBracket]}}, {145, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 1/3} \[RightAngleBracket]}}, 
 {146, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {2/3, 1/3, 1/3} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {1/3, 2/3, 2/3} \[RightAngleBracket]}}, 
 {147, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket]}}, 
 {148, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IEe, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {2/3, 1/3, 1/3} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3zi, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/3, 2/3, 2/3} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3z, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {1/3, 2/3, 2/3} \[RightAngleBracket]}}, {149, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 0} \[RightAngleBracket]}}, {150, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket]}}, {151, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 2/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 0} \[RightAngleBracket]}}, {152, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 1/3} \[RightAngleBracket]}}, {153, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 1/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 0} \[RightAngleBracket]}}, {154, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 2/3} \[RightAngleBracket]}}, {155, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2D, {2/3, 1/3, 1/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2C, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {1/3, 2/3, 2/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2C, {1/3, 2/3, 2/3} \[RightAngleBracket]}}, 
 {156, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 0} \[RightAngleBracket]}}, 
 {157, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 0} \[RightAngleBracket]}}, 
 {158, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {159, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {160, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {2/3, 1/3, 1/3} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2A, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {1/3, 2/3, 2/3} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {1/3, 2/3, 2/3} \[RightAngleBracket]}}, {161, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {2/3, 1/3, 1/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C3zi, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {2/3, 1/3, 5/6} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {2/3, 1/3, 5/6} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {2/3, 1/3, 5/6} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/3, 2/3, 2/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C3z, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {1/3, 2/3, 1/6} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/3, 2/3, 1/6} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {1/3, 2/3, 1/6} \[RightAngleBracket]}}, 
 {162, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 0} \[RightAngleBracket]}}, 
 {163, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {164, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 0} \[RightAngleBracket]}}, 
 {165, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {166, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2D, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2y, {2/3, 1/3, 1/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C2C, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IEe, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {2/3, 1/3, 1/3} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {1/3, 2/3, 2/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2C, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {1/3, 2/3, 2/3} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3zi, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {1/3, 2/3, 2/3} \[RightAngleBracket]}}, 
 {167, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2D, {2/3, 1/3, 5/6} \[RightAngleBracket], \[LeftAngleBracket] C2y, {2/3, 1/3, 5/6} \[RightAngleBracket], 
   \[LeftAngleBracket] C2C, {2/3, 1/3, 5/6} \[RightAngleBracket], \[LeftAngleBracket] IEe, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {2/3, 1/3, 1/3} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {2/3, 1/3, 5/6} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {2/3, 1/3, 5/6} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {2/3, 1/3, 5/6} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3z, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {1/3, 2/3, 2/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {1/3, 2/3, 1/6} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/3, 2/3, 1/6} \[RightAngleBracket], \[LeftAngleBracket] C2C, {1/3, 2/3, 1/6} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {1/3, 2/3, 2/3} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3zi, {1/3, 2/3, 2/3} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {1/3, 2/3, 1/6} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/3, 2/3, 1/6} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {1/3, 2/3, 1/6} \[RightAngleBracket]}}, 
 {168, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 0} \[RightAngleBracket]}}, 
 {169, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 5/6} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/6} \[RightAngleBracket]}}, 
 {170, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/6} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 5/6} \[RightAngleBracket]}}, 
 {171, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/3} \[RightAngleBracket]}}, 
 {172, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 2/3} \[RightAngleBracket]}}, 
 {173, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {174, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 0} \[RightAngleBracket]}}, 
 {175, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 0} \[RightAngleBracket]}}, 
 {176, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {177, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 0} \[RightAngleBracket]}}, 
 {178, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 5/6} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/6} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 5/6} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 1/6} \[RightAngleBracket]}}, 
 {179, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/6} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 5/6} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 1/6} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 5/6} \[RightAngleBracket]}}, 
 {180, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 1/3} \[RightAngleBracket]}}, 
 {181, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 2/3} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 2/3} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 1/3} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 2/3} \[RightAngleBracket]}}, 
 {182, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {183, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 0} \[RightAngleBracket]}}, 
 {184, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {185, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 0} \[RightAngleBracket]}}, 
 {186, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {187, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 0} \[RightAngleBracket]}}, 
 {188, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 0} \[RightAngleBracket]}}, 
 {189, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 0} \[RightAngleBracket]}}, 
 {190, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {191, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2B, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 0} \[RightAngleBracket]}}, 
 {192, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {193, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 0} \[RightAngleBracket]}}, 
 {194, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2D, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2C, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2B, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2A, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC6z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2B, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2A, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2D, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2C, {0, 0, 1/2} \[RightAngleBracket]}}, 
 {195, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket]}}, 
 {196, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {197, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {198, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket]}}, 
 {199, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 1/2, 0} \[RightAngleBracket]}}, {200, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket]}}, {201, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {202, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket]}}, {203, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {3/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 0, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/4, 0, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma], {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha]i, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {3/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/4, 1/2, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/4, 1/2, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2z, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma]i, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta], {3/4, 3/4, 0} \[RightAngleBracket]}}, {204, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {205, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket]}}, {206, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 1/2, 0} \[RightAngleBracket]}}, {207, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket]}}, {208, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {209, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2b, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {210, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {3/4, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2e, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {3/4, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2b, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2e, {3/4, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2f, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {3/4, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2d, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/4, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2b, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2f, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2c, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2d, {3/4, 3/4, 1/4} \[RightAngleBracket]}}, {211, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {212, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {3/4, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2e, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/4, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/4, 1/4, 1/4} \[RightAngleBracket]}}, 
 {213, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/4, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {3/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {3/4, 3/4, 3/4} \[RightAngleBracket]}}, 
 {214, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/4, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {3/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2b, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2e, {3/4, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2f, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2c, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {3/4, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2d, {1/4, 1/4, 1/4} \[RightAngleBracket]}}, {215, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {0, 0, 0} \[RightAngleBracket]}}, {216, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2d, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2f, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4xi, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2c, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/2, 1/2, 0} \[RightAngleBracket]}}, {217, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2c, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {218, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {219, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2e, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4xi, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2f, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2c, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 0, 1/2} \[RightAngleBracket]}}, {220, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/4, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/4, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2e, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/4, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {3/4, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4xi, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/4, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2c, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {3/4, 1/4, 1/4} \[RightAngleBracket]}}, {221, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket]}}, {222, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, {223, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {224, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2e, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket]}}, 
 {225, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2e, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2e, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 0} \[RightAngleBracket]}}, {226, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4xi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4x, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2e, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2d, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2f, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 1/2} \[RightAngleBracket]}}, {227, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/4, 1/2, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/4, 1/2, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 1/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {3/4, 3/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha]i, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2e, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {3/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2x, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2b, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2f, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2d, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma], {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4z, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4xi, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {3/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {3/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4yi, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma]i, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta], {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4x, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {3/4, 1/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2c, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 0} \[RightAngleBracket]}}, 
 {228, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4zi, {0, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {3/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4yi, {3/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {0, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Gamma]i, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta], {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4x, {1/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/4, 3/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2c, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] Ee, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2y, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha]i, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2a, {3/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {3/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C2e, {0, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {3/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4y, {0, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/4, 0, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2x, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/4, 0, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/4, 1/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2b, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2f, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {0, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2d, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {1/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {3/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {3/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {3/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {3/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {3/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/2, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {3/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/4, 0, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {3/4, 3/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2z, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {1/4, 0, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {3/4, 1/4, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4z, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {1/4, 3/4, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/4, 3/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {3/4, 1/2, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IEe, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha]i, {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {3/4, 0, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {0, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/4, 3/4, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {3/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {3/4, 1/4, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2e, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/4, 1/2, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {3/4, 1/4, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/4, 1/2, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {1/2, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 1/2} \[RightAngleBracket]}}, {229, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4x, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {0, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/2, 1/2, 1/2} \[RightAngleBracket]}}, 
 {230, {\[LeftAngleBracket] Ee, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2z, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2a, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4z, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {1/4, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4x, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2e, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4xi, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {3/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C2c, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IEe, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2y, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta]i, {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha]i, {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {0, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {0, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha], {1/2, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {1/2, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2a, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4zi, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {1/4, 3/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC2e, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4yi, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {3/4, 3/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4y, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {1/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] Ee, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C2x, {0, 0, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta]i, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Gamma]i, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha]i, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] C3\[Gamma], {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] C3\[Alpha], {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] C3\[Beta], {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] C2a, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2b, {1/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4z, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4zi, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C4x, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2e, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2f, {1/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] C4xi, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4yi, {1/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2c, {3/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] C4y, {3/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] C2d, {1/4, 1/4, 1/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IEe, {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC2z, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2y, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2x, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta]i, {1/2, 1/2, 1/2} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Beta]i, {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma]i, {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Alpha]i, {1/2, 0, 0} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Delta], {1/2, 1/2, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Gamma], {1/2, 0, 0} \[RightAngleBracket], 
   \[LeftAngleBracket] IC3\[Alpha], {0, 0, 1/2} \[RightAngleBracket], \[LeftAngleBracket] IC3\[Beta], {0, 1/2, 0} \[RightAngleBracket], \[LeftAngleBracket] IC2a, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2b, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4z, {1/4, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4zi, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC4x, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2e, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2f, {3/4, 3/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4xi, {1/4, 1/4, 3/4} \[RightAngleBracket], 
   \[LeftAngleBracket] IC4yi, {3/4, 1/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2c, {1/4, 1/4, 3/4} \[RightAngleBracket], \[LeftAngleBracket] IC4y, {1/4, 3/4, 1/4} \[RightAngleBracket], \[LeftAngleBracket] IC2d, {3/4, 3/4, 3/4} \[RightAngleBracket]}}}

Protect[gtspacegroupelements]
(*--- end space group elements *)

(*--- Colors for the tables ---*)
Unprotect[GTBackGroundColor1, GTBackGroundColor2, GTCornerColor, GTDividerColor1, GTDividerColor2];
GTBackGroundColor1 = Pink;
GTBackGroundColor2 = LightGray;
GTCornerColor      = Yellow;
GTDividerColor1    = Red;
GTDividerColor2    = Blue;
Protect[GTBackGroundColor1, GTBackGroundColor2, GTCornerColor, GTDividerColor1, GTDividerColor2];
(*--- End colors for the tables ---*)


On[General::"spell1"]; On[General::"spell"];

Begin["`Private`"]
(*InitSymbols:=Module[{}];*)

End[ ]

EndPackage[ ]
