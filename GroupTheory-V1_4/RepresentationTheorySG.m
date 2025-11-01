(* ::Package:: *)

(****z* /RepresentationTheorySG.m
!
! NAME
!  RepresentationTheorySG.m
! AUTHOR
!  M. Geilhufe, W. Hergert
! MODIFICATION HISTORY
!  10/27/18 : initial documentation 
! USAGE
!  Contains modules to calculate space group representations
! 
! ERROR MESSAGES
!  
! GTPack MODULES

! GTPack NOTEBOOKS 
!  GTPack_Dev/Modules_New/Space_Group_Representations/SGRep.nb
! DESCRIPTION
!  RepresentationTheorySG.m contains all modules to calculate space group representations.
!
! LITERATURE
! 
! TODO
!  i)   Documentation pages
!  ii)  Embed into GTPack
!  iii) Cross check against bilbao crystallographic server
! PROBLEMS
!
!--------------------------------------------------------------------------------
!
***)

BeginPackage["GroupTheory`RepresentationTheorySG`",{"GroupTheory`Basic`","GroupTheory`RepresentationTheory`","GroupTheory`Auxiliary`","GroupTheory`Symbols`","GroupTheory`Lattice`"}]
  GTCharacterTableOfK ::usage = "GTCharacterTableOfK[\*StyleBox[\"space group, k-point, reciprocal basis\", \"TI\"]] gives the character table of a \*StyleBox[\"space group\", \"TI\"] at the wave vector k. (The command was replaced by the new command GTSGCharacterTable)"
  GTSGCharacterTable ::usage = "GTSGCharacterTable[\*StyleBox[\"space group, k-point, reciprocal basis\", \"TI\"]] gives the character table of a \*StyleBox[\"space group\", \"TI\"] at the wave vector k"  
  GTSGCosetRepresentative ::usage = "GTSGCosetRepresentative[\*StyleBox[\"space group, subgroup, basis\", \"TI\"]] determines a coset representative for a left coset decomposition of the factor group G/T (G being a \*StyleBox[\"space group\", \"TI\"]], T the group of pure lattice translations) into a normal \*StyleBox[\"sub group\", \"TI\"]] of index 2 or 3"  
  GTSGClasses ::usage = "GTSGClasses[\*StyleBox[\"space group, basis\", \"TI\"]] gives the classes of a factor group of the factor group G/T (G being a \*StyleBox[\"space group\", \"TI\"], T the group of pure lattice translations)."
  GTSGgmat ::usage = "GTSGgmat[\*StyleBox[\"element1, element2, basis\", \"TI\"]] gives the product of two coset representatives of the factor group G/T (G being a \*StyleBox[\"space group\", \"TI\"], T the group of pure lattice translations)of a space group"
  GTSGGetInvSubGroup ::usage = "GTSGGetInvSubgroup[\*StyleBox[\"space group, classes, basis, index\", \"TI\"]] gives an invariant subgroup of the factor group G/T (G being a \*StyleBox[\"space group\", \"TI\"], T the group of pure lattice translations) with a specified index"
  GTSGGetIreps ::usage = "GTSGGetIreps[\*StyleBox[\"space group, k-point, reciprocal basis\", \"TI\"]] gives the character table and the representation matrices of a \*StyleBox[\"space group\", \"TI\"] at a given wave vector k."
  GTSGLeftCosets  ::usage  = "GTSGLeftCosets[\*StyleBox[\"space group, subgroup, basis\", \"TI\"]] gives the left cosets of the factor group G/T (G being a \*StyleBox[\"space group\", \"TI\"], T the group of pure lattice translations) with respect to a subgroup"
  GTSGOrderOfElement ::usage ="GTSGOrderOfElement[\*StyleBox[\"element, basis\", \"TI\"]] gives the order of a space group \*StyleBox[\"element\", \"TI\"]"
  GTSGRightCosets ::usage = "GTSGRightCosets[\*StyleBox[\"space group, subgroup, basis\", \"TI\"]] gives the right cosets of the factor group G/T (G being a \*StyleBox[\"space group\", \"TI\"], T the group of pure lattice translations) with respect to a subgroup"
  
  GTStarOfK                 ::usage = "GTStarOfK[\*StyleBox[\"group,\", \"TI\"]\*StyleBox[\"k\",\"TI\",FontWeight->\"Bold\"]\*StyleBox[\"vector, \", \"TI\"]\*StyleBox[\"reciprocal basis,\", \"TI\"]] gives the star of the wave vector \*StyleBox[\"k\", \"TI\",FontWeight->\"Bold\"] for a given \*StyleBox[\"group\", \"TI\"]."

(*  GTSGInduceOrbitLengthOne ::usage = "to be removed (internal)"
  GTSGInduceOrbitLengthTwoThree ::ssage = "to be removed (internal)"
  GTSGGetUnitaryTrafoMatrix ::ssage = "to be removed (internal)" *)

  Options[GTCharacterTableOfK]  = {GOVerbose->True,GOGeneralPositionIreps->True}
  Options[GTSGCharacterTable]  = {GOVerbose->True,GOGeneralPositionIreps->True}  
  Options[GTSGgmat]      = {GOTakeMod->True}
  Options[GTSGGetInvSubGroup]      ={MaxIterations -> 1000, GOVerbose -> True};
  Options[GTSGGetIreps]  = {GOVerbose->True,GOGeneralPositionIreps->True}
  Options[GTStarOfK]         = {GOFast->GOFastValue,GOLeftCosets->False}
Begin["`Private`"] (* Begin Private Context *) 

(****z* /GTCharacterTableOfK
! NAME
!  GTCharacterTableOfK
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   old version in 2016, disappeared in the current version, but was erroneous anyhow
!   11/02/2018: replaced by GTSGCharacterTable
! USAGE
!   GTCharacterTableOfK gives the character table of a space group at the wave vector k. (The command was replaced by the new command GTSGCharacterTable)
! INPUT
!   space group
!   k-vector
!   reciprocal basis vectors
! OUTPUT
!  {classes of G(k), characters, names of the ireps}
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTSGCharacterTable
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
  GTCharacterTableOfK[grp_,kvec_,recbas_,OptionsPattern[]] := Module[{},
   Print["GTCharacterTableOfK was replaced by the new command GTSGCharacterTable"];
   GTSGCharacterTable[grp,kvec,recbas,GOVerbose->OptionValue[GOVerbose],GOGeneralPositionIreps->OptionValue[GOGeneralPositionIreps]]
  ]
(*
***)

(****z* /GTSGCharacterTable
! NAME
!  GTSGCharacterTable
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   12/06/2018: first version
! USAGE
!   GTSGCharacterTable gives the character table of a space group at the wave vector k
! INPUT
!   space group
!   k-vector
!   reciprocal basis vectors
! OUTPUT
!  {classes of G(k), characters, names of the ireps}
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTSGGetIreps, GTSGPrintGeneral, GTSGPrintChars
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
  GTSGCharacterTable[grp_,kvec_,recbas_,OptionsPattern[]] := Module[{classes,chars,names,grpk,ireps,gr,genireps},
  If[OptionValue[GOGeneralPositionIreps],
   {{classes,chars,names},{grpk,ireps},{gr,genireps}}=GTSGGetIreps[grp,kvec,recbas,GOVerbose->False];
   If[OptionValue[GOVerbose],GTSGPrintGeneral[grpk,kvec,recbas]];
   If[OptionValue[GOVerbose],GTSGPrintChars[grpk,classes,chars,OptionValue[GOGeneralPositionIreps]]],
   {{classes,chars,names},{grpk,ireps}}=GTSGGetIreps[grp,kvec,recbas,GOVerbose->False,GOGeneralPositionIreps->False];
   If[OptionValue[GOVerbose],GTSGPrintGeneral[grpk,kvec,recbas]];
   If[OptionValue[GOVerbose],GTSGPrintChars[grp,classes,chars,OptionValue[GOGeneralPositionIreps]]];
   ];
   Return[{classes,chars,names}]
  ]
(*
***)

(****z* /GTSGCosetRepresentative
! NAME
!  GTSGCosetRepresentative
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   11/02/2018: first version
! USAGE
!   GTSGCosetRepresentative determines a coset representative for a left coset decomposition of a factor group of a space group into a sub group of index 2 or 3
! INPUT
!   factor group of a space group
!   a subgroup of index 2 or 3
!   reals space basis vectors
! OUTPUT
!  a space group element q which serves as coset representative according to G = H + q*H (+ q^2*H)
!  permutation for the map between G and G = H + q H (+ q^2 H)
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTSGOrderOfElement, GTSGgmat
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
GTSGCosetRepresentative[grp_, grpn_, basis_] :=  Module[{n, gln, comp, qel, reordgroup, perm, qelsq,Tm,groupfound},
  (*Calculate Index and Length of group*)
  n = Length[grp]/Length[grpn];
  gln = Length[grpn];
  (*Estimate the complement of group G and subgroup H to write G = 
  H + qH (+q^2 H); 
  *)

  comp = Complement[grp, grpn];
  reordgroup={};
  groupfound=False;
  Tm=0;
  While[Not[groupfound] && Tm<gln,
  Tm=Tm+1;
  qel=comp[[Tm]];
  qelsq = GTSGgmat[qel, qel, basis];
  (*write G = H + qH + (q^2H);*)
  reordgroup =
   Which[
    (*Index 2*)
    n == 2,
    Flatten[{grpn, 
      Table[GTSGgmat[qel, grpn[[T]], basis], {T, 1, gln}]}, 1],
    (*Index 3*)
    n == 3,
    Flatten[{grpn, 
      Table[GTSGgmat[qel, grpn[[T]], basis], {T, 1, gln}], 
      Table[GTSGgmat[qelsq, grpn[[T]], basis], {T, 1, gln}]}, 1]
    ];
    groupfound=Length[Complement[reordgroup, grp]] === 0;
    ];
    If[Tm>gln,Print["Error: GTSGCosetRepresentative: Could not identify coset representative!"];Abort[]];
  (*Find the corresponding Permutation of the reordered group*)
  perm = FindPermutation[reordgroup, grp];
  Return[{qel, perm}];
  ]
(*
***)

(****z* /GTSGClasses
! NAME
!  GTSGClasses
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   12/05/2018: first version
! USAGE
!   GTSGClasses gives the classes of a factor group of the factor group G/T (G being a space group, T the group of pure lattice translations).
! INPUT
!   factor group of a space group
!   reals space basis vectors
! OUTPUT
!   list of classes
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTInverseElement, GTSGgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTSGGetSubGroup to determine invariant subgroups from complete classes
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
GTSGClasses[grp_,basis_]:=Module[{g,tb1,grpi,classes,grpleft,tmp,possym,cltmp,el,pos},
g=Length[grp];
grpi=Table[GTInverseElement[grp[[i]],GOLattice->basis],{i,1,g}];
tb1={};
grpleft=grp;
While[Length[grpleft]>0,
tmp=Union[
Table[
el=GTSGgmat[GTSGgmat[grp[[i]],grpleft[[1]],basis],grpi[[i]],basis];
pos=Position[grp[[;;,1]],el[[1]]][[1,1]];
grp[[pos]]
,{i,1,g}]];
cltmp={};
Do[
pos = Position[grpleft[[;; , 1]], tmp[[i, 1]]];
If[Length[pos]>0,
cltmp=Append[cltmp,tmp[[i]]];
grpleft=Delete[grpleft,pos[[1,1]]];
]
,{i,1,Length[tmp]}];
tb1=Append[tb1,cltmp];
];
classes=Map[Union,{Map[Union,tb1]}][[1]];
possym=Position[classes,\[LeftAngleBracket]Ee,{0,0,0}\[RightAngleBracket]];
If[Length[possym]==0,possym=Position[classes,\[LeftAngleBracket]IdentityMatrix[Length[GTGetMatrix[grp[[1,1]]]]],{0,0,0}\[RightAngleBracket]]];
If[First@First@possym>1,classes=RotateLeft[classes,First@First@possym-1]];
Return[classes]
]
(*
***)

(****z* /GTSGgmat
! NAME
!  GTSGgmat
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   11/02/2018: first version
! USAGE
!   GTSGgmat gives the product of two coset representatives of a space group
! INPUT
!   two space group elements and a real space basis
! OUTPUT
!   a space group element (product of the input)
! GTPack OPTIONS
!   GOTakeMod: determines if the real product is taken or the product of the coset representatives of the factor group
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTSGGetIreps for the induction of two representations from an irreducible representationwith orbit 1
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

GTSGgmat[g1_, g2_, latt_,OptionsPattern[]] := Module[{t1, t2, r1, g, trans, c, list, lt, sol, dc},
 	
  (*Estimate translation vectors in Cartesian coordinates*)
  t1 = g1[[2]].latt;
  t2 = g2[[2]].latt;
  (*Calculate final translation and transform back into direct coordinates*)
  r1 = GTGetRotationMatrix[g1[[1]]];
  trans = r1.t2 + t1;
  lt = Length[trans];
  list = Sum[c[i] latt[[i]], {i, 1, lt}];
  sol = Solve[Table[list[[i]] == trans[[i]], {i, 1, lt}],Table[c[i], {i, 1, lt}]];
  dc = Table[c[i], {i, 1, lt}] /. sol[[1]];
  (*Calculate rotational part*)
  g = GTgmat[g1[[1]], g2[[1]]];
  (*result, take translation mod 1 to keep vectors like (0,1/2,
  0) but get rid of vectors like (0,1,
  0)*)
  If[OptionValue[GOTakeMod],
  Return[\[LeftAngleBracket]g, Mod[dc, 1] Sign[dc]\[RightAngleBracket]],
  Return[\[LeftAngleBracket]g, dc\[RightAngleBracket]]]
  ]

(*
***)


(****z* /GTSGGetInvSubGroup
! NAME
!  GTSGGetInvSubGroup
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   11/02/2018: first version
! USAGE
!   GTSGGetInvSubgroup gives an invariant subgroup of a factor group of G/T (G being a space group, T the group of pure lattice translations) with a specified index
! INPUT
!   factor group G/T, classes, real space basis, index
! OUTPUT
!   an invariant subgroup of index n of the given factor group
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTSGgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTSGGetIreps for the construction of a subgroup of index 2 or 3
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


GTSGGetInvSubGroup[grp_, classes_, basis_, order_, OptionsPattern[]] := Module[{lg, pos, reducedclasses, reducedclassmult, lcl, numb, ln, 
   orderpos, lorderpos, nn, llfin, eq, eqint, const, opt, erg, poscl, C, it, cc, case, sets, choice, compare, grpout, groupfound},
  lg = Length[grp];
  
  (*remove the class containing the identity element*)
  pos = {Position[classes, \[LeftAngleBracket] Ee,{0,0,0}\[RightAngleBracket]], 
     Position[classes, \[LeftAngleBracket]IdentityMatrix[Length[GTGetMatrix[classes[[1, 1,1]]]]],{0,0,0}\[RightAngleBracket]], 
     Position[classes, \[LeftAngleBracket]GTGetEulerAngles[Ee],{0,0,0}\[RightAngleBracket]], Position[classes, \[LeftAngleBracket]GTGetQuaternion[Ee],{0,0,0}\[RightAngleBracket]]};
    
  reducedclasses = Delete[classes, pos[[First@First@Position[Map[Length, pos], 1], 1, 1]]];
  (*estimate the class lengths, positions of classes with a specified length, and numbers of classes with specified length*)
  
  lcl = Map[Length, reducedclasses];
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
     poscl = Position[reducedclasses[[;;,;;,1]], GTSGgmat[reducedclasses[[i, m]], reducedclasses[[j, n]],basis][[1]]];
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
GTSGGetInvSubGroup[grp_,cl_,basis_,n_]:=Module[{nel,possym,cls,lcl,cllist,clrule,groupfound,mm,subs,pos,ii,rule,grptest,testgroup,totallist,ng,multtable,el},
(*Final number of elements*)
ng=Length[grp];
nel=ng/n;
(*Search for the identity, remove it from the list of classes*)
possym=Position[cl,\[LeftAngleBracket]Ee,{0,0,0}\[RightAngleBracket]];
If[Length[possym]==0,possym=Position[cl,\[LeftAngleBracket]IdentityMatrix[Length[GTGetMatrix[grp[[1,1]]]]],{0,0,0}\[RightAngleBracket]]];
cls=Delete[cl,possym[[1,1]]];
(*If nel = 1 return identity element as subgroup*)
If[nel==1,Return[cl[[possym[[1,1]]]]]];
(*Create dummy lists to work on, representing the classes*)
lcl=Length[cls];
cllist=Table[c[i],{i,1,lcl}];
clrule=Table[c[i]->Length[cls[[i]]],{i,1,lcl}];
groupfound=False;
(*form combinations of classes and querry for the requested length nel, check if the result is a group*)
mm=1;
totallist={0};
multtable=Table[
el=GTSGgmat[grp[[i]],grp[[j]],basis];
pos=Position[grp[[;;,1]],el[[1]]][[1,1]];
grp[[pos]]
,{i,1,ng},{j,1,ng}];
While[Not[groupfound]&&mm<=lcl&&Min[totallist]<=nel,
 subs=Subsets[cllist,{mm}];
 pos=Position[totallist=Map[Total,subs]/.clrule,nel-1];
 If[Length[pos]>0,
   ii=1;
   While[ii<=Length[pos]&&Not[groupfound],
     rule=Table[c[k]->cls[[k]],{k,1,lcl}];
     grptest=Flatten[Prepend[subs[[pos[[ii,1]]]]/.rule,cl[[possym[[1,1]]]]]];
     If[AllTrue[Divisible[nel,GTSGOrderOfElement[#,basis]]&/@grptest,TrueQ],
        testgroup=Union@Flatten@Table[multtable[[First@First@Position[grp,grptest[[i]]],First@First@Position[grp,grptest[[j]]]]],{i,1,nel},{j,1,nel}];
        groupfound=(Length[Intersection[testgroup,grptest]]==Length[grptest])&&Length[testgroup]==Length[grptest],
        groupfound=False];
     ii=ii+1;
   ];
  ,mm=mm+1]
];
If[groupfound,Return[grptest],Return[{}]];
]
*)

(*
***)

(****z* /GTSGGetIreps
! NAME
!  GTSGGetIreps
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   12/06/2018: first version
! USAGE
!   GTSGGetIreps gives the character table and the representation matrices of a space group at a given wave vector k.
! INPUT
!   space group, kvector, reciprocal basis
! OUTPUT
!   {{classes of G(k)/T, characters, names of ireps},{G(k),representation matrices }}
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTReciprocalBasis, GTGroupOfK, GTSGPrintGeneral, GTSGPrintChars, GTSGPrintIreps, GTSGgmat, GTSGClasses, GTSGGetInvSubGroup, GTSGCosetRepresentative, GTSGGetIreps, GTInverseElement, GTSGInduceOrbitLengthOne, GTSGInduceOrbitLengthTwoThree
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
GTSGGetIreps[group_,kvec_,recbas_,OptionsPattern[]]:=Module[{grpsym,cttable,lengthone,grp,classes,basis,ireps,irepsn,n,grpn,gl,gln,el,orbitq,orbitqsq,pos,qel,perm,qelsq, qelqu, qelinv,qelsqinv,irep,irepqsq,irepq,expfacq,expfacqsq,chars,charscomp,irepsfinal,charsfinal,clt,star, leftcos, genireps},
(*Initial information*)
basis=GTReciprocalBasis[recbas];
(*calculate group of k*)
grp=GTGroupOfK[group,kvec.recbas,recbas,GOFast->True];
{star,leftcos}=GTStarOfK[group,kvec.recbas,recbas,GOLeftCosets->True];
(**********************************)
(*Irreducible representation of SG1*)
(**********************************)
If[Length[grp]==1,
ireps={{{{1}}}};
charsfinal={{1}};
classes={grp};
(*General position ireps*)
If[OptionValue[GOGeneralPositionIreps],
genireps=GTGetGeneralPosIreps[group,grp,leftcos,ireps,basis,recbas,kvec];
charsfinal=Table[Table[Tr[genireps[[h,T]]],{T,1,Length[group]}],{h,1,Length[genireps]}];
classes=GTSGClasses[group,basis],
genireps={}
];
If[OptionValue[GOVerbose],GTSGPrintGeneral[grp,kvec,recbas,star]];
If[OptionValue[GOVerbose],GTSGPrintIreps[grp,ireps,genireps,group]];
If[OptionValue[GOVerbose],GTSGPrintChars[grp,classes,charsfinal,OptionValue[GOGeneralPositionIreps]]];
If[OptionValue[GOGeneralPositionIreps],
Return[{{classes,charsfinal,{Superscript["\[CapitalGamma]", ToString[1]]}},{grp,ireps},{group,genireps}}],
Return[{{classes,charsfinal,{Superscript["\[CapitalGamma]", ToString[1]]}},{grp,ireps}}]]];

(*****************************************)
(*Irreducible representation of other SGs *)
(*****************************************)

classes=GTSGClasses[grp,basis];

(*Check if the group is symmorphic*)
(*-------------------------*)

If[Total[Map[Norm,grp[[;;,2]]]] == 0,
	(*use conventional method for generating representation matrices*)
	If[OptionValue[GOVerbose],Print["Symmorphic Group"]];
	grpsym = grp[[;;,1]];
	gl = Length[grpsym];
	cttable = GTCharacterTable[grpsym,GOFast->True,GOVerbose->False];
	ireps = Table[GTGetIrep[grpsym,ir,cttable,GOFast->True],{ir,1,Length[cttable[[1]]]}];
	,
	
	(*If not symmorphic, continue from here*)
	
	If[OptionValue[GOVerbose],Print["Nonsymmorphic Group"]];
	
	(*Start with subgroup chain*)
	(*-------------------------*)
	(*check if group has a symorphic invariant subgroup of order 2 or 3*)
	grpsym=Select[grp, Norm[#[[2]]] == 0 &];
	grpn={};
	If[Length[grp]/Length[grpsym]==2||Length[grp]/Length[grpsym]==3,
		If[GTGroupQ[grpsym[[;;,1]]],
		(*Print["Symmorphic invariant subgroup"];*)
		n=Length[grp]/Length[grpsym];
		If[Length[grpsym]==Length[Union[Flatten[Table[classes[[Position[classes, T][[1, 1]]]], {T, grpsym}]]]],
			grpn=grpsym]];
	];
	(*Alternatively, check if divisible by 3 or 2; 3 is preferred to keep the group order small*)
	If[Length[grpn]==0,
		n=If[Divisible[Length[grp],2],2,3];
		grpn=GTSGGetInvSubGroup[grp,classes,basis,n,GOVerbose->False];
		If[Length[grpn]==0&&n==2,
		n=3;
		grpn=GTSGGetInvSubGroup[grp,classes,basis,n,GOVerbose->False]];
	];
	(*Calculate ireps of the subgroup*)
	irepsn=Last@Last@GTSGGetIreps[grpn,kvec,recbas,GOVerbose->False,GOGeneralPositionIreps->False];
	(*Calculate Coset decomposition*)
	{qel,perm}=GTSGCosetRepresentative[grp,grpn,basis];
	qelsq=GTSGgmat[qel,qel,basis,GOTakeMod->False];
	qelqu=GTSGgmat[qel,qelsq,basis,GOTakeMod->False];
	qelinv=GTInverseElement[qel,GOLattice->basis];
	qelsqinv=GTInverseElement[qelsq,GOLattice->basis];
	
	(*Induce representations*)
	(*----------------------*)
	(*Calculate orbits*)
	gl=Length[grp];
	gln=Length[grpn];
    expfacq={};
	orbitq=Table[
		el=GTSGgmat[GTSGgmat[qel,grpn[[j]],basis,GOTakeMod->False],qelinv,basis,GOTakeMod->False];
		pos=First@First@Position[grpn,el[[1]]];
		expfacq=Append[expfacq,Exp[-I ((el[[2]]-grpn[[pos,2]]).basis).(kvec.recbas)]];
		pos
	,{j,1,gln}];
	expfacqsq={};
	orbitqsq=Table[
		el=GTSGgmat[GTSGgmat[qelsq,grpn[[j]],basis,GOTakeMod->False],qelsqinv,basis,GOTakeMod->False];
		pos=First@First@Position[grpn,el[[1]]];
		expfacqsq=Append[expfacqsq,Exp[-I ((el[[2]]-grpn[[pos,2]]).basis).(kvec.recbas)]];
		pos
	,{j,1,gln}];
	(*Scan through all ireps, check length of orbits and induce representations*)
	ireps=FullSimplify[Flatten[Table[
	irep=irepsn[[h]];
	irepq=Table[irep[[orbitq[[T]]]]*expfacq[[T]],{T,1,gln}];
	irepqsq=Table[irep[[orbitqsq[[T]]]]*expfacqsq[[T]],{T,1,gln}];
	lengthone=If[n==2,
	Length[Union[{Map[Tr,irep],Map[Tr,irepq]}]]==1,
	Length[Union[{Map[Tr,irep],Map[Tr,irepq],Map[Tr,irepqsq]}]]==1
	];
	If[lengthone,
	(*Length of orbit is 1*)
	GTSGInduceOrbitLengthOne[grp,grpn,qel,qelsq,qelqu,perm,irep,irepq,basis,kvec,recbas],
	(*Length of orbit is 2 or 3*)
	{GTSGInduceOrbitLengthTwoThree[grp,grpn,qel,qelsq,qelqu,perm,irep,irepq,irepqsq,basis,kvec,recbas]}]
	,{h,1,Length[irepsn]}],1]];
];
chars=Table[Table[Tr[ireps[[h,T]]],{T,1,gl}],{h,1,Length[ireps]}]//FullSimplify;
charscomp=SortBy[Union[chars],First];
(*Final Ireps*)
irepsfinal=Table[
pos=First@First@Position[chars,charscomp[[i]]];
ireps[[pos]],{i,1,Length[charscomp]}];
(*General position ireps*)
If[OptionValue[GOGeneralPositionIreps],
genireps=GTGetGeneralPosIreps[group,grp,leftcos,irepsfinal,basis,recbas,kvec];
chars=Table[Table[Tr[genireps[[h,T]]],{T,1,Length[group]}],{h,1,Length[genireps]}];
classes=GTSGClasses[group,basis];
(*Final chars full group*)
clt=Transpose[chars];
charsfinal=FullSimplify[Transpose[Table[pos=First@First@Position[group,classes[[i,1]]];
clt[[pos]]
,{i,1,Length[classes]}]]];
,
genireps={};
(*Final chars little group*)
clt=Transpose[charscomp];
charsfinal=FullSimplify[Transpose[Table[pos=First@First@Position[grp,classes[[i,1]]];
clt[[pos]]
,{i,1,Length[classes]}]]];
];
If[OptionValue[GOVerbose],GTSGPrintGeneral[grp,kvec,recbas,star]];
If[OptionValue[GOVerbose],GTSGPrintIreps[grp,irepsfinal,genireps,group]];
If[OptionValue[GOVerbose],GTSGPrintChars[grp,classes,charsfinal,OptionValue[GOGeneralPositionIreps]]];
If[OptionValue[GOGeneralPositionIreps],
Return[{{classes,charsfinal,Table[Superscript["\[CapitalGamma]", ToString[i]],{i,1,Length[charsfinal]}]},{grp,irepsfinal},{group,genireps}}],
Return[{{classes,charsfinal,Table[Superscript["\[CapitalGamma]", ToString[i]],{i,1,Length[charsfinal]}]},{grp,irepsfinal}}]]
]

(*
***)

(****z* /GTGetGeneralPosIreps
! NAME
!  GTGetGeneralPosIreps
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   12/12/2018: first version
! USAGE
!   GTGetGeneralPosIreps induces the irreducible representations for the space group from a group of k
! INPUT
!   space group, group of k, left cosets, ireps, basis
! OUTPUT
!   representation matrices
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTgmat, GTInverseElement
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTSGGetIreps
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


GTGetGeneralPosIreps[grp_,grpk_,leftcos_,ireps_,basis_,recbas_,kvec_]:=Module[{lc,nm,nireps,pos,indmat,irtmp,el,lci,lcj,expf},
(*Definitions*)
lc=Length[leftcos];
nm=Table[0,{i,1,lc},{j,1,lc}];
nireps=Length[ireps];
(*Induce matrices*)
Transpose@Table[
Sum[
Do[
lci=GTInverseElement[leftcos[[i,1]],GOLattice->basis];
lcj=leftcos[[j,1]];
el=GTSGgmat[GTSGgmat[lci,grp[[T]],basis,GOTakeMod->False],lcj,basis,GOTakeMod->False];
pos=Position[grp[[;;,1]],el[[1]]][[1,1]];
expf=Exp[-I ((el[[2]]-grp[[pos,2]]).basis).(kvec.recbas)];
pos=Position[grpk[[;;,1]],el[[1]]];
If[Length[pos]>0,
indmat=nm;
indmat[[i,j]]=1;
irtmp=Table[KroneckerProduct[indmat,expf*ireps[[h,First@First@pos]]],{h,1,nireps}]];
(*irtmp=Table[KroneckerProduct[indmat,ireps[[h,First@First@pos]]],{h,1,nireps}]];*)
,{j,1,lc}];
irtmp,
{i,1,lc}]
,{T,1,Length[grp]}]
]

(*
***)

(****z* /GTSGLeftCosets
! NAME
!  GTSGLeftCosets
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   11/02/2018: first version
!   09/09/2020: enforce ordering of the cosets with respect to the order of the element to avoid nontrivial coset representatives (taken as first element of the list) for trivial coset decompositions G = E G
! USAGE
!   GTSGLeftCosets gives the left cosets of the factor group G/T (G being a space group, T the group of pure lattice translations) with respect to a subgroup
! INPUT
!   two factor groups and a real space basis
! OUTPUT
!   List of cosets
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTStarOfK
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

GTSGLeftCosets[group_, grpk_, lattice_] := Module[{el, pos, lc},
  lc = Union[Map[Sort[#] &, Table[el = GTSGgmat[group[[i]], grpk[[j]], lattice];
      pos = Position[group[[;; , 1]], el[[1]]];
      group[[pos[[1, 1]]]], {i, 1, Length[group]}, {j, 1, Length[grpk]}]]];
  Table[SortBy[lc[[i]], GTSGOrderOfElement[#, lattice] &], {i, 1, Length[lc]}]
  ]


(*
***)


(****z* /GTSGInduceOrbitLengthOne
! NAME
!  GTSGInduceOrbitLengthOne
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   11/02/2018: first version
! USAGE
!   GTSGInduceOrbitLengthOne induces two (three) irreducible representations from an irreducible representation of orbit length 1
! INPUT
!   factor group G/T, subgroup of index n, coset representative (G/T = H + q*H (+q^2*H)), permutation between initial G/T and G/T given by coset decomposition,
!   the irreducible representation, real space basis
! OUTPUT
!   induced representations
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTSGgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTSGGetIreps for the induction of 2 irreducible representations from an irep with an orbit of length 1
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


GTSGInduceOrbitLengthOne[grp_, grpn_, qel_,qelsq_,qelqu_, perm_, irepsn_,irepsq_,basis_,kvec_,recbas_] :=
  Module[{n, gln, el, posbel, pos, Amat, Bmat, Cmat, Um,  Umsq, erg, ind, elf,expfacqqu, expfacq,dim,elfound},
  (*Calculate Index and Length of group*)
  n = Length[grp]/Length[grpn];
  gln = Length[grpn];  
  (*--------------------*)
  (*Generate the matrix U*)
  (*--------------------*)
  (*Pick an element of the original group, here, the last one and calculate the conjugate element with q*)
  dim = Length[irepsn[[1]]];
  posbel=0;
  elfound=False;
  While[Not[elfound]&&posbel<gln,
   posbel=posbel+1;
   elfound = Not[irepsn[[posbel]] === irepsq[[posbel]]] && Not[irepsn[[posbel]] === IdentityMatrix[dim]] && Not[irepsq[[posbel]] === IdentityMatrix[dim]]];
  (*Determine U from D(q b q^-1) = U D(b) U^-1*)
  
  Amat = irepsn[[posbel]];(* Exp[I ((qibelqf[[2]] - qibelq[[2]]).basis).(kvec.recbas)];*)
  Bmat = irepsq[[posbel]];
  (*Determine U*)
  pos = Position[grp,qelsq[[1]]][[1,1]];
  expfacq=Exp[-I ((qelsq[[2]] - grp[[pos,2]]).basis).(kvec.recbas)];
  pos = Position[grp,qelqu[[1]]][[1,1]];
  expfacqqu=Exp[-I ((qelqu[[2]] - grp[[pos,2]]).basis).(kvec.recbas)];
  Which[
   (*Index 2*)
   n == 2,
     Cmat=irepsn[[First@First@Position[grpn[[;;,1]], qelsq[[1]]]]]*expfacq;
     Um=GTGetUnitaryTrafoMatrix[Amat,Bmat,Cmat,2];,
  (*Index 3*)
   n==3,
     Cmat=irepsn[[First@First@Position[grpn[[;;,1]], qelqu[[1]]]]]*expfacqqu;
     Um=GTGetUnitaryTrafoMatrix[Amat,Bmat,Cmat,3];
  ];
  Umsq=Um.Um;
  Which[
   (*Index 2*)
   n == 2,
   ind = Table[
      Flatten[{irepsn, Table[
         el = GTSGgmat[qel, grpn[[T]], basis];
         elf = GTSGgmat[qel, grpn[[T]], basis, GOTakeMod -> False];
        Exp[i I \[Pi]]*Exp[I ((elf[[2]] - el[[2]]).basis).(kvec.recbas)]*Um.irepsn[[T]], {T, 1, gln}]}, 1], {i, 0, 1}];
   ,
   (*Index 3*)
   n == 3,
   ind = Table[
      Flatten[{irepsn, 
        Table[
      	        el  = GTSGgmat[qel, grpn[[T]], basis];
                elf = GTSGgmat[qel, grpn[[T]], basis, GOTakeMod -> False];
                Exp[-2 \[Pi] I i/3]*Exp[I ((elf[[2]] - el[[2]]).basis).(kvec.recbas)]*Um.irepsn[[T]], {T, 1, gln}],
        Table[
        	    el = GTSGgmat[qelsq, grpn[[T]], basis];
                elf = GTSGgmat[qelsq, grpn[[T]], basis, GOTakeMod -> False];
                Exp[-4 \[Pi] I i/3]*Exp[I ((elf[[2]] - el[[2]]).basis).(kvec.recbas)]*Umsq.irepsn[[T]], {T, 1, gln}]},         
       1], {i, 0, 2}];];
  erg = Table[Permute[ind[[i]], perm], {i, 1, n}];
  Return[erg]
  ]

(*
***)


(****z* /GTSGInduceOrbitLengthTwoThree
! NAME
!  GTSGInduceOrbitLengthOne
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   11/02/2018: first version
! USAGE
!   GTSGInduceOrbitLengthTwoThree induces an irreducible representations from an irreducible representation of orbit length 2 or 3
! INPUT
!   factor group G/T, subgroup of index n, coset representative q (G/T = H + q*H (+q^2*H)), q^2, q^3, permutation between initial G/T and G/T given by coset decomposition,
!   the irreducible representation, the irreducible representation conjugated with q, the irreducible representation conjugated with q^2, real space basis
!   kvector, reciprocal basis
! OUTPUT
!   induced representation
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTSGgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTSGGetIreps for the construction of a subgroup of index 2 or 3
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

GTSGInduceOrbitLengthTwoThree[grp_, grpn_, qel_, qelsq_, qelqu_,perm_, irep_, irepq_, irepqsq_, basis_, kvec_, recbas_] := 
 Module[{n, gln, grpnrep, pos, qmat, ind, rep, qmatsq, el, qelsqexp,qelquexp},
  n = Length[grp]/Length[grpn];
  gln = Length[grpn];
  Which[
   (*index 2*)
   n == 2,
   grpnrep = Table[KroneckerProduct[{{1, 0}, {0, 0}}, irep[[T]]] + KroneckerProduct[{{0, 0}, {0, 1}}, irepq[[T]]], {T, 1, gln}]; 
   pos = First@First@Position[grpn[[;;,1]], qelsq[[1]]];
   qelsqexp=Exp[-I ((qelsq[[2]] - grpn[[pos,2]]).basis).(kvec.recbas)];
   qmat = KroneckerProduct[{{0, 1}, {0, 0}}, IdentityMatrix[Length[irep[[pos]]]]] + KroneckerProduct[{{0, 0}, {1, 0}}, qelsqexp*irep[[pos]]];
   ind = Flatten[{grpnrep,
      Table[
       el = GTSGgmat[qel, grpn[[T]], basis, GOTakeMod -> False];
       pos = Position[grp[[;;,1]],el[[1]]][[1,1]];
       Exp[I ((el[[2]] - grp[[pos,2]]).basis).(kvec.recbas)]*qmat.grpnrep[[T]], {T, 1, gln}]}, 1],
   (*index 3*)
   n == 3,
   grpnrep = 
    Table[KroneckerProduct[{{1, 0, 0}, {0, 0, 0}, {0, 0, 0}}, irep[[T]]] + KroneckerProduct[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}, irepq[[T]]] + 
      KroneckerProduct[{{0, 0, 0}, {0, 0, 0}, {0, 0, 1}},irepqsq[[T]]], {T, 1, gln}];
   pos = First@First@Position[grpn[[;;,1]], qelqu[[1]]];
   qelquexp=Exp[-I ((qelqu[[2]] - grpn[[pos,2]]).basis).(kvec.recbas)];
   qmat = KroneckerProduct[{{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},IdentityMatrix[Length[irep[[pos]]]]] + KroneckerProduct[{{0, 0, 0}, {0, 0, 0}, {1, 0, 0}},qelquexp*irep[[pos]]];
   qmatsq = qmat.qmat // Simplify;
   ind = Flatten[{grpnrep,
      Table[
       el = GTSGgmat[qel, grpn[[T]], basis, GOTakeMod -> False];
       pos = Position[grp[[;;,1]],el[[1]]][[1,1]];
       Exp[I ((el[[2]] - grp[[pos,2]]).basis).(kvec.recbas)]*qmat.grpnrep[[T]], {T, 1, gln}],
      Table[
       el = GTSGgmat[qelsq, grpn[[T]], basis, GOTakeMod -> False];
       pos = Position[grp[[;;,1]],el[[1]]][[1,1]];
       Exp[I ((el[[2]] -grp[[pos,2]]).basis).(kvec.recbas)]*qmatsq.grpnrep[[T]], {T,1, gln}]}, 1]
   ];
  rep = Permute[ind, perm];
  Return[rep];
  ]
(*
***)



(****z* /GTSGOrderOfElement
! NAME
!  GTSGOrderOfElement
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   11/02/2018: first version
! USAGE
!   GTSGOrderOfElement gives the order of a space group element
! INPUT
!   space group element, real space basis
! OUTPUT
!   order of a space group element
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTSGgmat
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTSGCosetRepresentative for the construction of an element q with the poperty G = H + q*H (+q^2*H)
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

GTSGOrderOfElement[el_, basis_] := Module[{ord, stop, eln},
  ord = 0;
  stop = False;
  eln = el;
  While[Not[stop],
   eln = GTSGgmat[el, eln, basis];
   stop = eln === el;
   ord = ord + 1];
  Return[ord]
  ]

(*
***)

(****z* /GTSGPrintGeneral
! NAME
!  GTSGPrintGeneral
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   12/06/2018: first version
! USAGE
!   Prints the k-vector and the group of k
! INPUT
!   group of k, k-vector, reciprocal basis
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
!  needed by GTSGGetIreps for the final output messages
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

GTSGPrintGeneral[grpk_,kvec_,recbas_,star_]:=Module[{Pr},
Print["Specified k-vector (Cartesian): ",kvec.recbas];
Print["Specified k-vector (Direct)   : ",kvec];
Print["Group of the k-vector:"];
Print[grpk];
If[Length[star]==1,
Print["Star of the k-vector: (1 arm)"];
Print[star];
,
Print["Star of the k-vector: ("<>ToString[Length[star]]<>" arms)"];
Print[star];
];
]

(*
***)

(****z* /GTSGPrintIreps
! NAME
!  GTSGPrintIreps
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   12/06/2018: first version
! USAGE
!   Prints the representation matrices of a space group at a specific k-point
! INPUT
!   group of k, ireps
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
!  needed by GTSGGetIreps for the final output messages
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
GTSGPrintIreps[grpk_,ireps_,genireps_,grp_]:=Module[{Pr},
Print["Representation matrices of the group of k:"];
Pr=Table[Table[MatrixForm[ireps[[h,T]]],{T,1,Length[grpk]}],{h,1,Length[ireps]}];
Pr=Prepend[Pr,grpk];
Pr=Table[Prepend[Pr[[i]],Superscript["\[CapitalGamma]", ToString[i-1]]],{i,1,Length[Pr]}];
Pr[[1,1]]="";
Print[Grid[Pr,Frame->All,ItemSize->Full,Dividers->{{2->GTDividerColor1},{2->GTDividerColor1}},Background->{{1->GTBackGroundColor1},{1->GTBackGroundColor1},{1,1}->GTCornerColor}]];
If[Length[genireps]>0,
Print["Representation matrices of the space group:"];
Pr=Table[Table[MatrixForm[genireps[[h,T]]],{T,1,Length[grp]}],{h,1,Length[ireps]}];
Pr=Prepend[Pr,grp];
Pr=Table[Prepend[Pr[[i]],Superscript["\[CapitalGamma]", ToString[i-1]]],{i,1,Length[Pr]}];
Pr[[1,1]]="";
Print[Grid[Pr,Frame->All,ItemSize->Full,Dividers->{{2->GTDividerColor1},{2->GTDividerColor1}},Background->{{1->GTBackGroundColor1},{1->GTBackGroundColor1},{1,1}->GTCornerColor}]];
]
]
(*
***)
(****z* /GTSGPrintChars
! NAME
!  GTSGPrintChars
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   12/06/2018: first version
! USAGE
!   Prints the character table of a space group at a specific k-point
! INPUT
!   group of k, classes, characters
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
!  needed by GTSGGetIreps for the final output messages
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

GTSGPrintChars[grpk_,classes_,chars_,genireptrue_]:=Module[{Pr,labels},
If[genireptrue,
Print["Classes of G/T:"],
Print["Classes of G(k)/T:"]];
Table[Print[Subscript["C", i], " = ",classes[[i]]] ,{i,1,Length[classes]}];
Pr=chars;
labels=Table[Length[classes[[i]]]classes[[i,1]],{i,1,Length[classes]}];
Pr=Prepend[Pr,labels];
Pr=Table[Prepend[Pr[[i]],Superscript["\[CapitalGamma]", ToString[i-1]]],{i,1,Length[Pr]}];
Pr[[1,1]]="";
Print[Grid[Pr,Frame->All,ItemSize->Full,Dividers->{{2->GTDividerColor1},{2->GTDividerColor1}},Background->{{1->GTBackGroundColor1},{1->GTBackGroundColor1},{1,1}->GTCornerColor}]];
]

(*
***)

(****z* /GTSGRightCosets
! NAME
!  GTSGRightCosets
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   11/02/2018: first version
!   09/09/2020: enforce ordering of the cosets with respect to the order of the element to avoid nontrivial coset representatives (taken as first element of the list) for trivial coset decompositions G = E G
! USAGE
!   GTSGRightCosets gives the right cosets of the factor group G/T (G being a space group, T the group of pure lattice translations) with respect to a subgroup
! INPUT
!   two factor groups and a real space basis
! OUTPUT
!   List of cosets
! GTPack OPTIONS
!   
! STANDARD OPTIONS
!
! GTPack MODULES
!   GTgmat
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

GTSGRightCosets[group_, grpk_, lattice_] := Module[{el, pos, lc},
  lc = Union[Map[Sort[#] &, Table[el = GTSGgmat[grpk[[j]],group[[i]], lattice];
        pos = Position[group[[;; , 1]], el[[1]]];
        group[[pos[[1, 1]]]], {i, 1, Length[group]}, {j, 1, Length[grpk]}]]];
  Table[SortBy[lc[[i]], GTSGOrderOfElement[#, lattice] &], {i, 1, Length[lc]}]
  ]


(*
***)

(****z* /GTStarOfK
! NAME
!  GTStarOfK
! AUTHOR
!  M. Geilhufe, W. Hergert
! PACKAGE
!   Lattice.m
! MODIFICATION HISTORY
!  * 29.01.2014 - first version
!  * 19.08.2016 - Star of K was not correct. See Cornwell for definition.
!  * 23.08.2016 - Simplify command, allow for 2D lattices, GTGetRotationMatrix instead of switching standard representation
!  * 24.05.2018 - check header an documentation
!  * 03.07.2018 - change towards double groups
!  * 12.12.2018 - moved to RepresentationTheorySG, made compatible with space group elements
! USAGE
!  GTStarOfK[group,kvect,recbas] gives the star of the wave vector k for a given group.
! INPUT
!   * group    -  the point group
!   * kvec     - the wavevector  
!   * recbas   - reciprocal basis
! OUTPUT
!   star of the wave vector
! GTPack OPTIONS
!   GOFast \[Rule] not working right now!
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  * GTGroupQ
!  * GTGetRotationMatrix
!  * GTGroupOfK
!  * GTLeftCosets
!  * GTSGLeftCosets
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  Needs to be extended to double groups
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTStarOfK[group_, kvec0_, recbas_, OptionsPattern[]] := Module[{i,gstar,star,grpk,leftcos,kvec,kvt,list,sol,dc,lattice},
   (*Check input - GOFast*)  
   (*If[OptionValue[GOFast],None,
   	If[GTGroupQ[group0],
         None,
         Print[$ErrorNoGroup]; Abort[]
      ];
   ];*)
  lattice=GTReciprocalBasis[recbas];
  If[Not[Length[kvec0]==Length[lattice]],Print[$ErrorDimensionKvec];Abort[]];
  (*Calculate Group of K and left cosets with respect to full group*)
  grpk=GTGroupOfK[group,kvec0,recbas,GOFast->True];
  If[GTWhichInput[group[[1]]]<5,
       leftcos = GTLeftCosets[group,grpk,GOFast->True],
       leftcos =GTSGLeftCosets[group,grpk,lattice]
   ];
   (*Apply one coset representative each to the k-vector -> Star of K*)
   gstar = leftcos[[;;,1]];
   kvec=If[Length[kvec0]==3,kvec0,Append[kvec0,0]];
   star = Table[
            kvt=GTGetRotationMatrix[gstar[[i]]].kvec;
            list = Sum[c[i] recbas[[i]], {i, 1, Length[recbas]}];
            sol = Solve[Table[list[[i]] == kvt[[i]], {i, 1, Length[recbas]}],Table[c[i], {i, 1, Length[recbas]}]];
            dc = Table[c[i], {i, 1, Length[recbas]}] /. sol[[1]];
            Mod[dc, 1].recbas
          , {i, 1, Length[gstar]}] // Union;
   If[OptionValue[GOLeftCosets],Return[{star[[;;,1;;Length[kvec0]]],leftcos}],Return[star[[;;,1;;Length[kvec0]]]]];
]

(*
***)

(* Attributes[GTStarOfK]={Protected, ReadProtected}*)


End[]


EndPackage[]
