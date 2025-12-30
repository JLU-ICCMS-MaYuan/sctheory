(****m* /CrystalField.m
!
! NAME
!  CrystalField.m
! AUTHOR
!  M. Geilhufe
! MODIFICATION HISTORY
!  1/27/16 : initially created and documented  
! USAGE
!  all modules connected to crystal field theory
!
! ERROR MESSAGES
!  
! GTPack MODULES
!
! GTBSTOperator calculates the matrix elements of the Buckmaster-Smith-Thornley operator equivalent
! GTBSTOperatorElement calculates a single matrix element of the Buckmaster-Smith-Thornley operator equivalent
! GTCFDatabaseInfo gives information about the crystal field expectation value parameter sets available in a database
! GTCFDatabaseRetrieve loads a crystal field expectation value parameter set from a given database
! GTCFDatabaseUpdate reads the database, adds a new parameter set and stores the database
! GTCrystalField gives the crystal field Hamiltonian for a certain group up to a maximal angular momentum quantum number lmax
! GTCrystalFieldParameter calculates the crystal field parameters using the point charge model
! GTCrystalFieldSplitting calculates the decomposition of irreducible representations if the symmetry is lowered from group1 to group2
! GTStevensOperator calculates the matrix elements of the Stevens operator equivalent
! GTStevensOperatorElement calculates a single matrix element of a Stevens operator equivalent
! GTStevensTheta gives the specified numerical prefactor for operator equivalents within the crystal field expansion
!
! GTPack NOTEBOOKS 
!  none in the moment
!
! DESCRIPTION
!  CrystalField.m contains all modules connected to crystal field theory.
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



BeginPackage["GroupTheory`CrystalField`",{"GroupTheory`Basic`","GroupTheory`RepresentationTheory`","GroupTheory`Symbols`","GroupTheory`Auxiliary`","GroupTheory`AngularMomentum`"}]

 GTBSTOperator  	      ::usage = "GTBSTOperator[\*StyleBox[\"l,m,J\",\"TI\"]] calculates the matrix elements \*Cell[BoxData[FormBox[RowBox[{\"\[LeftAngleBracket]\",RowBox[{\"J\",\",\",RowBox[{SubscriptBox[\"m\",\"1\"],\"\[VerticalSeparator]\",SubsuperscriptBox[OverscriptBox[\"O\",\"^\"],\"l\",\"m\"],\"\[VerticalSeparator]\",\"J\"}],\",\",SubscriptBox[\"m\",\"2\"]}],\"\[RightAngleBracket]\"}],TraditionalForm]],\"InlineMath\"], where \*SubsuperscriptBox[OverscriptBox[\"O\",\"^\"],\"l\",\"m\"] is the Buckmaster-Smith-Thornley operator equivalent."
 GTBSTOperatorElement     ::usage = "GTBSTOperatorElement[\*StyleBox[\"l,m,J,m1,m2\",\"TI\"]] calculates a single matrix element \*Cell[BoxData[FormBox[RowBox[{\"\[LeftAngleBracket]\",RowBox[{\"J\",\",\",RowBox[{SubscriptBox[\"m\",\"1\"],\"\[VerticalSeparator]\",SubsuperscriptBox[OverscriptBox[\"O\",\"^\"],\"l\",\"m\"],\"\[VerticalSeparator]\",\"J\"}],\",\",SubscriptBox[\"m\",\"2\"]}],\"\[RightAngleBracket]\"}],TraditionalForm]],\"InlineMath\"], where \*SubsuperscriptBox[OverscriptBox[\"O\",\"^\"],\"l\",\"m\"] is the Buckmaster-Smith-Thornley operator equivalent."
 GTCFDatabaseInfo         ::usage = "GTCFDatabaseInfo[\*StyleBox[\"database\",\"TI\"]] gives information about a parameter set with radial expectation values available in a \*StyleBox[\"database\",\"TI\"] and to be used within a crystal field expansion."
 GTCFDatabaseRetrieve     ::usage = "GTCFDatabaseRetrieve[\*StyleBox[\"database,parameter set\",\"TI\"]] loads a \*StyleBox[\"parameter set\",\"TI\"] containing radial expectation values from a given \*StyleBox[\"database\",\"TI\"]."
 GTCFDatabaseUpdate       ::usage = "GTCFDatabaseUpdate[\*StyleBox[\"database\",\"TI\"]] reads the \*StyleBox[\"database\",\"TI\"], adds a new parameter set and stores the \*StyleBox[\"database\",\"TI\"]."
 GTCrystalField           ::usage = "GTCrystalField[\*StyleBox[\"group,lmax\",\"TI\"]] gives the crystal field Hamiltonian for a certain \*StyleBox[\"group\",\"TI\"] up to a maximal angular momentum quantum number \*StyleBox[\"lmax\",\"TI\"]."
 GTCrystalFieldParameter  ::usage = "GTCrystalFieldParameter[\*StyleBox[\"l,m,positions,charges\",\"TI\"]] calculates the crystal field parameters \*Cell[BoxData[FormBox[SubsuperscriptBox[\"A\",\"l\",\"m\"],TraditionalForm]],\"InlineMath\"] using the point charge model."
 GTCrystalFieldSplitting  ::usage = "GTCrystalFieldSplitting[\*StyleBox[\"group1, group2, ct1, ct2\",\"TI\"]] calculates the decomposition of irreducible representations if the symmetry is lowered from \*StyleBox[\"group1\", \"TI\"] to \*StyleBox[\"group2\", \"TI\"]."
 GTStevensOperator	      ::usage = "GTStevensOperator[\*StyleBox[\"l,m,J\",\"TI\"]] calculates the matrix elements \*Cell[BoxData[FormBox[RowBox[{\"\[LeftAngleBracket]\",RowBox[{\"J\",\",\",RowBox[{SubscriptBox[\"m\",\"1\"],\"\[VerticalSeparator]\",SubsuperscriptBox[OverscriptBox[\"O\",\"^\"],\"l\",\"m\"],\"\[VerticalSeparator]\",\"J\"}],\",\",SubscriptBox[\"m\",\"2\"]}],\"\[RightAngleBracket]\"}],TraditionalForm]],\"InlineMath\"], where \*SubsuperscriptBox[OverscriptBox[\"O\",\"^\"],\"l\",\"m\"] is the Stevens operator equivalent."
 GTStevensOperatorElement ::usage = "GTStevensOperatorElement[\*StyleBox[\"l,m,J,m1,m2\",\"TI\"]] calculates a single matrix element \*Cell[BoxData[FormBox[RowBox[{\"\[LeftAngleBracket]\",RowBox[{\"J\",\",\",RowBox[{SubscriptBox[\"m\",\"1\"],\"\[VerticalSeparator]\",SubsuperscriptBox[OverscriptBox[\"O\",\"^\"],\"l\",\"m\"],\"\[VerticalSeparator]\",\"J\"}],\",\",SubscriptBox[\"m\",\"2\"]}],\"\[RightAngleBracket]\"}],TraditionalForm]],\"InlineMath\"], where \*SubsuperscriptBox[OverscriptBox[\"O\",\"^\"],\"l\",\"m\"] is the Stevens operator equivalent."
 GTStevensTheta           ::usage = "GTStevensTheta[\*StyleBox[\"element,l\",\"TI\"]] gives the specified numerical prefactor for operator equivalents within the crystal field expansion."

 Options[GTCrystalField]             = {GOHarmonics->"Complex", GOFast->True}
 Options[GTCrystalFieldParameter]    = {GOHarmonics->"Complex"}
 Options[GTCrystalFieldSplitting]  = {GOFast -> GOFastValue, GOIrepNotation -> "Mulliken"}
 
Begin["`Private`"] (* Begin Private Context *) 

(****d* /GTBSTOperator
! NAME
!  GTBSTOperator
! AUTHOR
!  M. Geilhufe, W. Hergert
! PACKAGE
!  test.m -> Auxiliary.m
! MODIFICATION HISTORY
!   02/18/2016: first version
! USAGE
!   GTBSTOperator[l,m,J] calculates the matrix elements \[LeftAngleBracket]J,Subscript[m, 1]\[VerticalSeparator]Subsuperscript[Overscript[O, ^], l, m]\[VerticalSeparator]J,Subscript[m, 2]\[RightAngleBracket], where Subsuperscript[Overscript[O, ^], l, m] is the so called Buckmaster-Smith-Thornley operator equivalent.
!   matrix form of dimension (2*lmax+1) x (2*lmax+1)
! INPUT
!   l,m, lmax
! OPTIONS
!  GOVerbose
! OUTPUT
!  
! GTPack OPTIONS
!  
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTBSTOperatorElement
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
 GTBSTOperator[k_, q_, j_] := Table[GTBSTOperatorElement[k, q, j, ms, m], {ms, -j, j}, {m, -j, j}]
  (*
***)

(****d* /GTBSTOperatorElement
! NAME
!  GTBSTOperatorElement
! AUTHOR
!  M. Geilhufe, W. Hergert
! PACKAGE
!  test.m -> Auxiliary.m
! MODIFICATION HISTORY
!   02/18/2016: first version
! USAGE
!   GTBSTOperatorElement[l,m,J,m1,m2] calculates a single matrix element \[LeftAngleBracket]J,Subscript[m, 1]\[VerticalSeparator]Subsuperscript[Overscript[O, ^], l, m]\[VerticalSeparator]J,Subscript[m, 2]\[RightAngleBracket], where Subsuperscript[Overscript[O, ^], l, m] is the so called Buckmaster-Smith-Thornley operator equivalent.
!   matrix form of dimension (2*lmax+1) x (2*lmax+1)
! INPUT
!   l,m, lmax
! OUTPUT
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
 GTBSTRedMatElement[k_, j_] := Sqrt[(2*j + k + 1)!/(2*j - k)!]/2^k
 GTBSTOperatorElement[k_, q_, j_, ms_,m_] := (-1)^(j - ms) Quiet[ThreeJSymbol[{j, -ms}, {k, q}, {j, m}]] GTBSTRedMatElement[k, j]
  (*
***)
(****d* /GTCFDatabaseInfo
! NAME
!  GTTbDatabaseInfo
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  CrystalField.m 
! MODIFICATION HISTORY
!  * February 2016 : 1st version 
! USAGE
!  GTCFDatabaseInfo[database] gives information about the crystal field expectation value parameter sets available in database.
! INPUT
!  database  - name of the database  
! OUTPUT
!  database and printout of the information
! GTPack OPTIONS
!  
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTReadFromFile
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
! - 
! RELEASE
!
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)
 GTCFDatabaseInfo[db_] := Module[{dbname, data, tab},
  dbname = StringJoin[db, ".parm"];
  data = GTReadFromFile[dbname];
  tab = Table[Flatten[{data[[i, 1 ;; 2]], r^2, r^4, r^6, r^-3}] /. data[[i, 3]], {i, 1, Length[data]}];
  tab = Prepend[tab, {"", "Reference", "<\!\(\*SuperscriptBox[\(r\), \(2\)]\)>", "<\!\(\*SuperscriptBox[\(r\), \(4\)]\)>", "<\!\(\*SuperscriptBox[\(r\), \(6\)]\)>", "<\!\(\*SuperscriptBox[\(r\), \(-3\)]\)>"}];
  Print[Grid[Sort@tab, Frame -> All, Dividers -> {{2 ->GTDividerColor1}, {2 -> GTDividerColor1}}, 
                 Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} ->GTCornerColor}]];
  ]
  (*
***)

(****d* /GTCFDatabaseRetrieve
! NAME
!  GTCFDatabaseRetrieve
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  CrystalField.m 
! MODIFICATION HISTORY
!  * February 2016 : 1st version 
! USAGE
!  GTCFDatabaseRetrieve[db, name] loads a crystal field expectation value parameter set from a given database.
! INPUT
!  * db   - name of the database
!  * name - name of the parameter set
! OUTPUT
!  parameter set
! GTPack OPTIONS
!  
! STANDARD OPTIONS
!
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)
GTCFDatabaseRetrieve[db_, name_] := Module[{dbp, dbname, nset, pos}, 
  dbname = StringJoin[db, ".parm"];
  (*--- Check if database and parameters exist ---*)
  If[Not[FileExistsQ[dbname]], Print["Error: the database " <> dbname <>" was not found in the working directory!"]; Abort[], None];
  (*--- retrieve the parameter set If[Position[dbp,name]\[Equal]{},Print[name," -> parameters not in database"];Abort[],None]; ---*)
  dbp = ReadList[dbname][[1]];
  nset = Length[dbp];
  pos = Position[dbp, name];
  If[Length[pos] == 0, Print["Error: The requested parameter set is not found in the file " <> dbname <> "."]; Abort[], Return[dbp[[First@First@pos, 3]]]]]
  (*
***)

(****d* /GTCFDatabaseUpdate
! NAME
!  GTCFDatabaseUpdate
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  CrystalField.m 
! MODIFICATION HISTORY
!  * February 2016 : 1st version 
! USAGE
!  GTCFDatabaseUpdate[db, name] reads the database, adds a new parameter set and stores the database.
! INPUT
!  * db   - name of the database
!  * name - name of the parameter set
! OUTPUT
!  parameter set
! GTPack OPTIONS
!  
! STANDARD OPTIONS
!
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The input of data is organized interactively, using the routine INPUT 
!  of Mathematica. The input requires the following steps:
!     Input                 Type                Example 
!     -----------------------------------------------------------------------------------------------  
!     Name of the entry     String              "Ho_SIC"   
!     Reference             String              "J. Forstreuter, et al., Phys. Rev. B, 55, 15, 1997"
!     r^2                   Real                 0.756
!     r^4                   Real                 1.43
!     r^6                   Real                 5.63
!     r^-3                  Real                10.7
!
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)
GTCFDatabaseUpdate[db_] := Module[{dbname, dbp, name, src, r2, r4, r6, r3},
  dbname = StringJoin[db, ".parm"];
  (*--- Check if database exists ---*)
  If[Flatten[Position[FileNames[], dbname]] == {}, Print[dbname, " -> database not in working directory"]; Abort[],  None];
  (*--- Read database ---*)
  dbp = ReadList[dbname][[1]];
  (*--- input of general information ---*)
  name = ToString@Input["Name of the database entry (String)"];
  src = ToString@Input["Reference (String)"];
  r2 = Input["Parameter: r^2"];
  r4 = Input["Parameter: r^4"];
  r6 = Input["Parameter: r^6"];
  r3 = Input["Parameter: r^-3"];
  dbp = Append[dbp, {name, src, {r^2 -> r2, r^4 -> r4, r^6 -> r6, r^-3 -> r3}}];
  Put[dbp, dbname];
  Return[dbp]]
  (*
***)

(****d* /GTCrystalField
! NAME
!  GTCrystalField
! AUTHOR
!  M. Geilhufe
! PACKAGE
!   
! MODIFICATION HISTORY
!  02/05/2021: new version
!  12/28/2013 : first version
! USAGE
!  GTCrystalField[group,lmax] gives the crystal field Hamiltonian for a certain group up to a maximal angular momentum quantum number lmax.
! INPUT
!  group, lmax
! OUTPUT
!  
! GTPack OPTIONS
!  GOHarmonics, GOFast
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
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTCrystalField[grp_, lmax_,OptionsPattern[]] := 
 Module[{gg, cproj, alm, ylm, ca, solution},
  If[Not[OptionValue[GOFast]],If[GTGroupQ[grp],None,Print[$ErrorNoGroup];Return[]]];
  (*--- Length of the group ---*)
  gg = Length[grp];
  Sum[
  (*--- Coefficient vectors ---*)
  alm = Table[A[ll, m], {m, -ll, ll}];
  If[OptionValue[GOHarmonics]=="Complex",
  			ylm = Table[Y[ll, m], {m, -ll, ll}],
  			ylm = Table[S[ll, m], {m, -ll, ll}]];
  (*--- Character projection operator for AngularMomentum representation ---*)
  cproj = Simplify[Total[GTAngularMomentumRep[grp, ll, GOFast -> True, GOHarmonics->OptionValue[GOHarmonics]]]]/Length[grp];
  (*--- Transformation of parameters ---*)
  ca = cproj.alm;
  (*--- allowed parameters ---*)
  solution = Flatten[Quiet[Solve[Table[alm[[i]] == ca[[i]], {i, 1, Length[ca]}], alm]]];
  r^ll alm.ylm /. solution, {ll, 0, lmax,2}]
  ] 
  
  (*
***)



(****d* /GTCrystalFieldSplitting
! NAME
!  GTCrystalFieldSplitting
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  CrystalField.m
! MODIFICATION HISTORY
!   January 2016: first version
!   April 2016: table changed, last head "splitting in group 2" also with yellow background
! USAGE
!   calculates the decomposition of irreducible representations if the symmetry is lowered from group1 to group2.
! INPUT
!   rotation matrix
! OPTIONS

! OUTPUT
!  index
! GTPack OPTIONS
!  GOFast, GOIrepNotation
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTSubGroupQ, GTCharacterTable, GTIrep
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  calculates the decomposition of irreducible representations if the symmetry is lowered from group1 to group2.
! LITERATURE
!  -
! TODO
! -
! RELEASE
!
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTCrystalFieldSplitting[grp1_, grp2_, ct1in_: {}, ct2in_: {},OptionsPattern[]] := 
 Module[{rep,gf,ct1, ct2, cl1, cl2, postable, charsys, ireps, names, np, tb,tb2, nrop, l, tb3,tabf},
  If[Not[OptionValue[GOFast]],If[Not[GTSubGroupQ[grp1,grp2]],Print[$ErrorGTCrystalFieldSubGroup];Abort[]]];
  rep=OptionValue[GOIrepNotation];
  gf = OptionValue[GOFast];
  ct1 = If[Length[ct1in] == 0, GTCharacterTable[grp1, GOIrepNotation -> rep, GOVerbose -> False, GOFast -> gf], ct1in];
  ct2 = If[Length[ct2in] == 0, GTCharacterTable[grp2, GOIrepNotation -> rep, GOVerbose -> False, GOFast -> gf], ct2in];
  
  cl1 = ct1[[1]];
  cl2 = ct2[[1]];
  postable = Table[pos = First@First@Position[cl1, cl2[[i, 1]]], {i, 1,Length[cl2]}];
  charsys =  Table[Table[ct1[[2, i, postable[[j]]]], {j, 1, Length[cl2]}], {i,1, Length[cl1]}];
  ireps =    Table[GTIrep[charsys[[i]], ct2, GOVerbose -> False], {i, 1,Length[cl1]}];
  names = Table[
    np = ireps[[j]];
    tb = Flatten[Table[Which[np[[i]] == 0, "a", np[[i]] == 1, {ct2[[3, i]], "\[CirclePlus]"}, np[[i]] > 1, {np[[i]] ct2[[3, i]], "\[CirclePlus]"}], {i, 1,Length[cl2]}]];
    tb2 = Delete[tb, Position[tb, "a"]];
    nrop = Length[Position[tb2, "\[CirclePlus]"]];
    l = Length[tb2] - Length[Position[tb2, "\[CirclePlus]"]];
    tb3 = If[l == nrop,Delete[tb2, Position[tb2, "\[CirclePlus]"][[nrop]]], None], {j,1, Length[ireps]}];
    tabf=Prepend[Table[Append[Prepend[charsys[[i]], ct1[[3, i]]],Row[names[[i]]]], {i, 1, Length[charsys]}],Flatten@{"Repr. in group 1", cl2[[;; , 1]], "Splitting in group 2"}];
    Grid[tabf,Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> 
    GTDividerColor1}}, Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {{1, 1} -> GTCornerColor, {1, Length[tabf] + 1} -> GTCornerColor}}]
  ]

(*
***)

(****d* /GTCrystalFieldParameter
! NAME
!  GTCrystalFieldParameter
! AUTHOR
!  M. Geilhufe
! PACKAGE
!   
! MODIFICATION HISTORY
!  06/02/2014 : first version
!  02/19/2016 : New version
! USAGE
!  GTCrystalFieldParameter[l,m,positions,charges] calculates the crystal field parameters A^m_l using the point charge model.
! INPUT
!  l,m, list with vectors, list with charges
! OUTPUT
!  parameter
! GTPack OPTIONS
!  GOHarmonics
! GTPack MODULES
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  GTTesseralHarmonicY
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

GTCrystalFieldParameter[l_, m_, Rlist_, qlist_,OptionsPattern[]] := Module[{R, th, ph,Ylm,i}, 
	1/(2*l + 1)*Sum[
       R = Rlist[[i]];
       th = ArcTan[R[[3]], Norm[{R[[1]], R[[2]]}]];
       If[Tan[th] == 0 || Norm[{R[[1]], R[[2]]}] == 0, ph = 0, ph = ArcTan[R[[1]], R[[2]]]];
       Which[	OptionValue[GOHarmonics]=="Real",
  					Ylm=GTTesseralHarmonicY[l, m, th, ph],
  		        OptionValue[GOHarmonics]=="Complex",
  			        Ylm=Conjugate[SphericalHarmonicY[l, m, th, ph]]];
       qlist[[i]]*Ylm/Norm[R]^(l + 1)
   ,{i, 1, Length[Rlist]}]]
 (*
***)

(****d* /GTStevensOperator
! NAME
!  GTStevensOperator
! AUTHOR
!  M. Geilhufe, W. Hergert
! PACKAGE
!  test.m -> Auxiliary.m
! MODIFICATION HISTORY
!   04/22/2014: first version
! USAGE
!   GTStevensOperator[l,m,J] calculates the matrix elements \[LeftAngleBracket]J,Subscript[m, 1]\[VerticalSeparator]Subsuperscript[Overscript[O, ^], l, m]\[VerticalSeparator]J,Subscript[m, 2]\[RightAngleBracket], where Subsuperscript[Overscript[O, ^], l, m] is the so called Stevens operator equivalent.
!   matrix form of dimension (2*lmax+1) x (2*lmax+1)
! INPUT
!   l,m, lmax
! OPTIONS
!  GOVerbose
! OUTPUT
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
!  not fully tested
!  documentation
! PROBLEMS
!
! 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTPow[a_, b_] := If[Abs[a] + Abs[b] == 0, 1, a^b];
GTfn[k_, qi_, j_, M1_, M2_, m_, s_] := 
 KroneckerDelta[M1, M2 + qi]*
   Sqrt[((j - M2)! (j + M2 + qi)!)/((j + M2)! (j - M2 - qi)!)] + 
  s*(-1)^(k - qi - m) KroneckerDelta[M1, M2 - qi]*
   Sqrt[((j + M2)! (j - M2 + qi)!)/((j - M2)! (j + M2 - qi)!)]

GTStevensOperator[k_, qii_, j_] := Module[{qi, sq, a, f1, f2, s, flist, Fkq, afkq, alpha, x},
  qi = Abs[qii];
  If[qii == 0, sq = 1, sq = Sign[qii]];
  (*--- according to Ryabov (Journal of Magnetic Resonance140, 141-145 (1999)) we define the coefficients a ---*)
  a = Table[Table[0, {k + 1}], {k + 1}];
  a[[k + 1, 1]] = 1;
  Do[q = q1 + 1;
   Do[Clear[f1, f2];
    If[q + m - 1 <= k && m >= 1, f1 = a[[q + 1, m]], f1 = 0];
    If[q + m <= k && m >= 0, f2 = a[[q + 1, m + 1]], f2 = 0];
    s = (2*q + m - 1)*f1 + (q (q - 1) - m (m + 1)/2)*f2;
    s = s + 
      If[k - q - m >= 1, 
       Sum[(-1)^
          n (Binomial[m + n, n]*x - Binomial[m + n, m - 1] - 
           Binomial[m + n, m - 2])*
         If[q + m + n <= k && m + n > 0, a[[q + 1, m + n + 1]], 
          0], {n, 1, k - q - m}], 0];
    a[[q1 + 1, m + 1]] = s, {m, 0, k - q1}], {q1, k - 1, 0, -1}];
  (*--- we need to normalize these coefficients by finding the greatest common divisor F[q,g] ---*)
  flist = Map[First,Map[First, Map[FactorList, Factor[a], {2}], {2}], {2}];
  Fkq = Table[GCD @@ flist[[i]], {i, 1, k + 1}];
  afkq = Factor[Table[If[Fkq[[i]] > 0, a[[i]]/Fkq[[i]], 0], {i, 1, k + 1}]];
  (*--- furthermore the parameter \[Alpha] is introduced. alpha= 1 for all q if k is an odd integer and alpha=1/2 and alpha=1 for odd and even q. ---*)
  alpha = If[OddQ[k], 1, If[OddQ[qi], 1/2, 1]];
  If[sq < 0,
   Table[alpha*Sum[GTfn[k, qi, j, M1, M2, m, sq]*afkq[[qi + 1, m + 1]]*GTPow[M2, m], {m, 0, k - qi}]/(2*I) /. x -> j*(j + 1), {M1, -j,j}, {M2, -j, j}],
   Table[alpha*Sum[GTfn[k, qi, j, M1, M2, m, sq]*afkq[[qi + 1, m + 1]]*GTPow[M2, m], {m, 0, k - qi}]/2 /. x -> j*(j + 1), {M1, -j, j}, {M2, -j, j}]]]

(*
***)

(****d* /GTStevensOperatorElement
! NAME
!  GTStevensOperator
! AUTHOR
!  M. Geilhufe, W. Hergert
! PACKAGE
!  test.m -> Auxiliary.m
! MODIFICATION HISTORY
!   04/22/2014: first version
! USAGE
!   GTStevensOperator[l,m,J,m1,m2] calculates a single matrix element \[LeftAngleBracket]J,Subscript[m, 1]\[VerticalSeparator]Subsuperscript[Overscript[O, ^], l, m]\[VerticalSeparator]J,Subscript[m, 2]\[RightAngleBracket], where Subsuperscript[Overscript[O, ^], l, m] is the so called Stevens operator equivalent.
!   matrix form of dimension (2*lmax+1) x (2*lmax+1)
! INPUT
!   l,m, lmax
! OUTPUT
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
!  not fully tested
!  documentation
! PROBLEMS
!
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTStevensOperatorElement[k_, qii_, j_,M1_ ,M2_] := Module[{qi, sq, a, f1, f2, s, flist, Fkq, afkq, alpha, x},
  qi = Abs[qii];
  If[qii == 0, sq = 1, sq = Sign[qii]];
  (*--- according to Ryabov (Journal of Magnetic Resonance140, 141-145 (1999)) we define the coefficients a ---*)
  a = Table[Table[0, {k + 1}], {k + 1}];
  a[[k + 1, 1]] = 1;
  Do[q = q1 + 1;
   Do[Clear[f1, f2];
    If[q + m - 1 <= k && m >= 1, f1 = a[[q + 1, m]], f1 = 0];
    If[q + m <= k && m >= 0, f2 = a[[q + 1, m + 1]], f2 = 0];
    s = (2*q + m - 1)*f1 + (q (q - 1) - m (m + 1)/2)*f2;
    s = s + 
      If[k - q - m >= 1, 
       Sum[(-1)^
          n (Binomial[m + n, n]*x - Binomial[m + n, m - 1] - 
           Binomial[m + n, m - 2])*
         If[q + m + n <= k && m + n > 0, a[[q + 1, m + n + 1]], 
          0], {n, 1, k - q - m}], 0];
    a[[q1 + 1, m + 1]] = s, {m, 0, k - q1}], {q1, k - 1, 0, -1}];
  (*--- we need to normalize these coefficients by finding the greatest common divisor F[q,g] ---*)
  flist = Map[First,Map[First, Map[FactorList, Factor[a], {2}], {2}], {2}];
  Fkq = Table[GCD @@ flist[[i]], {i, 1, k + 1}];
  afkq = Factor[Table[If[Fkq[[i]] > 0, a[[i]]/Fkq[[i]], 0], {i, 1, k + 1}]];
  (*--- furthermore the parameter \[Alpha] is introduced. alpha= 1 for all q if k is an odd integer and alpha=1/2 and alpha=1 for odd and even q. ---*)
  alpha = If[OddQ[k], 1, If[OddQ[qi], 1/2, 1]];
  If[sq < 0,
   alpha*Sum[GTfn[k, qi, j, M1, M2, m, sq]*afkq[[qi + 1, m + 1]]*GTPow[M2, m], {m, 0, k - qi}]/(2*I) /. x -> j*(j + 1),
   alpha*Sum[GTfn[k, qi, j, M1, M2, m, sq]*afkq[[qi + 1, m + 1]]*GTPow[M2, m], {m, 0, k - qi}]/2 /. x -> j*(j + 1)]]

(*
***)

(****d* /GTStevensTheta
! NAME
!  GTStevensTheta
! AUTHOR
!  M. Geilhufe, W. Hergert
! PACKAGE
!  
! MODIFICATION HISTORY
!   02/19/2016: first version
! USAGE
!   GTStevensTheta[element,l] gives the specified numerical prefactor for operator equivalents within the crystal field expansion. 
! INPUT
!   l,m, lmax
! OPTIONS
!  GOVerbose
! OUTPUT
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
!  not fully tested
!  documentation
! PROBLEMS
!
! 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTStevensTheta[element_, l_] := Module[{list, pos, erglist},
  list = {{"Ce", -2/35, 2/(7 45),0}, 
  	      {"Pr", -52/(11 15^2), -4/(55 33 3),17 16/(7 11^2 13 5 3^4)}, 
  	      {"Nd", -7/(33 33), -8 17/(11 11 13 297), -17 19 5/(13^2 11^3 3^3 7)},
  	      {"Pm", 14/(11 11 15),952/(13 3^3 11^3 5), 2584/(11^2 13^2 3 63)}, 
  	      {"Sm", 13/(7 45),26/(33 7 45), 0}, 
  	      {"Eu", 0, 0, 0}, 
  	      {"Gd", 0, 0, 0},
  	      {"Tb", -1/(99), 2/(11 1485), -1/(13 33 2079)},
  	      {"Dy", -2/(9 35), -8/(11 45 273),  4/(11^2 13^2 3^3 7)},
  	      {"Ho", -1/(30 15), -1/(11 2730), -5/(13 33 9009)}, 
  	      {"Er", 4/(45 35), 2/(11 15 273), 8/(13^2 11^2 3^3 7)},
  	      {"Tu", 1/99, 8/(3 11 1485), -5/(13 33 2079)}, 
  	      {"Yb", 2/63, -2/(77 15),4/(13 33 63)}};
  pos = Position[list, ToString[element]];
  erglist = If[Length[pos] > 0, list[[First@First@pos]], Print[$ErrorGTStevensTheta]; Abort[]];
  Which[l == 0, Return[1], l == 2, Return[erglist[[2]]], l == 4, Return[erglist[[3]]], l == 6, Return[erglist[[4]]]];
  Print[$ErrorGTStevensThetaLval];
 ]
(*
***)
 Attributes[GTBSTOperator]={Protected, ReadProtected}
 Attributes[GTBSTOperatorElement]={Protected, ReadProtected}
 Attributes[GTCFDatabaseInfo]={Protected, ReadProtected}
 Attributes[GTCFDatabaseUpdate]={Protected, ReadProtected}
 Attributes[GTCFDatabaseRetrieve]={Protected, ReadProtected}
 Attributes[GTCrystalField]={Protected, ReadProtected}
 Attributes[GTCrystalFieldParameter]={Protected, ReadProtected}
 Attributes[GTCrystalFieldSplitting]={Protected, ReadProtected}
 Attributes[GTStevensOperator]={Protected, ReadProtected}
End[]



(*/*-------------------------------------------- Commands for GTPack2.0 ---------------------------------------------*)


(*---------------- Usage ---------------*) 
 

(*---------------- Options -------------*)

(*---------------- Error messages ------*)



Begin["`Private`"] (* Begin Private Context  GTPack2.0*) 

(****y* 
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