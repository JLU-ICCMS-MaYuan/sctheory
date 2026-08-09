(****p* /Test.m
!
! NAME
!  Test.m
! AUTHOR
!  W. Hergert, M. Geilhufe
! MODIFICATION HISTORY
! 10/26/12 : initially created and documented  
! USAGE
!  Moduls under test , not incorporated in the production packages
!
! ERROR MESSAGES
!  
! GTPack MODULES
!  GTGenerators           - Find the generators from a group
!  GTGroupFromGenerators  - Constructs the group from generators
!
! GTPack NOTEBOOKS 
!  Wolfram GT_Basic.nb
! DESCRIPTION
!   The following modules have been implemented somewhere else:
!
!         GTSymmetryElementQ      -> CrystalStructure.m
!         GTIcosaederAxes[]       -> Install.m
!         GTGetIrepMatrix         -> Auxiliary.m
!         GTSymmetryProductChars  -> RepresentationTheory.m
!         GTVectorRep             -> RepresentationTheeory.m
!
!
! LITERATURE
! 
! TODO
!  My (Wolfram) moduls are perhaps far from being elegant. A reformulation could help.
!
! PROBLEMS
!
!--------------------------------------------------------------------------------
!
***)

BeginPackage["GroupTheory`Test`",{"GroupTheory`RepresentationTheory`","GroupTheory`Basic`","GroupTheory`Symbols`","GroupTheory`Install`","GroupTheory`Auxiliary`","GroupTheory`Lattice`","GroupTheory`TightBinding`",
	"GroupTheory`CrystalStructure`","GroupTheory`PseudoPotential`","GroupTheory`Photonics`","GroupTheory`ElectronicStructure`"}]
 
 
  GTReadAbinitBands           ::usage=  "GTReadAbinitBands[SCF-File, Bands-File] reads k-mesh, energy bands and Fermi-energy from Abinit-Output" 
  GTAbinitConvertWaveFunction ::usage=  "GTAbinitConvertWaveFunction[binary wave function file, k-mesh, number of bands, lattice] transforms an Abinit wavefunction and writes wave function files for further analysis."
  GTTbConvertWaveFunction     ::usage=  "GTTbConvertWaveFunction[wav-file, k-mesh, number of bands, basis, ll, geometry] transforms a Tight Binding wavefunction and writes wave function files for further analysis."
  GTAnalyzeBandStructure      ::usage=  "GTAnalyzeBandStructure[g0, lattice, kpt, mbands] analyzes the wave funcions corresponding to a band structure calculated using Abinit. "
  GTReadWaveFunction          ::usage=  "GTReadWaveFunction[kp,band] reads a converted wave function."
  GTTransformFunction         ::usage=  "GTTransformFunction[sym,f,{x,y,z}] transforms the scalar function f[x,y,z] with respect to the symmetry elements sym."
  GTCharacterTablePrint       ::usage=  "GTCharacterTablePrint[group] prints a character table if group elments are given only in matrix form."
  GTSymmetryElementQ          ::usage=  "GTSymmetryElementQ[sym,cluster] checks if sym is symmetry element of cluster." 
  GTIcosahedronAxes           ::usage=  "GTIcosahedronAxes[] installs axes for icosahedron group."
  GTInstallColorGroup         ::usage=  "GTInstallColorGroup[group, invariant subgroup] installs Shubnikov groups of first or second kind." 
  GTSymmetrizedProductChars   ::usage=  "GTSymmetrizedProductChars[character table, irepn, n , symmetry] gives the characters of the nth power of the irepn-th IREP of a group with character table."
  GTExportVASP                ::usage=  "GTExportVASP[structure] exports a structure to a POSCAR file."
  GTCharacters                ::usage=  "GTCharacters[matrices] calculates characters of a representation matrix or a list of represtation matrices."
  GTVectorRep                 ::usage=  "GTVectorRep[character table] gives the vector representation for a group with character table."
  GTGetIrepMatrix           ::usage=  "GTPrintIrepMatrix[group,irep,sym] prints matrices of irreducible representations in a nice form." 
  GTIrepMatrices              ::usage=  "GTIrepMatrices[index] gives the hard wired Irep matrices for Oh according to Cornwell for G1,G12,G25' and G15."
   
   	
Options[GTAnalyzeBandStructure]   = {GOMethod->"TightBinding", GOMaxPoints->10,GOIrepNotation -> "Mulliken"}
Options[GTTransformFunction]      = {GOVerbose ->True}
Options[GTCharacterTablePrint]    = {GOFast -> False, GOVerbose -> True}
Options[GTSymmetryElementQ]       = {GOVerbose -> False,GOPrecision->10^(-7)}
Options[GTInstallColorGroup]      = {GOVerbose -> False,GOFast->False}
Options[GTIcosahedronAxes]        = {GOVerbose -> True}
Options[GTInstallColorGroup]      = {GOVerbose -> False,GOFast->True}
Options[GTSymmetrizedProductChars]= {GOVerbose->True}
Options[GTCharacters]             = {GOClasses -> False}
Options[GTVectorRep]              = {GOVerbose -> False}
Options[GTGetIrepMatrix]        = {GOVerbose -> True}

Begin["`Private`"] 


(****g* /GTIrepMatrices
! NAME
!  GTIrepMatrices
! AUTHOR
!  W. Hergert
! PACKAGE
!   test.m 
! MODIFICATION HISTORY
!  08/08/2016 : first version
! USAGE
!  GTIrepMatrices[index] gives the hard wired Irep matrices for Oh according to Cornwell for G1,G12,G25' and G15.
! INPUT
!  
! OUTPUT
!  vector representation
! ERROR MESSAGES
!  -
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
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTIrepMatrices[index_] := Module[{oh, c2, c3, rep, sym, pos},
  oh = GTInstallGroup[Oh, GOVerbose -> False];
  c2 = {{{Ee, C2x, C2y, C2z, IEe, IC2x, IC2y, IC2z}, {{1, 0}, {0, 1}}},
        {{C3\[Alpha], C3\[Beta], C3\[Gamma], C3\[Delta], 
      IC3\[Alpha], IC3\[Beta], IC3\[Gamma], 
      IC3\[Delta]}, {{-1/2, -Sqrt[3]/2}, {Sqrt[3]/2, -1/2}}},
           {{C3\[Alpha]i, C3\[Beta]i, C3\[Gamma]i, C3\[Delta]i, 
      IC3\[Alpha]i, IC3\[Beta]i, IC3\[Gamma]i, 
      IC3\[Delta]i}, {{-1/2, Sqrt[3]/2}, {-Sqrt[3]/2, -1/2}}},
            {{C4x, C4xi, C2e, C2f, IC4x, IC4xi, IC2e, 
      IC2f}, {{1/2, -Sqrt[3]/2}, {-Sqrt[3]/2, -1/2}}},
    {{C4y, C4yi, C2c, C2d, IC4y, IC4yi, IC2c, 
      IC2d}, {{1/2, Sqrt[3]/2}, {Sqrt[3]/2, -1/2}}},
    {{C4z, C4zi, C2a, C2b, IC4z, IC4zi, IC2a, 
      IC2b}, {{-1, 0}, {0, 1}}}};
  c3 = {{{Ee, IEe}, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}},
           {{C3\[Alpha], 
      IC3\[Alpha]}, {{0, -1, 0}, {0, 0, 1}, {-1, 0, 0}}},
           {{C3\[Beta], 
      IC3\[Beta]}, {{0, 1, 0}, {0, 0, -1}, {-1, 0, 0}}},
           {{C3\[Gamma], 
      IC3\[Gamma]}, {{0, -1, 0}, {0, 0, -1}, {1, 0, 0}}},
           {{C3\[Delta], 
      IC3\[Delta]}, {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}}, 
           {{C3\[Alpha]i, 
      IC3\[Alpha]i}, {{0, 0, -1}, {-1, 0, 0}, {0, 1, 0}}},
            {{C3\[Beta]i, 
      IC3\[Beta]i}, {{0, 0, -1}, {1, 0, 0}, {0, -1, 0}}},
     {{C3\[Gamma]i, 
      IC3\[Gamma]i}, {{0, 0, 1}, {-1, 0, 0}, {0, -1, 0}}},
      {{C3\[Delta]i, IC3\[Delta]i}, {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}}},
    {{C2x, IC2x}, {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}}},
    {{C2y, IC2y}, {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}},
    {{C2z, IC2z}, {{1, 0, 0}, {0, -1, 0}, {0, 0, -1}}},
    {{C4x, IC4x}, {{0, 0, 1}, {0, -1, 0}, {-1, 0, 0}}},
    {{C4y, IC4y}, {{0, -1, 0}, {1, 0, 0}, {0, 0, -1}}},
    {{C4z, IC4z}, {{-1, 0, 0}, {0, 0, -1}, {0, 1, 0}}},
    {{C4xi, IC4xi}, {{0, 0, -1}, {0, -1, 0}, {1, 0, 0}}},
    {{C4yi, IC4yi}, {{0, 1, 0}, {-1, 0, 0}, {0, 0, -1}}},
    {{C4zi, IC4zi}, {{-1, 0, 0}, {0, 0, 1}, {0, -1, 0}}},
    {{C2a, IC2a}, {{1, 0, 0}, {0, 0, -1}, {0, -1, 0}}},
    {{C2b, IC2b}, {{1, 0, 0}, {0, 0, 1}, {0, 1, 0}}},
    {{C2c, IC2c}, {{0, -1, 0}, {-1, 0, 0}, {0, 0, 1}}},
    {{C2d, IC2d}, {{0, 1, 0}, {1, 0, 0}, {0, 0, 1}}},
    {{C2e, IC2e}, {{0, 0, -1}, {0, 1, 0}, {-1, 0, 0}}},
    {{C2f, IC2f}, {{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}}
           };
  (*--- selection of matrices ---*)
    If[index == 0,
       rep = Table[{{1}}, {48}],
       None
    ];
    If[index == 1,
       rep = {};
       Do[sym = oh[[i]];
          pos = Flatten[Position[c2, sym]][[1]];
          rep = Append[rep, c2[[pos, 2]]]
       , {i, 1, 48}],
        None
    ];
    If[index == 2,
       rep = {};
       Do[sym = oh[[i]];
          pos = Flatten[Position[c3, sym]][[1]];
          rep = Append[rep, c3[[pos, 2]]]
       , {i, 1, 48}],
       None
    ];
    If[index == 3,
       rep = GTGetMatrix[#] & /@ oh,
       None
    ];
    Return[rep]
]

(*
***)


(****g* /GTGetIrepMatrix
! NAME
!  GTGetIrepMatrix
! AUTHOR
!  W. Hergert
! PACKAGE
!   test.m 
! MODIFICATION HISTORY
!  08/08/2016 : first version
!  23.02.2017 : changes from Notebook in Return included
! USAGE
!  GTGetIrepMatrix[group,irep,sym] prints matrices of irreducible representations in a nice form."
! INPUT
!  
! OUTPUT
!  vector representation
! ERROR MESSAGES
!  -
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
!
!--------------------------------------------------------------------------------
! SOURCE
*)


GTGetIrepMatrix[grp_, irep_, sym_,OptionsPattern[]] := Module[{pos,lg,ld,imat,tab,tb,verb},
  (*--- options ---*)
   verb=OptionValue[GOVerbose];	
   If[Head[sym] === List,
      lg   = Length[sym]; 
      ld   = IntegerPart[lg/10];
      imat = {};
      Do[
         pos = Flatten[Position[grp, sym[[i]]]][[1]];
        imat = Append[imat, irep[[pos]]]
      , {i, 1, lg}];
      tab = Join[{sym}, {MatrixForm[#] & /@ imat}] // Transpose;
      If[verb,
         Do[
            tb = Take[tab, {(i - 1)*10 + 1, i*10}] // Transpose;
            Print[ Grid[tb, Frame -> All, 
                        Dividers   -> {Black, {2 -> GTDividerColor1}},
                        Background -> {None, {1 -> GTBackGroundColor1}}
                       ]
                 ]
         , {i, 1, ld}];
         tb = Take[tab, {ld*10 + 1, lg}] // Transpose;
         Print[Grid[tb, Frame -> All, 
                    Dividers -> {Black, {2 -> GTDividerColor1}}, 
                    Background -> {None, {1 -> GTBackGroundColor1}}
                   ]          
              ],
         None
      ];
      Return[ Join[{sym},  {imat}]//Transpose],     
      pos = Flatten[Position[grp, sym]][[1]];
      tab = {{sym}, {MatrixForm[irep[[pos]]]}};
      If[verb,
         Print[Grid[tab, Frame -> All, 
                    Dividers -> {Black, {2 -> GTDividerColor1}}, 
                    Background -> {None, {1 -> GTBackGroundColor1}}
                   ]
              ],
         None
      ];     
      Return[{{sym,irep[[pos]]}}]
   ];  
]


(*
***)


(****g* /GTVectorRep
! NAME
!  GTVectorRep
! AUTHOR
!  W. Hergert
! PACKAGE
!   test.m -> RepresentationTheory.m
! MODIFICATION HISTORY
!  08/07/2016 : first version
! USAGE
!  GTVectorReps[character table] calculates the vector representation to a group with character table.
! INPUT
!  
! OUTPUT
!  vector representation
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTGetMatrix
!  GTCharacters
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

(****g* /GTCharacters
! NAME
!  GTCharacters
! AUTHOR
!  W. Hergert
! PACKAGE
!   test.m -> RepresentationTheory.m
! MODIFICATION HISTORY
!  08/01/2016 : first version
! USAGE
!  GTCharacters[matrices] calculates the characters of a representation matrix or a list of represtation matrices.
! INPUT
!  
! OUTPUT
!  charcters
! ERROR MESSAGES
!  -
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
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTCharacters[elements_, OptionsPattern[]] := Module[{cls,elm1,cc,lc,elm2,i,chars},
   cls = OptionValue[GOClasses];
   elm1=elements;
   If[Head[elements] === List,
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

(****g* /GTTransformFunction
! NAME
!  GTReadAbinitBands
! AUTHOR
!  W. Hergert
! PACKAGE
!   test.m -> RepresentationTheory.m
! MODIFICATION HISTORY
!  05/13/2013 : first version
! USAGE
!  Transformation of a scalar function
! INPUT
!  sym - a single symmetry element or alist of them
!  f   - function to transform
!  xyz - lsit of arguments of f
! OUTPUT
!  transformed functions and a table
! ERROR MESSAGES
!  -
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  GTTransFormFunction.nb in Reference_Symbols
! DESCRIPTION
!  
! LITERATURE
!  
! TODO
!  may be it can be extended to space group elements or vector fields
! PROBLEMS
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTransformFunction[sym_, f_, xyz_,OptionsPattern[]] := Module[{ln,tf,arg,i,verb},
  verb = OptionValue[GOVerbose];
  If[Head[sym] === List,
     ln = Length[sym]; tf = {};
     Do[
        arg = Inverse[GTGetMatrix[sym[[i]]]].xyz;    
        tf = Append[tf, f[arg[[1]], arg[[2]], arg[[3]]]]
     , {i, 1, ln}];tf=tf//Expand;
     If[verb,
        Print[
     	  Grid[{sym, tf}, Frame -> All,Alignment -> Center,
          Dividers -> {2 -> Black, {2 -> GTDividerColor1}}, 
          Background -> {None, {1 -> GTBackGroundColor1}}]
        ]
     ];Return[tf],
     arg = Inverse[GTGetMatrix[sym]].xyz;
     Return[ f[arg[[1]], arg[[2]], arg[[3]]]//Expand]
   ]
]


(*
***)

(****g* /GTReadAbinitBands
! NAME
!  GTReadAbinitBands
! AUTHOR
!  S. Thomas, M. Geilhufe
! PACKAGE
!   test.m -> ElectronicStructure.m
! MODIFICATION HISTORY
!  10/22/2013 : first version
! USAGE
!  reads k-mesh, energy bands and Fermi-energy from Abinit-Output
! INPUT
!  Eigenvalue file from self-consisten calculation, Eigenvalue file from band structure calculation
! OUTPUT
!  k-mesh, energy bands, Fermi-energy
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
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTReadAbinitBands[fermifile_, eigfile_] := Module[{fermi, bdat, nkpt, bands, band, kpt, nl},
  fermi = ToExpression[StringCases[ReadList[fermifile, Record][[1]],RegularExpression["[^\\s]+\\.\\d+"]][[1]]];
  bdat = ReadList[eigfile, Record];
  nkpt = StringCases[bdat[[1]], RegularExpression["[0-9]+"]][[1]] // ToExpression;
  bands = ConstantArray[{}, nkpt];
  band = 0;
  kpt = {};
  For[nl = 2, nl <= Length[bdat], nl++,
   If[Length[StringPosition[bdat[[nl]], "kpt"]] != 0,
   			AppendTo[kpt,StringCases[bdat[[nl]], RegularExpression[".\\.\\d+"]][[2 ;; 4]] // ToExpression]; band++,
   			
    		bands[[band]] = bands[[band]]~Join~ToExpression[StringCases[bdat[[nl]], RegularExpression["[^\\s]+\\.\\d+"]]]
    ]
   ];
  bands -= fermi;
  Return[{kpt, bands, fermi}]]

(*
***) 

(****g* /GTConvertAbinitWaveFunction
! NAME
!  GTConvertAbinitWaveFunction
! AUTHOR
!  S. Thomas, M. Geilhufe
! PACKAGE
!   test.m -> ElectronicStructure.m
! MODIFICATION HISTORY
!  10/22/2013 : first version
! USAGE
!  reads wave functions from Abinit output and produces readable files.
! INPUT
!  binary wave function file, k-mesh, number of bands, lattice vectors
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
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTAbinitConvertWaveFunction[wfkfile_, kpt_, mbands_, lattice_] := 
 Module[{spin,file, nkpt, k1, k2, k3, band, kp, infile, indat, dat, dim, f,tdat, a, b, c, aa, bb, cc, r1, r2, r3, na, nb, nc, printout, s,nx, ny, nz, res, cmin, cmax, cd, k,x,y,z},
  (*--- Parameters are given explicitely up to now ---*)
  spin = 1;
  res = 11;
  cmin = -1/2;
  cmax = 1/2;
  cd = (cmax - cmin)/(res - 1);
  {r1, r2, r3} = lattice;
  {k1, k2, k3} = GTReciprocalBasis[lattice]/(2*\[Pi]);
  nkpt = Length[kpt];
  file = wfkfile;
  Do[
  	k = kpt[[kp, 1]]*k1 + kpt[[kp, 2]]*k2 + kpt[[kp, 3]]*k3;
  	Do[
  		NotebookDelete[printout]; 
  		printout = PrintTemporary["k-Point: " <> ToString[kp] <> " / " <> ToString[nkpt] <> " Band: " <> ToString[band] <> " / " <> ToString[mbands]];
  		(*If[nkpt==1,	
  			OutputForm[file <> "\n1\n0\n" <> ToString[band] <> "\n0\n0\n4\nwave\n0\n"] >> "!cut3d";,
  			OutputForm[file <> "\n1\n0\n" <> ToString[kp] <> "\n" <> ToString[band] <> "\n0\n0\n4\nwave\n0\n"] >> "!cut3d";];*)
  		RunProcess[$SystemShell, "StandardOutput", "\n(echo " <> file <> "; echo 1; echo 0; echo " <> ToString[kp] <>"; echo " <> ToString[band] <>"; echo 0; echo 0; echo 1; echo wave; echo 0;) | cut3d\nexit\n"];
  		infile = "wave_k" <> ToString[kp] <> "_b" <> ToString[band] <> "_s1";
		indat = ReadList[infile, Array[Number &, 5]];
		(*<< ("!rm " <> infile);*)
		dat = Map[#[[4]] + I #[[5]] &, indat];
		dim = Length[dat]^(1/3);
		dat = Transpose[Partition[dat, dim^2]];
		dat = Transpose[Partition[dat, dim]];
		f = ListInterpolation[dat, {{0, 1 - 1/dim}, {0, 1 - 1/dim}, {0, 1 - 1/dim}}];
		tdat = Table[{a, b, c} = {aa, bb, cc} /.Flatten[Solve[aa r1 + bb r2 + cc r3 == {x, y, z}]];
   		{na, nb, nc} = Floor[{a, b, c}];
   		{a, b, c} = {a, b, c} - {na, nb, nc};(*r=na r1+nb r2+ nc r3;*)
   		{x, y, z,Quiet[f[a, b, c]*Exp[I*2*\[Pi]*k.(a r1 + b r2 + c r3)]*Exp[I*2*\[Pi]*k.(na r1 + nb r2 + nc r3)],InterpolatingFunction::dmval]}, {x, cmin, cmax, cd}, {y, cmin,cmax, cd}, {z, cmin, cmax, cd}];
		s = OpenWrite["wave_" <> ToString[kp] <> "_" <> ToString[band]];
		WriteString[s,ToString[spin] <> "\n" <> ToString[res] <> "\n" <> ToString[res] <> "\n" <> ToString[res] <> "\n"];
		For[nx = 1, nx <= res, nx++,
  			For[ny = 1, ny <= res, ny++,
   				For[nz = 1, nz <= res, nz++,
    				WriteString[s,ToString[tdat[[nx, ny, nz, 4]], InputForm] <> "\n"]
    			]
   			]
  		];
  		Close[s]
  	,{band,1,mbands}]  	
  	,{kp,1,nkpt}];
	NotebookDelete[printout];]
(*
***) 


(****g* /GTAnalyzeBandStructure
! NAME
!  GTAnalyzeBandStructure
! AUTHOR
!  S. Thomas, M. Geilhufe
! PACKAGE
!   test.m -> ElectronicStructure.m
! MODIFICATION HISTORY
!  10/22/2013 : first version
!  12/16/2013 : faster version
! USAGE
!  analyzes the wave funcions corresponding to a band structure calculated using Abinit. 
! INPUT
!  point group, lattice vectors, k-mesh, number of bands
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
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTAnalyzeBandStructure[grp_, lattice_, kpt0_, mbands_,OptionsPattern[]] := 
 Module[{kpoint,nkpt, rlat, GoKList, Groups, ng, GoKPos, CTGroups, v, p,CTProj, xval, yval, zval, ft, ftt, g, sh, irep, printoutk,irepname, spin,out,q, ct, clline, printtable, ireps, irtable,cl,kpt,kptnames},
  kpt = kpt0[[1]];
  kptnames = kpt0[[2]];
  NotebookDelete[printoutk]; 
  printoutk = PrintTemporary["Group of k-points"];
  nkpt = Length[kpt];
  Which[	OptionValue[GOMethod]=="Abinit",
  							rlat = GTReciprocalBasis[lattice]/(2 \[Pi]);
  							GoKList = Table[GTGetMatrix@GTGroupOfK[grp, kpt[[i]].rlat, rlat, GOFast -> True], {i, 1,nkpt}];,
  			OptionValue[GOMethod]=="TightBinding",
  							rlat = GTReciprocalBasis[lattice];
  							GoKList = Table[GTGetMatrix@GTGroupOfK[grp, kpt[[i]]*\[Pi], rlat, GOFast -> True], {i, 1,nkpt}];];
  Groups = Intersection[GoKList, GoKList];
  ng = Length[Groups];
  GoKPos = Table[Flatten[Position[Groups, GoKList[[i]]]][[1]], {i, 1, nkpt}];
  NotebookDelete[printoutk]; 
  printoutk = PrintTemporary["Calculation of character tables"];
  CTGroups = Table[GTCharacterTable[Groups[[i]], GOIrepNotation -> OptionValue[GOIrepNotation],GOVerbose -> "False", GOFast -> True], {i, 1, ng}];
  Which[grpdgrp==1,  v[x_, y_, z_] = p[x, y, z];spin=1,
  		grpdgrp==2,  v[x_, y_, z_] = {p[x, y, z],q[x,y,z]};spin=2,
  		grpdgrp==3,  v[x_, y_, z_] = p[x, y, z];spin=1,
  		grpdgrp==5,  v[x_, y_, z_] = {p[x, y, z],q[x,y,z]};spin=2];
  NotebookDelete[printoutk]; 
  printoutk = PrintTemporary["Application of the character projection operator"];
  CTProj = Table[GTCharProjectionOperator[CTGroups[[i, 1]], CTGroups[[i, 2]],v, {x, y, z}, GOFast -> True], {i, 1, ng}]//N;
  If[spin==1,
  (*--- Without Spin ---*)
  out=Table[
  	NotebookDelete[printoutk];
   	printoutk = 
    PrintTemporary["Analyze bands, k-Point: " <> ToString[kp] <> "/" <> ToString[nkpt]];
   	Table[
    	{{xval, yval, zval}, ft} = GTReadWaveFunction[kp, band];
    	(*g[x_, y_, z_] = CTProj[[GoKPos[[kp]]]] /. p -> ft;*)
    	sh = Table[
    		g[x_, y_, z_] = CTProj[[GoKPos[[kp]], i]] /. p -> ft;
      			Quiet[NIntegrate[Abs[g[x, y, z] - ft[x, y, z]], {x, -xval, xval}, {y, -yval,yval}, {z, -zval, zval},Method -> {"MonteCarlo", "MaxPoints" -> OptionValue[GOMaxPoints]}]]
      			,{i, 1, Length[CTGroups[[GoKPos[[kp]], 1]]]}];
		irep = First[Flatten[Position[sh, Min[sh]]]];
    	irepname = CTGroups[[GoKPos[[kp]]]][[3, irep]]
    ,{band, 1,mbands}]
  ,{kp, 1, nkpt}];,
  
  (*--- Spin ---*)
   out= Table[
  	NotebookDelete[printoutk];
   	printoutk = 
    PrintTemporary["Analyze bands, k-Point: " <> ToString[kp] <> "/" <> ToString[nkpt]];
   	Table[
    	{{xval, yval, zval}, ft} = GTReadWaveFunction[kp, band];
    	ftt[x_, y_, z_] = {ft[[1]][x, y, z], ft[[2]][x, y, z]};
    	(*g[x_, y_, z_] = CTProj[[GoKPos[[kp]]]] /. p -> ft;*)
    	sh = Table[
    		g[x_, y_, z_] = CTProj[[GoKPos[[kp]], i]] /. p -> ft[[1]]/. q -> ft[[2]];
      			Norm[Quiet[NIntegrate[Abs[g[x, y, z] - ftt[x, y, z]], {x, -xval, xval}, {y, -yval,yval}, {z, -zval, zval},Method -> {"MonteCarlo", "MaxPoints" -> OptionValue[GOMaxPoints]}]]]
      			,{i, 1, Length[CTGroups[[GoKPos[[kp]], 1]]]}];
		irep = First[Flatten[Position[sh, Min[sh]]]];
    	irepname = CTGroups[[GoKPos[[kp]]]][[3, irep]]
    ,{band, 1,mbands}]
  ,{kp, 1, nkpt}];];
  (*--- Print Output ---*)
  Do[
  ct = CTGroups[[GoKPos[[kp]]]];
  Print["--------------------------------------------------------------------------"];
  kpoint = If[OptionValue[GOMethod]=="TightBinding",kpt[[kp]],kpt[[kp]].rlat];
  Print["k-point:     "<>ToString[kptnames[[kp]]]<>", ", kpoint];
  Print["Point group:"];
  Print[GTGetSymbol@GoKList[[kp]]];
  (*Print["Classes:"];*)
  cl=GTGetSymbol@ct[[1]];
  (*Print[cl];*)
  ireps=Union[out[[kp]]];
  Print["Symmetry"];
  (*clline = Prepend[Table[Subscript["C", ToString[i]], {i, 1, Length[ct[[1]]]}], ""];*)
  clline = Prepend[Table[Length[cl[[i]]] cl[[i,1]], {i, 1, Length[ct[[1]]]}], ""];
  Do[ Print["Band ",i,":",out[[kp,i]]],{i,1,mbands}];  
  irtable=Table[  		
  	Prepend[ct[[2, First@Flatten[Position[ct[[3]], ireps[[i]]]]]],ireps[[i]]]
  	,{i,1,Length[ireps]}];
  printtable=Prepend[irtable,clline];
  Print[DisplayForm[GridBox[printtable, RowLines -> {1, 0}, ColumnLines -> {1, 0}]]];
    ,{kp,1,nkpt}];  
  Return[{kpt,kptnames,out}]
  ]
(*
***) 

(****g* /GTReadWaveFunction
! NAME
!  GTReadWaveFunction
! AUTHOR
!  M. Geilhufe
! PACKAGE
!   test.m -> ElectronicStructure.m
! MODIFICATION HISTORY
!  12/16/2013 : first version
! USAGE
!  read a wave funcions corresponding to a band structure. The wave function has to be converted using GTTbConvertWaveFunction or GTConvertAbinitWaveFunction.
! INPUT
!  k-point (index), band
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
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTReadWaveFunction[kp_, band_] := Module[{filename, rawdat, xval, yval, zval, f, data, spin, dataup, datadown, fup, fdown},
  filename = "wave_" <> ToString[kp] <> "_" <> ToString[band];
  rawdat = ReadList[filename, Expression];
  {spin,xval, yval, zval} = rawdat[[1 ;; 4]];
  data = rawdat[[5 ;;]];
  (*--- Convert data to 3d array ---*)
  data = Partition[data, zval];
  data = Partition[data, yval];
  If[spin==1,
  	f = ListInterpolation[data, {{-xval, xval}, {-yval, yval}, {-zval, zval}}];
  	,
  	dataup   = Table[data[[x,y,z,1]],{x,1,xval},{y,1,xval},{z,1,xval}];
  	datadown = Table[data[[x,y,z,2]],{x,1,xval},{y,1,xval},{z,1,xval}];
  	fup = ListInterpolation[dataup, {{-xval, xval}, {-yval, yval}, {-zval, zval}}];
  	fdown = ListInterpolation[datadown, {{-xval, xval}, {-yval, yval}, {-zval, zval}}];
  	f = {fup,fdown};
  	];
  Return[{{xval, yval, zval}, f}]]
(*
***) 

(****g* /GTTbConvertWaveFunction
! NAME
!  GTTbConvertWaveFunction
! AUTHOR
!  M. Geilhufe
! PACKAGE
!   test.m -> TightBinding.m
! MODIFICATION HISTORY
!  12/16/2013 : first version
! USAGE
!  produces readable files for GTAnalyizeBandStructure.
! INPUT
!  wav-file, k-mesh, number of bands, basis, ll, geometry
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
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTTbConvertWaveFunction[file_, kpt_, mbands_, bas_, ll_, geo_,spin_:1] := 
 Module[{kp, band, s, nkpt, printout, field},
  nkpt = Length[kpt];
  Do[
   	Do[
    	NotebookDelete[printout];
        printout = PrintTemporary["k-Point: " <> ToString[kp] <> " / " <> ToString[nkpt] <> " Band: " <> ToString[band] <> " / " <> ToString[mbands]];
    	s = OpenWrite["wave_" <> ToString[kp] <> "_" <> ToString[band]];
        WriteString[s,ToString[spin] <> "\n" <> ToString[geo[[3]]] <> "\n" <> ToString[geo[[6]]] <> "\n" <>ToString[geo[[9]]] <> "\n"];
    	field = GTTbWaveFunction[file, kpt[[kp]], band, bas, ll, geo,spin];
    	Do[WriteString[s, ToString[field[[i]], InputForm] <> "\n"], {i, 1,Length[field]}];
    	Close[s]
    ,{band, 1, mbands}]
  ,{kp, 1, nkpt}]]

(*
***)   


(*
*)
(****g* /GTCharacterTableM
! NAME
!  GTCharacterTableM
! AUTHOR
!  W. Hergert
! PACKAGE
!  test.m -> RepresentationTheory.m
! MODIFICATION HISTORY
!   24/05/2016 : first version
! USAGE
!  GTCharacterTablePrint[group] prints a character table if group elments are given only in matrix form.
! INPUT
!  grp - group in matrix form or symbolic, symbolic will be converted in matrix form
!
! OUTPUT
! character table
! ERROR MESSAGES
!  
! GTPack MODULES
!  GTCharacterTable
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
*)item

GTCharacterTablePrint[grp_, OptionsPattern[]] := 
 Module[{fast,verb,group,classes,chars,names,nc,chtmp,chpr,i},
  (*--- options ---*)
  fast  = OptionValue[GOFast];
  verb  = OptionValue[GOVerbose];
  (*--- calculate character table --*) 
  group = Map[GTGetMatrix[#] &, grp];
  {classes, chars, names} ==  GTCharacterTable[group, GOFast -> fast, GOVerbose -> False];
  (*--- construct print ---*)
  nc    = Map[Length[#] &, classes];
  chtmp = Table[Prepend[chars[[i]], Superscript["\[CapitalGamma]", ToString[i]]], {i, 1, Length[chars]}];
  chpr  = Prepend[chtmp, 
  	               Flatten[{"",Table[
  	               	                 Subscript[ToString[nc[[i]]] <> "C", ToString[i]]
  	               	          , {i, 1, Length[chars]}]
  	                        }
  	               ]
  	      ];
  If[verb,
     Print[
     	   Grid[chpr, Frame -> All, 
                Dividers   -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
                Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} -> GTCornerColor}
           ]
      ],
      Return[{classes, chars, names}]
  ]
]


(*
***)


(****g* /GTSymmetryElementQ
! NAME
!  GTSymmetryElementQ
! AUTHOR
!  W. Hergert
! PACKAGE
!  test.m -> CrystalStructure.m
! MODIFICATION HISTORY
!   29/05/2016 : first version
! USAGE
!  GTSymmetryElementQ[symmetry,cluster] checks if the symmetry elements leave the cluster invariant.
! INPUT
!  symmetry - one symmetry element, or a list of symmetry elements
!  cluster  - ensemble of atoms to check with respect to symmetry 
! OUTPUT
!  logical value True or False
! ERROR MESSAGES
!  Compare: lists have to have same length
!  List of those symmetry elements which are not elements of the cluster.
! GTPack MODULES
!  GTGetMatrix
! GTPack NOTEBOOKS 
!  GTSymmetryElementQ in Wolfram_Devel/4_Discrete_Symmetry_Groups/New_Commands
! DESCRIPTION
!  The need for such a command appeared in connection with the description of the nanotubes. It is possible to show if 
!  the Dn groups are point groups of the nanotubes.
! LITERATURE
!
! TODO
!  more tests
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)


GTSymmetryElementQ[sym_, cluster_, OptionsPattern[]] := Module[{verb,ns,mats,clv,qtf,cl1,i,test,qtf1,log,pp,prec},
   (*--- options ---*)
   verb = OptionValue[GOVerbose];
   prec =OptionValue[GOPrecision];
   (*--- matrix form of symmmetry elements ---*)
   If[Head[sym] === List,
      ns = Length[sym]; mats = Map[GTGetMatrix[#] &, sym],
      ns = 1; mats = {GTGetMatrix[sym]}
   ];
   (*--- evaluate test ---*)
   clv = cluster;qtf = {};
   Do[
      cl1      = Transpose[clv]; 
      cl1[[1]] = Map[mats[[i]].# &, cl1[[1]]];
      cl1      = Transpose[cl1];
      test     = Compare[cl1, clv,GOPrecission->prec];
      If[test == {},
         qtf = Append[qtf, True],
         qtf = Append[qtf, False]
      ]
   , {i, 1, ns}];
   (*--- result of test ---*)
   qtf1 = Union[qtf];
   If[Length[qtf1] == 1 && qtf1[[1]] == True,
      log = True;pp={},
      log = False;
      If[verb,
         pp = Position[qtf, False] // Flatten;
         Print["Elements at positions ", pp," are not symmetry elements"]
      ]
   ]; 
   If[verb,
   	  Return[{log,pp}],
      Return[log]
   ]
]


Compare[cl1_, cl2_,OptionsPattern[]] := Module[{eps, cmp, cmp1, i, ll, lc, nn, dropc,vec, j}, 
	eps = 10^(-7);
    If[Length[cl1] == Length[cl2],
       ll = Length[cl1]; cmp = cl1; lc = ll,
       Print["Error: Lists do not have equal length"]; Abort[]
    ];
    Do[dropc = 0; cmp1 = cmp;
       Do[vec = cmp[[i, 1]] - cl2[[j, 1]];
          nn  = Norm[vec];
          If[nn <= eps && cmp[[i, 2]] == cl2[[j, 2]],
             cmp1 = Drop[cmp, {i}]; dropc = dropc + 1,
             None
          ]
       , {i, 1, lc}];
       lc = Length[cmp1];
       If[dropc > 1,
          Print["More than 1 element dropped  ", dropc]; Abort[],
          None
       ];
       cmp = cmp1
     , {j, 1, ll}];
     Return[cmp]
]

(*
***)





(****g* /GTSIcosahedronAxes
! NAME
!  GTIcosahedronAxes
! AUTHOR
!  M. Geilhufe, W. Hergert
! PACKAGE
!  test.m -> Install.m
! MODIFICATION HISTORY
!   06/08/2016 : first version
! USAGE
!  GTIcosahedronAxes[] installs axes for icosahedron group.
! INPUT
!  no input argument
! OUTPUT
!  no direct output, but all the axes are installed an can be used symbollically.
! OPTIONS
!  GOVerbose
! ERROR MESSAGES
!  Complains due to precision, but this does not matter. Is switched of oif GOVerbose=False
! GTPack MODULES
!  GTInstallAxes
! GTPack NOTEBOOKS 
!  GTIcosaeder_new.nb in Wolfram_Devel/4_Discrete_Symmetry_Groups/Nanotubes_Buckyballs.
! DESCRIPTION
! 
! LITERATURE
!  Altman, HErzig
! TODO
! 
! PROBLEMS
!  ! Two axes are only defined and not installed !!
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


(****g* /GTInstallColorGroup
! NAME
!  GTInstallColorgroups
! AUTHOR
!  M. Geilhufe, W. Hergert
! PACKAGE
!  test.m -> Install.m
! MODIFICATION HISTORY
!   06/27/2016 : first version
!   05/20/2021 : there was a typo in the line cmpl = Complement[grp1, 1 invs1]; (Before: "cmpl = Complement[grp, 1 invs1];")
! USAGE
!  GTIInstallColorGroup installs Shubnikov Groups of firsr and second kind
! INPUT
!  grp  - group used to construct the color group 
!  invs - {}   - Shubnikov group of first kind is installed
!         invs - invariant subgroup of index 2 of group  
! OUTPUT
!  the corresponding color group
! OPTIONS
!  GOVerbose
!  GOFast
! ERROR MESSAGES
!  If GOFast->False it will be checked if grp and invs are groups and if invs is an invariant subgroup of index 2 of 
!  group. If the checks are negative, the commant will stop with Abort[]
! GTPack MODULES
!  GTMagnetic, GTInstallGroup, GTGroupQ, GTInvSubGroupQ
! GTPack NOTEBOOKS 
!  GTIcosaeder_new.nb in Wolfram_Devel/4_Discrete_Symmetry_Groups/NMagnetic_Groups_Lattices/Colour_Groups
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

GTInstallColorGroup[grp_, invs_, OptionsPattern[]] := Module[{verb, fast, mgrp, cmpl,grp1,invs1},
 	(*---options---*)
     verb = OptionValue[GOVerbose];
     fast = OptionValue[GOFast];
    (*--- check magnetic ---*)
     If[GTMagneticQ[] == True,
        None,
        GTMagnetic[True];
        If[verb,          
           Print["Set internal variable to allow magnetic groups"], 
           None
        ]
     ];
    (*--- installation ---*)
     If[Head[grp] === Symbol,
        grp1 = GTInstallGroup[grp],
        grp1 = grp
     ];
     If[Head[invs] === Symbol,
        invs1 = GTInstallGroup[invs],
        invs1 = invs
     ];
     (*---tests---*)
     If[fast,
     	None,
        If[GTGroupQ[grp1],
             None,
             Print[grp1, " ist not a group"]; Abort[];
        ];
        If[invs == {},
           None,
           If[GTGroupQ[invs1],
              None,
              Print[invs1, " ist not a group"]; Abort[];
           ];
           If[GTInvSubGroupQ[grp1, invs1],
              If[GTGroupOrder[invs1] == GTGroupOrder[grp1]/2,
                 None,
                 Print["Error: Not an invariant subgroup of index 2"]; 
                 Abort[]
              ],
              Print["Error: ", invs1, " is not an invariant subgroup!"]; Abort[]
           ]
        ]
     ];
     If[invs1 == {},
    (*---Shubnikoc group of first kind---*)
        mgrp = Join[grp1, Map[(Ee)'\[SmallCircle]# &, grp1]];
        If[verb,
           Print["Shubniov group of first kind constructed"];
           Print["Order group ", Length[grp], " Order magnetic group ", Length[mgrp]],
           None
        ],
     (*---Shubnikov group of second kind" ---*)
        cmpl = Complement[grp1, 1 invs1];
        mgrp = Join[invs1, Map[(Ee)'\[SmallCircle]# &, cmpl]];
        If[verb,
           Print["Shubniov group of second kind constructed"];
           Print["Order group ", Length[grp], " Order magnetic group ", Length[mgrp]],
           None
        ]
     ];
     Return[mgrp]
]

(*
***)



(****g* /GTSymmetrizeProductChars
! NAME
!  GTSymmetrizeProductChars
! AUTHOR
!  W. Hergert
! PACKAGE
!  test.m -> RepresentationTheory.m
! MODIFICATION HISTORY
!   07/15/2016 : first version
! USAGE
!  GTSymmetrizedProductChars[ct,irep,power,symm] gives the characters of the nth power of an irreducible 
!  represantation of a given group
! INPUT
!  ct    - character table of the group under consideration
!  irep  - number of the irreducibel representation  
!  power - the power of the Irep to be investigated
!  sym   - decide if symmetri or antisymmetric power is considered, can be only "sym" od "asymm" 
! OUTPUT
!  
! OPTIONS
!  GOVerbose
!  
! ERROR MESSAGES
!  If GOFast->False it will be checked if grp and invs are groups and if invs is an invariant subgroup of index 2 of 
!  group. If the checks are negative, the commant will stop with Abort[]
! GTPack MODULES
!  GTMagnetic, GTInstallGroup, GTGroupQ, GTInvSubGroupQ
! GTPack NOTEBOOKS 
!  GTIcosaeder_new.nb in Wolfram_Devel/4_Discrete_Symmetry_Groups/NMagnetic_Groups_Lattices/Colour_Groups
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
(*
GTSymmetrizedProductCharsA[ct_, irep_, power_, sign_,OptionsPattern[]] := Module[
	{classes, chars, names, fac, pch, ncl, head, pchar, g, gg, c1,
     c2, c3, c4, c5, pos, i, tab, tn, no,verb},
  (*--- options, setup calculation ---*)  
   verb=OptionValue[GOVerbose]; 
   classes = ct[[1]]; chars = ct[[2]]; names = ct[[3]];
   ncl = Length[classes];
   head = {}; pchar = {};
  (*--- test of input ---*)
   If[sign == "symm" || sign == "asymm",
      None,
      Print["Error: ", sign, " not valid, has to be symm or asymm"]
   ];
   If[power > 5||power<2,
     Print["Error: only cases for powers from 2 to 5 are implemented"],
     None
   ];
   If[sign == "symm",
      fac = +1; pch = "[T]"^power,
      fac = -1; pch = "{T}"^power
   ]; 
  (*---  powers of the Irep ---*)
   If[power == 2,
      Do[
         g     = classes[[i, 1]]; head = Append[head, g];
         c1    = chars[[irep, i]]; 
         gg    = g\[SmallCircle]g; 
         pos   = Flatten[Position[classes, gg]][[1]];
         c2    = chars[[irep, pos]];
         pchar = Append[pchar, (c1^2 + fac* c2)/2]
      , {i, 1, ncl}],
      None
   ];
   If[power == 3,
      Do[
         g     = classes[[i, 1]]; head = Append[head, g];
         c1    = chars[[irep, i]]; 
         gg    = g\[SmallCircle]g; 
         pos   = Flatten[Position[classes, gg]][[1]];
         c2    = chars[[irep, pos]];
         gg    = g\[SmallCircle]gg; 
         pos   = Flatten[Position[classes, gg]][[1]];
         c3    = chars[[irep, pos]];
         pchar = Append[pchar, c1^3/6 + fac c2*c1/2 + c3/3]
      , {i, 1, ncl}],
      None
   ];
   If[power == 4,
      Do[
         g     = classes[[i, 1]]; head = Append[head, g];
         c1    = chars[[irep, i]]; 
         gg    = g\[SmallCircle]g; 
         pos   = Flatten[Position[classes, gg]][[1]];
         c2    = chars[[irep, pos]];
         gg    = g\[SmallCircle]gg; 
         pos   = Flatten[Position[classes, gg]][[1]];
         c3    = chars[[irep, pos]];
         gg    = g\[SmallCircle]gg; 
         pos   = Flatten[Position[classes, gg]][[1]];
         c4    = chars[[irep, pos]];    
         pchar = Append[pchar, 
                        fac*c4/4 + c3*c1/3 + c2^2/8 + fac  c2 c1^2/4 + c1^4/24]
      , {i, 1, ncl}],
      None
   ];
   If[power == 5,
      Do[
         g     = classes[[i, 1]]; head = Append[head, g];
         c1    = chars[[irep, i]]; 
         gg    = g\[SmallCircle]g; 
         pos   = Flatten[Position[classes, gg]][[1]];
         c2    = chars[[irep, pos]];
         gg    = g\[SmallCircle]gg; 
         pos   = Flatten[Position[classes, gg]][[1]];
         c3    = chars[[irep, pos]];
         gg    = g\[SmallCircle]gg; 
         pos   = Flatten[Position[classes, gg]][[1]];
         c4    = chars[[irep, pos]];
         gg    = g\[SmallCircle]gg; 
         pos   = Flatten[Position[classes, gg]][[1]];
         c5    = chars[[irep, pos]];    
         pchar = Append[pchar, 
                        c5/5 + fac c4*c1/4 + fac c3 c2 /6  + 
                        c3 c1^2/6 + c2^2 c1/8 + fac c2 c1^3/12 + c1^5/120]
      , {i, 1, ncl}],
      None
   ];
  (*--- construct table ---*)
  If[verb,
     tab = Append[Prepend[chars, head], pchar] // Transpose;
     tn  = {" ", names, pch} // Flatten; no = Length[tn];
     tab = Prepend[tab, tn] // Transpose;
     Print[Grid[tab, Frame      -> All, 
                     Dividers   -> {{2 -> GTDividerColor1}, 
                     	            {2 -> GTDividerColor1,no -> GTDividerColor1}
                     	           }, 
                     Background -> {{1 -> GTBackGroundColor1}, 
                     	            {1 -> GTBackGroundColor1, no -> GTCornerColor}, 
                     	            {1, 1} - GTCornerColor
                     	           }
                ]
           ],
      None
    ]; 
    Return[pchar]  
  ]
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

(****g* /GTExportVASP
! NAME
!  GTExportVASP
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  test.m -> CrystalStructure.m
! MODIFICATION HISTORY
!   07/14/2016 : first version
! USAGE
!  GTExportVASP exports a structure to a POSCAR file.
! INPUT
!  structure
! OUTPUT
!  POSCAR file
! OPTIONS
!  
! ERROR MESSAGES
!  
! GTPack MODULES
! 
! GTPack NOTEBOOKS 
!  GTExportVASP.nb in Matthias/Commands
! DESCRIPTION
! 
! LITERATURE
!  
! TODO
! 
! PROBLEMS
!   Documentation file missing
!--------------------------------------------------------------------------------
! SOURCE
*)


GTExportVASP[strin_] := Module[{str, atompos,atoms,chars,nn,atomposnew,atomsdiff,position,coordinates,poscartable,mm,i},
str = strin /.strin[[8]];

atompos = str[[7]];
atoms = Table[chars = Characters[atompos[[nn, 2]]]; If[Length[chars] > 1 && Not[NumberQ[ToExpression[chars[[2]]]]],chars[[1]] <> chars[[2]], chars[[1]]], {nn, 1, Length[atompos]}];
atomposnew = Table[{atompos[[nn, 1]], atoms[[nn]]}, {nn, 1, Length[atoms]}];
atomsdiff = Union[atoms];
position = Table[Position[atomposnew, atomsdiff[[nn]]], {nn, 1,Length[atomsdiff]}];
coordinates = Flatten[Table[Table[atomposnew[[position[[mm, nn, 1]], 1]], {nn, 1,Length[position[[mm]]]}], {mm, 1, Length[position]}], 1];
poscartable = Flatten[{{"Created by GTPack"}, {"1.0"}, N@str[[6]], {ToString@Row[Table[ToString[atomsdiff[[i]]] <> " ", {i, 1,Length[atomsdiff]}]]}, {ToString@Row[Table[ToString[Length[position[[nn]]]] <> " ", {nn, 1,Length[position]}]]}, {"Cartesian"}, N@coordinates}, 1];
Export["POSCAR", poscartable, "Table"]        
]

(*
***)

Attributes[GTSymmetryElementQ]={Protected, ReadProtected}
End[]
	 
EndPackage[]
