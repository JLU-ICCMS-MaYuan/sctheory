(****m* /Vibrations.m
!
! NAME
!  Vibrations.m
! AUTHOR
!  W. Hergert
! MODIFICATION HISTORY
!  * 07/20/2016  : initially created and documented  
!  * 28.12.2017  : check headers and documentation
! USAGE
!  Moduls for vibrational analysis of molecules and solids
!
! GTPack MODULES
!
! --- vibration modes of molecules ---
!  * GTVibDisplacementRep   - gives the displacement representation
!  * GTVibModeSymmetry      - Constructs the group from generators
!  * GTVibSpectroscopy      - Auftreten von IR und Raman Moden
!  
! --- vibrational modes of solids ---
!  * GTVibDynamicalMatrix   - dynamical matrix from a structure
!  * GTVibSetParameters     - sets spring constants and masses in a dynamical matrix
!  * GTVibTbToPhonon        - transforms a TB-p-Hamiltonian into a dynamical matrix
!  * GTVibTbToPhononRule    - rules to transdform the TB-parameters
!  * GTVibLatticeModes      - vibrational modes in a lattice
!
! --- internal procedures ---
!  * GTVibForceConstant     - calculation of force constants
!  * GTVibPotential         - set up of potential
!
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
! -
! PROBLEMS
! -
!--------------------------------------------------------------------------------
!
***)

BeginPackage["GroupTheory`Vibrations`",{"GroupTheory`RepresentationTheory`","GroupTheory`Basic`","GroupTheory`Symbols`","GroupTheory`Install`","GroupTheory`Auxiliary`","GroupTheory`Lattice`","GroupTheory`TightBinding`",
	"GroupTheory`CrystalStructure`","GroupTheory`PseudoPotential`","GroupTheory`Photonics`","GroupTheory`Molecules`"}]
 
(*------------ vibrational modes of molecules  -----------------------*)
 GTVibDisplacementRep    ::usage=  "GTVibDisplacementRep[\*StyleBox[\"molecule data\",\"TI\"]] gives the displacement representation of a molecule, defined by \*StyleBox[\"molecule data\",\"TI\"]." 
 GTVibModeSymmetry       ::usage=  "GTVibModeSymmetry[\*StyleBox[\"molecule data\",\"TI\"]] gives the vibrational modes of a molecule."
 GTVibSpectroscopy       ::usage=   "GTVibSpectroscopy[\*StyleBox[\"vibrational modes, character table,index\",\"TI\"]] gives seleection rules for IR and Raman spectroscopy."
  
(*------------ vibrational modes of solids  -----------------------*)
 GTVibDynamicalMatrix    ::usage=  "GTVibDynamicalMatrix[\*StyleBox[\"basis, cut off, nshell\",\"TI\"]] gives the dynamical matrix from a structure with \*StyleBox[\"basis\",\"TI\"] and one atom per unit cell. A cluster up to the \*StyleBox[\"cut off\",\"TI\"] radius is calculated. \*StyleBox[\"Nshell\",\"TI\"] shells are taken into account."
 GTVibSetParameters      ::usage=  "GTVibSetParmaeters[\*StyleBox[\"dynnamical matrix, spring constants, masses\",\"TI\"]] substitutes \*StyleBox[\"spring constants\",\"TI\"] and \*StyleBox[\"masses\",\"TI\"] in a \*StyleBox[\"dynnamical matrix\",\"TI\"]."
 GTVibTbToPhonon         ::usage=  "GTVibTbToPhonon[\*StyleBox[\"p-Hamiltonian, shells\",\"TI\"]] transforms a tight-binding \*StyleBox[\"p-Hamiltonian\",\"TI\"] in a dynamical matrix, taking into account \*StyleBox[\"shells\",\"TI\"]."
 GTVibTbToPhononRule     ::usage=  "GTVibTbToPhononRules[\*StyleBox[\"shells\",\"TI\"]] gives rules to transform a tight-binding p-Hamiltonian into a dynamical matrix."
 GTVibLatticeModes       ::usage=  "Lattice modes in a solid."

(*--------------------------- Options --------------------------------*)
 Options[GTVibDisplacementRep] = {GOFast -> GOFastValue,GOVerbose -> True}
 Options[GTVibModeSymmetry]    = {GOFast -> GOFastValue, GOVerbose -> True, GOIrepNotation -> "Bethe"}
 Options[GTVibTbToPhononRule]  = {GOTbBasis -> 0}
 Options[GTVibSetParameters]   = {GOTbBasis -> 0}
 Options[GTVibTbToPhonon]      = {GOVerbose -> True, GOTbBasis -> 0}
 Options[GTVibDynamicalMatrix] = {GOVerbose-> False}
 Options[GTVibLatticeModes]    = {GOIrepNotation -> "Bethe", GOVerbose -> True, GOTbEquivalence -> True, GOLattice -> {}}
 Options[GTVibSpectroscopy]     = {GOVerbose -> False, GOSpecMode -> "Raman"}


pA        ::usage = "Atom A"
pB        ::usage = "Atom B"
xyz1      ::usage = "Dynamical Matrix"
vx        ::usage = "Dynamical Matrix"
vy        ::usage = "Dynamical Matrix"
vz        ::usage = "Dynamical Matrix"


Begin["`Private`"] 


(****u* /GTVibSpectroscopy
! AUTHOR
!  W. Hergert
! PACKAGE
!   Molecules.m
! MODIFICATION HISTORY
!  * 31.12.2016 : first version
!  * 28.12.2017 : check header and documentation
!  * 24.06.2018 : check header and documentation
! USAGE
! GTVibSpectroscopy[vibrational modes,character table] calculates
! selection rules for IR and Raman spectroscopy
!
! INPUT
!  vibrational modes - the modes from the normal mode analysis
!  character table   - character table of the point group of the molecule
!
! OUTPUT
! Table with allowed transitions
!
! GTPack OPTIONS
!  * GOVerbose:  
!
!      - True  - additional information
!      - False - no additional information (standard)
!  * GOSpecMode:
!
!      - "Raman" - Raman active modes (Standard)
!      - "IR"    - IR active modes
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTSymmetryBasisFunctio
ns, 
!  GTDirectProductChars, 
!  GTIrep
! GTPack NOTEBOOKS
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  It works only for point groups, i.e. symmorphic space groups. A general extension to space groups would be helpful.
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)

 GTVibSpectroscopy[vibmodes_, ctab_, OptionsPattern[]] := Module[
 	{verb,mode,func,head,nf,names,chars,basis,test,pos,i,tabp,irs,
     polt,ndp,poln,tab,id,c1,n1,iv,k,c2,n2,iract,xyzl,l,nvib,xw,yw,zw,
     rule},
     rule = {xw -> Global`x, yw -> Global`y, zw -> Global`z};
  (*--- options ---*)
     verb = OptionValue[GOVerbose];
     mode = OptionValue[GOSpecMode];
     If[Intersection[{mode}, {"Raman", "InfraRed"}] == {},
        Print["Error : Type not allowed "]; Abort[],
        None
     ];
     If[mode == "Raman",
        func = {xw^2, yw^2, zw^2, xw*zw, xw*yw, yw*zw} /.rule;
        head = {"\[CapitalGamma](\!\(\*SubscriptBox[\(\[Alpha]\), \(ij\)]\))", 
        	    "component", "allowed mode"},
        func = {xw, yw, zw} /. rule;
        head = {"\[CapitalGamma](\!\(\*SubscriptBox[\(\[Mu]\), \(k\)]\))", "component", "allowed mode"}
     ];
     nf    = Length[func];
     nvib  = Length[vibmodes];
     names = ctab[[3]];
     chars = ctab[[2]];  
  (*--- find IRs of dipole moment or polarizability ---*)
     basis = GTSymmetryBasisFunctions[ctab, func, GOVerbose -> False];
     test  = Drop[basis[[2 ;; nf + 1]] // Transpose, 1] // Transpose;
     tabp  = {}; 
     Do[
     	pos = Position[test[[i]], func[[i]]];
        If[pos == {},
           None,
           irs  = Take[pos // Transpose, 1] // Flatten;
           tabp = Append[tabp, irs]
        ]
     ,{i, 1, nf}];
     polt = tabp // Flatten // Union; 
     ndp  = Length[polt];
     poln = names[[#]] & /@ polt;
     If[mode == "Raman",
        Print["Polarizability in Ireps : ", poln],
        Print["Dipole in Ireps : ", poln]
     ];
  (*--- calculate selection rules ---*)
     tab = {};
     Do[ 
        id = polt[[i]];
        c1 = chars[[id]]; 
        n1 = names[[id]];
        Do[
            iv = vibmodes[[k]]; 
            c2 = chars[[iv]]; 
            n2 = names[[iv]];
            iract = GTIrep[GTDirectProductChars[c1, c2], ctab, GOVerbose -> verb];
            If[iract[[1]] > 0,
               xyzl = {};
               Do[
                  pos = Position[tabp[[l]], id] // Flatten;
                  If[pos == {},
                     None,
                     xyzl = Append[xyzl, func[[l]]/. {x->Global`x,y->Global`y,z->Global`z}]
                  ]
               , {l, 1, nf}];
               tab = Append[tab, {n1, xyzl, n2}],
               None
            ]
        , {k, 1, nvib}]
     , {i, 1, ndp}];
     If[tab == {},
        Print["No active modes found!"],
        tab = Prepend[tab, head];
        Print[Grid[tab, Frame -> All, Background -> {None, 1 -> GTBackGroundColor1}]]
     ]
 ]
  
(*
***)

(****u* /GTVibDisplacementRep
! NAME
!  GTVibDisplacementRep
! AUTHOR
!  W. Hergert
! PACKAGE
!   Vibrations.m
! MODIFICATION HISTORY
!  * 08/01/2016  : first version
!  * 28.12.2017  : check headers and documentation
!  * 24.06.2018  : check headers and documentation
! USAGE
!  GTVibDisplacementRep[molecule data] calculates the displacement representation for a molecule.
! INPUT
!  * molecular data -  data from the database
! OUTPUT
!  * displacement representation
! GTPack OPTIONS
!  * GOFast - $GOFastValue used to avoid additionl tests
!  * GOVerbose:  
!
!      - True  - additional information
!      - False - no additional information (standard)
! STANDARD OPTIONS
!  -
! GTPack MODULES
!   GTMolToCluster,
!   GTSymmetryElementQ,
!   GTGetMatrix
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
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTVibDisplacementRep[moldat_, OptionsPattern[]] := Module[
    {fast,verb, grp, gord,clu,log,vert,lmol,grp1,cc,k,nl,no,eps,j,part,i,rep,drep},
    eps=10^(-7);
  (*--- options ---*)
    fast = OptionValue[GOFast];
    verb = OptionValue[GOVerbose];
  (*--- install group ---*)
    grp  = moldat[[4]];
    grp  = GTInstallGroup[grp, GOVerbose -> False]; 
    gord = Length[grp];
  (*-- cluster from molecule data ---*)  
    clu  = GTMolToCluster[moldat, GOMolPlot -> False, GOVerbose -> False];
  (*--- check of data with respect to symmetry ---*)
    If[fast,
       Do[
          log = GTSymmetryElementQ[grp[[i]], clu];
          If[log,
             None,
             Print["Error: Check symmetry!"]; Abort[]
          ]
       , {i, 1, gord}];
       If[verb,
             Print["correct symmetry properties."],
             None
       ]
    ];
  (*--- cluster data ---*)  
    vert = Transpose[clu][[1]]; lmol = Length[vert];
  (*--- construct displacement representation ---*)
    grp1 = GTGetMatrix[grp];
    cc   = Flatten[Table[ Map[grp1[[k]].# &, vert] , {k, 1, gord}], 1];
    nl   = {}; 
    Do[
       Do[
          no = Norm[vert[[k]] - cc[[j]]];
          If[no < eps,
             nl = Append[nl, k],
             None
          ]
       , {k, 1, lmol}]
    , {j, 1, Length[cc]}];
    part = Partition[nl, lmol];
    rep = Table[0, {k, 1, gord}, {i, 1, lmol}, {j, 1, lmol}];
    Do[
       rep[[k, part[[k, i]], i]] = 1, {i, 1, lmol}
    , {k, 1, gord}];
    drep = Table[0, {gord}];
    Do[
       drep[[i]] = KroneckerProduct[rep[[i]], grp1[[i]]]
    , {i, 1, gord}];
    Return[drep]
]

(*
***)



(****u* /GTVibModeSymmetry
! NAME
!  GTVibDisplacementRep
! AUTHOR
!  W. Hergert
! PACKAGE
!   Vibrations.m
! MODIFICATION HISTORY
!  * 08/01/2016  : first version
!  * 28.12.2017  : check headers and documentation
!  * 28.12.2017  : check headers and documentation
! USAGE
!  GTVibModeSymmetry[moldat] gives the vibrational modes of a molecule.
! INPUT
!  molecular data -  data from the database
! OUTPUT
!  vibrational modes
! GTPack OPTIONS
!  * GOIrepNotation - select the notation for the ireducible representations
!  * GOFast - $GOFastValue used to avoid additionl tests
!  * GOVerbose:  
!
!      - True  - additional information
!      - False - no additional information (standard)
! STANDARD OPTIONS
!  -
! GTPack MODULES
!   GTVibDisplacementRep,
!   GTCharacters,
!   GTInstallGroup,
!   GTCharacterTable,
!   GTIrep,
!   GTGetMatrix
! GTPack NOTEBOOKS
!  -
! DESCRIPTION
!  - 
! LITERATURE
!  -
! TODO
!  1.0.0
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)


GTVibModeSymmetry[moldat_, OptionsPattern[]] := Module[
	{verb,fast,inot,drep,dim,chdrep,symg,grp,ct,classes,names,ncl,red,cht,chr,symel,i,dt,ca,chrot,
     chtrans,np,tb,tb2,nrop,l,tb3},
  (*--- options ---*)
    verb = OptionValue[GOVerbose];
    fast = OptionValue[GOFast];
    inot = OptionValue[GOIrepNotation];
  (*--- construction of displacement representation ---*)  
    drep = GTVibDisplacementRep[moldat, GOVerbose -> verb, GOFast -> fast];
    If[verb,
       Print["displacment representation:"];
       dim = Length[drep[[1]]];
       Print[Length[drep], " matrices of dimension ", dim, " x ",dim],
       None
    ];
    chdrep = GTCharacters[drep, GOClasses -> True];
    symg   = moldat[[4]];
    If[verb,
       Print["characters disp. rep. : ", chdrep];
       Print["symmetry group        : ", symg],
       None
    ];
  (*--- reduction of displacement rep ---*)  
    grp = GTInstallGroup[symg, GOVerbose -> False];
    ct  = GTCharacterTable[grp, GOIrepNotation -> inot, GOVerbose -> verb];
    classes = ct[[1]]; names = ct[[3]]; ncl = Length[classes];
    If[verb,
       Print["Reduction of direct product representation:"],
       None
    ];
    red = GTIrep[chdrep, ct, GOVerbose -> verb];
  (*--- characters of translation and rotation rep ---*)  
    cht = {}; chr = {}; 
    Do [
    	symel = classes[[i, 1]]; 
    	dt    = Det[GTGetMatrix[symel]];
        ca    = GTGetQuaternion[symel][[1]]; ca = 4 ca^2 - 1;
        chr   = Append[chr, ca];
        If[dt == 1,
           cht = Append[cht, ca],
           cht = Append[cht, -ca]
        ]
    , {i, 1, ncl}];
    If[verb,
       Print["Ireps of rotations:"],
       None
    ];
    chrot = GTIrep[chr, ct, GOVerbose -> verb];
    If[verb,
       Print["Ireps of translations:"],
       None
    ];
    chtrans = GTIrep[cht, ct, GOVerbose -> verb];
    np = red - chrot - chtrans; 
  (*--- prepare print of result (Matthias Irep) ---*)  
    tb = Flatten[
                 Table[
                	   Which[np[[i]] == 0, "a", 
                             np[[i]] == 1, {ct[[3, i]], "\[CirclePlus]"}, 
                             np[[i]]  > 1, {np[[i]] ct[[3, i]], "\[CirclePlus]"}]
                 , {i, 1, Length[np]}]
               ];
    tb2  = Delete[tb, Position[tb, "a"]];
    nrop = Length[Position[tb2, "\[CirclePlus]"]];
    l    = Length[tb2] - Length[Position[tb2, "\[CirclePlus]"]];
    tb3  = If[l == nrop, 
                   Delete[tb2, Position[tb2, "\[CirclePlus]"][[nrop]]], 
                   None
           ];
    Print[Row[tb3]]; 
    Return[np]
  ]


(*
***)





(****u* /GTVibTbToPhononRule
! NAME
!  GTVibTbToPhononRule
! AUTHOR
!  W. Hergert
! PACKAGE
!   Vibrations.m
! MODIFICATION HISTORY
!  * 08/03/2016  : first version
!  * 28.12.2017  : check headers and documentation
!  * 24.06.2017  : check headers and documentation
! USAGE
!  GTVibTbToPhononRules[shells] gives rules to transform a tight-pinding p-Hamiltonian into a dynamical matrix
! INPUT
!  shells -  information about the shells used in the construcution of the tight-binding Hamiltonian
! OUTPUT
!  list of rules
! GTPack OPTIONS
!  GOTbBasis
! STANDARD OPTIONS
!  -
! GTPack MODULES
!   GTTbSymbol2C
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
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTVibTbToPhononRule[shells_, OptionsPattern[]] := Module[
	{bas,nb,baslist,bas1,rule,sym,i,j,mass},
  (*--- option ---*)
     bas   = OptionValue[GOTbBasis];
  (*--- check input  ---*)
     nb    = Length[shells];
     If[nb == 1 || nb == 3,
        None,
        Print["Error: Length of shell has to be 1 or 3!"]; Abort[]
     ];
     If[bas == 0,
        baslist = {},
        bas1 = StringSplit[bas, ","];
        If[Length[bas1] > 2,
           Print["Error: To many atoms!"]; Abort[],
           None
        ];
        bas1    = Sort[bas1];
        baslist = {bas1[[1]] <> "," <> bas1[[1]], 
                   bas1[[2]] <> "," <> bas1[[2]], 
                   bas1[[1]] <> "," <> bas1[[2]]
                  }
    ];
  (*--- construction of the substitution rules , one atom ---*)
    rule = {};
    If[baslist == {},
       Do[
  (*--- (ppp) is set to 0 ---*)
          sym  = GTTbSymbol2C[1, 1, 1, i];
          rule = Append[rule, sym -> 0];
  (*--- (pps) is substituted by k/m ---*)        
          sym  = GTTbSymbol2C[1, 1, 0, i];
          rule = Append[rule, sym -> - Subscript["k", i]/"m"];
      , {i, 1, shells[[1]]}],
  (*--- construction of the substitution rules , two atoms ---*)   
      Do[
         Do[
  (*--- (ppp) is set to 0 ---*)                              
            sym  = GTTbSymbol2C[1, 1, 1, i, GOTbBasis ->  baslist[[j]]];
            rule = Append[rule, sym -> 0];
  (*--- (pps) is usbstituted by k/m ---*)   
            sym  = GTTbSymbol2C[1, 1, 0, i, GOTbBasis ->  baslist[[j]]];
            If[j < 3,
               mass = Superscript["m", bas1[[j]]],
               mass = Sqrt[Superscript["m", bas1[[1]]] Superscript["m", bas1[[2]]] ]
            ];
            rule = Append[rule, sym -> - Subsuperscript["k", i, baslist[[j]]]/mass]
         , {i, 1, shells[[j]]}]
     , {j, 1, nb}]
    ]; 
    Return[rule]
  ]


(*
***)


(****u* /GTVibTbToPhonon
! NAME
!  GTVibTbToPhonon
! AUTHOR
!  W. Hergert
! PACKAGE
!   Vibrations.m
! MODIFICATION HISTORY
!  * 08/03/2016  : first version
!  * 28.12.2017  : check headers and documentation
!  * 24.06.2017  : check headers and documentation
! USAGE
!  GTVibTbToPhonon[ham,shells] transforms a TB p-Hamiltonian in a dynamical matrix
! INPUT
!  * hamiltonian - a TB Hamiltonian
!  * shells      -  information about the shells used in the construcution of 
!                 the tight-binding Hamiltonian
! OUTPUT
!  list of rules
! GTPack OPTIONS
!  * GOTbBasis - 
!  * GOVerbose:  
!
!      - True  - additional information
!      - False - no additional information (standard)
!  GOVerbose
!  GOTbBasis
! STANDARD OPTIONS
!  -
! GTPack MODULES
!   GTVibTbToPhononRule,
!   GTTbSymbol2C
!  
! GTPack NOTEBOOKS
   -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  - 
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)

 GTVibTbToPhonon[ham_, shells_, OptionsPattern[]] := Module[
 	 {verb,bas,bas1,rule,ham1,htest,glg,sym,glg1,sol,ham2,hamf},
  (*--- options ---*)
      verb = OptionValue[GOVerbose];
      bas  = OptionValue[GOTbBasis];
  (*--- interpretation of basis ---*)    
      If[bas == 0,
         None,
         bas1 = StringSplit[bas, ","];
         If[Length[bas1] > 2,
            Print["Error: To many atoms!"]; Abort[],
            None
         ]
      ];
  (*--- generation of substitution rules ---*)
      rule  = GTVibTbToPhononRule[shells, GOTbBasis -> bas];
      ham1  = ham /. rule;
      htest = ham1 /. {\[Xi] -> 0, \[Eta] -> 0, \[Zeta] -> 0};
      If[verb,
      	 Print["Matrix at the \[CapitalGamma] point:"];
      	 Print[htest//MatrixForm],
      	 None
      ];	 
      If[bas == 0,
  (*--- one atom in basis ---*)
         glg = Map[Plus @@ # &, htest] // Union;
         sym = GTTbSymbol2C[1, 1, -1, 0];
         glg1=glg[[1]] == 0 /. {sym -> pA};
         sol = Solve[glg1, pA];sol=Flatten[sol];
         If[verb,
            Print[{{"equation onsite-element"," : ",glg1},
                   {"solution"," : ",sol}
                  }//TableForm],
            None
         ];         
         ham2 = ham1 /. sym -> pA; 
         hamf = ham2 /. sol,
  (*--- two atoms in basis ---*)         
         glg  = Map[Plus @@ # &, htest] // Union;
         glg  = Map[# == 0 &, glg]; 
         rule = {GTTbSymbol2C[1, 1, -1, 0, GOTbBasis -> bas1[[1]]] -> pA, 
                 GTTbSymbol2C[1, 1, -1, 0, GOTbBasis -> bas1[[2]]] -> pB
                }; 
         glg1 = glg /.rule;
         sol  = Solve[glg1, {pA, pB}];
         sol  = Flatten[sol];
         If[verb,
            Print[{{"equations onsite-element"," : ",glg1},
                   {"solution"," : ",sol}
                  }//TableForm],
            None
         ];  
         ham2 = ham1 /. rule; 
         hamf = ham2 /. sol
      ];
      Return[hamf]
  ]


(*
***)




(****u* /GTVibSetParameters
! NAME
!  GTVibSetPArameters
! AUTHOR
!  W. Hergert
! PACKAGE
!   Vibrations.m
! MODIFICATION HISTORY
!  * 08/04/2016 : first version
!  * 28.12.2017  : check headers and documentation
!  * 24.06.2017  : check headers and documentation
!
! USAGE
!   GTVibSetParmaeters[dmat,spring,masses] substitutes parameter values in dynamical matrix. 
! INPUT
!  * dmat   - dynamical matrix
!  * spring - spring constants
!  * masses - masses of the atoms
! OUTPUT
!  parametrized dynamical matrix
! GTPack OPTIONS
!  * GOTbBasis - necessary if basis appears
! STANDARD OPTIONS
!  -
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
!  1.0.0
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)


GTVibSetParameters[dmat_, spring_, masses_, OptionsPattern[]] := Module[
	{bas,bas1,bas2,baslist,ma,dmat1,ls,rule,i,j,k},
  (*--- options ---*)
     bas = OptionValue[GOTbBasis];
  (*--- interpretation of basis ---*)
     If[bas == 0,
        None,
        bas1 = StringSplit[bas, ","];
        bas2 = bas1 // Sort;
        If[Length[bas1] > 2,
           Print["Error: To many atoms!"]; Abort[],
           baslist = {bas1[[1]] <> "," <> bas1[[1]], 
                      bas1[[2]] <> "," <> bas1[[2]], 
                      bas2[[1]] <> "," <> bas2[[2]]
                     }
        ]
     ];
  (*--- one atom in basis ---*)
     If[bas == 0,
        ma    = masses[[1]];
        dmat1 = dmat /. {"m" -> ma};
        ls    = Length[spring];
        rule  = {};
        Do[
           rule = Append[rule, Subscript["k", i] -> spring[[i]]]
        , {i, 1, ls}],
  (*--- two atoms in basis --*)
        dmat1 = dmat /. {Superscript["m", bas1[[1]]] -> masses[[1]], 
                         Superscript["m", bas1[[2]]] -> masses[[2]]
                        };
        rule  = {};
        Do[ 
           ls = Length[spring[[k]]];
           Do[
           	  rule = Append[rule, Subsuperscript["k", j, baslist[[k]]] -> spring[[k, j]]]
           , {j, 1, ls}]
        , {k, 1, 3}]
     ];
     dmat1 = dmat1 /. rule;
     Return[dmat1]
 ]
 
(*
***)


(****u* /GTVibDynamicalMatrix
! NAME
!  GTVibDynamicalMatrix
! AUTHOR
!  W. Hergert
! PACKAGE
!   Vibrations.m
! MODIFICATION HISTORY
!  * 08/04/2016  : first version
!  * 28.12.2017  : check headers and documentation
!  * 24.06.2017  : check headers and documentation
! USAGE
!   GTVibDynamicalMatrix[basis, cutoff, nshell] constructs the dynamical matrix from a given strucutre
! INPUT
!  * basis  - basis of the lattice
!  * cutoff -  cutoff in cluster construction
!  * nshell - number of shells taken into account
!  masses - masses of the atoms
! OUTPUT
!   dynamical matrix
! GTPack OPTIONS
!  * GOVerbose:  
!
!      - True  - additional information
!      - False - no additional information (standard)
! STANDARD OPTIONS
!  -
! GTPack MODULES
!   GTLatCluster,
!   GTLatShells,
!   GTVibPotential,
!   GTVibForceConstant
! GTPack NOTEBOOKS
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  Extension to more atoms per cell?
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTVibDynamicalMatrix[basis_, rcut_, ns_,OptionsPattern[] ]:= Module[
	{cl,sh0,ls,sh1,pot,const,na,rv,k,j,l,dynmat,i,verb},
  (*--- option ---*) 
    verb  = OptionValue[GOVerbose];	 
    cl    = GTLatCluster[basis, rcut,GOVerbose->verb]; 
    sh0   = GTLatShells[cl]; 
    ls    = Length[sh0];
    Print["Number of shells : ", ls];
    sh1   = Take[sh0, {2, ns + 2}];
    pot   = GTVibPotential[sh1]; 
    const = Table[0, {3}];
    Do[
       Do[
       	  na = Length[sh1[[l]]];
          Do[
          	 rv         = sh1[[l, k]];
             const[[j]] = const[[j]] + GTVibForceConstant[j, j, rv, pot, sh1]
          , {k, 1, na}]
       , {l, 1, ns}]
    , {j, 1, 3}];
    dynmat = Table[0, {3}, {3}]; 
    Do[
       dynmat[[i, i]] = -const[[i]]
    , {i, 1, 3}];
    If[verb,
       Print["diagonal terms"];
       Print[dynmat//MatrixForm],
       None
    ];
    Do[
       Do[
          Do[
             na = Length[sh1[[l]]];
             Do[
                rv             = sh1[[l, k]];
                dynmat[[i, j]] = dynmat[[i, j]] + 
                                 GTVibForceConstant[i, j, rv, pot, sh1]*
                                 Exp[2 \[Pi] I rv.{\[Xi], \[Eta], \[Zeta]}];
             , {k, 1, na}]
          , {l, 1, ns}];
       , {i, 1, 3}]
    , {j, 1, 3}];
    Return[dynmat // FullSimplify]
]

GTVibForceConstant[i_, j_, rv_, pot_, shells_] := Module[
	{vj, pos, sh, num, ui,fc}, 
	xyz1 = {"x", "y", "z"}; 
    vj  = Symbol["v" <> xyz1[[j]]]; 
    pos = Flatten[Position[shells, rv]];
    If[pos == {},
       Print["geht nicht1"]; Abort[],
       {sh, num} = pos;
       ui        = Symbol["u" <> ToString[sh] <> ToString[num] <> xyz1[[i]]]; 
       fc        = D[pot, ui, vj]
    ]; 
    Return[fc]
]

GTVibPotential[shells_] := Module[
	{s, ns, i, na, sat, nv, u, ia, is,spring} , 
	 s      = 0; 
     ns     = Length[shells]; 
     spring = Table[Subscript["k" ,i]/"m", {i, 1, ns}];
     Do[
        na  = Length[shells[[is]]];
        sat = shells[[is]]; nv = Norm[sat[[1]]];
        u   = Table[{Symbol["u" <> ToString[is] <> ToString[i] <> "x"], 
                     Symbol["u" <> ToString[is] <> ToString[i] <> "y"], 
                     Symbol["u" <> ToString[is] <> ToString[i] <> "z"]
                    }
              ,{i, 1,na}];
        Do[
     	   s = s + spring[[is]]/2*(sat[[ia]].(u[[ia]] - {vx, vy, vz})/nv)^2
        , {ia, 1, na}]
     , {is, 1, ns}]; 
     Return[s]
  ]

(*
***)

(****u* /GTVibLatticeModes
! NAME
!  GTVibLatticeModes
! AUTHOR
!  W. Hergert
! PACKAGE
!   Vibrations.m
! MODIFICATION HISTORY
!  * 02.01.2017 : first version
!  * 28.12.2017 : check headers and documentation
!  * 29.09.2021 : input of group as a list allowed 
! USAGE
!   GTVibLatticeModes[structure, cutoff, point group] gives vibration modes of a lattice
! INPUT
!  * structure     - structure file from database
!  * cutoff        - cutoff radius from structure construction
!  * point group   - point group of the system
! OUTPUT
!  * Vibrational modes of the lattices at zone center
! GTPack OPTIONS
!  * GOTbEquivalence
!  * GOLattice
!  * GOIrepNotation
!  * GOVerbose:  
!
!      - True  - additional information
!      - False - no additional information (standard)
! Standard OPTIONS
!  -
! GTPack MODULES
!   GTCluster,
!   GTInstallGroup,
!   GTCharacterTable,
!   GTGetMatrix,
!   GTIrep,
!   GTVectorRep
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
!  - 
!
!--------------------------------------------------------------------------------
! SOURCE
*)


GTVibLatticeModes[struc_, rcut_, group_, OptionsPattern[]] := Module[
   {verb, equi, latt, inot,bas,tbas,atoms,mult,nat,sub,clu,clt,g,pos,cc,la,coord,i,j,t,
   grp,ctab,lc,selm,mats,chars1,test,k,eqrep,vrep,efac,epos,vfac,vpos,np,chars,neq,nvec,
   ie,jv,ctest,rep,tb,tb2,tb3,nrop},
  (*--- options ---*)
    verb = OptionValue[GOVerbose];
    equi = OptionValue[GOTbEquivalence];
    latt = OptionValue[GOLattice];
    inot = OptionValue[GOIrepNotation];
  (*--- discussion of the  basis ---*)
    bas = struc[[7]] /. latt;
    If[equi,
       tbas  = Transpose[bas]; 
       atoms = tbas[[2]],
       tbas  = Transpose[bas];
       mult  = Tally[tbas[[2]]]; 
       nat   = Length[mult];
       atoms = {}; 
       sub   = {};
       Do[
          t     = Table[mult[[i, 1]] <> ToString[j], {j, 1, mult[[i, 2]]}];
          atoms = Flatten[Append[atoms, t]]
       , {i, 1, nat}];
       tbas[[2]] = atoms; 
       bas = Transpose[tbas]
    ];
    If[verb,
       Print["Atoms in Basis:"]; 
       g=Grid[{tbas[[2]], tbas[[1]]}, Frame -> All, Dividers -> {Black, {2 -> GTDividerColor1}}, 
             Background -> {None, {1 -> GTBackGroundColor1}}
            ]; 
       Print[g],
       None
    ];
  (*--- construct cluster,  sort atoms ---*)
    If[verb,
       Print["Cluster construction:"],
       None
    ];
    clu   = GTCluster[struc, rcut, GOLattice -> latt, GOTbEquivalence -> equi];
    clt   = Transpose[clu];
    la    = Length[atoms]; 
    coord = {};
    Do[
   	   pos = Position[clt[[2]], atoms[[i]]] // Flatten;
       cc = Map[clt[[1, #]] &, pos]; 
       coord = Append[coord, cc]
    , {i, 1, la}];
  (* *)
  (*coord=Flatten[coord,1];*)
  (* *)
  (*--- install group and character table ---*)
    If[verb,
       Print["Installation of group and character table :"],
       None
    ];
    grp  = GTInstallGroup[group, GOVerbose -> verb];
    ctab = GTCharacterTable[grp, GOVerbose -> verb, GOIrepNotation -> inot];
    lc   = Length[ctab[[1]]];
  (*--- representant of class as matrix ---*) 
    selm = Map[#[[1]] &, ctab[[1]]]; 
    mats = Map[GTGetMatrix[#] &, selm];
  (*--- find characters ---*)
    chars1 = Table[0, {lc}];
    If[equi,
       coord = Flatten[coord, 1];
       Do[
          Do[
          	 test = mats[[i]].tbas[[1, k]];
          	 If[Intersection[{test}, coord] == {},
                None, 
                chars1[[i]] = chars1[[i]] + 1
             ],
             {k, 1, la}]
       , {i, 1, lc}],
       Do[
          Do[
          	 test = mats[[i]].tbas[[1, k]];
             If[Intersection[{test}, coord[[k]]] == {},
                None, 
                chars1[[i]] = chars1[[i]] + 1
             ]
          , {k, 1, la}]
       , {i, 1, lc}]
    ];
    Print["Equivalence representation :"];
    eqrep = GTIrep[chars1, ctab, GOVerbose -> True];
  (*--- vector representation ---*)
    Print["Vector representation :"];
    vrep = GTVectorRep[ctab, GOVerbose -> True];
    Print["Lattice Modes :"];
    efac = Select[eqrep, # != 0 &]; 
    neq = Length[efac];
    epos = Position[eqrep, #] & /@ efac // Flatten;
    vfac = Select[vrep, # != 0 &];
    nvec = Length[vfac];
    vpos = Position[vrep, #] & /@ vfac // Flatten;
    np = Table[0, {lc}];
    chars = ctab[[2]];
    Do[
       Do[
          ie    = epos[[i]]; 
          jv    = vpos[[j]];
          ctest = chars[[ie]]*chars[[jv]];
          rep   = efac[[i]]*vfac[[j]]*GTIrep[ctest, ctab, GOVerbose -> False];
          np = np + rep
    , {i, 1, neq}]
   , {j, 1, nvec}];
  (*---prepare print of result (Matthias Irep)---*)
   tb = Flatten[
           Table[Which[np[[i]] == 0, "a", 
                 np[[i]] == 1, {ctab[[3, i]], "\[CirclePlus]"}, 
                 np[[i]] > 1, {np[[i]] ctab[[3, i]], "\[CirclePlus]"}]
           , {i, 1, Length[np]}]
        ];
   tb2  = Delete[tb, Position[tb, "a"]];
   nrop = Length[Position[tb2, "\[CirclePlus]"]];
  l = Length[tb2] - Length[Position[tb2, "\[CirclePlus]"]];
  tb3 = If[l == nrop, 
    Delete[tb2, Position[tb2, "\[CirclePlus]"][[nrop]]], None];
  Print[Row[tb3]];
  Return[np]

(*
***)



   ]

End[]
	 
EndPackage[]
