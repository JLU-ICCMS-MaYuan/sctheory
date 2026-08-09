(*
!
! List of Changes:
!  29.10.2021   Use as a package in GTPack
!  09.11.2021   name of global variable wyckoff
!
*)

(****m* /Wyckoff.m  
!
! NAME
!  Wyckoff.m
! AUTHOR
!  W. Hergert
!
! MODIFICATION HISTORY
! * 29.10.2021 : intitial documentation
!
! USAGE
!   Install, read and use Wyckoff positions for generation of structural data.
!
! GTPack MODULES
!
! --- Load, Save and Modify Wyckoff position data ---
!
!  * GTInstallWyckoff   - installation of a new space group
!  * GTAllWyckoff       - list of all installed Wyckoff positions
!  * GTGetWyckoff       - get information for a space group
!  * GTCoordFromWyckoff - generate corresponding part of GTPack structure file from Wyckoff positions
!  * GTSaveWyckoff      - save Wyckoff positions to database
!  * GTClearWyckoff     - clear Protected Variable gtwyckoff
!  * GTLoadWyckoff      - loads database with Wyckoff Positions
!
!      
! DESCRIPTION
!  It is not the aim of this part of GTPack to substitute all the sophisticated programs available in the community in this
!  field. There are also a lot of online resources available. This part was developed to be independent up to a certain level from 
!  external sources. Some commannds provide basic information about point groups and space groups. 
!
! LITERATURE
!
!  Interesting sources online are:
!
!	* http://www.cryst.ehu.es		Bilbao server
!	* http://it.iucr.org/A/		International Tables for Crystallography A: Space-group symmetry
!	* http://www.crystallography.net	Crystallography Open Database
!	* http://icsd.fiz-karlsruhe.de/	ICSD Web (Karlsruhe)
***)


BeginPackage["GroupTheory`Wyckoff`",{"GroupTheory`Symbols`","GroupTheory`Auxiliary`","GroupTheory`Lattice`","GroupTheory`Basic`","GroupTheory`RepresentationTheory`","GroupTheory`CrystalStructure`"}]

(*--------- Load, Save and Modify Wyckoff Data ------------*)

 GTClearWyckoff  	      ::usage = "GTClearWyckoff[] removes all currently installed Wckoff positions of the corresponding space groups."
 GTAllWyckoff    	      ::usage = "GTAllWyckoff[] lists all currently installed Wyckoff positions."
 GTLoadWyckoff   	      ::usage = "GTLoadWyckoff[\*StyleBox[\"database\",\"TI\"]] loads a list of installed Wyckoff positions from a \*StyleBox[\"database\",\"TI\"]."
 GTInstallWyckoff	      ::usage = "GTInstallWyckoff[\*StyleBox[\"Wyckoff Positions\",\"TI\"]] adds a new set of \*StyleBox[\"Wyckoff positions\",\"TI\"] to the global variable \*ButtonBox[\"wyckoff\", BaseStyle->\"Link\", ButtonData->\"paclet:GroupTheory/ref/wyckoff\"]."
 GTGetWyckoff   	      ::usage = "GTGetWyckoff[\*StyleBox[\"number, setting\",\"TI\"]] gives an set of Wyckoff positions to \*StyleBox[\"number\",\"TI\"] in the database for special \*StyleBox[\"settimg\",\"TI\"]."
 GTSaveWyckoff		      ::usage = "GTSaveWyckoff[\*StyleBox[\"database\",\"TI\"]] saves the content of the global variable \*ButtonBox[\"gtwyckoff\", BaseStyle->\"Link\", ButtonData->\"paclet:GroupTheory/ref/wyckoff\"] to a \*StyleBox[\"database\",\"TI\"]."
 GTCoordFromWyckoff       ::usage = "GTCoordFromWyckoff[\*StyleBox[\"num,setting,atoms\",\"TI\"]] generates a part of the GTPack structure file from Wyckoff positions from space group \*StyleBox[\"num\",\"TI\"] and \*StyleBox[\"setting\",\"TI\"] and \*StyleBox[\"atoms\",\"TI\"], a list of definitions  for atoms." 
 

(*--------------------------- Options ----------------------------*)

Options[GTGetWyckoff]            = {GOVerbose -> True}
Options[GTLoadWyckoff]           = {GOVerbose -> False}
Options[GTCoordFromWyckoff]      = {GOVerbose -> False}


Begin["`Private`"] (* Begin Private Context *) 

(*--------------------------- Modules ----------------------------*)

(****e* /GTClearWyckoff
! NAME
!  GTClearWyckoff
! AUTHOR
!  W.Hergert
! PACKAGE
!  -> CrystalStructure.m
! MODIFICATION HISTORY
!  11.06.2021 : first complete version, ready for tests
!  01.11.2021 : set up of documentation page
! USAGE
!  GTClearWyckoff[] clears the global variable wyckoff.
! INPUT
!  -
! OUTPUT
!  internal variable wyckoff={}
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  - 
! GTPack NOTEBOOKS 
!  Test_Wyckoff.nb - Notebook demonstrates the use of the different functions.
! DESCRIPTION
!   This is only a small service command. See also GTClearStructures.
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  -
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)
 
GTClearWyckoff[] := Module[{}, 
	Unprotect[wyckoff]; 
	Clear[wyckoff]; 
	wyckoff = {}; 
    Protect[wyckoff]; 
    Return[wyckoff];
] 
 
 
 (*
 ***)

(****e* /GTAllWyckoff
! NAME
!  GTAllWyckoff
! AUTHOR
!  Wolfram Hergert
! PACKAGE
!  -> CrystalStructure.m 
! MODIFICATION HISTORY
!  11.06.2021 : first complete version, ready for tests
!  01.11.2021 : documentation page setup
!  
! USAGE
! GTAllWyckoff[] lists all space groups and  Wyckoff positions actually implemented in wyckoff.
! INPUT
!  o no input necessary, internal variable wyckoff is used.
! OUTPUT
!  Table of implemented space goups, sorted with respect to space group number 
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack NOTEBOOKS 
!  Test_Wyckoff.nb
! DESCRIPTION
!  This command gives an overview on all Wyckoff positions implemented in a database.
! LITERATURE
!  International tables for Crystallography
! TODO
!  -
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTAllWyckoff::data = "No data in wyckoff."

GTAllWyckoff[] := Module[{ns, w1, tab0, tab, w2, tmp, tt, t1, t2},
    ns  = Length[wyckoff];
    If[ns==0,
       Message[GTAllWyckoff::data]; 
       Return[],
       None
    ];    	
  (*--- sort with respect to space group number ---*)
    w1   = Sort[wyckoff, #1[[1]] < #2[[1]] &] // Transpose;
    tab0 = w1[[1 ;; 2]] // Transpose;
    tab  = {};
    w2   = w1[[4]];
    Do[
       tmp = Take[#, 1 ;; 2] & /@ w2[[k]];
       tt = {}; 
       Do[
       	 t1  = ToString[tmp[[i, 1]]]; t2 = tmp[[i, 2]]; 
         tt  = Append[tt, StringJoin[t1, t2]]
       , {i, 1, Length[tmp]}];
       tt  = {tab0[[k, 1]], tab0[[k, 2]], tt};
       tab = Append[tab, tt]
    , {k, 1, ns}];
    tab = Prepend[tab, {"Space group", "Setting", "Wyckoff positions"}];
    Print[Grid[tab, Frame -> All, Background -> {None, 1 -> GTBackGroundColor1}, 
          Dividers -> {{3 -> GTDividerColor1}, {2 -> GTDividerColor2}}]
    ]
 ]
 
 (*
 ***)
 (****e* /GTLoadWyckoff
! NAME
!  GTLoadWyckoff
! AUTHOR
!  W. Hergert
! PACKAGE
!  -> CrystalStructure.m 
! MODIFICATION HISTORY
!   11.06.2021 : first complete version, ready for tests
!   01.11.2021 : documentation page setup
! USAGE
!  GTLoadWyckoff[database] loads a list of space groups with installed Wyckoff positions from a database.
! INPUT
!  database - must have the extension ".wyck"
! OUTPUT
!  the internal variable gtwyckoff will be set
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTReadFromFile
! GTPack NOTEBOOKS 
!  Test_Wyckoff.nb - Notebook demonstrates the use of the different functions.
! DESCRIPTION
!   The Wyckoff positions to the space groups are stored in one or more databases.
!   This command is used to load the Wckoff positions into the internal varable
!   wyckoff.
!   See also GTLoadStructures
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  - 
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*) 
GTLoadWyckoff::data = "Database `1` not in working directory `2`." 

GTLoadWyckoff[file0_,OptionsPattern[]] := Module[{file,dir,verb,ls},
	verb=OptionValue[GOVerbose];
	file=StringJoin[file0,".wyck"];
	(*--- check if database exists ---*)
	If[Flatten[Position[FileNames[], file]] == {},
	   dir=Directory[];
	   Message[GTLoadWyckoff::data,file,dir];
	   Return[],
       None
    ];
	(*--- load the structures ---*)
	Unprotect[wyckoff];
	wyckoff=GTReadFromFile[file];
	Protect[wyckoff];
	ls=Length[wyckoff];
	Print["Wyckoff positions of ",ls," space groups implemented in wyckoff."];
	If[verb,
	   Return[wyckoff],
	   Return[];
	]   
]

(****e* /GTInstallWyckoff
! NAME
!  GTInstallWyckoff
! AUTHOR
!  Wolfram Hergert
! PACKAGE
!  -> CrystalStructure.m 
! MODIFICATION HISTORY
!  11.06.2021 : first complete version, ready for tests
!  02.11.2021 : setup of documentation page
! USAGE
! GTInstallWyckoff[newwyc_] installs a set of Wyckoff positions for a space group into wyckoff.
! INPUT
!  o newwyc   - Wyckoff positions of a space group
! OUTPUT
!  internal list gtwyckoff
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack NOTEBOOKS 
!  Test_Wyckoff.nb
! DESCRIPTION
!  The information from The "International Tables for Crystallography"  are put in an appropriate
!  form. 
! LITERATURE
!  International tables for Crystallography
! TODO
!  Step by step more space groups should be installed.
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTInstallWyckoff::error = "Wyckoff positions of space group `1` with setting `2` is already installed."

GTInstallWyckoff[newwyc_] := Module[{sel, spgn, w1}, 
	Unprotect[wyckoff];
    w1   = Take[#, 1 ;; 2] & /@ wyckoff;
    spgn = newwyc[[1 ;; 2]];
    If[MemberQ[w1, spgn],
       Message[GTInstallWyckoff::error, spgn[[1]], spgn[[2]]]; 
       Return[],
       None
    ];
    wyckoff = Append[wyckoff, newwyc];
    sel = Length[wyckoff];
    Print["Wyckoff positions of ", Length[wyckoff], 
          " space groups implemented in wyckoff."];
    Protect[wyckoff];
    Return[];
]

(*
***)


 
 
 (****e* /GTGetWyckoff
! NAME
!  GTGetWyckoff
! AUTHOR
!  Wolfram Hergert
! PACKAGE
!  -> CrystalStructure.m 
! MODIFICATION HISTORY
!  11.06.2021 : first complete version, ready for tests
!  01.11.2021 : setup of documentation page
! USAGE
! GTGetWyckoff[num_,set_] retrieves Wyckoff positions of a space group with a special setting from wyckoff.
! INPUT
!  o num  - space group number
!  o set  - setting
! OUTPUT
!  Wyckoff positions of a space group with a special setting
! GTPack OPTIONS
!   * GOVerbose:
!
!     -  True   - additional information (nice table, standard)
!     -  False  - no additional information 

! Standard OPTIONS
!  -
! GTPack NOTEBOOKS 
!  Test_Wyckoff.nb
! DESCRIPTION
!  The Wyckoff positions extracted from a database are mainy used to construct structural information 
!  in the GTPack format.
!  It might be helpful to run GTAllWyckoff before to see, which settings are implemented.
! LITERATURE
!  International tables for Crystallography
! TODO
!  -
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTGetWyckoff::data = "Space group `1` With setting `2` is not installed."

GTGetWyckoff[num_, setting_, OptionsPattern[]] := Module[{verb, w1, px, dat}, 
  (*--- Clear variables used for Wyckoff definition ---*)	
    Clear[x, y, z];
  (*--- set options ---*)
    verb = OptionValue[GOVerbose];
  (* get the space group, but there might be different settings *) 
    w1   = Take[#, 1 ;; 2] & /@ wyckoff;
    If[MemberQ[w1, {num, setting}],
       w1 = Transpose[w1][[1]];
       px = Flatten[Position[w1, num]][[1]];
       dat = wyckoff[[px]],
       Message[GTGetWyckoff::data, num, setting]; 
       Return[]
    ];
    If[verb,
       GTSpaceGroups[num];
       Print["Setting  : ", dat[[2]]];
       Print["Basis    : ", dat[[3]]];
       tab = Prepend[dat[[4]], {"Multiplicity", "Wyckoff", "Symmetry", "Coordinates"}];
       Print[Grid[tab, Frame -> All, Background -> {None, 1 -> GTBackGroundColor1}, 
                  Dividers -> {{4 -> GTDividerColor1}, {2 -> GTDividerColor2}}
                 ]
       ]
    ];
    Return[dat]
]


(*
***)

(****e* /GTCoordFromWyckoff
! NAME
!  GTCoordFromWyckoff
! AUTHOR
!  Wolfram Hergert
! PACKAGE
!  -> CrystalStructure.m 
! MODIFICATION HISTORY
!  11.06.2021 : first complete version, ready for tests
!  02.11.2021 :setup of documentation page
! USAGE
! GTFromWyckoffWyckoff[num_,settings_,atoms_] generates a part of the GTPack structure file from Wyckoff positions
! INPUT
!  o num   - space group number
!  o set   - setting
!  o atoms - list of definitions  for atoms  corresponding to the  Wyckoff positions
!            atom : {Wyckoff position, atom type, substionen rules} : {"1b","O",{z-> zBo}}
! OUTPUT
!  data set for definition of structure files in GTPack
! GTPack OPTIONS
!   * GOVerbose:
!
!     -  True   - additional information (nice table, standard)
!     -  False  - no additional information 

! Standard OPTIONS
!  -
! GTPack NOTEBOOKS 
!  Test_Wyckoff.nb
! DESCRIPTION
!  This is the central command of all the commands concerning the Wyckoff positions. The information on the crystal
!  structures on more complicated system are in literature not given in the form, that we use in the package. By means 
!  of the Wyckoff positions it is possible to get the structural information in the form that is used in GTPack.
! LITERATURE
!  International tables for Crystallography
! TODO
!  -
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTCoordFromWyckoff::atoms = "Number of atoms does not agree with number of Wyckoff positions."

GTCoordFromWyckoff::wyckoff = "Wyckoff position `1` not in this space group."

GTCoordFromWyckoff[num_, setting_, atoms_, OptionsPattern[]] := 
  Module[{verb, coord, nb, pos, wyck, symb, wyp, np, nata, nat, atcoord, sort, rule, type, p, wy, co, j, k, l},
  (*--- set options ---*)
    verb  = OptionValue[GOVerbose];
  (*--- read Wyckoff positions ---*)  
    wyck  = GTGetWyckoff[num, setting, GOVerbose -> verb];
  (*--- construct the file ---*)  
    coord = wyck[[3]];
    nb    = Length[coord];
    pos   = wyck[[4]];
    symb  = Transpose[pos][[2]];
    wyp   = Transpose[pos][[4]];
    np    = Length[pos];
    nata  = Length[atoms];
    nat   = Plus @@ Transpose[pos][[1]];
    If[verb,
       Print[Grid[{{"Space group", "Basis atoms", "Wyckoff positions","Number of atoms"}, 
       	          {wyck[[1]], nb, np, nat}}, 
       	          Frame -> All,Background -> {None, 1 -> GTBackGroundColor1}
       	         ]
       ],
       None
    ];
    atcoord = {};
    Do[
       sort = atoms[[l, 2]];
       rule = atoms[[l, 3]];
       type = Last[StringPartition[atoms[[l, 1]], 1]];
       p    = Flatten[Position[symb, type]];
       If[p == {},
          Message[GTCoordFromWyckoff::wyckoff, type]; 
          Return,
          p = p[[1]]
       ]; 
       wy = wyp[[p]];
       Do[
          Do[
             co      = coord[[k]] + wy[[j]] /. rule;
             co      = co {a, b, c};
             atcoord = Append[atcoord, {co, sort}]
          , {j, 1, Length[wy]}]
       , {k, 1, nb}]
    , {l, 1, nata}];
    atcoord = atcoord /. {x -> x1, y -> y1, z -> z1};
    If[verb,
       Print[Grid[Reverse[#] & /@ atcoord, Frame -> All,Background -> { 1 -> GTBackGroundColor1, None}
       	         ]
       	   ],
       None
   ];
   Return[atcoord]
]

(*
***)


(****e* /GTSaveWyckoff
! NAME
!  GTSaveWyckoff
! AUTHOR
! W. Hergert
! PACKAGE
!  -> CrystalStructure.m 
! MODIFICATION HISTORY
!  11.06.2021 : first complete version, ready for tests
!  01.11.2021 : setup of documentation page
! USAGE
!  GTSaveWyckoff[database] saves the content of the global variable wyckoff to a database.
!  GTWriteToFile allows to use the standard name filename.wyck for the database 
! INPUT
!  database - file name
! OUTPUT
!  database to hard disc
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTWriteToFile
! GTPack NOTEBOOKS 
!  Test_Wyckoff.nb - Notebook demonstrates the use of the different functions.
! DESCRIPTION
!   In principle all is working like for the structures itself. By means of GTInstallWyckoff
!   a new structure iss installed in wyckoff. This variable is saved as the new database content
!   with GTSaveWyckoff.
! LITERATURE
!  International tables for Crystallography
! TODO
!  -
! RELEASE
! -
! PROBLEMS
! -  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSaveWyckoff[file0_] := Module[{file},
	file=StringJoin[file0,".wyck"];
    GTWriteToFile[wyckoff,file]
] 

(*
***)






  
  
(*-------------------------- Attributes ------------------------------*)

Attributes[GTLoadWyckoff]={Protected, ReadProtected}
Attributes[GTSaveWyckoff]={Protected, ReadProtected}
Attributes[GTInstallSWyckoff]={Protected, ReadProtected}
Attributes[GTAllWyckoff]={Protected, ReadProtected}
Attributes[GTGetWyckoff]={Protected, ReadProtected}
Attributes[GTCoordFromWyckoff]={Protected, ReadProtected}
Attributes[GTGroupNotation]={Listable,Protected, ReadProtected}

End[]  (* End Private Context *)



EndPackage[]
