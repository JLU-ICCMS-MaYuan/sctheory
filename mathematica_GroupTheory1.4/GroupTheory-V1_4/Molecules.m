(****m* /Molecules.m
!
! NAME
!  Molecules.m
! AUTHOR
!  W. Hergert
! MODIFICATION HISTORY
!  * 31.07.2016  - initially created and documented  
!  * 21.2.2017   - check of the package
!  * 28.12.2017  - check headers and documentation
!  * 25.05.2018  - check headers and documentation
! USAGE
! Contains only a few moduls to handle molecules.
! Molecules are mainly considered in the discussion of vibrations.
!
!
!  GTPack MODULES
! 
!  --- Parameter sets and libraries for molecules ---
!
!   * GTMolChemicalData      - Prints data of a molecule from Mathematica data base in  GTPack fashion.
!   * GTMolDatabaseInfo      - Information about molecules in database.
!   * GTMolDatabaseUpdate    - Add another molecule to database.
!   * GTMolGetMolecule       - Gives information about a molecule from a molecular database. 
!   * GTMolToCluster         - Transforms a molecular data set to GTPack cluster.
! 
! --- representation theory for molecules ---
!
!   * GTMolPermutationRep    - Calculats the Permutation representation
! DESCRIPTION
!   -
! LITERATURE
!   -
! TODO
!  -
! PROBLEMS
! -
!--------------------------------------------------------------------------------
!
***)

BeginPackage["GroupTheory`Molecules`",{"GroupTheory`RepresentationTheory`","GroupTheory`Basic`","GroupTheory`Symbols`","GroupTheory`Install`","GroupTheory`Auxiliary`","GroupTheory`Lattice`","GroupTheory`TightBinding`",
	"GroupTheory`CrystalStructure`","GroupTheory`PseudoPotential`","GroupTheory`Photonics`"}]
 

(*---------- Parameter sets and libraries for molecules ----------*)
 GTMolChemicalData   ::usage=   "GTMolChemicalData[\*StyleBox[\"molecule\",\"TI\"]] gives data of \*StyleBox[\"molecule\",\"TI\"] from the Mathematica database in the form used in \*StyleBox[\"GTPack\",FontWeight->\"Bold\"]."
 GTMolDatabaseInfo   ::usage=   "GTMolDatabaseInfo[\*StyleBox[\"database\",\"TI\"]] gives information about \*StyleBox[\"database\",\"TI\"], containing data of molecules."
 GTMolDatabaseUpdate ::usage=   "GTMolDatabaseUpdate[\*StyleBox[\"database\",\"TI\"]] reads \*StyleBox[\"database\",\"TI\"], adds data of another molecule and stores \*StyleBox[\"database\",\"TI\"]."
 GTMolGetMolecule    ::usage=   "GTMolGetMolecule[\*StyleBox[\"database,information\",\"TI\"]] gives \*StyleBox[\"information\",\"TI\"] of a molecule characterized from a molecular \*StyleBox[\"database\",\"TI\"]."
 GTMolToCluster      ::usage=   "GTMolToCluster[\*StyleBox[\"database,index\",\"TI\"]] transforms a molecular data set characterized by \*StyleBox[\"index\",\"TI\"] from \*StyleBox[\"database\",\"TI\"] into a cluster of atoms."
 
(*---------- representation theory for molecules -----------------*)
 GTMolPermutationRep ::usage=   "GTMolPermutationRep[data set] gives the permutation representation for a molecule defined by data set."
  
(*---------- Options ---------------------------------------------*) 	
 Options[GTMolChemicalData]    = {GOMolPlot -> True, GOVerbose -> True}
 Options[GTMolGetMolecule]     = {GOVerbose -> False}
 Options[GTMolToCluster]       = {GOMolPlot -> True, GOVerbose -> True}
 Options[GTMolPermutationRep]  = {GOFast -> GOFastValue, GOVerbose -> False}
 
Begin["`Private`"]

(****j* /GTMolPermutationRep
!  
! AUTHOR
!  W. Hergert
! NAME
!  GTMolPermutationRep
! PACKAGE
!   Molecules.m
! MODIFICATION HISTORY
!  * 13.08.2016 - first version
!  * 28.12.2017 - check headers and documentation
!  * 25.05.2018 - check headers and documentation
! USAGE
!   GTMolPermutationRep[moldat] gives the permutation representation for a molecule defined by moldat.
! INPUT
!  moldat - molecular data set
! OUTPUT
!  permutation representation
! GTPack OPTIONS
!  * GOFast  -  changes the internal varable GOFastValue 
!  * GOVerbose:
!
!      - True - permutation representation as a table 
!      - False - pretty print suppressed (Standard)
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  * GTInstallGroup
!  * GTMolToCluster
!  * GTSymmetryElementQ
!  * GTGetMatrix
! 
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


GTMolPermutationRep[moldat_, OptionsPattern[]] := 
   Module[{fast, verb, grp, gord, clu, log, vert,lmol, grp1, cc, 
           k, nl, no, eps, j, part, i, rep, tab,ld,lp,tb}, 
  (*--- parameter ---*)
     eps  = 10^(-7);lp=6;
  (*---options---*)
    fast  = OptionValue[GOFast];
     verb = OptionValue[GOVerbose];
  (*---install group---*)grp = moldat[[4]];
     grp  = GTInstallGroup[grp, GOVerbose -> False]; 
     gord = Length[grp];
  (*--cluster from molecule data---*)
     clu  = GTMolToCluster[moldat, GOMolPlot -> False, GOVerbose -> False];
  (*---check of data with respect to symmetry---*)
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
  (*---cluster data---*)
    vert = Transpose[clu][[1]]; 
    lmol = Length[vert];
  (*---construct displacement representation---*)
    grp1 = GTGetMatrix[grp];
    cc   = Flatten[Table[Map[grp1[[k]].# &, vert], {k, 1, gord}], 1];
    nl = {};
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
    rep  = Table[0, {k, 1, gord}, {i, 1, lmol}, {j, 1, lmol}]; 
    Do[
       Do[
       	  rep[[k, part[[k, i]], i]] = 1
       , {i, 1, lmol}]
    , {k, 1, gord}];
    If[verb,
       tab = Join[{grp}, {MatrixForm[#] & /@ rep}] // Transpose;
       ld = IntegerPart[gord/lp];
       If[ld==0,
       	  Print[Grid[tab//Transpose, Frame -> All,   
                                     Dividers -> {Black, {2 -> GTDividerColor1}}, 
                                     Background -> {None, {1 -> GTBackGroundColor1}}
                     ]
          ],
          Do[
          	tb = Take[tab, {(i - 1)*lp + 1, i*lp}] // Transpose;        
            Print[Grid[tb, Frame -> All,   
                            Dividers -> {Black, {2 -> GTDividerColor1}}, 
                            Background -> {None, {1 -> GTBackGroundColor1}}
                      ]
                 ]
         , {i, 1, ld}];
         tb = Take[tab, {ld*lp + 1, gord}];
         If[tb=={},
         	None,
            Print[Grid[tb//Transpose, Frame -> All, 
                        Dividers -> {Black, {2 -> GTDividerColor1}}, 
                        Background -> {None, {1 -> GTBackGroundColor1}}
                   ]
              ]
         ]
       ], 
       None
   ];
   Return[rep]
]

(*
***)



(****j* /GTMolChemicalData
!  
! AUTHOR
!  W. Hergert
! NAME
!  GTMolChemicalData
! PACKAGE
!   Molecules.m
! MODIFICATION HISTORY
!  * 31.07.2016  - first version
!  * 28.12.2017  - check headers and documentation
!  * 25.05.2018  - check headers and documentation
!  * 23.03.2021  - Molecules are sometimes represented in Mathematica without a complete set of
!                  properties. Comand is extented to handle smootly such cases.
!  * 24.03.2021  - Mathematica contains isomers: the first one is taken automatically now
!
! USAGE
!  GTMolChemicalData[molecule] gives information about a molecule from the Mathematica database.
! INPUT
!   molecule - name of a molecule or some information Mathematica can indentify.
! OUTPUT
!  general information and position of the atoms
! GTPack OPTIONS
!  * GOMolPlot:
!
!      - True  - plot of molecule (Standard)
!      - False - no plot of molecule
!  * GOVerbose:
!
!      - True  - complete information (Standard)
!      - False - information suppressed
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
!  
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTMolChemicalData::pos = "Atomic positions are not available. All set to zero.";
GTMolChemicalData::plt = "Atomic positions are not available. Plot not possible.";
GTMolChemicalData::iso = "There are isomers in Mathematica's database. The first is taken.";

GTMolChemicalData[mol_, OptionsPattern[]] := Module[{plt,verb,data,gg,ttt,pltt,data1,i},
   (*--- options ---*)
    pltt = False;
    plt  = OptionValue[GOMolPlot];
    verb = OptionValue[GOVerbose];
   (*--- data from Mathematica ---*) 
    data = { mol, 
             ChemicalData[mol, "Formula"] , 
             ChemicalData[mol, "FormulaString"] , 
             ChemicalData[mol, "BeilsteinNumber"], 
             ChemicalData[mol, "VertexTypes"], 
             ChemicalData[mol, "AtomPositions"], 
             ChemicalData[mol, "AlternateNames"]
           };
   (* check for isomers *)
    If[Length[data[[3]]]>1,
       Message[GTMolChemicalData::iso];
       data1 = Table[0, {Length[data]}]; 
       Do[
          If[Length[data[[i]]] > 1,
             data1[[i]] = data[[i, 1]], 
             data1[[i]] = data[[i]]]
       ,{i, 1, Length[data]}];
       data = data1,
       None
    ];    
   (*--- check for mising data ---*)
    data = If[Head[#] === Missing, " ", #] & /@ data;  
    If[data[[6]] == " ",
       pltt = True;
       Message[GTMolChemicalData::pos];
       data[[6]] = Table[{0, 0, 0}, {data[[5]]}],
       None
    ];      
   (*--- pretty print ---*)        
    If[verb,
       Print[{{"Name", " : ", data[[1]]},
              {"Formula", " : ", data[[2]]},
              {"FormulaString", " : ", data[[3]]},
              {"Alternate names", " : ", data[[7]]},
              {"BeilsteinNumber", " : ", data[[4]]}
             } // TableForm];
       gg  = Join[{data[[5]]}, data[[6]] // Transpose] // Transpose; 
       ttt = Join[{{"Atom", "x", "y", "z"}}, gg];
       Print[Grid[ttt, Frame -> All, 
                       Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
                       Background -> {{1 -> GTBackGroundColor1}, {1 ->  GTBackGroundColor1}, 
                       	             {1, 1} -> GTCornerColor}
                 ]
       ],
       None
    ];
    If[plt,
       If[pltt,
          Message[GTMolChemicalData::plt]; Return[],
          None
       ];
       Show[ChemicalData[mol, "MoleculePlot"]],
       pos = Table[QuantityMagnitude[data[[6, i, j]]], {i, 1, Length[data[[6]]]}, {j, 1, 3}];
       data1 = {data[[1]], data[[2]], data[[3]], " ", data[[4]], " ", data[[5]], pos};
        Return[data1]
    ]
]

(*
***)


(****j* /GTMolDatabaseInfo
!  
! AUTHOR
!  W. Hergert
! NAME
!  GTMolDatabaseInfo
! PACKAGE
!   Molecules.m
! MODIFICATION HISTORY
!  * 31.07.2016  : first version
!  * 28.12.2017  : check headers and documentation
!  * 28.12.2017  : check headers and documentation
! USAGE
!  GTMolDatabaseInfo[database] gives information about molecules in database
! INPUT
!   database - name of the database
! OUTPUT
!  general information about database entries
! GTPack OPTIONS 
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  * GTReadFromFile
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
!--------------------------------------------------------------------------------
! SOURCE
*)

GTMolDatabaseInfo[db_] := Module[{dbname,head,dbp,info,nr,tt}, 
    dbname = StringJoin[db, ".parm"];
    head = {"Number", "Name", "Formula", "FormulaString", "Symmetry", "Reference",  "Remarks"};
  (*---Check if database exists---*)
    If[Flatten[Position[FileNames[], dbname]] == {}, 
       Print["Error : database not in working directory"]; Abort[], 
       None
    ];
  (*---print info---*)
    dbp  = GTReadFromFile[dbname];
    info = Transpose[Take[Transpose[dbp], 6]];
     nr  = Table[i, {i, 1, Length[info]}];
     tt  = Join[{nr}, info // Transpose] // Transpose;
     tt  = Join[{head}, tt];
     Print[Grid[tt, Frame -> All, 
                    Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
                    Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, 
                    	           {1, 1} -> GTCornerColor
                    	          }
                ]
          ];
     Return[dbp]
 ]
 
(*
***)



(****j* /GTMolGetMolecule
!  
! AUTHOR
!  W. Hergert
! NAME
!  GTMolGetMolecule
! PACKAGE
!   Molecules.m
! MODIFICATION HISTORY
!  * 31.07.2016  : first version
!  * 28.12.2017  : check headers and documentation
!  * 25.05.2018  : check headers and documentation
! USAGE
!  GTMolGetMolecule[database,index] gives the data of a molecule from database.
! INPUT
!   * database - name of the database
!   * index    - data item to indentify the molecule in the database 
! OUTPUT
!  information about the molecule
! GTPack OPTIONS
!  * GOVerbose
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  * GTReadFromFile
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
!--------------------------------------------------------------------------------
! SOURCE
*)

GTMolGetMolecule[db_, index_, OptionsPattern[]] := Module[
	{verb,dbname,dbp,nset,mols,pset,pos,gg,ttt},
	(*--- options ---*)
     verb = OptionValue[GOVerbose]; 
    (*---Check if database and parameters exist---*)
     dbname = StringJoin[db, ".parm"];
     If[Flatten[Position[FileNames[], dbname]] == {}, 
        Print[dbname, " -> database not in working directory"]; Abort[], 
        None
     ];
   (*---read the data---*)
     dbp  = GTReadFromFile[dbname];
     nset = Length[dbp];
     mols = Take[Transpose[dbp], 1] // Flatten;
  (*---prepare the print---*)
     If[Head[index] === Integer, 
        If[index <= nset, 
           pset = dbp[[index]], 
           Print["Error: The requested parameter set is not found in the file " <> dbname <> "."]; Abort[]
        ],
        pos = Position[mols, index];
        If[Length[pos] == 0, 
           Print["Error: The requested parameter set is not found in the file " <> dbname <> "."]; Abort[], 
                 pset = dbp[[First@First@pos]]
        ]
     ];
     If[verb,
        Print[{{"Name"  , " : ", pset[[1]]},
               {"Formula", " : ", pset[[2]]},
               {"FormulaString", " : ", pset[[3]]},
               {"Symmetry", " : ", pset[[4]]},
               {"Reference", " : ", pset[[5]]},
               {"Remarks", " : ", pset[[6]]}
              } // TableForm];
       gg  = Join[{pset[[7]]}, pset[[8]] // Transpose] // Transpose;
       ttt = Join[{{"Atom", "x", "y", "z"}}, gg]; 
       Print[Grid[ttt, Frame -> All, 
                       Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
                       Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, 
                                      {1, 1} -> GTCornerColor
                                     }
                  ]
            ]
  ];
  Return[pset]
]

(*
***)


(****j* /GTMolToCluster
!  
! AUTHOR
!  W. Hergert
! NAME
!  GTMolToCluster
! PACKAGE
!   Molecules.m
! MODIFICATION HISTORY
!  * 31.07.2016  - first version
!  * 28.12.2017  - check headers and documentation
!  * 28.12.2017  - check headers and documentation
! USAGE
! GTMolToCLuster[database,index] transforms a molecular data set characterized by index from database into a cluster of atoms.
! INPUT
!  * database - molecular data base
!  * index    - data item to indentify the molecule in the database   
! OUTPUT
!   molecule as cluster data
! GTPack OPTIONS
!  * GOMolPlot:
!
!      - True  - plot of molecule (Standard)
!      - False - no plot of molecule
!  * GOVerbose:
!
!      - True  - information about min and max distances (Standard)
!      - False - information suppressed
! Standard OPTIONS
!  -
! GTPack MODULES
!  GTPlotCluster
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
!--------------------------------------------------------------------------------
! SOURCE
*)

GTMolToCluster[molecule_, OptionsPattern[]]:= Module[{plt,verb,names,pos,dist,i,j,cl,midi,madi,scale},
	(*--- options ---*)
     plt  = OptionValue[GOMolPlot];
     verb = OptionValue[GOVerbose];
    (*--- calculate distances ---*)
     names = molecule[[7]]; pos = molecule[[8]]; dist = {};
     Do[
        Do[dist = Append[dist, Norm[pos[[i]] - pos[[j]]]]
        ,{i, 1, j - 1}]
     ,{j, 1, Length[pos]}];
     cl = {pos, names} // Transpose;
     midi = Min[dist]; madi = Max[dist];
     scale = IntegerPart[madi 3];
     If[verb,
       Print[{{"minimal distance", " : ", midi},
              {"maximal distance", " : ", madi}
             } // TableForm],
       None
     ];
     If[plt,
        Return[GTPlotCluster[cl, madi + 2, scale, GOSphere -> False]],
        Return[cl]
     ]
  ]

(*
***)


(****j* /GTMolDatabaseUpdate
!  
! AUTHOR
!  W. Hergert
! NAME 
!  GTMolDatabaseUpdate
! PACKAGE
!   Molecules.m
! MODIFICATION HISTORY
!  * 31.07.2016   - first version
!  * 28.12.2017  - check headers and documentation
!  * 25.05.2018  - check headers and documentation
! USAGE
! GTMolDatabaseUpdate[database] adds information about a molecule to a database
! INPUT
!   database -  name of a molecular database
! OUTPUT
!   database with added data
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  * GTReadFromFile
!  * GTMolChemicalData
!  * GTWriteToFile
!  
! DESCRIPTION
!   The data can be found from the Mathematica database, from a file or by interactive input.
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)

GTMolDatabaseUpdate[db_] := Module[{dbname,dbp,log,mol,dataset,fname,fo,fs,sy,ref,rem,nat,csym,pos,
	                                atom},
  (*---Check if database exists---*)
   dbname = StringJoin[db, ".parm"];
   If[Flatten[Position[FileNames[], dbname]] == {}, 
      Print[dbname, " -> database not in working directory"]; Abort[], 
      None
   ];
  (*---Read database---*)
   dbp = GTReadFromFile[dbname];
   log = Input["Molecule from Mathematica (M), file (F), direct (D)? (String)"];
   If[Head[log] === String,
      None,
      Print["Error: Input is not a string"]; Abort[]
   ];
  (*--- chemical data from Mathematica ---*)
   If[log == "M",
   	  mol = Input["Molecule (String)"];
   	  If[Head[mol] === String,
         None,
         Print["Error: Name of molecule is not a string."]; Abort[]
   	  ];
      dataset = GTMolChemicalData[mol, GOMolPlot -> False, GOVerbose -> False],
      None
   ];
  (*--- input from file ---*)
   If[log == "F",
      fname = Input["Filename (String)"];
      dataset = GTReadFromFile[fname],
      None
   ];
  (*--- direct input ---*) 
   If[log == "D",
      mol  = Input["Molecule (String)"];
      fo   = Input["Formula  (String)"];
      fs   = Input["Formula string (String)"];
      sy   = Input["Symmetry(String)"];
      ref  = Input["Reference(String)"];
      rem  = Input["Remarks(String)"];
      nat  = Input["Number of atoms"]; 
      csym = {}; 
      pos  = {};
      Do[
      	 atom = Input["{chemical symbol,{x,y,z}}"];
         csym = Append[csym, atom[[1]]];
         pos  = Append[pos , atom[[2]]]
      , {i, 1, nat}];
      dataset = {mol, fo, fs, sy, ref, rem, csym, pos},
      None
   ];
  (*--- Output of database ---*)
   dbp = Append[dbp, dataset];
   GTWriteToFile[dbp, dbname];
   Return[dbp]
]

(*
***)


End[]
	 
EndPackage[]
