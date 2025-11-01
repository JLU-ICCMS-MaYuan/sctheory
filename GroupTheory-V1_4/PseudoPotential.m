(****m* /Pseudopotential.m
!
! NAME
!  Pseudopotential.m
! AUTHOR
!  W. Hergert
! MODIFICATION HISTORY
!  * 02.02.2014 : first version
!  * 01.07.2018 : check mowhole package
! USAGE
!  Modules for the calculation of electronic structures by means of model pseudopotentials.
!  
! GTPack MODULES
!  
! --- Parameter Sets and Libraries ---
!
!   * GTPwDatabaseInfo		 - information on data in the pseudopotential database
!   * GTPwDatabaseRetrieve	 - retrieves a pseudopotential dataset
!   * GTPwDatabaseUpdate      - update the database, i.e. a new dataset is introduced
!   * GTPwPrintParmSet		 - pretty print of a parameter set
!
! --- Model-Pseudopotential Hamiltonians ---
!  
!   * GTPwHamiltonian		 - constructs a pseudopotential Hamiltonian
!   * GTPwDielectricF		 - screening function
!   * GTPwModelPotential     - constructs several model potentials
!   
! --- Plane Waves and Group Theory ---
!
!   * GTPwEmptyLatticeIrep   - provides Ireps for an empty lattice band structure
!   * GTPwSymmetrizePW       - finds symmetrized plane waves at a certain k-point
!   
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  A lot of development in the TB and the Photonics package apply also to the 
!  Pseudopotential case. Therefore the extension to this case is  natural. The specific
!  commands for pseudopotential theory are summarized here.
!  The commands start with GTPw#. 
! LITERATURE
! 
!
***)

BeginPackage["GroupTheory`PseudoPotential`",{"GroupTheory`RepresentationTheory`","GroupTheory`Basic`","GroupTheory`Symbols`","GroupTheory`Install`","GroupTheory`Auxiliary`","GroupTheory`Lattice`","GroupTheory`TightBinding`","GroupTheory`CrystalStructure`"}]


(* ---------- Parameter sets and libraries -----------------------------------*)
 GTPwDatabaseInfo     ::usage = "GTPwDatabaseInfo[\*StyleBox[\"database\", \"TI\"]] gives information about the pseudopotential parameter sets available in \*StyleBox[\"database\", \"TI\"]."
 GTPwDatabaseRetrieve ::usage = "GTPwDatabaseRetrieve[\*StyleBox[\"database,parameter set\", \"TI\"]] loads a pseudopotential \*StyleBox[\"parameter set\", \"TI\"] from a given \*StyleBox[\"database\", \"TI\"]."
 GTPwDatabaseUpdate   ::usage = "GTPwDatabaseUpdate[\*StyleBox[\"database\", \"TI\"]] reads a \*StyleBox[\"database\", \"TI\"], adds a new parameter set and stores the \*StyleBox[\"database\", \"TI\"]."
 GTPwPrintParmSet     ::usage = "GTPwPrintParmSet[\*StyleBox[\"database,parameter set\", \"TI\"]] prints the \*StyleBox[\"parameter set\", \"TI\"] form a pseudopotential \*StyleBox[\"database\", \"TI\"]."

(* ---------- Model-Pseudopotential Hamiltonians -----------------------------*)
 GTPwDielectricF       ::usage =  "GTPwDielectricF[\*StyleBox[\"wave vector,epsilon,Fermi wave vector,scaling\",\"TI\"]] defines a screening function."
 GTPwHamiltonian       ::usage =  "GTPwHamiltonian[\*StyleBox[\"space group,gcut,parameters\",\"TI\"]] constructs a Hamiltonian based on pseudopotential theory for a system belonging to \*StyleBox[\"space group\",\"TI\"] with a set of \*StyleBox[\"parameters\",\"TI\"] defining the pseudopotential. Reciprocal lattice vectors up to a length \*StyleBox[\"gcut\",\"TI\"] are used."
 GTPwModelPotential    ::usage =  "GTPwModelPotential[q,ionp,scrp,lfac,vfac] constructs several model potentials." 

(* ---------- Plane Waves and Group Theory -----------------------------------*)
 GTPwEmptyLatticeIrep  ::usage = "GTPwEmptyLatticeIrep[\*StyleBox[\"point group,reciprocal lattice vectors,gcut,number of bands\",\"TI\"]] determines the irreducible representations of the empty lattice band structure to a structure with \*StyleBox[\"point group\",\"TI\"]. The \*StyleBox[\"reciprocal lattice vectors\",\"TI\"] and their maximum length \*StyleBox[\"gcut\",\"TI\"] have to be provided. The analysis is done for a restricted \*StyleBox[\"number of bands\",\"TI\"]."                     
 GTPwSymmetrizePW      ::usage=  "GTPwSymmetrizePW[group,kvec] finds symmetrized plane waves" 
 
(*--------------------------- Options ----------------------------*)
Options[GTPwEmptyLatticeIrep]   = {GOVerbose -> False, GOIrepNotation -> "Mulliken",GOPhotonics->False,GOVectorField -> True}
Options[GTPwHamiltonian]        = {GOLattice -> {}, GOBondCharges -> False}
Options[GTPwSymmetrizePW]       = {GOVerbose -> True, GOProjection -> {}, GOIrepMat -> "Automatic", GOIrepNotation -> "Mulliken", GOlmax -> 15, GOMethod -> "Flodmark"}

	

	
Begin["`Private`"] (* Begin Private Context *) 


(****l* /GTPwDatabaseRetrieve
! NAME
!  GTPwDatabaseRetrieve
! AUTHOR
!  W. Hergert
! PACKAGE
!  Pseudopotential.m 
! MODIFICATION HISTORY
!  1st version August 2014
!  3.9.2014 change due to change in GTPwDatabseInfo (5th part is parameterset)
!  March 2015 check of existance of database and parmset
!  * August 2016 : check of database name and parmset fixed for subdirectories
!  * August 2016 : data base retrieve with name or index
!  * 19.05.2018  : check docu
! USAGE
!  GTPwDatabaseRetrieve[database,parameter set] loads a pseudopotential parameter set from a given database.
! INPUT
!  database      - name of the database (without ".parm")
!  parameter set - name of the parameter set
! OUTPUT
!  parameter set
! ERROR MESSAGES
!  
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  it is a copy of GTTbDatabaseRetrieve
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

GTPwDatabaseRetrieve[db_,index_]:=Module[{dbp,dbname,nset,i,pset},
	(*--- Append file extension, if necessary and check if database and parameters exist ---*)
	dbname=StringTrim[db, ".parm"] <> ".parm";
	If[!FileExistsQ[dbname],
       Print[db <> ".parm", " -> database not in working directory"];Abort[]
    ];
    (*--- Retrieve the data ---*)
	dbp=ReadList[dbname][[1]];
    nset=Length[dbp];
	(*--- Test, if name is part of data ---*)
	If[Head[index] === Integer,
       If[index > nset,         
          Print["Index = ", index, " -> only ", nset," datasets in database"]; Abort[],
          None
       ];
       pset = {dbp[[index, 2]], dbp[[index, 5]]},
       If[Position[dbp, index] == {},
          Print[index, " -> parameters not in database"]; Abort[],
          None
       ];
       Do[
          If[dbp[[i, 1]] == index,
             pset = {dbp[[i, 2]], dbp[[i, 5]]},
             None]
       , {i, 1, nset}]
    ];
    Return[pset]
]


(*
***) 


(****l* /GTPwDatabaseInfo
! NAME
!  GTPwDatabaseInfo
! AUTHOR
!  W. Hergert
! PACKAGE
!  Pseudopotential.m 
! MODIFICATION HISTORY
!  1st version August 2014
!  3.9.2014    : references added
!  1.8.2016    : index can be a number or name, column with numbers included, improved tabular print
!  March 20015 : check existance of database
!  August 2016 : check of database name and parmset fixed for subdirectories
!  19.05.2018  : ceck docu
! USAGE
!  GTPwDatabaseInfo[database] gives information about the pseudopotential parameter sets available in database.
! INPUT
!  database  - the database has to have the extension ".parm"
!  
! OUTPUT
!  database and printout of the information
! ERROR MESSAGES
!  
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  it is a copy of GTTbDatabaseInfo
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
GTPwDatabaseInfo[db_] := Module[
	{dbp, dbname,info,i,head,tab},
	(*---Append file extension, if necessary and check if database and parameters exist ---*)
     dbname = StringTrim[db, ".parm"] <> ".parm";
     If[! FileExistsQ[dbname], 
        Print["Error: the database " <> dbname <> " was not found in the working directory!"]; Abort[]
     ];
   (*---Print dtabase info---*)
     dbp  = ReadList[dbname][[1]];
     info = Transpose[Take[Transpose[dbp], 4]];
     Do[
        info[[i]] = Prepend[info[[i]], i]
     , {i, 1, Length[info]}];
     head = {"Number", "Name", "PP Type  ", "Authors", "Reference"};
     tab  = Prepend[info, head];
     Print[Grid[tab, Frame -> All, 
                     Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
                     Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, 
                     {1, 1} -> GTCornerColor}
               ]
          ];
     Return[dbp]
 ]


(*
***) 

(****l* /GTPwPrintParmSet
! NAME
!  GTPwPrintParmSet
! AUTHOR
!  W. Hergert
! PACKAGE
!  Pseudopotential.m 
! MODIFICATION HISTORY
!  1st version August 2014
! March 2015 : check database name and parmset
!  * August 2016  : check of database name and parmset fixed for subdirectories (S)
!  * August 2016  : address data with name or dataset number (W)
!  * January 2017 : the standard table form is introduced
!  * 19.05.2018   : check docu
! USAGE
!  GTPwPrintParmSet[database,parameter set] prints the parameter set form a pseudopotential database.
! INPUT
!  database        - name odf database (without ".parm")
!  parameter set   - name of the parameter set
! OUTPUT
!  printout
! ERROR MESSAGES
!  
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  formulated analogously to GTTbPrintParmSet
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

 GTPwPrintParmSet[db_, index_] := Module[
 	{dbname,dbp,nset,pset,tab,np,i,pp,tt,hl},
  (*---Append file extension,if necessary and check if database and parameters exist---*)
    dbname = StringTrim[db, ".parm"] <> ".parm";
    If[! FileExistsQ[dbname], 
       Print["Error: the database " <> dbname <> " was not found in the working directory!"]; Abort[]
       ,None
    ];
  (*---Prepare print---*)
    dbp  = ReadList[dbname][[1]];
    nset = Length[dbp];
  (*---extract data from database---*)
    If[Head[index] === Integer, 
       If[index > nset, 
          Print["Index = ", index, " -> only ", nset, " datasets in database"]; Abort[], 
          None
       ];
       pset = dbp[[index]], 
       If[Position[dbp, index] == {}, 
          Print[index, " -> parameters not in database"]; Abort[], 
          None
       ];
       Do[
       	  If[dbp[[i, 1]] == index, 
       	  	 pset = dbp[[i]], 
       	  	 None
       	  ]
       , {i, 1, nset}]
    ];
  (*---print the data---*)  
    tab = {{"Name", pset[[1]]}, 
    	   {"PP type", pset[[2]]}, 
    	   {"Source", pset[[3]]}, 
    	   {"Reference  ", pset[[4]]}
    	  };
    Print[Grid[tab, Frame -> All,Background -> {1 -> GTBackGroundColor1}]];
    If[pset[[2]] == "PPCoeff", 
       Print["reciprocal lattice vectors:"];
       np  = Length[pset[[5, 2]]];
       tab = {Table[i, {i, 1, np}], pset[[5, 2]]};
       Print[Grid[tab, Frame -> All, Background -> {None, 1 -> GTBackGroundColor1}]];
       pp  = Map[Flatten[#] &, pset[[5, 1]]];
       Print["Pseudopotential coefficients:"];
       tab = Prepend[pp, {"Atom ", Table[i, {i, 1, np}]} // Flatten];
       Print[Grid[tab, Frame -> All, Background -> {1 -> GTBackGroundColor1, 
                       1 -> GTBackGroundColor1, {1, 1} -> GTCornerColor}
                 ]
           ];
       Print["Units : ", pset[[5, 3]]], 
       None
    ];
    If[pset[[2]] === "EmptyCore" || pset[[2]] === "Heine-Abarenkov", 
       tab = {{"Dielectric constant :", pset[[5, 3]]}, 
       	      {"Volume              :", pset[[5, 4]]}, 
       	      {"Fermi vector        :", pset[[5, 5]]}
       	     };
       Print[Grid[tab, Frame -> All, Background -> {1 -> GTBackGroundColor1}]];
       pp = Map[Flatten[#] &, pset[[5, 1]]];
       Print["Parameter potential:"];
       tab = Prepend[pp, {"Atom ", "z_i", "r_i", "u"}];
       Print[Grid[tab, Frame -> All, Background -> {1 -> GTBackGroundColor1, 
                       1 -> GTBackGroundColor1, {1, 1} -> GTCornerColor}
                 ]
           ];
       Print["bond charge positions:"]; 
       np  = Length[pset[[5, 2]]];
       tt  = Prepend[pset[[5, 2]], {"\!\(\*SubscriptBox[
             StyleBox[\"\[Tau]\",\nFontWeight->\"Bold\"], \(x\)]\)", 
             "\!\(\*SubscriptBox[StyleBox[\"\[Tau]\",\nFontWeight->\"Bold\"], \(y\)]\)", 
             "\!\(\*SubscriptBox[StyleBox[\"\[Tau]\",\nFontWeight->\"Bold\"], \(z\)]\)"}];
       hl  = {" ", Table[i, {i, 1, 4}]} // Flatten;
       tab = Prepend[tt // Transpose, hl];
       Print[Grid[tab, Frame -> All, Background -> {1 -> GTBackGroundColor2, 
             1 -> GTBackGroundColor1, {1, 1} -> GTCornerColor}]]
      
    ]
  ]
     
     
(* ersetzt am 27.1.2017
  GTPwPrintParmSet[db_, index_] := Module[{dbname,dbp,nset,i,pset,np,pp}, 
	(*--- Append file extension, if necessary and check if database and parameters exist ---*)
	dbname=StringTrim[db, ".parm"] <> ".parm";
    If[!FileExistsQ[dbname], 
       Print["Error: the database " <> dbname <> " was not found in the working directory!"]; Abort[]
    ];
   (*--- Prepare print ---*)
    dbp = ReadList[dbname][[1]]; 
    nset = Length[dbp];
  (*--- extract data from database ---*) 
   If[Head[index] === Integer,
      If[index > nset,
         Print["Index = ", index, " -> only ", nset," datasets in database"]; Abort[],
         None
      ];
      pset = dbp[[index]],
      If[Position[dbp, index] == {},
         Print[index, " -> parameters not in database"]; Abort[],
         None
      ];
      Do[
         If[dbp[[i, 1]] == index,
            pset = dbp[[i]],
            None
         ]
      , {i, 1, nset}]
    ];
   (*--- print the data ---*)
   Print[{{"Name       :", pset[[1]]}, {"PP type    :", pset[[2]]}, {"Source     :", pset[[3]]},{"Reference  :", pset[[4]]}} // TableForm];
   If[pset[[2]] == "PPCoeff",
      Print["reciprocal lattice vectors:"]; 
      np = Length[pset[[5, 2]]];
      Print[TableForm[{pset[[5, 2]]}, TableHeadings -> {None, Table[i, {i, 1, np}]}]];
      pp = Map[Flatten[#] &, pset[[5, 1]]];
      Print["Pseudopotential coefficients:"];
      Print[TableForm[pp, TableHeadings -> {None, {"Atom ", Table[i, {i, 1, np}]} // Flatten}]];
      Print["Units : ", pset[[5, 3]]],
      None
   ];
   If[pset[[2]] === "EmptyCore" || pset[[2]] === "Heine-Abarenkov",
      Print[{{"Dielectric constant :", pset[[5, 3]]}, {"Volume              :", pset[[5,4]]},{"Fermi vector        :", pset[[5,5]]}} // TableForm];
      pp = Map[Flatten[#] &, pset[[5, 1]]];
      Print["Parameter potential:"];
      Print[TableForm[pp, TableHeadings -> {None, {"Atom ", "z_i", "r_i", "u"}}]];
      Print["bond charge positions:"]; np = Length[pset[[5, 2]]];
      Print[TableForm[{pset[[5, 2]]}, TableHeadings -> {None, Table[i, {i, 1, np}]}]],
      None
    ];
]
*)  

(*
***) 


(****l* /GTPwDatabaseUpdate
! NAME
!  GTPwDatabaseUpdate
! AUTHOR
!  W. Hergert
! PACKAGE
!  Pseudopotential.m 
! MODIFICATION HISTORY
!  * 1st version : August 2014
!  * March 2015  : check if database exists
!  * August 2016 : check of database name and parmset fixed for subdirectories
!  * 19.05.2015  : check docu
! USAGE
!  GTPwDatabaseUpdate[database] reads a database, adds a new parameter set and stores the database.
! INPUT
!  database  - name of the database (without ".parm")
!
!  The input of data is organized interactively, using the routine INPUT 
!  of Mathematica. The input requires the following steps:
!  Input           Type                Example
!  ----------------------------------------
!  Name           String               "GaAs"
!  PP Type        String               "PPCoeff"
!  Source         String               "Cohen&Bergstresser""
!  Input          Number               1/2
! 
! The other input depends on the PP type.
!
! OUTPUT
!  database stored and in output
! ERROR MESSAGES
!  
! GTPack MODULES
! 
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  it is buildt in the same way like GTTbDatabaseUpdate
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

GTPwDatabaseUpdate[db_] := Module[{dbname, dbp, name, str,pset, auth,file,src,it,basis,gset,ng,nam,atset,i,unit,models}, 
  models={"EmptyCore","Heine-Abarenkov","PPCoeff"};
	(*--- Append file extension, if necessary and check if database and parameters exist ---*)
	dbname=StringTrim[db, ".parm"] <> ".parm";
    If[!FileExistsQ[dbname], 
       Print["Error: the database " <> dbname <> " was not found in the working directory!"]; Abort[]
    ];
  (*--- Read database ---*)
  dbp = ReadList[dbname][[1]];
  (*--- input of general information ---*) 
  name = Input["Name of Element or Compound (String)"];
  str = Input["PP type: EmptyCore, Heine-Abarenov, PPCoeff (String)"];
  If[Intersection[models,{str}]=={},
  	 Print["Error GTPwDatabaseUpdate : Type not implemented!"];Abort[],
     None
  ];   
  auth= Input["Authors (String)"];
  src = Input["Source of data (String)"];
  it = Input["Input parameters (1), read from file (2)"];
  If[it == 1, 
     basis = Input["number of basis atoms (=1 if the atoms in a basis are equal, see Si)"];
     If[str == "PPCoeff", pset = {}; 
        gset = Input["List gvec contributions"]; ng = Length[gset];
        Do[
           nam = Input["Name of atom"]; 
           atset = Input["potential parameters (list)"];
           If[Length[atset] == ng, 
           	  None, 
              Print["Error GTPwDatabaseUpdate : Number of potential parameters and number of recip. lattice vectors do not aggree!"]; 
              Abort[]
           ];
           pset = Append[pset, {nam, atset}]
        , {i, 1, basis}]; 
        pset = {pset};
        unit = Input["Unit of parameters (eV/Ryd) (String)"]; 
        pset = Append[pset, atset];pset = Append[pset, unit],
        None
     ];
     If[str == "EmptyCore" || str == "Heine-Abarenkov", 
     	pset = {};
        Do[
           nam = Input["Name of atom"]; 
           atset = Input["z_i, r_i (a.u.), u (list)"];
           pset = Append[pset, {nam, atset}]
        , {i, 1, basis}]; 
        pset = {pset};
        atset = Input["List of bond charge positions"]; 
        pset = Append[pset, atset];
        unit = Input["Dielectric constant"]; pset = Append[pset, unit];
        unit = Input["Volume"]; pset = Append[pset, unit];
        unit = Input["Fermi vector (2Pi/a)"]; pset = Append[pset, unit],
        None
     ],
     None
  ];
  If[it == 2, 
  	file = Input["file name"];pset = GTReadFromFile[file], 
  	None
  ];
  dbp = Append[dbp, {name, str, auth, src, pset}];
  Put[dbp, dbname];
  Return[dbp]
]

(*
***) 


(****l* /GTPwDielectricF
! NAME
!  GTPwDielelctricF
! AUTHOR
!  W. Hergert
! PACKAGE
!  PseudoPotential.m
! MODIFICATION HISTORY
!   1.8.2014 : first version
! USAGE
!  GTPwDielectricF[wave vector,epsilon,Fermi wave vector,scaling] defines a screening function. (according to Hybertsen and Louie)
! INPUT
!  q    - wave vector q in units (2Pi/a)
!  kf   - Fermi energy (?) in units (2Pi/a)
!  eps0 - dielectric constant at zero fequency
!  lfac - scaling function
! OUTPUT
!  value of the screening function
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The  bare ionic potential is screened. The screening function describes this effect in q-space.
! LITERATURE
!   Hybertsen, Louie, Phys. Rev. B 37, 2733 (1988) ;
!   Levine, Louie, Phys. Rev. B 25,6310 (1987) ;
!   see also: Boehm and Unger, phys. stat. sol (b) 216, 961 (1999)
! TODO
!  not fully tested
!  documentation
! PROBLEMS
!
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTPwDielectricF[q_, kf_, eps0_, lfac_] := Module[{z, yp, ym, lamb, l2qf2, f},
  (*--- q,kf in units 2Pi/a ---*)
      z = q/kf;
      l2qf2 = Pi^2  q^2 lfac/2/kf;
      yp = z (2 + z); ym = z (2 - z);
      lamb = 2/Sqrt[3 Pi^2 lfac kf (eps0 - 1)/2];
      f = 1. + (1/2 - lamb (ArcTan[yp/lamb] + 
           ArcTan[ym/lamb])/(4 z) + (lamb^2/(16 z^3) + 1/(4 z) - 
          z/16) Log[(lamb^2 + yp^2)/(lamb^2 + ym^2)])/(l2qf2);
      Return[f]
      
  ]
            
 
(*
***)


(****l* /GTPwModelPotential
! NAME
!  GTPwDModelPotential
! AUTHOR
!  W. Hergert
! PACKAGE
!  PseudoPotential.m
! MODIFICATION HISTORY
!   4.8.2014 : first version
! USAGE
!  GTPwModelPotential defines model potential for the PP method. Implemented are the empty core 
!  model and the Heine-Abarenkov model.
! INPUT
!  q    - wave vector q in units (2Pi/a)
!  ionp - (zi,ri,u)
!  scrp - (kf,eps0)
!  lfac -  scaling function (a0/a)
!  OUTPUT
!   v(q)/eps(q)
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTPwDielectricF
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The  bare ionic potential is screened. The screening function describes this effect in q-space.
! LITERATURE
! potential:
!   Boehm, Unger, pss(b), 216, 961 (1999)
!   Fritsch et al. Phys. Rev. B 69, 165204 (2004)
!   Fritsch et al. Phys. Rev. B 67, 235205 (2003)
! screening:
!   Hybertsen, Louie, Phys. Rev. B 37, 2733 (1988) ;
!   Levine, Louie, Phys. Rev. B 25,6310 (1987) ;
!   see also: Boehm and Unger, phys. stat. sol (b) 216, 961 (1999)
! TODO
!  not fully tested
!  documentation
! PROBLEMS
!
! 
!--------------------------------------------------------------------------------
! SOURCE
*) 

GTPwModelPotential[q_, ionp_, scrp_, lfac_, vfac_, modell_] := Module[{kf, eps0, zi, ri, vg, u},
        kf = scrp[[1]]; eps0 = scrp[[2]];
        zi = ionp[[1]]; ri = ionp[[2]]*lfac; u = ionp[[3]];
        If[modell === "EmptyCore",
           vg = -2 zi*lfac*Cos[2*Pi*q*ri]/q^2/vfac,
           None
        ];
        If[modell === "Heine-Abarenkov",
           vg = -2 zi*lfac*((1 - u) Cos[2*Pi*q*ri] + u  Sin[2*Pi*q*ri]/(2 Pi q ri)) /q^2/vfac, 
           None
        ];
        vg = vg /GTPwDielectricF[q, kf, eps0, lfac]*13.6067;
        Return[vg]
]

 
(*
***)
          
            
(****l* /GTPwHamiltonian
! NAME
!  GTPwHamiltonian
! AUTHOR
!  W. Hergert
! PACKAGE
!  PseudoPotential.m
! MODIFICATION HISTORY
!   * 05.02.2014  : first version
!   * 04.08.2014  : implenetation of al parts
!   * 08.03.2015  : Help and Empty lattice work without full paramter list
!   * 17.03.2016  : units for empty lattice test changed, corresponds to Cornwell old   
!   * 01.07.2018  : check header and documentation
!                   bug fix GTPwModelPotential was called with an option, but it does not have an option
!                   option in GTPwModelPotential call removed.
! USAGE
!  GTPwHamiltonian[space group,gcut,parameters] constructs a Hamiltonian based on pseudopotential theory for a system belonging to space group with a set of parameters defining the pseudopotential. Reciprocal lattice vectors up to a length gcut are used.
! INPUT
!  * spc   - a file from the crystal strucutre database 
!  * gcut  - cutoff is the norm of the reciprocal lattice vector in units (2 Pi /a)
!  * ppdef - definition of the pseudopotentials, the list depends on the kind of emprical pseudopotentials used
! OUTPUT
!  pseudopotential Hamiltonian
! GTPack OPTIONS
!  o GOLattice  - rescales the lattice constant
!
!  o GOVerbose:
!     
!     - True   - additional information
!     - False  - no additional information
! GTPack MODULES
!  GTReciprocalBasis, 
!  GTLatCluster, 
!  GTPwModelPotential
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  
! LITERATURE
! o potential
!
!   -  Boehm, Unger, pss(b), 216, 961 (1999)
!   -  Fritsch et al. Phys. Rev. B 69, 165204 (2004)
!   -  Fritsch et al. Phys. Rev. B 67, 235205 (2003)
!
! TODO
!
! PROBLEMS
!  full meaning of kf,
!  data sets for model potentials : no need for energy units
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTPwHamiltonian[spc_, gcut_, parms_, OptionsPattern[]] := 
       Module[{ptype, lconst, combu, ufac, rlat, bas, nb, baspos, atoms, rule, gvec, blat, np, ham, gkv, i, dist, pcoeff, dg, gn, pos, k, j,
               formfac,n,models,bcharge,scrp,vfac,lfac,bondpos,nbond,rb,zb,ionp,vq,sfac,err,msg1,msg2,msg3,ppdef},
       rule = spc[[8]]; lconst = a /. rule;
       models = {"EmptyLattice", "PPCoeff", "EmptyCore", "Heine-Abarenkov"};
       (*--- error messages ---*)
       err  ="Error GTPwHamiltonian : ";
       msg1 ="PP-type not implemented, use Help!";
       msg2 =" Basis and PP definition do not aggree.";
       msg3 ="Bond charges for nb>2 not implemented!";
       ptype  = parms[[1]];
        (*--- test if model implemented ---*)
       If[ptype == "Help",
          Print["Implemented are : ", models]; Return[],
          If[Intersection[models, {ptype}] == {},
             Print[err,msg1]; Abort[],
             None
          ]
       ];
       If[ptype=="EmptyLattice",
       	 None,
         ppdef  = parms[[2]];
         bcharge = OptionValue[GOBondCharges]
       ];
       (*--- lattice information ---*) 
       If[OptionValue[GOLattice] == {},
          rule = spc[[8]],
          rule = OptionValue[GOLattice]
       ];
       rlat   = spc[[6]]; 
       bas    = spc[[7]]; nb = Length[bas]; 
       baspos = Transpose[bas][[1]];
       atoms  = Transpose[bas][[2]]; 
       If[ptype=="EmptyLattice",
       	 None, 
         If[atoms == Transpose[ppdef[[1]]][[1]],
            pcoeff = Transpose[ppdef[[1]]][[2]],
            Print[err,msg2]; Abort[]
         ];
         If[ppdef[[3]] == "eV",
            ufac = 1,  
            ufac = 13.6067
         ];
         combu = 0.511 10^6*2.4263^2/lconst^2
       ];
       rlat = rlat /. rule; baspos = baspos /. rule;
       (*--- reciprocal lattice in units 2pi/a ---*)
       blat = GTReciprocalBasis[rlat]/(2 \[Pi]);
       gvec = GTLatCluster[blat, gcut];np = Length[gvec]; 
       (*--- start with empty Hamitonian ---*)
       ham = Table[0, {np}, {np}];
  (*-----------------------------------------------*)
  (*---           Empty lattice test            ---*)    
  (*---  units : m a^2/hbar^2/pi^2              ---*)
  (*-----------------------------------------------*)
       If[ptype == "EmptyLattice",
          gkv = Map[2*(Norm[xez + #]^2) &, gvec];
          Do[
          	 ham[[i, i]] = gkv[[i]]
          , {i, 1, np}];
          Return[ham], 
          None
       ];
  (*-----------------------------------------------*)
  (*---  Pseudopotential formfactors given      ---*)
  (*-----------------------------------------------*)
       If[ptype == "PPCoeff",  
          dist = ppdef[[2]]; 
          Do[
             Do[
                dg = gvec[[i]] - gvec[[j]]; gn = Norm[dg];
                formfac = 0; pos = 0;
                Do[
                   If[Abs[dist[[n]] - gn] <= 10^(-6),
                      pos = n,
                      None
                   ]
                , {n, 1, Length[dist]}];
                If[pos == 0,
                   None,
                   Do[
                      formfac = formfac + Exp[-I 2 \[Pi] dg.baspos[[k]]]*pcoeff[[k, pos]]
                   , {k, 1, nb}]
                ];
                ham[[i, j]] = combu (gvec[[i]] + xez).(gvec[[i]] + xez) KroneckerDelta[i, j]/2 + formfac*ufac
            , {i, 1, np}]
          , {j, 1, np}];
          Return[ham],
   (*-----------------------------------------------*)
   (*---      Model potentials                   ---*)
   (*-----------------------------------------------*) 
          scrp = {ppdef[[5]], ppdef[[3]]}; vfac = ppdef[[4]]; 
          lfac = 52.9177/lconst;
          If[bcharge == True,
             If[nb == 2,
                None,
                Print[err,msg3];Abort[]
             ];
             bondpos = ppdef[[2]]; nbond = Length[bondpos];
             rb = Sqrt[(pcoeff[[1, 1]]*pcoeff[[1, 2]]^2 + 
                  pcoeff[[2, 1]]*pcoeff[[2, 2]]^2)/(pcoeff[[1, 1]] + 
                  pcoeff[[2, 1]])]*lfac;
             zb = (pcoeff[[1, 1]] + pcoeff[[2, 1]])/(4 ppdef[[3]]),
             None
          ];
          Do[
          	 ham[[i, i]] = combu (gvec[[i]] + xez).(gvec[[i]] + xez)/2
          ,{i, 1, np}];
          Do[
             Do[
                dg = gvec[[i]] - gvec[[j]];
                gn = Norm[dg]; formfac = 0;
                Do[
                   ionp = pcoeff[[k]]; 
                   vq = GTPwModelPotential[gn, ionp, scrp, lfac, vfac,ptype];
                   formfac = formfac + Exp[-I 2 Pi dg.baspos[[k]]]*vq
                , {k, 1, nb}];
                 If[bcharge,
                    ionp = {zb, rb, 0};
                    vq = GTPwModelPotential[gn, ionp, scrp, lfac, vfac, ptype];
                    sfac = 0; 
                    Do[
                       sfac = sfac + Exp[-I 2 Pi dg.bondpos[[k]]]
                    , {k, 1, nbond}]; 
                    formfac = formfac + sfac*vq*13.606, 
                    None
                 ];
                 ham[[i, j]] = formfac;
                 ham[[j, i]] = Conjugate[ham[[i, j]]]
             , {i, j + 1, np}]
          , {j, 1, np}];
          Return[ham]
      ]
  ]


(*
***)

(****l* /GTPwEmptyLatticeIrep
! NAME
!  GTPwEmptyLatticeIrep
! AUTHOR
!  W. Hergert
! PACKAGE
!  PseudoPotential.m  (Photonics.m)
! MODIFICATION HISTORY
!   March 2015: first complete version
!   03/17/2016: extension to photonic case, by option and corresponding calculation of
!   energies
!   08/25/2016: correction such that two-dimensional case are handeled appropriately
! USAGE
!  GTPwEmptyLatticeIrep[point group,reciprocal lattice vectors,gcut,number of bands] determines the irreducible representations of the empty lattice band structure
!   to a structure with point group. The reciprocal lattice vectors and their maximum length gcut have to be provided. The analysis is done for a restricted number of bands.
! INPUT
!  group  - the real space point group
!  blat   - reciprocal lattice vectors
!  kvec   - k-vector
!  gcut   - cutoff is the norm of the reciprocal lattice vector in units (2 Pi /a)
!  nenerg - number of energies to  consider
! OPTIONS
!  GOVerbose
!  GOIrepNotation
! OUTPUT
!  table with energies in the empty lattice, the degeneracy and the Ireps
! ERROR MESSAGES
!  
! GTPack MODULES
!  GTCharacterTable
!  GTLatCluster
!  GTIrep
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

GTPwEmptyLatticeIrep[group_, blat_, kvec_, gcut_, nenerg_, OptionsPattern[]] :=
     Module[{gonot,gov,phot,field,grpk,classes,characters,names,ncl,gvec,kv,el,sh,energy,mult,i,ns,chi,ll,l,
     	     tab,rep,f,test,j,k,vec,dim,mat,s,vec1,posvec,pos,charep,irp,head,ptpt,xw,yw,zw,hh},
  (*---check of input data---*)
    If[Head[group] === List, 
       None, 
       Print[group, " -> not a list of group elements"]; Abort[]
    ];
    If[Head[blat]  === List, 
       None, 
       Print[blat, " -> not a list of reciprocal lattice vectors"]; Abort[]
    ];
    If[Head[kvec]   === List, 
       None, 
       Print[kvec, " -> not a k-vector"]; Abort[]
    ];
  (*---Options---*)
    gonot = OptionValue[GOIrepNotation];
    gov   = OptionValue[GOVerbose];
    phot  = OptionValue[GOPhotonics];
    field = OptionValue[GOVectorField];
  (*---G(k) and character table---*)
    grpk = GTGroupOfK[group, kvec, blat];
    If[gov,
       Print["Group of k : ", grpk],
       None
    ];
    {classes, characters, names} =  GTCharacterTable[grpk, GOIrepNotation -> gonot, GOVerbose -> gov];
    ncl = Length[classes];
  (*---take one element of each class in element list el---*)
    el = {};
    Do[
       el = Append[el, GTGetMatrix[classes[[i, 1]]]]
   , {i, 1, ncl}];
  (*---cluster of vectors G+k ordered in shells around k---*)
   gvec = GTLatCluster[blat, gcut, GOVerbose -> gov];
   If[Length[kvec] == 3,
      kv = kvec,
      kv = Append[kvec, 0]
   ];
   gvec = Map[kv + # &, gvec];
   sh   = Take[GTLatShells[gvec], {1, nenerg}];
 (*---energies and multiplicity---*)
   mult   = Map[Length[#] &, sh];
   energy = {};
   Do[
   	  If[phot,
         energy = Append[energy, Norm[sh[[ns, 1]]]],   
         energy = Append[energy, 2*sh[[ns, 1]].sh[[ns, 1]]]
      ]
   , {ns, 1, nenerg}];
   If[gov,
   	  If[phot,
   	  	 Print["Frequenc     : ", energy],
         Print["Energy       : ", energy]
   	  ];
      Print["Multiplicity : ", mult],
      None
   ];
 (*--- discussion for scalar fields or vector fields ---*)
  tab = {};
  If[field,
 (*--- vector fields ---*)
     chi = Table[0, {ncl}];
     ll  = Map[GTGetQuaternion[#] &, el];
     ll  = Map[2 ArcCos[#[[1]]] &, ll];
     ll  = Map[2 Cos[#] &, ll];
     Do[
        If[Det[el[[l]]] == 1,
           chi[[l]] = ll[[l]],
           chi[[l]] = 0
        ]
     , {l, 1, ncl}];
  (*--- loop shells ---*) 
     Do[
        rep = Table[0, {ncl}];  
  (*--- special case k={0,0,0} ---*)       
       If[Length[sh[[i]]] == 1 && sh[[i, 1]] == {0, 0, 0},
          f[xw_, yw_, zw_] := {1, 0, 0};        
          test = GTCharProjectionOperator[classes, characters, f, {xw, yw, zw}];
          rep=If[#== {0, 0, 0},0,1]& /@ test;
          tab = Append[tab, Flatten[{energy[[i]], mult[[i]], sh[[i, 1]], rep}]],
     (*--- standard case ---*)
          Do[
             test = Map[#.sh[[i, j]] - sh[[i, j]] &, el];
             Do[
                If[test[[k]] == {0, 0, 0},
                   rep[[k]] = rep[[k]] + 1,
                   None
                ]
            , {k, 1, ncl}];
         , {j, 1, Length[sh[[i]]]}];
         rep =rep*chi;
         rep = GTIrep[rep, {classes, characters, names}, GOVerbose -> gov];
         tab = Append[tab, Flatten[{energy[[i]], mult[[i]], sh[[i, 1]], rep}]]
       ]
    , {i, 1, nenerg}],
 (*--- scalar fields ---*)
    Do[
       rep = {}; vec = sh[[i]]; dim = mult[[i]];
       Do[
       	  mat = Table[0, {dim}, {dim}]; 
       	  s = el[[k]];
          Do[
          	 vec1 = Map[s.# &, vec];
             posvec = Position[vec1, vec[[j]]]; 
             posvec = Flatten[posvec];pos = posvec[[1]]; 
             mat[[j, pos]] = 1
          , {j, 1, dim}];
          rep = Append[rep, mat]
       , {k, 1, Length[classes]}];
       charep = Map[Tr[#] &, rep];
       irp    = GTIrep[charep, {classes, characters, names}, GOVerbose -> gov];
       tab    = Append[tab, Flatten[{energy[[i]], mult[[i]], vec[[1]], irp}]]
    , {i, 1,nenerg}];
  ];
 (*--- output of the results ---*)
  If[phot,
     hh="\[Omega]",
     hh="E"
  ];
  head = Flatten[{{hh, "deg", "(k+G\!\(\*SubscriptBox[\()\), \(x\)]\)", 
                               "(k+G\!\(\*SubscriptBox[\()\), \(y\)]\)", 
                               "(k+G\!\(\*SubscriptBox[\()\), \(z\)]\)"},
                    names}];
  tab = Join[{head}, tab];
 (*---Pretty print of table---*)
  ptpt = Grid[tab, Dividers -> {{2 -> GTDividerColor1, 3 -> GTDividerColor2, 
                   6 -> GTDividerColor2}, {2 -> GTDividerColor1}}, Frame -> True, 
                   Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, None}, {{{2, nenerg + 1}, {2, 2}} -> 
                   GTCornerColor, {{2, nenerg + 1}, {3, 5}} -> GTBackGroundColor2}}
             ];
  Print[ptpt];
  Return[tab]
  ]

(* substituted 8.11.2016
GTPwEmptyLatticeIrep[group_, blat_, kvec_, gcut_, nenerg_, OptionsPattern[]] := 
   Module[{gonot,gov,grpk,classes,characters,names,el,i,gvec,np,gtsh,mult,energy,ns,tab,rep,vec,dim,
   	       mat,s,k,vec1,posvec,j,pos,charep,irp,head,ptpt,phot,kv},
   (*--- check of input data ---*) 
   If[Head[group] === List, 
   	  None, 
      Print[group, " -> not a list of group elements"]; Abort[]
   ];
   If[Head[blat] === List, 
   	  None, 
      Print[blat, " -> not a list of reciprocal lattice vectors"]; Abort[]
   ];
   If[Head[kvec] === List, 
   	  None, 
   	  Print[kvec, " -> not a k-vector"]; Abort[]
   ];
   (*--- Options ---*)
  gonot = OptionValue[GOIrepNotation];
  gov   = OptionValue[GOVerbose];
  phot  = OptionValue[GOPhotonics];
  (*--- Group of K character table ---*)
  grpk = GTGroupOfK[group, kvec, blat];
  If[gov,
     Print["Group of k : ", grpk],
     None
  ]; 
  {classes, characters, names} =GTCharacterTable[grpk, GOIrepNotation -> gonot, GOVerbose -> gov];
  (*--- take one element of each class in element list el ---*)
  el = {};
  Do[
     el = Append[el, GTGetMatrix[classes[[i, 1]]]]
  , {i, 1, Length[classes]}];
  (*--- cluster of  vectors G+k ordered in shells around k ---*)
  gvec = GTLatCluster[blat, gcut];np = Length[gvec];
  If[Length[kvec]==3,
  	 kv=kvec,
  	 kv=Append[kvec,0]
  ];
  gvec = Map[kv + # &, gvec];
  gtsh = GTLatShells[gvec];
  (*--- energies and multiplicity ---*)
  mult = Map[Length[#] &, gtsh];
  energy = {};
  Do[
  	 If[phot,
  	 	 energy = Append[energy, Norm[gtsh[[ns, 1]]]],
  	     energy = Append[energy, 2*gtsh[[ns, 1]].gtsh[[ns, 1]]]
  	 ]    
  , {ns, 1,nenerg}];
  If[gov,
     Print["Shells : ", mult];
     Print["Energy : ", energy],
     None
  ];
  (*--- construct of representations and reduction ---*)
  tab = {};
  Do[
      rep = {}; vec = gtsh[[i]]; dim = mult[[i]];
      Do[
         mat = Table[0, {dim}, {dim}]; s = el[[k]];
         Do[
            vec1 = Map[s.# &, vec];
            posvec = Position[vec1, vec[[j]]];posvec = Flatten[posvec];
            pos = posvec[[1]];mat[[j, pos]] = 1
         , {j, 1, dim}];
         rep = Append[rep, mat]
      , {k, 1, Length[classes]}];
     charep = Map[Tr[#] &, rep];
     irp = GTIrep[charep, {classes, characters, names}, GOVerbose -> gov];
     tab = Append[tab, Flatten[{energy[[i]], mult[[i]], vec[[1]], irp}]]
  , {i, 1, nenerg}]; 
  head = Flatten[{{"E", "deg", 
      "(k+G\!\(\*SubscriptBox[\()\), \(x\)]\)", 
      "(k+G\!\(\*SubscriptBox[\()\), \(y\)]\)", 
      "(k+G\!\(\*SubscriptBox[\()\), \(z\)]\)"}, names}];
  tab = Join[{head}, tab];
  (*--- Pretty print of table ---*)
  ptpt = 
  Grid[tab, 
    Dividers -> {{2 -> GTDividerColor1, 3 -> GTDividerColor2, 6 -> GTDividerColor2}, {2 ->GTDividerColor1}}, 
    Frame -> True, 
    Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, None}, {{{2, nenerg + 1}, {2, 2}} -> 
                  GTCornerColor, {{2, nenerg + 1}, {3, 5}} -> GTBackGroundColor2}}
  ]; Print[ptpt];
 (*--- return vectors in shells of k without k ---*)
  gtsh = Take[gtsh, {1, nenerg}]; 
  Do[
  	 Do[
  	 	gtsh[[i, j]] = gtsh[[i, j]] - kv
  	 , {j, 1, Length[gtsh[[i]]]}]
  , {i, 1, nenerg}];
  Return[gtsh]
]

*)
(*
***) 


(****l* /GTPwSymmetrizePW
! NAME
!  GTPwSymmetrizePW
! AUTHOR
!  W. Hergert
! PACKAGE
!  Pseudopotential.m (Photonics.m)
! MODIFICATION HISTORY
!   March 2015: first complete version
! USAGE
!  GTPwSymmetrizeWF constructs symmetrized WF for Ireps of the group of k
! INPUT
!  grpk   - group of the wave vector
!  blat   - reciprocal lattice vectors
!  kvec   - k-vector
!  shells - reciprocal lattice vectors to use for the contruction of the 
!          the Ireps
! OPTIONS
!  GOVerbose       -    False
!  GOIrepNotation  -    Mulliken
!  GOIreps         -    ivestigate all representations
!                  -    list : only listed reps will be studied
! OUTPUT
!  table with symmetrized wave functions
! ERROR MESSAGES
!  all input will be investigated if it has the correct form, otherwise the Module stops.
! GTPack MODULES
!  GTCharacterTable
!  GTGetIrep
!  GTCharacterProjectionOperator
!  GTProjectionOperator
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
(*
GTPwSymmetrizePW[grpk_, blat_, kvec_, shells_,OptionsPattern[]] :=
  Block[{gtshrep,shnum,gov,gon,rpn,nirep,characters,classes,names,rpnum,i,irmat,kn,j,vecf,repnames,results,
  	      deg,repdim,erg,f,ik,spw,tab,rv},
  rv={x,y,z};	      
  (*--- check input ---*)
  If[Head[grpk] === List,
     None,
     Print[grpk, " -> not a list of group elements"]; Abort[]
  ];
  If[Head[blat] === List,
     None, 
     Print[blat, " -> not a list of reciprocal lattice vectors"];Abort[]
  ];
  If[Head[kvec] === List,
     None, 
     Print[kvec, " -> not a vector"]; Abort[]
  ];
  If[Head[shells] === List,
     None, 
     Print[shells, " -> not a list of G-vectors"]; Abort[]
  ];
  gtshrep = Map[kvec + # &, shells]; shnum = Length[gtshrep];
  (*--- interprete the options ---*)
  gov = OptionValue[GOVerbose];
  gon = OptionValue[GOIrepNotation];
  rpn = OptionValue[GOIreps];
  (*--- calculate character table, select Ireps ---*)
  {classes, characters, names} = 
  GTCharacterTable[grpk, GOIrepNotation -> gon, GOVerbose -> gov];
  nirep = Length[characters];
  If[rpn == {},
     rpnum = Table[i, {i, 1, nirep}],
     rpnum = rpn; nirep = Length[rpn]
  ];
  (*--- calculate representation matrices ---*)
  If[gov,
     Print["Calculate the representation matrices. May be this takes some time."],
     None
  ];
  irmat = Table[0, {nirep}];
  Do[kn = rpnum[[j]];
     irmat[[j]] = GTGetIrep[grpk, kn, {classes, characters, names}]
  , {j, 1, nirep}];
  (*--- create the trial function ---*) 
  vecf = Map[Exp[I 2 \[Pi] #.{x, y, z}] &, gtshrep];
  (*--- symmetrize wavefunctions ---*) 
  Do[repnames = {}; results = {}; deg = {};
      Do[kn = rpnum[[j]]; repdim = characters[[kn, 1]];
       (*--- use character projection if possible ---*)
\
           If[repdim == 1,
                f[x_, y_, z_] = vecf[[i]]; 
     erg = GTCharProjectionOperator[classes, characters[[kn]], 
       f, rv]; c = FactorTermsList[erg][[1]];
            erg = erg*Sqrt[c]/c // FullSimplify;
               results = Append[results, erg]; deg = Append[deg, 1];
               repnames = Append[repnames, names[[kn]]],
               None
          ];
    (*--- use full projection in other cases ---*)
  \
         If[repdim > 1, 
                Do[
                    f[x_, y_, z_] = vecf[[i]]; 
      erg = GTProjectionOperator[grpk, irmat[[j]], ik, ik, 
        f, {x, y, z}] ; c = FactorTermsList[erg][[1]];
             erg = erg*Sqrt[c]/c // FullSimplify;
                    repnames = Append[repnames, names[[kn]]]; 
      results = Append[results, erg]; deg = Append[deg, repdim];
                   , {ik, 1, characters[[kn, 1]]}],
              None
        ];
    , {j, 1, nirep}];
     Print["Vector : ", shells[[i]]];
     spw = {repnames, deg, results} // Transpose;
   spw = Join[{{"Irep", "n", "Symmetrized linear combinations"}}, spw];
   tab = Grid[spw, Dividers -> {{2 -> Red, 3 -> Black}, {2 -> Red}}, 
     Frame -> True, 
     Background -> {{Pink, None}, {Pink, 
        None}, {{{2, Length[spw + 1]}, {2, 2}} -> Yellow}}];
   Print[tab]
   , {i, 1, shnum}]
  ]
*)



(****g* /GTPwSymmetrizePW
! NAME
!  GTPwSymmetrizePW
! AUTHOR
!  W. Hergert
! PACKAGE
!  PseudoPotential.m
!
! MODIFICATION HISTORY
!  *  16.2.2014 : first version
!  * 08.09.2016 : voellig ueberarbeitete Version
!  * 17.01.2017 : still some bugs, slightly changed, hopfully this helps
!  * 19.05.2018 : documentation missing , test necessary
!  * 01.07.2018 : start docuemntation, some tests
! USAGE
!  GTPwSymmetrizePW creates a symmetrized plane wave  to a given k-vector
! INPUT
!  chartab  - character table (former version group of k)
!  kv       - k-vector
!  irepn    - number  of IREP of group of k (see chartab)
!
! OUTPUT
!  table of the vectors in linear combination, used to present the result in a symbolic form
!  linear combinations corresponding to the Ireps of the group of the wave vector
! ERROR MESSAGES
!  
! GTPack MODULES
!  GTCharacterTable,
!  GTGetIrep,
!  GTProjectionOperator,
!  GTCharProjectionOperator,
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  Character projection operator and full projection operators are used
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
!  
! PROBLEMS
!   not fully tested, 
!  problem with representatin matrices? 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTPwSymmetrizePW[chartab_, kvec_, irepn_, OptionsPattern[]] := Block[
	{classes,characters,names,grk,verb,irepmat,not,plist,method,lmax,lg,
	 kv,tmat,gt,i,lg1,names1,tab,dlist,tab1,rv,subst,rxyz,dim,swf0,f,fl,swf,
	 norm,spw,spw1,fac,ll,spw2,log,pl,start,final,k,fu},
    {classes, characters, names} = chartab;
     grk     = classes // Flatten;
     lg      = GTGroupOrder[grk,GOVerbose->False];
  (*---options---*)
     verb    = OptionValue[GOVerbose];
     irepmat = OptionValue[GOIrepMat];
     not     = OptionValue[GOIrepNotation];
     plist   = OptionValue[GOProjection];
     method  = OptionValue[GOMethod];
     lmax    = OptionValue[GOlmax];
    
  (*---reciprocal vectors in symmetrized combinations---*)
     kv   = kvec;
     tmat = GTGetMatrix[Flatten[classes]];
     gt   = Table[tmat[[i]].kv, {i, 1, lg}];
     gt   = Union[gt]; 
     lg1  = Length[gt];
  (*---construct table k-vectors---*)
     If[verb, 
        names1 = {" ", "\!\(\*SubscriptBox[\(k\), \(x\)]\)", 
        "\!\(\*SubscriptBox[\(k\), \(y\)]\)", 
        "\!\(\*SubscriptBox[\(k\), \(z\)]\)"};
        tab   = Table[{Subscript["\[Psi]", i], gt[[i, 1]], gt[[i, 2]], gt[[i, 3]]}
                , {i, 1, lg1}];
        dlist = Table[{2 + i -> Black}, {i, 1, lg1}];
        dlist = Join[{{2 -> GTDividerColor1}}, dlist] // Flatten;
        tab1  = Prepend[tab, names1] // Transpose;
        tab   = Grid[tab1, Frame -> True, 
                     Dividers -> {dlist, {2 -> GTDividerColor1}}, 
                     Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, None}, 
                     {1, 1} -> GTCornerColor}
                    ];
        Print[tab],
        None
     ];
  (*---substitution rules---*)
     rv    = {xswf, yswf, zswf};
     subst = Table[Simplify[Exp[I 2 \[Pi] gt[[i]].rv]] -> Subscript["\[Psi]", i]
     	     , {i, 1, lg1}];
     rxyz  = {xswf -> Global`x, yswf -> Global`y, zswf -> Global`z};
  (*---construct symmetrized PW---*)
     dim = characters[[irepn, 1]];
     If[dim == 1,
        f[xswf_, yswf_, zswf_] = Exp[I 2 \[Pi] kv.rv];
        swf0 = GTCharProjectionOperator[classes, characters[[irepn]], 
        f, {xswf, yswf, zswf}] // Expand;
   (*---normalization,standard expression---*)
        fl   = Transpose[FactorTermsList[#] & /@ List @@ swf0][[1]];
        norm = Sqrt[Plus @@ Map[#^2 &, fl]];
        If[norm == 0,
           swf = swf0,
           swf = swf0/norm   
        ];      
        spw1 = swf /. rxyz // FullSimplify;
        spw  = FactorTermsList[swf];
        fac  = spw[[1]];
        ll   = Map[Simplify[#] &, List @@ spw[[2]]];
        spw2 = fac*(Plus @@ ll) /. subst // Expand;
        spw  = {{names[[irepn]], dim, spw2}, {" ", " ", spw1}};
        spw  = Join[{{"Irep", "n", "Symmetrized linear combinations"}}, spw];
        tab  = Grid[spw, Frame -> True, Dividers -> {{2 -> GTDividerColor1, 
                    3 -> Black}, {2 -> GTDividerColor1}}, Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, 
                    None}, {{{2, Length[spw + 1]}, {2, 2}} -> GTCornerColor}}
               ];
        Print[tab],
        None
     ];
     If[dim > 1,
        If[irepmat == "Automatic" || Length[irepmat] > 0,
           None,
           Print["Error: Dimension of Irep is ", dim, " but matrices are not provided"]; Abort[]
        ];
        If[irepmat == "Automatic",           
           irepmat = GTGetIrep[grk, irepn, {classes, characters, names}, GOlmax -> lmax, GOMethod -> method],
           None
        ];
   (*--construction of partners---*)
     If[plist == {},
        log = True; pl = Table[{i, i}, {i, 1, dim}],
        log = False; pl = plist
     ];
     If[Length[pl] == dim,
        None,
        Print["Error: number of projection operators has to be dimension of Irep"]; Abort[]
     ];
     swf = Table[0, {dim}];
     Do[
        start = pl[[k, 1]];
        final = pl[[k, 2]];
        If[log,
           f[xswf_, yswf_, zswf_] = Exp[I 2 \[Pi] kv.rv],   
           If[k == 1,
              f[xswf_, yswf_, zswf_] = Exp[I 2 \[Pi] kv.rv],  
              f[xswf_, yswf_, zswf_] = swf[[start]] // Expand
           ]
        ];
        swf[[final]] = GTProjectionOperator[grk, irepmat, start, final,f, {xswf, yswf, zswf}] // Expand
     , {k, 1, dim}];
  (*---normalization---*)
     Do[
     	fu = swf[[k]] // Expand;
        swf[[k]] = NormF[fu]
     , {k, 1, dim}];
  (*---construction of table---*)
     spw1 = swf[[1]] /. rxyz // FullSimplify;
     spw2 = swf[[1]] /. subst // Expand;
     spw  = {{names[[irepn]], 1, spw2}, {" ", " ", spw1}};
     Do[
     	spw1 = swf[[i]] /. rxyz // FullSimplify;
        spw2 = swf[[i]] /. subst // Expand;
        spw  = Append[spw, {" ", i, spw2}];
        spw  = Append[spw, {" ", " ", spw1}]
     , {i, 2, dim}];
     spw  = Join[{{"Irep", "n", "Symmetrized linear combinations"}}, spw];
     tab = Grid[spw, Frame -> True, Dividers -> {{2 -> GTDividerColor1, 3 -> Black}, {2 -> GTDividerColor1}}, 
                     Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, None}, 
                     {{{2, Length[spw + 1]}, {2, 2}} -> GTCornerColor}}
               ];
     Print[tab], None];
     Return[swf /. rxyz]
  ]



(* ersetzt am 17.1.2017

GTPwSymmetrizePW[grk_, kvec_, irepn_, OptionsPattern[]] := 
      Block[{verb, lg, kv, tmat, gt, names1, tab, dlist, tab1, classes, characters, names, pl, i, ik, ij, spw1, spw2, dim,
             irepmat, k, spwf, kl, spw,not,plist,method,lmax,lg1,log,swf0,fl,norm,rv,fac,ll,swf,subst,fu},
  (*--- options ---*)
    verb    = OptionValue[GOVerbose];
    irepmat = OptionValue[GOIrepMat];
    not     = OptionValue[GOIrepNotation];
    plist   = OptionValue[GOProjection];
    method  = OptionValue[GOMethod];
    lmax    = OptionValue[GOlmax];
    lg      = GTGroupOrder[grk, GOVerbose -> verb];
  (*--- reciprocal vectors in symmetrized combinations ---*)
    kv      = kvec;
    tmat    = GTGetMatrix[grk];
    gt      = Table[tmat[[i]].kv, {i, 1, lg}];
    gt      = Union[gt]; lg1 = Length[gt];
  (*--- construct table ---*)
    names1  = {" ", "\!\(\*SubscriptBox[\(k\), \(x\)]\)", 
                    "\!\(\*SubscriptBox[\(k\), \(y\)]\)", 
                    "\!\(\*SubscriptBox[\(k\), \(z\)]\)"
              };
    tab    = Table[{Subscript["\[Psi]", i], gt[[i, 1]], gt[[i, 2]], gt[[i, 3]]}, {i, 1, lg1}];
    dlist  = Table[{2 + i -> Black}, {i, 1, lg1}];
    dlist  = Join[{{2 -> Red}}, dlist] // Flatten;
    tab1   = Prepend[tab, names1] // Transpose;
    tab    = Grid[tab1, Frame -> True, Dividers -> {dlist, {2 -> GTDividerColor1}}, 
                        Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, None}, {1, 1} -> GTCornerColor}
                 ];
    Print[tab];
  (*--- substitution rules---*)
    rv={xswf,yswf,zswf};
    subst = Table[Simplify[Exp[I 2 \[Pi] gt[[i]].rv]] -> Subscript["\[Psi]", i], {i, 1,lg1}];
    rxyz = {xswf->Global`x,yswf->Global`y, zswf->Global`z};
  (*--- Calculate character table ---*) 
    If[verb, 
       Print["Character table group of k"], 
       None
    ];
    {classes, characters, names} = GTCharacterTable[grk, GOIrepNotation -> not, GOVerbose -> verb];
  (*--- construct symmetrized PW ---*)
    dim = characters[[irepn, 1]];
    If[dim == 1,
       f[xswf_, yswf_, zswf_] = Exp[I 2 \[Pi] kv.rv];
       swf0 = GTCharProjectionOperator[classes, characters[[irepn]], f, {xswf, yswf, zswf}]//Expand;
  (*--- normalization, standard expression ---*)
       fl   = Transpose[FactorTermsList[#] & /@ List @@ swf0][[1]]; 
       norm = Sqrt[Plus @@ Map[#^2 &, fl]]; 
       If[norm==0,
       	  swf=swf0,
       	  swf=swf0/norm
       ];
       spw1 = swf  /.rxyz // FullSimplify;
       spw=FactorTermsList[swf];fac=spw[[1]];ll=Map[Simplify[#]&,List@@spw[[2]]];
       spw2=fac*(Plus@@ll) /. subst//Expand;       
       spw  = {{names[[irepn]], dim, spw2}, {" ", " ", spw1}};
       spw  = Join[{{"Irep", "n", "Symmetrized linear combinations"}}, spw];
       tab  = Grid[spw, Frame -> True, Dividers -> {{2 -> GTDividerColor1, 3 -> Black}, {2 -> GTDividerColor1}}, 
                        Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, None}, {{{2, Length[spw + 1]}, {2, 2}} -> GTCornerColor}}
                  ]; 
       Print[tab],
       None
     ];
     If[dim > 1,
       If[irepmat == "Automatic" || Length[irepmat] > 0,
          None,
          Print["Error: Dimension of Irep is ", dim, " but matrices are not provided"]; Abort[]
       ];
       If[irepmat == "Automatic",
          irepmat = GTGetIrep[grk, irepn, {classes, characters, names}, GOlmax -> lmax, GOMethod -> method],
          None
       ];
   (*-- construction of partners ---*)
       If[plist == {},
           log = True; pl = Table[{i, i}, {i, 1, dim}],
           log = False; pl = plist
       ];
       If[Length[pl] == dim,
         None,
         Print["Error: wrong number of projection operators"]; Abort[]
       ];
       swf = Table[0, {dim}];
       Do[ 
       	  ik = pl[[k, 1]]; ij = pl[[k, 2]];
          If[log,
             kl = k; f[xswf_, yswf_, zswf_] = Exp[I 2 \[Pi] kv.rv],
             kl = ik;
             If[k == 1,
                f[xswf_, yswf_, zswf_] = Exp[I 2 \[Pi] kv.rv],
                f[xswf_, yswf_, zswf_] = swf[[ij]] // Expand
             ]
          ];          
          swf[[kl]] = GTProjectionOperator[grk, irepmat, ik, ij, f, {xswf, yswf, zswf}]//Expand;
       , {k, 1, dim}];    
   (*--- normalization ---*)
       Do[
         fu = swf[[k]] // Expand;
         swf[[k]] = NormF[fu]
       ,{k, 1, dim}];          
   (*--- construction of table ---*)
       spw1 = swf[[1]]/. rxyz // FullSimplify;
       spw2 = swf[[1]] /. subst // Expand;
       spw  = {{names[[irepn]], 1, spw2}, {" ", " ", spw1}};
       Do[
          spw1 = swf[[i]] /. rxyz // FullSimplify;
          spw2 = swf[[i]] /. subst // Expand;
          spw  = Append[spw, {" ", i, spw2}];
          spw  = Append[spw, {" ", " ", spw1}]
       , {i, 2, dim}];
       spw = Join[{{"Irep", "n", "Symmetrized linear combinations"}}, spw];
       tab = Grid[spw, Frame -> True, Dividers -> {{2 -> GTDividerColor1, 3 -> Black}, {2 -> GTDividerColor1}}, 
                                      Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, None}, {{{2, Length[spw + 1]}, {2, 2}} -> GTCornerColor}}
                 ]; 
       Print[tab],
       None
     ];
  Return[swf /.rxyz]
  ]

*)

NormF[fu_] := Module[{f1, nsu, fac, ln, ft, func, i, j, norm, nfu,fc},
      f1  = Map[Simplify[#] &, (List @@ fu)];
      f1  = List @@ # & /@ f1;
      nsu = Length[f1];
      fac = {}; func = {};
      Do[
      	 ft = f1[[i]]; ln = Length[ft];
         If[ln == 2 && ft[[1]] == E,
            fac = Append[fac, 1]; 
            func = Append[func, ft[[2]]],
            fc = 1; 
            Do[
               fc = fc*ft[[j]]
            , {j, 1, ln - 1}];
           fac = Append[fac, fc]; 
           func = Append[func, ft[[ln]]]
         ]
      , {i, 1, nsu}];
      norm = Sqrt[Plus @@ Map[#^2 &, fac]];
      nfu  = fac.func;
      If[norm == 0, 
      	Return[nfu], 
      	Return[nfu/norm // Expand]
      ]
  ]






(*
***)


End[] (* End Private Context *)

EndPackage[]

(*
***)
