(*
!
! List of Changes:
!  11.06.2021   GTLoadStructure    error message changed
!  18.06.2021   GTSymmetryElementQ check of arguments at the beginning.
!  29.06.2021   GTAllStructures    change of error message
!
*)


(****m* /CrystalStructure.m  
!
! NAME
!  CrystalStructure.m
! AUTHOR
!  W. Hergert, M. Geilhufe, S. Schenk
!
! MODIFICATION HISTORY
! * 21.08.2012 : intitial documentation
! * 27.06.2018 : check of headers and documentation before release fo version 1.0.0
! USAGE
!  Basic information about point groups, space groups and crystallographic data.
!  
! GTPack MODULES
!
! --- Load, Save and Modify Crystal Structure Libraries ---
!
!  * GTLoadStructures	    - loads a list of structures
!  * GTSaveStructures	    - saves a list of structures
!  * GTClearStructures      - remove all structures from spcgrp
!  * GTInstallStructure	    - a new structure is added to the list of structures
!  * GTAllStructures        - list of all implemented structures
!  * GTGetStructure         - gives structure data file
!  * GTImportCIF            - import a cif file and change to GTPack-Format
!  * GTChangeCoord          - change of coordinates from units of lattice vectors to Cartesian
!  
! --- Manipulate and Plot Crystal Structures ---
!
!  *  GTPlotStructure	    - plots the crystal structure
!  *  GTPlotStructure2D     - plots a 2D structure
!  *  GTPlotCluster         - plots a 3D cluster
!  *  GTClusterManipulate   - manipulate structural data of a cluster
!  *  GTAtomsInCell         - finds atoms in a region of a cluster
!
! --- Lattices, Point and Space Group Information ---
!
!  *  GTBravaisLattice	    - gives the lattice vectors of the crystal system
!  *  GTCellPlot            - plots the unit cell in conventional or primitive setting
!  *  GTCrystalData         - provides information about crystal structures.
!  *  GTCrystalSystem       - gives the crystl system of a grop and vice versa
!  *  GTGroupNotation       - changes between notations
!  *  GTLatticeVectors      - generates lattice vectors
!  *  GTShowSymmetryElements- shows symmetry elemenets 
!  *  GTSymmetryElementQ    - checks if a cluster is invariant under a symmetry operation
!  *  GTSpaceGroups         - gives the nomenclature of the 230 space groups 
! 
! --- Carbon nanostructures ---
! 
!  *  GTBuckyBall           - generates a cluster of C-atoms 
!  *  GTTubeParameters      - calculatesd geometric properties of singe wall carbon nanotubes
!  *  GTTubeStructure       - creates a singel wall carbon nanotube
! 
! --- internal modules ---
!  
!  *  GTPoly                -  graphical representation of rotation axis
!  *  GTPlane               -  graphical representations of reflection planes
!  *  GTDelInv              -  delete inverse elements from a list
!  *  GTDelRots             -  selects highest order of rotation and deletes the others
!  *  GTSelectElements      -  selects minimum number of symmetry elements for graphical representation
!  *  
!      
! DESCRIPTION
!  It is not the aim of this part of GTPack to substitute all the sophisticated programs available in the community in this
!  field. There are also a lot of online resources available. This part was developed to be independent up to a certain level from 
!  external sources. Some commannds provide basic information about point groups and space groups. 
!
!  Additionally a small database is provided containing information about crystal structures which appear in connection to the electronic
!  structure calculations in the book. This database might be extended easily by the user. To store this information we use our one dataset
!  structure.For each crystal structure
!  the following information is given:
!
!	Prototype: chemical compound, name of structure
!	Pearson Symbol
!	Strukturbericht Designation
!	Space Group symbol
!	Space group number
!	Lattice Vectors
!	Basis Vectors
!	Lattice constant
!
! LITERATURE
! 
!
!  Interesting sources online are:
!
!	* http://www.cryst.ehu.es		Bilbao server
!	* http://it.iucr.org/A/		International Tables for Crystallography A: Space-group symmetry
!	* http://www.crystallography.net	Crystallography Open Database
!	* http://icsd.fiz-karlsruhe.de/	ICSD Web (Karlsruhe)
***)

BeginPackage["GroupTheory`CrystalStructure`",{"GroupTheory`Symbols`","GroupTheory`Auxiliary`","GroupTheory`Lattice`","GroupTheory`Basic`","GroupTheory`RepresentationTheory`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
 GTAllStructures	      ::usage = "GTAllStructures[] prints all currently installed structures."
 GTClearStructures	      ::usage = "GTClearStructures[] removes all currently installed structures."
 GTExportXSF              ::usage = "GTExportXSF[\*StyleBox[\"filename,data\",\"TI\"]] exports structure \*StyleBox[\"data\",\"TI\"] to .xsf format for visualisation with e.g. XCrySDen."
 GTGetStructure		      ::usage = "GTGetStructure[\*StyleBox[\"information\",\"TI\"]] gives an installed crystal structure from a certain \*StyleBox[\"information\",\"TI\"] (e.g. Pearson symbol, Strukturbericht designation, name)."
 GTInstallStructure	      ::usage = "GTInstallStructure[\*StyleBox[\"crystal structure\",\"TI\"]] adds a new \*StyleBox[\"crystal structure\",\"TI\"] to the global variable \*ButtonBox[\"spcgrp\", BaseStyle->\"Link\", ButtonData->\"paclet:GroupTheory/ref/spcgrp\"]."
 GTImportCIF		      ::usage = "GTImportCIF[\*StyleBox[\"filename\",\"TI\"]] opens .cif file and imports relevant data to generate structure data for GTInstallStructure."
 GTLoadStructures	      ::usage = "GTLoadStructures[\*StyleBox[\"database\",\"TI\"]] loads a list of installed crystal structures from a \*StyleBox[\"database\",\"TI\"]."
 GTSaveStructures	      ::usage = "GTSaveStructures[\*StyleBox[\"database\",\"TI\"]] saves the content of the list global variable \*ButtonBox[\"spcgrp\", BaseStyle->\"Link\", ButtonData->\"paclet:GroupTheory/ref/spcgrp\"] to a \*StyleBox[\"database\",\"TI\"]."
 GTChangeCoord	          ::usage = "GTChangeCoord[\*StyleBox[\"crystal structure\",\"TI\"]] changes positon of basisatoms from units of lattice vectors to Cartesian."

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)
 GTClusterManipulate	   ::usage=  "GTClusterManipulate[\*StyleBox[\"cluster\",\"TI\"]] manipulations a \*StyleBox[\"cluster\",\"TI\"] in different ways. The manipulation is controlled by the option GOMethod."
 GTPlotCluster             ::usage=  "GTPlotCluster[\*StyleBox[\"cluster,distance,scaling\",\"TI\"]] creates a plot or a graphics object of \*StyleBox[\"cluster\",\"TI\"]. \*StyleBox[\"Distance\",\"TI\"] is the maximal distance between atoms to draw a bond. \*StyleBox[\"Scaling\",\"TI\"] can be used to scale the radius of the spheres, representing the atoms. If \*StyleBox[\"cluster\",\"TI\"] is three-dimensional a Graphics3D object is generated. For a plane cluster a 2D graphic object will be constructed."
 GTPlotStructure	       ::usage = "GTPlotStructure[\*StyleBox[\"structure,Min,Max, Scale\",\"TI\"]] plots a crystal \*StyleBox[\"structure\", \"TI\"] within a certain plot range, defined by \*StyleBox[\"Min\", \"TI\"] and \*StyleBox[\"Max\", \"TI\"]. \*StyleBox[\"scale\", \"TI\"] scales the atomic radii. "
 GTPlotStructure2D         ::usage = "GTPlotStructure2D[\*StyleBox[\"struc,rc\",\"TI\"]] plots a 2D structure defined by the structure file \*StyleBox[\"struc\",\"TI\"] inside a radius \*StyleBox[\"rc\",\"TI\"]."
 GTAtomsInCell             ::usage = "GTAtomsInCell[\*StyleBox[\"object,cluster,plotdata\",\"TI\"]] finds all atoms in a region defined by \*StyleBox[\"object\",\"TI\"] in \*StyleBox[\"cluster\",\"TI\"].  \*StyleBox[\"plotdata\",\"TI\"] is used to customize the plot."
 GTCellPlot                ::usage = "GTCellPlot[\*StyleBox[\"basis, atoms, plotdata\",\"TI\"]] plots unit cells as conventional or primitive cell. \*StyleBox[\"basis\",\"TI\"] contains the lattice vectors and \*StyleBox[\"atoms\",\"TI\"] the position of the atoms in conventional setting. \*StyleBox[\"plotdata\",\"TI\"] is used to customize the plot."
(*--------- Point and Space Group Information ---------------------------*)
 GTBravaisLattice  	       ::usage = "GTBravaisLattice[\*StyleBox[\"crystal system\",\"TI\"]] gives the lattice vectors of the \*StyleBox[\"crystal system\",\"TI\"]."
 GTCrystalData             ::usage = "GTCrystalData[\*StyleBox[\"input\",\"TI\"]] provides information about crystal structures. \*StyleBox[\"Input\",\"TI\"] can be a Strukturbericht classification, or a molecular formula."
 GTCrystalSystem	       ::usage = "GTCrystalSystem[\*StyleBox[\"input\",\"TI\"]] If \*StyleBox[\"input\",\"TI\"] is a crystal system all related point groups of the system are given. If \*StyleBox[\"input\",\"TI\"] is a point group, the crystal system is determined."
 GTGroupNotation	       ::usage = "GTGroupNotation[\*StyleBox[\"group\",\"TI\"]] changes the notation of \*StyleBox[\"group\",\"TI\"] between Schoenflies and Hermann-Mauguin." 
 GTLatticeVectors          ::usage = "GTLatticeVectors[\*StyleBox[\"list\",\"TI\"]] gives the lattice vectors from a \*StyleBox[\"list\",\"TI\"] of lattice constants and lattice angles."
 GTShowSymmetryElements    ::usage=  "GTShowSymmetryElements[\*StyleBox[\"point group\",\"TI\"]] produces a graphical representation of the symmetry elements of \*StyleBox[\"point group\",\"TI\"]."
 GTSymmetryElementQ        ::usage=  "GTSymmetryElementQ[\*StyleBox[\"symmetry element, cluster\",\"TI\"]] gives True, if  \*StyleBox[\"symmetry element\",\"TI\"] leaves \*StyleBox[\"cluster\",\"TI\"] invariant and False otherwise." 
 GTSpaceGroups  	       ::usage = "GTSpaceGroups[\*StyleBox[\"argument\",\"TI\"]] gives the nomenclature of the 230 space groups." 

(*--------- Carbon nano structures ---------------------------*)
 GTBuckyBall               ::usage=  "GTBuckyBall[\*StyleBox[\"r5,r6\",\"TI\"]] generates a cluster of C atoms representing buckminsterfullerene. \*StyleBox[\"r5\",\"TI\"] denotes the bond lengths in the pentagons and \*StyleBox[\"r6\",\"TI\"] represents the length of the bonds connecting pentagons."    
 GTTubeParameters          ::usage=  "GTtubeParameters[\*StyleBox[\"n,m\",\"TI\"]] calculates the important geometric properties of single wall carbon nanotubes (\*StyleBox[\"n,m\",\"TI\"])."
 GTTubeStructure           ::usage=  "GTTubeStructure[\*StyleBox[\"n,m,cluster,ncell,ndistance\",\"TI\"]] creates a (\*StyleBox[\"n,m\",\"TI\"]) single wall carbon nanotube. \*StyleBox[\"cluster\",\"TI\"] has to be a sufficiently large part of a graphene sheet, which will be rolled up. The finite tube contains \*StyleBox[\"ncell\",\"TI\"] unit cells. \*StyleBox[\"ndistance\",\"TI\"] is an integer that sets the \*StyleBox[\"n\",\"TI\"]th distance that appears in the tube to be used as maximum bond length."
 
 (*--------------------------------------------------*)
 (*-------------------------- Internal --------------*)
 (*--------------------------------------------------*)
 
 
 (*
a        ::usage = "a Lattice constant"
b        ::usage = "b Lattice constant"
c        ::usage = "c Lattice constant"
\[Alpha] ::usage = "alpha Lattice angle"
\[Beta]  ::usage = "beta Lattice angle"
\[Gamma] ::usage = "gamma Lattice angle"

*)
Protect[a,b,c,\[Alpha],\[Beta],\[Gamma]]
(*--------------------------- Options ----------------------------*)
Options[GTAllStructures]        = {GOVerbose->False}
Options[GTChangeCoord]          = {GOVerbose->False}
Options[GTBravaisLattice]       = {GOBravais -> "Conventional", GOData -> "Data", FontSize -> 20, GOCoordinateSystem -> True, 
                                   GOImage -> {a -> 1, b -> 1.5,  c -> 2, \[Alpha] -> \[Pi]/4, \[Beta] -> \[Pi]/3, \[Gamma] -> 2 \[Pi]/6}}                                
Options[GTBuckyBall]            = {GOPlot->True,GOColorScheme -> "ElementData", GOPlot -> True, GOSphere -> True,GOCoordinateSystem->False,GODirection->{0,0,0}}
Options[GTClusterManipulate]    = {GOMethod -> {"Translate", {0,0,0}}}
Options[GTCrystalData]          = {GOTable -> False, GOStructures -> False, GOVerbose -> False}
Options[GTCrystalSystem]        = {GONotation->"SFL"}
Options[GTExportXSF]            = {GOVerbose->True}
Options[GTGetStructure]         = {GOVerbose->True}
Options[GTGroupNotation]        = {GOVerbose->True}
Options[GTImportCIF]            = {GOFractional->False,GOCorrectLabels->True,GOTolerance->10^-6}
Options[GTPlotCluster]          = {GOColorScheme->"ElementData",GOPlot->True,GOSphere->True,GOCoordinateSystem->False,GODirection->{0, 0, 0},GODimension->3,Boxed->True,GOBondColor->{}}
Options[GTPlotStructure2D]      = {PlotLabel -> {" "}, PlotStyle -> Automatic, GOBonds -> True,Frame -> False,GOLattice->{},GOTbEquivalence->True}
Options[GTShowSymmetryElements] = {GOSelectSymmetry->"All"}
Options[GTSpaceGroups]          = {GOVerbose -> True, GOTable -> False}
Options[GTTubeParameters]       = {GOLattice -> {}, GOVerbose -> False}
Options[GTTubeStructure]        = {GOVerbose -> False,GOBonds -> True,GOPlot->True}
Options[GTSymmetryElementQ]     = {GOVerbose -> False,GOPrecision->10^(-7)}
Options[GTAtomsInCell]          = {GOVerbose -> True, GOPlot -> True}
Options[GTCellPlot]             = {FontSize -> 14, GOCoordinateSystem -> True}

(*--------------------------- external Modules -------------------*)
(* only non public functions must have this 
GTOverScript = GroupTheory`Auxiliary`Private`GTOverScript*)
Begin["`Private`"] (* Begin Private Context *) 

(*--------------------------- Modules ----------------------------*)

(****e* /GTLoadStructures
! NAME
!  GTLoadStructures
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  CrystalStructure.m 
! MODIFICATION HISTORY
!  * 21.08.2012 : intitial documentation
!  * 17.02.2014 : GTReadFromFile introduced , W.
!  * 08.83.2015 : check if database in working directory
!  * 26.06.2018 : check header and documentation
!  * 11.06.2021 : error message changed
! USAGE
!  GTLoadStructures[database] loads a list of installed crystal structures from a database.
! INPUT
!  database - must have the extension ".struc"
! OUTPUT
!  the internal variabl spcgrp will be set
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTReadFromFile
! GTPack NOTEBOOKS 
!  GTSpaceGroup.nb - Notebook demonstrates the use of the different functions.
! DESCRIPTION
!  - 
! LITERATURE
!  -
! TODO
!  work on databases
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*) 

GTLoadStructures::data = "Database `1` not in working directory `2`." 

GTLoadStructures[file0_] := Module[{file,dir},
	file=StringJoin[file0,".struc"];
	(*--- check if database exists ---*)
	If[Flatten[Position[FileNames[], file]] == {},
       dir=Directory[];
	   Message[GTLoadStructures::data,file,dir];
	   Return[],
       None
    ];
	(*--- load the structures ---*)
	Unprotect[spcgrp];
	spcgrp=GTReadFromFile[file];
	Protect[spcgrp];
	Return[spcgrp]
]

(*
***)

(****e* /GTSaveStructures
! NAME
!  GTSaveStructures
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  CrystalStructure.m 
! MODIFICATION HISTORY
!  * 21.08.2012: initial documentation
!  * 17.02.2014 : use GTWriteToFile 
!  * 26.06.2018 : check header and documentation
! USAGE
!  GTSaveStructures[database] saves the content of the global variable spcgrp to a database.
!  GTWriteToFile allows to use the standard name GTPack.struc for the database 
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
!  GTSpaceGroup.nb - Notebook demonstrates the use of the different functions.
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
! -  
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSaveStructures[file0_] := Module[{file},
	file=StringJoin[file0,".struc"];
    GTWriteToFile[spcgrp,file]
] 

(*
***)


(****e* /GTClearStructures
! NAME
!  GTClearStructures
! AUTHOR
!  W.Hergert
! PACKAGE
!  SpaceGroup.m 
! MODIFICATION HISTORY
!  * 26.06.2018 : initial description and documentation
! USAGE
!  GTClearStructures[] clears the global variable spcgrp.
! INPUT
!  -
! OUTPUT
!  internal variable spcgrp={}
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  - 
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The idea came up during the ceck  of the documentation of GTSaveStructures. The documentation
!  example works better, if spcgrp is cleaned up before, otherwise it might be, that diamond is already
!  installed and one can hardly see, what happens.
! 
!  It is also useful inpther applications to clean up the structure database if a new problem will be 
!  investigated with another set of structures. Of course one could also quit the Kernel...
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
 
GTClearStructures[] := Module[{}, 
	Unprotect[spcgrp]; 
	Clear[spcgrp]; 
	spcgrp = {}; 
    Protect[spcgrp]; 
    Return[spcgrp]
] 
 
 
 (*
 ***)
 
(****e* /GTInstallStructure
! NAME
!  GTInstallStructure
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  CrystalStructure.m 
! MODIFICATION HISTORY
!  21.08.2012 : intitial documentation
!  17.02.2014 : GTGetStructure cancelled, W.
!  26.06.2018 : Check header and documentation
! USAGE
!  GTInstallStructure[crystal structure] adds a new crystal structure to the global variable spcgrp.
! INPUT
!  crystal structure - a list containing the structure information
! OUTPUT
!  spcgrp changed
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  GTSpaceGroup.nb - Notebook demonstrates the use of the different functions.
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
 


GTInstallStructure[newspc_]:=Module[{sel},
	Unprotect[spcgrp];
    If[Intersection[spcgrp,{newspc}]=={newspc},
       Print["Already implemented"],
       spcgrp=Append[spcgrp,newspc];sel=Length[spcgrp];
       If[sel==1,
          Print[Length[spcgrp]," structure is implemented in spcgrp."],
          Print[Length[spcgrp]," structures are implemented in spcgrp."]
       ]   
    ];
    Protect[spcgrp];
]

(*
***)

(****e* /GTGetStructure
! NAME
!  GTGetStructure
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  CrystalStructure.m 
! MODIFICATION HISTORY
!  * 21.08.2012 : intitial documentation
!  * 17.02.2014 : Correction of list in output not {{..}}, but{..}
!  * 01.09.2014 : Option GOVerbose introduced
!  * 26.06.2018 : check header and documentation, pretty print of table (W.)
!  * 27.06.2018 : restrict print of the first 8 entries to avoid problems with structures imported from cif-files (W.)
! USAGE
!  GTGetStructure[information] retrieves an installed crystal structure by means of an information like  Pearson symbol, Strukturbericht designation, name.
! INPUT
!  information - to get the crystal strcture
! OUTPUT
!  list containing the informatio to the crystal strucutre
! GTPack OPTIONS
!   * GOVerbose:
!
!      - True   - additional information (nice table, standard)
!     -  False  - no additional information 
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  GTSpaceGroup.nb - Notebook demonstrates the use of the different functions.
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
 
GTGetStructure[spc_, OptionsPattern[]] := 
 Module[{test, sel, selspc, actspc, prt,i},
 	prt = TrueQ[OptionValue[GOVerbose]];
  If[IntegerQ[spc] && spc <= Length[spcgrp],
    Print[" 1 structure found in spcgrp."];
    actspc=spcgrp[[spc]];
    If[prt, 
       	    Print[
       		    Grid[{{"Prototype         : ",
                       "Pearson Symbol    : ",
                       "Strukturbericht   : ",
                       "Space Group       : ",
                       "Number            : ",
                       "Lattice Vectors   : ",
                       "Basis Vectors     : ",
                       "Lattice Constants : "}, Take[actspc,{1,8}]} // Transpose, 
                Frame -> All, Background -> {{GTBackGroundColor1}, None}
           ]
        ],
        None
    ];
    Return[actspc],
    test = {spc};
    selspc = Select[spcgrp, Intersection[Flatten[#], test] == test &];
    sel = Length[selspc];
    If[sel == 0, Print[" structure not found in spcgrp."]];
    If[sel == 1, Print[" 1 structure found in spcgrp."]];
    If[sel > 1, Print[sel," structures found in spcgrp. Returning first entry."]];
    Print["----------------------------------------------------------------"];
    If[prt, 
    	      Do[actspc = selspc[[i]];
       	    Print[
       		    Grid[{{"Prototype         : ",
                       "Pearson Symbol    : ",
                       "Strukturbericht   : ",
                       "Space Group       : ",
                       "Number            : ",
                       "Lattice Vectors   : ",
                       "Basis Vectors     : ",
                       "Lattice Constants : "}, Take[actspc,{1,8}]} // Transpose, 
                Frame -> All, Background -> {{GTBackGroundColor1}, None}
           ]
        ]
       , {i, 1, sel}], None];
    If[sel == 0, Abort[], Return[selspc[[1]]]]
 ]
 ]

(*
***)

(****e* /GTAllStructures
! NAME
!  GTAllStructures
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  CrystalStructure.m 
! MODIFICATION HISTORY
!  * 21.08.2012 : intitial documentation
!  * 13.10.2015 : option GOVerbose implemente, Table of implemented structures
!  * 10.02.2016 : output of None if vb=True suppressed,
!  * 26.06.2018 : cheack headers and documentation
!  * 29.06.2021 : change of error messages
! USAGE
!  GTAllStructures prints all currently installed structures.
! INPUT
!   -
! OUTPUT
!  information about the implemented structures
! GTPack OPTIONS
!   * GOVerbose:
!
!      - True   - additional information (nice table, standard)
!     -  False  - no additional information 
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  GTSpaceGroup.nb - Notebook demonstrates the use of the different functions.
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
GTAllStructures::nostruc = "No structures implemented in spcgrp."

GTAllStructures[OptionsPattern[]] := 
   Module[{i, actspc, sel, vb,tab,grd}, sel = Length[spcgrp]; 
   tab = {{" No ", " Prototype ", " Structure ", " Strukturbericht"}};
   If[sel==0, 
   	  Message[GTAllStructures::nostruc];Return[],
   	  None
   ];
   vb = OptionValue[GOVerbose];
   If[vb,
   	  Print[sel," structures are implemented in spcgrp."],
   	  None
   ];
   Do[
   	  actspc = spcgrp[[i]];
      If[vb, 
      	 Print[
       		    Grid[{{"Prototype         : ",
                       "Pearson Symbol    : ",
                       "Strukturbericht   : ",
                       "Space Group       : ",
                       "Number            : ",
                       "Lattice Vectors   : ",
                       "Basis Vectors     : ",
                       "Lattice Constants : "}, actspc} // Transpose, 
                     Frame -> All, Background -> {{GTBackGroundColor1}, None}
                   ]
              ];
         Print["--------------------------------"], 
         tab = Append[tab, Flatten[{i, actspc[[1]], actspc[[3]]}]]
      ]
   , {i, 1, sel}];
   grd = Grid[tab, Background -> {None, {GTBackGroundColor1}}, Frame -> All,Dividers -> {{2 -> GTDividerColor1, 3 -> GTDividerColor2, 4 -> GTDividerColor2}}];
   If[!vb&&sel>0, 
   	  Print[grd]
   ]
  
  ]
  
(*
***)

(****e* /GTPlotStructure
! NAME
!  GTPlotStructure
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  CrystalStructure.m 
! MODIFICATION HISTORY
!  * 21.08.2012 : intitial documentation
!  * 18.04.2016 : take first element from Elementdata[element,"Atomicradius"] to skip the units
!  * 08.07.2016 : test implemented, Abort[] if radii incompatible with dimension of cell,
!                 scale factor for radii allows to scale the size of the atoms
!                 the correspondence of the colors to the elements is also shown.
! USAGE
!  GTPlotStructure[structure,Min,Max,scale] plots a crystal structure within a certain plot range, defined by Min and Max. Atomic radii can be scaled.
! INPUT
!  * structure - index of the structure
!  * min, max  - plot range
! OUTPUT
!  graphics object
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  - 
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  GTSpaceGroup.nb - Notebook demonstrates the use of the different functions.
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
!--------------------------------------------------------------------------------
! SOURCE
*)

GTPlotStructure[spc0_, min0_, max0_, scale_] := 
    Module[{spc, min = min0, max = max0, ers, lat, bastemp, basnr,bas, baselm$, gg, val, minval, maxval, GVec, 
    	    atomvec, k, h, i, l,test,tm,bdif},
    spc     = GTChangeCoord[spc0];
    ers     = spc[[8]];
    lat     = spc[[6]] /. ers;
    bastemp = spc[[7]] /. ers;
    basnr   = Length[bastemp];
    bas     = Table[bastemp[[i, 1]], {i, 1, basnr}];
    baselm$ = Table[bastemp[[i, 2]], {i, 1, basnr}];
  (*--- test ---*)
    test = Table[ElementData[baselm$[[i]], "AtomicRadius"][[1]], {i, 1, basnr}]; 
    tm = Min[test];
    If[tm > max0[[3]],
       Print["Lattice constant and boundaries incompatibel with atomic radii!"]; Abort[],
       None
    ];
    bdif = Union[baselm$];
    Do[
       Print[bdif[[i]], "  ", ElementData[bdif[[i]], "IconColor"]]
    , {i, 1, Length[bdif]}];
    gg     = {};
    val    = Sqrt[lat[[1]].lat[[1]]] // Simplify;
    minval = Table[3*(IntegerPart[min[[i]]/val] - 1), {i, 1, 3}];
    maxval = Table[3*(IntegerPart[max[[i]]/val] + 1), {i, 1, 3}];
    Do[
       Do[
       	  Do[
       	  	 GVec = h*lat[[1]] + k*lat[[2]] + l*lat[[3]];
             Do[
             	atomvec = GVec + bas[[i]]; 
                If[atomvec[[1]] <= max[[1]] && atomvec[[1]] >= min[[1]] && 
                   atomvec[[2]] <= max[[2]] && atomvec[[2]] >= min[[2]] && 
                   atomvec[[3]] <= max[[3]] && atomvec[[3]] >= min[[3]], 
                   None;
                   gg = Append[gg, ElementData[baselm$[[i]], "IconColor"]]; 
                   gg = Append[gg,Sphere[atomvec, ElementData[baselm$[[i]], "AtomicRadius"][[1]]*scale]], 
                   None
                ]
             , {i, 1, basnr}]
          , {l, minval[[1]], maxval[[1]]}]
       , {k, minval[[2]], maxval[[2]]}]
    , {h, minval[[3]], maxval[[3]]}];
    Graphics3D[gg]
]

(*
***)

(****e* /GTPlotStructure2D
! NAME
!  GTPlotStructure2D
! AUTHOR
!  W. Hergert
! PACKAGE
!   CrystalStructure.m
! MODIFICATION HISTORY
!  * 29.01.2015 : Moved to CrystalStructure.m
!  * 29.01.2014 : first version
!  * 28.09.2014 : option GOTbEquivalence introduced
!  * 26.06.2018 : check header and documentation  
! USAGE
!  GTPlotStructure2D[struc,rc] plots a 2D structure defined by the structure file struc inside a radius rc.
! INPUT
!  * struc    -  the usual structure file
!  *  rc       -  a cutoff radius for the cluster
! OUTPUT
!   graphics object
! GTPack OPTIONS
!
!  * GOBonds			    -	True
!  *  GOLattice	      	-	{}
!  *  GOTbEquivalence	-	True
!
! STANDARD OPTIONS   
!  * PlotLabel          -    "  "
!  * PlotStyle          -   Automatic 
!  * Frame              -   False
! GTPack MODULES
!  GTCluster
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  - 
! TODO
!  check if options do work correctly, i.e. that bonds are plotted in any case when they should plotted.
!  Seems that with GOEquivalence=True bonds are not plotted correctly.
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTPlotStructure2D[struc_, rc_, OptionsPattern[]] := 
     Module[{cl, frm, i, j, style, bas, nb, c0, c1, cl1, c, bnd, pl, d, np, g, g1, lop, equi, tbas, mult, nat, tb, sub, t, bas1},
     (*--- options ---*)	
       lop   = OptionValue[GOLattice];
       equi  = OptionValue[GOTbEquivalence];
       bnd   = OptionValue[GOBonds];
       pl    = OptionValue[PlotLabel];
       frm   = OptionValue[Frame];
       style = OptionValue[PlotStyle];
    (*--- recalculate basis atoms ---*)
       str1  = GTChangeCoord[struc];
    (*--- check for nonequivalent basis atoms ---*)
       bas = str1[[7]]; nb = Length[bas];
       If[equi === True,
          bas1=bas,
          tbas = Transpose[bas];
          mult = Tally[tbas[[2]]]; nat = Length[mult];
          tb = {}; sub = {};
          Do[
             t = Table[mult[[i, 1]] <> ToString[j], {j, 1, mult[[i, 2]]}];
             tb = Append[tb, t] // Flatten
          , {i, 1, nat}];
          tbas[[2]] = tb; bas1 = Transpose[tbas]
       ];
     (*--- construction of the cluster ---*)  
       cl = GTCluster[str1, rc, GOLattice -> lop, GOTbEquivalence -> equi];
     (*--- prepare the plot ---*)
       bas = Transpose[bas1][[2]]//Union;
       nb = Length[bas];
       c = {};
       Do[
       	  c0 = Transpose[Select[cl, #[[2]] == bas[[i]] &]][[1]];
          c1 = Take[Transpose[c0], {1, 2}] // Transpose;
          c = Append[c, c1]
       , {i, 1, nb}];
       If[bnd, 
       	  cl1 = Flatten[c, 1]; np = Length[cl1];
          d = Table[Table[Norm[cl1[[i]] - cl1[[j]]], {i, 1, j - 1}], {j, 1, np}];
          d=d //Flatten // Union // Min; g = {};
          Do[
          	 Do[
          	 	If[Norm[cl1[[i]] - cl1[[j]]] <= d + 0.01, 
                   g = Append[g, Line[{cl1[[i]], cl1[[j]]}]], 
                   None
                ]
             , {i, 1, j - 1}]
          , {j, 1, np}], 
          None
        ];
     (*--- plot ---*)
        g1 = ListPlot[c, Prolog -> PointSize[0.05], PlotStyle -> style, 
                         AspectRatio -> 1, Frame -> frm, PlotLabel -> pl, 
                        PlotLegends -> 
                        PointLegend[Automatic, bas, LegendFunction -> "Frame", 
                        LegendLabel -> "Atoms"]
             ];
         If[bnd, 
         	Show[g1, Graphics[g]], 
         	Show[g1]
         ]
   ]

(*
***)


(****e* /GTCrystalSystem
! NAME
!  GTCrystalSystem
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!  * 22.01.2013 : first version
!  * 04.05.2013 : revised version, with options 
!  * 30.05.2015 : generalization of input with GTSwitchNotationIntern 
!  * 15.01.2016 : reformulation with respect to new namenclature
!  * 23.02.2016 : GTSwitchNotation -> GTGroupNotation replaced
!  * 27.06.2016 : check header and documentation
! USAGE
!  GTCrystalSystem[input] returns  all related point groups, if input is a crystal system If
!  input is a point group, the crystal system is determined.
! INPUT
!  input - can be either a crystal system or a point group (SFL or HM). Pointgroups have to be in symbolic form.
! OUTPUT
!  * groups of a crystal system
!  * crystal system to which the group belongs 
! GTPack OPTIONS
!
!   GONotation : 
!
!     "SFL"   - Schoenfliess notation
!      "HM"   - Herman-Maughuin notation
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTGroupNotation
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  It is based on lists, i.e. no complex alogithm necessary.
! LITERATURE
!  Mathiak/Stingl , Ludwig/Falter 
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

GTCrystalSystem[grp0_, OptionsPattern[{GONotation -> "SFL"}]] := 
 Module[{x, nt, cs, i, t, j,namesdbl,grp}, 
 	nt = OptionValue[GONotation];
    cs = {{"Triclinic", "Monoclinic", "Orthorhombic", "Tetragonal",  "Trigonal", "Hexagonal",  "Cubic"}, 
    	  {{C1, Ci}, {C2, Cs, C2h}, {D2, C2v, D2h}, {C4, S4, C4h,
            D4, C4v, D2d, D4h}, {C3, C3i, D3, C3v, D3d}, {C6, C3h, C6h, D6,
            C6v, D3h, D6h}, {T, Th, Td, O, Oh}}
         };
    namesdbl = {{C1h, S2, S6, S1, V, Vd, Vh}, {Cs, Ci, C3i, Cs, D2, D2d, D2h}};
    (*--- resolve double notation ---*)
    grp = grp0;
    If[Head[grp] == Symbol,
       Do[
          If[grp === namesdbl[[1, i]],
             y = namesdbl[[2, i]]; Print[grp, " is ", y];
             grp = y,
             None
          ]
       , {i, 1, 7}]
    ];
    (*--- switch HM to SFL ---*)
    If[Head[grp] === BracketingBar,
       x = GTGroupNotation[grp, GOVerbose -> False],
       (*x = GTSwitchNotation[grp], *)
       x = grp
    ];
    (*--- pointgroups to crystal system ---*)
    If[Head[grp] === String,
       If[Position[cs[[1]], grp] == {},
          Print["Error: GTCrystalSystem: Wrong sytem name or pointgroup not a symbol"]; Abort[],
         None
       ];
       Do[
       	  If[x == cs[[1, i]], 
             If[nt == "HM",
                t = Map[GTGroupNotation[#, GOVerbose -> False] &, cs[[2, i]]],
                t = cs[[2, i]]
             ],
             None
          ]
       , {i, 1, 7}];
       Print["Crystal system - ", x, " : ", t],
       Do[
          Do[
             If[x == cs[[2, i, j]],
                If[nt == "HM",
                   Print[GTGroupNotation[x, GOVerbose -> False], "(", x,   
                   ") belongs to the crystal system: ", cs[[1, i]]],
                   Print[x, "(", GTGroupNotation[x, GOVerbose -> False],   
                   ") belongs to the crystal system: ", cs[[1, i]]]
                ],
                None 
             ]
         , {j, 1, Length[cs[[2, i]]]}]
        , {i, 1, 7}]
    ]
  ]

(*
***)


(****e* /GTGroupNotation
! NAME
!  GTGroupNotation
! AUTHOR
!  W. Hergert
! PACKAGE
!   CrystalStructure.m
! MODIFICATION HISTORY
!  * 22.01.2013 : first version
!  * 04.05.2013 : modified, use of options, better HM notation
!  * 23.06.2015 : HM mit bracketingbar
!  * 15.01.2016 : all notations corresponding to the new scheme
!  * 27.06.2018 : check header
! USAGE
!  GTGroupNotation[group] changes the notation of group between Schoenflies and Hermann-Mauguin.
! INPUT
!  group name
! OUTPUT
!  group name
! ERROR MESSAGES
! GTPack OPTIONS
!  GOVerbose:
!
!    - True   - additional information if two names for the same group can be used (standard)
!    - False  - no additional information
! STANDARD OPTIONS 
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!   -
! DESCRIPTION
!  It is based on lists, i.e. no complex alogithm necessary.
! LITERATURE
!  Mathiak/Stingl , Ludwig/Falter 
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

GTGroupNotation[grp_, OptionsPattern[]] := 
 Module[{x, y, namesdbl, names, i, sw1, sw2},
 	(*--- Definition of the names and Symbols ---*)
   y="nn";
  namesdbl = {{C1h, S2, S6, S1, V, Vd, Vh}, {Cs, Ci, C3i, Cs, D2, D2d,
      D2h}};
  names = {{Oh, O, Th, Td, D6h, D4h, T, C6h, C6v, D6, D3h, D3d, C4h, 
     C4v, D4, D2h, D2d, C6, C3i, C3h, C3v, D3, C4, S4, C2h, C2v, D2, 
     Cs, Ci, C2, C1, 
     C3}, {
     \[LeftBracketingBar]m-3m\[RightBracketingBar],
     \[LeftBracketingBar]432\[RightBracketingBar],
     \[LeftBracketingBar]m-3\[RightBracketingBar],
     \[LeftBracketingBar]-43m\[RightBracketingBar],
     \[LeftBracketingBar]6/mmm\[RightBracketingBar],
     \[LeftBracketingBar]4/mmm\[RightBracketingBar],
     \[LeftBracketingBar]23\[RightBracketingBar],
     \[LeftBracketingBar]6/m\[RightBracketingBar],
     \[LeftBracketingBar]6 mm\[RightBracketingBar],
     \[LeftBracketingBar]622\[RightBracketingBar],
     \[LeftBracketingBar]-6m2\[RightBracketingBar],
     \[LeftBracketingBar]-3m\[RightBracketingBar],
     \[LeftBracketingBar]4/m\[RightBracketingBar],
     \[LeftBracketingBar]4 mm\[RightBracketingBar], 
     \[LeftBracketingBar]422\[RightBracketingBar],
     \[LeftBracketingBar]mmm\[RightBracketingBar],
     \[LeftBracketingBar]-42m\[RightBracketingBar],
     \[LeftBracketingBar]6\[RightBracketingBar],
     \[LeftBracketingBar]-3\[RightBracketingBar], 
     \[LeftBracketingBar]-6\[RightBracketingBar],
     \[LeftBracketingBar]3m\[RightBracketingBar],
     \[LeftBracketingBar]32\[RightBracketingBar],
     \[LeftBracketingBar]4\[RightBracketingBar],
     \[LeftBracketingBar]-4\[RightBracketingBar],
     \[LeftBracketingBar]2/m\[RightBracketingBar],
     \[LeftBracketingBar]mm2\[RightBracketingBar],
     \[LeftBracketingBar]222\[RightBracketingBar],
     \[LeftBracketingBar]m\[RightBracketingBar],
     \[LeftBracketingBar]-1\[RightBracketingBar],
     \[LeftBracketingBar]2\[RightBracketingBar],
     \[LeftBracketingBar]1\[RightBracketingBar],
     \[LeftBracketingBar]3\[RightBracketingBar]}};
  If[Head[grp] === Symbol,
       sw1 = 1; sw2 = 2,
       If[Head[grp] === BracketingBar,
            sw1 = 2; sw2 = 1,
            Print["Error: GTGroupNotation: Symbol not found!"]; Abort[]
     ]
   ];
  x = grp;
  If[Head[grp] === Symbol,
              Do[
                   If[x === namesdbl[[sw1, i]],
                        y = namesdbl[[sw2, i]];
                        If[OptionValue[GOVerbose],
                             Print[x, " is ", y],
                            None
                       ];
                        x = y; y = "nn",
                        None
                ]
           , {i, 1, 7}]
   ];
  Do[
       If[x === names[[sw1, i]],
             y = names[[sw2, i]],
             None
     ]
   , {i, 1, 32}];
  Return[y]
  ]

(*
***)

(****e* /GTCrystalData
! NAME
!  GTCrystalData
! AUTHOR
!  W. Hergert
! PACKAGE
!   CrystalStructure.m
! MODIFICATION HISTORY
!  * 01.02.2014 : first version
!  * 15.12.2015 : totally new version
!  * 13.01.2016 : version completed, data filled in
!  * 27.06.2018 : check header and documentation
! USAGE
!  GTCrystalData[input] provides information about crystal structures. argument can be a Strukturbericht 
!  nomenclature, or a molecular formula
! INPUT
!  input -  Strukturbericht nomenclature, molecular formula 
!  
! OUTPUT
!  - if a single space group is in the INPUT all the different names and an information about
!    symmorphic/nonsymmorphic are printed and also given in  list.
!  - if the crystal system is INPUT, all space groups which belog to the crystal system are printed.
! ERROR MESSAGES
!
! GTPack OPTIONS
!  GOTable: 
!    - False (standard) - input interpreted als Struktubericht information, molecular formula or spacegroup number
!    - True             - Table with implemented information is printed
!  GOStructures: 
!    - False (standard) - input is Strukturbericht information         
!    - True             - input is molecular formula 
!  GOVerbose: 
!    - False (Standard)  - no additional information
!    - True              - additional information about the space group, and structural examples
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTSpaceGroups
! GTPack NOTEBOOKS 
!  Wolfram_Devel/4_Discrete/Symmetry_Groups/GTCrystalData.nb
! DESCRIPTION
!  In principle it is a database consisting of two lists. The list table contains information to 
!  the Strukturbericht classification. The list structures contains examples of crystal structures.
! LITERATURE
!  * Landolt-Boernstein 1. Band, Atom- und Molekuelphysik, 4. Teil Kristalle, Springer 1955
!  * W.B. Pearson, A Handbook of Lattice Spacings and Structures of Metals and Alloys, Pergamon Press,
!  * G.S. Rohrer, Structure and Bonding in Crystalline Materials, Cambridge University Press,2001
!  * F. Cardelli, Materials Handbook, A Concise Desktop Reference, Springer, 2000
! TODO
!  Has to be checked carefully !!
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTCrystalData::struc = "Classification or structure not found.";

GTCrystalData[struc_, OptionsPattern[]] := Module[{head, tab, row, pos,table,comp,verb,structures,info,spc,ttt},
  head = {"No.", "Strukturbericht", "Citation", "Prototype", 
    "Pearson symbol", "Space group", "GTPack.struc"} ; 
    
  (*--- options ---*)
  comp = OptionValue[GOStructures];
  verb = OptionValue[GOVerbose];
  tab = OptionValue[GOTable];
  If[verb,
  	Print["Citation (v,p) is related to Strukturbericht, volume v and page p."],
  	None
  ];
  (*--- Strukturbericht table ---*)
  table = {
    {1, "A1", "(1,13)", "Cu", "cF4", "225", "*"},
    {2, "A2", "(1,15)", "W", "cI2", "229", " "},
    {3, "A3", "(1,16)", "Mg", "hP2", "194", " "},
    {4, "A4", "(1,19)", "C (Diamond)", "cF8", "227", " "},
    {5, "A5", "(1,21)", "Sn (White Tin)", "tI4", "141", " "},
    {6, "A6", "(1,23)", "In", "tI2", "139", " "},
    {7, "A7", "(1,25)", "As", "hR2", "166", " "},
    {8, "A8", "(1,27)", "Se", "hP3", "152,154", " "},
    {9, "A9", "(1,28)", "C (Graphite)", "hP4", "194", " "},
    {10, "A10", "(1,737)", "Hg", "hR1", "166", " "},
    {11, "A11", "(2,1;3,3)", "Ga", "oC8", "64", " "},
    {12, "A12", "(2,2)", "\[Alpha]-Mn type", "cI58", "217", " "},
    {13, "A13", "(2,3)", "\[Beta]-Mn", "cP20", 213, " "},
    {14, "A14", "(2,5)", "J", "cP8", "64", " "},
    {15, "A15", "(2,6)", "\!\(\*SubscriptBox[\(W\), \(3\)]\)O", "cP8",
      223, " "},
    {16, "A16", "(3,4)", "\[Alpha]-S", "oF128", "70", " "},
    {17, "A17", "(3,6)", "P (black)", "oC8", "64", " "},
    {18, "A18", "(4,4)", "Cl", " ", "138", " "},
    {19, "A19", "(4,4)", "Po", " ", "5", " "},
    {20, "A20", "(6,3)", "\[Alpha]-U", "oC4" , "63", " "},
    {21, "B1", "(1,72)", "NaCl", "cF8", "225", "*"},
    {22, "B2", "(1,74)", "CsCl", "cP2", "221", "*"},
    {23, "B3", "(1,76)", "ZnS (Sphalerite)", "cF8", "216", "*"},
    {24, "B4", "(1,78)", "ZnO (Wurtzite)", "hP4", "186", "*"},
    {25, "\!\(\*SubscriptBox[\(B8\), \(1\)]\)", "(1,84)", "NiAs", 
     "hP4", "194", " "},
    {26, "\!\(\*SubscriptBox[\(B8\), \(2\)]\)", "(1,84)", 
     "\!\(\*SubscriptBox[\(Ni\), \(2\)]\)In", "hP6", "194", " "},
    {27, "B9", "(1,87)", "HgS (Cinnabar)", "hP6", "152,154", " "},
    {28, "B10", "(1,89)", "PbO", "tP4", "129", " "},
    {29, "B11", "(1,94)", "\[Gamma]-CuTi", "tP4", "129", " "},
    {30, "B12", "(1,95)", "BN", "hP4", "194"},
    {31, "B13", "(2,6)", "NiS (Millerite)", "hR6", "160", " "},
    {32, "B16", " ", "GeS", "oP8", "62", " "},
    {33, "B17", "(2,9)", "PtS (Cooperite)", "tP4", "131" , " "},
    {34, "B18", "(2,10)", "CuS (Covellite)", "hP12", "194", " "},
    {35, "B19", "(2,11)", "AuCd", "oP4 ", "51", " "},
    {36, "B20", "(2,13)", "FeSi", "cP8", "198", " "},
    {37, "B27", "(3,12)", "FeB", "oP8", "62", " "},
    {38, "B29", "(3,14)", "SnS", "oP8", "62", " "},
    {39, "B31", "(3,17)", "MnP", "oP8", "62", " "},
    {40, "B32", "(3,19)", "NaTl", "cF16", "227", " "},
    {41, "B34", " " , "PdS", "tP16", "84", " "},
    {42, "B35", "(6,4)", "CoSn", "hP6", "191", " "},
    {43, "B37", " ", "TlSe", "tI16", "140", " "},
    {44, "\!\(\*SubscriptBox[\(B\), \(e\)]\)", " ", "CdSb", "oP16", 
     "61", " "},
    {45, "\!\(\*SubscriptBox[\(B\), \(f\)]\)", " ", "\[Zeta]-CrB", 
     "oC8", "63", " "},
    {46, "\!\(\*SubscriptBox[\(B\), \(g\)]\)", " ", "MoB", "tI16", 
     "141", " "},
    {47, "\!\(\*SubscriptBox[\(B\), \(h\)]\)", " ", "WC", "hP2", 
     "191", " "},
    {48, "\!\(\*SubscriptBox[\(B\), \(i\)]\)", " ", "\[Gamma]'-MoC", 
     "hP8", "194", " "},
    {49, "C1", "(6,4)", 
     "\!\(\*SubscriptBox[\(CaF\), \(2\)]\) (Fluorite)", "cF12", "225",
      " "},
    {50, "\!\(\*SubscriptBox[\(C1\), \(b\)]\)", "(6,4)", "MgAgAs", 
     "cF12", "216", " "},
    {51, "C2", "(1,150)", 
     "\!\(\*SubscriptBox[\(FeS\), \(2\)]\) (Pyrite)", "cF12", "205", 
     " "},
    {52, "C3", "(1,153)", 
     "\!\(\*SubscriptBox[\(Cu\), \(2\)]\)O (Cuprite)", "cP6", "224", 
     " "},
    {53, "C4", "(1,155)", 
     "\!\(\*SubscriptBox[\(TiO\), \(2\)]\) (Rutile)", "tP6", "136", 
     " "},
    {54, "C5", "(1,158)", 
     "\!\(\*SubscriptBox[\(TiO\), \(2\)]\)(Anatase)", "tI6", "141" , 
     " "},
    {55, "C6", "(1,161)", "\!\(\*SubscriptBox[\(CdJ\), \(2\)]\)", 
     "hP3", "156", " "},
    {56, "C7", "(1,164)", 
     "\!\(\*SubscriptBox[\(MoS\), \(2\)]\) (Molybdenite)", "hP6", 
     "194", " "},
    {57, "C8", "(1,166)", 
     "\!\(\*SubscriptBox[\(SiO\), \(2\)]\) (Quartz)", "hP9", "181", 
     " "},
    {58, "C9", "(1,169)", 
     "\!\(\*SubscriptBox[\(SiO\), \(2\)]\) (\[Beta]-Cristoballite)", 
     "cP12", "227", " "},
    {59, "C10", "(1,171)", 
     "\!\(\*SubscriptBox[\(SiO\), \(2\)]\) (\[Beta]-Tridymite)", 
     "hP12", "194", " "},
    {60, "\!\(\*SubscriptBox[\(C11\), \(a\)]\)", "(1,740)", 
     "\!\(\*SubscriptBox[\(CaC\), \(2\)]\)", "tI6", "139", " "},
    {61, "\!\(\*SubscriptBox[\(C11\), \(b\)]\)", "(1,740)", 
     "\!\(\*SubscriptBox[\(MoSi\), \(2\)]\)", "tI6", "139", " "},
    {62, "C12", "(1,175)", "\!\(\*SubscriptBox[\(CaSi\), \(2\)]\)", 
     "hR6", "166", " "},
    {63, "C14", "(1,180)", 
     "\!\(\*SubscriptBox[\(MgZn\), \(2\)]\) (Laves)", "hP12", "194", 
     " "},
    {64, "C15", "(1,490)", 
     "\!\(\*SubscriptBox[\(MgCu\), \(2\)]\) (Laves)", "cF24", "227", 
     " "},
    {65, "\!\(\*SubscriptBox[\(C15\), \(b\)]\)", "(1,490)", 
     "\!\(\*SubscriptBox[\(AuBe\), \(5\)]\)", "cF24", "216,196", 
     " "},
    {66, "C16", "(1,491)", "\!\(\*SubscriptBox[\(Al\), \(2\)]\)Cu", 
     "tI12", "140", " "},
    {67, "C18", "(1,495)", 
     "\!\(\*SubscriptBox[\(FeS\), \(2\)]\) (Marcasite)", "oP6", "58", 
     " "},
    {68, "C19", "(1,742)", "\!\(\*SubscriptBox[\(CdCl\), \(2\)]\)", 
     "hR3", "166", " "},
    {69, "C21", "(2,14)", 
     "\!\(\*SubscriptBox[\(TiO\), \(2\)]\) (Brookite)", "oP24", "61", 
     " "},
    {70, "C22", "(2,15)", "\!\(\*SubscriptBox[\(Fe\), \(2\)]\)P", 
     "hP9", "150", " "},
    {71, "C23", "(2,16)", "\!\(\*SubscriptBox[\(PBCl\), \(2\)]\)", 
     "oP12", "62", " "},
    {72, "C32", " ", "\!\(\*SubscriptBox[\(AlB\), \(2\)]\)", " hP3", 
     "191", " "},
    {73, "C33", " ", 
     "\!\(\*SubscriptBox[\(Bi\), \(2\)]\)\!\(\*SubscriptBox[\(Te\), \
\(2\)]\)S", "hR5" , "166", " "},
    {74, "C34", " ", "\!\(\*SubscriptBox[\(AuTe\), \(2\)]\)", "mC6", 
     "12", " "},
    {75, "C36", " ", "\!\(\*SubscriptBox[\(MgNi\), \(2\)]\)", "hP24", 
     "194", " "},
    {76, "C38", "(3,33)", "\!\(\*SubscriptBox[\(Cu\), \(2\)]\)Sb", 
     "tP6", "129", " "},
              {77, "C40", "(3,35)", 
     "\!\(\*SubscriptBox[\(CrSi\), \(5\)]\)", "hP9", "180", " "},
    {78, "C44", " ", "\!\(\*SubscriptBox[\(GeS\), \(2\)]\)", "oF72", 
     "43", " "},
    {79, "C46", " ", "\!\(\*SubscriptBox[\(AuTe\), \(2\)]\)", "oP24", 
     "28", " "},
    {80, "C49", " ", "\!\(\*SubscriptBox[\(ZrSi\), \(2\)]\)", "oC12", 
     "63", " "},
             {81, "C54", " ", "\!\(\*SubscriptBox[\(TiSi\), \(2\)]\)",
      "oF24", "70", " "},
             {82, "\!\(\*SubscriptBox[\(C\), \(c\)]\)", " ", 
     "\!\(\*SubscriptBox[\(ThSi\), \(2\)]\)", "tI12", "141", " "},
    {83, "\!\(\*SubscriptBox[\(C\), \(e\)]\)", " ", 
     "\!\(\*SubscriptBox[\(CoGe\), \(2\)]\)", "oA24", "41", " "},
    {84, "\!\(\*SubscriptBox[\(D0\), \(2\)]\)", " ", 
     "\!\(\*SubscriptBox[\(CoAs\), \(3\)]\) (Skutterudite)", "cI32", 
     "204", " "},
    {85, "\!\(\*SubscriptBox[\(D0\), \(3\)]\)", "(2,22)", 
     "\!\(\*SubscriptBox[\(BiF\), \(3\)]\)", "cF16", "225", " "},
    {86, "\!\(\*SubscriptBox[\(D0\), \(9\)]\)", "(2,31)", 
     "\!\(\*SubscriptBox[\(ReO\), \(3\)]\)", "cP4", "221", " "},
    {87, "\!\(\*SubscriptBox[\(D0\), \(11\)]\)", "(2,33)", 
     "\!\(\*SubscriptBox[\(Fe\), \(3\)]\)C", "oP16", "62", " "},
    {88, "\!\(\*SubscriptBox[\(D0\), \(18\)]\)", "(5,6)", 
     "\!\(\*SubscriptBox[\(Na\), \(3\)]\)As", "hP8", "194", " "},
    {89, "\!\(\*SubscriptBox[\(D0\), \(19\)]\)", "(5,7)", 
     "\!\(\*SubscriptBox[\(Ni\), \(3\)]\)Sn", "hP8", "194", " "},
    {90, "\!\(\*SubscriptBox[\(D0\), \(20\)]\)", " ", 
     "\!\(\*SubscriptBox[\(NiAl\), \(3\)]\)", "oP16" , "62", " "},
    {91, "\!\(\*SubscriptBox[\(D0\), \(21\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Cu\), \(3\)]\)P", "hP24", "165", " "},
    {92, "\!\(\*SubscriptBox[\(D0\), \(22\)]\)", "(7,13)", 
     "\!\(\*SubscriptBox[\(TiAl\), \(3\)]\)", "tI8", "139", " "},
    {93, "\!\(\*SubscriptBox[\(D0\), \(23\)]\)", "(7,14)", 
     "\!\(\*SubscriptBox[\(ZrAl\), \(3\)]\)", "tI16", "139", " "},
    {94, "\!\(\*SubscriptBox[\(D0\), \(24\)]\)", " ", 
     "\!\(\*SubscriptBox[\(TiNi\), \(3\)]\)", "hP16", "194", " "},
    {95, "\!\(\*SubscriptBox[\(D0\), \(c\)]\)", " ", 
     "\!\(\*SubscriptBox[\(U\), \(3\)]\)Si", "tI16", "140", " "},
    {96, "\!\(\*SubscriptBox[\(D0\), \(e\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Ni\), \(3\)]\)P", "tI32", "82", " "},
    {97, "\!\(\*SubscriptBox[\(D1\), \(3\)]\)", "(3,45)", 
     "\!\(\*SubscriptBox[\(BaAl\), \(4\)]\)", "tI10", "139", " "},
    {98, "\!\(\*SubscriptBox[\(D1\), \(a\)]\)", " ", 
     "\!\(\*SubscriptBox[\(MoNi\), \(4\)]\)", "tI10", "87", " "},
    {99, "\!\(\*SubscriptBox[\(D1\), \(b\)]\)", " ", 
     "\!\(\*SubscriptBox[\(UAl\), \(4\)]\)", "oI20", "74", " "},
    {100, "\!\(\*SubscriptBox[\(D1\), \(c\)]\)", " ", 
     "\!\(\*SubscriptBox[\(PtSn\), \(4\)]\)", "oC20", "41", " "},
    {101, "\!\(\*SubscriptBox[\(D1\), \(e\)]\)", " ", 
     "\!\(\*SubscriptBox[\(ThB\), \(4\)]\)", "tP20", "127", " "},
    {102, "\!\(\*SubscriptBox[\(D1\), \(f\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Mn\), \(4\)]\)B", "oF40", "70", " "},
    {103, "\!\(\*SubscriptBox[\(D2\), \(1\)]\)", " ", 
     "\!\(\*SubscriptBox[\(CaB\), \(6\)]\)", "cP7", "221", " "},
    {104, "\!\(\*SubscriptBox[\(D2\), \(3\)]\)", "(6,8)", 
     "\!\(\*SubscriptBox[\(NaZn\), \(13\)]\)", "cF112", "226", " "},
    {105, "\!\(\*SubscriptBox[\(D2\), \(b\)]\)", " ", 
     "\!\(\*SubscriptBox[\(ThMn\), \(12\)]\)", "tI26", "140", " "},
    {106, "\!\(\*SubscriptBox[\(D2\), \(c\)]\)", " ", 
     "\!\(\*SubscriptBox[\(U\), \(6\)]\)Mn", "tI28 ", "140", " "},
    {107, "\!\(\*SubscriptBox[\(D2\), \(d\)]\)", " ", 
     "\!\(\*SubscriptBox[\(CaCu\), \(5\)]\)", "hC6", "191", " "},
    {108, "\!\(\*SubscriptBox[\(D2\), \(f\)]\)", " ", 
     "\!\(\*SubscriptBox[\(UB\), \(12\)]\)", "cF52", "225", " "},
    {109, "\!\(\*SubscriptBox[\(D2\), \(h\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Al\), \(6\)]\)Mn", "oC28", "63", " "},
    {110, "\!\(\*SubscriptBox[\(D5\), \(1\)]\)", "(1,240)", 
     "\[Alpha]-\!\(\*SubscriptBox[\(Al\), \(2\)]\)\!\(\*SubscriptBox[\
\(O\), \(3\)]\) (Corundum)", "hR10", "167", " "},
    {111, "\!\(\*SubscriptBox[\(D5\), \(2\)]\)", "(1,744)", 
     "\!\(\*SubscriptBox[\(La\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", "hP5", "164", " "},
    {112, "\!\(\*SubscriptBox[\(D5\), \(3\)]\)", "(2,38)", 
     "\!\(\*SubscriptBox[\(Mn\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", "cI80", "206", " "},
    {113, "\!\(\*SubscriptBox[\(D5\), \(8\)]\)", "(3,49)", 
     "\!\(\*SubscriptBox[\(Sb\), \(2\)]\)\!\(\*SubscriptBox[\(S\), \
\(3\)]\)", "oP20", "62", " "},
    {114, "\!\(\*SubscriptBox[\(D5\), \(9\)]\)", "(3,51)", 
     "\!\(\*SubscriptBox[\(Zn\), \(3\)]\)\!\(\*SubscriptBox[\(P\), \
\(2\)]\)", "tP40", "137", " "},
              {115, "\!\(\*SubscriptBox[\(D5\), \(10\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Cr\), \(3\)]\)\!\(\*SubscriptBox[\(C\), \
\(2\)]\)", "oP20", "62", " "},
    {115, "\!\(\*SubscriptBox[\(D5\), \(13\)]\)", "(5,10)", 
     "\!\(\*SubscriptBox[\(Ni\), \(2\)]\)\!\(\*SubscriptBox[\(Al\), \
\(3\)]\)", "hC5", "164", " "},
    {117, "\!\(\*SubscriptBox[\(D5\), \(a\)]\)", " ", 
     "\!\(\*SubscriptBox[\(U\), \(3\)]\)\!\(\*SubscriptBox[\(Si\), \
\(2\)]\)", "tP10", "127", " "},
    {118, "\!\(\*SubscriptBox[\(D5\), \(c\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Pu\), \(2\)]\)\!\(\*SubscriptBox[\(C\), \
\(3\)]\)", "cI40", "220", " "},
    {119, "\!\(\*SubscriptBox[\(D7\), \(1\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Al\), \(4\)]\)\!\(\*SubscriptBox[\(C\), \
\(3\)]\)", "hR7", "166", " "},
    {120, "\!\(\*SubscriptBox[\(D7\), \(2\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Co\), \(3\)]\)\!\(\*SubscriptBox[\(S\), \
\(4\)]\) (Linaeite)", "cF56", "227", " "},
    {121, "\!\(\*SubscriptBox[\(D7\), \(3\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Th\), \(3\)]\)\!\(\*SubscriptBox[\(P\), \
\(4\)]\)", "cI26", "220", " "},
    {122, "\!\(\*SubscriptBox[\(D7\), \(b\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Ta\), \(3\)]\)\!\(\*SubscriptBox[\(B\), \
\(4\)]\)", "oI14", "71", " "},
    {123, "\!\(\*SubscriptBox[\(D8\), \(1\)]\)", "(1,497)", 
     "\!\(\*SubscriptBox[\(Fe\), \(3\)]\)\!\(\*SubscriptBox[\(Zn\), \
\(10\)]\)", " ", "229", " "},
    {124, "\!\(\*SubscriptBox[\(D8\), \(2\)]\)", "(1,497)", 
     "\!\(\*SubscriptBox[\(Cu\), \(5\)]\)\!\(\*SubscriptBox[\(Zn\), \
\(8\)]\) (\[Gamma]-brass)", "cI52", "217", " "},
    {125, "\!\(\*SubscriptBox[\(D8\), \(3\)]\)", "(1,497)", 
     "\!\(\*SubscriptBox[\(Cu\), \(9\)]\)\!\(\*SubscriptBox[\(Al\), \
\(4\)]\)", "cP52", "215", " "},
    {126, "\!\(\*SubscriptBox[\(D8\), \(4\)]\)", "(3,59)", 
     "\!\(\*SubscriptBox[\(Cr\), \(23\)]\)\!\(\*SubscriptBox[\(C\), \
\(6\)]\)", "cF116", "225", " "},
    {127, "\!\(\*SubscriptBox[\(D8\), \(5\)]\)", "(3,61;6,177)", 
     "\!\(\*SubscriptBox[\(Fe\), \(7\)]\)\!\(\*SubscriptBox[\(W\), \
\(6\)]\)", "hR13", "166", " "},
    {128, "\!\(\*SubscriptBox[\(D8\), \(6\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Cu\), \(15\)]\)\!\(\*SubscriptBox[\(Si\), \
\(4\)]\)", "cI76", "220", " "},
    {129, "\!\(\*SubscriptBox[\(D8\), \(8\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Mn\), \(5\)]\)\!\(\*SubscriptBox[\(Si\), \
\(3\)]\)", "hP16", "193", " "},
    {130, "\!\(\*SubscriptBox[\(D8\), \(9\)]\)", "(4,26)", 
     "\!\(\*SubscriptBox[\(Co\), \(9\)]\)\!\(\*SubscriptBox[\(S\), \
\(8\)]\)", "cF68", "225", " "},
    {131, "\!\(\*SubscriptBox[\(D8\), \(10\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Cr\), \(5\)]\)\!\(\*SubscriptBox[\(Al\), \
\(8\)]\)", "hR26", "160", " "},
    {132, "\!\(\*SubscriptBox[\(D8\), \(11\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Co\), \(2\)]\)\!\(\*SubscriptBox[\(Al\), \
\(5\)]\)", "hP28", "194", " "},
    {133, "\!\(\*SubscriptBox[\(D8\), \(a\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Th\), \(6\)]\)\!\(\*SubscriptBox[\(Mn\), \
\(23\)]\)", "cF116", "225", " "},
    {134, "\!\(\*SubscriptBox[\(D8\), \(b\)]\)", " ", 
     "FeCr (\[Sigma]-phase)", "tP30", "136", " "},
    {135, "\!\(\*SubscriptBox[\(D8\), \(f\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Ir\), \(3\)]\)\!\(\*SubscriptBox[\(Sn\), \
\(7\)]\)", "cI40", "229", " "},
    {136, "\!\(\*SubscriptBox[\(D8\), \(i\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Mo\), \(2\)]\)\!\(\*SubscriptBox[\(B\), \
\(5\)]\)", "hR7", "166", " "},
    {137, "\!\(\*SubscriptBox[\(D8\), \(h\)]\)", " ", 
     "\!\(\*SubscriptBox[\(W\), \(2\)]\)\!\(\*SubscriptBox[\(B\), \(5\
\)]\)", "hP14", "194", " "},
    {138, "\!\(\*SubscriptBox[\(D8\), \(l\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Cr\), \(5\)]\)\!\(\*SubscriptBox[\(B\), \
\(3\)]\)", "tI32", "140", " "},
    {139, "\!\(\*SubscriptBox[\(D8\), \(m\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Si\), \(3\)]\)\!\(\*SubscriptBox[\(W\), \
\(5\)]\)", "tI32", "140", " "},
    {140, "\!\(\*SubscriptBox[\(D10\), \(1\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Cr\), \(7\)]\)\!\(\*SubscriptBox[\(C\), \
\(3\)]\)", "oP40", "159", " "},
    {141, "\!\(\*SubscriptBox[\(D10\), \(2\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Fe\), \(3\)]\)\!\(\*SubscriptBox[\(Th\), \
\(7\)]\)", "hP20", "186", " "},
    {142, "\!\(\*SubscriptBox[\(E0\), \(1\)]\)", "(2,45)", "PbClF", 
     "tP6", "129", " "},
    {143, "\!\(\*SubscriptBox[\(E1\), \(1\)]\)", "(2,48)", 
     "\!\(\*SubscriptBox[\(CuFeS\), \(2\)]\) (Chalcopyrite)", " tI16",
      "122", " "},
    {144, "\!\(\*SubscriptBox[\(E2\), \(1\)]\)", "(1,300;3,68)", 
     "\!\(\*SubscriptBox[\(CaTiO\), \(3\)]\)(Perovskite)", "cP5", 
     "221", "*"},
    {145, "E3", " ", 
     "\!\(\*SubscriptBox[\(Al\), \(2\)]\)\!\(\*SubscriptBox[\(CdS\), \
\(4\)]\)", "tI14", "82", " "},
    {146, "\!\(\*SubscriptBox[\(E9\), \(3\)]\)", "(3,71)", 
     "\!\(\*SubscriptBox[\(Fe\), \(3\)]\)\!\(\*SubscriptBox[\(W\), \
\(3\)]\)C", "cF112", "227", " "},
    {147, "\!\(\*SubscriptBox[\(E9\), \(a\)]\)", " " , 
     "\!\(\*SubscriptBox[\(Al\), \(7\)]\)\!\(\*SubscriptBox[\(Cu\), \
\(2\)]\)Fe", "tP40", "128" , " "},
    {148, "\!\(\*SubscriptBox[\(E9\), \(b\)]\)", " ", 
     "\!\(\*SubscriptBox[\(AlLi\), \(3\)]\)\!\(\*SubscriptBox[\(N\), \
\(2\)]\)", "hP18", "189", " "},
    {149, "\!\(\*SubscriptBox[\(F0\), \(1\)]\)", " ", 
     "CoAsS (Cobaltite)", "cP12", "198", " "},
    {150, "\!\(\*SubscriptBox[\(F5\), \(1\)]\)", "(1,271)", 
     "\!\(\*SubscriptBox[\(CrNaS\), \(2\)]\)", "hR4", "160", " "},
    {151, "\!\(\*SubscriptBox[\(H1\), \(1\)]\)", "(1,350;4,170)", 
     "\!\(\*SubscriptBox[\(MgAl\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(\(4\)\(\\\ \)\)]\) (Spinel)", "cF56", "227", " "},
    {152, "\!\(\*SubscriptBox[\(H2\), \(4\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Cu\), \(3\)]\)\!\(\*SubscriptBox[\(S\), \
\(4\)]\)V (Sulvanite)", "cP8" , "215", " "},
    {153, "\!\(\*SubscriptBox[\(L1\), \(0\)]\)", "(1,484)", "CuAu", 
     "tP4", "123", " "},
    {154, "\!\(\*SubscriptBox[\(L1\), \(2\)]\)", "(1,486)", 
     "\!\(\*SubscriptBox[\(Cu\), \(3\)]\)Au", "cP4", "221", " "},
    {155, "\!\(\*SubscriptBox[\(L2\), \(1\)]\)", "(1,488)", 
     "\!\(\*SubscriptBox[\(AlCu\), \(2\)]\)Mn (Heusler)", "cF16", 
     "225", " "},
    {156, "L'\!\(\*SubscriptBox[\(2\), \(b\)]\)", " ", 
     "\!\(\*SubscriptBox[\(ThH\), \(2\)]\)", "tI16", "139", " "},
    {157, "L\!\(\*SubscriptBox[\('\), \(3\)]\)", " ", 
     "\!\(\*SubscriptBox[\(Fe\), \(2\)]\)N", "hP3", "194", " "},
    {158, "\!\(\*SubscriptBox[\(L6\), \(0\)]\)", " ", 
     "\!\(\*SubscriptBox[\(CuTi\), \(3\)]\)", "tP4", "123", " "}
    };
  structures = {
    (*--- A1 ---*)
    {"Cu", "\[Alpha]-Sr", "Al", "Sc", "Yb", 
     "\[Alpha]-Th", "Ni", "Rh", "Pd", "Ir", "Pt", "Ag", "Au", "Pb", 
     "Ne", "Ar", "Kr", "Xe", "\[Beta]-La", "Ce", "\[Beta]-Pr", 
     "\[Gamma]-Mn", "\[Gamma]-Fe", "Co"},
    (*--- A2 ---*)
    {"W", "Li", "Na", "K", "Rb", "Cs", "V", "Nb", "Ta", 
     "Mo", "Cr", "\[Beta]-Ti", "\[Beta]-Zr", "\[Gamma]-U"},
    (*--- A3 ---*)
    {"Mg", "Be", "\[Beta]-Sr", "Y", "Gd", "Tb", "Dy", 
     "Ho", "Er", "Tm", "Tc", "Re", "Ru", "Os", "Zn", "Cd", "Ce", 
     "\[Alpha]-Pr", "\[Alpha]-Ti", "\[Alpha]-Zr"},
    (*--- A4 ---*)
    {"C", "Si", "\[Alpha]-Sn", "Ge"},
    (*--- A5 ---*)
    {"\[Beta]-Zn"},
    (*--- A6 ---*)
    {"In", "\[Gamma]-Mn", 
     "\!\(\*SubscriptBox[\(InPd\), \(3\)]\)"},
    (*--- A7 ---*)
    {"As", "Sb", "Bi"},
    (*--- A8 ---*)
    {"Se", "Te"},
    (*--- A9 ---*)
    {"C"},
    (*--- A10 ---*)
    {"Hg"},
    (*--- A11 ---*)
    {"Ga"},
    (*--- A12 ---*)
    {"\[Alpha]-Mn", 
     "\!\(\*SubscriptBox[\(MoRe\), \(4\)]\)"},
    (*--- A13 ---*)
    {"\[Beta]-Mn", 
     "\[Gamma]-\!\(\*SubscriptBox[\(Cu\), \(5\)]\)Si", 
     "\[Beta]'-\!\(\*SubscriptBox[\(Au\), \(3\)]\)Al"},
    (*--- A14 ---*)
    {"J", "Br"},
    (*--- A15 ---*)
    {"\!\(\*SubscriptBox[\(V\), \(3\)]\)Si", 
     "\!\(\*SubscriptBox[\(Cr\), \(3\)]\)Si", 
     "\!\(\*SubscriptBox[\(Mo\), \(3\)]\)Si", 
     "\!\(\*SubscriptBox[\(Ge\), \(3\)]\)V", 
     "\!\(\*SubscriptBox[\(GeCr\), \(3\)]\)"},
    (*--- A16 ---*)
    {"S"},
    (*--- A17 ---*)
    {"P"},
    (*--- A18 ---*)
    {"Cl"},
    (*--- A19 ---*)
    {"Po"},
    (*--- A20 ---*)
    {"\[Alpha]-U"},
    (*--- B1 ---*)
    {"NaCl", "CaO", "BaO", "EuSe", "TiO", "TiC", "VO", 
     "CoO", "MnO", "FeO"},
    (*--- B2 ---*)
    {"CsCl", "CsBr", "RuSi", "AlFe", "AuMg", "CoFe", 
     "TiZn"},
    (*--- B3 ---*)
    {"ZnS", "AlAs", "ZnSe", "ZnTe", "GaP", "GaAs", "InP",
      "InSb"},
    (*--- B4 ---*)
    {"ZnO", "BeO", "MgTe", "AlN", "\[Gamma]-MnS", "ZnS", 
     "CdSe", "GaN", "InN"},
    (*--- B8_1 ---*)
    {"NiAs", "TiSe", "TiSb", "CrTe", "PdTe", "CoSb", 
     "NiSb", "AsMn", "BiNi"},
    (*--- B8_2 ---*)
    {"\!\(\*SubscriptBox[\(InNi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(SnMn\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(GeFe\), \(2\)]\)", "BeSiZr", "CoFeGe", 
     "FeGeNi"},
    (*--- B9 ---*)
    {"HgS"},
    (*--- B10 ---*)
    {"PbO", "LiOH", 
     "\!\(\*SubscriptBox[\(NH\), \(4\)]\)SH"},
    (*--- B11 ---*) 
    {"PbO", "\[Beta]-AgZr", "SZr"},
    (*--- B12 ---*)
    {"BN"},
    (*--- B13 ---*)
    {"NiS", "NiSe"},
    (*--- B16 ---*)
    {"GeS"},
    (*--- B17 ---*)
    {"PtS"},
    (*--- B18 ---*)
    {"CuS", "CuSe"},
    (*--- B19 ---*)
    {"\[Beta]'-AuCd", "CdMg"},
    (*--- B20 ---*)
    {"FeSi", "CrSi", "MnSi", "CoSi", "GaPd", "GaPt", 
     "GeCr", "SnRh"},
    (*--- B27 ---*)
    {"FeB", "USi", "MnB", "CoB", "FeB"},
    (*--- B29 ---*)
    {"\!\(\*SubscriptBox[\(PbS\), \(2\)]\)Sn", "SnS", 
     "SeSn"},
    (*--- B31 ---*)
    {"MnP", "CrP", "CrAs", "MnP", "MnAs", "FeP", "FeAs",
      "CoP", "CoAs", "SnPd"},
    (*--- B32 ---*)
    {"NaTl", "AlLi", "ZnLi", "CdLi", "GaLi", "InLi", 
     "InNa"},
    (*--- B34 ---*)
    {"PdS"},
    (*--- B35 ---*)
    {"CoSn", "InNi", "TlPt", "\[Beta]-SnFe"},
    (*--- B37 ---*)
    {"FeB", "InTe", "STl", "SeTl"},
    (*--- B_e ---*)
    {"CdSb", "SbZn"},
    (*--- B_f ---*)
    {"\[Zeta]-CrB", "BNb", "BTa", "BNi"},
    (*--- B_g ---*)
    {"\[Delta]-MoB", "\[Delta]-BW"},
    (*--- B_h ---*)
    {"WC", "MoP"},
    (*--- B_i ---*)
    {"\[Gamma]'-MoC", "\[Alpha]-AsTi", "PTi"},
    (*--- C1 ---*)
    {"\!\(\*SubscriptBox[\(CaF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(EuF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(UO\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(PuO\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CuF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CdF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(PtAl\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Ga\), \(2\)]\)Au", 
     "\!\(\*SubscriptBox[\(Ga\), \(2\)]\)Pt", 
     "\!\(\*SubscriptBox[\(In\), \(2\)]\)Au", 
     "\!\(\*SubscriptBox[\(Li\), \(2\)]\)O"},
    (*--- C1_b ---*)
    {"AgAsMg", "ZnNLi", "CuMgSb", "CdSbCu"},
    (*--- C2 ---*)
    {"\!\(\*SubscriptBox[\(FeS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(MnS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(MnSe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(MnTe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CoS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CoSe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(NiS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(NiSe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(RuS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(RuSe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(RuTe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(PdAs\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(PdSb\), \(2\)]\)"},
    (*--- C3 ---*)
    {"\!\(\*SubscriptBox[\(Cu\), \(2\)]\)O", 
     "\!\(\*SubscriptBox[\(Ag\), \(2\)]\)O"},
    (*--- C4 ---*)
    {"\!\(\*SubscriptBox[\(TiO\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(MgF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(VO\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CrO\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(MnF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(FeF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CoF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(NiF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(RuO\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(PdF\), \(2\)]\)" \
"\!\(\*SubscriptBox[\(OsO\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(IrO\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(ZnF\), \(2\)]\)"},
    (*--- C5 ---*)
    {"\!\(\*SubscriptBox[\(TiO\), \(2\)]\)"},
    (*--- C6 ---*)
    {"\!\(\*SubscriptBox[\(CdJ\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(MgJ\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(MgBr\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(TiS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(TiSe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(TiTe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(NiTe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(PTSe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(S\), \(2\)]\)Zr"},
    (*--- C7 ---*)
    {"\!\(\*SubscriptBox[\(MoS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(WS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(WSe\), \(2\)]\)"},
    (*--- C8 ---*)
    {"\!\(\*SubscriptBox[\(BeF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(GeO\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(SiO\), \(2\)]\)"},
    (*--- C9 ---*)
    {"\!\(\*SubscriptBox[\(BeF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(SiO\), \(2\)]\)"},
    (*--- C10 ---*)
    {"\!\(\*SubscriptBox[\(SiO\), \(2\)]\)"},
    (*--- C11_a ---*)
    {"\!\(\*SubscriptBox[\(CaC\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(SrC\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(BaC\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(LaC\), \(2\)]\)"},
    (*--- C11_b ---*)
    {"\!\(\*SubscriptBox[\(MoSi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(WSi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(ReSi\), \(2\)]\)"},
    (*--- C12 ---*)
    {"\!\(\*SubscriptBox[\(CaSi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Ge\), \(2\)]\)Ca"},
    (*--- C14 ---*)
    {"\!\(\*SubscriptBox[\(Cu\), \
\(3\)]\)\!\(\*SubscriptBox[\(Mg\), \(2\)]\)Si", 
     "\!\(\*SubscriptBox[\(MgZn\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(KNa\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CaLi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CaMg\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(SrMg\), \(2\)]\)"},
    (*--- C15 ---*)
    {"\!\(\*SubscriptBox[\(Cu\), \(2\)]\)Mg", 
     "\!\(\*SubscriptBox[\(Al\), \(2\)]\)Ca", 
     "\!\(\*SubscriptBox[\(Co\), \(2\)]\)Zr", 
     "\!\(\*SubscriptBox[\(Co\), \(2\)]\)Nb", 
     "\!\(\*SubscriptBox[\(Fe\), \(2\)]\)Zr", 
     "\!\(\*SubscriptBox[\(Fe\), \(2\)]\)Gd"},
    (*--- C15_b ---*)
    {"\!\(\*SubscriptBox[\(AuBe\), \(5\)]\)", 
     "\!\(\*SubscriptBox[\(FeBe\), \(5\)]\)", 
     "\!\(\*SubscriptBox[\(PdBe\), \(5\)]\)"},
    (*--- C16 ---*)
    {"\!\(\*SubscriptBox[\(CuAl\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Co\), \(2\)]\)B", 
     "\!\(\*SubscriptBox[\(Mn\), \(2\)]\)B", 
     "\!\(\*SubscriptBox[\(Ni\), \(2\)]\)B", 
     "\!\(\*SubscriptBox[\(Ge\), \(2\)]\)Fe", 
     "\!\(\*SubscriptBox[\(Sn\), \(2\)]\)Mn", 
     "\!\(\*SubscriptBox[\(Sn\), \(2\)]\)Fe", 
     "\!\(\*SubscriptBox[\(Sn\), \(2\)]\)Co"},
    (*--- C18 ---*)
    {"\!\(\*SubscriptBox[\(FeS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CrSb\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(FeSe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(FeTe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(FeP\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CosB\), \(2\)]\)"},
    (*--- C19 ---*)
    {"\!\(\*SubscriptBox[\(CdCl\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CoCl\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(NiCl\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(NiJ\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(RuCl\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(ZnBr\), \(2\)]\)"},
    (*--- C21 ---*)
    {"\!\(\*SubscriptBox[\(TiO\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(TeO\), \(2\)]\)"},
    (*--- C22 ---*)
    {"\!\(\*SubscriptBox[\(Mn\), \(2\)]\)P", 
     "\!\(\*SubscriptBox[\(Fe\), \(2\)]\)P", 
     "\!\(\*SubscriptBox[\(Ni\), \(2\)]\)P", 
     "\!\(\*SubscriptBox[\(Pd\), \(2\)]\)Si", 
     "\!\(\*SubscriptBox[\(GePd\), \(2\)]\)"},
    (*--- C23 ---*)
    {"\!\(\*SubscriptBox[\(PbCl\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(BaCl\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(BaBr\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(EuCl\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Co\), \(2\)]\)P", 
     "\!\(\*SubscriptBox[\(PbBr\), \(2\)]\)"},
    (*--- C32 ---*)
    {"\!\(\*SubscriptBox[\(AlB\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Al\), \(2\)]\)Th", 
     "\!\(\*SubscriptBox[\(CeGa\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Ni\), \(2\)]\)Th", 
     "\!\(\*SubscriptBox[\(Be\), \(2\)]\)Zr"},
    (*--- C33 ---*)
    {"\!\(\*SubscriptBox[\(Bi\), \
\(2\)]\)\!\(\*SubscriptBox[\(STe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Bi\), \(2\)]\)\!\(\*SubscriptBox[\(Se\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Bi\), \(2\)]\)\!\(\*SubscriptBox[\(Te\), \
\(3\)]\)" "\!\(\*SubscriptBox[\(Bi\), \(2\)]\)\!\(\*SubscriptBox[\(Se\
\), \(3\)]\)"},
    (*--- C34 ---*)
    {"\!\(\*SubscriptBox[\(AuTe\), \(2\)]\)"},
    (*--- C36 ---*)
    {"\!\(\*SubscriptBox[\(MgNi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Co\), \(2\)]\)Ti", 
     "\!\(\*SubscriptBox[\(Cr\), \(2\)]\)Hf", 
     "\!\(\*SubscriptBox[\(HfMn\), \(2\)]\)"},
    (*--- C38 ---*)
    {"\!\(\*SubscriptBox[\(Cu\), \(2\)]\)Sb", 
     "\!\(\*SubscriptBox[\(Mn\), \(2\)]\)As", 
     "\!\(\*SubscriptBox[\(Mn\), \(2\)]\)Sb", 
     "\!\(\*SubscriptBox[\(Cr\), \(2\)]\)As"},
    (*--- C40 ---*)
    {"\!\(\*SubscriptBox[\(Cr\), \(2\)]\)Si", 
     "\!\(\*SubscriptBox[\(VSi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(NbSi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(TaSi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Ge\), \(2\)]\)Nb"},
    (*--- C44 ---*)
    {"\!\(\*SubscriptBox[\(GeS\), \(2\)]\)"},
    (*--- C46 ---*)
    {"\!\(\*SubscriptBox[\(AuTe\), \(2\)]\)"},
    (*--- C49 ---*)
    {"\!\(\*SubscriptBox[\(ZrSi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Ge\), \(2\)]\)Zr", 
     "\!\(\*SubscriptBox[\(HfSi\), \(2\)]\)"},
    (*--- C54 ---*)
    {"\!\(\*SubscriptBox[\(TiSi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Ge\), \(2\)]\)Ti", 
     "\!\(\*SubscriptBox[\(Sn\), \(2\)]\)Zr"},
    (*--- C_c ---*)
    {"\!\(\*SubscriptBox[\(ThSi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CeSi\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Ge\), \(2\)]\)Pu", 
     "\!\(\*SubscriptBox[\(LaSi\), \(2\)]\)"},
    (*--- C_e ---*)
    {"\!\(\*SubscriptBox[\(CoGe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(PdSn\), \(2\)]\)"},
    (*--- D0_2 ---*)
    {"\!\(\*SubscriptBox[\(CoAs\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(CoSb\), \(3\)]\)"},
    (*--- D0_3 ---*)
    {"\[Alpha]-\!\(\*SubscriptBox[\(BiF\), \(3\)]\)", 
     "\[Alpha]-\!\(\*SubscriptBox[\(Fe\), \(3\)]\)Si", 
     "\[Beta]-\!\(\*SubscriptBox[\(Cu\), \(3\)]\)Sb", 
     "\!\(\*SubscriptBox[\(LaMg\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(CeMg\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(PrMg\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Fe\), \(3\)]\)Al" \
"\!\(\*SubscriptBox[\(BiLi\), \(3\)]\)"},
    (*--- D0_9 ---*)
    {"\!\(\*SubscriptBox[\(ReO\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(ScF\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(TaF\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(MoF\), \(3\)]\)"},
    (*--- D0_11 ---*)
    {"\!\(\*SubscriptBox[\(Fe\), \(3\)]\)C", 
     "\!\(\*SubscriptBox[\(Co\), \(3\)]\)C", 
     "(Fe,Mn\!\(\*SubscriptBox[\()\), \(3\)]\)C"},
    (*--- D0_18 ---*)
               {"\!\(\*SubscriptBox[\(Na\), \
\(3\)]\)As", "\!\(\*SubscriptBox[\(Li\), \(3\)]\)P", 
     "\!\(\*SubscriptBox[\(Li\), \(3\)]\)As", 
     "\!\(\*SubscriptBox[\(Na\), \(3\)]\)P", 
     "\!\(\*SubscriptBox[\(Na\), \(3\)]\)Sb", 
     "\!\(\*SubscriptBox[\(BiNa\), \(3\)]\)"},
    (*--- D0_19 ---*)
    {"\!\(\*SubscriptBox[\(Ni\), \(3\)]\)Sn", 
     "\!\(\*SubscriptBox[\(Co\), \(3\)]\)N", 
     "\!\(\*SubscriptBox[\(Cd\), \(3\)]\)Mg", 
     "\!\(\*SubscriptBox[\(CdMg\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Hg\), \(3\)]\)Li"},
    (*--- D0_20 ---*)
    {"\!\(\*SubscriptBox[\(Al\), \(3\)]\)Ni"},
    (*--- D0_21 ---*)
    {"\!\(\*SubscriptBox[\(Cu\), \(3\)]\)P", 
     "\!\(\*SubscriptBox[\(AsCu\), \(3\)]\)"},
    (*--- D0_22 ---*)
    {"\!\(\*SubscriptBox[\(Al\), \(3\)]\)Ti", 
     "\!\(\*SubscriptBox[\(Al\), \(3\)]\)V", 
     "\!\(\*SubscriptBox[\(Al\), \(3\)]\)Nb", 
     "\!\(\*SubscriptBox[\(Al\), \(3\)]\)Ta", 
     "\!\(\*SubscriptBox[\(Ga\), \(3\)]\)Ti"},
    (*--- D0_23 ---*)
    {"\!\(\*SubscriptBox[\(Al\), \(3\)]\)Zr"},
    (*--- D0_24 ---*)
    {"\!\(\*SubscriptBox[\(Ni\), \(3\)]\)Ti", 
     "\!\(\*SubscriptBox[\(Pd\), \(3\)]\)Ti", 
     "\!\(\*SubscriptBox[\(Pd\), \(3\)]\)U", 
     "\!\(\*SubscriptBox[\(Pt\), \(3\)]\)Zr"},
    (*--- D0_c ---*)
    {"\!\(\*SubscriptBox[\(SiU\), \(3\)]\)"},
    (*--- D0_e ---*)
    {"\!\(\*SubscriptBox[\(Ni\), \(3\)]\)P", 
     "\!\(\*SubscriptBox[\(Cr\), \(3\)]\)P", 
     "\!\(\*SubscriptBox[\(Mn\), \(3\)]\)P", 
     "\!\(\*SubscriptBox[\(Fe\), \(3\)]\)P", 
     "\!\(\*SubscriptBox[\(Mo\), \(3\)]\)P"},
    (*--- D1_3 ---*)
    {"\!\(\*SubscriptBox[\(Al\), \(4\)]\)Ba", 
     "\!\(\*SubscriptBox[\(Al\), \(4\)]\)Sr", 
     "\!\(\*SubscriptBox[\(Al\), \(4\)]\)Ca", 
     "\!\(\*SubscriptBox[\(LaAl\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(ceAl\), \(4\)]\)"},
    (*--- D1_a ---*)
    {"\!\(\*SubscriptBox[\(MoNi\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(Ni\), \(4\)]\)W"},
    (*--- D1_b ---*)
    {"\!\(\*SubscriptBox[\(Al\), \(4\)]\)U", 
     "\!\(\*SubscriptBox[\(Al\), \(4\)]\)Np", 
     "\!\(\*SubscriptBox[\(Al\), \(4\)]\)Pu"},
    (*--- D1_c ---*)
    {"\!\(\*SubscriptBox[\(PtSn\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(AuSn\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(PdSn\), \(4\)]\)"},
    (*--- D1_e ---*)
    {"\!\(\*SubscriptBox[\(B\), \(4\)]\)Th", 
     "\!\(\*SubscriptBox[\(B\), \(4\)]\)Ce", 
     "\!\(\*SubscriptBox[\(B\), \(4\)]\)U"},
    (*--- D1_f ---*)
    {"\!\(\*SubscriptBox[\(Mn\), \(4\)]\)B", 
     "\!\(\*SubscriptBox[\(Mn\), \(4\)]\)Cr"},
    (*--- D2_1 ---*)
    {"\!\(\*SubscriptBox[\(CaB\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(SrB\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(BaB\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(YB\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(LaB\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(CeB\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(GdB\), \(6\)]\)"},
    (*--- D2_3 ---*)
    {"\!\(\*SubscriptBox[\(NaZn\), \(13\)]\)", 
     "\!\(\*SubscriptBox[\(CeBe\), \(13\)]\)", 
     "\!\(\*SubscriptBox[\(ZrBe\), \(13\)]\)", 
     "\!\(\*SubscriptBox[\(Zn\), \(13\)]\)Na", 
     "\!\(\*SubscriptBox[\(Zn\), \(13\)]\)Ca", 
     "\!\(\*SubscriptBox[\(Zn\), \(13\)]\)Sr", 
     "\!\(\*SubscriptBox[\(Zn\), \(13\)]\)Ba", 
     "\!\(\*SubscriptBox[\(Cd\), \(13\)]\)Rb"},
    (*--- D2_b ---*)
    {"\!\(\*SubscriptBox[\(Mn\), \(12\)]\)Th", 
     "\!\(\*SubscriptBox[\(Be\), \(12\)]\)Cr", 
     "\!\(\*SubscriptBox[\(Be\), \(12\)]\)Nb", 
     "\!\(\*SubscriptBox[\(Be\), \(12\)]\)Mo", 
     "\!\(\*SubscriptBox[\(Be\), \(12\)]\)V"},
    (*--- D2_c ---*)
    {"\!\(\*SubscriptBox[\(MnU\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(FeU\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(CoU\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(NiU\), \(6\)]\)"},
    (*--- D2_d ---*)
    {"\!\(\*SubscriptBox[\(CaZn\), \(5\)]\)", 
     "\!\(\*SubscriptBox[\(Ag\), \(5\)]\)Ba", 
     "\!\(\*SubscriptBox[\(Ag\), \(5\)]\)Sr", 
     "\!\(\*SubscriptBox[\(Au\), \(5\)]\)Ba", 
     "\!\(\*SubscriptBox[\(CaCu\), \(5\)]\)", 
     "\!\(\*SubscriptBox[\(CaNi\), \(5\)]\)"},
    (*--- D2_f ---*)
    {"\!\(\*SubscriptBox[\(UB\), \(12\)]\)", 
     "\!\(\*SubscriptBox[\(ZrB\), \(12\)]\)"},
    (*--- D2_h ---*)
    {"\!\(\*SubscriptBox[\(Al\), \(6\)]\)Mn"},
    (*--- D5_1 ---*)
    {"\[Alpha]-\!\(\*SubscriptBox[\(Al\), \(2\)]\)\!\(\
\*SubscriptBox[\(O\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Ti\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(V\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \(3\
\)]\)", "\!\(\*SubscriptBox[\(Cr\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", 
     "\[Alpha]-\!\(\*SubscriptBox[\(Fe\), \(2\)]\)\!\(\*SubscriptBox[\
\(O\), \(\(3\)\(\\\ \)\)]\)", 
     "\!\(\*SubscriptBox[\(Rh\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", 
     "\[Alpha]-\!\(\*SubscriptBox[\(Ga\), \(2\)]\)\!\(\*SubscriptBox[\
\(O\), \(3\)]\)"},
    (*--- D5_2 ---*)
    {"\!\(\*SubscriptBox[\(La\), \
\(2\)]\)\!\(\*SubscriptBox[\(O\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Ce\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Pr\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Nd\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Th\), \(2\)]\)\!\(\*SubscriptBox[\(N\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Pu\), \(2\)]\)\!\(\*SubscriptBox[\(SO\), \
\(3\)]\)"},
    (*--- D5_3 ---*)
    {"\!\(\*SubscriptBox[\(Mn\), \
\(2\)]\)\!\(\*SubscriptBox[\(O\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Be\), \(3\)]\)\!\(\*SubscriptBox[\(N\), \
\(2\)]\)", 
     "\!\(\*SubscriptBox[\(Mg\), \(3\)]\)\!\(\*SubscriptBox[\(N\), \
\(2\)]\)", 
     "\!\(\*SubscriptBox[\(Mg\), \(3\)]\)\!\(\*SubscriptBox[\(P\), \
\(2\)]\)", 
     "\!\(\*SubscriptBox[\(Sc\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Nd\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Sm\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Ho\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(3\)]\)"},
    (*--- D5_8 ---*)
    {"\!\(\*SubscriptBox[\(Sb\), \
\(2\)]\)\!\(\*SubscriptBox[\(S\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Th\), \(2\)]\)\!\(\*SubscriptBox[\(S\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(U\), \(2\)]\)\!\(\*SubscriptBox[\(S\), \(3\
\)]\)", "\!\(\*SubscriptBox[\(Bi\), \(2\)]\)\!\(\*SubscriptBox[\(S\), \
\(3\)]\)"},
    (*--- D5_9 ---*)
    {"\!\(\*SubscriptBox[\(Zn\), \
\(3\)]\)\!\(\*SubscriptBox[\(P\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Zn\), \(3\)]\)\!\(\*SubscriptBox[\(As\), \
\(2\)]\)", 
     "\!\(\*SubscriptBox[\(Cd\), \(3\)]\)\!\(\*SubscriptBox[\(P\), \
\(2\)]\)", 
     "\!\(\*SubscriptBox[\(Cd\), \(3\)]\)\!\(\*SubscriptBox[\(As\), \
\(2\)]\)"},
    (*--- D5_10 ---*)
    {"\!\(\*SubscriptBox[\(C\), \
\(2\)]\)\!\(\*SubscriptBox[\(Cr\), \(3\)]\)"},
    (*--- D5_13 ---*)
    {"\!\(\*SubscriptBox[\(Ni\), \
\(2\)]\)\!\(\*SubscriptBox[\(Al\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Ga\), \(3\)]\)\!\(\*SubscriptBox[\(Pt\), \
\(2\)]\)", 
     "\!\(\*SubscriptBox[\(In\), \(3\)]\)\!\(\*SubscriptBox[\(Pd\), \
\(2\)]\)", 
     "\!\(\*SubscriptBox[\(In\), \(3\)]\)\!\(\*SubscriptBox[\(Pt\), \
\(2\)]\)", 
     "\!\(\*SubscriptBox[\(Sn\), \(3\)]\)\!\(\*SubscriptBox[\(Pt\), \
\(2\)]\)"},
    (*--- D5_a ---*)
    {"\!\(\*SubscriptBox[\(Si\), \
\(2\)]\)\!\(\*SubscriptBox[\(U\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Al\), \(2\)]\)\!\(\*SubscriptBox[\(Th\), \
\(3\)]\)"},
    (*--- D5_c ---*)
    {"\!\(\*SubscriptBox[\(Pu\), \
\(2\)]\)\!\(\*SubscriptBox[\(C\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(3\)]\)\!\(\*SubscriptBox[\(U\), \(2\
\)]\)"},
    (*--- D7_1 ---*)
    {"\!\(\*SubscriptBox[\(Al\), \
\(4\)]\)\!\(\*SubscriptBox[\(C\), \(3\)]\)"},
    (*--- D7_2 ---*)
    {"\!\(\*SubscriptBox[\(Co\), \
\(3\)]\)\!\(\*SubscriptBox[\(S\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(FeNi\), \(2\)]\)\!\(\*SubscriptBox[\(S\), \
\(4\)]\)"},
    (*--- D7_3 ---*)
    {"\!\(\*SubscriptBox[\(Th\), \
\(3\)]\)\!\(\*SubscriptBox[\(P\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(As\), \(4\)]\)\!\(\*SubscriptBox[\(Th\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Bi\), \(4\)]\)\!\(\*SubscriptBox[\(U\), \
\(3\)]\)"},
    (*--- D7_b ---*)
    {"\!\(\*SubscriptBox[\(Ta\), \
\(3\)]\)\!\(\*SubscriptBox[\(B\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(B\), \(4\)]\)\!\(\*SubscriptBox[\(Cr\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(B\), \(4\)]\)\!\(\*SubscriptBox[\(Mn\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(B\), \(4\)]\)\!\(\*SubscriptBox[\(FeMo\), \
\(2\)]\)"},
    (*--- D8_1 ---*)
    {"\!\(\*SubscriptBox[\(Fe\), \
\(3\)]\)\!\(\*SubscriptBox[\(Zn\), \(10\)]\)"},
    (*--- D8_2 ---*)
    {"\[Gamma]-\!\(\*SubscriptBox[\(Cu\), \(5\)]\)\!\(\
\*SubscriptBox[\(Zn\), \(8\)]\)", 
     "\!\(\*SubscriptBox[\(Al\), \(8\)]\)\!\(\*SubscriptBox[\(V\), \
\(5\)]\)"},
    (*--- D8_3 ---*)
    {"\[Gamma]-\!\(\*SubscriptBox[\(Cu\), \(9\)]\)\!\(\
\*SubscriptBox[\(Al\), \(4\)]\)"},
    (*--- D8_4 ---*)
    {"\!\(\*SubscriptBox[\(Cr\), \
\(23\)]\)\!\(\*SubscriptBox[\(C\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(W\), \(2\)]\)\!\(\*SubscriptBox[\(C\), \(6\
\)]\)\!\(\*SubscriptBox[\(Cr\), \(21\)]\)", 
     "\!\(\*SubscriptBox[\(Mn\), \(23\)]\)\!\(\*SubscriptBox[\(C\), \
\(6\)]\)", 
     "\!\(\*SubscriptBox[\(Fe\), \(21\)]\)\!\(\*SubscriptBox[\(Mo\), \
\(2\)]\)\!\(\*SubscriptBox[\(C\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(Fe\), \(21\)]\)\!\(\*SubscriptBox[\(W\), \
\(2\)]\)\!\(\*SubscriptBox[\(C\), \(6\)]\)"},
    (*--- D8_5 ---*)
    {"\!\(\*SubscriptBox[\(Fe\), \
\(7\)]\)\!\(\*SubscriptBox[\(W\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(Fe\), \(7\)]\)\!\(\*SubscriptBox[\(Mo\), \
\(6\)]\)", 
     "\!\(\*SubscriptBox[\(Fe\), \(7\)]\)\!\(\*SubscriptBox[\(W\), \
\(6\)]\)" "\!\(\*SubscriptBox[\(Co\), \(7\)]\)\!\(\*SubscriptBox[\(Mo\
\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(Co\), \(7\)]\)\!\(\*SubscriptBox[\(W\), \
\(6\)]\)"},
    (*--- D8_6 ---*)
    {"\!\(\*SubscriptBox[\(Cu\), \
\(15\)]\)\!\(\*SubscriptBox[\(Si\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(Na\), \(15\)]\)\!\(\*SubscriptBox[\(Pb\), \
\(4\)]\)"},
    (*--- D8_8 ---*)
    {"\!\(\*SubscriptBox[\(Mn\), \
\(5\)]\)\!\(\*SubscriptBox[\(Si\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Ge\), \(3\)]\)\!\(\*SubscriptBox[\(Mn\), \
\(5\)]\)", 
     "\!\(\*SubscriptBox[\(Ge\), \(3\)]\)\!\(\*SubscriptBox[\(Ti\), \
\(5\)]\)", 
     "\!\(\*SubscriptBox[\(Pb\), \(3\)]\)\!\(\*SubscriptBox[\(Zr\), \
\(5\)]\)"},
    (*--- D8_9 ---*)
    {"\!\(\*SubscriptBox[\(Co\), \
\(9\)]\)\!\(\*SubscriptBox[\(S\), \(8\)]\)", 
     "(Ni,Fe\!\(\*SubscriptBox[\()\), \
\(9\)]\)\!\(\*SubscriptBox[\(S\), \(8\)]\)"},
    (*--- D8_10 ---*)
    {"\!\(\*SubscriptBox[\(Al\), \
\(8\)]\)\!\(\*SubscriptBox[\(Cr\), \(5\)]\)"},
    (*--- D8_11 ---*)
    {"\!\(\*SubscriptBox[\(Al\), \
\(5\)]\)\!\(\*SubscriptBox[\(Co\), \(2\)]\)"},
    (*--- D8_a ---*)
    {"\!\(\*SubscriptBox[\(Mn\), \
\(23\)]\)\!\(\*SubscriptBox[\(Th\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(Cu\), \(16\)]\)\!\(\*SubscriptBox[\(Mg\), \
\(6\)]\)\!\(\*SubscriptBox[\(Si\), \(7\)]\)", 
     "\!\(\*SubscriptBox[\(Mn\), \(23\)]\)\!\(\*SubscriptBox[\(Th\), \
\(6\)]\)"},
    (*--- D8_b ---*)
    {"\[Sigma]-FeCr", "CoCr", "CoV", "FeMo", "FeV", 
     "MnTi", "NbRe"},
    (*--- D8_f ---*)
    {"\!\(\*SubscriptBox[\(Ge\), \
\(7\)]\)\!\(\*SubscriptBox[\(Ir\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Ga\), \(7\)]\)\!\(\*SubscriptBox[\(Pd\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Ga\), \(7\)]\)\!\(\*SubscriptBox[\(Pt\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(In\), \(7\)]\)\!\(\*SubscriptBox[\(Pt\), \
\(3\)]\)"},
    (*--- D8_i ---*)
    {"\!\(\*SubscriptBox[\(Mo\), \
\(2\)]\)\!\(\*SubscriptBox[\(B\), \(5\)]\)"},
    (*--- D8_h ---*)
    {"\!\(\*SubscriptBox[\(W\), \
\(2\)]\)\!\(\*SubscriptBox[\(B\), \(5\)]\)", 
     "\!\(\*SubscriptBox[\(Ti\), \(2\)]\)\!\(\*SubscriptBox[\(B\), \
\(5\)]\)"},
    (*--- D8_l ---*)
    {"\!\(\*SubscriptBox[\(Cr\), \
\(5\)]\)\!\(\*SubscriptBox[\(B\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Nb\), \(5\)]\)\!\(\*SubscriptBox[\(Si\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Si\), \(3\)]\)\!\(\*SubscriptBox[\(Ta\), \
\(5\)]\)"},
    (*--- D8_m ---*)
    {"\!\(\*SubscriptBox[\(Si\), \
\(3\)]\)\!\(\*SubscriptBox[\(W\), \(5\)]\)", 
     "\!\(\*SubscriptBox[\(Cr\), \(5\)]\)\!\(\*SubscriptBox[\(Ge\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Cr\), \(5\)]\)\!\(\*SubscriptBox[\(Si\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Mo\), \(5\)]\)\!\(\*SubscriptBox[\(Si\), \
\(3\)]\)", 
     "\!\(\*SubscriptBox[\(Si\), \(3\)]\)\!\(\*SubscriptBox[\(V\), \
\(5\)]\)", 
     "\!\(\*SubscriptBox[\(Si\), \(3\)]\)\!\(\*SubscriptBox[\(W\), \
\(5\)]\)"},
    (*--- D10_1 ---*)
    {"\!\(\*SubscriptBox[\(Cr\), \
\(7\)]\)\!\(\*SubscriptBox[\(C\), \(3\)]\)"},
    (*--- D10_2 ---*)
    {"\!\(\*SubscriptBox[\(Fe\), \
\(3\)]\)\!\(\*SubscriptBox[\(Th\), \(7\)]\)"},
    (*--- E0_1 ---*)
    {"PbFCl", "YoCl", "LaOCl", "LaOBr", "LaOJ", 
     "PrOCl", "BiFO", "BiClO"},
    (*--- E1_1 ---*)
    {"\!\(\*SubscriptBox[\(CuFeS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(AgAlS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(AlCuSe\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CuGaS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CuInS\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CuGaSe\), \(2\)]\)"},
    (*--- E2_1 ---*)
    {"\!\(\*SubscriptBox[\(CaTiO\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(PbTiO\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(PbZrO\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(BaTiO\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(PbSnO\), \(3\)]\)" \
"\!\(\*SubscriptBox[\(BiAlO\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(CdTiO\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(SrTiO\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(SrZnO\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(CaZrO\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(BaZrO\), \(3\)]\)"},
    (*--- E3 ---*)
    {"\!\(\*SubscriptBox[\(Al\), \
\(2\)]\)\!\(\*SubscriptBox[\(CdS\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(Al\), \(2\)]\)\!\(\*SubscriptBox[\(CdSe\), \
\(4\)]\)", 
     "\!\(\*SubscriptBox[\(CdGa\), \(2\)]\)\!\(\*SubscriptBox[\(Se\), \
\(4\)]\)" "\!\(\*SubscriptBox[\(Ga\), \(2\)]\)\!\(\*SubscriptBox[\(Te\
\), \(4\)]\)Zn"
     },
    (*--- E9_3 ---*)
    {"\!\(\*SubscriptBox[\(Fe\), \
\(3\)]\)\!\(\*SubscriptBox[\(W\), \(3\)]\)C", 
     "\!\(\*SubscriptBox[\(Fe\), \(3\)]\)\!\(\*SubscriptBox[\(Mo\), \
\(3\)]\)C", 
     "\!\(\*SubscriptBox[\(Co\), \(3\)]\)\!\(\*SubscriptBox[\(W\), \
\(3\)]\)C"},
    (*--- E9_a ---*)
    {"\!\(\*SubscriptBox[\(Al\), \
\(7\)]\)\!\(\*SubscriptBox[\(Cu\), \(2\)]\)Fe"},
    (*--- E9_b ---*)
    {"\!\(\*SubscriptBox[\(Al\), \
\(8\)]\)\!\(\*SubscriptBox[\(FeMg\), \
\(3\)]\)\!\(\*SubscriptBox[\(Si\), \(6\)]\)"},
    (*--- F0_1 ---*)
    {"CoAsS"},
    (*--- F5_1 ---*)
    {"\!\(\*SubscriptBox[\(NaHF\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(NaN\), \(3\)]\)", "NaCNO", 
     "\!\(\*SubscriptBox[\(CsJCl\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CaCN\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(CrS\), \(2\)]\)Na", 
     "\!\(\*SubscriptBox[\(CrS\), \(2\)]\)K", 
     "\!\(\*SubscriptBox[\(FeO\), \(2\)]\)Na"},
    (*--- H1_1 ---*)
    {"\!\(\*SubscriptBox[\(MgAl\), \
\(2\)]\)\!\(\*SubscriptBox[\(O\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(Mg\), \(2\)]\)\!\(\*SubscriptBox[\(TiO\), \
\(4\)]\)", 
     "\!\(\*SubscriptBox[\(Mg\), \(2\)]\)\!\(\*SubscriptBox[\(VO\), \
\(4\)]\)", 
     "\!\(\*SubscriptBox[\(MoO\), \(4\)]\)\!\(\*SubscriptBox[\(Na\), \
\(2\)]\)", 
     "\!\(\*SubscriptBox[\(WO\), \(4\)]\)\!\(\*SubscriptBox[\(Na\), \
\(2\)]\)", 
     "\!\(\*SubscriptBox[\(Mn\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \
\(4\)]\)Ti", 
     "\!\(\*SubscriptBox[\(Fe\), \(3\)]\)\!\(\*SubscriptBox[\(O\), \
\(4\)]\)"},
    (*--- H2_4 ---*)
    {"\!\(\*SubscriptBox[\(Cu\), \
\(3\)]\)\!\(\*SubscriptBox[\(VS\), \(4\)]\)"},
    (*--- L1_0 ---*)
    {"AuCu", "PdFe", "PtFe", "PtNi", "ZnPt", "CdPd", 
     "CdPt", "PtCo"},
    (*--- L1_2 ---*)
    {"\!\(\*SubscriptBox[\(Cu\), \(3\)]\)Au", 
     "\!\(\*SubscriptBox[\(Pt\), \(3\)]\)Ti", 
     "\!\(\*SubscriptBox[\(PtNi\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(Cu\), \(3\)]\)Pt", 
     "\!\(\*SubscriptBox[\(Zn\), \(3\)]\)Mn", 
     "\!\(\*SubscriptBox[\(Ni\), \(3\)]\)Mn", 
     "\!\(\*SubscriptBox[\(AuTi\), \(3\)]\)"},
    (*--- L2_1 ---*)
    {"\!\(\*SubscriptBox[\(Cu\), \(2\)]\)AlMn", 
     "\!\(\*SubscriptBox[\(Ni\), \(2\)]\)MgSb", 
     "\!\(\*SubscriptBox[\(InCu\), \(2\)]\)Mn", 
     "\!\(\*SubscriptBox[\(SnNi\), \(2\)]\)Mg", 
     "\!\(\*SubscriptBox[\(SnCu\), \(2\)]\)Mn"},
    (*--- L'2_b ---*)
    {"\!\(\*SubscriptBox[\(ThH\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(H\), \(2\)]\)Hf", 
     "\!\(\*SubscriptBox[\(H\), \(2\)]\)Zr"},
    (*--- L'_3 ---*)
    {"\!\(\*SubscriptBox[\(Fe\), \(2\)]\)N", 
     "\!\(\*SubscriptBox[\(W\), \(2\)]\)C", 
     "\!\(\*SubscriptBox[\(Co\), \(2\)]\)N", 
     "\!\(\*SubscriptBox[\(CV\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(Mn\), \(2\)]\)N"},
    (*--- L6_0 ---*)
    {"\!\(\*SubscriptBox[\(CuTi\), \(3\)]\)"}
    };
  If[tab,
     Print[Grid[Join[{head}, table], Background -> {None, {GTBackGroundColor1}}, Frame -> All]],           
     If[Head[struc] === Integer, 
  	    row = struc,  
        If[comp, 
     	   pos = Position[structures, struc] // Flatten, 
           pos = Position[table, struc] // Flatten];
           If[pos == {}, 
              Message[GTCrystalData::struc]; Return[],
              row = pos[[1]]
           ]
         ];
         info = {{"Strukturbericht", table[[row, 2]], "Prototype", 
                 table[[row, 4]]}, {"Citation", table[[row, 3]], 
                 "Pearson symbol", table[[row, 5]]}, {"GTPack.struc", 
                 table[[row, 7]], "Space group no.", table[[row, 6]]}
                };
         Print[Grid[info, Frame -> All, Background -> {{1 -> GTBackGroundColor1, 3 -> GTBackGroundColor1}, None}]]
  ];
  If[verb,
     Print["Structures: ", structures[[row]]];
     spc = Read[StringToStream[table[[row, 6]]], Number];
     ttt = GTSpaceGroups[spc, GOVerbose -> True]
  ]
]

(*
***)

(****e* /GTSpaceGroups
! NAME
!  GTSpaceGroups
! AUTHOR
!  W. Hergert
! PACKAGE
!   CrystalStructure.m
! MODIFICATION HISTORY
!  * 02.01.2014 : first version
!  * 15.01.2016 : new version corresponding to the new nomenclature
!  * 27.06.2018 : check header and documentation
! USAGE
!  GTSpaceGroups[input] gives the nomenclature of the 230 space groups. 
! INPUT
!  argument - space group number, international, HM or crystal system
!  
! OUTPUT
!  * if a single space group is in the INPUT all the different names and an information about
!    symmorphic/nonsymmorphic are printed and also given in  list.
!  * if the crystal system is INPUT, all space groups which belog to the crystal system are printed.

! GTPack OPTIONS
!  GOVerbose: 
!     - True  (standard) - print and list in output 
!     - False            - print only if crystal system in input
!  GOTable: 
!     - True             - complete table
!     - False            - no table
! STANDARD OPTIONS
!  -
! GTPack MODULES
! -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The names are implemented as symbols in Symbols.m. For the international nomenclature
!  bracketing bars are used.
! LITERATURE
!  * R.C. Powell : Symmetry, Group Theory and the Physical Properties of Crystals
!  * International Tables of Crystallography
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


GTSpaceGroups[arg_, 
  OptionsPattern[]] := 
 Module[{prt, i, j, symm,
   triclinicSF, monoclinicSF, orthorhombicSF, tetragonalSF, 
   trigonalSF, hexagonalSF, cubicSF, triclinicHM, monoclinicHM, 
   orthorhombicHM, tetragonalHM, trigonalHM, hexagonalHM, cubicHM,
   head,head1,system,sf,hm,spg,syst,num,sv,tab,type,index}, 
  prt = OptionValue[GOVerbose];
  head = {"Number" , "Schoenfliess", "International", "symmorphic"};
  head1 = {"Number" , "System", "Schoenfliess", "International", 
    "symmorphic"};
    
  (*--- Crystal systems ---*)
  system = {{"Triclinic", 1, 2}, {"Monoclinic", 3, 
     15}, {"Orthorhombic", 16, 74}, {"Tetragonal", 75, 
     142}, {"Trigonal", 143, 167}, {"Hexagonal", 168, 194}, {"Cubic", 
     195, 230}};
     
  (*--- symmorphic space groups ---*)
  symm = {1, 2, 3, 5, 6, 8, 10, 12, 16, 21, 22, 23, 25, 35, 38, 42, 
    44, 47, 65, 69, 71, 75, 79, 81, 82, 83, 87, 89, 97, 99, 107, 111, 
    115, 119, 121, 123, 139, 143, 146, 147, 148, 149, 150, 155, 156, 
    157, 160, 162, 164, 166, 168, 174, 175, 177, 183, 187, 189, 191, 
    195, 196, 197, 200, 202, 204, 207, 209, 211, 215, 216, 217, 221, 
    225, 229};

  (*--- definition of space groups symbols SF ---*)
  triclinicSF = {C11, Ci1};
  monoclinicSF = {C21, C22, C23, Cs1, Cs2, Cs3, Cs4, C2h1, C2h2, C2h3,
     C2h4, C2h5, C2h6};
  orthorhombicSF = {D21, D22, D23, D24, D25, D26, D27, D28, D29, 
    C2v1,
    C2v2, C2v3, C2v4, C2v5, C2v6, C2v7, C2v8, C2v9, C2v10, C2v11, 
    C2v12, C2v13, C2v14, C2v15, C2v16, C2v17, C2v18, C2v19, C2v20, 
    C2v21, C2v22, D2h1, D2h2, D2h3, D2h4, D2h5, D2h6, D2h7, D2h8, 
    D2h9, D2h10,
    D2h11, D2h12, D2h13, D2h14, D2h15, D2h16, D2h17, D2h18, D2h19, 
    D2h20, D2h21, D2h22, D2h23, D2h24, D2h25, D2h26, D2h27, D2h28};
  tetragonalSF = {C41, C42, C43, C44, C45, C46, S41, S42, C4h1, C4h2,
    C4h3, C4h4, C4h5, C4h6, D41, D42, D43, D44, D45, D46, D47, D48, 
    D49, D410, C4v1, C4v2, C4v3, C4v4, C4v5, C4v6, C4v7, C4v8, C4v9, 
    C4v10, C4v11, C4v12, D2d1, D2d2, D2d3, D2d4,
    D2d5, D2d6, D2d7, D2d8, D2d9, D2d10, D2d11, D2d12, D4h1, D4h2, 
    D4h3, D4h4, D4h5, D4h6, D4h7, D4h8, D4h9, D4h10, D4h11, D4h12, 
    D4h13, D4h14, D4h15, D4h16, D4h17, D4h18, D4h19, D4h20};
  trigonalSF = {C31, C32, C33, C34, C3i1, C3i2, D31, D32, D33, D34, 
    D35, D36, D37, C3v1, C3v2, C3v3, C3v4, C3v5, C3v6, D3d1, D3d2, 
    D3d3, D3d4, D3d5, D3d6};
  hexagonalSF = {C61, C62, C63, C64, C65, C66, C3h1, C6h1, C6h2, D61, 
    D62, D63, D64, D65, D66, C6v1, C6v2, C6v3, C6v4, D3h1, D3h2, D3h3,
     D3h4, D6h1, D6h2, D6h3, D6h4};
  cubicSF = {T1, T2, T3, T4, T5, Th1, Th2, Th3, Th4, Th5, Th6, Th7, 
    O1, O2, O3, O4, O5, O6, O7, O8, Td1, Td2, Td3, Td4, Td5, Td6, Oh1,
     Oh2, Oh3, Oh4, Oh5, Oh6, Oh7, Oh8, Oh9, Oh10};

  (*--- definition of space groups symbols HM ---*)
  triclinicHM = {
  	\[LeftBracketingBar]P1   \[RightBracketingBar],
  	\[LeftBracketingBar]P - 1\[RightBracketingBar]
  };
  monoclinicHM = {
  	\[LeftBracketingBar]P2   \[RightBracketingBar],
  	\[LeftBracketingBar]P2_ 1\[RightBracketingBar],
  	\[LeftBracketingBar]C2   \[RightBracketingBar],
  	\[LeftBracketingBar]Pm   \[RightBracketingBar],
  	\[LeftBracketingBar]Pc   \[RightBracketingBar],
  	\[LeftBracketingBar]Cm   \[RightBracketingBar],
  	\[LeftBracketingBar]Cc   \[RightBracketingBar],
  	\[LeftBracketingBar]P2/ m\[RightBracketingBar],
  	\[LeftBracketingBar]P2_1/m\[RightBracketingBar],
  	\[LeftBracketingBar]C2/ m\[RightBracketingBar],
  	\[LeftBracketingBar]P2/ c\[RightBracketingBar],
  	\[LeftBracketingBar]P2_ 1/c\[RightBracketingBar],
  	\[LeftBracketingBar]C2/ c\[RightBracketingBar]
  };
  orthorhombicHM = {
  	\[LeftBracketingBar]P222\[RightBracketingBar],
  	\[LeftBracketingBar]P222_ 1\[RightBracketingBar],
  	\[LeftBracketingBar]P2_ 12 _ 12\[RightBracketingBar],
  	\[LeftBracketingBar]P2_ 12 _ 12 _ 1\[RightBracketingBar],
  	\[LeftBracketingBar]C222_ 1\[RightBracketingBar],
  	\[LeftBracketingBar]C222\[RightBracketingBar],
  	\[LeftBracketingBar]F222\[RightBracketingBar],
  	\[LeftBracketingBar]I222\[RightBracketingBar],
  	\[LeftBracketingBar]I2_ 12 _ 12 _ 1\[RightBracketingBar],
  	\[LeftBracketingBar]Pmm2\[RightBracketingBar],
    \[LeftBracketingBar]Pmc2_ 1\[RightBracketingBar],
    \[LeftBracketingBar]Pcc2\[RightBracketingBar],
    \[LeftBracketingBar]Pma2\[RightBracketingBar], 
    \[LeftBracketingBar]Pca2_ 1\[RightBracketingBar],
    \[LeftBracketingBar]Pnc2\[RightBracketingBar],
    \[LeftBracketingBar]Pmn2_ 1\[RightBracketingBar],
    \[LeftBracketingBar]Pba2\[RightBracketingBar],
    \[LeftBracketingBar]Pna2_ 1\[RightBracketingBar],
    \[LeftBracketingBar]Pnn2\[RightBracketingBar],
    \[LeftBracketingBar]Cmm2\[RightBracketingBar],
    \[LeftBracketingBar]Cmc2_ 1\[RightBracketingBar],
    \[LeftBracketingBar]Ccc2\[RightBracketingBar],
    \[LeftBracketingBar]Amm2\[RightBracketingBar],
    \[LeftBracketingBar]Aem2\[RightBracketingBar],
    \[LeftBracketingBar]Ama2\[RightBracketingBar],
    \[LeftBracketingBar]Aea2\[RightBracketingBar],
    \[LeftBracketingBar]Fmm2\[RightBracketingBar],
    \[LeftBracketingBar]Fdd2\[RightBracketingBar],
    \[LeftBracketingBar]Imm2\[RightBracketingBar],
    \[LeftBracketingBar]Iba2\[RightBracketingBar],
    \[LeftBracketingBar]Ima2\[RightBracketingBar],
    \[LeftBracketingBar]Pmmm\[RightBracketingBar],
    \[LeftBracketingBar]Pnnn\[RightBracketingBar],
    \[LeftBracketingBar]Pccm\[RightBracketingBar],
    \[LeftBracketingBar]Pban\[RightBracketingBar],
    \[LeftBracketingBar]Pmma\[RightBracketingBar],
    \[LeftBracketingBar]Pnna\[RightBracketingBar],
    \[LeftBracketingBar]Pmna\[RightBracketingBar],
    \[LeftBracketingBar]Pcca\[RightBracketingBar],
    \[LeftBracketingBar]Pbam\[RightBracketingBar],
    \[LeftBracketingBar]Pccn\[RightBracketingBar],
    \[LeftBracketingBar]Pbcm\[RightBracketingBar],
    \[LeftBracketingBar]Pnnm\[RightBracketingBar],
    \[LeftBracketingBar]Pmmn\[RightBracketingBar],
    \[LeftBracketingBar]Pbcn\[RightBracketingBar],
    \[LeftBracketingBar]Pbca\[RightBracketingBar],
    \[LeftBracketingBar]Pnma\[RightBracketingBar],
    \[LeftBracketingBar]Cmcm\[RightBracketingBar],
    \[LeftBracketingBar]Cmce\[RightBracketingBar],
    \[LeftBracketingBar]Cmmm\[RightBracketingBar],
    \[LeftBracketingBar]Cccm\[RightBracketingBar],
    \[LeftBracketingBar]Cmme\[RightBracketingBar],
    \[LeftBracketingBar]Ccce\[RightBracketingBar],
    \[LeftBracketingBar]Fmmm\[RightBracketingBar],
    \[LeftBracketingBar]Fddd\[RightBracketingBar],
    \[LeftBracketingBar]Immm\[RightBracketingBar],
    \[LeftBracketingBar]Ibam\[RightBracketingBar],
    \[LeftBracketingBar]Ibca\[RightBracketingBar],
    \[LeftBracketingBar]Imma\[RightBracketingBar]
  };
  tetragonalHM = {
  	\[LeftBracketingBar]P4\[RightBracketingBar],
  	\[LeftBracketingBar]P4_ 1\[RightBracketingBar],
  	\[LeftBracketingBar]P4_ 2\[RightBracketingBar],
  	\[LeftBracketingBar]P4_ 3\[RightBracketingBar],
  	\[LeftBracketingBar]I4\[RightBracketingBar],
  	\[LeftBracketingBar]I4_ 1\[RightBracketingBar],
  	\[LeftBracketingBar]P - 4\[RightBracketingBar],
  	\[LeftBracketingBar]I - 4\[RightBracketingBar],
  	\[LeftBracketingBar]P4/ m\[RightBracketingBar],
  	\[LeftBracketingBar]P4_ 2/ m\[RightBracketingBar],
    \[LeftBracketingBar]P4/ n\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2/ n\[RightBracketingBar],
    \[LeftBracketingBar]I4/ m\[RightBracketingBar],
    \[LeftBracketingBar]I4_ 1/ a\[RightBracketingBar],
    \[LeftBracketingBar]P422\[RightBracketingBar],
    \[LeftBracketingBar]P42_ 12\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 122\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 12 _ 12\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 222\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 22 _ 12\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 322\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 32 _ 12\[RightBracketingBar],
    \[LeftBracketingBar]I422\[RightBracketingBar],
    \[LeftBracketingBar]I4_ 122\[RightBracketingBar],
    \[LeftBracketingBar]P4mm\[RightBracketingBar],
    \[LeftBracketingBar]P4bm\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2 cm\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2 nm\[RightBracketingBar],
    \[LeftBracketingBar]P4cc\[RightBracketingBar],
    \[LeftBracketingBar]P4nc\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2 mc\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2 bc\[RightBracketingBar],
    \[LeftBracketingBar]I4mm\[RightBracketingBar],
    \[LeftBracketingBar]I4cm\[RightBracketingBar],
    \[LeftBracketingBar]I4_ 1 md\[RightBracketingBar],
    \[LeftBracketingBar]I4_ 1 cd\[RightBracketingBar],
    \[LeftBracketingBar]P -42 m\[RightBracketingBar],
    \[LeftBracketingBar]P -42 c\[RightBracketingBar],
    \[LeftBracketingBar]P -42 _ 1 m\[RightBracketingBar],
    \[LeftBracketingBar]P -42 _ 1 c\[RightBracketingBar],
    \[LeftBracketingBar]P -4 m2\[RightBracketingBar],
    \[LeftBracketingBar]P -4 c2\[RightBracketingBar],
    \[LeftBracketingBar]P -4 b2\[RightBracketingBar],
    \[LeftBracketingBar]P -4 n2\[RightBracketingBar],
    \[LeftBracketingBar]I -4 m2\[RightBracketingBar],
    \[LeftBracketingBar]I -4 c2\[RightBracketingBar],
    \[LeftBracketingBar]I -42 m\[RightBracketingBar],
    \[LeftBracketingBar]I -42 d\[RightBracketingBar],
    \[LeftBracketingBar]P4/mmm\[RightBracketingBar],
    \[LeftBracketingBar]P4/mcc\[RightBracketingBar],
    \[LeftBracketingBar]P4/nbm\[RightBracketingBar],
    \[LeftBracketingBar]P4/nnc\[RightBracketingBar],
    \[LeftBracketingBar]P4/mbm\[RightBracketingBar],
    \[LeftBracketingBar]P4/mnc\[RightBracketingBar],
    \[LeftBracketingBar]P4/nmm\[RightBracketingBar],
    \[LeftBracketingBar]P4/ncc\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2/mmc\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2/mcm\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2/nbc\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2/nnm\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2/mbc\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2/mnm\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2/nmc\[RightBracketingBar],
    \[LeftBracketingBar]P4_ 2/ncm\[RightBracketingBar],
    \[LeftBracketingBar]I4/  mmmm\[RightBracketingBar],
    \[LeftBracketingBar]I4/   mcm\[RightBracketingBar],
    \[LeftBracketingBar]I4_ 1/amd\[RightBracketingBar],
    \[LeftBracketingBar]I4_ 1/acd\[RightBracketingBar]
  };
  trigonalHM = {
  	\[LeftBracketingBar]P3\[RightBracketingBar],
  	\[LeftBracketingBar]P3_ 1\[RightBracketingBar],
  	\[LeftBracketingBar]P3_ 2\[RightBracketingBar],
  	\[LeftBracketingBar]R3\[RightBracketingBar],
  	\[LeftBracketingBar]P - 3\[RightBracketingBar],
  	\[LeftBracketingBar]R - 3\[RightBracketingBar],
  	\[LeftBracketingBar]P312\[RightBracketingBar],
  	\[LeftBracketingBar]P321\[RightBracketingBar],
  	\[LeftBracketingBar]P3_ 112\[RightBracketingBar],
  	\[LeftBracketingBar]P3_ 121\[RightBracketingBar],
  	\[LeftBracketingBar]P3_ 212\[RightBracketingBar],
  	\[LeftBracketingBar]P3_ 221\[RightBracketingBar],
  	\[LeftBracketingBar]R32\[RightBracketingBar],
  	\[LeftBracketingBar]P3m1\[RightBracketingBar],
  	\[LeftBracketingBar]P31m\[RightBracketingBar],
  	\[LeftBracketingBar]P3c1\[RightBracketingBar],
  	\[LeftBracketingBar]P31c\[RightBracketingBar],
  	\[LeftBracketingBar]R3m\[RightBracketingBar],
  	\[LeftBracketingBar]R3c\[RightBracketingBar],
  	\[LeftBracketingBar]P -31 m\[RightBracketingBar],
  	\[LeftBracketingBar]P -31 c\[RightBracketingBar],
  	\[LeftBracketingBar]P -3 m1\[RightBracketingBar],
  	\[LeftBracketingBar]P -3 c1\[RightBracketingBar],
  	\[LeftBracketingBar]R -3 m\[RightBracketingBar],
  	\[LeftBracketingBar]R -3 c\[RightBracketingBar]
  };
  hexagonalHM = {
  	\[LeftBracketingBar]P6\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 1\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 5\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 2\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 4\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 3\[RightBracketingBar],
  	\[LeftBracketingBar]P - 6\[RightBracketingBar],
  	\[LeftBracketingBar]P6/ m\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 3/m\[RightBracketingBar],
  	\[LeftBracketingBar]P622\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 122\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 522\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 222\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 422\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 322\[RightBracketingBar],
  	\[LeftBracketingBar]P6mm\[RightBracketingBar],
  	\[LeftBracketingBar]P6cc\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 3 cm\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 3 mc\[RightBracketingBar],
  	\[LeftBracketingBar]P - 6 m2\[RightBracketingBar],
  	\[LeftBracketingBar]P - 6 c2\[RightBracketingBar],
  	\[LeftBracketingBar]P - 62 m\[RightBracketingBar],
  	\[LeftBracketingBar]P - 62 c\[RightBracketingBar],
  	\[LeftBracketingBar]P6/ mmm\[RightBracketingBar],
  	\[LeftBracketingBar]P6/ mcc\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 3/ mcm\[RightBracketingBar],
  	\[LeftBracketingBar]P6_ 3/ mmc\[RightBracketingBar]
  };
  cubicHM = {
  	\[LeftBracketingBar]P23\[RightBracketingBar],
  	\[LeftBracketingBar]F23\[RightBracketingBar],
  	\[LeftBracketingBar]I23\[RightBracketingBar],
  	\[LeftBracketingBar]P2_ 13\[RightBracketingBar],
  	\[LeftBracketingBar]I2_ 13\[RightBracketingBar],
  	\[LeftBracketingBar]Pm -3\[RightBracketingBar],
  	\[LeftBracketingBar]Pn -3\[RightBracketingBar],
  	\[LeftBracketingBar]Fm -3\[RightBracketingBar],
  	\[LeftBracketingBar]Fd -3\[RightBracketingBar],
  	\[LeftBracketingBar]Im -3\[RightBracketingBar],
  	\[LeftBracketingBar]Pa -3\[RightBracketingBar],
  	\[LeftBracketingBar]Ia -3\[RightBracketingBar],
  	\[LeftBracketingBar]P432\[RightBracketingBar],
  	\[LeftBracketingBar]P4_ 232\[RightBracketingBar],
  	\[LeftBracketingBar]F432\[RightBracketingBar],
  	\[LeftBracketingBar]F4_ 132\[RightBracketingBar],
  	\[LeftBracketingBar]I432\[RightBracketingBar],
  	\[LeftBracketingBar]P4_ 332\[RightBracketingBar],
  	\[LeftBracketingBar]P4_ 132\[RightBracketingBar],
  	\[LeftBracketingBar]I4_ 132\[RightBracketingBar],
  	\[LeftBracketingBar]P -43 m\[RightBracketingBar],
  	\[LeftBracketingBar]F -43 m\[RightBracketingBar],
  	\[LeftBracketingBar]I -43 m\[RightBracketingBar],
  	\[LeftBracketingBar]P -43 n\[RightBracketingBar],
  	\[LeftBracketingBar]F -43 c\[RightBracketingBar],
  	\[LeftBracketingBar]I -43 d\[RightBracketingBar],
  	\[LeftBracketingBar]Pm -3 m\[RightBracketingBar],
  	\[LeftBracketingBar]Pn -3 n\[RightBracketingBar],
  	\[LeftBracketingBar]Pm -3 n\[RightBracketingBar],
  	\[LeftBracketingBar]Pn -3 m\[RightBracketingBar],
  	\[LeftBracketingBar]Fm -3 m\[RightBracketingBar],
  	\[LeftBracketingBar]Fm -3 c\[RightBracketingBar],
  	\[LeftBracketingBar]Fd -3 m\[RightBracketingBar],
  	\[LeftBracketingBar]Fd -3 c\[RightBracketingBar],
  	\[LeftBracketingBar]Im -3 m\[RightBracketingBar],
  	\[LeftBracketingBar]Ia -3 d\[RightBracketingBar]
  };
  (*--- combine all the notations ---*)
  sf = {triclinicSF, monoclinicSF, orthorhombicSF, tetragonalSF, 
     trigonalSF, hexagonalSF, cubicSF} // Flatten;
  hm = {triclinicHM, monoclinicHM, orthorhombicHM, tetragonalHM, 
     trigonalHM, hexagonalHM, cubicHM} // Flatten;
  spg = Transpose[{sf, hm}];syst = {};
  (*--- List of crystal systems is constructed ---*)
  Do[
     Do[
        syst = Append[syst, system[[i, 1]]]
     , {j, system[[i, 2]], system[[i, 3]]}]
  , {i, 1, 7}];
  (*--- Output of a complete table ---*)
  If[OptionValue[GOTable],
     num = Table[i, {i, 1, 230}]; 
     sv = Table["no", {230}]; 
     Do[
     	sv[[symm[[i]]]] = "yes"
     , {i, 1, Length[symm]}];
     spg = Transpose[{num, syst, sf, hm, sv}]; 
     Print[ Grid[Join[{head1}, spg], Background -> {None, {GTBackGroundColor1}}, Frame -> All]]; Return[],
     None
  ];
  (*--- input is crystal system ---*)
  If[Head[arg] === String,  
     tab = {};     
     Do[
        If[arg === system[[i, 1]],
           Do[
              If[Intersection[{j}, symm] === {},
                 tab = Append[tab, {j, spg[[j, 1]], spg[[j, 2]], "no"}],
                 tab = Append[tab, {j, spg[[j, 1]], spg[[j, 2]], "yes"}]
              ]   
            , {j, system[[i, 2]], system[[i, 3]]}]
            , None
        ]
     , {i, 1, 7}];
     If[Length[tab]>0,Print[ Grid[Join[{head}, tab], Background -> {None, {GTBackGroundColor1}}, Frame -> All]]];
     Return[],
     None
  ];
  (*--- input is space group number ---*)
  If[Head[arg] === Integer,
     If[Intersection[{arg}, symm] == {},
        type = "no",
        type = "yes"
     ];
     tab = {arg, syst[[arg]], spg[[arg, 1]], spg[[arg, 2]], type};
     If[prt,
       Print[ Grid[Join[{head1}, {tab}], Background -> {None, {GTBackGroundColor1}}, Frame -> All]],
       None
     ];
     Return[tab]
  ];
  (*--- input is symbol or bracketingbar ---*)
  If[Head[arg] === Symbol,
     index = 1,
     If[Head[arg] === BracketingBar,
        index = 2,
        Print["GTSpaceGroups Error : Wrong Head in notation"]; 
        Abort[]
     ]
  ];
  Do[
     If[arg === spg[[i, index]],
       If[Intersection[{i}, symm] == {},
          type = "no",
          type = "yes"
       ];
       tab = {i, syst[[i]], spg[[i, 1]], spg[[i, 2]], type};
       If[prt,
          Print[ Grid[Join[{head1}, {tab}], Background -> {None, {GTBackGroundColor1}}, Frame -> All]],
          None
       ];
       Return[tab],
       None
    ]
  , {i, 1, 230}]
  ]
  
 
(*
***)


(****e* /GTBravaisLattice
! NAME
!  GTBravaiLattice
! AUTHOR
!  W. Hergert
! PACKAGE
!   CrystalStructure.m
! MODIFICATION HISTORY
!  * 01.02.2014 : first version
!  * 27.06.2018 : check header and documentation
!  * 08.03.2023
! USAGE
!  GTBravaisLattice[crystal system] gives the lattice vectors of the crystal system.
!   
! INPUT
!  crystal system - in form {name,type} like : {"Cubic","F"}
!  
! OUTPUT
!  - image of the Bravais lattice
!  - Lattice vectors and basis
! 
! GTPack OPTIONS
!  GOBravais:  
!     - "Conventional"   - conventional cell with basis (standard)
!     - "Primitive"      - primitive cells with one atom per cell
!  GOData:  
!     - "Image"          - Image of the conventional cell
!     - "Basis"          - Lattice vectors (standard)
!  GOImage:              
     -  {a->1,b->1.5,c->2,alph->pi/3,bet->pi/4,gam->2 pi/3} some standard values for GTCellPlot
! STANDARD OPTIONS
!  FontSize: 20          - fontsize coordinate system
!  Spheresize
! GTPack MODULES
!   GTLatticeVectors, GTCellPlot
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The difficulty here is to guarantee, that the names will be understood correctly. May be it is a bit 
!  complicated now, but it should work.
! LITERATURE
!  -
! TODO
!  extend to the twodimensional lattices
!  tests!
!  documentation page
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTBravaisLattice::inp  = "Structure not correct or old input style for GTLatticeVectors not allowed";
GTBravaisLattice::cell = "GOBravais has to be Conventional or Primitive.";
GTBravaisLattice::out  = "Output option has to be Data or Image.";


GTBravaisLattice[grp_, pldat_, OptionsPattern[]] := Module[
  {gbr,plt,fnt,cs,rule,latt,v1,v2,v3,bas,basc,grp1,t,t1,p,vec,bas1},
   Clear[a, b, c, \[Alpha], \[Beta], \[Gamma], v1, v2, v3];
 (*--- options ---*)
   gbr  = OptionValue[GOBravais];
   plt  = OptionValue[GOData];
   fnt  = OptionValue[FontSize];
   cs   = OptionValue[GOCoordinateSystem];
   rule = OptionValue[GOImage];
 (*---Check the options---*)
   If[gbr == "Conventional" || gbr == "Primitive",
      None,
      Message[GTBravaisLattice::cell]; Return[]
   ];
   If[plt == "Data" || plt == "Image",
      None,
      Message[GTBravaisLattice::out]; Return[]
   ];
 (*--- lattice types ---*)
   latt = {{"Triclinic",    "aP"}, 
  	      {"Monoclinic",   "mP"}, 
  	      {"Monoclinic",   "mC"}, 
  	      {"Orthorhombic", "oP"}, 
  	      {"Orthorhombic", "oC"}, 
  	      {"Orthorhombic", "oI"}, 
  	      {"Orthorhombic", "oF"}, 
  	      {"Tetragonal",   "tP"}, 
  	      {"Tetragonal",   "tI"}, 
  	      {"Hexagonal",    "hP"}, 
  	      {"Rhombohedral", "hR"}, 
  	      {"Cubic",        "cP"}, 
  	      {"Cubic",        "cF"}, 
  	      {"Cubic",        "cI"}
  	     };
 (*--- information for conventional cell ---*)
   basc = {{}, 
   	       {}, 
   	       {(v1 + v2)/2, (v1 + v2)/2 + v3}, 
   	       {}, 
   	       {(v1 + v2)/2, (v1 + v2)/2 + v3},
   	       {(v1 + v2 + v3)/2}, 
   	       {(v1 + v2)/2, (v1 + v3)/2, (v2 + v3)/2, (v1 + v2)/2 + v3, (v2 + v3)/2 + 
            v1, (v1 + v3)/2 + v2}, {}, {(v1 + v2 + v3)/2}, 
           {}, 
           {}, 
           {}, 
           {(v1 + v2)/2, (v1 + v3)/2, (v2 + v3)/2, (v1 + v2)/2 + v3, (v2 + v3)/2 + v1, (v1 + v3)/2 + v2}, 
           {(v1 + v2 + v3)/2}
          };
 (*---Test input---*)
   If[Intersection[latt, {grp}] == {},
      Message[GTBravaisLattice::inp]; Return[],
      None
   ];
 (*--- Construct the cells ---*)  
   If[gbr == "Primitive",
      {v1, v2, v3} = GTLatticeVectors[grp];
       bas   = {},
       grp1  = grp;
       t     = Characters[grp1[[2]]];
       If[grp1[[1]] == "Rhombohedral",
          t1 = "R",
          t1 = "P"
       ];
       grp1[[2]]    = StringJoin[{t[[1]], t1}];
       {v1, v2, v3} = GTLatticeVectors[grp1];
   p = Flatten[Position[latt, grp]][[1]];
   bas = basc[[p]]
   ];
 (*--- outpu vectors or plot ---*)  
   If[plt == "Data",
      Return[{{v1, v2, v3}, bas}],
      vec = {v1, v2, v3} /. rule;
      bas1 = bas /. rule;
      Print[GTCellPlot[vec, bas1, pldat, FontSize -> fnt, GOCoordinateSystem -> cs]]
   ]
]
  
(*
GTBravaisLattice::inp  = "Structure not correct, old input style for GTLatticeVectors not allowed";
GTBravaisLattice::cell = "GOBravais has to be Conventional or Primitive.";
GTBravaisLattice::conv = "Primitive lattice not in list.";

GTBravaisLatticeT[grp_, OptionsPattern[]] := Module[
  {v1, v2, v3, grp1, gbr, latt, bas, p},
  (* Clear a,b,c,\[Alpha],\[Beta],\[Gamma] *)
  
  Clear[a, b, c, \[Alpha], \[Beta], \[Gamma], v1, v2, v3];
  gbr  = OptionValue[GOBravais];
  plt  = OptionValue[GOPlot];
  fnt  = OptionValue[FontSize];
  sph  = OptionValue[SphereSize];
  cs   = OptionValue[CoordSystem];
  rule = OptionValue[GOImage];   
  (*---Check the options---*)
  
  If[gbr == "Conventional" || gbr == "Primitive",
       None,
       Message[GTBravaisLattice::cell]; Return[]
   ];
  (* lattice types and the corresponding space groups *)
  
  latt = {{"Triclinic", "aP"}, {"Monoclinic", "mP"}, {"Monoclinic", 
     "mC"}, {"Orthorhombic", "oP"}, {"Orthorhombic", 
     "oC"}, {"Orthorhombic", "oI"}, {"Orthorhombic", 
     "oF"}, {"Tetragonal", "tP"}, {"Tetragonal", "tI"}, {"Hexagonal", 
     "hP"}, {"Rhombohedral", "hR"}, {"Cubic", "cP"}, {"Cubic", 
     "cF"}, {"Cubic", "cI"}};
  (* information for conventional cell *)
  
  basc = {{}, {}, {(v1 + v2)/2, (v1 + v2)/2 + v3}, {}, {(v1 + v2)/
      2, (v1 + v2)/2 + v3}, {(v1 + v2 + v3)/2}, {(v1 + v2)/
      2, (v1 + v3)/2, (v2 + v3)/2, (v1 + v2)/2 + v3, (v2 + v3)/2 + 
      v1, (v1 + v3)/2 + v2}, {}, {(v1 + v2 + v3)/
      2}, {}, {}, {}, {(v1 + v2)/2, (v1 + v3)/2, (v2 + v3)/
      2, (v1 + v2)/2 + v3, (v2 + v3)/2 + v1, (v1 + v3)/2 + 
      v2}, {(v1 + v2 + v3)/2}};
  (*---Test input---*)
     If[Intersection[latt, {grp}] == {},
          Message[GTBravaisLattice::inp]; Return[],
         None
     ];
  If[gbr == "Primitive",
   {v1, v2, v3} = GTLatticeVectors[grp];
   bas = {},
   grp1 = grp;
   t = Characters[grp1[[2]]];
   If[grp1[[1]] == "Rhombohedral",
    t1 = "R",
    t1 = "P"
    ];
   grp1[[2]] = StringJoin[{t[[1]], t1}];
   {v1, v2, v3} = GTLatticeVectors[grp1];
   p = Flatten[Position[latt, grp]][[1]]; Print[p];
   bas = basc[[p]]; Print[bas]
   ];
  Return[{{v1, v2, v3}, bas}]
  ]
(*
GTBravaisLattice::inp = "Structure not correct. Allowed are `1` "
GTBravaisLattice::gbr = "GOBravais has to be Conventional or Primitive."
GTBravaisLattice::out = "Output type has to be Image or Basis."

GTBravaisLattice[grp_,OptionsPattern[]] := Module[{pc,gbr,system,im,basc,basp,nn,i,ltype,grp1,num,syst,ln,tt,gc,tb,bvec,lat,lb,bv,bn,tl},
	                                           (*    a,b,c,\[Alpha], \[Beta], \[Gamma]}, *)
	                                           Print["under revision"];
	 pc = OptionValue[GOData]; gbr = OptionValue[GOBravais];
     (*--- Check the options ---*)
     If[gbr == "Conventional" || gbr == "Primitive", 
     	None, 
        Message[GTBravaisLattice::gbr]; Return[]
     ];
     If[pc == "Image" || pc == "Basis", 
     	None, 
        Message[GTBravaisLattice::out]; Return[]
     ]; 
     (*--- data base ---*)
     system = {{"Triclinic", 1, 2}, {"Monoclinic", 3,15}, {"Orthorhombic", 16, 74}, {"Tetragonal", 75, 142}, 
     	       {"Trigonal", 143, 167}, {"Hexagonal", 168, 194}, {"Cubic", 195, 230}};
     im = {{{"P", "SimpleTriclinic"}}, 
     	   {{"P", "SimpleMonoclinic"}, {"B", "BaseCenteredMonoclinic"}}, 
     	   {{"P", "SimpleOrthorhombic"}, {"I", "BodyCenteredOrthorhombic"}, 
     	   	{"C", "BaseCenteredOrthorhombic"}, {"F","FaceCenteredOrthorhombic"}}, 
     	   {{"P", "SimpleTetragonal"}, {"I","CenteredTetragonal"}}, 
     	   {{"R", "SimpleTrigonal"}}, 
     	   {{"P", "SimpleHexagonal"}}, 
     	   {
     	   	{"P", "SimpleCubic"}, 
     	    {"I","BodyCenteredCubic"}, 
     	    {"F", "FaceCenteredCubic"}
     	   }
     	  };
     basc = {{{"P", {{a, b, c, \[Alpha], \[Beta], \[Gamma]},{{0,0,0}}}}},
             {{"P", {{a, b, c, \[Pi]/2, \[Beta], \[Pi]/2},{{0,0,0}}}}, {"B", {{a, b, c, \[Pi]/2, \[Beta], \[Pi]/2}, {{0, 0, 0}, {1/2, 1/2, 0}}}}},
             {{"P", {{a, b, c, \[Pi]/2, \[Pi]/2, \[Pi]/2},{{0,0,0}}}}, {"C", {{a, b, c, \[Pi]/2, \[Pi]/2, \[Pi]/2}, {{0, 0, 0}, {1/2, 1/2, 0}}}}, 
              {"I", {{a, b, c, \[Pi]/2, \[Pi]/2, \[Pi]/2}, {{0, 0, 0}, {1/2, 1/2, 1/2}}}},
              {"F", {{a, b, c, \[Pi]/2, \[Pi]/2, \[Pi]/2}, {{0, 0, 0}, {1/2, 1/2, 0}, {0, 1/2, 1/2}, {1/2, 0, 1/2}}}}},
             {{"P", {{a, a, c, \[Pi]/2, \[Pi]/2, \[Pi]/2},{{0,0,0}}}}, {"I", {{a, a, c, \[Pi]/2, \[Pi]/2, \[Pi]/2}, {{0, 0, 0}, {1/2, 1/2, 1/2}}}}},
             {{"R", {{a, a, a, \[Alpha], \[Alpha], \[Alpha]},{{0,0,0}}}}},
             {{"P", {{a, a, c, \[Pi]/2, \[Pi]/2, 2 \[Pi]/3},{{0,0,0}}}}},
             {
              {"P", {{a, a, a, \[Pi]/2, \[Pi]/2, \[Pi]/2},{{0,0,0}}}}, 
              {"I", {{a, a, a, \[Pi]/2, \[Pi]/2, \[Pi]/2}, {{0, 0, 0}, {1/2, 1/2,  1/2}}}}, 
              {"F", {{a, a, a, \[Pi]/2, \[Pi]/2, \[Pi]/2}, {{0, 0, 0}, {1/2, 1/2, 0},{0, 1/2, 1/2}, {1/2, 0, 1/2}}}}
             }
            };
     basp = {{{"P", {a, b, c, \[Alpha], \[Beta], \[Gamma]}}},
             {{"P", {a, b, c, \[Pi]/2, \[Beta], \[Pi]/2}}, 
              {"B", {{a/2, b/2, 0}, {-a/2, b/2, 0}, {a c Cos[\[Beta]]/Sqrt[a^2 + b^2], b c/Sqrt[a^2 + b^2], c Sin[\[Beta]]}}}},
             {{"P", {a, b, c, \[Pi]/2, \[Pi]/2, \[Pi]/2}}, {"F", {{0, b/2, c/2}, {a/2, 0, c/2}, {a/2, b/2, 0}}}, 
              {"I", {{-a/2, b/2, c/2}, {a/2, -b/2, c/2}, {a/2, b/2, -c/2}}},
              {"C", {{a/2, b/2, 0}, {-a/2, b/2, 0}, {0, 0, c}}}},
             {{"P", {a, a, c, \[Pi]/2, \[Pi]/2, \[Pi]/2}}, {"I", {{-a/2, a/2, c/2}, {a/2, -a/2, c/2}, {a/2, a/2, -c/2}}}},
             {{"R", {a, a, a, \[Alpha], \[Alpha], \[Alpha]}}},
             {{"P", {a, a, c, \[Pi]/2, \[Pi]/2, 2 \[Pi]/3}}},
             {
              {"P", {a, a, a, \[Pi]/2, \[Pi]/2, \[Pi]/2}}, 
              {"I", {{-a/2, a/2, a/2}, {a/2, -a/2, a/2}, {a/2, a/2, -a/2}}}, 
              {"F", {{0, a/2, a/2}, {a/2, 0, a/2}, {a/2, a/2, 0}}}} 
             };
  
  If[Head[grp] === List,
    (*--- Input is crystal system and lattice type ---*)
    (*--- Test input ---*)
    tl= {{"Triclinic","P"}, {"Monoclinic","P"}, {"Monoclinic","B"}, {"Orthorhombic","P"}, 
    	{"Orthorhombic","C"}, {"Orthorhombic","I"}, {"Orthorhombic","F"}, {"Tetragonal","P"}, 
    	{"Tetragonal","I"}, {"Trigonal","R"}, {"Hexagonal","P"}, {"Cubic","P"}, {"Cubic","F"}, 
    	{"Cubic","I"}};    	
    If[Intersection[tl,{grp}]=={},
        Message[GTBravaisLattice::inp,tl]; Return[],
    	None
    ];	     	
    Do[
    	If[grp[[1]] == system[[i, 1]],
    	  nn = i, 
    	  None
       ]
   , {i, 1, 7}]; 
   ltype = grp[[2]],
   (*--- Input is space group spezification ---*)
   grp1 = GTSpaceGroups[grp, GOVerbose -> True]; num = grp1[[1]];
   Do[
   	  If[num \[GreaterSlantEqual] system[[i, 2]] && num \[LessSlantEqual] system[[i, 3]],
        syst = system[[i, 1]];nn = i,
        None
      ]
   , {i, 1, 7}];
    ltype = Characters[grp1[[3]]][[1]]
  ];
  (*--- Picture of the Bravais lattice ---*)
  If[pc == "Image",
     ln = Length[im[[nn]]]; 
     Do[If[ltype == im[[nn, i, 1]], 
     	   tt = im[[nn, i, 2]], 
     	   None
     	]
     , {i, 1, ln}];
    gc = LatticeData[tt, "Image"], 
    None
  ];
  (*--- Calculate the basis vectors ---*)
  If[pc == "Basis",
  (*--- in the conventional basis the lattice vectors and basis are calculated in cartesian coordinates ---*) 
    If[gbr == "Conventional",
    (*--- Conventional unit cell, non-primitive cell ---*)
       ln = Length[basc[[nn]]];
       Do[
          If[ltype == basc[[nn, i, 1]], 
          	tb = basc[[nn, i, 2]],
            None
          ]
       , {i, 1, ln}];
       bvec = GTLatticeVectors[tb[[1]]];
       lb=Length[tb[[2]]];
       If[ lb== 1,
          (*--- primitive lattice ---*)
          lat = bvec,
          (*--- non-primitive lattice ---*)
          bv = {}; 
          Do[bn = bvec[[1]] tb[[2,i,1]] + bvec[[2]] tb[[2, i, 2]] + bvec[[3]] tb[[2, i, 3]]; 
        	  bv = Append[bv, bn] 
          , {i, 1, lb}];
          lat = {bvec, bv}
       ],
       (*--- Primitive unite cells ---*) 
       ln = Length[basp[[nn]]];
       Do[If[ltype == basp[[nn, i, 1]],
             If[ltype === "P",
                lat = GTLatticeVectors[basp[[nn, i, 2]]],
                lat = basp[[nn, i, 2]]
             ],
             None
          ];
       , {i, 1, ln}];
     ],
     None
  ];
  If[pc == "Image", Show[gc],Return[lat]]
  ]
*)

*)

(*
***)


(****e* /GTLatticeVectors
! NAME
!  GTLatticeVectors
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!  * 01.03.2014 : first version
!  * 15.12.2015 : correction of an error in formula 
!  * 07.03.2023 : new module with definitions corresponding to AFLOW
! USAGE
!  GTLatticeVectors[list] gives the lattice vectors from a list of lattice constants and lattice angles.
! INPUT
!  list  - {a,b,c,alph,beta,gamma}          - old verson
!  list  - {crystal system, Pearson symbol} - new version
! OUTPUT
!  the lattice vectors
!
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  .
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  this was the idea of the previous version:
!  all primitive lattice vectors can be calculated from a general formula, a is chosen in x direction, b lies in the xy plane
!  and z is constructed to enclose the correct angles with a and b. In principle this formula corresponds to the defintion of 
!  the triclinic system.
!  In principle the idea is fine, but the result doesn't correspond in all cases to the standard definitions.
!
!  In the new version all AFLOW definitions are hard coded. The old version is included.
!
! LITERATURE
!  see Crystal info in https://www.aflowlib.org/prototype-encyclopedia/
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

GTLatticeVectors::inp  = "Lattice not correct. Allowed are `1` ";
GTLatticeVectors::imp  = "Not implemented! ";

GTLatticeVectors[lattice_] := Module[
	{a1, b1, c1, \[Alpha]1, \[Beta]1, \[Gamma]1, bt, latt, lat1, tr, tab, cx, cy,vecs, rule,dat},
      Clear[a, b, c, \[Alpha], \[Beta], \[Gamma],a1,b1,c1];
(* Bravais lattices *)     
   latt = {{"Triclinic", "aP"}, {"Monoclinic", "mP"}, {"Monoclinic", 
          "mC"}, {"Orthorhombic", "oP"}, {"Orthorhombic", 
          "oC"}, {"Orthorhombic", "oI"}, {"Orthorhombic", 
          "oF"}, {"Tetragonal", "tP"}, {"Tetragonal", "tI"}, {"Hexagonal", 
          "hP"}, {"Rhombohedral", "hR"}, {"Cubic", "cP"}, {"Cubic", 
          "cF"}, {"Cubic", "cI"}};
(* substitution rule *)          
   rule = {a1 -> a, b1 -> b, c1 -> c, \[Alpha]1 -> \[Alpha], \[Beta]1 -> \[Beta], \[Gamma]1 -> \[Gamma]};
(* old version of command *)
   If[Length[lattice] == 6,
      dat = lattice;
      a1  = dat[[1]]; 
      b1  = dat[[2]]; 
      c1  = dat[[3]]; 
      \[Alpha]1 = dat[[4]];
      \[Beta]1  = dat[[5]]; 
      \[Gamma]1 = dat[[6]];
      bt        = {{a1, 0, 0}, 
      	           {b1 Cos[\[Gamma]1], b1 Sin[\[Gamma]1], 0}, 
      	           {c1 Cos[\[Beta]1], c1 (Cos[\[Alpha]1] - Cos[\[Beta]1] Cos[\[Gamma]1])Sin[\[Gamma]1], 
                    c1 Sqrt[1 - (Cos[\[Alpha]1]^2 + Cos[\[Beta]1]^2 + Cos[\[Gamma]1]^2) + 2 Cos[\[Alpha]1] \
                    Cos[\[Beta]1] Cos[\[Gamma]1]] /Sin[\[Gamma]1]}
                  };
      Return[bt],
      None
   ];   
(* Output information *)
   If[Head[lattice] === List,
      None,
      If[lattice === "Help",
         lat1 = lattice,
         lat1 = "Help"
      ]
   ]; 
   If[lat1 == "Help",
      tr  = latt // Transpose;
      tab = {Table[i, {i, 1, 14}], tr[[1]], tr[[2]]} // Transpose;
      tab = Prepend[tab, {"#","System", "Pearson"}];
      Print[Grid[tab, Frame -> All,Background -> {1 -> LightBlue, 1 -> GTBackGroundColor1, {1, 1} -> GTCornerColor}]]; 
      Return[], 
      None
   ];
(* Wrong input *)
   If[Intersection[{lattice}, latt] == {},
      Message[GTLatticeVectors::inp, latt]; 
      Return[],
      None
   ];
(* Definition of lattices *)
   Switch[lattice,
      {"Triclinic", "aP"},
        cx   = c1 Cos[\[Beta]1]; 
        cy   = c1 (Cos[\[Alpha]1] - Cos[\[Beta]1] Cos[\[Gamma]1])/Sin[\[Gamma]1];
        vecs = {{a1, 0, 0}, {b1 Cos[\[Gamma]1] ,b1 Sin[\[Gamma]1], 0}, {cx, cy, Sqrt[c1^2 - cx^2 - cy^2]}},
      {"Monoclinic", "mP"},
        vecs = {{a1, 0, 0}, {0, b1, 0}, {c1 Cos[\[Beta]1], 0, c1 Sin[\[Beta]1]}},
      {"Monoclinic", "mC"},
        vecs = {{a1/2, -b1/2, 0}, {a1/2, b1/2, 0}, {c1 Cos[\[Beta]1], 0, c1 Sin[\[Beta]1]}},
      {"Orthorhombic", "oP"},
        vecs = {{a1, 0, 0}, {0, b1, 0}, {0, 0, c1}},
      {"Orthorhombic", "oC"},
        vecs = {{a1/2, -b1/2, 0}, {a1/2, b1/2, 0}, {0, 0, c1}},
      {"Orthorhombic", "oI"},
        vecs = {{-a1/2, b1/2, c1/2}, {a1/2, -b1/2, c1/2}, {a1/2, b1/2, -c1/2}},
      {"Orthorhombic", "oF"},
        vecs = {{0, b1/2, c1/2}, {a1/2, 0, c1/2}, {a1/2, b1/2, 0}},
      {"Tetragonal", "tP"},
        vecs = {{a1, 0, 0}, {0, a1, 0}, {0, 0, c1}},
      {"Tetragonal", "tI"},
        vecs = {{-a1/2, a1/2, c1/2}, {a1/2, -a1/2, c1/2}, {a1/2, a1/2, -c1/2}},
      {"Hexagonal", "hP"},
        vecs = {{a1/2, -Sqrt[3] b1/2, 0}, {a1/2, +Sqrt[3] b1/2, 0}, {0, 0, c1}},
      {"Rhombohedral", "hR"},
        vecs = {{a1/2, -a1/(2 Sqrt[3]) , c1/3}, {0, a1/Sqrt[3], c1/3}, {-a1/2, -a1/(2 Sqrt[3]) , c1/3}},
      {"Cubic", "cP"},
        vecs = {{a1, 0, 0}, {0, a1, 0}, {0, 0, a1}},
      {"Cubic", "cF"},
        vecs = {{0, a1/2, a1/2}, {a1/2, 0, a1/2}, {a1/2, a1/2, 0}},
      {"Cubic", "cI"},
        vecs = {{-a1/2, a1/2, a1/2}, {a1/2, -a1/2, a1/2}, {a1/2, a1/2, -a1/2}},
      _,
        Message[GTLatticeVectors::imp]; Return[]
   ];
   (* output *)
   vecs = vecs /. rule;
   Return[vecs]
   ]
 
(*
***) 


(****e* /GTImportCIF
! NAME
!  GTImportCIF
! AUTHOR
!  S. Schenk
! PACKAGE
!   CrystalStructure.m
! MODIFICATION HISTORY
!  * 0.1.10.2015 : first version
!  * 21.10.2016 : update fractional coordinates
!  * 26.01.2016 : option GOCorrectLabels added, tolerance added to find equivalent atoms
!  * 27.06.2018 : chek header and documentation
! USAGE
!  GTImportCIF[filename] opens .cif file and imports relevant data to generate structure data for GTInstallStructure.
! INPUT
!  filename - path to file
! OUTPUT
!  structure data 
!
! GTPack OPTIONS
!  * GOFractional - Switches output format of atom coordinates between cartesian and fractional
!  * GOCorrectLabels
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  .cif is a common data type for crystal structure data, so the package can simply be extended via import of these files
!
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

GTImportCIF[file_,OptionsPattern[]] := Module[{CifData, CifFlags, xyzName, ElemSym, ElemName, SpgName, 
  xyzStrData, xyzExpData, xyzTranslation, Atoms, EveryAtomPos, LatticeVectors, LatticeScales, sol, i, j, k, failed,
  x,y,z, SpgNumber, aa,bb,cc,alpha,beta,gamma,FracOutText},

    If[!FileExistsQ[file], Print["Error:" <> file <> " was not found in the working directory!"]; Abort[] ];
	CifData = Import[file, "string"] <> "\n"; (*--- fix for some strange behaviour ---*)
	CifData = ImportString[CifData, "CIF"];
	CifFlags = Table[Level[CifData[[i]], 1][[1]], {i, Length[CifData]}];
	ElemSym = If[Position[CifFlags, "_chemical_formula_sum"] != {},"_chemical_formula_sum" /. CifData,""];
	ElemName = If[Position[CifFlags, "_chemical_name_mineral"] != {},"_chemical_name_mineral" /. CifData,""];
	
	(*--- get lattice ---*)
	aa = If[Position[CifFlags, "_cell_length_a"] != {},Rationalize["_cell_length_a" /. CifData], Print["lattice error"]; Abort[]];
	bb = If[Position[CifFlags, "_cell_length_b"] != {},Rationalize["_cell_length_b" /. CifData], Print["lattice error"]; Abort[]];
	cc = If[Position[CifFlags, "_cell_length_c"] != {},Rationalize["_cell_length_c" /. CifData], Print["lattice error"]; Abort[]];
	alpha = If[Position[CifFlags, "_cell_angle_alpha"] != {},Rationalize["_cell_angle_alpha" /. CifData] \[Degree], Print["lattice error"]; Abort[]];
	beta = If[Position[CifFlags, "_cell_angle_beta"] != {},Rationalize["_cell_angle_beta" /. CifData] \[Degree], Print["lattice error"]; Abort[]];
	gamma = If[Position[CifFlags, "_cell_angle_gamma"] != {},Rationalize["_cell_angle_gamma" /. CifData] \[Degree], Print["lattice error"]; Abort[]];
	LatticeVectors = GTLatticeVectors[{aa,bb,cc,alpha,beta,gamma}];
	LatticeScales = Table[ToExpression[FromCharacterCode[96+i]] -> Simplify@Norm[LatticeVectors[[i]]] ,{i,Length[LatticeVectors]}];
	LatticeVectors = Table[ToExpression[FromCharacterCode[96+i]]* Simplify@Normalize[LatticeVectors[[i]]] ,{i,Length[LatticeVectors]}];
	
	If[Position[CifFlags, "_atom_site_fract_x"|"_atom_site_fract_y"|"_atom_site_fract_z"] != {},
		Atoms = Table[{
			Rationalize[(
				("_atom_site_fract_x" /. CifData)[[i]]*LatticeVectors[[1]]+
				("_atom_site_fract_y" /. CifData)[[i]]*LatticeVectors[[2]]+
				("_atom_site_fract_z" /. CifData)[[i]]*LatticeVectors[[3]] )/.LatticeScales
			],
			("_atom_site_label" /. CifData)[[i]]
		}, {i,Length[("_atom_site_label" /. CifData)]}];
	,
		If[Position[CifFlags, "_atom_site_Cartn_x"|"_atom_site_Cartn_y"|"_atom_site_Cartn_z"] != {},
			Print["Error: cartesian coords found, which are currently not implemented."];
			Abort[]
		,
			Print["Error: no atom_site_fract entry found!"];
			Abort[]
		];	
	];

	(*--- get multiplicity and transform to cartesian ---*)
	xyzName    = Part[CifFlags, Flatten[Position[StringCount[CifFlags, "xyz"], 1]]];
	xyzStrData = Flatten[xyzName /. CifData];
	xyzExpData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i,Length[xyzStrData]}];
	xyzTranslation = xyzExpData  /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0} ;
	xyzExpData = xyzExpData - xyzTranslation;
	xyzTranslation = Table[ Total[ Table[xyzTranslation[[i,j]]*LatticeVectors[[j]] ,{j,3}] ], {i,Length[xyzTranslation]} ] /.LatticeScales;
	xyzExpData = xyzExpData + xyzTranslation;

	(*--- apply multiplicity ---*)
	For[i = 1, i <= Length[Atoms], i++,
		EveryAtomPos = 
			Sort[
				Union[
					xyzExpData  /. {ToExpression["x"] -> Atoms[[i, 1, 1]], ToExpression["y"] -> Atoms[[i, 1, 2]], ToExpression["z"] -> Atoms[[i, 1, 3]]} 
				],
			Total[#1] > Total[#2] &];
	   
		For[j = 1, j <= Length[EveryAtomPos] - 1, j++,
		For[k = j + 1, k <= Length[EveryAtomPos], k++,
			sol = {x, y, z} /. 
				Flatten[
					Solve[EveryAtomPos[[j]] == (EveryAtomPos[[k]] + x*LatticeVectors[[1]] + y*LatticeVectors[[2]] + z*LatticeVectors[[3]] /. LatticeScales ), {x, y, z}]
				, 1];
			sol = Rationalize[Table[Round[sol[[i]], OptionValue[GOTolerance]], {i, 1, Length[sol]}]];
			If[ IntegerQ[sol[[1]]] && IntegerQ[sol[[2]]] && IntegerQ[sol[[3]]]
				, EveryAtomPos = Delete[EveryAtomPos, k]; k--]
		]];
		Atoms[[i, 1]] = Table[ (EveryAtomPos[[j]]/{a,b,c} /. LatticeScales )*{a,b,c}  ,{j, Length[EveryAtomPos]}];
	 
	(*--- correct the element label ---*)
	    If[OptionValue[GOCorrectLabels],
		sol = Characters[Atoms[[i, 2]]];
		sol = ToUpperCase[sol[[1]]] <> If[Length[sol] > 1 && LetterQ[sol[[2]]], sol[[2]], ""];
		Atoms[[i, 2]] = 
			Quiet[Check[
				ElementData[sol, "Abbreviation"],
				If[! ValueQ[failed], failed = "Possible element name failure, please check manually."];
	    		ElementData[Characters[sol][[1]], "Abbreviation"]
			]];];
	]; (* end for loop *)
	
	Atoms = Flatten[Table[{Atoms[[i, 1, j]], Atoms[[i, 2]]}, {i, 1, Length[Atoms]}, {j, 1, Length[Atoms[[i, 1]]]}], 1];
	If[ValueQ[failed], Print[failed]; Clear[failed]];

	(*--- filter equal length ---*)
	{x,y,z} = {a,b,c} /. LatticeScales;
	If[ PossibleZeroQ[y-z], 
		LatticeScales = Delete[LatticeScales,3]; 
		Atoms = Atoms /.{ c -> b};
		LatticeVectors = LatticeVectors /.{c->b},
	  If[ PossibleZeroQ[x-z],
		LatticeScales = Delete[LatticeScales,3];
		Atoms = Atoms /.{c->a};
		LatticeVectors = LatticeVectors /.{c->a}]
	];
	If[ PossibleZeroQ[x-y],
		LatticeScales = Delete[LatticeScales,2];
		Atoms = Atoms /.{b->a};
		LatticeVectors = LatticeVectors /.{b->a}];
	Clear[x,y,z];
		
	(*--- space group names ---*)
	SpgNumber = If[Position[CifFlags, "_space_group_IT_number"] != {},"_space_group_IT_number" /. CifData,0];
	If[SpgNumber!=0 && IntegerQ[SpgNumber],
		SpgName = GTSpaceGroups[SpgNumber, GOVerbose -> False][[4]]
	,
		If[Position[CifFlags, "_symmetry_space_group_name_H-M"] != {},
			SpgName = "_symmetry_space_group_name_H-M" /. CifData;
			SpgName = StringJoin[DeleteCases[Characters[SpgName], " "]];
			SpgName = StringReplace[SpgName, {"m3m" -> "m-3m"}];
			SpgName = BracketingBar[ToExpression[SpgName]],
			Print["space group name error"];
			SpgName = ""
		]
	];
	
	If[TrueQ[OptionValue[GOFractional]],
		Atoms = Table[{
			Flatten[{x,y,z} /. Solve[Transpose[LatticeVectors].{x,y,z} == Atoms[[i, 1]] /. LatticeScales,{x,y,z}],1]
			, Atoms[[i,2]]}
			, {i, Length[Atoms]}
		];
		FracOutText="fractional";		
	,
		FracOutText="cartesian"
	];

	If[ValueQ[failed], Print[failed]];
	{{ElemSym, ElemName}, "","", SpgName, SpgNumber, LatticeVectors, Atoms, LatticeScales, FracOutText}
]
 
(*
***) 

(****e* /GTExportXSF
! NAME
!  GTExportXSF
! AUTHOR
!  S. Schenk
! PACKAGE
!   CrystalStructure.m
! MODIFICATION HISTORY
!  * 23.06.2016 : first version
!  * 27.06.2018 : check header and documentation 
! USAGE
!  GTExportXSF[filename,data] saves data from e.g. GTGetStructure to filename.
! INPUT
!  filename - path to file
!  data - structure data
! OUTPUT
!  file with structure data
!
! ERROR MESSAGES
!
! GTPack OPTIONS
!  GOVerbose:
!     - True   -  prints the file content in notebook
!     - False  -  no additional information 
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  .xsf is a simple data type, which can be used to visualise the structure in XCrySDen
!
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

 GTExportXSF[filename_, data_, OptionsPattern[]] := 
 Module[{file, atoms, atomstr},
  file = OpenWrite[filename];
  WriteString[file, "CRYSTAL\n"];
  WriteString[file, "PRIMVEC\n"];
  WriteString[file, StringJoin[Table[StringJoin[{Table[
        ToString[
         PaddedForm[N[data[[6, j, i]] /. data[[8]]], {14, 10}]], {i, 
         1, 3}], "\n"}], {j, 1, 3}]]];
  (*WriteString[file,"CONVVEC\n"];
  WriteString[file,"ToDo\n"];*)
  WriteString[file, "PRIMCOORD\n"];
  atoms = 
   Table[
     Join[{data[[7, i, 2]]}, Table[ data[[7, i, 1, j]] /. data[[8]] , {j, Length[data[[7, i, 1]]]}] ] , {i,Length[data[[7]]]}];
  WriteString[file, StringJoin[{ToString[Length[atoms]], " 1\n"}]];
  atomstr = 
   StringJoin[
    Table[ToString[ElementData[atoms[[i, 1]], "AtomicNumber"]] <> 
      ToString[PaddedForm[N[atoms[[i, 2]]], {14, 10}]] <> 
      ToString[PaddedForm[N[atoms[[i, 3]]], {14, 10}]] <> 
      ToString[PaddedForm[N[atoms[[i, 4]]], {14, 10}]] <> "\n", {i, 
      Length[atoms]}]];
  WriteString[file, atomstr];
  Close[file];
  If[TrueQ[OptionValue[GOVerbose]],
   FilePrint[filename]];
  ]
  
(*
***)


(****e* /GTPoly
! NAME
!  GTPoly
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!   * 15.08.2015  : first version
!   * 05.04.2016  : from Test.m to CrystalStructure.m
!   * 27.06.2016  :  check header
! USAGE
!  constructs objects for graphical representation of rotational axes
! INPUT
!  * p     - order of axes
!  * scale - scaling factor
!  * vec   - direction of rotational axes
! OUTPUT
!  a graphical object
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  - 
! GTPack MODULES
! -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  the roational axis will be represented by a small cylinder in the direction of the axes
!  and a handle at both ends to symbolize the order of the axes:
!  2 -small bar, 3- triangle, 4 -square, 6 -hexagon
!
!  It is an internal module.  
!
! LITERATURE
!  
! TODO
!  not fully tested
!  
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTPoly[p_, scale_, vec_] := Module[{k, object, vo, mat, pg, el, v1},
  If[Cross[{0, 0, 1}, vec] == {0, 0, 0},
     mat = IdentityMatrix[3],
     mat = RotationMatrix[{{0, 0, 1}, vec}]
  ];
  vo = mat.{0, 0, 1}; v1 = mat.{.1, 0, 0};
  If[p > 2,
     pg = Table[scale {Cos[2 Pi k /p], Sin[2 Pi k /p], 0}, {k, 1, p}];
     pg = Map[mat.# &, pg];
     pg = Polygon[pg];
     object = {Translate[pg, vo], Translate[pg, -vo], Cylinder[{vo, -vo}, .01]},
     el = Cylinder[{-v1, v1}, 0.01];
     object = {Translate[el, -vo], Translate[el, vo], Cylinder[{-vo, vo}, .01]}
  ];
  Return[object]
]

(*
***)


(****e* /GTPlane
! NAME
!  GTPlane
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!   * 15.08.2015  : first version
!   * 10.04.2016  : from Test.m to CrystalStructure.m
!   * 27.06.2016  : check header
! USAGE
!  Constructs reflection planes for a graphical representation of symmetry elements.
! INPUT
!  sym - symmetry element
! OUTPUT
!  a square, representing the reflection plane
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  - 
! GTPack MODULES
!  GTGetQuaternion
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The vector part of the quaternion representation of the symmetry element gives the normal of the
!  reflection plane. The normal is used to rotate a square, lying in the x-y-plane in the correct position.
!
! It is an internal module.
!
! LITERATURE
!  
! TODO
!  not fully tested
!  
! PROBLEMS
!  may be not really a problem: Due to the rotations the position of the squares are in principle correct, but sometimes
!  rotated.
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTPlane[sym_] := Module[{pg, v, m, object}, 
  pg = .5 {{-1, -1, 0}, {-1, 1, 0}, {1, 1, 0}, {1, -1, 0}};
  v = GTGetQuaternion[sym][[2]]; 
  If[Cross[{0, 0, 1}, v] == {0, 0, 0},
      m = IdentityMatrix[3],
      m = RotationMatrix[{{0, 0, 1}, v}]
  ]; 
  object = {Opacity[0.4], Polygon[Map[m.# &, pg]]};
  Return[object]
]

(*
***)


(****e* /GTDelInv
! NAME
!  GTDelInv
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!   * 10.08.2015 : first version
!   * 08.042016  : from Test.m to CrystalStructure.m
!   * 27.06.2018 : check header
! USAGE
!  deletes the inverse elements in a list of elements
! INPUT
!  elm - list of symmetry elements
! OUTPUT
!  list of elements without inverse elements
!
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTInverseElement
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  For the graphical representation of symmetry elements only one element of the rotational axis is necessary.
!  Therefore, i.e. C3z and C3zi are equivalent in the sense that both represent the same rotational axes.
!
!  This is an internal module.
!
! LITERATURE
!  -
! TODO
!  not fully tested
! RELEASE
!  1.0.0 
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTDelInv[elm_] := Module[{l, t, i},
   l = {}; t = Sort[elm]; 
   Do[
   	  l = Append[l, t[[1]]]; 
      t = Complement[t, {t[[1]], GTInverseElement[t[[1]]]}]
   , {i, 1, Length[elm]/2}]; 
   Return[l]
]


(*
***)


(****e* /GTSelectElements
! NAME
!  GTSelectElements
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!   * 10.08.2015  : first version
!   * 20.04.2016  : from Test.m to CrystalStructure.m
!   * 27.06.2018  : check header
! USAGE
!  Selects the minimum number of elements necessary for the representation of the symmetry elements
! INPUT
!  list of symmetry elements (proper rotations)
! OUTPUT
!  list of lists: {r2,r3,r4,r6}, containing the symmetry elements to demonstrate 2fold (r2), 3fold (r3)
!  4fold (r4) and 6fold (r6) rotations
!
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTOrderOfElement, GTDelInv, GTDelRots
! GTPack NOTEBOOKS 
!  -  
! DESCRIPTION
!  All inverse elements are cancelled.
!  If an rotational axes exists in a certain direction, it is necessary to keep only the rotation of highest order.
!  The module selects the minimum number of rotations. If, for exapmle, the z-axis is a 6fold rotational axis, 3fold and 
!  2fold rotations with this axis are cancelled.
! LITERATURE
!  -
! TODO
!  not fully tested
! RELEASE
!  1.0.0 
! PROBLEMS
! -
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSelectElements[elm_] := Module[{r2, r3, r4, r6, r21, r31},
  r2 = Select[elm, GTOrderOfElement[#] == 2 &];
  r3 = Select[elm, GTOrderOfElement[#] == 3 &]; r3 = GTDelInv[r3];
  r4 = Select[elm, GTOrderOfElement[#] == 4 &]; r4 = GTDelInv[r4];
  r6 = Select[elm, GTOrderOfElement[#] == 6 &]; r6 = GTDelInv[r6];
  (*--- Select elements for graphics ---*)
  r21 = GTDelRots[r2, r4];
  r21 = GTDelRots[r21, r6];
  r31 = GTDelRots[r3, r6];
  Return[{r21, r31, r4, r6}]
  ]

(*
***)


(****e* /GTDelRots
! NAME
!  GTDelRots
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrytalStructure.m
! MODIFICATION HISTORY
!   * 13.08.2015  : first version
!   * 20.04.2016  : from Test.m to CrystalStructure.m
!   * 27.06.2018  : check header
! USAGE
!  Delets rotations of lower order
! INPUT
!  * elm1 - list of elements (lower order)
!  * elm2 - list of elements (higher order)
! OUTPUT
!  minimum number of rotations
!
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTGetQuaternion
! GTPack NOTEBOOKS 
!  -  
! DESCRIPTION
!  If an rotational axis exists in a certain direction, it is necessary to keep only the rotation of highest order.
!  The module selects the minimum number of rotations. If, for example, elm1 contains 2fold rotations and elm2 contains 6fold rotations,
!  all 2fold rotations in directions of 6fold rotations are cancelled.
! LITERATURE
!  
! TODO
!  not fully tested
!  
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTDelRots[elm1_, elm2_] := Module[{t1, t2, l1, l2, v1, v2, arg, arg1, i, k},
   If[elm2 == {},
      Return[elm1],
      t1 = elm1; l1 = Length[t1];
      t2 = elm2; l2 = Length[t2];
      Do[
         v2 = GTGetQuaternion[t2[[k]]][[2]]; v2 = v2/Norm[v2];
         Do[
         	arg = t1[[i]];
            v1 = GTGetQuaternion[arg][[2]]; v1 = v1/Norm[v1];
            If[v1 == v2,
               arg1 = arg,
               None
            ]
         , {i, 1, l1}];
         t1 = Complement[t1, {arg}]; l1 = Length[t1];
      , {k, 1, l2}];
      Return[Union[t1]]
   ]
]

(*
***)


(****e* /GTShowSymmetryElements
! NAME
!  GTShowSymmetryElements
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!  * 13.08.2015 : first version
!  * 14.04.2016 : from Test.m to CrystalStructure.m
!  * 27.06.2018 : check header and documentation  
! USAGE
!  Graphical representation of point group symmetries
! INPUT
!  group - one of the 32 point groups 
! OUTPUT
!  visualization of the poing group symmetries
!  GTPack OPTIONS
!  GOSelectSymmetry 
!                   - "All"             - all symmetry elements taken nto account (standard)
!                   - "Reflections"     - only reflections taen into account
!                   - "Rotations"       - only rotations, but all taken into account
!                   - {"Rotations",{4}} - only 4fold rotations are shown   
! STANDARD OPTIONS
!   -
! GTPack MODULES
!  GTGetMatrix, GTSelectElements, GTGetQuaternion, GTPoly, GTOrderOfElement, GTPlane
! GTPack NOTEBOOKS 
!  * Reference_Symbols : GTShowSymmetryElements.nb 
!  * Wolfram_Book_Boxes : GTShowSymmetryElements.nb
! DESCRIPTION
!  The symmetry elements of a point group are visualized. A yellow sphere at the origin is used to 
!  demonstrate the inversion. Rotation axes are presented by a small cylinder in the coresponding 
!  direction. Handels at both ends demonstrate the order (2,3,4,6). Rotation planes are given by squares.
!  Large groups contain a lot of symmetry elements. The options allows to reduce the amout of presented
!  symmetry elements. 
! LITERATURE
!  
! TODO
!  not fully tested
!  improper symmetry elements have to be checked.
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTShowSymmetryElements[group_, OptionsPattern[]] := 
	 Module[{grp,ge,rots,ss,st,tre,trr,inv,prop,improp,r2,r3,r4,r6,rot,ref,i},
     grp = group; ge = {}; rots = {2, 3, 4, 5, 6};
 (*--- Interpretation of option ---*)
     ss = OptionValue[GOSelectSymmetry];
     If[Length[ss] == 2,
        st = ss[[1]]; rots = ss[[2]],
        st = ss
     ];
     If[st == "All",
        tre = True; trr = True,
        tre = False; trr = False;
        If[ss === "Reflections",
           tre = True,
           tre = False
        ];
        If[st === "Rotations",
           trr = True,
           trr = False
        ]
     ];
(*--- Test for inversion ---*)
     inv = {Specularity[White, 30], Orange, Sphere[{0, 0, 0}, .1]};
     If[Intersection[grp, {IEe}] == {},
        None,
        ge = Append[ge, inv]; grp = Complement[group, {IEe}]
     ];
(*--- Proper and improper operations ---*)
     prop     = Select[grp, Det[GTGetMatrix[#]] == 1 &];
     improp = Complement[grp, prop];
(*--- Select elements for graphics ---*)
     {r2, r3, r4, r6} = GTSelectElements[prop];
(*--- Investigation proper rotations ---*)
     If[trr,
(*--- 3fold rotations ---*)
        If[r3 == {},
           None,
           If[Intersection[rots, {3}] == {},
              None,
              rot = Transpose[Map[GTGetQuaternion[GTGetMatrix[#]] &, r3]][[2]];
              ge = Append[ge, Map[GTPoly[3, .1, #] &, rot]]
           ]
        ];
(*--- 4fold rotations ---*)
        If[r4 == {},
           None,
           If[Intersection[rots, {4}] == {},
              None,   
              rot = Transpose[Map[GTGetQuaternion[GTGetMatrix[#]] &, r4]][[2]];
              ge = Append[ge, Map[GTPoly[4, .1, #] &, rot]]
           ]
        ];
(*--- 2fold rotations ---*)
        If[r2 == {},
           None,
           If[Intersection[rots, {2}] == {},
              None,  
              rot = Transpose[Map[GTGetQuaternion[GTGetMatrix[#]] &, r2]][[2]];
              ge = Append[ge, Map[GTPoly[2, .1, #] &, rot]]
           ]
        ];
 (*--- 6fold rotations ---*)
        If[r6 == {},
           None,
           If[Intersection[rots, {6}] == {},
              None,
              rot = Transpose[Map[GTGetQuaternion[GTGetMatrix[#]] &, r6]][[2]];
              ge = Append[ge, Map[GTPoly[6, .1, #] &, rot]]
           ]
        ],
        None
     ];
(*--- Reflections ---*)
     If[tre,
        ref = Select[improp, GTOrderOfElement[#] == 2 &];
        Do[
           ge = Append[ge, GTPlane[ref[[i]]]]
        , {i, 1, Length[ref]}],
        None
     ];
     Graphics3D[ge, PlotRange -> {{-1.2, 1.2}, {-1, 1.2}, {-1, 1.2}}]
]

(*
***)


(****e* /GTPlotCluster
! NAME
!  GTPlotCluster
! AUTHOR
!  W. Hergert
! PACKAGE
!   CrystalStructure.m
! MODIFICATION HISTORY
!  * 23.08.2015 : first version
!  * 11.11.2015 : GTPlotCluster and GTPlotStructure2D unified
!  * 04.04.2016 : from Test.m to CrystalStructure.m
!  * 09.11.2021 : coloring of bonds implemented, documentation extended accordinglingly
! USAGE
!  Plot of manipulated clusters as a 2D or 3D plot
! INPUT
!  * cluster - a cluster of atoms as a result of GTCluster or GTClusterManipulate
!  * dist    - maximum length for the bonds
!  * scale   _ scaling factor for the size of the atoms
! OUTPUT
!  plot of the cluster (2D or 3D) or a graphics object
! GTPack OPTIONS
!  GOColorScheme
!	o "ElementData" - The color scheme for the atoms stored in ElementData will be used. (Standard)
!	o List of colors for the basis atoms.
!  GOPlot 
!	o True  - output is a plot of the cluster (Standard)
!	o False - output is a graphics object   
!  GOSphere
!	o True  - a circumscribing sphere around the cluster is plotted (Standard)
!	o False - no sphere 
!  GOCoordinateSystem
!	o No coordinate system will be draw (standard)
!	o a coordinate system will be given
!  GODimension
!	o =2 2D cluster
!	o =3 3D cluster (standard)
!  GODirection
!   o {0,0,0}  nothig additional plotted 
!   o {a,b,c}  defines a direction that will be indicated by an arrow 
!  GOBondColor
!	o {}                    standard setting
!   o {{at1,at2,},{col}}    define color for special bonds, example : {{{"O", "O"}, {Thick, Red}}}
!  
! STANDARD OPTIONS
!   - 
! GTPack MODULES
!   -
! GTPack NOTEBOOKS 
!  Reference_Symbols : GTPlotCluster.nb, 
! DESCRIPTION
!  A simple programm to get an overview on a cluster of atoms generated by GTCluster and
!  eventually manipulated by GTClusterManipulate. If ClusterManipulate is used, a tool is
!  helpful that allows also to control the manipulations done.
! LITERATURE
!  -
! TODO
!  not fully tested
! RELEASE
!  1.0.0 
! PROBLEMS
!  The new version includes also the capabilities of the command GTPlotStructure2D.
!  The new version was buildt because it was not optimal that a pure 2D cluster was
!  constructed as a 3D graphics object in GTPlotCluster.
!--------------------------------------------------------------------------------
! SOURCE
*)
GTPlotCluster::color =  "Number of atoms in basis and number of colors do not correspond";


GTPlotCluster[cluster_, dist_, scale_, OptionsPattern[]] := Module[
	{atoms, basisat, opc, oplt, osp, rule, aml, i, j, lc, glines, nom, gsph, no, col, cosys, la,
     codir, odim, coord,gstruc,gcos,box,bonds,bcl,p,bcols,bond},
           atoms   = Transpose[cluster][[2]]; 
           lc      = Length[cluster]; 
           basisat = Union[atoms];
(*--- options ---*)
           opc   = OptionValue[GOColorScheme];
           oplt  = OptionValue[GOPlot];
           osp   = OptionValue[GOSphere];
           cosys = OptionValue[GOCoordinateSystem];
           codir = OptionValue[GODirection];
           odim  = OptionValue[GODimension];
           bcl   = OptionValue[GOBondColor];
           box   = OptionValue[Boxed];
           bcl   = OptionValue[GOBondColor];
           coord = Transpose[cluster][[1]];
(*--- prepare bond colors ---*)
           If[bcl=={},
           	  None,
           	  {bonds, bcols} = Transpose[bcl];
           	  bonds=Sort[#] & /@ bonds
           ];	  
(*--- coordinates 2D ---*)
           If[odim == 2,
              coord = Take[Transpose[coord], {1, 2}] // Transpose,
              None
           ];
(*--- Color from Element Data ---*)
           If[opc === "ElementData",
              aml = {};
              Do[
                 aml = Append[aml, ElementData[cluster[[i, 2]], "IconColor"]]
              , {i, 1, lc}];
              Print[TableForm[{basisat, Map[ElementData[#, "IconColor"] &, basisat]} //Transpose]],
(*--- Color from Input ----*)
              If[Length[basisat] == Length[opc],
                 None,
                 Message[GTPlotCluster::color];Abort[]
              ];
              rule = {};
              Do[
                 rule = Append[rule, basisat[[i]] -> opc[[i]]]
              , {i, 1, Length[basisat]}];
              aml = Transpose[cluster][[2]] /. rule;
              Print[TableForm[{basisat, opc} // Transpose]]
           ];
(*--- generate bonds and spheres/discs/ coloring of bonds ---*)
           lc     = Length[cluster];
           glines = {}; 
           gsph   = {};
           nom    = Map[Norm[#] &, coord]; 
           nom    = Max[nom];
           la = 1.4*nom;
           Do[
           	  Do[
           	  	 no = Norm[cluster[[i, 1]] - cluster[[j, 1]]];
                 If[no <= dist, 
                 	If[bcl=={},
                       glines = Append[glines, Line[{coord[[i]], coord[[j]]}]], 
                       bond=Sort[{cluster[[i, 2]],cluster[[j, 2]]}];                    
                       p = Flatten[Position[bonds, bond]];   
                       If[p=={},
                          glines = Append[glines,{Gray, Line[{coord[[i]], coord[[j]]}]}],                
                          glines = Append[glines,Flatten[{bcols[[p[[1]]]], Line[{coord[[i]], coord[[j]]}]}]]
                       ]
                 	],
                    None
                 ]
              , {j, i + 1, lc}];
              If[opc === "ElementData", 
                 col = ElementData[cluster[[i, 2]], "IconColor"], 
                 col = aml[[i]]
              ];
              If[odim == 2,
                 gsph = Append[gsph, {aml[[i]], Disk[coord[[i]], 0.03*scale]}],
                 gsph = Append[gsph, {aml[[i]], Sphere[coord[[i]], 0.03*scale]}]
              ]
           , {i, 1, lc}];
(*--- the full structure ---*)
           gstruc = {gsph, glines};
           If[osp,
              If[odim == 3,
                 gstruc = Append[gstruc, {Opacity[0.4], Sphere[{0, 0, 0}, nom]}],
                 gstruc = Append[gstruc, {Opacity[0.4], Disk[{0, 0}, nom]}]
              ]
           ];
(*--- coordinate system ---*)
           If[odim == 3,
              If[cosys,
                 gcos = {{Red,   Arrowheads[0.05], Arrow[Tube[{{0, 0, 0}, {la, 0, 0}}]]}, 
                 	     {Green, Arrowheads[0.05], Arrow[Tube[{{0, 0, 0}, {0, la, 0}}]]}, 
                 	     {Blue,  Arrowheads[0.05], Arrow[Tube[{{0, 0, 0}, {0, 0, la}}]]}
                 	    };
                 gstruc = Append[gstruc, gcos];
                 Print["x-axis: red, y-axis: green, z-axis: blue"], 
                 None
              ], 
              None
           ];
           If[odim == 2,
              If[cosys,
                 gcos = {{Red,   Arrow[{{0, 0}, {la, 0}}]}, 
                 	     {Green, Arrow[{{0, 0}, {0, la}}]}
                 	    };
                 gstruc = Append[gstruc, gcos];
                 Print["x-axis: red, y-axis: green"], 
                 None
              ]
           ];
(*--- plot direction ---*)
           If[codir == {0, 0, 0}, 
           	  None, 
           	  codir = codir*la*nom;
              If[odim == 2,
                 gstruc = Append[gstruc, {Magenta, Arrow[{{0, 0}, codir}]}],
                 gstruc = Append[gstruc, {Magenta, Arrowheads[0.05],Arrow[Tube[{{0, 0, 0}, codir}]]}]
              ]
           ];
(*--- export either graph or list ---*)
           If[oplt,
              If[odim == 2, 
              	 Return[Graphics[gstruc]], 
                 Return[Graphics3D[gstruc,Boxed-> box]]
              ], 
              Return[gstruc]
            ]
]

(* old version of GTPlotCluster without coloring of bonds 

GTPlotCluster[cluster_, dist_, scale_, OptionsPattern[]] := 
   Module[{atoms, basisat, opc, oplt, osp, rule, aml, i, j, lc, glines, nom, gsph, no, col, cosys, la,
           codir, odim, coord,gstruc,gcos,box},
           atoms   = Transpose[cluster][[2]]; lc = Length[cluster]; 
           basisat = Union[atoms];
(*--- options ---*)
           opc   = OptionValue[GOColorScheme];
           oplt  = OptionValue[GOPlot];
           osp   = OptionValue[GOSphere];
           cosys = OptionValue[GOCoordinateSystem];
           codir = OptionValue[GODirection];
           odim  = OptionValue[GODimension];
           box   = OptionValue[Boxed];
           coord = Transpose[cluster][[1]];
(*--- coordinates 2D ---*)
           If[odim == 2,
              coord = Take[Transpose[coord], {1, 2}] // Transpose,
              None
           ];
(*--- Color from Element Data ---*)
           If[opc === "ElementData",
              aml = {};
              Do[
                 aml = Append[aml, ElementData[cluster[[i, 2]], "IconColor"]]
              , {i, 1, lc}];
              Print[TableForm[{basisat, Map[ElementData[#, "IconColor"] &, basisat]} //Transpose]],
(*--- Color from Input ----*)
              If[Length[basisat] == Length[opc],
                 None,
                 Print["Error: Number of atoms in basis and number of colors do not correspond"]; Abort[]
              ];
              rule = {};
              Do[
                 rule = Append[rule, basisat[[i]] -> opc[[i]]]
              , {i, 1, Length[basisat]}];
              aml = Transpose[cluster][[2]] /. rule;
              Print[TableForm[{basisat, opc} // Transpose]]
           ];
(*--- generate bonds and spheres/discs ---*)
           lc = Length[cluster];
           glines = {}; gsph = {};
           nom = Map[Norm[#] &, coord]; nom = Max[nom];
           la = 1.4*nom;
           Do[
           	  Do[
           	  	 no = Norm[cluster[[i, 1]] - cluster[[j, 1]]];
                 If[no <= dist, 
                    glines = Append[glines, Line[{coord[[i]], coord[[j]]}]], 
                    None
                 ]
              , {j, i + 1, lc}];
              If[opc === "ElementData", 
                 col = ElementData[cluster[[i, 2]], "IconColor"], 
                 col = aml[[i]]
              ];
              If[odim == 2,
                 gsph = Append[gsph, {aml[[i]], Disk[coord[[i]], 0.03*scale]}],
                 gsph = Append[gsph, {aml[[i]], Sphere[coord[[i]], 0.03*scale]}]
              ]
           , {i, 1, lc}];
(*--- the full structure ---*)
           gstruc = {gsph, glines};
           If[osp,
              If[odim == 3,
                 gstruc = Append[gstruc, {Opacity[0.4], Sphere[{0, 0, 0}, nom]}],
                 gstruc = Append[gstruc, {Opacity[0.4], Disk[{0, 0}, nom]}]
              ]
           ];
(*--- coordinate system ---*)
           If[odim == 3,
              If[cosys,
                 gcos = {{Red,   Arrowheads[0.05], Arrow[Tube[{{0, 0, 0}, {la, 0, 0}}]]}, 
                 	     {Green, Arrowheads[0.05], Arrow[Tube[{{0, 0, 0}, {0, la, 0}}]]}, 
                 	     {Blue,  Arrowheads[0.05], Arrow[Tube[{{0, 0, 0}, {0, 0, la}}]]}
                 	    };
                 gstruc = Append[gstruc, gcos];
                 Print["x-axis: red, y-axis: green, z-axis: blue"], 
                 None
              ], 
              None
           ];
           If[odim == 2,
              If[cosys,
                 gcos = {{Red,   Arrow[{{0, 0}, {la, 0}}]}, 
                 	     {Green, Arrow[{{0, 0}, {0, la}}]}
                 	    };
                 gstruc = Append[gstruc, gcos];
                 Print["x-axis: red, y-axis: green"], 
                 None
              ]
           ];
(*--- plot direction ---*)
           If[codir == {0, 0, 0}, 
           	  None, 
           	  codir = codir*la*nom;
              If[odim == 2,
                 gstruc = Append[gstruc, {Magenta, Arrow[{{0, 0}, codir}]}],
                 gstruc = Append[gstruc, {Magenta, Arrowheads[0.05],Arrow[Tube[{{0, 0, 0}, codir}]]}]
              ]
           ];
(*--- export either graph or list ---*)
           If[oplt,
              If[odim == 2, 
              	 Return[Graphics[gstruc]], 
                 Return[Graphics3D[gstruc,Boxed-> box]]
              ], 
              Return[gstruc]
            ]
]

*)


(*
***)
  
  
 
(****e* /GTClusterManipulate
! NAME
!  GTClusterManipulate
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!  * 15.10.2015 : first version
!  * 04.04.2016 : from Test.m to CrystalStructure.m
!  * 13.07.2016 : application of symmetry element included
!  * 27.06.2018 : check header and documentation
!  * 19.09.2018 : Method Impurities implemented
! USAGE
!  Manipulate a given cluster 
! INPUT
!  cluster -  a cluster of atoms as a result of GTCluster or GTClusterManipulate
! OUTPUT
!  the manipulated cluster
! GTPack OPTIONS
!  GOMethod
!  The option GOMethod can be a keyword (String) only or a list consisting of a keyword and an argument.
!   Keyword only:
!	o "Help"      - lists all possibilities of GOMethod
!	o "Surface"   - keeps atoms with z-coordinates z<=0
!	o "Project"   - projects atoms in the plane z=0
!	Keyword and argument:
!   o {"Translate",vector}         - translates the whole cluster with vector 
!	o {"AddAtoms",listofatoms}     - append the list of atoms to the cluster
!   o SymmetryElement",element}	   - applies symmetry element
!	o {"RemoveAtoms",listofatoms}  - remove the atoms in the list from the cluster
!	o {"SelectType",type}          - select atoms of a given type from the cluster
!   o {"RemoveType",type}          - remove atoms of a given type from the cluster
!   o {"CutPlane",{normal,origin}} - cuts out of the cluster aplane define by the normal of the plan and a point in the plane
!	o {"Slab",number"}             - the first number layers near the surface are cut out
!	o {"Rotate",rotation}          - the whole cluster will be rotated , the rotation can be described in different ways.
!									  theta, w -> counterclockwise rotation with theta around vector w
!									  {u,w} -> gives the matrix, that rotatesthe vector u to vector w
!   o {"CutSphere",{position, radius}} 
!                                   - a sphere of radius will be cut out of the thecluster around position.
!   o {"Impuities,{{sort1,sort2},rad,conc} 
!                                   - in a sphere of radius rad atoms of sort1 will be substtuded by atoms of sort2 randomly.
!                                     The concentration is conc.
!                        
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  Reference_Symbols : GTClusterManipulate.nb, GOMethod.nb
! DESCRIPTION
!  Sometimes it is necessary to modify the cluster generated by means of GTCluster, to represent finally
!  special crystal structures. In this modul a series of methods to manipulate the cluster are implemented. 
!  New methods can be easily implemented.
! LITERATURE
!  -
! TODO
!  not fully tested, especially rotation
! RELEASE
!  1.0.0 
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)

GTClusterManipulate[cluster_, OptionsPattern[]] := Module[
    {opt, task, arg, coord, atoms, cll, lc, i, pos, opts,test,normal,origin,cp,zco,mat,atom,subst,rc,conc,pos1,p,nc,rpos}, 
     opts= {{"Translate", "Translates the whole cluster"},
            {"AddAtoms", "Adds additional atoms to the cluster"},
            {"RemoveAtoms", "Removes atoms from th cluster"},
            {"SelectType", "Selects atoms of a certain type"},
            {"RemoveType", "Removes atoms of a certain type from the cluster"},
            {"Rotate", "Rotates the whole cluster"},
            {"Slab", "Cuts a slab out of the cluster"},
            {"Surface", "Generates a surface"},      
            {"CutPlane", "Cuts a plane out of a cluster"},
            {"SymmetryElement", "Applies pointgroup element"},
            {"CutSphere", "Cuts a sphere out of cluster around a given position"},
            {"Project", "Projects the atoms in the cluster to plane z=0"},
            {"Impurities","Substitutes atoms in Cluster randomly"}
           };
     coord = Transpose[cluster][[1]]; atoms = Transpose[cluster][[2]]; lc = Length[cluster];cll = {}; 
     (*--- Evaluate options ----*)
  opt = OptionValue[GOMethod];
  If[Head[opt] === List,
     task = OptionValue[GOMethod][[1]];
     arg = OptionValue[GOMethod][[2]];
     test = Intersection[Transpose[opts][[1]], {task}];
     If[test == {},
        Print["Error: Option ", task, " not implemented!"]; Abort[],
        None
     ],
     If[opt == "Help", 
        Print["Options are:"]; 
        Print[TableForm[opts]];Return[],
        task = opt
     ]
  ];  
  If[task == "CutPlane",
     normal = arg[[1]]; origin = arg[[2]];
     Do[cp = cluster[[i, 1]];
        If[(cp - origin).normal == 0,
           cll = Append[cll, cluster[[i]]],
           None
        ] 
     , {i, 1, lc}],
     None
  ]; 
  If[task == "Surface",
     Do[
        If[cluster[[i, 1, 3]] <= 0,
           cll = Append[cll, cluster[[i]]],
           None
        ] 
      , {i, 1, lc}],
      None
   ];
  If[task == "Translate", 
     coord = Map[arg + # &, coord]; 
     cll = Transpose[{coord, atoms}],
     None
  ];
  If[task == "AddAtoms",
     cll = Join[cluster, arg],
     None
  ];
  If[task == "RemoveAtoms",
     cll = Complement[cluster, arg],
     None
  ];
  If[task == "SelectType",
     pos = Position[Transpose[cluster][[2]], arg] // Flatten;
     Do[
     	cll = Append[cll, cluster[[pos[[i]]]]]
     , {i, 1, Length[pos]}],
     None
  ];
  If[task == "RemoveType",
     pos = Position[Transpose[cluster][[2]], arg] // Flatten;
     pos = Partition[pos, 1];
     cll = Delete[cluster, pos],
     None
  ];
  If[task == "Slab",
     pos = Sort[Transpose[Transpose[cluster][[1]]][[3]] // Union, Greater];
     zco = pos[[arg]]; 
     Do[
        If[cluster[[i, 1, 3]] >= zco,
           cll = Append[cll, cluster[[i]]],
           None
        ]
     , {i, 1, lc}],
     None
  ];
  If[task == "Rotate",
     If[Head[arg[[1]]] === List,
        mat = RotationMatrix[arg], 
        mat = RotationMatrix[arg[[1]], arg[[2]]]
     ];
     coord = Map[mat.# &, coord]; 
     cll = Transpose[{coord, atoms}],
     None
  ];
  If[task == "Project",  
     Do[
        coord[[i, 3]] = 0
     , {i, 1, lc}]; cll = Transpose[{coord, atoms}],
     None
  ];
  If[task == "CutSphere", 
     cll = Select[cluster, Norm[#[[1]] - arg[[1]]] <= arg[[2]] &],
     None
  ];
  If[task == "SymmetryElement",
     mat = GTGetMatrix[arg]//N;
     Do[
     	coord[[i]] = mat.coord[[i]]
     , {i, 1, lc}]; 
     cll = Transpose[{coord, atoms}]
     , None
  ];
  If[task == "Impurities",
     atom  = arg[[1, 1]];
     subst = arg[[1, 2]];
     rc    = arg[[2]];
     conc  = arg[[3]];
     pos   = Position[atoms, atom] // Flatten;
     pos1 = {};
     Do[
        p = pos[[i]];
        If[Norm[coord[[p]]] < rc,
           pos1 = Append[pos1, p],
           None
        ]
     , {i, 1, Length[pos]}];
     nc    = Round[conc*Length[pos1]];
     Print[
     Grid[
     	{{"Number of atoms", atom, Length[pos]}, 
     	 {"Number selected positions", atom, Length[pos1]}, 
     	 {"Number substituted atoms", subst, nc}, 
     	 {"actual concentration (%)", subst, nc/Length[pos1]*100.}}, 
     	 Frame -> All,Background->{{1->GTBackGroundColor1,2->GTCornerColor},None}
     	 ]
     ];
     rpos = RandomSample[pos1, nc];
     Do[
        atoms[[rpos[[i]]]] = subst
     , {i, 1, nc}];
     cll = {coord, atoms} // Transpose,
     None
  ];
  Return[cll]
]

(*
***)

(****e* /GTBuckyBall
! NAME
!  GTBuckyBall
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!   * 10.10.2015 : first version
!   * 18.04.2016 : from Test.m to CrystalStructure.m
!   * 18.06.2016 : correction of bugs, rotation of the cluster such, that it fits to the definition of the
!                  icosahedron group, bond length  is now much simpler
! USAGE
!  Constructs a C60 molecule (buckyball)
! INPUT
!  o r5 - bond lengths in pentagons
!  o r6 - lengths of bonds connecting pentagons
! OUTPUT
!  cluster, describing buckminsterfullerene
! GTPAck OPTIONS
!  GOColorScheme
!	o "ElementData" - The color scheme for the atoms stored in ElementData will be used. (Standard)
!	o List of colors for the basis atoms.
!  GOPlot 
!	o True  - output is a plot of the cluster (Standard)
!	o False - output is a graphics object   
!  GOSphere
!	o True  - a circumscribing sphere around the cluster is plotted (Standard)
!	o False - no sphere 
!  GOCoordinateSystem
!	o No coordinate system will be draw (standard)
!	o a coordinate system will be given
!  GODirection
!   o {0,0,0}  nothig additional plotted 
!   o {a,b,c}  defines a direction that will be indicated by an arrow 
!  
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  Reference_Symbols : GTBuckyBall.nb
! DESCRIPTION
!  The coordinates of the buckyball are generated. The output is 
!  a list in the form also generated by GTCluster. Therefore this 
!  form fits well for further use in GTPack.
!
!  The method of Peter Senn is used to generate the coordinates.
! LITERATURE
!  o Peter Senn, Comutation of the Cartesian Coordinates of Buckminsterfullere,
!    Journal of Chemical Education, 72, 302 (1995)
!  o W.F. David et al., Crystal structure and bonding of ordered C_60, Nature, 353, 147 (1991).
! TODO
!  It is not totally clear if this is the best method to generate the cluster. It could be also
!  good to implement the ideal truncated isosaedron.
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)

GTBuckyBall[r5_, r6_, OptionsPattern[]] := 
   Module[{\[Zeta],\[Lambda],\[Delta], iko, numbs, fak, k, v, n, x, y, z, j, 
           cs, sp, cos, r12, cll, co, cn, pts, dir, rmat, con, blength},
   (*---numbering scheme---*)
    numbs = {{11, 1}, {1, 11}, {1, 9}, {9, 1}, {9, 11}, {11, 9}, {11, 7}, {11, 5}, {11, 3}, {3, 11}, 
    	     {3, 1}, {1, 3}, {1, 2}, {1, 10}, {10, 1}, {10, 9}, {9, 10}, {9, 8}, {9, 7}, {7, 9},
    	     {7, 11}, {7, 5}, {5, 7}, {5, 11}, {5, 3}, {3, 5}, {3, 4}, {3, 2}, {2,3}, {2, 1}, {2, 10}, 
    	     {10, 2}, {10, 12}, {10, 8}, {8, 10}, {8, 9}, {8, 7}, {7, 8}, {7, 6}, {6, 7}, {6, 5}, 
    	     {5, 6}, {5, 4}, {4, 5}, {4, 3}, {4, 2}, {2, 4}, {2, 12}, {12, 2}, {12, 10}, {12, 8}, 
    	     {8, 12}, {8, 6}, {6, 8}, {6, 12}, {6, 4}, {4, 6}, {4, 12}, {12, 4}, {12, 6}
    	    };
  (*---parameters---*)
   \[Zeta]   = r5 + r6/2;
   \[Lambda] = r5/(2 r5 + r6);
   \[Delta]  = \[Pi]/5;
   r12       = \[Zeta]  Sqrt[(5 + Sqrt[5])/2]; fak = r12/Sqrt[5];
  (*---vertices icosahedron---*)
   iko = Table[0, {12}];
   Do[
      v = {2 fak Cos[k \[Delta]], 2 fak Sin[k \[Delta]], (-1)^k fak};
      iko[[k]] = v
   , {k, 1, 10}];
   iko[[11]] = {0, 0, -r12}; 
   iko[[12]] = {0, 0,  r12};
  (*---generate buckyball---*)
   cll = Table[0, {60}];
   Do[
  	  j = numbs[[n, 1]]; k = numbs[[n, 2]];
      x = iko[[j, 1]] + \[Lambda]*(iko[[k, 1]] - iko[[j, 1]]);
      y = iko[[j, 2]] + \[Lambda]*(iko[[k, 2]] - iko[[j, 2]]);
      z = iko[[j, 3]] + \[Lambda]*(iko[[k, 3]] - iko[[j, 3]]);
      cll[[n]] = {{x, y, z}, "C"}
   , {n, 1, 60}];
  (*---rotate the cluster---*)
   {co, cn} = Transpose[cll];
   pts  = Select[co, Abs[#[[2]]] < 0.01 &];
   dir  = pts[[3]] + pts[[4]]; dir = dir/Norm[dir];
   rmat = RotationMatrix[{dir, {0, 0, 1}}];
   con = rmat.# & /@ co; cll = Transpose[{con // Chop, cn}];
  (*--- prepare the plot ---*) 
   If[OptionValue[GOPlot], 
   	  blength = Max[r5,r6]+ 0.001;   
      cs  = OptionValue[GOColorScheme];
      sp  = OptionValue[GOSphere];
      cos = OptionValue[GOCoordinateSystem];
      dir = OptionValue[GODirection];
      GTPlotCluster[cll, blength, 7, GOColorScheme -> cs, GOPlot -> True,
                    GOSphere -> sp, GOCoordinateSystem -> cos, GODirection -> dir
                   ], 
      Return[cll]
   ]
 ]
   


(*
***)


(****e* /GTTubeParameters
! NAME
!  GTBTubeParameters
! AUTHOR
!  W. Hergert
! PACKAGE
!   CrystalStructure.m
! MODIFICATION HISTORY
!   * 10.10.2015 : first version
!   * 11.11.2015 : revision, GOLattice changed to rule
!   * 04.04.2016 : from Test.m to CrystalStructure.m
! USAGE
!  Calculates main geometric properties of carbon nanotubes (n,m)
! INPUT
!  n,m - characteristics of the nanotubes
! OUTPUT
!  list of the characteristic parameters
! GTPack OPTIONS
!  GOVerbose 
!     - True : a table with all the data is printed
!     - False : no table (Standard)
!  GOLattice 
!     - set lattice constant to a value (standard =1)
! STANDARD OPTIONS
!   - 
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  Reference_Symbols : GTTubeParameters.nb
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
! There was problem with GOLattice. In the first version it was used in the form GOLattice->value.
! Usually it is used in the form GOLattice-> substitution rule. In the present form it was changed 
! accordingly. A simple application of the rule was not possible, because the internally used "a" is
! different from the "a" in the substitution rule {a->2.49}. Instead of applying the rule, the number 
! is read out of the list {a->2.49}.
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTubeParameters[n_, m_, OptionsPattern[] ]:=Module[
	{a, a1, a2, ch, chl, dt, chla, t, \[Theta], dr, t1, t2, tl, test, np, nq, i, j, test1, n60, mn,
	 \[Tau], \[Psi], res, nu, npl, nql,prt,rule},
  (*--- lattice vectors ---*)
  rule=OptionValue[GOLattice];
  If[rule==={},
  	 a=1,
  	 a=rule[[1,2]]
  ];
  a1 = {Sqrt[3] a/2, a/2}; a2 = {Sqrt[3] a/2, -a/2};
  (*--- chiral vector ---*)
  ch = n*a1 + m*a2; chl = Norm[ch]; 
  dt = chl/\[Pi]; chla = Sqrt[n^2 + m^2 + n m];
  (*--- chiral angle ---*)
  t = (2 n + m)/2/chla; \[Theta] = ArcCos[t];
  While[Abs[\[Theta]] >= \[Pi]/6, 
   If[\[Theta] < 
     0, \[Theta] = \[Theta] + \[Pi]/6, \[Theta] = \[Theta] - \[Pi]/6]];
  (*--- translational vector ---*)
  dr = GCD[2 n + m, 2 m + n];
  t1 = (2 m + n)/dr; t2 = -(2 n + m)/dr;
  tl = Sqrt[3] chla/dr;
  (*--- number of hexagons ---*)
  nu = 2 (n^2 + m^2 + n m)/dr;
  (*--- symmetry vector ---*)
  npl = {}; nql = {};
  If[t1 == 0, n60 = 1, n60 = t1];
  Do[
     Do[
        test = t1*j - t2*i;
        If[test == 1,
           np = i; nq = j; test1 = m*np - n*nq,
           test1 = -1
        ];
        If[test1 > 0 && test1 <= nu,
           npl = Append[npl, np]; nql = Append[nql, nq],
           None
        ]
    , {j, -Abs[t2], Abs[t2]}]
   , {i, -Abs[n60], Abs[n60]}];
   If[Length[npl] == 1, 
   	  np = npl[[1]]; nq = nql[[1]], 
      Print["Error: ", npl, " ", nql]
   ];
   (*--- number of T in NR ---*)
   mn = m*np - n*nq;
   (*--- pitch of R ---*)
   \[Tau] = mn/nu*tl;
   (*--- rotation angle of R ---*)
   \[Psi] = 2 \[Pi]/nu;
   (*--- nanotube data ---*)
   res = {a, {n, m}, chla, dt, \[Theta], {t1, t2}, tl, nu, {np, nq}, 
   \[Tau], \[Psi], mn};
   (*--- Print data ---*)
  If[OptionValue[GOVerbose] == True,
   prt = {{"Property", 
      "Value"}, {"lattice constant" <> If[a === 1, " ", "(in \[Angstrom])"], 
      a}, {"Chiral vector \!\(\*SubscriptBox[\(C\), \(h\)]\) (in \
units \!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \
\(2\)]\))", {n, 
       m}}, {"Length of \!\(\*SubscriptBox[\(C\), \(h\)]\) (units of \
a)", chla}, {"Diameter " <> If[a === 1, "(units of a)", "(in \[Angstrom])"], 
      dt},
     {"chiral angle \[Theta]", \[Theta]},
     {"Translational vector T (in units \!\(\*SubscriptBox[\(a\), \(1\
\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\))", {t1, t2}},
     {"Length of T (units of a)", tl},
     {"Number of hexagons N in unit cell", nu},
     {"Symmetry vector R (in units \!\(\*SubscriptBox[\(a\), \
\(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\))", {np, nq}},
     {"Pitch \[Tau] of R", \[Tau]}, {"Rotation angle \[Psi]", \[Psi]},
     {"Number M of T in N\[CenterDot]R", mn}};
   Print[Grid[prt, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> 
    GTDividerColor1}},Background -> {None, {GTBackGroundColor1}}]], 
   None];
  Return[res]
  ]
  
  
(*
***)

 



(****e* /GTTubeStructure
! NAME
!  GTTubeStructure
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!  * 07.10.2015 : first version
!  * 18.04.2016 : from Test.m to CrystalStructure.m
!  * 27.06.2018 : check header and documentation
! USAGE
!  Out of the graphene layer a part is cut which is rolled to form a unit cell
!  of the nanotube. Several unit cells can be plotted 
! INPUT
!  o nt,mt   - characteristics of the nanotube
!  o cluster - a cluster containing the C-atoms of the carbon sheet
!  o ncells  - number of unitcells to be constructed
!  o dist    - shell number to be used for the plot of bonds
! OUTPUT
!  plot, graphics object, list of coordinates
! OPTIONS
!  GOVerbose 
!     - False : no additional information (standard)
!     - True  : additional information about tube parameters and construction
!  o GOBonds  
!     -  True : bonds are plotted (standard)
!     - False : no bonds
!  o GOPlot   
!     -  True  : the nanotube will be plotted (standard)
!     -  False : list with atomic coordinates and graphic object
! STANDARD OPTIONS
!   - 
! GTPack MODULES
!  GTTubeStructure, GTPointInSquareQ
! GTPack NOTEBOOKS 
!  GTNanotubes
! DESCRIPTION
!  -
! LITERATURE
!  R. Saito, G. Dresselhaus, M.S: Dresselhaus, Physical Properties of Carbon Nanotubes, Imperial College Press 2007
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!-
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTubeStructure[nt_, mt_, cluster_, ncells_, dist_, OptionsPattern[]] := 
  Module[{a, a1, a2, v1, v2, v3, v4, n, m, rad,lch, t1, t2, lt, cl1,  pt, i, phi, nat, x, y, z, alph, l, atoms, bonds,
  	      verb,bopt,oplot,parms,vert,uc1,uc,new,j,tube,blist,test,blength,graph},
  (*--- eveluation of options ---*) 
  a=1;
  verb = OptionValue[GOVerbose];
  bopt  = OptionValue[GOBonds];
  oplot = OptionValue[GOPlot];
  (*--- find the parameters ---*) 
  parms = GTTubeParameters[nt, mt, GOVerbose -> verb, GOLattice -> {}];
  (*--- construction of lattice ---*)
  a1 = {Sqrt[3] a/2, a/2}; 
  a2 = {Sqrt[3] a/2, -a/2};
  (*--- parameters from list ---*)
  n = parms[[2, 1]];
  m = parms[[2, 2]];
  lch = parms[[3]];
  rad = lch/(2 \[Pi]);
  t1 = parms[[6, 1]];
  t2 = parms[[6, 2]];
  lt = parms[[7]];
  (*--- vertices for cut ---*)
  v1 = {0., 0.};
  v2 = n*a1 + m*a2;
  v4 = t1*a1 + t2*a2;
  v3 = v2 + v4;
  vert = {v1, v2, v3, v4};
  (*--- coordinates in plane ---*) 
  cl1 = Take[Transpose[Transpose[cluster][[1]]], {1, 2}] // 
    Transpose;
  (*--- cut unit cell ---*)
  uc = {};
  Do[
     pt = cl1[[i]];
     If[GTPointInSquareQ[vert, pt],
        uc = Append[uc, pt],
        None
     ]
  , {i, 1, Length[cl1]}];
  (*--- multiple cells ---*)
  If[ncells > 1,
     Do[
        Do[new = uc[[i]] + j*v4;
           uc = Append[uc, new]
        , {i, 1, Length[uc]}]
     , {j, 1, ncells - 1}],
     None
  ];
  (*--- normal positions ---*)
  phi = ArcCos[v2.{1, 0}/lch];
  m = {{Cos[phi], Sin[phi]}, {-Sin[phi], Cos[phi]}};
  uc1 = Map[m.# &, uc];
  (*--- construct the tube ---*)
  nat = Length[uc1];
  tube = {};bonds = {};
  Do[
  	 pt = uc1[[i]];
     z = pt.{0, 1};
     l = pt.{1, 0};
     alph = 2*\[Pi]*l/lch;
     x = rad*Cos[alph];
     y = rad*Sin[alph];
     tube = Append[tube, {x, y, z}//Chop]
  , {i, 1, nat}];
  atoms = Map[Sphere[#, 0.1] &, tube];
  (*--- find the bonds ---*)
  If[bopt,
     If[dist > 0,
  (*--- calculate bond length ---*)
        blist = {};
        Do[
           Do[
              test = Norm[tube[[i]] - tube[[j]]];
              blist = Append[blist, test];
           , {j, 1, i - 1}]
        , {i, 1, Length[tube]}];
   (*--- filter bond lengths ---*)   
        blist = Select[Floor[blist*100000]/100000 // Union, # > 0 &];
        If[verb,
           Print["bond lengths: " , Take[blist*1., {1, 10}]],
           None
        ];
        blength = blist[[dist]] + 0.001,
        blength = Abs[dist]
     ];
   (*--- construct all bonds ---*)
     Do[
        Do[test=Norm[tube[[i]] - tube[[j]]];
           If[test <= blength&&test>0,
              bonds = Append[bonds, Cylinder[{tube[[i]], tube[[j]]}, 0.02]], 
              None
           ]
        , {i, 1, j-1}]
      , {j, 1, nat}],
      None
  ];
  (*--- decide output ---*)
  graph = Join[atoms, bonds];
  If[oplot,
     Show[Graphics3D[graph]],
     Return[{tube,graph}]
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
!  CrystalStructure.m
! MODIFICATION HISTORY
!   * 29.05.2016 : first version
!   * 10.10.2018 : implementation in CrystalStructure.m, check of Header
!   * 18.06.2021 : check of arguments at the beginning.
! USAGE
!  GTSymmetryElementQ[symmetry,cluster] checks if the symmetry elements leaves the cluster invariant.
! INPUT
!  symmetry - one symmetry element, or a list of symmetry elements
!  cluster  - ensemble of atoms to check with respect to symmetry 
! OUTPUT
!  logical value True or False
! GTPack Options
!  o GOVerbose 
!     - False : no additional information (standard)
!     - True  : additional information about tube parameters and construction
!  o GOPrecision  
!     -  10^(-7) - epsilon for numerical comparison
!  STANDARD OPTONS
!   -
! GTPack MODULES
!  GTGetMatrix
! GTPack NOTEBOOKS 
!  GSymmetryElementQ in Wolfram_Devel/4_Discrete_Symmetry_Groups/New_Commands
! DESCRIPTION
!  The need for such a command appeared in connection with the description of the nanotubes. It is possible to show that
!  the Dn groups are point groups of the nanotubes.
! LITERATURE
!
! TODO
!  more tests
! RELEASE
!  1.0.1
! PROBLEMS
!  It might be that it was not included in the first release, because Test.m wasn't included.
!--------------------------------------------------------------------------------
! SOURCE
*)

GTSymmetryElementQ::error = "Wrong order of arguments: symmetry element first, than cluster."

GTSymmetryElementQ[sym_, cluster_, OptionsPattern[]] := Module[{verb,ns,mats,clv,qtf,cl1,i,test,qtf1,log,pp,prec},
   (*--- options ---*)
   verb = OptionValue[GOVerbose];
   prec =OptionValue[GOPrecision];
   (*--- matrix form of symmmetry elements ---*)
   If[Head[sym] === List&&Head[cluster]===Symbol,
   	  Message[GTSymmetryElementQ::error],
   	  None
   ];
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
      test     = Compare[cl1, clv,GOPrecision->prec];
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




(****e* /GTChangeCoord
! NAME
!  GTChangeCoord
! AUTHOR
!  W. Hergert
! PACKAGE
!  CrystalStructure.m
! MODIFICATION HISTORY
!   * 11.02.2023 : first version
! USAGE
!  GTChangeCoord[structure] recalculates, if necessary positions of basis atoms to cartesian coordinates. 
! INPUT
!  structure - a GTPack crystal strucutre file
! OUTPUT
!  crystal structure file with atomic positions in Cartesian Coordinates
! GTPack Options
!  o GOVerbose 
!     - False : no additional information (standard)
!     - True  : additional information about tube parameters and construction
!  STANDARD OPTONS
!   -
! GTPack MODULES
!   -
! GTPack NOTEBOOKS 
!   -
! DESCRIPTION
!  In the original database GTPack.struc the positions of the basis atoms are given in Cartesian coordinates.
!  A common version is to use fractions of the lattice vectors for those positions. In a new database this version
!  is used. To be able to use old database notation and the new style, the command performs the calculation to Cartesian
!  coordinates.
!  - Length[structure]=8   -> no recalculation necessary
!  - Length[structure]=9 && structure[[9]]="C"   -> no recalculation necessary
!  - Length[structure]=9 && structure[[9]]="R"   -> ecalculation is done
!
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.33
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)

GTChangeCoord::coord =  "Type of coordinates is `1`. It has to be R or C."
                         

GTChangeCoord[struc_, OptionsPattern[]] := Module[
	{ln, st, struc1, pos, gv, nb, apos, npos, apos1, verb}, 
     verb   = OptionValue[GOVerbose];
     struc1 = struc;
  (*---check type of data set---*)ln = Length[struc];
     If[ln == 9,
        st = struc1[[9]];
        If[Intersection[{st}, {"R", "C"}] == {},
           Message[GTChangeCoord::coord, st]; Abort[],
           None        ]
     ];
  (*---recalculate dataset if necessary---*)
     If[ln == 8,
        If[verb,
           Print["No recalculation to Cartesian coordinates."];
           struc1 = Append[struc1, "C"],
           None
        ],
        If[ln == 9 && st == "C",
           If[verb,
              Print["No recalculation to Cartesian coordinates."],
              None
           ],
           pos   = struc1[[7]];
           gv    = struc1[[6]];
           nb    = Length[pos];
           {apos, npos} = Transpose[pos];
           apos1 = apos.gv;
           struc1[[7]] = Transpose[{apos1, npos}];
           struc1[[9]] = "C";
           If[verb, 
              Print["Recalculation to Cartesian coordinates performed."], 
              None
           ]
        ]
     ];
     Return[struc1]
]

                        

(*
***)
 

(****e* /GTAtomsInCell
! NAME
!  GTAtomsInCell
! AUTHOR
!  W. Hergert
! PACKAGE
! CrystalStructure.m
! MODIFICATION HISTORY
!  * 02.12.2021 - 1st version 
!  * 12.02.2023 - implemented in package
! USAGE
!  GTAtomsInCell[object,cluster,plot] finds all atoms in a region of the cluster defined by object.
! INPUT
!  * object    - graphics object for region or list for definition of a parallelepiped
!  * cluster   - the cluster of atoms
!  * plot      - (distance, size} necessary if result is plotted optional argument set to {0.5,2}
! OUTPUT
!   graphics output if GOPlot->True,  list of atoms in the region
!  
! GTPack OPTIONS
!  * GOVerbose:
!
!     - False - no additional information
!     - True  - detailed information 
!  * GOPlot: 
!
!     - False - no plot
!     - True  - plot of cluster nad region
!
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTPlotCluster
!  
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  If a larger cluster is constructed it is sometimes hard to inspect the local arrangement of the atoms. In the construction 
!  of the cluster all atoms inside a sphere are found. Here we can cut out of  this cluster all atoms in a parallelepiped.  
!  The parallelepiped can be defined before using GTAtomsInCell or the parallelepiped is constructed inside.
! LITERATURE
!-
! TODO
!  In the moment the region, that can be used is a parallelepiped. It should be possible to extend to other regions too. 
!  I have problems with the RegionMember function in other cases than parallelepiped. If this problem is solved, the 
!  extension is more or less trivial.
! RELEASE
!  1.33
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTAtomsInCell::error      = "Vectors are not linear independent.";
GTAtomsInCell::definition = "This is neither a definition of a graphics object nor a list to define a parallelepiped.";

GTAtomsInCell[object_, cluster_, plc_ : {0.5, 2}, OptionsPattern[]] := Module[
	{cl1, sel, mf, cl2, nat, verb, plt,object1,cent,vec,rot,rotm,lu,vecs,cube,g2},
(*---  options ---*)
     verb = OptionValue[GOVerbose];
     plt  = OptionValue[GOPlot];
(*--- prepare the region ---*)
     Switch[Head[object],
            List,
            If[verb,
               Print["Construction of parallelepiped"],
               None
            ];
            cent = object[[1]];
            vec  = object[[2]];
            rot  = object[[3]];
            lu   = Det[vec];
            If[lu == 0.,
               Message[GTAtomsInCell::error];Return[],
               None
            ];
            If[rot == {},
               vecs = vec,
               rotm = RotationMatrix[rot];
               vecs = rotm.# & /@ vec
            ];
            object1 = Parallelepiped[cent, vecs],           
            Parallelepiped,
            object1 = object,            
            _,
            Message[GTAtomsInCell::definition]; Return[]
      ];
      nat = Length[cluster];
      cl1 = Transpose[cluster][[1]];
      mf  = RegionMember[object1];
      sel = mf[cl1]; 
      cl2 = {}; 
      Do[
   	     If[sel[[i]] == True, 
   	  	    cl2 = Append[cl2, cluster[[i]]], 
   	  	    None
   	     ]
      ,{i,1, nat}];
      If[verb,
         Print[Length[cl2], " atoms selected from ", nat];
         Print[Grid[Transpose[cl2][[2]] // Tally // Transpose, Frame -> All,Background -> {None, 1 -> GTBackGroundColor1}]],
         None
      ];
      If[plt,
         cube = Graphics3D[{EdgeForm[{Thick, Black}], FaceForm[{LightRed, Opacity[0.2]}], object1}, 
                AspectRatio -> 1];
         g2   = GTPlotCluster[cluster, plc[[1]], plc[[2]], GOSphere -> False,GOCoordinateSystem -> True];
         Print[Show[cube, g2]],
         None
      ];
      Return[cl2]
  ]

(*
***)  

(****e* /GTAtomsInCell
! NAME
!  GTAtomsInCell
! AUTHOR
!  W. Hergert
! PACKAGE
! CrystalStructure.m
! MODIFICATION HISTORY
!  * 02.12.2021 - 1st version 
!  * 12.02.2023 - implemented in package
! USAGE
!  GTAtomsInCell[object,cluster,plot] finds all atoms in a region of the cluster defined by object.
! INPUT
!  * object    - graphics object for region or list for definition of a parallelepiped
!  * cluster   - the cluster of atoms
!  * plot      - (distance, size} necessary if result is plotted optional argument set to {0.5,2}
! OUTPUT
!   graphics output if GOPlot->True,  list of atoms in the region
!  
! GTPack OPTIONS
!  * GOVerbose:
!
!     - False - no additional information
!     - True  - detailed information 
!  * GOPlot: 
!
!     - False - no plot
!     - True  - plot of cluster nad region
!
! STANDARD OPTIONS
!  
! GTPack MODULES
!  GTPlotCluster
!  
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  If a larger cluster is constructed it is sometimes hard to inspect the local arrangement of the atoms. In the construction 
!  of the cluster all atoms inside a sphere are found. Here we can cut out of  this cluster all atoms in a parallelepiped.  
!  The parallelepiped can be defined before using GTAtomsInCell or the parallelepiped is constructed inside.
! LITERATURE
!-
! TODO
!  In the moment the region, that can be used is a parallelepiped. It should be possible to extend to other regions too. 
!  I have problems with the RegionMember function in other cases than parallelepiped. If this problem is solved, the 
!  extension is more or less trivial.
! RELEASE
!  1.33
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTCellPlot[vec_, bas_, pldat_, OptionsPattern[]] := Module[
	{ts,ah,sc,p,s1,si,co,v1,v2,v3,ko,kv,koord,kop,t1,t2,t3,po,pl1,ba,gb,gpo,go,obj},
 (*--- standard definitions for plot ---*)
     ts = 0.01; 
     ah = 0.05; 
     sc = 1.05;
 (*--- options ---*)	
     p  = OptionValue[FontSize];
     co = OptionValue[GOCoordinateSystem];
 (*--- plot data: scaling, sphere size  ---*)   
     s1 = pldat[[1]];
     si = pldat[[2]];
 (*--- lattice vactors ko and coordinate system  kv ---*)
     ko = kv = Table[0, {3}]; {v1, v2, v3} = vec;
     ko[[1]] = {Red, Arrowheads[ah], Arrow[Tube[{{0, 0, 0}, v1}, ts]], {Black, 
                Text[Style["\!\(\*SubscriptBox[\(a\), \(1\)]\)", FontSize -> p], v1/2]}};
     ko[[2]] = {Red, Arrowheads[ah], Arrow[Tube[{{0, 0, 0}, v2}, ts]], {Black, 
                Text[Style["\!\(\*SubscriptBox[\(a\), \(2\)]\)", FontSize -> p], v2/2]}};
     ko[[3]] = {Red, Arrowheads[ah], Arrow[Tube[{{0, 0, 0}, v3}, ts]], {Black, 
                Text[Style["\!\(\*SubscriptBox[\(a\), \(3\)]\)", FontSize -> p], v3/2]}};
     kv[[1]] = {Blue, Arrowheads[ah], Arrow[Tube[{{0, 0, 0}, s1 {1, 0, 0}}, ts]], {Blue, 
                Text[Style["x", FontSize -> p], s1 sc {1, 0, 0}]}}; 
     kv[[2]] = {Blue, Arrowheads[ah], Arrow[Tube[{{0, 0, 0}, s1 {0, 1, 0}}, ts]], {Blue, 
                Text[Style["y", FontSize -> p], s1 sc {0, 1, 0}]}}; 
     kv[[3]] = {Blue, Arrowheads[ah], Arrow[Tube[{{0, 0, 0}, s1 {0, 0, 1}}, ts]], {Blue, 
                Text[Style["z", FontSize -> p], s1 sc {0, 0, 1}]}};
 (*--- Plots --- *)
  koord = Graphics3D[kv, AspectRatio -> 1, AxesOrigin -> {0, 0, 0}];
  kop   = Graphics3D[ko, AspectRatio -> 1, AxesOrigin -> {0, 0, 0}];
 (*--- Plot the cell ---*)
  pl1   = Graphics3D[{Opacity[0.1], LightGray, 
     Parallelepiped[{0, 0, 0}, vec]}, AspectRatio -> 1];
 (*--- Plot atoms in primitive cell ---*)
  t1  = vec[[1]];
  t2  = vec[[2]]; 
  t3  = vec[[3]];
  po  = {{0, 0, 0}, t1, t2, t3, t1 + t2, t1 + t3, t2 + t3, t1 + t2 + t3};
  gpo = {Map[{RGBColor[.0, .5, .5], Sphere[#, si]} &, po]};
  go  = Graphics3D[gpo, AspectRatio -> 1, BoxStyle -> Directive[Opacity[0.9]], AxesOrigin -> {0, 0, 0}];
 (*---  basis atoms if any ---*)
  If[bas == {},
     None,
     ba = Map[Sphere[#, si] &, bas];
     gb = Graphics3D[ba, AspectRatio -> 1, BoxStyle -> Directive[Opacity[0.9]], AxesOrigin -> {0, 0, 0}];
  ];
  If[bas == {},
     If[co,
        obj = {koord, pl1, kop, go},
        obj = {pl1, kop, go}
     ],
     If[co,
        obj = {koord, pl1, kop, go, gb},
        obj = {pl1, kop, go, gb}
     ]
  ];
  Show[obj, Axes -> True, Ticks -> None]
]

(*
***)  

(*-------------------------- Attributes ------------------------------*)

Attributes[GTLoadStructures]={Protected, ReadProtected}
Attributes[GTGroupNotation]={Listable,Protected, ReadProtected}
Attributes[GTSaveStructures]={Protected, ReadProtected}
Attributes[GTInstallStructure]={Protected, ReadProtected}
Attributes[GTAllStructures]={Protected, ReadProtected}
Attributes[GTGetStructure]={Protected, ReadProtected}
Attributes[GTPlotStructure]={Protected, ReadProtected}
Attributes[GTPlotStructure2D]={Protected, ReadProtected}

End[]  (* End Private Context *)



EndPackage[]
