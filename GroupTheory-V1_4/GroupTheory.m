(****m* /GroupTheory.m
!
! NAME
!  GroupTheory.m
! AUTHOR
!  W. Hergert, M. Geilhufe, S. Schenk
! MODIFICATION HISTORY
!  27.12.2017 : header of package
! USAGE
!  Group theory package to be used with the corresonding book
!
! DESCRIPTION
!  The package was developed over a series of years. This is the final form
!  which corresponds to the description of the book.
!  During the years we have learned that several things can be done better.
!  This willbe considered in the next verson of th package.
! LITERATURE
!  W. Hergert, M. Geilhufe 
!  Group Theory in Solid State Physics and Photonics. - Problem Solving
!  with Mathematica 
! TODO
!  Find all the errors :-))
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
!
***)

BeginPackage["GroupTheory`",
	        {"GroupTheory`Install`",
	         "GroupTheory`Auxiliary`",
	         "GroupTheory`Additional`",
	         "GroupTheory`AngularMomentum`",
	         "GroupTheory`CrystalField`",
	         "GroupTheory`Symbols`",
	         "GroupTheory`Basic`", 
	         "GroupTheory`RepresentationTheory`",
	         "GroupTheory`RepresentationTheorySG`",
	         "GroupTheory`CrystalStructure`",
	         "GroupTheory`TightBinding`",
	         "GroupTheory`ElectronicStructure`",
	         "GroupTheory`Photonics`",
	         "GroupTheory`PseudoPotential`",
	         "GroupTheory`Molecules`",
	         "GroupTheory`Vibrations`",
	         "GroupTheory`Lattice`",
	         "GroupTheory`SimPack`",
	         "GroupTheory`LandauTheory`",
(*--- Third party packages ---*)	         
	         "GroupTheory`ThirdParty`",
	       
(*--- Test of new routines ---*)
                 "GroupTheory`Wyckoff`"
            (*     "GroupTheory`Test1`"	   *)     
	        }
	        ]

GTGOFast		::usage = "GTGOFast[\*StyleBox[\"logical\", \"TI\"]] switches between standard (False) and expert (True) mode. In expert mode, input of modules is not testet for errors! \[NonBreakingSpace]\*ButtonBox[StyleBox[\"\[RightSkeleton]\", \"SR\"], BaseStyle->\"Link\", ButtonData->\"paclet:GroupTheory/ref/GTGOFast\"]"

(*--- Installation of axes ---*)
tmp1 = PrintTemporary["The Mathematica Group Theory Package GTPack (Version 1.4)"];
tmp2 = PrintTemporary["-------------------------------------------------------------------"];
tmp1a = PrintTemporary["Version release date: April 5th 2023"];
tmp2b = PrintTemporary["Copyright W. Hergert, R.M. Geilhufe"];
tmp2c = PrintTemporary[""];
tmp2d = PrintTemporary["More information and updates can be found at: https://gtpack.org"];
tmp2e = PrintTemporary["GTPack is an academic project and free to use. To support the development of GTPack,"];
tmp2f = PrintTemporary["we ask you to cite the following publications:"];
tmp2g = PrintTemporary["1) W. Hergert, R. M. Geilhufe, Group Theory in Solid State Physics and Photonics: Problem Solving with Mathematica, Wiley-VCH, 2018"];
tmp2h = PrintTemporary["2) R. M. Geilhufe, W. Hergert, Frontiers in Physics, 6:86, 2018"];
tmp2i = PrintTemporary[""];
tmp1b = PrintTemporary["Installation of symmetry elements for "<>If[gtactvpsv>0,"passive","active"]<>" rotations."];
tmp2 = PrintTemporary["-------------------------------------------------------------------"];
tmp3 = PrintTemporary["install symmetry axes:"];

tmp4 = PrintTemporary["x-Axis"];
GTInstallAxis["x", ex,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["y-Axis"];
GTInstallAxis["y", ey,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["z-Axis"];
GTInstallAxis["z", ez,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["a-Axis"];
GTInstallAxis["a", ea,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["b-Axis"];
GTInstallAxis["b", eb,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["c-Axis"];
GTInstallAxis["c", ec,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["d-Axis"];
GTInstallAxis["d", ed,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["e-Axis"];
GTInstallAxis["e", ee,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["f-Axis"];
GTInstallAxis["f", ef,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["\[Alpha]-Axis"];
GTInstallAxis["\[Alpha]", e\[Alpha],GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["\[Beta]-Axis"];
GTInstallAxis["\[Beta]", e\[Beta],GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["\[Gamma]-Axis"];
GTInstallAxis["\[Gamma]", e\[Gamma],GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["\[Delta]-Axis"];
GTInstallAxis["\[Delta]", e\[Delta],GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["A-Axis"];
GTInstallAxis["A", eA,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["B-Axis"];
GTInstallAxis["B", eB,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["C-Axis"];
GTInstallAxis["C", eC,GOVerbose->None]
NotebookDelete[tmp4];

tmp4 = PrintTemporary["D-Axis"];
GTInstallAxis["D", eD,GOVerbose->None]
NotebookDelete[tmp4];
NotebookDelete[tmp3];

tmp3 = PrintTemporary["Extracting archives:"];
tmp4 = PrintTemporary["datasets"];
GTExtractDatasets[];
NotebookDelete[tmp4];
NotebookDelete[tmp3];

NotebookDelete[tmp2];
NotebookDelete[tmp2a];
NotebookDelete[tmp2b];
NotebookDelete[tmp2c];
NotebookDelete[tmp2d];
NotebookDelete[tmp2e];
NotebookDelete[tmp2f];
NotebookDelete[tmp2g];
NotebookDelete[tmp2h];
NotebookDelete[tmp2i];
NotebookDelete[tmp1b];
NotebookDelete[tmp1a];
NotebookDelete[tmp1];

(****v* /GTGOFast
! NAME
!  GTGOFast
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  GroupTheory.m
! MODIFICATION HISTORY
!  05/03/2013 : first version
!  28.12.2017   : check header and documentation
! USAGE

!  "GTGOFast[logical] switches between standard (False) and expert (True) mode. In expert mode, 
!  input of modules is not testet for errors!
! INPUT
!  logical
! OUTPUT
!  
! GTPack OPTIONS
!  
! STANDARD OPTIONS
! 
! GTPack MODULES
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

GTGOFast[log_]:= Module[{},
  	Unprotect[GOFastValue];
	GOFastValue=log;
	Protect[GOFastValue];
	
	(*--- Options Basic.m  ---*)

	SetOptions[GTGenerators,GOFast->GOFastValue];
	SetOptions[GTGetSubGroups,GOFast->GOFastValue];
	SetOptions[GTGroupOrder,GOFast->GOFastValue];
	SetOptions[GTMultTable,GOFast->GOFastValue];
	
	(*--- Options RepresentationTheory.m ---*)
	
	SetOptions[GTClasses,GOFast->GOFastValue];
	SetOptions[GTClassMult,GOFast->GOFastValue];
	SetOptions[GTClassMultTable,GOFast->GOFastValue];
	
	(*--- characters ---*)
	
	SetOptions[GTCharacterTable,GOFast->GOFastValue];
	SetOptions[GTCharProjectionOperator,GOFast->GOFastValue];
	SetOptions[GTProjectionOperator,GOFast->GOFastValue];
	SetOptions[GTWignerProjectionOperator,GOFast->GOFastValue];
	
	(*--- Ireps ---*)
	
	SetOptions[GTClebschGordanSum,GOFast->GOFastValue];
	SetOptions[GTGetIrep,GOFast->GOFastValue];
	SetOptions[GTIrep,GOFast->GOFastValue];
	
	(*--- Direct Product ---*)
	
	SetOptions[GTClebschGordanCoefficients,GOFast->GOFastValue];
	SetOptions[GTDirectProductRep,GOFast->GOFastValue];
	
	(*--- Options Test.m ---*)
	SetOptions[GTGroupOfK,GOFast->GOFastValue];
]


(*-------------------------- Options usages --------------------------*)
GOAccuracy     ::usage=""
GOBands        ::usage="is an option to select a number of bands for the calculation of the density of states."
GOBondCharges  ::usage="is an option to decide if bond charges should be used in the construction of the plane wave Hamiltonian."
GOBonds        ::usage="is an option to decide, whether bonds between the atoms should be drawn or not."
GOBravais      ::usage="is an option to decide if the lattice should be presented by a primitive cell or alternatively by a convention cell with basis."
GOColorScheme  ::usage="is an option to define a color scheme used for the plot of the structure, i.e. which colors are used for the different atoms."
GOCompact      ::usage="is an option which fixes the output format of the matrix."
GOData         ::usage="is an option to decide if the lattice vectors or a picture of the cell is required."
GODimension    ::usage="is an option which specifies the dimension of space used in the algorithm."
GOEigenvectors ::usage="is an option to control whether the eigenvectors will be calculated."
GOFast         ::usage="is an option which specifies if the input is validated."
GOFastValue    ::usage="is a global variable defining the default option value of GOFast."
GOGeneralPositionIreps ::usage="is an option which controls if the character table and representation matrices of the full space group or the character table of the space group of k is printed."
GOGroupOrder   ::usage="is an option which specifies the used group order in the algorithm."
GOHarmonics    ::usage="is an option to decide whether real or complex spherical harmonics are used."
GOIrepNotation ::usage="is an option which specifies whether the notation of the irreducible representations is given in Bethe, Mulliken or Bouckaert notation."
GOLattice      ::usage="is an option to provide rules for the rescaling of lattice constants."
GOLeftCosets   ::usage="is an option which controls if the left cosets are given as an additional output."
GOMatrixType   ::usage="is an option which provides information about the structure of a matrix."
GOMethod       ::usage=""
GONames        ::usage="is an option to control the names of irreducible representations."
GONotation     ::usage="is an option which specifies whether the group notation is given in Schoenfliess or Herman-Mauguin notation."
GOPhotonic     ::usage="is an option which is used to distinguish between the calculation of an electronic or a photonic bandstructure."
GOPhPol        ::usage="is an option to define the polarization for band structure calculations of two-dimensional photonic crystals."
GOPlot         ::usage="is an option to decide, if the adjacency matrix is plotted."
GOPlotBands    ::usage="is an option to control the calculation of energy bands in the Brillouin zone."
GOPlotDos      ::usage="is an option to decide whether the density of states or the integrated density of states will be plotted."
GOPlotStyle    ::usage="is an option which specifies whether the connection graph is embedded or not."
GOPosition     ::usage="is an option during the construction of a cluster or shells."
GORepresentation ::usage="is an option to change the standard representation during the installation of point groups."
GOSelectCase   ::usage="is an option which specifies whether connections are shown only between the groups of the list or to all connected groups."
GOShift        ::usage="is an option to shift the energies of a calculated band structure by a fixed value."
GOSort         ::usage="is an option to define if a canonical order of the output of GTAdjacencyMatrix or GTLatShells should be performed."
GOStore        ::usage="is an option which controls the output of data to files."
GOTakeMod      ::usage="is an option which controls if the full translation vector of the product of element1 and element2 is given or if the translation vector is given mod the basis vectors."
GOTbBasis      ::usage="is an option to control the construction of symbols for the two-center tight-binding parameters."
GOTbEquivalence ::usage="is an option which decides whether equivalent atoms in a basis should be renamed."
GOTbLattice    ::usage="is an option to provide detailed information for the construction of tight-binding Hamiltonians from the shell construction, if more than one atom is present in the basis."
GOTbOrthogonal ::usage="is an option of GTBandStructure used to define if orthogonal basis sets are assumed in the tight-binding calculations or not."
GOTbRule       ::usage="is an option to apply substitution rules in the transformation of tight-binding matrix elements in two-center form."
GOTolerance    ::usage="defines the maximal allowed deviation to regard two numerical values as equal."
GOVerbose      ::usage="is an option which controls the output of additional information."


(*
***)


Begin["`Private`"]


End[]

EndPackage[]

