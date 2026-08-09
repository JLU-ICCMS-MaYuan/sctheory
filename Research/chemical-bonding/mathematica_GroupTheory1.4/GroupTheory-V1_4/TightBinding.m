(****m* /TightBinding.m
!       
! NAME
!  TightBinding.m
! AUTHOR
!  W. Hergert
! MODIFICATION HISTORY
!  * 9/05/13     : initially created and documented  
!  * August 2015 : Check of documentation of TightBinding.m
! USAGE
!  Tight-binding models for electronic structure calculations
!  
! GTPack MODULES
!
! --- Construction of two-center tight-binding Hamiltonians --- 
!
!  * GTTbSpinMatrix            - matrices to construct Hamiltonian with spin
!  * GTTbMatrixElement         - matrix elements of a TB-Hamiltonian
!  * GTTbSymbol2C              - construction of the symbols for two-center tight-binding parameters
!  * GTTbHamiltonianElement    - constructs an element of the TB-Hamiltonian in k-space
!  * GTTbHamiltonian           - full Hamiltonian in two-center form
!  * GTTbSpinOrbit             - inclusion of the spin-orbit interaction
!  * GTTbHamiltonianOnsite 	   - correction of the on-site energies (crystal field splitting)
!
! --- Construction of three-center tight-binding Hamiltonians ---
!
!  * GTSymmetryBasisFunctions 	- evaluates to which Irep a functio is basis function
!  * GTTbSymbol3C               - symbols for three-center tight-binding parameters
!  * GTTbMatrixElement3C          - gives part of matrix element from a certain shell in 3-center form
!  * GTTbIntegralRules          - finds the rules to minimize the number of parameters
!  * GTIrepInfo                 - extracts the information from representation matrices necessary for Hamiltonian construction
!  * GTTbNumberOfIntegrals      - number of independent parameters in the HAmiltonian
!  * GTTbVarList                - ?
!
! --- Real space formulation of tight-binding models ---
! 
!  * GTTbRealSpaceMatrix        - constructs the interaction of two atoms in a certain neigbour shell
!  * GTTbHamiltonianRS          - constructs TB Hamiltonian in real space
!  * GTFindStateNumbers         - finds numbers of eingenvalues of real space Hamiltonian in a given energy interval
!  * GTPlotStateWeights         - graphical representation of the weight of eigenstates
!
! --- Parameter sets and libraries for tight-binding models ---
!
!  * GTTbParameterNames	        - gives the parameter names up to a certain lmax and shell number
!  * GTTbParmToRule             - constructes a substitution rule to substitute hopping integrals by numerical values
!  * GTTbGetParameter           - get a certain parameter out of a given parameter set
!  * GTTbSetParameter           - sets the parameter sym of set to value 
!  * GTTbDatabaseInfo           - information about the parameter sets in the library
!  * GTTbDatabaseUpdate         - add a new parameter set to the library 		 
!  * GTTbDatabaseRetrieve       - retrieve a parameter set from the library 
!  * GTTbPrintParmSet           - pretty print of a parameter set
!
! ---  Wave functions ---
!
!  * GTTbAtomicWaveFunction     - constructs an atomic-like wavefunction
!  * GTTbBlochFunction          - builds a Bloch function from atomic like wavefunctions
!  * GTTbWaveFunction           - constructs the total wavefunction for a special band at a certain k-point
!
! ---  Symmetry Analysis of Band Structures ---
!
!  * GTTbSymmetrySingleBand     - symmetry analysis of a single band
!  * GTTbSymmetryPoint          - symmetry analysis at a certain k-point
!  * GTTbSymmetryBands          - symmetry analysis of a band structure
!
! --- Output of Hamiltonians as FORTRAN modules ---
!
!  * GTTbToFortran              - transformation of a tight-binding Hamiltonian in a FORTRAN module
!  * GTTbRSToFortran			- transformation of a real space tight-binding Hamiltonian in a FORTRAN module
!  * GTTbToFortranList          - prints a tight-binding Hamiltonian as FORTRAN code
!
!  --- Miscellaneous ---
!
!  * GTHamiltonianPlot           - plots the structure of the Hamiltonian
!  * GTHamiltonianList           - lists the elements of a tight-binding Hamiltonian
!  * GTWaveFunctionPlot          - plots a wavefunction
!  * GTTbOrbitalsFromBasis       - creates a list of orbital names from a basis description
!
! --- NanoTubes ---
!
!  * GTTbTubeBands               - tight-binding Hamiltonian for nanotubes
!
!  --- Wannier90 ---
!
!  * GTTbWannier90               - reads information from Wannier90 package
!  * GTTbWannier90H0amiltonian   - constructs tight-binding Hamiltonian from Wannier90 package
!
! --- internal modules ---
!
!  * GTTbSimplify                - rules to simplify two-center matrix elements    
!  * GTTbParmToFortran           - transformation of paramter names to a parameter field  for FORTRAN Export
!  * TbCompose                   - glue parts of Hamiltonian together (spin-dependent calculations)
!  * HydrogenRadial              - radial part of a Hydogen-like wave function
!
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The package contains all modules corresponding to tight-binding models. The modules start usually 
!  with GTTb#. TB-Hamiltonians can be automatically constructed in two-center form. Spin-Orbit interaction
!  can be included. The consideration of tight-binding models in three-center form are also possible.
!
! LITERATURE
!  Some basic literature of this part of GTPack is:
!  - J.C. Slater, G.F. Koster, Simplified LCAO method for the periodic potential problem, Phys. Rev. 94, 1498 (1954)
!  - M. Miasek, Tight-Binding Method for Hexagonal Close-Packed Structure, Phys. Rev. 107, 92 (1957)
!  - A.V. Podolskiy, P. Vogl, Compact expressions for the angular dependendce of tight binding Hamiltonian matrix elements,
!    Phys. Rev. B 69, 233101 (2004)
!  - M.D. Stiles, Generalized Slater-Koster method for fitting band structures, Phys. Rev. B 55, 4168 (1997)
!  - R.F. Egorov, B.I. Reser, and V.P. Shirkovskii, Consistent Treatment of Symmetry in the Tight Binding Approximation,
!    phys. stat. sol. 26, 391 (1968)
!  - V.N. Kuznetsov and A.N. Men, Symmetry and Analytical Expresssions of Matrix Components int the Tight-Binding Method, 
!    phys. stat. sol (b) 85, 95 (1978)
! 
!
***)

BeginPackage["GroupTheory`TightBinding`",{"GroupTheory`ElectronicStructure`","GroupTheory`Symbols`","GroupTheory`Lattice`","GroupTheory`Auxiliary`","GroupTheory`RepresentationTheory`","GroupTheory`Basic`","GroupTheory`Photonics`","GroupTheory`CrystalStructure`","GroupTheory`AngularMomentum`"}]

(*---------- Construction of two-center tight-binding Hamiltonians ----------*)
 GTTbHamiltonian           ::usage = "GTTbHamiltonian[\*StyleBox[\"basis,lattice\",\"TI\"]] constructs the \*StyleBox[\"k\",FontWeight->\"Bold\"]-dependent Hamiltonian from the information about \*StyleBox[\"basis\",\"TI\"] and \*StyleBox[\"lattice\",\"TI\"]."
 GTTbHamiltonianElement    ::usage = "GTTbHamiltonianElement[\*StyleBox[\"l1,m1,l2,m2,shell,shell vectors\",\"TI\"]] constructs the \*StyleBox[\"k\",FontWeight->\"Bold\"]-dependent contribution of \*StyleBox[\"shell\",\"TI\"] characterized by the \*StyleBox[\"shell vectors\",\"TI\"] to the tight-binding matrix element between functions of symmetry \*Cell[BoxData[FormBox[RowBox[{\"(\",RowBox[{SubscriptBox[\"l\",\"1\"],\",\",SubscriptBox[\"n\",\"1\"]}],\")\"}], TraditionalForm]],\"InlineMath\"] and \*Cell[BoxData[FormBox[RowBox[{\"(\",RowBox[{SubscriptBox[\"l\",\"2\"],\",\",SubscriptBox[\"n\",\"2\"]}],\")\"}], TraditionalForm]],\"InlineMath\"]."
 GTTbMatrixElement         ::usage = "GTTbMatrixElement[\*StyleBox[\"l1,m1,l2,m2,shell\",\"TI\"]] gives the decomposition of the tight-binding three-center integral between atom1 \*Cell[BoxData[FormBox[RowBox[{\"(\",RowBox[{SubscriptBox[\"l\",\"1\"],\",\",SubscriptBox[\"m\",\"1\"]}],\")\"}], TraditionalForm]],\"InlineMath\"] and atom2 \*Cell[BoxData[FormBox[RowBox[{\"(\",RowBox[{SubscriptBox[\"l\",\"2\"],\",\",SubscriptBox[\"m\",\"2\"]}],\")\"}], TraditionalForm]],\"InlineMath\"], when atom2 belongs to the neighborhood shell and atom2 is located in direction relative to atom1."
 GTTbSpinMatrix            ::usage = "GTTbSpinMatrix[\*StyleBox[\"angular momentum,s1,s2\",\"TI\"]] gives elementary spin matrices for tight-binding Hamiltonians."
 GTTbSpinOrbit             ::usage = "GTTbSpinOrbit[\*StyleBox[\"hamiltonian,spin-orbit interaction\",\"TI\"]] adds spin-orbit coupling to a given tight-binding \*StyleBox[\"Hamiltonian\",\"TI\"] due to a specified \*StyleBox[\"spin-orbit interaction\",\"TI\"]."
 GTTbSymbol2C              ::usage = "GTTbSymbols2C[\*StyleBox[RowBox[{SubscriptBox[\"l\",\"<\"],\",\",SubscriptBox[\"l\",\">\"],\",m,shell\"}],\"TI\"]] gives a symbol according to the nomenclature of the two-center approximation."
 GTTbHamiltonianOnsite     ::usage = "GTTbHamiltonianOnsite[\*StyleBox[\"Hamiltonian, structure\",\"TI\"]] corrcts the on-site elements for fcc,bcc, or hcp strucutre in the Hamiltonian."
(*---------- Construction of three-center tight-binding Hamiltonians  -------*)
 GTIrepInfo                ::usage = "GTIrepInfo[pg,pgct,ireps] information on representation matrices"
 GTSymmetryBasisFunctions  ::usage = "GTSymmetryBasisFunctions[\*StyleBox[\"character table,wave functions\",\"TI\"]] calculates to which irreducible representations the \*StyleBox[\"wave functions\",\"TI\"] are basis functions."
 GTTbNumberOfIntegrals     ::usage = "GTTbNumberOfIntegrals[\*StyleBox[\"point group,character tables,irreducible representations\",\"TI\"]] calculates the number of independent tight-binding integrals for functions transforming like irreducible representations of \*StyleBox[\"point group\",\"TI\"]. The \*StyleBox[\"character tables\",\"TI\"] of the groups \*Cell[BoxData[FormBox[SubsuperscriptBox[\"G\",\"p\",\"l\"], TraditionalForm]],\"InlineMath\"] have to be provided."
 GTTbSymbol3C              ::usage = "GTTBSymbol3C[\*StyleBox[\"irep1,row1,irep2,row2,qlp\",\"TI\"]] changes a formal name for the three-center paramer in the Slater&Koster form."
 GTTbMatrixElement3C       ::usage = "GTTbMatrixElement[repmat1,repmat2,mu,mup,vars,qlp,kv] gives part of matrix element from a certain shell in 3-center form."
 GTTbIntegralRules         ::usage = "GTTbIntegralRules[ireps1, ireps2] finds the rules to minimize the number of parameters."
 GTTbVarList               ::usage = "GTTbVarlist[]"
(*---------- Tight-binding real space formulation ---------------------------*)
 GTFindStateNumbers        ::usage = "GTFindStateNumbers[\*StyleBox[\"eigenvalues,emin,emax\",\"TI\"]] findes the number of states in the list \*StyleBox[\"eigenvalues\",\"TI\"] lying in the energy interval [\*StyleBox[\"emin,emax\",\"TI\"]]."
 GTPlotStateWeights        ::usage = "GTPlotStateWeights[\*StyleBox[\"cluster,distance,basis,wave function,scale\",\"TI\"]] demonstrates at which atoms of a \*StyleBox[\"cluster\",\"TI\"], the \*StyleBox[\"wave function\",\"TI\"] (a solution of a real space tight-binding Hamiltonian) has the heighest weights."
 GTTbHamiltonianRS         ::usage = "GTTbHamiltonianRS[\*StyleBox[\"cluster,adjacency matrix,basis\",\"TI\"]] constructs a tight-binding Hamiltonian in real space."
 GTTbRealSpaceMatrix       ::usage = "GTTbRealSpaceMatrix[\*StyleBox[\"atom1,atom2,shell\",\"TI\"]] constructs the interaction of \*StyleBox[\"atom1\",\"TI\"] and \*StyleBox[\"atom2\",\"TI\"] in a tight-binding Hamiltonian in a certain neighbor \*StyleBox[\"shell\",\"TI\"]."
 
(*---------- Parameter sets and libraries for tight-binding models ----------*)
 GTTbDatabaseInfo          ::usage = "GTTbDatabaseInfo[\*StyleBox[\"database\",\"TI\"]] gives information about the tight-binding parameter sets available in \*StyleBox[\"database\",\"TI\"]."
 GTTbDatabaseRetrieve      ::usage = "GTTbDatabaseRetrieve[\*StyleBox[\"database,parameter set\",\"TI\"]] loads a tight-binding \*StyleBox[\"parameter set\",\"TI\"] from a given \*StyleBox[\"database\",\"TI\"]." 
 GTTbDatabaseUpdate        ::usage = "GTTbDatabaseUpdate[\*StyleBox[\"database\",\"TI\"]] reads the \*StyleBox[\"database\",\"TI\"], adds a new parameter set and stores the \*StyleBox[\"database\",\"TI\"]."
 GTTbGetParameter          ::usage = "GTTbGetParameter[\*StyleBox[\"parameter set,parameter\",\"TI\"]] gives the value of a \*StyleBox[\"parameter\",\"TI\"] from a \*StyleBox[\"parameter set\",\"TI\"]."
 GTTbSetParameter          ::usage = "GTTbSetParameter[\*StyleBox[\"parameter set,parameter, value\",\"TI\"]] sets \*StyleBox[\"parameter\",\"TI\"] from \*StyleBox[\"parameter set\",\"TI\"] to \*StyleBox[\"value\",\"TI\"]  ."
 GTTbParameterNames        ::usage = "GTTbParameterNames[\*StyleBox[\"maximum angular momentum,maximum shell number\",\"TI\"]] creates a set of parameter names up to \*StyleBox[\"maximum angular momentum\",\"TI\"] and \*StyleBox[\"maximum shell number\",\"TI\"]."
 GTTbParmToRule            ::usage = "GTTbParmToRule[\*StyleBox[\"parameter set\",\"TI\"]] gives a rule to replace tight-binding symbols using a given \*StyleBox[\"parameter set\",\"TI\"]."
 GTTbPrintParmSet          ::usage = "GTTbPrintParmSet[\*StyleBox[\"database,parameter set\",\"TI\"]] prints a tight-binding \*StyleBox[\"parameter set\",\"TI\"] from a \*StyleBox[\"database\",\"TI\"]."
 
(*--------- Wave functions --------------------------------------------------*)
 GTTbAtomicWaveFunction    ::usage = "GTTbAtomicWaveFunction[\*StyleBox[\"atomic number,n,l,m,position\",\"TI\"]] gives the value of an atomic-like wave function of the \*StyleBox[\"atomic number\",\"TI\"] to the quantum numbers (\*StyleBox[\"n,l,m\",\"TI\"]) at \*StyleBox[\"position\",\"TI\"]."
 GTTbBlochFunction         ::usage = "GTTbBlochFunction[kvec,rvec,pos,atom] constructs a TB bloch function"
 GTTbWaveFunction          ::usage = "GTTbWaveFunction[kvec,rvec,pos,atom] constructs a TB bloch function"

(*--------- Symmetry Analysis of Band Structures ------------------------------*)
 GTTbSymmetrySingleBand    ::usage = "GTTbSymmetrySingleBand[\*StyleBox[\"coefficients, orbitals, charactertable\",\"TI\"]] performes the symmetry analysis for a band ."
 GTTbSymmetryPoint         ::usage = "GTTbSymmetryPoint[\*StyleBox[\"kpoint,coefficients, orbitals, minb,maxb,point group,reciprocal basis\",\"TI\"]] performes the symmetry analysis for at a certain k-point."
 GTTbSymmetryBands         ::usage = "GTTbSymmetryBandStructure[\*StyleBox[\"fileb,filew,klist,minb,maxb,point group,reciprocal basis, orbitals\",\"TI\"]] performes the symmetry analysis for a band structure." 
 
(*---------- Output of Hamiltonians as FORTRAN modules -----------------------*)
 GTTbToFortran             ::usage = "GTTbToFortran[\*StyleBox[\"Hamiltonian,maximal angular momentum,maximum shell number,file name\",\"TI\"]] transforms a \*StyleBox[\"k\",FontWeight->\"Bold\"]-space tight-binding Hamiltonian into a FORTRAN module."
 GTTbToFortranList         ::usage = "GTTbToFortranList[\*StyleBox[\"Hamiltonian,parameter set,database\",\"TI\"]] prints a Hamiltonian as FORTRAN code."

(*---------- Miscellaneous ---------------------------------------------------*)
 GTHamiltonianPlot         ::usage = "GTHamiltonianPlot[\*StyleBox[\"Hamiltonian,basis\",\"TI\"]] plots the structure of a \*StyleBox[\"hamiltonian\",\"TI\"], using information about the \*StyleBox[\"basis\",\"TI\"] employed in the construction of the \*StyleBox[\"hamiltonian\",\"TI\"]."
 GTHamiltonianList         ::usage = "GTHamiltonianList[\*StyleBox[\"Hamiltonian,names\",\"TI\"]] lists the elements of \*StyleBox[\"hamiltonian\",\"TI\"], using information about the \*StyleBox[\"names\",\"TI\"] of the orbitals in the basis."
 GTWaveFunctionPlot        ::usage = "GTWavefunctionPlot[\*StyleBox[\"wavefunctions, names\",\"TI\"]] plots the structure of  \*StyleBox[\"wavefunctions\",\"TI\"] (coefficients), using information about the \*StyleBox[\"names\",\"TI\"] of the orbitals in the basis."

 GTTbOrbitalsFromBasis     ::usage = "GTTbOrbitalsFromBasis[bas_] gives a list of orbital names from a  standard basis description."
 
(*----------- NanoTubes ------------------------------------------------------*)
GTTbTubeBands              ::usage  = "GTTubebands[\*StyleBox[\"Hamiltonian, n,m,kpoints\",\"TI\"]] calculates the band structure of an (n,m) nanotube."

(*----------- Connection to Wannier90 package --------------------------------*)
GTTbReadWannier90          ::usage  = "GTTbReadWannier90[\*StyleBox[\"file name\",\"TI\"]] reads information from Wannier90 package."
GTTbWannier90Hamiltonian   ::usage  =  "GTTbWannier90Hamiltonian[\*StyleBox[\"data,basis\",\"TI\"]] constructs Hamiltonian from Wanier90 data."

(*---------- internal modules  -----------------------------------------------*)
(* GTTbSimplify            ::usage = "GTTbSimplify[mat,ll,lu,dist,rl] simplifies matrix elements with rule rl" *)
(* GTTbParameterSet        ::usage = "GTTbParameterSet[lmax,dmax] is an input module to generate two-center parameter sets" *)
(* GTTbParmToFortran       ::usage = "GTTbParamToFortran[set] constructes a substitution rule for the symbols of the hopping parameters for FORTRAN output"*)

(*----------- Options ---------------------------------------------------------*)
 Options[GTPlotStateWeights]       = {GOColorScheme -> "ElementData", GOPlot -> True,Boxed->True}
 Options[GTSymmetryBasisFunctions] = {GONames->{}, GOVerbose->True}
 Options[GTTbAtomicWaveFunction]   = {GOHarmonics->"Complex"}
 Options[GTTbBlochFunction]        = {GOHarmonics->"Complex"}
 Options[GTTbHamiltonianElement]   = {GOTbBasis -> 0, GOTbRealSpace -> False}
 Options[GTTbHamiltonianRS]        = {GOTbBasis -> 0, GOVerbose -> False}
 Options[GTTbMatrixElement]        = {GOTbRule->1,GOTbBasis->0}
 Options[GTTbNumberOfIntegrals]    = {GOVerbose -> False,GONames->{}}
 Options[GTTbParameterNames]       = {GOTbBasis->0,GOVerbose->False}
 Options[GTTbParameterSet]         = {GOTbBasis->0}
 Options[GTTbRealSpaceMatrix]      = {GOTbBasis->0}
 Options[GTTbSimplify]             = {GOTbBasis->0}
 Options[GTTbSymbol2C]             = {GOTbBasis->0}
 Options[GTTbSymmetrySingleBand]   = {GOVerbose -> False, GOSpinOrbit -> False}
 Options[GTTbSymmetryBands]        = {GOVerbose -> False, GOIrepNotation -> "Bouckaert", PlotStyle -> {{Thin, Black}}, Joined -> True, PlotRange -> All, GOLabelShift -> {0.05, 0.02},GOLabelStyle->{}, 
 	                                  GOPlot -> True, FrameLabel -> {"", "Energy (eV)"}, PlotLabel -> "Band structure", GOSpinOrbit -> False, GOOrbitalConstruction -> True,GOPrecision->4}
 Options[GTTbSymmetryPoint]        = {GOVerbose -> False, GOIrepNotation -> "Bouckaert", GOSpinOrbit -> False}
 Options[GTTbToFortran]            = {GOTbBasis->0,GOVerbose->False}
 Options[GTTbRSToFortran]          = {GOTbBasis -> 0, GOVerbose -> False,GOImpurity->False}
 Options[GTTbToFortranList]        = {GOVerbose -> False}
 Options[GTIrepInfo]               = {GOMethod->"Cornwell"} 
 Options[GTTbTubeBands]            = {GOVerbose->True,GOPlotBands->True,GOTbOrthogonal->True,GOStore->0}
 Options[GTTbHamiltonianOnsite]    = {GOVerbose -> False}
 Options[GTTbReadWannier90]        = {GOVerbose -> False}
 Options[GTTbWannier90Hamiltonian] = {GOCutOff  -> 0,GOVerbose->False}
 Options[GTTbMatrixElement3C]      = {GOVerbose -> False}
 Options[GTTbIntegralRules]        = {GOVerbose -> False}
 Options[GTTbIntegralRelations]    = {GOVerbose -> False}
 Options[GTTbSpinOrbit]            = {GOVerbose -> False}
 Options[GTHamiltonianPlot]        = {ItemSize-> All}
 Options[GTHamiltonianList]        = {GOList -> "All"}
 Options[GTFindStateNumbers]       = {GODecimals -> 0}
 
  Begin["`Private`"] (* Begin Private Context *) 

(****t* /GTWaveFunctionPlot
! NAME
!  GTTHamiltonianList
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  7.2.2017 : first version
! USAGE
!  GTHamiltonianList[ham_,names_] list the elements of an hamiltonian using names
! INPUT
!  a list that describes a basis
! OUTPUT
!  list of orbitals
! ERROR MESSAGES
!  
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  This is a part, what is used originally in GTHamiltonianplot. The list of orbital names
!  may be generated separately by this command.
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

 GTWaveFunctionPlot[wav_, bas_] := Module[
    {bez,nwf,ncof,hamp,i,j,blist,test},
     bez  = bas;
     nwf  = Length[wav] + 1;
     ncof = Length[wav[[1]]] + 1;
     hamp = Table[" ", {nwf}, {ncof}];
     If[bez == {},
        bez = Table[ToString[i], {i, 1, ncof - 1}],
        None
     ];
     hamp[[1, 2 ;; ncof]] = bez;
     hamp[[2 ;; nwf, 1]] = Table[ToString[i], {i, 1, nwf - 1}];
     blist = {};
     Do[
        Do[
           test = Abs[wav[[i, j]]];
           If[test > 0,
              blist = Append[blist, {{i + 1, j + 1} -> GTBackGroundColor1}],
              None
           ]
        , {j, 1, ncof - 1}]
     , {i, 1, nwf - 1}];
     blist = Flatten[blist, 1];
     tab=Grid[hamp, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
                Background -> {None, None, blist}, ItemSize -> All
             ];
     Print[tab]
 ]

(*
***)

(****t* /GTHamiltonianList
! NAME
!  GTHamiltonianList
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 07.02.2017 : first version
!  * 21.02.2017 : option GOList implemented 
!  * 25.06.2018 : check header and documentation
!  * 20.09.2022 : line and row of the matrix element  aew additionale indicated
!
! USAGE
!  GTHamiltonianList[ham_,names_] list the elements of a hamiltonian using names
! INPUT
!  o ham   - the Hamiltonian
!  o names - list of names created by means of GTTbOrbitalsFromBasis
! OUTPUT
!  Table of matrix elements
! GTPack OPTIONS
!  * GOList:
!
!    - "All"  - complete HAmiltonina
!    - "NonZero" - only nonzero elements
!    - "Lower" - lower part of the Hermitian matrix (Slater&Koster) style
!  
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The commands creates at the end of the day a list which is similar to that in the seminal paper of 
!  Slater&Koster.
! LITERATURE
!  Slater, Koster, Phys. Rev.  94, 1498 (1954)
! TODO
!  -
! RELEASE
!  1.0.0 
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTHamiltonianList[ham_, names_, OptionsPattern[]] :=
  Module[{verb, dim1, dim2, tab, ni, nj, i, j, bez, test}, 
  verb = OptionValue[GOList];
  If[Intersection[{"All", "NonZero", "Lower"}, {verb}] == {}, 
   Print["Error: Method not implemented"]; Abort[], None];
  dim1 = Length[ham];
  dim2 = Length[ham[[1]]];
  tab = {{"i", "j", "Element", " Analytic expression of Hamiltonian"}};
  If[verb == "Lower", Do[Do[If[names == {}, ni = i;
      nj = j;
      bez = "(" <> ToString[nj] <> "/" <> ToString[ni] <> ")", 
      ni = names[[1, i]];
      nj = names[[2, j]];
      bez = "(" <> nj <> "/" <> ni <> ")"];
     tab = Append[tab, {i, j, bez, ham[[j, i]]}], {i, j, dim1}], {j, 
     1, dim2}], None];
  If[verb == "NonZero" || verb == "All", 
   Do[Do[If[names == {}, ni = i;
      nj = j;
      bez = "(" <> ToString[nj] <> "/" <> ToString[ni] <> ")", 
      ni = names[[1, i]];
      nj = names[[2, j]];
      bez = "(" <> ni <> "/" <> nj <> ")"];
     test = ham[[i, j]];
     If[(Head[test] === Integer || Head[test] === Real || 
         Head[test] === Complex) && verb == "NonZero", None, 
      tab = Append[tab, {i, j, bez, test}]], {i, 1, dim1}], {j, 1, 
     dim2}], None];
  Print[Grid[tab, Frame -> All, 
    Background -> {{1 -> Yellow, 2 -> Yellow, 
       3 -> GTBackGroundColor1}, 
      1 -> GTBackGroundColor1, {1, 1} -> GTCornerColor}]]]
      
     
(*  GTHamiltonianList[ham_, names_, OptionsPattern[]] := Module[
  {verb,dim1,dim2,tab,ni,nj,i,j,bez,test},
   verb = OptionValue[GOList];
   If[Intersection[{"All", "NonZero", "Lower"}, {verb}] == {},
      Print["Error: Method not implemented"]; Abort[],
      None
   ];
   dim1 = Length[ham];
   dim2 = Length[ham[[1]]];
   tab  = {{"Element", " Analytic expression of Hamiltonian"}};
   If[verb == "Lower",
      Do[
      	 Do[
      	    If[names == {}, 
      	       ni  = i; 
      	       nj  = j;
               bez = "(" <> ToString[nj] <> "/" <> ToString[ni] <> ")", 
               ni  = names[[1, i]]; 
               nj  = names[[2, j]];
               bez = "(" <> nj <> "/" <> ni <> ")"
            ];
            tab = Append[tab, {bez, ham[[j, i]]}]
         , {i, j, dim1}]
      , {j, 1, dim2}],
      None
   ];
   If[verb == "NonZero" || verb == "All",
      Do[
         Do[
            If[names == {},
               ni  = i; 
               nj  = j;
               bez = "(" <> ToString[nj] <> "/" <> ToString[ni] <> ")",
               ni  = names[[1, i]]; 
               nj  = names[[2, j]];
               bez = "(" <> ni <> "/" <> nj <> ")"
            ];
            test = ham[[i, j]];
            If[(Head[test] === Integer||Head[test] === Real||Head[test] ===Complex) && verb == "NonZero",
               None,
               tab = Append[tab, {bez, test}]
            ]
         , {i, 1, dim1}]
      , {j, 1, dim2}],
      None
   ];
   Print[Grid[tab, Frame -> All, Background -> {1 -> GTBackGroundColor1, 
                   1 -> GTBackGroundColor1, {1, 1} -> GTCornerColor}
             ]
        ]
]
*)

(*
 ***)  
  
(****t* /GTTbOrbitalsFromBasis
! NAME
!  GTTborbitalsFromBasis
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 07.02.2017 : first version
!  * 25.06.2018 : check header and documentation
! USAGE
!  GTOrbitalsFromBasis[bas_] creates a list of orbital names from a basis description.
! INPUT
!  a list that describes a basis
! OUTPUT
!  list of orbitals
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  This is a part, what is used originally in GTHamiltonianplot. The list of orbital names
!  may be generated separately by this command.
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)

 GTTbOrbitalsFromBasis[bas_] := Module[
 	{nb,mom,defbez,bez,tbas,mult,angmom,nat,i,j,ang,la,l,pos,k,ll,s},
     nb     = Length[bas];
     mom    = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4};
     defbez = {"s", "\!\(\*SubscriptBox[\(p\), \(y\)]\)", 
               "\!\(\*SubscriptBox[\(p\), \(z\)]\)", 
               "\!\(\*SubscriptBox[\(p\), \(x\)]\)", 
               "\!\(\*SubscriptBox[\(d\), \(xy\)]\)", 
               "\!\(\*SubscriptBox[\(d\), \(yz\)]\)", 
               "\!\(\*SubscriptBox[\(d\), \(3 \*SuperscriptBox[\(z\),            \
                 \(2\)]\)]\)", "\!\(\*SubscriptBox[\(d\), \(xz\)]\)", 
               "\!\(\*SubscriptBox[\(d\), \(\*SuperscriptBox[\(x\), \(2\)] -     \
                 \*SuperscriptBox[\(y\), \(2\)]\)]\)"
              };
     bez    = {};
     tbas   = Transpose[bas];
     mult   = tbas[[2]]; 
     angmom = {};
     nat    = Plus @@ mult;
     If[nat == nb, 
     	angmom = tbas[[3]], 
        Do[
           Do[
           	  angmom = Append[angmom, tbas[[3, i]]]
           , {j, 1, mult[[i]]}]
        , {i, 1, nb}]
     ];
     Do[
        ang = angmom[[k]] /. mom; 
        la = Length[ang];
        Do[
           l   = ang[[ll]];
           pos = Sum[2*s + 1, {s, 0, l - 1}];
           If[ang[[ll]] == -1, 
              bez = Append[bez, "\!\(\*SuperscriptBox[\(s\), \(\[Star]\)]\)"], 
              None
           ];
           Do[
              bez = Append[bez, defbez[[pos + m]]]
           , {m, 1, 2 l + 1}]
        , {ll, 1,la}]
     , {k, 1, nat}];
     Return[bez]
 ]
  
 (*
 ***) 


(****t* /GTTbIntegralRelations
! NAME
!  GTTbIntegratalRelations
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  05.01.2017   : first version
! USAGE
!  GTTbIntegrationRules[ireps1, ireps2] finds the rules to minimize the number of parameters
! INPUT 
!   * ireps1  - list of representation matrices to the generators ogf Glp for irep of basis function1
!   * ireps2  - same for basis function 2
!
! OPTIONS
!  -
! OUTPUT
!  part of matrix element for a certain shell 
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  If large eigenvalue problems are under consideration, it is sometimes necessary to find the number of eigenstates in a certain
!  energy region for a detailed inspection of the wave function.
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

GTTbIntegralRelations[rep1_, rep2_, OptionsPattern[]] := Module[
	{verb, lr1, lr2, facs,leq,k,ir1,ir2},
  (*--- options ---*)
    verb = OptionValue[GOVerbose];
  (*--- ---*)
  lr1  = Length[rep1]; 
  lr2  = Length[rep2]; leq = lr1*lr2;
  facs = Table[0, {leq + 2}, {leq + 2}];
  facs[[1, 1]] = facs[[1, 2]] = facs[[2, 1]] = facs[[2, 2]] = " ";
  k    = 0;
  Do[
     Do[
        k                = k + 1;
        facs[[2, k + 2]] = facs[[k + 2, 2]] = ir1;
        facs[[1, k + 2]] = facs[[k + 2, 1]] = ir2
      , {ir1, 1, lr1}]
   , {ir2, 1, lr2}];
   (*--- coefficients from representations ---*)
  k = 0;
  Do[
       Do[k = k + 1; l = 0;
            Do[
               Do[
                     l = l + 1;
                     
      facs[[k + 2, l + 2]] = rep1[[nu, mu]] rep2[[nup, mup]]
        , {nu, 1, lr1}]
     , {nup, 1, lr2}]
      , {mu, 1, lr1}]
   , {mup, 1, lr2}];
  (*--- print of matrix ---*)
  If[verb,
       Print[Grid[facs, Frame -> True,
                              
     Dividers -> {{3 -> GTDividerColor1}, {3 -> GTDividerColor1}},
                              
     Background -> {{1 -> GTBackGroundColor1, 
        2 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1, 
        2 -> GTBackGroundColor1},
       {{1, 2}, {1, 2}} -> GTCornerColor}]],
   None];
  (*--- list of variables ---*)
  
  vars = Table[
    Symbol["I" <> ToString[facs[[1, i + 2]]] <> 
      ToString[facs[[2, i + 2]]]], {i, 1, leq}];
  (*--- set up the equations ---*)
  eqs = {}; varsi = {}; 
  varse = vars;
  Do[
   s = 0;
   Do[
        s = s + vars[[i]]*facs[[j + 2, i + 2]]
    , {i, 1, leq}];
   eqt = vars[[j]] == s;
   If[Head[eqt] === Symbol,
        varsi = Append[varsi, vars[[j]]];
        varse = Complement[varse, {vars[[j]]}],
         eqs = Append[eqs, eqt]
    ]
   , {j, 1, leq}];
  Return[{eqs, varsi, varse}]
  ]
  
 (*
 ***) 

(****t* /GTTbIntegralRules
! NAME
!  GTTbIntegrationRules
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  05.01.2017   : first version
! USAGE
!  GTTbIntegrationRules[ireps1, ireps2] finds the rules to minimize the number of parameters
! INPUT 
!   * ireps1  - list of representation matrices to the generators ogf Glp for irep of basis function1
!   * ireps2  - same for basis function 2
!
! OPTIONS
!  -
! OUTPUT
!  part of matrix element for a certain shell 
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  If large eigenvalue problems are under consideration, it is sometimes necessary to find the number of eigenstates in a certain
!  energy region for a detailed inspection of the wave function.
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

 GTTbIntegralRules[ireps1_, ireps2_, OptionsPattern[]] := Module[
   {verb, nrep, lr2, eqs, varsi, varse, eqt, vi, ve, irep, rule, reps1, reps2},
  (*--- options ---*)
    verb = OptionValue[GOVerbose];
  (*--- set initila data ---*)  
    nrep  = Length[ireps1];
    lr2   = Length[ireps2];
    reps1 = Transpose[ireps1][[2]];
    reps2 = Transpose[ireps2][[2]];
    If[nrep == lr2,
       None,
       Print["Error: same number of Ireps have to be provided!"];  Abort[]
    ];
    eqs   = {}; 
    varsi = {}; 
    varse = {};
    Do[
       {eqt, vi, ve} = GTTbIntegralRelations[reps1[[irep]], reps2[[irep]], GOVerbose -> verb];
        eqs   = Append[eqs, eqt];
        varsi = Append[varsi, vi];
        varse = Append[varse, ve];
    , {irep, 1, nrep}];
  (*--- Solve equations , get rules ---*)
    eqs   = Flatten[eqs] // Union;
    varse = Flatten[varse] // Union;
    varsi = Flatten[varsi] // Union;
    rule  = Solve[eqs, varse];
  Return[{rule, varsi}]
  ]

(*
***)


(****t* /GTTbMatrixElement3C
! NAME
!  GTTbMatrxiElement3C
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  05.01.2017   : first version
! USAGE
!  GTTbMatrixElement[repmat1,repmat2,mu,mup,vars,qlp,kv] part of matrix element from a certain shell
! INPUT 
!   * repmat1   - list of symmetry elements and rep matrices {{sym element1, matrix 1},{sym element2, matrix 2}..} 
!   * repmat2   - similar to repmat1 for another irreducible representation
!   * mu        - row to rep1
!   * mup       - row to rep2
!   * vars      - list is output of GTVarList
!   * qlp       - the vector Q_p^l
!   * kv        - k-vector, usually {xi,eta,zeta}
! OPTIONS
!  -
! OUTPUT
!  part of matrix element for a certain shell 
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  If large eigenvalue problems are under consideration, it is sometimes necessary to find the number of eigenstates in a certain
!  energy region for a detailed inspection of the wave function.
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

(*
GTTbMatrixElement3C[rep1_, rep2_, mu_, mup_, vars_, qlp_, kv_, OptionsPattern[]] := Module[
	{s1,s2,irep1,irep2,verb,sym,ns,kvec,lp,s,su,i,nu,nup,arg}, 
    {s1, irep1} = Transpose[rep1];
    {s2, irep2} = Transpose[rep2];
  (*---options---*)
    verb = OptionValue[GOVerbose];
  (*---tests---*)
    If[Complement[s1, s2] == {},
       sym = s1;
       ns  = Length[sym];
       If[verb,
          Print["Symmetry operations involved :", sym], 
          None
       ], 
       Print["Error: Inconsistency in input of Ireps!"]; Abort[]
    ];
    If[Length[#] & /@ qlp // Union == {1}, 
       None, 
       Print["Error: Not formulated for more than one Qlp yet!"]; Abort[]
    ];
    If[ns == Length[irep1] && ns == Length[irep2], 
       None, 
       Print["Error: Check input of Ireps!"]; Abort[]
    ];
  (*---calculate matrix element---*)
    s    = 0;
    kvec = kv;
    lp   = Length[irep1[[1]]]*Length[irep2[[1]]];
    Do[
       If[Head[vars[[2, i]]] === Symbol, 
       	  su = 0;
          Do[
          	 nu  = vars[[3, i]];
             nup = vars[[4, i]];
             arg = GTGetMatrix[sym[[j]]].kvec;
             su  = su + Exp[I arg.qlp]*Conjugate[irep1[[j, nu, mu]]]*irep2[[j, nup, mup]]
          , {j, 1, ns}];
          s = s + vars[[2, i]]*su, 
          None
       ]
    , {i, 1, lp}];
    s = s // ExpToTrig // TrigFactor // FullSimplify;
    Return[s]   
 ]
*)

GTTbMatrixElement3C[rep1_, rep2_, mu_, mup_, vars_, qlp_, kv_, OptionsPattern[]] := Module[
  {irep1, irep2, s1, s2, sym, lp, ns, s, su, kvec, nu, nup, arg, j, i,verb},
    {s1, irep1} = Transpose[rep1];
    {s2, irep2} = Transpose[rep2];
  (*--- options---*)
    verb = OptionValue[GOVerbose];
  (*--- tests ---*)
    If[Complement[s1, s2] == {},
       sym = s1; 
       ns  = Length[sym];
       If[verb,
          Print["Symmetry operations involved :", sym],
          None
       ],
       Print["Error: Inconsistency in input of Ireps!"]; Abort[]
    ];
    If[Length[#] & /@ qlp // Union=={1},
       None,
       Print["Error: Not formulated for more than one Qlp yet!"]; Abort[]
    ];
    If[ns == Length[irep1] && ns == Length[irep2],
       None,
       Print["Error: Check input of Ireps!"]; Abort[]
    ];
  (*--- calculate matrix element ---*)
    s    = 0; 
    kvec = kv;
    lp   = Length[irep1[[1]]]*Length[irep2[[1]]];
    Do[
       If[Head[vars[[2, i]]] === Symbol,
          su = 0;
          Do[
             nu  = vars[[3, i]];
             nup = vars[[4, i]];
             arg = GTGetMatrix[sym[[j]]].kvec;
             su = su + Exp[I arg.qlp]*Conjugate[irep1[[j, nu, mu]]]*irep2[[j, nup, mup]]
          , {j, 1, ns}];
          s = s + vars[[2, i]]*su,
          None
       ]
    , {i, 1, lp}];
    s = s // ExpToTrig // TrigFactor // FullSimplify;
    Return[s]
  ]

  
  (*
  ***)

(****t* /GTTbSymbol3C
! NAME
!  GTTBSymbol3C
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 05.01.2017 : first version
!  * 27.06.2018 : check header create documentation page 
! USAGE
!  GTTBSymbol3C[ir1_,row_1,ir_2,row_2,qlp] changes a formal name for the three-center parameter in the Slater&Koster form.
! INPUT
!   * ir1   - Irep_1
!   * row1  - row to Irep_1
!   * ir2   - Irep_2
!   * row2  - row to Irep2
!   * qlp   - vector Q_p^l
! OPTIONS
!  -
! OUTPUT
!  name of energy integral and substitution rule
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  In the first step of construction of the matrix elements formal names are used. In a last step those names are substituted by 
!  names which are used also in the slater-Koster approach. 
! LITERATURE
!  Egorov
! TODO
!  In the moment it is no variable. It is fixed to cubic structures. Perhaps it doesnt work well for the hcp example in Egorov's paper.
!  A genralization is necessary.
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)


GTTbSymbol3C[j_, mu_, jp_, mup_, qlp_] := Module[
	{p, d1, d2, s1, s2, sym, sym0, su, str},
     p  = {"x", "z", "y"}; 
     d1 = {"\!\(\*SuperscriptBox[\(x\), \
\(2\)]\)-\!\(\*SuperscriptBox[\(y\), \(2\)]\)", 
  "3\!\(\*SuperscriptBox[\(z\), \(2\)]\)-\!\(\*SuperscriptBox[\(r\), \
\(2\)]\)"};
     d2 = {"xz", "xy", "yz"};
     If[j == 1,
        s1 = "s",
        None
     ];
     If[j == 9,
       s1 = p[[mu]],
       None
     ];
     If[j == 6,
        s1 = d1[[mu]],
        None
     ];
     If[j == 8,
       s1 = d2[[mu]],
       None
     ];
     If[jp == 1,
        s2 = "s",
        None
     ];
     If[jp == 9,
       s2 = p[[mup]],
       None
     ];
     If[jp == 6,
       s2 = d1[[mup]],
       None
     ];
     If[jp == 8,
        s2 = d2[[mup]],
        None
     ];
     str  = "E(" <> ToString[qlp[[1]]] <> ToString[qlp[[2]]] <> ToString[qlp[[3]]] <> ")";
     su   = s1 <> "," <> s2;
     sym  = Subscript[str, su];
     sym0 = Symbol["I" <> ToString[mup] <> ToString[mu]];
     Return[{sym, sym0 -> sym}]
 ]

(*
***)

(****t* /GTTbVarList
! NAME
!  GTTbVarList
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  31.01.2017    : first version
!  23.02.2017    : j,j' are not necessary here, taken away from argment list, external function now
! USAGE
!  GTTbarlist[] findes the number of states in the list eigenvalues lying in the energy interval [emin,emax].
! INPUT
!   * eigenvalues - vector of eigenenergies
!   * emin        - energy minimum
!   * emax        - energy maximum
! OPTIONS
!  -
! OUTPUT
!  table of state numbers and energies
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  If large eigenvalue problems are under consideration, it is sometimes necessary to find the number of eigenstates in a certain
!  energy region for a detailed inspection of the wave function.
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

 GTTbVarList[mu_, mup_, rule_] := Module[
 	{m1, m2,  v1, k},
     v1 = Table[0, {4}, {mu*mup}];
     k  = 0;
     Do[
        Do[k = k + 1;          
           v1[[1, k]] = 
           v1[[2, k]] = Symbol["I" <> ToString[m2] <> ToString[m1]];
           v1[[3, k]] = m2;
           v1[[4, k]] = m1
        , {m2, 1, mup}]
     , {m1, 1, mu}];
     v1[[2]] = v1[[2]] /. Flatten[rule];
     Return[v1]
  ]

(* 23.2 Auskommentiert und Variante aus Notebook uebernommen
 GTVarList[j_, mu_, jp_, mup_, rule_] := Module[
 	 {v1,k,m1,m2},
     v1 = Table[0, {4}, {mu*mup}];
     k  = 0;
     Do[
        Do[
           k          = k + 1;
           v1[[1, k]] = 
           v1[[2, k]] = Symbol["I" <> ToString[m2] <> ToString[m1]];
           v1[[3, k]] = m2;
           v1[[4, k]] = m1
        , {m2, 1, mup}]
     , {m1, 1, mu}];
     v1[[2]] = v1[[2]] /. Flatten[rule];
     Return[v1]
  ]
*)

  
(*
***)  

(****t* /GTFindStateNumbers
! NAME
!  GTFindStateNumbers
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 15.07.2014 : first version
!  * 10.12.2016 : check if emin<emax, otherwise interchange emin and emax.
!  * 25.06.2018 : check header and documentation
!  * 13.02.2023 : improved version with better table and GODecimals
! USAGE
!  GTFindStateNumbers[eigenvalues,emin,emax] findes the number of states in the list eigenvalues lying in the energy interval [emin,emax].
! INPUT
!   * eigenvalues	-	vector of eigenenergies
!   * emin			-	energy minimum
!   * emax			-   energy maximum
! OUTPUT
!  table of state numbers and energies
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  If large eigenvalue problems are under consideration, it is sometimes necessary to find the number of eigenstates in a certain
!  energy region for a detailed inspection of the wave function.
! LITERATURE
!  -
! TODO
! -
! RELEASE
!  1.0.0
! PROBLEMS
! -
! 
!--------------------------------------------------------------------------------
! SOURCE
*)


GTFindStateNumbers[energ_, emin_, emax_, OptionsPattern[]] := Module[
	{nround,ln,pos,tit,emi,ema,en,ten,lnp,m,posn,pos1},
     nround = OptionValue[GODecimals];
     ln     = Length[energ];
     pos    = {}; 
     tit    = {"#", "state number", "energy"};
     If[emin <= emax, 
     	emi = emin; ema = emax, 
  	    emi = emax, ema = emin
     ];
     Do[en = energ[[i]];
        If[en > emi && en < ema,
           If[nround == 0,
              ten = en,  
              ten = ToString[PaddedForm[N[en], {nround, nround}]]
           ];
           pos = Append[pos, {i, ten}], 
           None
        ]
     , {i, 1, ln}];
     lnp  = Length[pos]; 
     m    = Floor[lnp/10];
     posn = Transpose[pos]; 
     posn = Prepend[posn, Table[i, {i, 1, lnp}]] // Transpose;
     Do[pos1 = posn[[(i - 1)*10 + 1 ;; i*10]];
        pos1 = Prepend[pos1, tit];
        Print[Grid[pos1 // Transpose, Frame -> All, Background -> 
        	       {None, {1 -> LightBlue, 2 -> GTBackGroundColor1}}]]
     , {i, 1, m}];
     If[10 m < lnp,
        pos1 = posn[[10 m + 1 ;; lnp]]; pos1 = Prepend[pos1, tit]; 
        Print[Grid[pos1 // Transpose, Frame -> All, Background ->
                    {None, {1 -> LightBlue,2 -> GTBackGroundColor1}}]],
        None
     ];
     Return[]
]

(*
***)

(****t* /GTPlotStateWeights
! NAME
!  GTPlotStateWeights
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBnding.m
! MODIFICATION HISTORY
!  * 20.07.2014 : first version
!  * 25.06.2018 : check header and documentation
! USAGE
!  GTPlotStateWeights[cluster,distance,basis,wave function,scale] demonstrates at which atoms of a cluster, 
!  the wave function (a solution of a real space tight-binding Hamiltonian) has the heighest weights.
! INPUT
!  * cluster       - cluster
!  * distance      - minimal distance between atoms (bond length)
!  * basis         - standard description of the basis
!  * wave function - eigenvector
!  * scale         - scale factor: radius of the sphere with maximal weight
 OUTPUT
!  graphical presentation of the  weights of the selected states
! GTPack OPTIONS
!  GOColorScheme 
!    o "ElementData" (standard) -> Colors are taken from ElementData
!    o list of colors, Length of list has to correspond to the length of basis
!
!  GOPlot      
!    o True (standard) - graph is plotted
!    o False           - graphics object is given {gstruc,gwf}
!                        gstruc - structural graph
!                        gwf - plot of wavefunction weights
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

GTPlotStateWeights[cl_, dist_, basis_, ev_, scale_, OptionsPattern[]] :=Module[
	{mom,vec,opc,rule,basisat,bas,bas1,basisam,i,aml,nom,lc,glines,gsph,j,no,col,gstruc,am,k,su,sc,gwf,box},
     mom = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4};
     opc = OptionValue[GOColorScheme];
     box = OptionValue[Boxed];
     (*--- discuss basis ---*)
      bas     = Transpose[basis]; 
      bas1    = Transpose[{bas[[1]], bas[[3]]}];
      basisat = Transpose[bas1][[1]];
      basisam = Transpose[bas1][[2]];
      If[opc === "ElementData", 
       	 None,
       	 If[Length[basisat]==Length[opc],
       	 	None,
       	 	Print["Error GTPlotStateWeights : Number of atoms in basis and number of colors do not correspond"];Abort[]
       	 ];	
         rule = {}; 
         Do[
           rule = Append[rule, basisat[[i]] -> opc[[i]]]
         , {i, 1, Length[basisat]}];
         aml = Transpose[cl][[2]] /. rule
      ];
      (*--- generate structural graph ---*)
      lc = Length[cl]; 
      glines = {}; 
      nom = Map[Norm[#] &, Transpose[cl][[1]]]; nom = Max[nom]; 
      gsph = {};
      Do[
         Do[
            no = Norm[cl[[i, 1]] - cl[[j, 1]]];
            If[no == dist, 
               glines = Append[glines, Line[{cl[[i, 1]], cl[[j, 1]]}]], 
               None
            ]
         , {j, i + 1, lc}];
         If[opc === "ElementData", 
            col = ElementData[cl[[i, 2]], "IconColor"],  
            col = aml[[i]]
         ];
         gsph = Append[gsph, {col, Sphere[cl[[i, 1]], 0.03]}]
      , {i, 1, lc}];
      gstruc = {gsph, glines, {Opacity[0.4], Sphere[{0, 0, 0}, nom]}};
      (*--- construct the plot of weights ---*)
      am = Map[Plus @@ # &, Map[(2 # + 1) &, basisam /. mom /. {-1 -> 0}]];
      rule = {}; 
      Do[
      	 rule = Append[rule, basisat[[i]] -> am[[i]]]
      , {i, 1, Length[basisat]}];
      aml = Transpose[cl][[2]] /. rule; 
      vec = Table[0, {lc}]; k = 0;
      Do[
         su = 0;
         Do[
         	k  = k + 1;
            su = su + ev[[k]]^2
         , {j, 1, aml[[i]]}];
         vec[[i]] = su
      , {i, 1, lc}];
      sc = scale/Max[vec]; gsph = {};
      Do[
         gsph = Append[gsph, Sphere[cl[[i, 1]], vec[[i]]*sc]]
      , {i, 1, lc}];
      gwf = {gsph, glines}; 
      If[OptionValue[GOPlot], 
      	 Return[Graphics3D[{{gstruc, gwf}},Boxed->box]],
         Return[{{gstruc, gwf}}]
      ]
  ]


(*
***)



(****t* /GTTbSpinMatrix
! NAME
!  GTTbSpinMatrix
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 03.07.2013 : first version
!  * 25.06.2018 :check header and documnetation
! USAGE
!  GTTbSpinMatrix[angular momentum,s1,s2] gives elementary spin matrices for tight-binding Hamiltonians.
! INPUT
!  * angular momentum l
!  * s1 - spin direction bra vector (+1,-1)
!  * s2 - spin direction ket vector (+1,-1)
! OUTPUT
!  spin matrix of dimension 2l+1
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTTbAuxVogl
! GTPack NOTEBOOKS 
!  GTPack_Notebooks/TightBinding/GT_Test_Vogl.nb for some tests
! DESCRIPTION
!  -
! LITERATURE
!  Podolskiy, Vogl, Phys. Rev. B 69, 233101 (2004)
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTbSpinMatrix[l_,s1_,s2_]:=Module[{\[Mu]},hso=Table[0,{2*l+1},{2*l+1}];
                If[Sign[s1]==Sign[s2],
                	Do[
                		Do[hso[[m1+l+1,m2+l+1]]=Sign[s1]*I*m2*KroneckerDelta[m1,-m2]/2,
                        {m1,-l,l}],
                   	{m2,-l,l}],
                   	Do[
                   		Do[\[Mu]=Conjugate[Tb\[Theta][m1]]*Tb\[Theta][m2](1+Tb\[Tau][-Abs[m1*m2]]);
                   			hso[[m1+l+1,m2+l+1]]= \[Mu](KroneckerDelta[Abs[m1],Abs[m2]+Sign[s1]]-(-1)^(Tb\[Tau][m1]+
                   			Tb\[Tau][m2])KroneckerDelta[Abs[m1],Abs[m2]-Sign[s1]])*Sqrt[(l+m1^2-Abs[m1*m2])(l+m2^2-Abs[m1*m2])]/2,
                   		{m1,-l,l}],
                   	{m2,-l,l}]];
Return[hso]
]
(*
***) 

(****t* /GTTbAuxVogl
! NAME
!  GTTbAuxVogl
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  July 2013 : first version
! USAGE
!  Collection of auxilary Functions for matrix element construction
! INPUT
!  -
! OUTPUT
!  -
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  GTPack_Notebooks/TightBinding/GT_Test_Vogl.nb for some tests
! DESCRIPTION
!  The auxilary functions are used for the construction of the TB matrix elements
!  in two-center form according to Podolskiy and Vogl:
!   * d-function   eq. (7)
!   * tau          eq. (15)
!   * a,b          eq. (21)
!   * s            eq. (24)
!   * t            eq. (25)
!   * theta        eq. (32)
! LITERATURE
!  Podolskiy, Vogl, Phys. Rev. B 69, 233101 (2004)
! TODO
!  -
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)

Tb\[Tau][m_]:=If[m>=0,1,0]

Tb\[Theta][m_]:=(1-KroneckerDelta[m,0])(Tb\[Tau][m]+I*Tb\[Tau][-m])/Sqrt[2] +KroneckerDelta[m,0]/2

Tbd[l_,m_,mp_]:=Module[{\[Beta],df,df1},df=WignerD[{l,m,mp},0,\[Beta],0]//TrigReduce//TrigExpand;
df1=df/.{Cos[\[Beta]]->DN,Sin[\[Beta]]->Sqrt[1-DN^2]}  //FullSimplify//Factor//FullSimplify;
df1=df1/.{ Sqrt[1-DN^2]->Sqrt[1-DN]Sqrt[1+DN]};Return[df1]]

TbSi[n_]:=Module[{upper,si,s,c},upper=IntegerPart[(n-1)/2];s=DM/Sqrt[1-DN^2];c=-DL/Sqrt[1-DN^2];
	              si=Sum[(-1)^k Binomial[n,2*k+1]s^(2*k+1)c^(n-2*k-1),{k,0,upper}];FullSimplify[si,{DL>0,DN>0}]
                ]

TbCo[n_]:=Module[{upper,co,s,c},upper=IntegerPart[n/2];s=DM/Sqrt[1-DN^2];c=-DL/Sqrt[1-DN^2];
	             co=Sum[(-1)^k Binomial[n,2*k]s^(2*k)c^(n-2*k),{k,0,upper}];FullSimplify[co,{DL>0,DN>0}]
                ]

TbA[m_]:=Module[{},If[m==0,1/Sqrt[2],FullSimplify[(-1)^Abs[m](Tb\[Tau][m]TbCo[Abs[m]]-Tb\[Tau][-m]TbSi[Abs[m]]),{DL>0,DM>0}]]]

TbB[m_]:=Module[{},FullSimplify[(-1)^Abs[m](Tb\[Tau][m]TbSi[Abs[m]]+Tb\[Tau][-m]TbCo[Abs[m]]),{DL>0,DM>0}]]

TbS[l_,m_,mp_]:=Module[{c},c=TbA[m]((-1)^Abs[mp]Tbd[l,Abs[m],Abs[mp]]+Tbd[l,Abs[m],-Abs[mp]]);Return[c//FullSimplify]]

TbT[l_,m_,mp_]:=Module[{c},c=(1-KroneckerDelta[m,0])TbB[m]((-1)^Abs[mp]Tbd[l,Abs[m],Abs[mp]]-Tbd[l,Abs[m],-Abs[mp]]);Return[c//FullSimplify]]

(*
***) 


(****t* /GTTbMatrixElement
! NAME
!  GTTbMatrixElement
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 15.07.2013 : first version
!  * 14.04.2014 : extension to s* orbital
!  * 25.06.2018 : check header and description
!  * 17.07.2018 : bug fix equal atoms
! USAGE
!  GTTbMatrixElement[l1,m1,l2,m2,shell] gives the decomposition of the tight-binding three-center integral between 
!  atom1 (angular momentum l1,m1) and atom2 (angular momentum l2,m2), when atom2 belongs to the neighborhood shell 
!  in two-center form in dependence on the two-center TB parameters and the direction cosines L,M,N.
! INPUT
!   * l1,l2  - angular momenta (l=-1 -> s* orbital)
!   * m1,m2  - magnetic quantum numbers
!   * shell  - shell number
! OUTPUT
!  Linear combination of two-center parameters in terms of direction cosines
! GTPack OPTIONS
!  GOTbBasis
!    o  0 -> no basis (standard)
!    o  basis "Cu,Zn" a list of chemical symbols 
!
!  GOTbRule
!    o -> 1 no rule applied
!    o -> 2 L^2 -> 1 - M^2 - N^2
!    o -> 3 M^2 -> 1 - L^2 - N^2
!    o -> 4 N^2 -> 1 - L^2 - M^2
!
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTTbAuxVogl, GTTbSymbol2C
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The application of transformation rules is sometimes necessary to get a simple form 
!  of the results. The option GOTbRule allows transformations.
!  
! LITERATURE
!  Podolskiy, Vogl, Phys. Rev. B 69, 233101 (2004)
! TODO
! -
! RELAESE
!  1.0.0
! PROBLEMS
!  -
! BUG FIXES
!   17.07.18 : Up to now it was only considered that GOTbBasis describes a case without a need of 
!   superscripts (monoatomic crystal). The other case sonsidered correctly was a basis of two DIFFERENT atoms.
!   If the atoms in the basis are equal ("Pb,Pb") in case of a real basis, superscripts are necessary but
!   (sps)^"Pb,Pb" and (pss)^"Pb,Pb" are the same. Hence, now only the name with lexicographic order of the
!   orbitals is generates
!
!	  
!--------------------------------------------------------------------------------
! SOURCE
*)
GTTbMatrixElement[l1_, m1_, l2_, m2_, dist_, OptionsPattern[]] := Module[
	{bas,l1s,l2s,ll,lu,rule,lls,sl,mp,c1,lt1,lt2,s1,c2,mel,cl,s,i,test,bas1}, 
	bas = OptionValue[GOTbBasis];
	bas1=bas;
	If[bas==0,
	   None,	
	   test=Length[StringSplit[bas,","]//Union];
	  If[test==1,
	  	 bas=0,
	  	 None
	  ];
	];  
	(*--- handle s* orbitals ---*)
	If[l1 == -1, 
		l1s = 0, 
		l1s = l1
	]; 
	If[l2 == -1, 
		l2s = 0, 
		l2s = l2
	];
    (*--- decide about the names:if a basis appears also (pss) is allowed ---*)
    ll = Min[l1, l2]; lu = Max[l1, l2]; 
    rule = {{}, {DL^2 -> 1 - DN^2 - DM^2}, {DM^2 -> 1 - DL^2 - DN^2}, {DN^2 -> 1 - DL^2 - DM^2}}; 
    If[ll == -1, 
    	lls = 0, 
    	lls = ll
    ];
    If[bas === 0, 
    	sl = Table[GTTbSymbol2C[ll, lu, mp, dist, GOTbBasis -> bas1], {mp, 0, lls}], 
        sl = Table[GTTbSymbol2C[l1, l2, mp, dist, GOTbBasis -> bas1], {mp, 0, lls}]
    ];
    c1 = 0;
    If[bas === 0, lt1 = ll; lt2 = lu, lt1 = l1; lt2 = l2];
    Do[s1 = (TbS[l1s, m1, Abs[mp]]*TbS[l2s, m2, Abs[mp]] + TbT[l1s, m1, Abs[mp]]*TbT[l2s, m2, Abs[mp]])*
       GTTbSymbol2C[lt1, lt2, Abs[mp], dist, GOTbBasis -> bas1];
       c1 = c1 + s1
    , {mp, 1, lls}];
    c2 = 2*TbA[m1]*TbA[m2]*Tbd[l1s, Abs[m1], 0]*Tbd[l2s, Abs[m2], 0]*GTTbSymbol2C[lt1, lt2, 0, dist, GOTbBasis -> bas1];
    mel = FullSimplify[(c1 + c2), {N > 0, M > 0, N > 0}]*(-1)^((l1s - l2s + Abs[l1s - l2s])/2);
    cl = 0;
    Do[
    	s = CoefficientList[mel, sl[[i]]]; 
        s = s[[2]]*sl[[i]] /. rule[[OptionValue[GOTbRule]]] // 
        FullSimplify; cl = cl + s
    , {i, 1, lls + 1}];
        cl = FullSimplify[cl, 
    Assumptions -> {DL^2 == 1 - DM^2 - DN^2, DM^2 == 1 - DN^2 - DL^2, 
      DN^2 == 1 - DL^2 - DM^2, 1 == (DL^2 + DM^2 + DN^2)}];
  cl = cl /. {DN -> Global`DN, DM -> Global`DM, DL -> Global`DL};
  Return[cl]]


(*
***) 

(****t* /GTTbSymbol2C
! NAME
!  GTTbSymbol2C
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 15.04.2014 : s* orbital is introduced
!  * 16.06.2014 : canonical form of the parameters
!  * 25.06.2018 : check header and description
! USAGE
!  GTTbSymbol2C[ll,lu,m,shell] gives a symbol according to the nomenclature of the two-center approximation.
! INPUT
!
!  * dist - shell numner (index for the two-center element)
!  * ll - lower angular momentum
!  * lu - upper angular momentum
!  * m  
! 
!    o larger than zero -> type of overlap
!      1,2,3,4     : sigma, pi, delta, fi
!    o less than zero -> construction of onsite names
!      -1,-2,-3,-4 : 0,1,2,0
!      s -> (ss0)  p -> (pp0)  d ->(dd0) d1 -> (dd1) d2 -> (dd2)  pd -> (pd0) s* -> (s*s*0)
!      according to the Slater&Koster tables
! OUTPUT
!  a symbol for a two-center element in the form of the original Slater-Koster paper
! GTPack OPTIONS
!  * GOTbBasis: 
!
!    - 0 -> no basis (standard)
!    - basis "Cu,Zn" a list of chemical symbols 
!
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  In real space calculations problems arise with respect to the names of the paramters.
!  In the construction parameters like (pds)^(Ga,As) AND (dps)^(As,Ga) might occur. The new version 
!  reformulates the parameters such, that the atoms are in lexicographic order. The angular
!  momenta are reordered accordingly.
! LITERATURE
!  - Podolskiy, Vogl, Phys. Rev. B 69, 233101 (2004)
!  - J.C. Slater, G.F. Koster, Simplified LCAO method for the periodic potential problem, Phys. Rev. 94, 1498 (1954)
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
! The s* orbital has angular momentum l=0 as a normal s-orbital too. Therfore l=-1
! is used for the s* orbital to distinguish the two cases.
!
!--------------------------------------------------------------------------------
! SOURCE
*)
GTTbSymbol2C[ll_, lu_, m_, dist_, OptionsPattern[]] := 
    Module[{am, type, bas, type1, tc, symb1, symb2, bas1, ll1, lu1, bas2,lx}, 
    am = {"s", "p", "d", "f", "g", "s*"};
    bas = OptionValue[GOTbBasis];
    (*--- canonical order for hopping parameter if a basis exists ---*)
    If[bas == 0,
       ll1 = ll; lu1 = lu,
       bas1 = StringSplit[bas, ","];
       If[Length[bas1] == 1,
          ll1 = ll; lu1 = lu,
          bas1 = Sort[bas1]; bas2 = bas1[[1]] <> "," <> bas1[[2]];
          If[bas2 == bas,
             ll1 = ll; lu1 = lu,
             lx = ll; ll1 = lu; lu1 = lx; bas = bas2
          ] 
       ]
    ];
    (*--- construction of the symbol ---*)
    If[ll1 == -1, 
       symb1 = am[[6]], 
       symb1 = am[[ll1 + 1]]
    ];
    If[lu1 == -1, 
       symb2 = am[[6]], 
       symb2 = am[[lu1 + 1]]
    ];
    type = {"\[Sigma]", "\[Pi]", "\[Delta]", "\[Phi]"};
    type1 = {"0)", "1)", "2)"};
    tc = "(" <> symb1 <> symb2;
    If[m >= 0, 
       tc = tc <> type[[m + 1]] <> ")";
       If[m >= 0 && bas === 0, 
       	  Return[Subscript[tc, dist]], 
          Return[Subsuperscript[tc, dist, bas]]
       ], 
       None
    ];
    If[bas === 0 && m < 0, 
       tc = tc <> type1[[Abs[m]]]; Return[tc], 
       tc = tc <> type1[[Abs[m]]];
       Return[Superscript[tc, bas]]
    ]
]

(*
***)

 




(****t* /GTTbSimplify
! NAME
!  GTTbSimplify
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  July 2013 : first version
! USAGE
!  Helps to simplify matrix elements after construction
! INPUT
!  * mat  - matrix element
!  * ll   - lower angular momentum
!  * lu   - upper angular momentum 
!  * dist - shell numner (index for the two-center element)
!  * rl   - a rule to express the direction cosines in a different way
! OUTPUT
!  simplified matrix element
! OPTIONS
!  GOTbBasis
!    o  0 -> no basis (standard)
!    o  basis "Cu,Zn" a list of chemical symbols 
!
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTTbSymbol2C
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  This is an internal function of the package TightBinding.m
! LITERATURE
!  Podolskiy, Vogl, Phys. Rev. B 69, 233101 (2004)
! TODO
!  -
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTbSimplify[mat_,ll_,lu_,dist_,rl_,OptionsPattern[]]:=Module[{sl,mp,cl,sx,i,s,bas},
	bas=OptionValue[GOTbBasis];
	sl=Table[GTTbSymbol2C[ll,lu,mp,dist,GOTbBasis->bas],{mp,0,ll}];
	cl=0;
    Do[
       s=CoefficientList[mat,sl[[i]]]//FullSimplify;sx=ReplaceAll[s[[2]],rl]*sl[[i]];
	   sx=FullSimplify[sx];cl=cl+sx
	,{i,1,ll+1}];
	Return[cl]
]

(*
***) 

(****t* /TbDcos
! NAME
!  TbDcos
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  JulY 2013 : first version
! USAGE
!  TBDcos[lat] calculates the direction cosines of a list of lattice vectors
! INPUT
!  lat - list of lattice vectors
! OUTPUT
!  direction cosines of list of lattice vectors
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  This is an internal function of the package TightBinding.m
! LITERATURE
!  Podolskiy, Vogl, Phys. Rev. B 69, 233101 (2004)
! TODO
!  -
! PROBLEMS
!  -
!
!--------------------------------------------------------------------------------
! SOURCE
*)

TbDcos[lat_]:=Module[{b,lat1,i},
	 b=Map[Sqrt[#.#]&,lat];lat1=lat/b;
	 Do[
	 	lat1[[i]]={{1,0,0}.lat1[[i]],{0,1,0}.lat1[[i]],{0,0,1}.lat1[[i]]}
	 ,{i,1,Length[lat]}];
	 Return[lat1]
]

(*
***) 


(****t* /GTTbHamiltonianElement
! NAME
!  GTTBHamiltonianElement
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 14.02.2014 : units (2 Pi/a) introduced
!  * 15.04.2014 : s* orbital 
!  * 25.06.2018 : check header and descriotion
! USAGE
!  GTTbHamiltonianElement[l1,m1,l2,m2,shell,shell vectors] constructs the k-dependent 
!  contribution of shell characterized by the shell vectors to the tight-binding matrix element 
!  between functions of symmetry l1,m1 and l2,m2.
!
! INPUT
!  * l1, m1    - angular momentum information 1st orbital  (if l=-1 s* orbital)
!  * l2, m2    - angular momentum information 2nd orbital  (if l=-1 s* orbital)
!  * dist      - distance, repective the shell number
!  * lat       - list of lattice vectors  in units of a
! OUTPUT
!  the corresponding analytical form of the tight-binding matrix element
! GTPack OPTIONS 
!  * GOTbBasis:
!
!    o  0 -> no basis (standard)
!    o  basis "Cu,Zn" a list of chemical symbols 
!
! STANDARD OPTIONS
!  -
! GTPack MODULES
!   GTTbMatrixElement, GTTbSymbol2C, TbDcos
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  Calculation of an element of the Hamiltonian matrix. The k-vector is
!  assumed to be given in units of 2 pi/a,  k={xi,eta,zeta}. 
!  The list of lattice vectors correspond to a certain shell. 
!  The shell information and the information about the basis is used to construct
!  the parameter symbols accordingly.
!
! LITERATURE
!  Podolskiy, Vogl, Phys. Rev. B 69, 233101 (2004)
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)
GTTbHamiltonianElement[l1_, m1_, l2_, m2_, dist_, lat_, OptionsPattern[]] := Module[{bas, hm, dc, mat, i, mm, mat1, l1s, l2s, rsp}, 
    bas = OptionValue[GOTbBasis]; rsp = OptionValue[GOTbRealSpace]; 
    dc = TbDcos[lat]; hm = 0;
    mat = GTTbMatrixElement[l1, m1, l2, m2, dist, GOTbBasis -> bas];
    Do[
      If[dc[[i, 3]] == 1, 
         If[bas == 0, 
            mat1 = GTTbSymbol2C[Min[l1, l2], Max[l1, l2], Abs[m1], dist, GOTbBasis -> bas], 
            mat1 = GTTbSymbol2C[l1, l2, Abs[m1], dist, GOTbBasis -> bas]
         ];
         If[l1 == -1, l1s = 0, l1s = l1];
         If[l2 == -1, l2s = 0, l2s = l2];
         mm = KroneckerDelta[m1, m2]*mat1*(-1)^((l1s - l2s + Abs[l1s - l2s])/2), 
         mm = mat /. {Global`DL -> dc[[i, 1]], Global`DM -> dc[[i, 2]], Global`DN -> dc[[i, 3]]}
      ];
      If[rsp, 
      	 hm = hm + mm, 
         hm = hm + Exp[I 2 \[Pi] xez.lat[[i]]]*mm
      ]
    , {i, 1, Length[lat]}];
    If[rsp, 
       None, 
       hm = ExpToTrig[hm] // FullSimplify
    ];
    Return[hm]
]

(*
***) 



(****t* /GTTbParameterSet
! NAME
!  GTTbParameterSet
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
! * July 2013  : first version
! * April 2014 : s* orbital April 2014
! USAGE
!  This module helps to form the parameter data sets, if the parameters are put in by hand.
! INPUT
!  * llst      - a list of angular momenta
!  * dmax      - maximum distance, i.e. shell number
! OPTIONS
!  GOTbBasis
!    o  0 -> no basis (standard)
!    o  basis "Cu,Zn" a list of chemical symbols 
!
! OUTPUT
!  The parameter set  contains the parameter symbols, constructed from the shell information
!  and the basis information, and the parameter values
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTTbSymbol2C
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  
! LITERATURE
!  Podolskiy, Vogl, Phys. Rev. B 69, 233101 (2004)
! TODO
!  
! PROBLEMS
!  
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTbParameterSet[llst_, dmax_, OptionsPattern[]] := Module[{p2c,bas,lmax,ind,ln,ln1,index1,tab,tab1,i,j,k,lt,lu,
	lu1,lu2,ll1,ll,ll2,bas1,symb,m}, 
	p2c = {};bas = OptionValue[GOTbBasis]; lmax = Length[llst];
   (*--- check if basis,list of basis atoms ---*)
   If[
    (*--- only one atom in basis ---*)
       bas == 0, ind = 0; ln = 1; ln1 = 1; index1 = {0},
    (*--- more atoms in basis ---*)
       tab = StringSplit[bas, ","]; ln = Length[tab]; tab1 = {}; 
       ind = 1; index1 = {};
    (*--- create a list of combinations of basis atoms ---*)
       Do[
          Do[
             tab1 = Append[tab1, StringJoin[{tab[[i]], ",", tab[[j]]}]];
             If[i == j, 
             	index1 = Append[index1, 0], 
             	index1 = Append[index1, 1]
             ]
          , {i, 1, j}]
       , {j, 1, ln}];
       ln1 = Length[tab1]; tab = tab // Flatten 
   ];
    (*--- if basis go over all combinations, otherwise only over all paramters ---*)
   Do[
      Do[
         Do[
         	(*--- decide if (sps) and (pss) are necessary ---*)
            If[index1[[k]] == 0,
               lt = lu,
               lt = lmax
            ]; 
            lu1 = llst[[lu]]; 
            If[llst[[lu]] == -1, 
               lu2 = 0, 
               lu2 = llst[[lu]]
            ];
            Do[
               ll1 = llst[[ll]]; 
               If[llst[[ll]] == -1, 
               	  ll2 = 0, 
               	  ll2 = llst[[ll]]
               ];
               If[ind == 1, 
               	  bas1 = tab1[[k]], 
               	  bas1 = 0
               ];
               Do[
                  symb = GTTbSymbol2C[ll1, lu1, m, i, GOTbBasis -> bas1];
                  p2c = Append[p2c, {symb, Input[symb]}]
               , {m, 0, Min[ll2, lu2]}]
            , {ll, 1, lt}]
         , {lu, 1, lmax}]
      , {i, 1, dmax}]
   , {k, 1, ln1}];
  (*--- append on-site paramters ---*)
   Do[
   	  If[ind == 1, 
   	  	 bas1 = tab[[k]], 
   	  	 bas1 = 0
   	  ];
      symb = GTTbSymbol2C[-1, -1, -1, 1, GOTbBasis -> bas1]; 
      p2c = Append[p2c, {symb, Input[symb]}];
      symb = GTTbSymbol2C[0, 0, -1, 1, GOTbBasis -> bas1]; 
      p2c = Append[p2c, {symb, Input[symb]}];
      symb = GTTbSymbol2C[1, 1, -1, 1, GOTbBasis -> bas1]; 
      p2c = Append[p2c, {symb, Input[symb]}];
      symb = GTTbSymbol2C[2, 2, -1, 1, GOTbBasis -> bas1]; 
      p2c = Append[p2c, {symb, Input[symb]}];
      symb = GTTbSymbol2C[2, 2, -2, 1, GOTbBasis -> bas1]; 
      p2c = Append[p2c, {symb, Input[symb]}];
      symb = GTTbSymbol2C[2, 2, -3, 1, GOTbBasis -> bas1]; 
      p2c = Append[p2c, {symb, Input[symb]}];
      symb = GTTbSymbol2C[1, 2, -1, 1, GOTbBasis -> bas1]; 
      p2c = Append[p2c, {symb, Input[symb]}]
   , {k, 1, ln}];
   Return[p2c]
]

(*
***) 


(****t* /GTTbGetParameter
! NAME
!  GTTbGetParameter
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 15.07.2013 : first version
!  * 13.06.2016 : error handling included 
!  * 25.06.2018 : check header and documentation
! USAGE
!  GTTbGetParameter[set, symb] gives the value of a parameter from a parameter set.
! INPUT
!  * set  -  the parameter set 
!  * symb -  symbol of the parameter
! OUTPUT
!  the corresponding parameter value of the set
! GTPack OPTIONS 
!  -
! STANDARD OPTIONS
!  - 
! GTPack MODULES
! - 
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
!--------------------------------------------------------------------------------
! SOURCE
*)
GTTbGetParameter[set_,sym_]:=Module[{pos,parm,pp},
	pos = Position[Transpose[set][[1]], sym] // Flatten;
    If[pos == {} || Length[pos] > 1, 
       Print["Symbol not in parameter set, or choice not unambiguous"]; 
       Return[],
       pp=pos[[1]];parm=set[[pp,2]];
       Return[parm]
    ]
 ]

(*
***) 


(****t* /GTTbSetParameter
! NAME
!  GTTbSetParameter
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  June 2016 : first version
! USAGE
!  GTTbSetParameter[set, sym, val] changes the value of parameter sym to val
! INPUT
!  * set  -  the parameter set 
!  * symb -  symbol of the parameter
!  * val  -  new value of the parameter
! OUTPUT
!  the corresponding parameter value of the set
! OPTIONS
!  - 
! ERROR MESSAGES
!  -
! GTPack MODULES
! - 
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  introduce GOVerbose for a pretty print of the modified parameter set.
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTbSetParameter[set_, sym_, val_] := Module[{set1,pos}, 
	set1 = set;
    pos  = Position[Transpose[set1][[1]], sym] // Flatten;
    If[pos == {} || Length[pos] > 1, 
       Print["Symbol not in parameter set, or choice not unambiguous"]; 
       Return[],
       set1[[pos, 2]] = val;Return[set1]
    ]
 ]
  
 (*
***)  

(****t* /GTTbDatabaseInfo
! NAME
!  GTTbDatabaseInfo
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 08.08.2013 : 1st version 
!  * 13.03.2015 : check of database name and parmset
!  * 06.07.2016 : nice print of th table consecutive numbering included
!  * 15.08.2016 : check of database name and parmset fixed for subdirectories
!  * 25.06.2018 : check header and documentation
!  * 21.02.2022 : check for errors included and Mathematica Messagesystem used.
! USAGE
!  GTTbDatabaseInfo[database] gives information about the tight-binding parameter sets available in database.
! INPUT
!  database  - name of the database
!  
! OUTPUT
!  database and printout of the information
! GTPack OPTIONS
!  -
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
! - 
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
! SOURCE
*)
GTTbDatabaseInfo::data =  "The database `1` was not found in the working directory." ;
GTTbDatabaseInfo::empty = "The database is empty." ;

GTTbDatabaseInfo[db_] := Module[{dbp, dbname,head,info,nr,i,tt}, 
	(*--- Append file extension, if necessary and check if database and parameters exist ---*)
	dbname=StringTrim[db, ".parm"] <> ".parm";
    If[!FileExistsQ[dbname], 
       Message[GTTbDatabaseInfo::data, dbname]; Abort[]
    ];
   head = {"Number", "Name", "Structure", "Authors", "Reference"};
  (*---print info---*)
  dbp = ReadList[dbname][[1]];
  If[dbp == {},
     Message[GTTbDatabaseInfo::empty]; Abort[],
     None
  ];
  info = Transpose[Take[Transpose[dbp], 4]]; 
  nr = Table[i, {i, 1, Length[info]}];
  tt = Join[{nr}, info // Transpose] // Transpose; 
  tt = Join[{head}, tt];
  Print[Grid[tt, Frame      -> All, 
                 Dividers   -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
                 Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} -> GTCornerColor}
            ]
  ];
  Return[dbp]
]

(*
GTTbDatabaseInfo[db_]:=Module[{dbp,dbname,info},
	   dbname=StringJoin[db,".parm"];
	   (*--- Check if database exist s ---*)
	   If[Flatten[Position[FileNames[], dbname]] == {},
         Print["Error : database not in working directory"];Abort[],
         None
       ];
	   (*--- print info ---*)
	   dbp=ReadList[dbname][[1]];
	   info=Transpose[Take[Transpose[dbp],4]];
       Print[TableForm[info,TableHeadings->{None,{"Name","Structure","Authors","Reference"}}]];
       Return[dbp]
]
*)

(*
***) 



(****t* /GTTbDatabaseUpdate
! NAME
!  GTTbDatabaseUpdate
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 08.08.2013 : first version
!  * 21.042014  : s* orbital introduced 
!  * 05.03.2015 : check if database exists
!  * 24.07.2016 : input of SO paramters implemented, better check of data
!  * 13.08.2016 : check of database name and parmset fixed for subdirectories
!  * 25.06.2018 : check header nad documentation
! USAGE
!  GTTbDatabaseUpdate[database] reads the database, adds a new parameter set and stores the database.
! INPUT
!  db  - name of the database
!
! OUTPUT
!  updated database stored and is available as output of the command
! GTPAck OPTIONS
!  -
! STANDARD OPTIONS
!  -
!
! GTPack MODULES
!  GTTbParameterSet, GTReadFromFile
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The input of data is organized interactively, using the routine INPUT 
!  of Mathematica. The input requires the following steps:
!     Input           Type                Example 
!     -------------------------------------------  
!     Name           String               "BTO"   
!     Structure      String               "Perovskite"
!     Source         String               "Book T. Wolfram"
!     Input          Number               1
!     Number atoms   Number               3
!     Names atoms    String               "Ba,Ti,O"
!     Angular mom.   List of Strings      {"s","p","d"} 
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  The angular momenta are defined for all atoms in the basis in the same way, i.e
!  the parameter set contains more entries than necessary.
!  The solution of this problem would require a totally new routine.
! 
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTbDatabaseUpdate[db_] := Module[
	{dbname, mom, dbp, name, str, src, aut, it, basis, llst, dmax,pset, file,atoms,so,namex,parm,l,sow,i,j},
  (*---prepare input---*)
   mom    = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4};
  (*--- Append file extension, if necessary and check if database and parameters exist ---*)
   dbname=StringTrim[db, ".parm"] <> ".parm";
   If[!FileExistsQ[dbname], 
     Print["Error: the database " <> dbname <> " was not found in the working directory!"]; Abort[]
   ];
  (*--- Read database ---*)
   dbp  = ReadList[dbname][[1]];
  (*--- input of general information with tests ---*)  
   name = Input["Name of Element or Compound (String)"];
   If[Head[name] === String,
      None,
      Print["Error: Name of Element is not a string."]; Abort[]
   ];
   str = Input["Structure (String)"];
   If[Head[str] === String,
      None,
      Print["Error: Structure is not a string."]; Abort[]
   ];
   aut = Input["Author (String)"];
   If[Head[aut] === String,
      None,
      Print["Error: Author is not a string."]; Abort[]
   ];
   src = Input["Reference (String)"];
   If[Head[src] === String,
      None,
      Print["Error: Reference is not a string."]; Abort[]
   ];
   it  = Input["Input parameters (1), read from file (2)"];
   If[Head[it] === Integer,
      None,
      Print["Error: Input selection is not an Integer"]; Abort[]
   ];
   If[it == 1,
      basis = Input["number of basis atoms (=1 if the atoms in a basis are equal, see Si)"];
      atoms = basis;
      If[Head[basis] === Integer,
         None,
         Print["Error: Number of basis atoms is not an Integer"]; Abort[]
      ];
      If[basis > 1,
         basis = Input["List of basis atoms (String, atom names separated by comma)"];
         If[Head[basis] === String,
            None,
            If[Head[basis] === List,
               Print["Error: List of atoms has to be a comma separated string NOT A LIST"]; Abort[],                     
               Print["Error: List of basis atoms is not a string"]; Abort[]
            ]
         ],
         basis = 0
      ];
      llst = Input["List of angular momenta , Strings {s,p,d}"];
      llst = llst /. mom;
      dmax = Input["Number of Shells"];
   (*--- Input of parameters ---*)
     pset = GTTbParameterSet[llst, dmax, GOTbBasis -> basis];
   (*--- SO parameter ---*)  
      so  = Input["Spin-Orbit coupling (True/False)"];
      If[so,
         Do[
            parm  = Input["angular momenta and SO-parameters {Atom,{l1,so1},...}"];
            namex = parm[[1]];
            Do[
               l   = parm[[i, 1]]; 
               sow = parm[[i, 2]];
               If[atoms == 1,
                  pset = Append[pset, {"\[Xi]", sow}],
                  pset = Append[pset,  {Subsuperscript["\[Xi]", 1, namex], sow}]
               ]
            , {i, 2, Length[parm]}]
       , {j, 1, atoms}],
       None
     ],
     None
   ];
   If[it == 2,
      file = Input["file name"];
      pset = GTReadFromFile[file],
      None
   ];
   dbp = Append[dbp, {name, str, aut, src, pset}];
   GTWriteToFile[dbp, dbname];
   Return[dbp]
]
 
(*
GTTbDatabaseUpdate[db_] := Module[{dbname,mom,dbp,name,str,src,aut,it,basis,llst,dmax,pset,file}, 
	 dbname = StringJoin[db, ".parm"];
	 (*--- Check if database exists ---*)
	 If[Flatten[Position[FileNames[], dbname]] == {},
        Print[dbname, " -> database not in working directory"];Abort[],
        None
	 ];
	 (*--- input of data ---*)
     mom = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4}; 
     (*--- Read database ---*)
     dbp = ReadList[dbname][[1]];
     (*--- input of general information ---*)
     name = Input["Name of Element or Compound (String)"];
     str  = Input["Structure (String)"];
     aut  = Input["Author (String)"];
     src  = Input["Reference (String)"];
     it   = Input["Input parameters (1), read from file (2)"];
     If[it == 1, 
        basis = Input[
        	"number of basis atoms (=1 if the atoms in a basis are equal, see Si)"];
        If[basis > 1, 
        basis = Input[
           "List of basis atoms (String, atom names separated by comma)"], 
        basis = 0];
        llst = Input["List of angular momenta , Strings {s,p,d}"];
        llst = llst /. mom;
        dmax = Input["Number of Shells"];
        (* Input of parameters *)
        pset = GTTbParameterSet[llst, dmax, GOTbBasis -> basis], 
        None
     ];
     If[it == 2, 
     	file = Input["file name"]; 
     	pset = GTReadFromFile[file], 
     None];
     dbp = Append[dbp, {name, str, aut, src, pset}];
     Put[dbp, dbname];
     Return[dbp]
]
  
*)
(*
***) 


(****t* /GTTbDatabaseRetrieve
! NAME
!  GTTbDatabaseRetrieve
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 08.08.2013 : 1st version 
!  * 14.03.2015 : check if database and parameter set exist
!  * 20.07.2016 : index to find the dataset can be the name of element or compund or the number /position in the dataset
!  * 08.08.2016 : check of database name and parmset fixed for subdirectories
!  * 25.06.2018 : check of header and documentation
! USAGE
!  GTTbDatabaseRetrieve[db, name] loads a tight-binding parameter set from a given database.
! INPUT
!  * db   - name of the database
!  * name - name of the parameter set
! OUTPUT
!  parameter set
! GTPack OPTIONS
!  -
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
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTbDatabaseRetrieve[db_, name_] := Module[{dbp, dbname, nset, pos, elm, pset}, 
	(*--- Append file extension, if necessary and check if database and parameters exist ---*)
	dbname=StringTrim[db, ".parm"] <> ".parm";
    If[!FileExistsQ[dbname], 
       Print["Error: the database " <> dbname <> " was not found in the working directory!"]; Abort[]
    ];
	(*---retrieve the parameter set---*)
    dbp  = ReadList[dbname][[1]];
    elm  = Take[Transpose[dbp], 1] // Flatten;
    nset = Length[dbp];
    If[Head[name] === Integer,
       If[name <= nset,
          pset = dbp[[name, 5]],
          Print["Error: The requested parameter set is not found in the file " <> dbname <> "."]   ; Abort[] 
       ],
       pos = Position[elm, name];
       If[Length[pos] == 0,
          Print["Error: The requested parameter set is not found in the file " <> dbname <> "."]; Abort[], 
          pset = dbp[[First@First@pos, 5]]
       ]
    ];
    Return[pset]
 ]

  
(*
GTTbDatabaseRetrieve[db_, name_] := Module[{dbp, dbname, nset, pos},
  dbname = StringJoin[db, ".parm"];
  (*--- Check if database and parameters exist ---*)
  If[Flatten[Position[FileNames[], dbname]] == {},
     Print["Error: the database "<>dbname<>" was not found in the working directory!"];Abort[],
     None
  ];
  (*--- retrieve the parameter set ---*)
(*  If[Position[dbp, name] == {},
     Print[name, " -> parameters not in database"]; Abort[], 
     None
  ];*)
  dbp = ReadList[dbname][[1]];
  nset = Length[dbp];
  pos = Position[dbp, name];
  If[Length[pos] == 0, 
  	Print["Error: The requested parameter set is not found in the file " <>dbname<>"."]; Abort[],
  	
  	Return[dbp[[First@First@pos, 5]]]
   ]
]

*)
(*
***) 


(****t* /GTTbPrintParmSet
! NAME
!  GTTbPrintParmSet
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 08.08.2013 : 1st version
!  * 15.03.2015 : check if database and parmset exist
!  * 09.07.2016 : it is now possible to take the name of element or compund or the number/position in the
!                 database to select the dataset. Print of table in standard form now.
!  * 10.08.2016 : check of database name and parmset fixed for subdirectories
! USAGE
!  GTTbPrintParmSet[db, name] prints a tight-binding parameter set from a database.
! INPUT
!  * db   - name of the database
!  * name - name of the parameter set
! OUTPUT
!  pretty print of the parameter set
! GTPack OPTIONS
!  -
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
!  -
! PROBLEMS
!  -
!-------------------------------------------------------------------------------
! SOURCE
*)


GTTbPrintParmSet[db_, name_] := Module[{ps, parm, names, no, n, i, dbname, nset, pset,dbp,elm,pos,plt}, 
	(*--- Append file extension, if necessary and check if database and parameters exist ---*)
	dbname=StringTrim[db, ".parm"] <> ".parm";
    If[!FileExistsQ[dbname], 
       Print["Error: the database " <> dbname <> " was not found in the working directory!"]; Abort[]
    ];
  (*---read the parmset---*)
    dbp  = ReadList[dbname][[1]]; 
    nset = Length[dbp];
    elm  = Take[Transpose[dbp], 1] // Flatten;
 (*--- prepare the print ---*)
    If[Head[name] === Integer,
       If[name <= nset,
          pset = dbp[[name]],
          Print["Error: The requested parameter set is not found in the file " <> dbname <> "."]; Abort[] 
       ],
       pos = Position[elm, name];
       If[Length[pos] == 0,
          Print["Error: The requested parameter set is not found in the file " <> dbname <> "."]; Abort[], 
          pset = dbp[[First@First@pos]]
       ]
    ];
    Print[{{"Name       :", pset[[1]]}, 
    	   {"Structure  :", pset[[2]]}, 
    	   {"Authors    :", pset[[3]]}, 
    	   {"Reference  :", pset[[4]]}} // TableForm];
     ps    = pset[[5]]; 
     parm  = Transpose[ps][[2]]; 
     names = Transpose[ps][[1]];
     n = Length[names]; no = IntegerPart[n/10]; n1 = n - no*10;
     Do[ 
     	plt = {Take[names, {(i - 1)*10 + 1, i*10}], Take[parm, {(i - 1)*10 + 1, i*10}]}; 
        Print[Grid[plt, Frame -> All, 
                        Dividers -> {Black, {2 -> GTDividerColor1}}, 
                        Background -> {None, {1 -> GTBackGroundColor1}}
                  ]
        ];
     , {i, 1, no}];
     plt = {Take[names, {no*10 + 1, n}], Take[parm, {no*10 + 1, n}]};
     Print[Grid[plt, Frame -> All, 
     	             Dividers -> {Black, {2 -> GTDividerColor1}},
                     Background -> {None, {1 -> GTBackGroundColor1}}
               ]
     ]
  ]
    
    
(*
GTTbPrintParmSet[db_,name_]:=Module[{ps,parm,names,no,n,n1,i,dbp,dbname,nset,pset},
	          dbname=StringJoin[db,".parm"];
	          (*--- Check if database and parameters exist ---*)
	          If[Flatten[Position[FileNames[], dbname]] == {},
                 Print[dbname, " -> database not in working directory"];Abort[],
                 None
              ];
	          (*--- print the parmset ---*)
	          dbp=ReadList[dbname][[1]];nset=Length[dbp];
              If[Position[dbp, name] == {},
                 Print[name, " -> parameters not in database"]; Abort[], 
                 None
              ];
              Do[If[dbp[[i,1]]==name,pset=dbp[[i]],None],{i,1,nset}];
              Print[{{"Name       :",pset[[1]]},{"Structure  :",pset[[2]]},{"Authors    :",pset[[3]]},{"Reference  :",pset[[4]]}}//TableForm];
              ps=pset[[5]];parm=Transpose[ps][[2]];names=Transpose[ps][[1]];
              n=Length[names];no=IntegerPart[n/10];n1=n-no*10;
              Do[
                 Print[TableForm[{Take[parm,{(i-1)*10+1,i*10}]},TableHeadings->{None,Take[names,{(i-1)*10+1,i*10}]}]],
              {i,1,no}];
              Print[TableForm[{Take[parm,{no*10+1,n}]},TableHeadings->{None,Take[names,{no*10+1,n}]}]]
]
*)
(*
***) 

(****t* /GTTbParmToRule
! NAME
!  GTTbParmToRule
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 08.08.2013 : first version
!  * 25.06.2018 : check header and documentation
! USAGE
!  GTTbParmToRule[set] gives a rule to replace tight-binding symbols using a given parameter set.
! INPUT
!  set - parameter set
! OUTPUT
!  list of substitution rules
! GTPack OPTIONS
!  -
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
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTbParmToRule[set_]:=Module[{ln,rl},ln=Length[set];rl={};
                         Do[rl=Append[rl,{set[[i,1]]->set[[i,2]]}],{i,1,ln}];
                         Return[Flatten[rl]]
]

(*
***) 



(****t* /GTTbParmToFortran
!  NAME
!  GTTbParmToFortran
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  September 2013 : 1st version 
! USAGE
!  Transformation of the parameter set into a list of substitution rules for FORTRAN output.
!  Additionally {xi,eta,zeta} will be transformed to {kx,ky,kz}
! INPUT
!  set - parameter set
! OUTPUT
!  list of substitution rules
! OPTIONS
!  -
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  For the export of the constructed Hamiltoninan in FORTRAN form the two-center parameters like (sss)_1 have to
! be substituted by a form which is accepted as a FORTRAN name. The paramters are given in a form of simple list pt.
!
!  This is an internal module of the package TightBinding.m
! LITERATURE
!  -
! TODO
!  -
! PROBLEMS
! - 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTTbParmToFortran[set_]:=Module[{ln,rl,pt,set1,i},
	       ln=Length[set];set1=set;rl={};
	       pt=Table[p[i],{i,1,Length[set]}];
           Do[
           	  rl=Append[rl,{set1[[i,1]]->pt[[i]]}]
           ,{i,1,ln}];
           rl={rl,{\[Xi]->kx,\[Eta]->ky,\[Zeta]->kz}};
           Return[Flatten[rl]]
]
                         
(*
***) 
                         

(****t* /GTTbToFortranList
! NAME
!  GTTbToFortranList
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
! * 13.09.2013 : first version
! * 15.07.2015 : 2nd versiom
!
!     o  renamed to GTTbToFortranList (former GTTbHamiltonianToFortran), 
!     o  Option GOVerbose included
!
! * 25.06.2018 : check header and documentation
! USAGE
!  GTTbToFortranList[ham,elem,database] prints a Hamiltonian as FORTRAN code.
! INPUT
!  * ham  - Hamiltonian
!  * elem - element of database
!  * database  - database of TB Parameters
! OUTPUT
!  print of Matrixelements in FORTRAN Form
! GTPack OPTIONS
!   GOVerbose 
!
!     o True  - ouput of the substitution rule of the parameters
!     o False - no output
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTTbDatabaseRetrieve, GTTbParmToFortran
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
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTbToFortranList[ham_, elem_, base_,OptionsPattern[]] := Module[{opt,pset,rule,ndim,hamf,i,j},
  opt=OptionValue[GOVerbose];
  pset = GTTbDatabaseRetrieve[base, elem];
  rule = GTTbParmToFortran[pset];
  ndim = Length[ham]; hamf = ham /. rule;
  Do[
   Do[Print[
     StringForm["h(``,``)=``", i, j, FortranForm[hamf[[i, j]]]]
     ]
    , {j, 1, ndim}]
   , {i, 1, ndim}];
   If[opt==True,
   	 Do[Print[rule[[i]]], {i, 1, Length[rule] - 3}],
   	 None
   ]	 
  ]

(*
***)



(****t* /GTTbHamiltonian
! NAME
!  GTTbHamiltonian
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * February 2014: 
!
!     Unification of cases -
!     Up to now the case wih one atom in the basis was different in comparison to the case with more atoms.
!     The list lat has to start with 0 if only one atom is present, i.e. the list from GTShells
!     could not be used directly. This is fixed. Also in the one atom per basis case the shell list 
!     can be used directly. 
!  * June 2014:
!
!     The situation becomes more complicate if substances like BTO are considered (equal atoms in
!     the basis). It was hard to do all the cases: one atom, several different atoms and true basis in 
!     one module. In the moment it is split in three different submoduls.
!  * September 2014:
!
!     TBHamiltonianTrueBasis generated Superscripts even if only one sort of atoms are in the basis.
!     This is changed. In such a case no superscripts will be used now.
!  * 23./24.12.2016:
! 
!	  TBHamiltonianBasis changed. If the atoms in the basis are equal, no indices are added at the paprameters.
!     There was a problem with the onsite elements after this change in the normal case. The onsite elements
!     got a double index. It is corrected.
!  * 25.06.2018 check header and description
! USAGE
!  GTTbHamiltonian[basis,lattice] constructs the k-dependent Hamiltonian from the information about basis and lattice.
! INPUT
!  * basis         - contains the information about the basis in the form
!                  {("atom1", natom, llist1},{"atom2",llist2},...}
!
!                  natom - numer of atoms of this type                   
!                  llist - angular momenta included in the construction of the Hamiltonian     
!  * lattice       - lattice information (The lattice is the corresponding output of GTShells
! OUTPUT
!  the Hamiltonian in terms of the two-center tight binding parameters in analytical form 
! GTPack OPTIONS
!  * GOVerbose:
!
!      - True  - additional information
!      - False - no addtional information (Standard)
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTTbHamiltonianElement, GTTbSymbol2C, 
!  TBHamiltonianOneAtom, TBHamiltonianBasis, TBHamiltonianTrueBasis
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  GTTbHamiltonian is a master module, that calls one of the following modules in accordance with the description of the basis:
!   * TbHamiltoninaOneAtom - one atom in the basis
!   * TbHamiltoninanBasis - several equivalent atoms in the basis (for example Si)
!   * TbHamiltonianTrueBasis - a real basis and multiple atoms of the same sort are also allowed.
! LITERATURE
!  Podolskiy, Vogl, Phys. Rev. B 69, 233101 (2004)
! TODO
!  Like it is now, it seems that it will be possible to unify all cases again. Anyway the advantage of the structured form is, that
! one has a better chance to check all in case of mistakes.
! RELEASE
!  1.0.0
! PROBLEMS
!  The whole procedure is not general enough. It is restricted tot he cubic case. Hexagonal structure dos not work.
!--------------------------------------------------------------------------------
! SOURCE
*)
(* hcp Test 24.12.
GTTbHamiltonian[basis_, lat_, OptionsPattern[]] :=
  Module[{nbas, tbas, nat, ham, verb},
  verb = OptionValue[GOVerbose];
  (*---select the correct case---*)
  nbas = Length[basis];
  tbas = Transpose[basis];
  nat = Plus @@ tbas[[2]];
  (*---one atom in basis---*)
  
  If[nat == 1, ham = TbHamiltonianOneAtom[basis, lat]; Return[ham], 
   None];
  If[nat == nbas, 
(*   ham = TbHamiltonianBasis[basis, lat, GOVerbose -> verb]; *)
  ham = TbHamiltonianBasis[basis, lat];
   Return[ham], None];
  If[nat > nbas, ham = TbHamiltonianTrueBasis[basis, lat]; 
   Return[ham], None]]
*)   
   

GTTbHamiltonian[basis_, lat_] := Module[{nbas,tbas,nat,ham},
  (*--- select the correct case ---*)
  nbas = Length[basis];
  tbas = Transpose[basis];
  nat = Plus @@ tbas[[2]];
  (*--- one Atom in the Basis ---*)
  If[nat == 1,
      ham = TbHamiltonianOneAtom[basis, lat]; Return[ham],
      None
   ];
  If [nat == nbas,
  	
      ham = TbHamiltonianBasis[basis, lat]; Return[ham],
      None
   ];
   If [nat> nbas,
   
      ham = TbHamiltonianTrueBasis[basis, lat]; Return[ham],
      None
   ]
  ]


TbHamiltonianOneAtom[bas_, lati_] := 
       Module[{latt, tbas, am, am1, nam, ldl, ndim, ham, l1, l2, il, jl, k, dmax, latt1, ro, co, su, i,mom,im,jm},
       mom = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4};  
       latt = lati; latt[[1, 1]] = 0;
       tbas = Transpose[bas];
       am = Flatten[tbas[[3]]] /. mom;
       am1 = am /. {-1 -> 0};
       nam = Length[am1];
       ldl = Map[(2 # + 1) &, am1]; ndim = Plus @@ ldl;
       ham = Table[0, {ndim}, {ndim}];
       Do[l1 = am1[[il]];
            Do[l2 = am1[[jl]];
                 Do[
                       Do[
      (*--- calulate the position of the matrix element ---*)                           
                          ro = Sum[ldl[[k]], {k, 1, il - 1}] + l1 + im + 1;
                          co = Sum[ldl[[k]], {k, 1, jl - 1}] + l2 + jm + 1;
      (*--- lattice information ---*)
                          dmax = latt[[1, 2]]; latt1 = latt[[1, 3]];
                          If[ro < co,
      (*--- add the contributions from the shells (nondiagonal elements) ---*)
                             su = 0;
                             Do[
                                su = su + GTTbHamiltonianElement[l1, im, l2, jm, i, latt1[[i]], GOTbBasis -> 0]
                             , {i, 1, dmax}];
                             ham[[ro, co]] = su;
                             ham[[co, ro]] = ComplexExpand[Conjugate[su]],
                             None
                          ];
      (*--- add contributions from shells+on - site (diagonal elements) ---*)
                          If[ro == co,
                             su = 0;
                             Do[
                             	su = su + GTTbHamiltonianElement[l1, im, l2, jm, i, latt1[[i]], GOTbBasis -> 0]
                             , {i, 1, dmax}];
                             su = su + GTTbSymbol2C[l1, l2, -1, 1, GOTbBasis -> 0];
                             ham[[ro, co]] = su,
                             None
                          ]
                     , {jm, -l2, l2}]
               , {im, -l1, l1}]
          , {jl, 1, nam}]
      , {il, 1, nam}];
     Return[ham]
]

(*
 TbHamiltonianBasis[basis_, lat_, OptionsPattern[{GOVerbose -> False}]] := Module[
   {is, nbas, am, nam, ndim, ldl, i, dbl, ham, atom1, atom2, js, 
    l1, l2, il, jl, ro, co, k, im, jm, bas1, dmax, latt1, su, bas2, 
    mom, am1, l1s, l2s, latt, tbas, names, sel, p1, verb,test,test1,mylog,bas3,bas4,
    onsite},
  verb = OptionValue[GOVerbose];
  nbas = Length[basis];
  mom = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4};
  tbas = Transpose[basis];
  names = tbas[[1]];
  (*--- check for suppression of indices ---*)
  
  test = Map[StringSplit[#, ","] &, Transpose[lat][[1]]] // Flatten //
      Tally // Transpose;
  test1 = Map[StringTake[#, {1, 2}] &, test[[1]]] // Tally;
  If[Length[test1] == 1,
        mylog = True,
        mylog = False
   ];
  If[verb,
        Print[test];
       Print[test1];
        If[Length[test1] == 1,
    test1 = test1 // Flatten;
             
    Print[test1[[2]], " equivalent ", test1[[1]], 
     " atoms in basis, indices at parameters suppressed"],
             None
       ],
   None];
  am = tbas[[3]] /. mom;
  am1 = am /. {-1 -> 0};
  nam = Map[Length[#] &, am1];
  ndim = Plus @@ Map[(2 # + 1) &, am1 // Flatten];
  latt = lat; sel = Transpose[latt][[1]];
  ldl = {};
  Do[ldl = Append[ldl, Map[(2 # + 1) &, am1[[i]]]], {i, 1, nbas}];
  dbl = Map[Plus @@ # &, ldl];
  ham = Table[0, {ndim}, {ndim}];
  (*---calculate all matrix elements---*)Do[atom1 = names[[is]];
   Do[atom2 = names[[js]];
    Do[l1 = am[[is, il]];
     If[l1 == -1, l1s = 0, l1s = l1];
     Do[l2 = am[[js, jl]];
      If[l2 == -1, l2s = 0, l2s = l2];
      Do[Do[(*---calculate the position of the matrix element---*)
   
             ro = 
         Sum[dbl[[k]], {k, 1, is - 1}] + 
          Sum[ldl[[is, k]], {k, 1, il - 1}] + l1s + im + 1;
        co = 
         Sum[dbl[[k]], {k, 1, js - 1}] + 
          Sum[ldl[[js, k]], {k, 1, jl - 1}] + l2s + jm + 1;
        (*---deccription of basis---*)
        
        bas1 = StringJoin[atom1, ",", atom2];
        bas2 = StringJoin[atom2, ",", atom1];
        If[mylog,
              bas3 = 0;
              bas4 = 0,
              bas3 = bas1;
              bas4 = bas2
         ];
        (*---lattice information---*)
         
        p1 = Flatten[Position[sel, bas1]];
        If[p1 == {},
             p = Flatten[Position[sel, bas2]][[1]], p = p1[[1]]];
              dmax = latt[[p, 2]]; latt1 = latt[[p, 3]];
        If[ro < co,
         (*---add the contributions from the shells (nondiagonal \
elements)---*)
         su = 0;
         Do[
               
          su = su + 
            GTTbHamiltonianElement[l1, im, l2, jm, i, latt1[[i]], 
             GOTbBasis -> bas3], {i, 1, dmax}];
                 ham[[ro, co]] = su;
                 ham[[co, ro]] = ComplexExpand[Conjugate[su]],
                 None
          ];
        (*---add contributions from shells+on-
        site (diagonal elements)---*)
        
        If[ro == co, 
           bas2 = atom1; su = 0;
          Do[
          	If[mylog,
          	   onsite=bas3,
          	   onsite=bas1
          	   ];	
          su = su + 
            GTTbHamiltonianElement[l1, im, l2, jm, i, latt1[[i]], 
             GOTbBasis -> onsite], {i, 1, dmax}];
             If[mylog,
          	   onsite=bas4,
          	   onsite=bas1
          	   ];	
         su = su + GTTbSymbol2C[l1, l2, -1, 1, GOTbBasis -> onsite];
         ham[[ro, co]] = su, None];
        If[verb,
         Print["angular momenta : (", jl, ",", jm, ")  (", il, ",", 
          im, ")"];
         Print["basis atoms ", js, " ", is],
         None]
                ,
        {jm, -l2s, l2s}], {im, -l1s, l1s}], {jl, 1, nam[[js]]}], {il, 
      1, nam[[is]]}], {js, 1, nbas}], {is, 1, nbas}];
  Return[ham]]
*)

(* 24.12. hcp Test*)
TbHamiltonianBasis[basis_, lat_] := 
     Module[{is, nbas, am, nam, ndim, ldl, i, dbl, ham, atom1, atom2, js, l1, l2, il, jl, ro, co, k, im, jm, bas1, 
     	     dmax, latt1, su, bas2, mom, am1, l1s, l2s, latt,tbas,names,sel,p1},
     nbas = Length[basis];
     mom = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4};
     tbas = Transpose[basis];
     names = tbas[[1]];
     am = tbas[[3]] /. mom;
     am1 = am /. {-1 -> 0};
     nam = Map[Length[#] &, am1];
     ndim = Plus @@ Map[(2 # + 1) &, am1 // Flatten];
     latt = lat; sel = Transpose[latt][[1]];
     ldl = {};
     Do[
     	ldl = Append[ldl, Map[(2 # + 1) &, am1[[i]]]]
     , {i, 1, nbas}];
     dbl = Map[Plus @@ # &, ldl];
     ham = Table[0, {ndim}, {ndim}];
     (*--- calculate all matrix elements ---*)
     Do[
     	atom1 = names[[is]];
        Do[
           atom2 = names[[js]];
           Do[
           	  l1 = am[[is, il]]; 
           	  If[l1 == -1, 
           	  	 l1s = 0, 
           	  	 l1s = l1
           	  ];
              Do[
              	 l2 = am[[js, jl]]; 
              	 If[l2 == -1, 
              	 	l2s = 0, 
              	 	l2s = l2
              	 ];
                 Do[
                 	Do[
     (*--- calulate the position of the matrix element ---*)
                       ro = Sum[dbl[[k]], {k, 1, is - 1}] + Sum[ldl[[is, k]], {k, 1, il - 1}] + l1s + im + 1;
                       co = Sum[dbl[[k]], {k, 1, js - 1}] + Sum[ldl[[js, k]], {k, 1, jl - 1}] + l2s + jm + 1;
     (*--- deccription of basis ---*)
                       bas1 = StringJoin[atom1, ",", atom2]; 
                       bas2 = StringJoin[atom2, ",", atom1];
     (*--- lattice information ---*)
                       p1 = Flatten[Position[sel, bas1]]; 
                       If[p1 == {}, 
                       	  p = Flatten[Position[sel, bas2]][[1]], 
                          p = p1[[1]]
                       ];
                       dmax = latt[[p, 2]]; latt1 = latt[[p, 3]];
                       If[ro < co,
     (*--- add the contributions from the shells (nondiagonal elements) ---*)
                          su = 0;                   
                          Do[
                             su = su + GTTbHamiltonianElement[l1, im, l2, jm, i, latt1[[i]],GOTbBasis -> bas1]
                          , {i, 1, dmax}];
                          ham[[ro, co]] = su;
                          ham[[co, ro]] = ComplexExpand[Conjugate[su]], 
                          None
                       ];
     (*--- add contributions from shells+on- site (diagonal elements) ---*)
                       If[ro == co, 
                       	  bas2 = atom1; su = 0;
                          Do[
                             su = su + GTTbHamiltonianElement[l1, im, l2, jm, i, latt1[[i]], GOTbBasis -> bas1]
                          , {i, 1, dmax}];
                          su = su + GTTbSymbol2C[l1, l2, -1, 1, GOTbBasis -> bas2];
                          ham[[ro, co]] = su, 
                          None
                       ]
                   , {jm, -l2s, l2s}]
                 , {im, -l1s, l1s}]
              , {jl, 1, nam[[js]]}]
           , {il, 1, nam[[is]]}]
        , {js, 1, nbas}]
     , {is, 1, nbas}];
     Return[ham]
]
(* *)
 
TbHamiltonianTrueBasis[basis_, lat_] := 
     Module[{is, nbas, am, nam, ndim, ldl, i, dbl, ham, js, l1, l2, il, jl, ro, co, k, im, jm, bas1, dmax, latt1, 
     	     su, bas2, mom, am1, l1s, l2s, latt,tbas, mult, n1, n2, aml, nat, names, atom1s, atom2s, p, p1, atom1n,
     	     atom2n,sel,j,bas1n},
     mom = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4};
     nbas = Length[basis];
     tbas = Transpose[basis];
     (*--- extension of basis ---*)
     mult = tbas[[2]]; n1 = {}; n2 = {}; aml = {};
     nat = Plus @@ mult;
     Do[
        Do[
           n1 = Append[n1, tbas[[1, i]]];      
           n2 = Append[n2, StringJoin[tbas[[1, i]] <> ToString[j]]];
           aml = Append[aml, tbas[[3, i]]]
        , {j, 1, mult[[i]]}]
     , {i, 1, nbas}];
     names = {n2, n1};
     am = aml /. mom;
     am1 = am /. {-1 -> 0};
     nbas = nat;
     nam = Map[Length[#] &, am1];
     ndim = Plus @@ Map[(2 # + 1) &, am1 // Flatten];
     latt = lat[[1]]; sel = lat[[2]];
     ldl = {};
     Do[
     	ldl = Append[ldl, Map[(2 # + 1) &, am1[[i]]]]
     , {i, 1, nbas}];
     dbl = Map[Plus @@ # &, ldl];
     ham = Table[0, {ndim}, {ndim}];
     (*--- calculate all matrix elements ---*)
     Do[
        atom1n = names[[2, is]]; atom1s = names[[1, is]];
        Do[
           atom2n = names[[2, js]]; atom2s = names[[1, js]];
           Do[
              l1 = am[[is, il]]; 
              If[l1 == -1, 
           	     l1s = 0, 
           	     l1s = l1
               ];
               Do[
           	      l2 = am[[js, jl]]; 
           	      If[l2 == -1, 
           	  	     l2s = 0, 
           	  	     l2s = l2
           	      ];
                  Do[
              	     Do[
     (*--- calulate the position of the matrix element ---*)
                        ro = Sum[dbl[[k]], {k, 1, is - 1}] + Sum[ldl[[is, k]], {k, 1, il - 1}] + l1s + im + 1;
                        co = Sum[dbl[[k]], {k, 1, js - 1}] + Sum[ldl[[js, k]], {k, 1, jl - 1}] + l2s + jm + 1;
     (*--- deccription of basis ---*)
                        bas1 = StringJoin[atom1s, ",", atom2s]; 
                        bas2 = StringJoin[atom2s, ",", atom1s]; 
                        bas1n = StringJoin[atom1n, ",", atom2n];
     (*--- lattice information ---*)
                        p1 = Flatten[Position[sel, bas1]]; 
                        If[p1 == {}, 
                           p = Flatten[Position[sel, bas2]][[1]], 
                           p = p1[[1]]
                        ];
                        dmax = latt[[p, 2]]; latt1 = latt[[p, 3]];
                        If[ro < co,
     (*--- add the contributions from the shells (nondiagonal elements) ---*)
                           su = 0;
     (*--- no suberscripts if basis, but only one sort ---*)                    
                          If[Length[basis]
                          	==1,bas1n=0,None];
     (* *) 
                           Do[
                              su = su + GTTbHamiltonianElement[l1, im, l2, jm, i, latt1[[i]], GOTbBasis -> bas1n]
                           , {i, 1, dmax}];
                           ham[[ro, co]] = su;
                           ham[[co, ro]] = ComplexExpand[Conjugate[su]], 
                           None
                        ];
     (*--- add contributions from shells+on- site (diagonal elements) ---*)
                        If[ro == co, 
                           bas2 = atom1n; su = 0;
     (*--- no suberscripts if basis, but only one sort ---*)                    
                          If[Length[basis]
                          	==1,bas1n=0;bas2=0,None];
     (* *) 
                           Do[
                              su = su + GTTbHamiltonianElement[l1, im, l2, jm, i, latt1[[i]], GOTbBasis -> bas1n]
                           , {i, 1, dmax}];
                           su = su + GTTbSymbol2C[l1, l2, -1, 1, GOTbBasis -> bas2];
                           ham[[ro, co]] = su, 
                           None
                        ]
                     , {jm, -l2s, l2s}]
                  , {im, -l1s, l1s}]
               , {jl, 1, nam[[js]]}]
            , {il, 1, nam[[is]]}]
         , {js, 1, nbas}]
     , {is, 1, nbas}];
     Return[ham]
]







(*
***)


(****t* /GTHamiltonianPlot
! NAME
!  GTHamiltonianPlot
! AUTHOR
!  W. Hergert
! PACKAGE
!   TightBinding.m
! MODIFICATION HISTORY
!  * 05.12.2014 : first version
!  * 01.12.2016 : BackgroudColor from Red to GTBackgroundColor1 -> a change is psossible
!  * 15.12.2016 : all boxes made of same size
!  * 25.06.2018 : check header and documentation
! USAGE
!  GTHamiltonianPlot[hamiltonian,basis] plots the structure of a hamiltionian, using information about the basis employed in the construction of the Hamiltonian.
! INPUT
!   * ham - Hamiltonian
!   * bas - Basis (in the form used for the construction of the Hamiltonian)
! OUTPUT
!   Plot of the structure of the Hamiltonian
! GTPack OPTIONS
!   -
! STANDARD OPTIONS   
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  This type of plot is very helpful during the development of a Hamiltonian. If the  TB parameters should be found by means of fitting at symmetry points,
!  the structure of the Hamiltonian can be easily inspected.
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
!
*) 

GTHamiltonianPlot[ham_, bas_,OptionsPattern[]] := 
 Module[{lham,defbez, bez, hamp, ang, k, la, l, ll, pos, s, m, i, blist, mom, j,tbas,angmom,mult,nb,nat,size}, 
 	(*--- options ---*)
 	size= OptionValue[ItemSize];
   lham = Length[ham] + 1;
   nb = Length[bas];
   mom = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4};
   defbez = {"s", "\!\(\*SubscriptBox[\(p\), \(y\)]\)", 
    "\!\(\*SubscriptBox[\(p\), \(z\)]\)", 
    "\!\(\*SubscriptBox[\(p\), \(x\)]\)", 
    "\!\(\*SubscriptBox[\(d\), \(xy\)]\)", 
    "\!\(\*SubscriptBox[\(d\), \(yz\)]\)", 
    "\!\(\*SubscriptBox[\(d\), \(3 \*SuperscriptBox[\(z\),            \
         \(2\)]\)]\)", "\!\(\*SubscriptBox[\(d\), \(xz\)]\)", 
    "\!\(\*SubscriptBox[\(d\), \(\*SuperscriptBox[\(x\), \(2\)] -     \
                    \*SuperscriptBox[\(y\), \(2\)]\)]\)"};
   bez = {};
   tbas = Transpose[bas];
   mult = tbas[[2]]; angmom = {};
   nat = Plus @@ mult;
   If[nat == nb,
      angmom = tbas[[3]],
      Do[
         Do[
            angmom = Append[angmom, tbas[[3, i]]]
         , {j, 1, mult[[i]]}]
       , {i, 1, nb}]
   ];
   hamp = Table[" ", {lham}, {lham}];
   Do[
      ang = angmom[[k]] /. mom; la = Length[ang];
      Do[
         l = ang[[ll]];
         pos = Sum[2*s + 1, {s, 0, l - 1}];
         If[ang[[ll]] == -1,                 
            bez = Append[bez, "\!\(\*SuperscriptBox[\(s\), \(\[Star]\)]\)"],
            None
         ];
         Do[
            bez = Append[bez, defbez[[pos + m]]]
         , {m, 1, 2 l + 1}]
       , {ll, 1, la}]
   , {k, 1, nat}]; 
   Do[
   	  hamp[[1, i]] = bez[[i - 1]];
      hamp[[i, 1]] = bez[[i - 1]]
   , {i, 2, lham}];
   blist = {};
   Do[
   	  Do[
   	  	 If[ham[[i - 1, j - 1]] === 0 || ham[[i - 1, j - 1]] === 0., 
            None, 
            blist = Append[blist, {{i, j} -> GTBackGroundColor1}]
         ]
       , {j, 2, lham}]
   , {i,2, lham}];
   blist = Flatten[blist, 1];
   Grid[hamp, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, Background -> {None, None, blist},ItemSize -> size]]

  
(*
***)

(****t* /TbCompose
! NAME
!  TbCompose
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  1st version August 2013
! USAGE
! TbCpompose[ham,pt,pos,ind] places or extracts the matrix pt in the Hamiltonian ham at a certain position pos.
! INPUT
!  * ham    - a given matrix, may be the Hamiltonian without SO
!  * pt     - matrix to insert
!  * pos    - position 
!
!         o    ind = -1 : {{row1,row2},{co1,co2}} dimension of the rectangular matrix to take
!         o    ind >= 0 : {row,column} the place where pt[[1,1]] is inserted
!  * ind    - index: 
!         o   -1 extract a matrix (pt not needed)
!         o    0 insert pt in ham
!         o    1 insert the sum of pt and ham in ham at the corresponding place 
! OUTPUT
!  * ind = -1 : exctracted matrix
!  * ind >= 0 : the modified matrix ham
! OPTIONS
!  - 
! ERROR MESSAGES
!  -
! GTPack MODULES
!  - 
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  To set up a Hamiltonian including the spin-orbit coupling, the Hamiltonian 
!  without spin will be constructed first. The inclusion of SO is than a new arrangement of
!  blocks of that Hamiltonian together with the SO matrices.
!  The module allows to extract parts of a matrix and to put in another matrix
!  at a certain place.
! LITERATURE
!  -
! TODO
!  -
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)

TbCompose[ham_,pt_,pos_,ind_]:=Module[{new,lzei,lsp,is,js,i,j},
     If[ind==-1,new=Take[ham,pos[[1]],pos[[2]]],None];
     If[ind>=0,
     	is=pos[[1]];js=pos[[2]];new=ham;lzei=Length[pt];lsp=Length[pt[[1]]];
     	Do[
     	   Do[
     	   	  new[[i-1+is,j-1+js]]=pt[[i,j]]+ind*ham[[i+is-1,j+js-1]]
           ,{i,1,lzei}]
        ,{j,1,lsp}]
        ,None
      ];
      Return[new]
]

(*
***)


(****t* /GTTbSpinOrbit
! NAME
!  GTTbSpinOrbit
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 01.09.2013 : first version
!  * 11.01.2017 : several l per atom possible, new structure of command
!  * 15.02.2017 : bug fix list soi
!  * 15.04.2017 : bug fix list soi, I should have fixed it 2 months ago?!
!  * 25.06.2018 : check header and documentation
! USAGE
!  GTTbSpinOrbit[hamiltonian,soi] adds spin-orbit coupling to a given tight-binding hamiltonian due to a specified spin-orbit interaction.
! INPUT
! 
!  ham  -  a list of matrices, either one or two. In this form different Hamltonians for spin up and down can be used,
!  information about SOI  - {atom1, atom2, ..}
!		o atom1 -> {"Ga",{{l1,pos1,xi1},{l2,pos2,xi2}}
!		o l1,l2,.. angular momenta
!		o pos1,pos2,.. position of the angular momentum block
!		o xi1,xi2,.. SOI parameters       
! OUTPUT
!  Hamiltonian with spin-orbit interaction included
! GTPack OPTIONS
!  * GOVerbose:
!
!     - True   - additional information
!     - False  - no additional information (Standard)
! STANDARD OPTIONS
!  -
! GTPack MODULES
!   TbCompose,
!   GTTbSpinMatrix
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)

 GTTbSpinOrbit[ham_, soi_,OptionsPattern[]] :=Module[
 	{verb,ndim,ndim2,ham0,smat,natoms,blk,nso,k,lv,pos,xi,ik,sup,sdn,sud,sdu,blist,hamp},
  (*--- options ---*)
     verb = OptionValue[GOVerbose];
  (*--- double the Hamiltonian ---*)
     If[Length[ham] == 2,
        If[verb,
           Print["Different Hamiltonians for the spin directions."],
           None
        ];
        ndim  = Length[ham[[1]]];
        ndim2 = Length[ham[[2]]];
        If[ndim == ndim2,
           None,         
           Print["Error: Hamiltonians for the two spin directions have different size!"]; 
           Abort[]
        ];
        ham0 = Table[0, {2*ndim}, {2*ndim}];
        smat = Table[0, {2*ndim}, {2*ndim}];
        ham0[[1 ;; ndim, 1 ;; ndim]] = ham[[1]];
        ham0[[1 + ndim ;; 2*ndim, 1 + ndim ;; 2*ndim]] = ham[[2]],
        If[Length[ham] == 1,
           ndim  = Length[ham[[1]]];
           ham0  = Table[0, {2*ndim}, {2*ndim}];
           smat  = Table[0, {2*ndim}, {2*ndim}];
           ham0[[1 ;; ndim, 1 ;; ndim]] = ham[[1]];
           ham0[[1 + ndim ;; 2*ndim, 1 + ndim ;; 2*ndim]] = ham[[1]],
           Print["Error: Input has to be list of form {Hamiltonian}."]; Abort[]
        ]
     ];
  (*--- print information about SOC interaction ---*)
      If[verb,
         Print[ Grid[soi, Frame -> All, Background -> {1 -> GTBackGroundColor1, None}]],
         None
      ];
  (*--- construct SOC matrix ---*)
      natoms = Length[soi];
      Do[
         blk = soi[[k, 2]];
         nso = Length[blk];
         Do[
            lv  = blk[[ik,1]];
            pos = blk[[ik,2]];
            xi  = blk[[ik,3]];
    (*--- diagonal blocks M---*)
            sup = xi*GTTbSpinMatrix[lv, 1, 1];
            sdn = xi*GTTbSpinMatrix[lv, -1, -1];
            smat[[pos ;; pos + 2*lv, pos ;; pos + 2*lv]] = sup;
            smat[[ndim + pos ;; ndim + pos + 2*lv, ndim + pos ;; ndim + pos + 2*lv]] = sdn;
    (*---non-digonal blocks N---*)
            sud = xi*GTTbSpinMatrix[lv, 1, -1];
            sdu = xi*GTTbSpinMatrix[lv, -1, 1];
            smat[[pos ;; pos + 2*lv, ndim + pos ;; ndim + pos + 2*lv]] = sud;
            smat[[ndim + pos ;; ndim + pos + 2*lv, pos ;; pos + 2*lv]] = sdu;
        , {ik, 1, nso}]
      , {k, 1, natoms}];
   (*-- print structure of SOC matrix ---*)
      If[verb,
      	Print["Structure of SOC matrix."];
         blist = {}; 
         hamp  = Table[" ", {2*ndim}, {2*ndim}];
         Do[
            Do[
               If[smat[[i, j]] === 0 || smat[[i, j]] === 0., 
               	  None, 
                  blist = Append[blist, {{i, j} -> GTBackGroundColor1}]
                  ]
            , {j, 1, 2*ndim}]
         , {i, 1, 2*ndim}];
         blist = Flatten[blist, 1];
        Print[Grid[hamp, Frame -> All, Background -> {None, None, blist}, ItemSize -> All]],
        None
      ]  ;
      Return[ham0+smat]
  ]


(*
***)

(****t* /HydrogenRadial
! NAME
! HydrogenRadial
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  1st version October 2013
! USAGE
!  Module calculates the radial part of a Hydrogen like wave function
!  
! INPUT
!   * Z   - atomic number
!   * n   - main quantum number
!   * l   - orbital momentum quantum number
!   * r   - distance from nucleus
!
! OUTPUT
!  value of the radial wave function
!  OPTIONS
!  -
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  This module is an internal module of the package TightBinding.m
! LITERATURE
!   P.Rennert, A. Chasse, w. Hergert, Einfuehrung in die Quantenphysik, Springer
!   
! TODO
!  -
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)

HydrogenRadial[Z_,n_,l_,r_]:=Module[{x},x=2 Z r/n; Z^(3/2) 2/n^2  Sqrt[(n-l-1)!/(n+l)!] Exp[-Z r/n] x^l LaguerreL[n-l-1,2l+1,x]]

(*
***)

(****t* /GTTbAtomicWaveFunction
! NAME
!  GTTbAtomicWaveFunction
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 1st version October 2013
!  * revised version November 2013
!  * Tesseral Harmonics Dezember 2013
! USAGE
!  GTTbAtomicWaveFunction[oz,n,l,m,rvec] gives the value of an atomic-like wave function of the atomic number to the quantum numbers (n,l,m) at position rvec
!
! INPUT
!   * oz  - atomic number
!   * n   - main quantum number
!   * l   - orbital momentum quantum number
!   * m   - magnetic quantum number
!   * r   - list of space coordinates (vector relative to center of the atom)
!
! OPTIONS
!   GOHarmonics  
!
!	o "Complex"  -> complex spherical harmonics 
!	o "Real"     -> real linear combinations
! OUTPUT
!  value of the atomiclike wave function
!  
! ERROR MESSAGES
!  "Error: Wrong type of Spherical harmonics!" if type of spherical harmonics is not "Complex" or "Real"
!
! GTPack MODULES
!   Hyrogenradial
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  If the type of the angular Part is "Complex" the Mathematica funtions SphericalHarmonicY are used,i.e.
!  for a given l, m runs from -l to l.
!  If the  type of the angular part is "Real" the real linear combinations are formed.
!  	m=0  Y_{l,0}(theta,phi)
!  	m>0  (Y_{l,m}(theta,phi)-Y_{l,-m}(theta,phi))/Sqrt[2]
!  	m<0  (Y_{l,-m}(theta,phi)-Y_{l,m}(theta,phi))/Sqrt[2]/I
! LITERATURE
!  see: Glaeske, Reinhold, Volkmer, Quantenchenmie V
! TODO
!  Check the definition of the real linear combinations. It gives now real linear combinations,
!  but doesn't agree with the literature.
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTTbAtomicWaveFunction[oz_,n_,l_,m_,rvec_,OptionsPattern[]]:=Module[{sph,x,y,z,r,z1,ang,rp,wf},
                   sph=OptionValue[GOHarmonics];
	               x=rvec[[1]];y=rvec[[2]];z=rvec[[3]];r=Sqrt[x^2+y^2+z^2];z1=Sqrt[x^2+y^2];
	               If[z1==0,
	               	  \[Theta]=0,
	               	  \[Theta]=ArcTan[z,z1]
	               ];
	               If[y==0,
	               	  \[Phi]=0,
	               	  \[Phi]=ArcTan[x,y]
	               ];
                   If[sph=="Complex",
                   	  			ang=SphericalHarmonicY[l,m,\[Theta],\[Phi]],
                       		If[sph=="Real",
                      			ang=GTTesseralHarmonicY[l,m,\[Theta],\[Phi]]
                        ,None
                     ]
                     ,Print["Error: Wrong type of Spherical harmonics!"]
                  ];
	              rp=HydrogenRadial[oz,n,l,r];
	              wf=rp*ang//N//Chop;
	              Return[wf]
]
(*
***)

(****t* /GTTbBlochFunction
! NAME
!  GTTbBlochFunction
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
! * 1st version November 2013
! * 5.10.2014 : fit to changes made in TB part in last year
! USAGE
!  The module provides the value of a Bloch function, constructed from atomic-like function to
!  a given wave vector kvec at the position rvec.
! INPUT
!   * kvec    - k-point under consideration
!   * rvec    - space coordinates
!   * pos     - positions of atoms in the problem (output of GTShells)
!   * atom     - information about the atom {symbol,n,l,m,,ns,tau,scale}
!
!             1   : symbol - chemical symbol
!             2-4 : n,l,m  - quantum numbers
!             5   : ns     - number of shells to include
!             6   : tau    - position in unit cell
!             7   : scale  -scaling factor for atomic number
!
! OPTIONS
!   GOHarmonics 
!	o "Complex"  -> complex spherical harmonics 
!	o "Real"     -> real linear combinations
!
! OUTPUT
!  value of the Bloch function
!  
! ERROR MESSAGES
!  -
! GTPack MODULES
!   GTTbAtomicWaveFunction
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!   up to now it is not tested, if the module works also in case of more atoms per unit cell.
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
*)

GTTbBlochFunction[kvec_,rvec_,pos_,atom_,OptionsPattern[]]:=Module[
	           {opt,asym,z,n,m,l,\[Tau],ntype,nshell,apos,shells,s1,s2,nl,asym1,i},
               (*--- options ---*)
               opt=OptionValue[GOHarmonics];
               (*--- set atomic data ---*)
               asym=atom[[1]];z=ElementData[asym,"AtomicNumber"]*atom[[7]];
               n=atom[[2]];l=atom[[3]];m=atom[[4]];shells=atom[[5]];\[Tau]=atom[[6]];
               (*--- positions to type asym ---*)
               ntype=Length[pos];
               asym1 = asym <> "," <> asym;
               Do[
                  If[pos[[i,1]]==asym1,
                     nshell=pos[[i,2]];apos=pos[[i,3]],
                     None
                  ]
               ,{i,1,ntype}];
               apos=Flatten[Take[apos,{1,shells}],1];
               (*--- append the central atom position ---*)
               apos=Append[apos,\[Tau]];nl=Length[apos];
               (*--- evaluate the Bloch function ---*) 
               s1=Map[Exp[I 2 \[Pi] kvec.#]&,apos];
               s2=Map[(rvec-\[Tau]-#)&,apos];
               s2=Map[GTTbAtomicWaveFunction[z,n,l,m,1.0*#,GOHarmonics->opt]&,s2];
               Return[s1.s2/Sqrt[nl]//Chop]
]  
          
(*
***)


(****t* /GTTbWaveFunction
! NAME
!  GTTbWaveFunction
! AUTHOR
!  W. Hergert
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * November 2013, 1st version 
!  * October 2014, spin is included 
!
! USAGE
!  The module calculates the wave function to a TB Hamiltonian with model atomic functions
!  for symmetry analysis
! INPUT
!  * file    - file name
!  * kpt     - k-point
!  * bdnr    - band number
!  * basis   - information about the basis set  {a1orb1,a1orb2,a1orbn1,a2orb1,...}  aiorbj correspond to the information used in GTTbBlochFunction
!  * pos     - information about all position as it comes from GTShells
!  * geo     - information about the region in which the data are calculated {xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz}. For example,
!            if the xy-plane should be calculated use: {xmin,xmax,nx,ymin,ymax,ny,0,0,1}
!  * spin
! OUTPUT
! List of values of the wavefunction
!  
! ERROR MESSAGES
!  -
! GTPack MODULES
!   GTTbBlochFunction
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
*)

GTTbWaveFunction[file_,kpt_,bdnr_,basis_,pos_,geo_,spin_:1]:=Module[
	                     {wav,nkpt,nx,xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz,ny,nz,field,ka,
	                      rvec,cvecs,i,j,k,l,nbas,s,atom,bf},
	                     (*--- read eigenvectors ---*)
                         wav=GTReadFromFile[file];nkpt=Length[wav];
                         Do[
                            If[wav[[i,1]]==kpt,
                               cvecs=wav[[i,2]],
                               None
                              ] 
                         ,{i,1,nkpt}];
                         (*--- check basis ---*)
                         nbas=Length[basis];
                         If[nbas==Length[cvecs[[1]]]/spin,
                         	None,
                         	Print["Error : basis and eigenvectors not compatible!"];Abort[]
                         ];
                         (*--- set mesh data ---*)
                         nx=geo[[3]];ny=geo[[6]];nz=geo[[9]];
                         xmax=geo[[1]];xmin=geo[[2]];
                         If[nx==1,
                         	dx=0,
                         	dx=(xmax-xmin)/(nx-1)
                         ];
                         ymax=geo[[4]];ymin=geo[[5]];
                         If[ny==1,
                         	dy=0,
                         	dy=(ymax-ymin)/(ny-1)
                         ];
                         zmax=geo[[7]];zmin=geo[[8]];
                         If[nz==1,
                         	dz=0,
                         	dz=(zmax-zmin)/(nz-1)
                         ];   
                         field={};ka=kpt;
                         If[spin==1,
                         (*--- Without Spin ---*)
	                         Do[             
	                            Do[
	                              Do[
	                                 rvec={xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(l-1)*dz}; 
	                                 s=0;         
	                                 Do[If[cvecs[[bdnr,k]]==0, 
	                                       None,
	                                       atom=basis[[k]];
	                                       bf=GTTbBlochFunction[ka,rvec,pos,atom,GOHarmonics->"Real"];
	                               	       s=s+cvecs[[bdnr,k]]*bf
	                               	    ]
	                                 ,{k,1,nbas}];
	                                 field=Append[field,s];
	                               ,{i,1,nx}]
	                            ,{j,1,ny}]
	                         ,{l,1,nz}];,
	                      (*--- Spin ---*)
	                         Do[             
	                            Do[
	                              Do[
	                                 rvec={xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(l-1)*dz}; 
	                                 s=0;         
	                                 Do[If[Abs[cvecs[[bdnr,k]]]+Abs[cvecs[[bdnr,k+nbas]]]==0, 
	                                       None,
	                                       atom=basis[[k]];
	                                       bf=GTTbBlochFunction[ka,rvec,pos,atom,GOHarmonics->"Real"];
	                               	       s=s+cvecs[[bdnr,k]]*{bf,0}+cvecs[[bdnr,k+nbas]]*{0,bf}
	                               	    ]
	                                 ,{k,1,nbas}];
	                                 field=Append[field,s];
	                               ,{i,1,nx}]
	                            ,{j,1,ny}]
	                         ,{l,1,nz}];
                         ];
                         Return[field]
]


(*
***)


(****t* /GTTbParameterNames
! NAME
!  GTTbParameterNames
! AUTHOR
!  S. Thomas
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 15.12.2013 : first version
!  * 15.01.2015 : GOVerbose is implemented, Extension to input of parameter sets 
!  * 25.06.2016 : check header and documentation
! USAGE
!  GTTbParameterNames[maximum angular momentum,maximum shell number] creates a set of parameter names up to maximum angular momentum and maximum shell number.
!
! INPUT
!  * lmax       - is either the maximum angular momentum or a parameter set from a database
!  * dmax       - maximum shell number
!
!  OUTPUT
!  a list of parameter names
!  
! GTPack OPTIONS
!  * GOTbBasis:
!
!      - 0       - no basis  (Standard)
!      - "Pb,Te" -basis for PbTe 
!  * GOVerbose:
!
!       - False  - no additional output (standard)
!       - True   - output of a table which shows to which parameter the Tb parameter
!                  will belong in the FORTRAN code.
!  
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTTbSymbol2C
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The module generates a a set of parameter names WITHOUT the interactive input mode.
!  This is especially helpful for the generation of substitution rules.
! LITERATURE
!  -
! TODO
!  It is not optimal in the moment. It could be better to have a list of angula momenta like
!  {"s*","s","p"} in the input. If one would like to construct a d Hamiltonian only, lmax=2 
!  will create a lot of unnecessary output.
!  More serious is, that "s*" can not be handled by an input of lmax only.
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)


GTTbParameterNames[lmax_, dmax_:1, OptionsPattern[]] :=
          Module[{p2c, ll, lu, m, i, symb, tab, bas, tab1, k, ln, ind, j, bas1, ln1, index1, lt, verb,pl,mom}, 
          	 p2c = {};mom = {"s*", "s", "p", "d", "f", "g"};
             bas = OptionValue[GOTbBasis];
             verb=OptionValue[GOVerbose];
             (*--- check wether list of parameters or lmax ---*)
             If[Head[lmax]===List, 
                ln = Length[lmax];
                Do[p2c = Append[p2c, lmax[[i, 1]]], {i, 1, ln}],
         (*--- check if basis,list of basis atoms ---*)
                If[bas == 0, 
                   ind = 0; ln = 1; ln1 = 1; index1 = {0}, 
                   tab = StringSplit[bas, ","]; ln = Length[tab]; 
                   tab1 = {}; ind = 1; index1 = {};
          (*--- create a list of combinations of basis atoms ---*)
                   Do[
                      Do[tab1 = Append[tab1, StringJoin[{tab[[i]], ",", tab[[j]]}]];
                         If[i == j, 
                       	    index1 = Append[index1, 0], 
                            index1 = Append[index1, 1]
                         ]
                      ,{i, 1, j}] 
                   ,{j, 1, ln}]; 
                   ln1 = Length[tab1]; tab = tab // Flatten
                 ];
           (*--- if basis go over all combinations, otherwise olny over all paramters ---*)
                 Do[
                    Do[
                       Do[
           (*--- decide if (sps) and (pss) are necessary ---*)
                            If[index1[[k]] == 0, 
                              lt = lu, 
                              lt = lmax
                            ];
                            Do[
                            	If[ind == 1, 
                            	   bas1 = tab1[[k]], 
                            	   bas1 = 0
                            	 ];
                                 Do[symb = GTTbSymbol2C[ll, lu, m, i, GOTbBasis -> bas1];
                                    p2c = Append[p2c, symb] 
                                 ,{m, 0, Min[ll, lu]}]
                            , {ll, 0, lt}]
                         ,{lu, 0, lmax}]
                      ,{i, 1, dmax}]
                   ,{k, 1, ln1}];
             (*--- append on-site paramters ---*)
                   Do[
                   	  If[ind == 1, 
                   	  	 bas1 = tab[[k]], 
                   	  	 bas1 = 0
                   	  ];
                      symb = GTTbSymbol2C[0, 0, -1, 1, GOTbBasis -> bas1]; 
                      p2c = Append[p2c, symb];
                      symb = GTTbSymbol2C[1, 1, -1, 1, GOTbBasis -> bas1]; 
                      p2c = Append[p2c, symb];
                      symb = GTTbSymbol2C[2, 2, -1, 1, GOTbBasis -> bas1]; 
                      p2c = Append[p2c, symb];
                      symb = GTTbSymbol2C[2, 2, -2, 1, GOTbBasis -> bas1]; 
                      p2c = Append[p2c, symb];
                      symb = GTTbSymbol2C[2, 2, -3, 1, GOTbBasis -> bas1]; 
                      p2c = Append[p2c, symb];
                      symb = GTTbSymbol2C[1, 2, -1, 1, GOTbBasis -> bas1]; 
                      p2c = Append[p2c, symb]
                    ,{k, 1, ln}]
             ];
             If[verb === True, 
                pl = Table["P(" <> ToString[i] <> ")", {i, 1, Length[p2c]}];
                Do[Print[p2c[[i]], "  ->  ", pl[[i]]], {i, 1, Length[p2c]}]; 
                Return[p2c], Return[p2c]
             ]
]

(*
***)


(****t* /GTTbToFortran
! NAME
!  GTTbToFortran
! AUTHOR
!  S. Thomas
! PACKAGE
!  TightBinding.m 
! MODIFICATION HISTORY
!  * 20.12.2013 : 1st version
!  * 13.01.2015 : as a results of the change in GTTbParameterNames lmax can be a parameter set 
!  * 25.06.2018 : check header and documentation             
! USAGE
!  The module is used to transform a tight-binding Hamiltonian to a FORTRAN program
! INPUT
!  * ham        - Hamiltonian
!  * lmax       - maximum angular momentum 
!  * dmax       - maximum shell number
!  * file       - file name for output
!
! OUTPUT
!  Hamitonian will be exported as a FORTRAN module into file
!
! GTPack OPTIONS
!  * GOTbBasis
!  
!	o  0       - no basis
!	o  "Pb,Te" - basis for PbTe 
!    
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTTbParameterNames
!
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! RELEAS
E
!  1.0.0
! PROBLEMS
! - 
!--------------------------------------------------------------------------------
! SOURCE
*)
GTTbToFortran[ham_,lmax_,dmax_,file_,OptionsPattern[]]:=
          Module[{bas,rule,pa,ham1,ham2,ham3,s,i,verb,j,ft},
          bas=OptionValue[GOTbBasis];verb=OptionValue[GOVerbose];
          pa=GTTbParameterNames[lmax,dmax,GOTbBasis->bas];
     (*--- select the acutally used parameters - drops some onsite elements ---*)
          rule=Select[pa,Not[FreeQ[ham,#]]&];
     (*--- construct the substitution rules ---*)
          rule=Table[rule[[i]]->p[i],{i,1,Length[rule]}];
          rule=rule~Join~{\[Xi]->kx,\[Eta]->ky,\[Zeta]->kz};ft=rule;
          ham1=Table[
          	         Table[
          	         	    "H("<>ToString[i]<>","<>ToString[j]<>") = "<>ToString[FortranForm[ham[[i,j]]/.ft/.x_Complex->Re[x]+ci*
          	         	     Im[x](*Replace Imaginary Unit by ci*)/.Table[p[i]->P[i],{i,1,Length[ft]-3}](*Replace p[i] by params[i]*)]]
          	          ,{i,1,Length[ham]}]
          	   ,{j,1,Length[ham]}]//Flatten;
      (*--- Convert Scientific Notation to normal product with exponent ---*)
          ham2=Map[StringReplace[#,RegularExpression["([0-9])e(-?)([0-9]+)"] -> "$1 * 10**($2$3)"]&,ham1];
      (*--- Convert Int to double ---*)
          ham2=Map[StringReplace[#,RegularExpression["([^PH][^0-9,.][0-9]+)([^0-9.])"]->"$1._8$2"]&,ham2];
      (*--- Convert float to double ---*)
          ham2=Map[StringReplace[#,RegularExpression["([.][0-9]*)([^_0-9])"]->"$1_8$2"]&,ham2];
          ham2=Map[StringReplace[#,RegularExpression["([0-9]+)[.]?([0-9]*)$"]->"$1.$2_8"]&,ham2];
      (*--- Manual optimization by replacing cos, sin and sqrt ---*)
          ham3=Map[StringReplace[#,RegularExpression["Cos\\(k(.)\\*Pi\\)"]->"ck$1"]&,ham2];
          ham3=Map[StringReplace[#,RegularExpression["Sin\\(k(.)\\*Pi\\)"]->"sk$1"]&,ham3];
          ham3=Map[StringReplace[#,RegularExpression["Cos\\(\\(k(.)\\*Pi\\)/2\\._8\\)"]->"ck$1d2"]&,ham3];
          ham3=Map[StringReplace[#,RegularExpression["Sin\\(\\(k(.)\\*Pi\\)/2\\._8\\)"]->"sk$1d2"]&,ham3];
          ham3=Map[StringReplace[#,RegularExpression["Sqrt\\(([23])\\._8\\)"]->"s$1"]&,ham3];
      (*--- Output into the file ---*)   
          s=OpenWrite[file];
(*--- Header ---*)
WriteString[s,"
! This module was automatically created with Mathematica!
MODULE Hamiltonian
  IMPLICIT NONE
  !
  COMPLEX(8), PARAMETER :: ci = CMPLX(0,1,KIND=8)
  REAL(8), PARAMETER :: pi = 3.14159265358979323846264_8
  ! REAL(8), PARAMETER :: E = 2.71828182845904523536028_8
  REAL(8), PARAMETER :: s2 = Sqrt(2._8), s3 = Sqrt(3._8)
  !
  CONTAINS
  ! 
  SUBROUTINE hamilton(H,kx,ky,kz,P)"
];
(*--- Variables ---*)
WriteString[s,"
  COMPLEX(8), DIMENSION("<>ToString[Length[ham]]<>","<>ToString[Length[ham]]<>") :: H"
];
WriteString[s,"
  REAL(8), DIMENSION("<>ToString[Length[ft]-3]<>") :: P"
];
WriteString[s,"
  REAL(8) :: kx, ky, kz
  !
  REAL(8) :: ckx, cky, ckz
  REAL(8) :: ckxd2, ckyd2, ckzd2
  REAL(8) :: skx, sky, skz
  REAL(8) :: skxd2, skyd2, skzd2
  !
  ckx = COS(kx*Pi)
  cky = COS(ky*Pi)
  ckz = COS(kz*Pi)
  ! 
  ckxd2 = COS((kx*Pi)/2._8)
  ckyd2 = COS((ky*Pi)/2._8)
  ckzd2 = COS((kz*Pi)/2._8)
  !
  skx = SIN(kx*Pi)
  sky = SIN(ky*Pi)
  skz = SIN(kz*Pi)
  ! 
  skxd2 = SIN((kx*Pi)/2._8)
  skyd2 = SIN((ky*Pi)/2._8)
  skzd2 = SIN((kz*Pi)/2._8)
  !
"
];
Do[
WriteString[s,
"  "<>ham3[[i]]<>"\n"
]
,{i,1,Length[ham3]}];
WriteString[s,"
  END SUBROUTINE
END MODULE"
];
If[verb==True,
   Do[
      Print["  "<>ham3[[i]]]
   ,{i,1,Length[ham3]}],
   None
]  ;
Close[s] 
]

(*
***)



(****t* /GTSymmetryBasisFunctions
! NAME
!  GTSymmetryBasisFunctions
! AUTHOR
!  W. Hergert
! PACKAGE
! TightBinding.m 
! MODIFICATION HISTORY
!  1st version August 2014
!  1. November 2016: Out put changed such, that in case of Verbose -> False the whole table is returned, i.e.
!                    all what is printed otherwise with Print[Grid[]]. This helps if 
! USAGE
!  GTSymmetryBasisFunctions[ct,wfu] calculates to which irreducible representations the functions are basis functions.
!
! INPUT
!   * ct  - character table of the point group (Output of GTCharacterTable)
!   * wfu - List of wave functions          
! OUTPUT
!  Table with the information
! OPTIONS
!  * GOVerbose
!
! 	True  : print of a table (standard) 
!	False : output f a list
!  * GONames   
!
!	 {}  : the names from the character table are used (standard)
!	list : the names in the list are used
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTCharProjectionOperator 
! GTPack NOTEBOOKS 
!  - 
! DESCRIPTION
!  A set of basis functions (polynoms in x,y,z) and a point group are given.
!  It is investigated to which irreducible representation the functions belong
!  as basis functions
! LITERATURE
!  R.F. Egorov et. al. , phys. stat. sol. 26, 391 (1968), p. 396
! TODO
!  -
! PROBLEMS
!  This command can be formulated only as a block. It is not clear to me why it is so.
! SOURCE
!--------------------------------------------------------------------------------
! 
*)
GTSymmetryBasisFunctions[ct_, wfu_, OptionsPattern[]] := Block[
	{classes, chars, ireps, nwf, nirep, tab1, i, j, f, tab, w1, names, prt}, 
    (*--- Preparation of data ---*)
	w1 = wfu; nwf = Length[w1]; classes = ct[[1]]; 
    chars = ct[[2]]; nirep = Length[classes];
    (*--- names ---*)
    names = OptionValue[GONames]; 
    If[Length[names] == 0, 
       ireps = ct[[3]], 
       ireps = names
    ];
    (*--- Calculation of the projections ---*)
    tab = {};
    Do[
       w = w1[[j]] /. {Global`x -> xswf, Global`y -> yswf, Global`z -> zswf};  
       f[xswf_, yswf_, zswf_] = w; tab1 = {};
       Do[
       	  tab1 = Append[tab1, GTCharProjectionOperator[classes, chars[[i]], f, {xswf, yswf, zswf}]]
       , {i, 1, nirep}]; 
       tab = Append[tab, tab1]
    , {j, 1, nwf}]; 
    tab = tab //. {xswf -> Global`x, yswf ->Global`y, zswf -> Global`z};
    (*--- Output of the results ---*)
    prt = OptionValue[GOVerbose];
    w1  = Prepend[w1, " "];
    tab = Join[{ireps}, tab];
    tab = Join[{w1}, Transpose[tab]] // Transpose;
    If[prt == True,
       Print[Grid[tab, Frame -> All, 
                  Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
                  Background -> {{1 -> GTBackGroundColor1}, {1 -> 
                  GTBackGroundColor1}, {1, 1} -> GTCornerColor}
                  ]
             ],
       Return[tab]
    ]
]


(*
***)



(****t* /GTIrepInfo
! NAME
!  GTIrepInfo
! AUTHOR
!  W. Hergert
! PACKAGE
! TightBinding.m 
! MODIFICATION HISTORY
!  1st version January 2014
! USAGE
!  GTIrepInfo[pg,pgct,ireps] exctracts those information from the representation matrices 
!  which are necessary to construct the three-center matrix elements
!
! INPUT
!  * pg - pointgroup
!  * pgct   - character table of a point group
!  * ireps  - List of Ireps
!               
! OUTPUT
!  List of matrix elements of the Ireps
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTGetIrep
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The matrix element appears in the list in correspondance to the position of a
!  symmetry transformation in the point group list
! LITERATURE
!  R.F. Egorov et. al. , phys. stat. sol. 26, 391 (1968), p. 396
! TODO
!  may be this should be a part of RepresentationTheory.m ?
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTIrepInfo[pg_,pgct_, ireps_,OptionsPattern[]] := Module[{tab,nir,ctab,ncl,nirep,repmat,dim,ir,i,j,meth}, tab = {}; nir = Length[ireps];
	  meth=OptionValue[GOMethod];
     (*--- character table ---*)
     ctab = pgct[[2]]; ncl = Length[ctab[[1]]];
     (*--- loop over the IREPs for basis ---*)
     Do[
        nirep = ireps[[ir]]; dim = ctab[[nirep, 1]];
        repmat = GTGetIrep[pg, nirep, pgct,GOMethod->meth];
       Do[  
          Do[      
             tab = Append[tab, {{nirep, i, j}, Map[#[[i, j]] &, repmat]}]
          , {i, 1, dim}]
      , {j, 1, dim}]
     , {ir, 1, nir}]; 
     Return[tab]
]

(*
***)



(****t* /GTTbNumberOfIntegrals
! NAME
!  GTTbNumberOfIntegrals
! AUTHOR
!  W. Hergert
! PACKAGE
! TightBinding.m 
! MODIFICATION HISTORY
!  1st version August 2014
!  01.02.2017 :slightly changed, input of subgroups fitted to other commands, can handle now several
!              Q_p^l per shell
! USAGE
!  GTTbNumberOfIntegrals[grp,sbgrps,ireps] calculates the number of independent tight-binding 
!  integrals for functions transforming like irreducible representations of point group. 
!  The character tables of the groups G_l^p have to be provided.
!
! INPUT
!  * grp - character table of group G_0 of the structure
!  * sbgrps - charcter tables of G_l^p for the different shells 
!  * ireps - list of numbers of Ireps. (Ireps to which the atomic-like 
!            functions are basis functions.) Number is the position of
!            the Irep in the character table.
! OPTIONS
!  * GOVerbose 
!
!	o False (standard) , output of intermediate results suppressed
!	o more information
!  * GONames  
! 
!	o {} (standard) , names of Ireps from charcter table are used  
!	o list of names for the Ireps      
! OUTPUT
!  Tables for the number of independent integrals for the shells
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTIrep
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  R.F. Egorov et. al. , phys. stat. sol. 26, 391 (1968), p. 396
! TODO
!  -
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTTbNumberOfIntegrals[grp_, subgrps_, ireps_, OptionsPattern[]] := Module[
	{verb,notat,nshell,nirep,clgrp,clct,ncl,names,nvec,cls,ncls,s,pos,z,cijs,rep,
	 k,nij,cij,cij1,j,tb,tab},
  (*--- options ---*)
     verb   = OptionValue[GOVerbose];
     notat  = OptionValue[GONames];
  (*---  ---*)
     nshell = Length[subgrps];
     nirep  = Length[ireps];
     tab    = {{"number of shells", "Vectors \!\(\*SubsuperscriptBox[\(Q\), \(p\), \(1\)]\) in shells"}, {nshell, Length[#] & /@ subgrps}} // Transpose;
     Print[Grid[tab, Frame -> All, Background -> {1 -> GTBackGroundColor1}]];
  (*--- Data of the group G0 ---*)
     clgrp  = grp[[1]];
     clct   = grp[[2]];
     ncl    = Length[clgrp];
     If[notat == {},
        names = Map[grp[[3, #]] &, ireps],
        names = notat
     ];
  (*--- loop over all shells ---*)
     Do[
        nvec = Length[subgrps[[s]]];
  (*--- loop over all subgroups G_l^p  in shell ---*)
        Do[
           cls  = subgrps[[s, p, 1]];
           ncls = Length[cls];
           If[verb,
              Print[cls],
              None
           ];
  (*--- Find the classes of the group G0, in which the subgroup classes appear ---*)
           pos = {};
           Do[
              z = Map[Intersection[cls[[i]], #] &, clgrp];
              If[verb,
                 Print[z],
                 None
              ];
              pos = Append[pos, Flatten[Position[z, Sort[cls[[i]]]]][[1]]]
           , {i, 1, ncls}];
           cijs = {};
           Do[
  (*--- charcters from G0 to the classes of the subgroups ---*)
              rep = Map[clct[[ireps[[k]], #]] &, pos];
  (*--- calculation of the cijs ---*)
              cijs = Append[cijs, GTIrep[rep, subgrps[[s, p]], GOVerbose -> False]]
           , {k, 1, nirep}];
  (*--- evaluation of the numbers of the parameters ---*)
           nij = Table[0, {nirep}, {nirep}];
           Do[
              cij = cijs[[k]]; 
              nij[[k, k]] = Plus @@ Map[# (# + 1)/2 &, cij]
           , {k, 1, nirep}];
           Do[
           	  cij = cijs[[k]];
              Do[
              	 cij1        = cijs[[j]];
                 nij[[j, k]] = cij.cij1; 
                 nij[[k, j]] = nij[[j, k]]
              , {j, k + 1, nirep}]
           , {k, 1, nirep}];
           tb  = Join[{Flatten[Join[{ToString[s] <> "/" <> ToString[p]}, names]]}, Transpose[Join[{names}, nij // Transpose]]];
           tab = Grid[tb, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, Frame -> True, 
                         Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, None}, {1, 1} -> GTCornerColor}
                     ];
           Print[tab]
        , {p, 1, nvec}]
     , {s, 1, nshell}]
  ]


(*-- ersetzt 1.2.2017 

 GTTbNumberOfIntegrals[grp_, subgrps_, ireps_, OptionsPattern[]] := Module[
	{nbsg,nirep,clgrp,clct,ncl,gv,gn,names},
  (*--- options ---*)	
     gv = OptionValue[GOVerbose];
     gn = OptionValue[GONames];
  (*---   ---*)    
     nbsg  = Length[subgrps];
     nirep = Length[ireps];
  (*--- Data of the group G0 ---*)
     clgrp = grp[[1]]; 
     clct  = grp[[2]]; 
     ncl   = Length[clgrp];
     If[gn=={},
       names=Map[grp[[3,#]] &, ireps],
       names=gn
     ];  
  
  (*--- loop over all subgroups G_l^p ---*)  
    Do[ 
       cls = subgrps[[s, 1]]; ncls = Length[cls];
       If[gv,
          Print[cls],
          None
       ];
  (*--- Find the classes of the group G0,in which the subgroup classes appear ---*)
       pos = {};
       Do[
          z = Map[Intersection[cls[[i]], #] &, clgrp];
          If[gv,
             Print[z],
             None
          ];
          pos = Append[pos, Flatten[Position[z, Sort[cls[[i]]]]][[1]]]
       , {i, 1, ncls}];
       cijs = {};
       Do[
  (*--- charcters from G0 to the classes of the subgroups ---*)
          rep = Map[clct[[ireps[[k]], #]] &, pos];
  (*--- calculation of the cijs ---*)
          cijs = Append[cijs, GTIrep[rep, subgrps[[s]], GOVerbose -> False]];
       , {k, 1, nirep}];
   (*--- evaluation of the numbers of the parameters ---*)
       nij = Table[0, {nirep}, {nirep}];
       Do[
          cij = cijs[[k]]; nij[[k, k]] = Plus @@ Map[# (# + 1)/2 &, cij]
       , {k, 1, nirep}];
       Do[
       	  cij = cijs[[k]];
          Do[
          	 cij1 = cijs[[j]];
             nij[[j, k]] = cij.cij1; nij[[k, j]] = nij[[j, k]]
          , {j, k + 1, nirep}]
       , {k, 1, nirep}];
       tb = Join[{Flatten[Join[{s}, names]]}, 
       Transpose[Join[{names}, nij // Transpose]]];
       tab = Grid[tb, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, Frame -> True,
       Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, None},{1,1}->GTCornerColor}];
       Print["Shell number ", s];
       Print[tab]
    , {s, 1, nbsg}]
  ]

*)

(*
***)



(****t* /GTTbRealSpaceMatrix
! NAME
!  GTTbRealSpaceMAtrix
! AUTHOR
!  W. Hergert
! PACKAGE
! TightBinding.m 
! MODIFICATION HISTORY
!  * 03.06.2014 : first version
!  * 20.06.2014 : revision
!  * 25.06.2018 : check header and documentation
! USAGE
!  GTTbRealSpaceMatrix[atom1,atom2,shell] constructs the interaction of atom1 and atom2 in a tight-binding Hamiltonian
!  in a certain neighbor shell.
! INPUT
!  * atom1      - description of atom 1, like {"Cu",{"s","p","p"}} (angular momentum "s*" is possible)
!  * atom2      - description of atom 2
!  * shell      - number of neighborhood shell 1,2,3
!  * dir        - distance vector between the two atoms               
! OUTPUT
!  TB interaction matrix 
! OPTIONS
!  GOTbBasis
!	o 0 no basis -> no superscripts at two-center parameters
!	o >0 basis -> superscripts constructed from the description of the atoms
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTTbHamiltonianElement
! GTPack NOTEBOOKS 
!  Test: GTTbRealSpaceMatrix.nb
! DESCRIPTION
!  GTTbHamiltonianElement has been modified to cover also the real space case.
!  If the angular momenta of the atoms atom1 und atom2 lead to the dimensions
!  nat1, nat2 an nat2xnat1 matrix is constructed.
! LITERATURE
!  J.C. Slater, G.F. Koster, Phys. Rev. 94, 1498 (1954)
! TODO
!  -
! PROBLEMS
!  In the first version the construction of the matrices was not stable. Instead of using Append during the
!  construction and a Partition at the end, the dimensions of the matrix are calculated first and the matrix elements
!  are calculated directly.
!
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTTbRealSpaceMatrix[atom1_, atom2_, sh_, OptionsPattern[]]:= 
      Module[{optb, mom, bas, hpart, mat, am1, am2, ang1, ang2, nd1, nd2, lz, ls, l1, l1s, l2, l2s, m1, m2, i, j}, 
      mom = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4};
      optb = OptionValue[GOTbBasis];
      If[optb == 0, 
      	 bas = 0, 
      	 bas = atom1[[1]] <> "," <> atom2[[1]]
      ];
      ang1 = atom1[[2]] /. mom; ang2 = atom2[[2]] /. mom;
      am1 = ang1 /. {-1 -> 0}; am2 = ang2 /. {-1 -> 0};
      nd1 = Plus @@ Map[(2 # + 1) &, am1 // Flatten];
      nd2 = Plus @@ Map[(2 # + 1) &, am2 // Flatten];
      hpart = Table[0, {nd1}, {nd2}];
      nd1 = Length[ang1]; nd2 = Length[ang2];
      lz = 0;
      Do[l1 = ang1[[i]];
         If[l1 == -1, 
         	l1s = 0, 
         	l1s = l1
          ];
          Do[
          	 lz = lz + 1; ls = 0; 
             Do[
             	l2 = ang2[[j]];
                If[l2 == -1, 
                   l2s = 0, 
                   l2s = l2
                ];
                Do[
                   ls = ls + 1;                        
                   mat = GTTbMatrixElement[l1, m1, l2, m2, sh, GOTbBasis -> bas];
                   hpart[[lz, ls]] = mat
                , {m2, -l2s, l2s}]
             , {j, 1, nd2}]
          , {m1, -l1s, l1s}]
       , {i, 1, nd1}];
       Return[hpart]
]

(*
***) 


(****t* /GTTbHamiltonianRS
! NAME
!  GTTbHamiltonianRS
! AUTHOR
!  W. Hergert
! PACKAGE
! TightBinding.m 
! MODIFICATION HISTORY
!  * 13.05.2014 : first version 
!  * 15.06.2014 : corrected version
!  * 20.07.2014 : test of the Hamiltonian
!  * 28.09.2014 : Problem with superscript if only one sort, but several atoms in unit cell solved
!  * 29.09.2014 : GTTbRealSpaceHamiltonian -> GTTbHamiltonianRS
!  * 25.06.2018 : check header and documentation
! USAGE
!  Constructs the real space Hamiltomian for a given crystal structure
! INPUT
!  * cl         - cluster of atoms corresponding to the strucutre
!  * admat      - adjacency matrix which belongs to the cluster cl
!  * basis      - basis information in the form {{atom A, number of A-atoms in cell, list of angular momenta},{ },{ },..}
!                                           
! OUTPUT
!  TB real space Hamiltonian matrix 
! GTPack OPTIONS
!  * GOTbBasis
!
!	o =0 no basis -> no superscripts at two-center parameters
!	o >0 basis -> superscripts constructed from the description of the atoms
!  * GOVerbose
!
!	o False -  minimal information (Standard)
!	o True  - information about intermediate steps.
! GTPack MODULES
!  GTRealSpaceMatrix, TbCompose, TbDCos 
! GTPack NOTEBOOKS 
!  Test: GTTbRealSpaceMatrix.nb
! DESCRIPTION
! The basis is formulated in the same way like for the k-space,  but information about the number of atoms in the basis is not really used.
! LITERATURE
!  -
! TODO
!  The Module seems to work also for semiconductures, but spurious peaks appear in the gap. What is the reason? It appears only for sp basis.
!  A transformation to hybrid orbitals does not help. -> the states are surface states of the finite cluster !
! RELEASE
!  1.0.0
! PROBLEMS
!  In the first version the interaction matrices have been calculated over and over again. This makes the whole story very slow.
!  Now the matrices are set up with  general direction cosines (Dl,DM,DN). The direction cosines are set to numerical values if necessary. 
!
! SOURCE
!--------------------------------------------------------------------------------
*)

  GTTbHamiltonianRS[cl_, admat_, basis_, OptionsPattern[]] := Module[
          {optb, verb, mom, cmd, msg1, adl, sh, ind, nat, bas, nb, sl, ndim, hp, i, ham, j, l, dist, dc, pos1, pos2, hpn, part, am, am1, k,
           bas1, basisam, basisat, mult, atomc, atomnc, test, tdim, sub, llist, iact, comp, nblock, ndblk, shn, ba, m, atom1, atom2, blk, 
           ww, pos,ll,l1,tab}, 
           mom = {"s*" -> -1, "s" -> 0, "p" -> 1, "d" -> 2, "f" -> 3, "g" -> 4};
  (*--- options ---*)
           optb = OptionValue[GOTbBasis];
           verb = OptionValue[GOVerbose];
  (*--- error messages ---*)
           cmd  = "GTTbHamiltonianRS : ";
           msg1 = "Cluster and adjacency matrix do not correspond.";
  (*--- cluster ---*)
           nat = Length[cl];
  (*--- adjacency matrix ---*)
           adl = admat[[1]];
           sh = Transpose[admat[[2]]][[1]]; ind = admat[[3]];
  (*--- correspondence of cluster and adjacency matrix ---*)         
           If[nat == adl,
              None,
              Print["Error :", cmd, msg1]; Abort[]
           ];
  (*---------------------------------------------------------------------------------*)         
  (*---                       one atom in basis                                   ---*)
  (*---------------------------------------------------------------------------------*)
           If[optb == 0,
              bas = {basis[[1, 1]], basis[[1, 3]]};
              nb = Plus @@ Map[(2 # + 1) &, bas[[2]] /. mom]; ndim = nb*nat;
              sl = Union[sh];
  (*--- print information ---*)
              tab={{"atoms in cluster","angular momenta","dimension Hamiltonian","shells"},
                   {" : "," : "," : "," : "},{nat bas[[1]],bas[[2]],ndim,sl}};
              Print[Grid[tab//Transpose,Frame->All,Background->{1->GTBackGroundColor1}]];
   (*           Print[nat, " ", bas[[1]], " atoms in cluster"];
              Print["angular momenta       : ", bas[[2]]];  
              Print["dimension Hamiltonian : ", ndim];
              Print["shells                : ", sl]; *)
              hp = {};
  (*--- construct parts of Hamionian ---*)
              Do[
                 hp = Append[hp,GTTbRealSpaceMatrix[bas, bas, sl[[i]], GOTbBasis -> 0]]
              , {i, 1, Length[sl]}];
  (*--- constructs full matrix ---*)
              ham = Table[0, {ndim}, {ndim}];
              Do[
  (*--- position in adjacency matrix ---*)
                 i    = IntegerPart[(ind[[l]] - 1)/adl] + 1;
                 j    = ind[[l]] - (i - 1)*adl;
                 dist = cl[[j, 1]] - cl[[i, 1]];
                 dc   = TbDcos[{dist}] // Flatten;
  (*--- position of block in Hamiltonian ---*)
                 pos1 = nb*(i - 1); 
                 pos2 = nb*(j - 1);
                 hpn  = Flatten[Position[sl, sh[[l]]]][[1]];
                 part = hp[[hpn]] //. {Global`DL -> dc[[1]], Global`DM -> dc[[2]], Global`DN -> dc[[3]]};
                 ham = TbCompose[ham, part, {pos1 + 1, pos2 + 1}, 0]
              , {l, 1, Length[sh]}];
              ham = ham + Transpose[ham];
   (*--- add diagonal blocks ---*)
              ind = 0;
              Do[
              	 am = bas[[2]] /. mom; am1 = am /. {-1 -> 0};
                 Do[
                 	l = am1[[j]];
                    Do[
                       ind = ind + 1;
                       ham[[ind, ind]] = GTTbSymbol2C[l, l, -1, 0, GOTbBasis -> 0]
                    , {k,1, 2*l + 1}]
                 , {j, 1, Length[am1]}]
              , {ll, 1, nat}]
              , None
           ];
  (*---------------------------------------------------------------------------------*)         
  (*---                 multiple atoms in basis                                   ---*)
  (*---------------------------------------------------------------------------------*)
           If[optb == 0,
              None,
  (*--- read informaion from basis ---*)             
              bas1 = basis // Sort;
              bas = Transpose[bas1]; 
              bas1 = Transpose[{bas[[1]], bas[[3]]}];
              basisat = Transpose[bas1][[1]]; 
              basisam = Transpose[bas1][[2]]; 
              mult = Tally[Transpose[cl][[2]]] // Sort;
   (*--- Print information ---*)      
              tab={{"atoms in basis","angular momenta","atoms in cluster"},
                   {" : "," : "," : "},{basisat,basisam,mult}};
              Print[Grid[tab//Transpose,Frame->All,Background->{1->GTBackGroundColor1}]];     
   (*           Print["Atoms in basis        : ", basisat];
              Print["Angular momenta       : ", basisam];
              Print["Atoms in cluster      : ", mult];*)
   (*---Test if basis and cluster correspond ---*)      
              atomc = Transpose[mult][[1]]; atomnc = Transpose[mult][[2]];
              test = Complement[atomc, basisat];
              If[test == {},
                 None,
                 Print["Error : ", cmd, msg1]; Abort[]
              ];
   (*--- Dimension of Hamiltonian and information for position of blocks ---*)          
              tdim = {basisat, Map[Plus @@ # &, Map[(2 # + 1) &, basisam /. mom /. {-1 -> 0}]]} // Transpose //Sort;
        ndim = atomnc.Transpose[tdim][[2]];
        Print["Dimension Hamiltonian : ", ndim];
   (*--- list for positions in full Hamiltonian ---*)
        sub = {};
        Do[
           sub = Append[sub, tdim[[j, 1]] -> tdim[[j, 2]]]
        , {j, 1, Length[tdim]}];
        llist = Transpose[cl][[2]] /. sub;
   (*--- construction of the nondiagonal blocks for later use ---*)
        iact = admat[[2]];
        comp = Union[iact]; nblock = Length[comp];
        If[verb,
           Print["Interactions (blocks) :", comp],
           None
        ];
        If[verb,      
           Print["Construction of ", nblock, " nondiagonal blocks"],
           None
        ];
        ndblk = {};
        Do[
           shn = comp[[i, 1]];bas = comp[[i, 2]];
           ba = StringSplit[bas, ","];
           m = Flatten[Position[basisat, ba[[1]]]][[1]]; atom1 = bas1[[m]];
           m = Flatten[Position[basisat, ba[[2]]]][[1]]; atom2 = bas1[[m]];
    (*--- no suberscripts if basis, but only one sort ---*)                    
         If[Length[basis]==1,
         	bas=0,
         	None
         ];
    (*--- ---*) 
           blk = GTTbRealSpaceMatrix[atom1, atom2, shn, GOTbBasis -> bas];
           If[verb,
              Print[atom1[[1]], "-", atom2[[1]], "-interaction"]; 
              Print[blk // MatrixForm],
              None
           ];
           ndblk = Append[ndblk, blk]
        , {i, 1, nblock}];
         Print["Construction of nondiagonal blocks finished!"];
   (*--- build the Hamiltonian ---*)
         ham = Table[0, {ndim}, {ndim}];
         Do[
   (*--- position in adjacency matrix ---*)
            i    = IntegerPart[(ind[[l]] - 1)/adl] + 1;
            j    = ind[[l]] - (i - 1)*adl;
            ww   = cl[[i, 2]] <> "," <> cl[[j, 2]]; 
            dist = cl[[j, 1]] - cl[[i, 1]];
    (*--- position in block list,get block ---*)
             pos = Flatten[Position[comp, ww]][[1]];
            blk = ndblk[[pos]];
            dc = TbDcos[{dist}] // Flatten;
    (*--- specialize block for interaction ---*) 
            part = blk //. {Global`DL -> dc[[1]], Global`DM -> dc[[2]], Global`DN -> dc[[3]]};
    (*--- position of block in Hamiltonian ---*)
            pos1 = Plus @@ Take[llist, i - 1];
            pos2 = Plus @@ Take[llist, j - 1];
            If[verb,
               Print["ind=", ind[[l]], " i=", i, " j=", j, " pos1=", pos1, " pos2=", pos2],
               None
            ];
    (*--- place the block ---*)       
            ham = TbCompose[ham, part, {pos1 + 1, pos2 + 1}, 0]
          , {l, 1, Length[ind]}];
          ham = ham + Transpose[ham];
          Print["All Nondiagonal parts are constructed!"];
    (*--- add diagonal elements ---*) 
          ind = 0;
          Do[
          	 m = Flatten[Position[basisat, cl[[ll, 2]]]][[1]];
          	 bas = basisat[[m]];
    (*--no suberscripts if basis, but only one sort---*)                    
             If[Length[basis]==1,
         	  bas=0,
              None
           ];
    (*--- ---*)
             am = basisam[[m]] /. mom ;
             Do[
             	If[am[[j]] == -1, 
             	   l1 = -1; l = 0, 
             	   l1 = am[[j]]; l = l1
             	];
                Do[
                   ind = ind + 1;
                   ham[[ind, ind]] = GTTbSymbol2C[l1, l1, -1, 0, GOTbBasis -> bas]
                , {k, 1, 2*l + 1}]
             , {j, 1, Length[am]}]
           ,{ll, 1, nat}]];
           Print["Full Hamitonian constructed!"];
           Return[ham // FullSimplify]
]


   
(*
***) 


(****t* /GTTbSymmetrySingleBand
! NAME
!  GTTbSymmetrySingleBand
! AUTHOR
!  W. Hergert
! PACKAGE
! TightBinding.m 
! MODIFICATION HISTORY
!  * 03/19/2016 : implementation
!  * 24.02.2017 : SOC included
!
! USAGE
!  GTTbSymmetryBand[coefficients,orbitals,chartab] performs a symmetry analysis for a certain band represented by the 
!  coefficients of the wave function. the orbitals are cartesian real spherical harmonics with respect to the basis set
!  of the ban structure calculation.
! INPUT
!  * coefficients - eigenvector of the corresponding band
!  * orbitals     - cartesian real spherical harmonics with respect to the basis set of the calculation
!  * chartab      - character table of the group of the wave vector
!             
! OUTPUT
! Irep
! OPTIONS
!  GOVerbose
!	o False -> no additional output (standard)
!	o True  -> output of additional information
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  Test: Symmetry_Analysis_Coefficients_TB.nb
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
! 
*)
   
GTTbSymmetrySingleBand[wave_, orbs_, chartab_, OptionsPattern[]] := Module[
    {verb, irep, wf2, pr, jj, so, nireps, lw},
  (*---options---*)
     verb   = OptionValue[GOVerbose];
     so     = OptionValue[GOSpinOrbit];
     nireps = Length[chartab[[2]]];
     Clear[wf2, xswf, yswf, zswf];
     irep   = {};
     If[so,
  (*---analysis with SO spinor WF ---*)
        lw = Length[wave];
        wf2[xswf_, yswf_, zswf_] := {orbs.wave[[1 ;; lw/2]], orbs.wave[[lw/2 + 1 ;; lw]]} /. 
        	{Global`x -> xswf, Global`y -> yswf, Global`z -> zswf},
   (*---analysis without SO---*)       
        wf2[xswf_, yswf_, zswf_] := orbs.wave /. 
        	{Global`x -> xswf, Global`y -> yswf, Global`z -> zswf}
     ];
     Do[
        pr = GTCharProjectionOperator[chartab[[1]], chartab[[2, jj]], wf2, {xswf, yswf, zswf}] // Chop;
        If[pr === 0 || pr === {0, 0},
           None,
           irep = Append[irep, chartab[[3, jj]]];
           If[verb,
              Print["IREP=",irep];
              Print[Grid[{{"Wave function",wf2[xswf, yswf, zswf] /. {xswf -> Global`x, yswf -> Global`y, zswf -> Global`z}},{"Projection",pr/. {xswf -> Global`x, yswf -> Global`y, zswf -> Global`z}}},Frame->All,Background -> {1 ->  GTBackGroundColor1}, 
 Dividers -> {2 -> GTDividerColor1}]],
              None
           ]
        ]
     , {jj, 1, nireps}];
     If[Length[irep] == 1,
        Return[irep[[1]]],
        Return[irep]
     ]
  ]
  
 
 
 (* auskommenitiert am 24.2 Variante aus symmetry.m uebernommen
GTTbSymmetrySingleBand[wave_, orbs_, chartab_, OptionsPattern[]] := Module[ 
	{verb,irep,wf2,pr,jj},	
  (*--- options ---*)
  verb = OptionValue[GOVerbose];
  (*--- analysis ---*)      
  Clear[wf2, xswf, yswf, zswf]; 
  wf2[xswf_, yswf_, zswf_] := orbs.wave /. {Global`x -> xswf, Global`y -> yswf, Global`z -> zswf};
  Do[
     pr = GTCharProjectionOperator[chartab[[1]], chartab[[2, jj]], wf2, {xswf, yswf, zswf}] // Chop;
     If[pr === 0, 
        None,
        irep = chartab[[3, jj]];
        If[verb,
           Print[irep, " ", wf2[xswf, yswf, zswf] /. {xswf -> Global`x, yswf -> Global`y, zswf -> Global`z}],
           None
        ]
     ]
  , {jj, 1, Length[chartab[[1]]]}];
  Return[irep]
 ]
*)

(*
***) 



(****t* /GTTbSymmetryPoint
! NAME
!  GTTbSymmetryPoint
! AUTHOR
!  W. Hergert
! PACKAGE
! TightBinding.m 
! MODIFICATION HISTORY
!  * 03/19/2016 : implementation
!  * 24.02.2017 : analysis with spin-orbit coupling included
!
! USAGE
!  GTTbSymmetryBand[kpoint,wave,orbs,minb,maxb,group,kbasis] performs a symmetry analysis at kpoint The bands are represented by the 
!  coefficients of the wave functions -wave. The orbitals are cartesian real spherical harmonics with respect to the basis set
!  of the ban structure calculation. The analysis will be performed for the bands [minb,maxb]. group is the point group of the crystal and
! INPUT
!  * kpoint       - k-point to analyse
!  * wave         - eigenvectors for this k-point
!  * orbs         - cartesian real spherical harmonics with respect to the basis set of the calculation
!  * [minb,maxb]  - bands to analyse
!  * group        - point group of the crystal
!  * kbasis       - basis vectors of the reciprocal lattice
!             
! OUTPUT
! list of Ireps
! OPTIONS
!  GOVerbose
!	o False -> no additional output (standard)
!	o True  -> output of additional information
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  Test: Symmetry_Analysis_Coefficients_TB.nb
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
! 
*)



GTTbSymmetryPoint[kpoint_, wave_, orbs_, minb_, maxb_, group_, kbasis_, 
	OptionsPattern[]] := Module[
	 {verb,so,reps,gofk,ctgofk,ireptab,irep,ib},
  (*---options---*)
    verb = OptionValue[GOVerbose];
    so   = OptionValue[GOSpinOrbit];
    reps = OptionValue[GOIrepNotation];
  (*---analysis---*)
    gofk   = GTGroupOfK[group, kpoint, kbasis];
    ctgofk = GTCharacterTable[gofk, GOVerbose -> verb, GOIrepNotation -> reps];
    If[so,
       ctgofk = GTExtraRepresentations[ctgofk, GOVerbose -> verb, GOIrepNotation -> reps],
       None
    ];
    ireptab={};
    Do[
       irep = GTTbSymmetrySingleBand[wave[[ib]], orbs, ctgofk, GOVerbose -> verb, GOSpinOrbit -> so];
     (*  Print[irep]; *)
       ireptab=Append[ireptab,irep]
    , {ib, minb, maxb}];
    Return[ireptab]
 ]

(* auskommentiert 24.2. variante aus symmetry.m uebernommen
GTTbSymmetryPoint[kpoint_, wave_, orbs_, minb_, maxb_, group_, kbasis_, OptionsPattern[]] := Module[
  {verb,repname,gofk,ctgofk,ireptab,irep,ib},
  (*--- options ---*)
  verb    = OptionValue[GOVerbose];
  repname = OptionValue[GOIrepNotation];
  (*--- analysis ---*)
  gofk    = GTGroupOfK[group, kpoint, kbasis];
  ctgofk  = GTCharacterTable[gofk, GOIrepNotation -> repname, GOVerbose -> verb];
  ireptab = {};
  Do[
     irep = GTTbSymmetrySingleBand[wave[[ib]], orbs, ctgofk, GOVerbose -> verb];
     ireptab = Append[ireptab, irep];
  , {ib, minb, maxb}];
  Return[ireptab]
]
*)


(*
***) 



(****t* /GTTbSymmetryBands
! NAME
!  GTTbSymmetryBands
! AUTHOR
!  W. Hergert
! PACKAGE
! TightBinding.m 
! MODIFICATION HISTORY
!  * 03/19/2016 : implementation
!  * 24.02.2017 : new version including SOC
!  * 10.03.2023 : correction of error:
!                 The output from GTBands into the bandstructure file was changed. The file contains now also the list of
!                 symmetry points. Thus, klist as argument is not necessary anymore. The code is changed accordingly.
! USAGE
!  GTTbSymmetryBands[bands,wave,minb,maxb,group,recip,basis] performs a symmetry analysis for the whole band structure.  
!  The analysis is based on the eigenvectors of the TB Hamiltonian and the orbitals represented by means of
!  cartesian real spherical harmonics with respect to the basis set of the band structure calculation. 
!
! INPUT
!  * bands        - band structure file
!  * wave         - file with eigenvectors 
!  * [minb,maxb]  - bands to analyse
!  * group        - point group of the crystal
!  * recip        - basis vectors of the reciprocal lattice
!  * basis
!             
! OUTPUT
!  list of Ireps and plot of band structure
! GTPack OPTIONS
!  GOIrepNotation - selects the required Irep notation
!  GOLabelStyle   - defines the style of the labels
!  GOPrecision    - 4 precision for numerical comparison
!  GOLabelShift   - global shift of labels
!  GOVerbose
!	o False -> no additional output (standard)
!	o True  -> output of additional information
!  GOPlot
!   o False -> no plot of band structure
!   o True  ->Plot of band structure
!  GOSpinOrbit
!   o True  -> with spin-orbit interaction
!   o False -> without spin-orbit interaction
!  GOOrbitalContruction
!   o True  -> internal construction of orbitals
!   o False -> provide list
! 
! Standard options
!   FrameLabel
!   PlotLabel
!   PlotStyle
!   Joined
!   PlotRange
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  Test: Symmetry_Analysis_Coefficients_TB.nb
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
! 
*)
GTTbSymmetryBands::band = "Band structure file `1` not in directory!";
GTTbSymmetryBands::wave = "Wave function  file `1` not in directory!";
GTTbSymmetryBands::orbs = "Number of orbitals and dimension of Hamiltonian are different.";

GTTbSymmetryBands[fileb_, filew_, minb_, maxb_, group_, kbasis_, basis_, OptionsPattern[]] := Module[
    {verb, repname, style, range, join, plot, del1, del2, nkp, files, bands, wav, orbs, norbs, nbands, i, j, irepl, ewl, nkpt, 
     kname, k, ew, kpoint, wave, ireps, enlist, ireplist, headline, tab,clist, plist, text, text1, arg, frame, plab, bplot, so, 
     nround,names, const,orbitals,lstyle,klist},
  (*---options---*)
     verb    = OptionValue[GOVerbose];
     frame   = OptionValue[FrameLabel];
     plab    = OptionValue[PlotLabel];
     repname = OptionValue[GOIrepNotation];
     style   = OptionValue[PlotStyle];
     range   = OptionValue[PlotRange];
     join    = OptionValue[Joined];
     plot    = OptionValue[GOPlot];
     nround  = OptionValue[GOPrecision];
     so      = OptionValue[GOSpinOrbit];
     const   = OptionValue[GOOrbitalConstruction];
     lstyle  = OptionValue[GOLabelStyle];
    {del1, del2} = OptionValue[GOLabelShift];
  (*---input bandstructure --*)
     If[Head[fileb] === String,
     	files = FileNames[];
        If[Intersection[{fileb}, files] == {},
           Message[GTTbSymmetryBands::band, fileb];   
           Return[],
           {bands, klist} = GTReadFromFile[fileb]
        ],
        {bands, klist} = fileb
     ];
(*---input wavefunction --*)
     If[Head[filew] == String,
       files = FileNames[];
       If[Intersection[{filew}, files] == {}, 
          Message[GTTbSymmetryBands::wave, filew];   
          Return[],
          wav = GTReadFromFile[filew]
       ],
       wav = filew
     ];
  (*--- number of k-points to be analysed ---*)     
     nkp = Length[klist];
  (*---construct the spherical harmonics for projection---*)
     orbs = {};
     orbitals={basis[[1,3]]};
     If[const,
        Do[
           Do[
              l    = orbitals[[i, j]];  
              orbs = Append[orbs, 
                     Flatten[Table[
                     GTCartesianTesseralHarmonicY[l, m, 
                     	Global`x, Global`y, Global`z]*(Global`x^2 + Global`y^2 + Global`z^2)^(l/2), {m, -l, l}]]]
           , {j, 1, Length[orbitals[[i]]]}]
        , {i, 1, Length[orbitals]}];
        orbs  = Flatten[orbs],
        orbs  = orbitals
     ];
     norbs = Length[orbs];
     nbands = Length[bands[[1,4]]];
     If[verb,
       Print["Basis functions: ", orbs];
       Print["Number of bands: ", nbands],
       None
     ];
  (*---consistency check---*)
  (*--- number of orbitals has to aggree with size of Hamiltonian ---*)
     If[norbs == nbands && so != True || 2*norbs == nbands && so,
        None,
        Message[GTTbSymmetryBands::orbs];Return[]
     ];
  (*---symmetry analysis---*)
     tab    = {};
     irepl  = {};
     ewl    = {};
     Do[
        nkpt   = klist[[k, 1]];
        kname  = klist[[k, 2]];
        kpoint = bands[[nkpt, 3]];
        Print[nkpt," ",kname," ", kpoint];
        ew     = Round[#*10.^nround]/10.^nround& /@bands[[nkpt, 4]];
        wave   = wav[[nkpt, 2]];
  (*---  use special names for extra representations ---*)
        If[so,
           If[Head[repname] === List,
              names = repname[[k]],
              names = {}
           ],
           names = repname
        ];
         ireps  = GTTbSymmetryPoint[kpoint, wave, orbs, minb, maxb, group, kbasis, 
                                     GOVerbose -> verb, GOIrepNotation -> names, GOSpinOrbit -> so
                                    ];
         irepl    = Append[irepl, ireps];
         ewl      = Append[ewl, ew];
         ireplist = Prepend[ireps, kname];
         enlist   = Prepend[ew, " "];
         tab      = Append[tab, {ireplist, enlist}]
     , {k, 1, nkp}];
     headline = Prepend[Table[i, {i, minb, maxb}], "Points"];
     tab      = Flatten[Prepend[tab, {headline}], 1];
     Print[Grid[tab, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, 
     	       {2 -> GTDividerColor1}},Background -> {{1 -> GTBackGroundColor1}, 
               {GTBackGroundColor1, {GTBackGroundColor2,White}},{1, 1} -> GTCornerColor}
              ]
     ];
     clist = Transpose[bands][[2]];
     plist = Map[clist[[#]] + del1 &, Transpose[klist][[1]]];
     text  = {};
     text1 = {};
     Do[
        Do[
           arg   = {irepl[[i, j]], {plist[[i]]+del1, ewl[[i, j]] + del2}};
           text1 = Append[text1, arg];
           If[lstyle=={},
              text  = Append[text, Text[irepl[[i, j]], {plist[[i]]+del1, ewl[[i, j]] + del2}]],
              text  = Append[text, Text[Style[irepl[[i, j]],lstyle], {plist[[i]]+del1, ewl[[i, j]] + del2}]]
           ];
        , {j, minb, maxb}]
     , {i, 1, nkp}];
     If[plot,          
        bplot = GTBandsPlot[{bands,klist}, maxb, PlotStyle -> style, Joined -> join, PlotRange -> range, FrameLabel -> frame, 
                         PlotLabel -> plab
                        ];
     Show[bplot, Graphics[text]],
     Return[{irepl, ewl, text1}]
     ]
  ]


(* auskommentiert 24.2.2017 neu aus symmetry.m uebernommen
GTTbSymmetryBands[fileb_, filew_, klist_, minb_, maxb_, group_, kbasis_, orbitals_,  OptionsPattern[]] := Module[
	{verb,repname,style,range,join,plot,del1,del2,nkp,files,bands,wav,orbs,norbs,nbands,i,j,irepl,ewl,nkpt,kname,k,
	 ew,kpoint,wave,ireps,enlist,ireplist,headline,tab,clist,plist,text,text1,arg,frame,plab,bplot},
  (*--- options ---*)
   verb    = OptionValue[GOVerbose];
   frame   = OptionValue[FrameLabel];
   plab    = OptionValue[PlotLabel];
   repname = OptionValue[GOIrepNotation];
   style   = OptionValue[PlotStyle];
   range   = OptionValue[PlotRange];
   join    = OptionValue[Joined];
   plot    = OptionValue[GOPlot];
   {del1, del2} = OptionValue[GOLabelShift];
    nkp    = Length[klist];
  (*--- input bandstructure and wavefunction ---*)  
   If[Head[fileb] == String,
     files = FileNames[];
     If[Intersection[{fileb}, files] == {},         
        Print["Error: band structure file ", fileb,  " not in the directory!"]; Abort[],
        bands = GTReadFromFile[fileb];
        bands = {bands, klist}
     ];
     If[Head[fileb] == List,
        bands = fileb,
        Print["Error: Input of band structure not correct!"]; Abort[]
     ]
   ];
   If[Head[filew] == String,
      files = FileNames[];
      If[Intersection[{filew}, files] == {},
         Print["Error: wavefunction file ", fileb, " not in the directory!"]; Abort[],
         wav = GTReadFromFile[filew]
      ];
      If[Head[fileb] == List,
         bands = fileb,
         Print["Error: Input of band structure not correct!"]; Abort[]
      ]
   ];
  (*--- construct the spherical harmonics for projection ---*) 
  orbs = {};
   Do[
      Do[
    	 l = orbitals[[i, j]];
         orbs = Append[orbs, Flatten[Table[GTCartesianTesseralHarmonicY[l, m, Global`x, Global`y, Global`z]*
         	                              (Global`x^2 + Global`y^2 + Global`z^2)^(l/2)
         	                         , {m, -l, l}]]]
      , {j, 1, Length[orbitals[[i]]]}]
   , {i, 1, Length[orbitals]}];
   orbs = Flatten[orbs]; norbs = Length[orbs];
   If[verb,
      Print["Basis functions: ", orbs],
      None
   ];
  (*--- consistency check ---*)
   nbands = Length[bands[[1, 1, 4]]];
   If[norbs == nbands,
      None,
      Print["Error: Number of orbitals and dimension of Hamiltonian are different."]; Abort[]
   ];
  (*--- symmetry analysis ---*)
   tab = {}; irepl = {}; ewl = {};
   Do[
      nkpt     = klist[[k, 1]]; 
      kname    = klist[[k, 2]];
      kpoint   = bands[[1, nkpt, 3]]; 
      ew       = bands[[1, nkpt, 4]];
      wave     = wav[[nkpt, 2]];
      ireps    = GTTbSymmetryPoint[kpoint, wave, orbs, minb, maxb, group, kbasis, GOVerbose -> verb,GOIrepNotation -> repname];
      irepl    = Append[irepl, ireps]; 
      ewl      = Append[ewl, ew];
      ireplist = Prepend[ireps, kname];
      enlist   = Prepend[ew, " "];
      tab      = Append[tab, {ireplist, enlist}]
   , {k, 1, nkp}];
   headline = Prepend[Table[i, {i, minb, maxb}], "Points"]; 
   tab = Flatten[Prepend[tab, {headline}], 1];
   Print[Grid[tab, Frame -> All,Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
              Background -> {{1 -> GTBackGroundColor1}, {GTBackGroundColor1, {GTBackGroundColor2, White}}, {1, 1} -> GTCornerColor}
             ]
   ];
   clist = Transpose[bands[[1]]][[2]];
   plist = Map[clist[[#]] + del1 &, Transpose[klist][[1]]];
   text = {};text1={};
   Do[
  	  Do[
  	  	 arg={irepl[[i, j]], {plist[[i]], ewl[[i, j]] + del2}};
  	  	 text1=Append[text1,arg];
         text = Append[text, Text[irepl[[i, j]], {plist[[i]], ewl[[i, j]] + del2}]];
      , {j, minb, maxb}]
   , {i, 1, nkp}];
   If[plot,
     bplot = GTBandsPlot[bands, maxb, PlotStyle -> style, Joined -> join, PlotRange -> range,
                        FrameLabel->frame,PlotLabel->plab];
     Show[bplot, Graphics[text]],
     Return[{irepl, ewl, text1}]
   ]
]
*)





(*
***) 





(****t* /GTTbTubeBands
! NAME
!  GTTbTubeBands
! AUTHOR
!  W. Hergert
! PACKAGE
! TightBinding.m 
! MODIFICATION HISTORY
!  * 21.12.2016 : implementation
!
! USAGE
!  GTTbTubeBands[hopin,n,m,kpath] calculates the band structure for a nanotube.
!
! INPUT
!  * hopin        - Hamiltonian or {Hamiltonian, Overlap
!  * n,m          - (n,m) nanotube
!  * kpath        - {{kmin,kmax,nkpt},{names of symmetry points}}
!             
! OUTPUT
! list of Ireps
! OPTIONS
!  GOVerbose
!	o False -> no additional output (standard)
!	o True  -> output of additional information
!  GOPlotBands
!   o True -> Data for band structur plot (GTBandsPlot) (Standard)
!   o False -> Data for DOS
!  GOTBOrthogonal
!   o True -> orthogonal Hamiltonian (Standard)
!   o False -> Nonorthogonal Hamiltonian
!  GTOStore
!   o 0 -> nothing stared (Standard)
!   o file -> stored to file
! ERROR MESSAGES
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  Test: Symmetry_Analysis_Coefficients_TB.nb
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTTbTubeBands[hopin_, n_, m_, kpath_, OptionsPattern[]] := Module[
	{nop,bnds,verb,store,hamop,overl,dim,b1,b2,a1,a2,dR,NR,t1,t2,tt,ttn,
	 K1,K2,ev,kmin,kmax,nkp,del,bands,ewall,kk,j,kp,mu,hop1,ew,over},
  (*---options---*)
    nop  = OptionValue[GOTbOrthogonal];
    bnds = OptionValue[GOPlotBands];
    verb = OptionValue[GOVerbose];
    store= OptionValue[GOStore];
    If[nop, 
       hamop = hopin, 
       hamop = hopin[[1]];
       overl = hopin[[2]]
    ];
    dim = Length[hamop];
  (*---reciprocal and direct lattice vectors---*)
    b1 = {1/Sqrt[3], -1/3, 0};
    b2 = {0, 2/3, 0};
    a1 = {Sqrt[3], 0, 0};
    a2 = {Sqrt[3]/2, 3/2, 0};
  (*---reciprocal vectors tube---*)
    dR   = GCD[2 n + m, 2 m + n];
    NR   = 2 (m^2 + n^2 + m n)/dR;
    t1   = (2 m + n)/dR;
    t2   = -(2 n + m)/dR;
    tt   = t1*a1 + t2*a2;
    ttn  = Norm[tt];
    K1   = (-t2 b1 + t1 b2)/NR;
    K2   = (m b1 - n b2)/NR;
    ev   = K2/Norm[K2];
    kmin = kpath[[1, 1]];
    kmax = kpath[[1, 2]];
    nkp = kpath[[1, 3]];
    del = (kmax - kmin)/(nkp - 1); 
    If[verb,
       Print[NR*dim," bands at ",nkp," k-points"],
       None
    ];    	
    bands = {};
    Do[
       ewall = {};
       Do[
       	  kk   = (kmin + (j - 1)*del);
          kp   = kk*ev/ttn + mu*K1;
          hop1 = hamop /. {\[Xi] -> kp[[1]], \[Eta] -> kp[[2]]};
          If[nop,
             ew   = Eigenvalues[hop1] // Chop // Sort,
             over = overl /. {\[Xi] -> kp[[1]], \[Eta] -> kp[[2]]};
             ew   = Eigenvalues[{hop1, over}] // Chop // Sort
          ];
          ewall= Append[ewall, ew]
       , {mu, 0, NR - 1}];
       If[bnds,
          bands = Append[bands, {j, Abs[kk - kmin], {kk, 0, 0}, Flatten[ewall]}],
          bands = Append[bands, {{kk, 0, 0}, Flatten[ewall]}]
       ]
    , {j, 1, nkp}];
  (*--- store data ---*)
    If[Head[store] == String, 
       Print["Write to file: ", store];
       If[bnds,
          GTWriteToFile[ {bands, {{1, kpath[[2, 1]]}, {nkp, kpath[[2, 2]]}}}, store],
          GTWriteToFile[bands, store]
       ],
       None
    ];  
    If[bnds,
       Return[{bands, {{1, kpath[[2, 1]]}, {nkp, kpath[[2, 2]]}}}],
       Return[bands]
    ]
 ]

(*
***)


(****t* /GTTbHamiltonianOnsite
! NAME
!  GTTbHamiltonianOnsite
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 21.12.2016 : first implementation
!  * 25.06.2018 : check header, create and link the documentation page
!
! USAGE
!  GTTbHamiltonianOnsite[ham,basis,struc] corrects the onsite energies with respect to the crystal field splitting.
!  
! INPUT
!  * ham         - tight-binding  Hamiltonian in analytic form
!  * basis       - basis used in the construction of the Hamiltonian
!  * structure   - structure can be fcc, bcc, hcp
! OUTPUT
!  Hamiltonian with splitting of on-site energies
! GTPack OPTIONS
!  * GOVerbose:       
!
!     - True  - additional information  
!     - False - no andditional information (standard)
!
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The module is not general enough. It is fixed to an spd-basis. An ERROR is printed, if not a complete basis is used.
! LITERATURE
! -
! TODO
!  Here is a need for a more flexible and general implementation.
! RELEASE
!  1.0.0
! PROBLEMS
!  ends with an error if basis is not spd
! SOURCE
!--------------------------------------------------------------------------------
!
*)


GTTbHamiltonianOnsite[ham_, basis_, struc_,OptionsPattern[]] := Module[
	{hop,verb,str,nbas,tbas,basat,nat,mom,tt,i,test,bas1},
	 hop  = ham;
  (*--- options ---*)
    verb = OptionValue[GOVerbose];
  (*--- check ---*)  
    str  = {"hcp", "fcc", "bcc"};
    If[Intersection[{{struc}, str}] == {},
       Print["Error: structure not implemented"]; Abort[],
       None
    ];
  (*--- analyse basis and parmeter set ---*)
    bas1  = basis /. {0 -> "s", 1 -> "p", 2 -> "d"};
    nbas  = Length[bas1];
    tbas  = Transpose[bas1];
    basat = tbas[[1]];
    nat   = Transpose[basis][[2]];
    mom   = tbas[[3]];
    If[verb,
       Do[
          If[nat[[i]] == 1,
             tt = " atom ",
             tt = " atoms "
          ];
          Print[nat[[i]], tt, "of type ", basat[[i]]]
       ,{i, 1, nbas}], 
       None
    ];
    If[Flatten[Union[mom]] == {"s", "p", "d"},
       None,
       Print["works only for complete basis sets in order s,p,d"]; Abort[]
    ];
   (*--- modify Hamiltonian  hcp---*)
    If[struc == "hcp" && nbas == 2,
       test = Map[StringTake[#, {1, 2}] &, basat] // Tally;
       If[Length[test] == 1,
          None,
          Print["Atoms are not identical ", test]; Abort[]
       ];
       hop[[6, 6]]   = ham[[6, 6]] /. {"(dd0)" -> "(dd1)"};
       hop[[8, 8]]   = ham[[8, 8]] /. {"(dd0)" -> "(dd1)"};
       hop[[7, 7]]   = ham[[7, 7]] /. {"(dd0)" -> "(dd2)"};
       hop[[15, 15]] = ham[[15, 15]] /. {"(dd0)" -> "(dd1)"};
       hop[[17, 17]] = ham[[17, 17]] /. {"(dd0)" -> "(dd1)"};
       hop[[16, 16]] = ham[[16, 16]] /. {"(dd0)" -> "(dd2)"};
       hop[[4, 5]]   = ham[[4, 5]] + "(pd0)";
       hop[[5, 4]]   = ham[[5, 4]] + "(pd0)";
       hop[[13, 14]] = ham[[13, 14]] + "(pd0)";
       hop[[14, 13]] = ham[[14, 13]] + "(pd0)";
       hop[[2, 9]]   = ham[[2, 9]] + "(pd0)";
       hop[[9, 2]]   = ham[[9, 2]] + "(pd0)";
       hop[[11, 18]] = ham[[11, 18]] + "(pd0)";
       hop[[18, 11]] = ham[[18, 11]] + "(pd0)";
       Return[hop],
       None
    ];
   (*--- modify Hamiltonian  fcc,bcc---*)
    If[struc == "bcc" && nbas == 1 || struc == "fcc" && nbas == 1,
       hop[[5, 5]] = ham[[5, 5]] /. {"(dd0)" -> "(dd1)"};
       hop[[6, 6]] = ham[[6, 6]] /. {"(dd0)" -> "(dd1)"};
       hop[[8, 8]] = ham[[8, 8]] /. {"(dd0)" -> "(dd1)"};
       hop[[7, 7]] = ham[[7, 7]] /. {"(dd0)" -> "(dd2)"};
       hop[[9, 9]] = ham[[9, 9]] /. {"(dd0)" -> "(dd2)"};
       Return[hop],
       Print["Conflict of data"]
    ]
 ]   
  
  
    
(*
***)


(****t* /GTTbReadWannier90
! NAME
!  GTTbReadWannier90
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 21.12.2016 : first implementation
!
! USAGE
!  GTFermiSurface[Hamiltonian, Fermi energy,  list of bands, ndel,  kbasis, clusterdata, kpath ] 
!  calculates the Fermi surface corresponding to a Hamiltonian and Fermi energy if the Fermi surface 
!  contains parts from list of bands. The electronic structure is calculated in a cube at ndel points 
!  per spatial dimension. kbasis is the basis of the reciprocal lattice. clusterdata contains the 
!  data for the lattice construction. The path used in electronic structure calculations can be given by BZpath.
! INPUT
!  ham         - parametrized Hamiltonian
!  fermi       - Fermi energy
!  nbands      - list bands where parts of the Fermi surface are expected
!  ndel        - number of points per spatial direction
!  kbasis      - basis vectors of the reciprocal lattice
!  cluster     - {cut,smin,smax} defines vectors for BZ construction
!  kpath       - Path in BZ
! OPTIONS
!  * GOVerbose        - additional information  (False)
!  
! OUTPUT
!  Plot of the Fermi surface with ListPlot3D
!  
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTVoronoiCell
! GTPack NOTEBOOKS 
!  Symmetry_Analysis-Coefficients_TB
! DESCRIPTION
!  
! LITERATURE
! 
! TODO
!  more tests
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)
 GTTbReadWannier90[file_, OptionsPattern[]] := Module[
 	{verb,fn,str,title,bands,vecs,notes,all,factors,tab,tab1,la,ttt},
  (*--- options ---*)
    verb = OptionValue[GOVerbose]; 
  (*---Check if database exists---*)
    fn = StringJoin[file, ".wann"];
    If[Flatten[Position[FileNames[], fn]] == {},
       Print["Error : database not in working directory"]; Abort[],
       None
    ];
    str     = OpenRead[fn];
    title = Read[str, String];
    bands = Read[str, Number];
    vecs  = Read[str, Number];
    notes = {{"Headline", "Number of bands", "Number of lattice vectors"}, {title, bands, vecs}};
    Print[Grid[notes // Transpose, Frame -> All, 
               Background -> {1 -> GTBackGroundColor1, None}
              ]
         ];
    all     = ReadList[str, Number]; la = Length[all];
    factors = Take[all, {1, vecs}];
    tab     = Partition[Take[all, {vecs + 1, la}], 7];
    If[verb,
       tab1 = {}; 
       Do[
       	  tab1 = Append[tab1, {i, tab[[i]]} // Flatten]
       , {i, 1, Length[tab]}]; 
       ttt = Join[{{"#", "a", "b", "c", "n", "m", "Re", "Im"}}, tab1];
       Print[Grid[ttt, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
                       Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, 
                       {1, 1} -> GTCornerColor}
                 ]
             ],
       None
    ];
    Return[{bands, vecs, factors, tab}]
  ]
  
 (*
 ***)
 
 
  
(****t* /GTTbWannier90Hamiltonian
! NAME
!  GTTbWannier90Hamiltonian
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 21.12.2016 : first implementation
!
! USAGE
!  GTFermiSurface[Hamiltonian, Fermi energy,  list of bands, ndel,  kbasis, clusterdata, kpath ] 
!  calculates the Fermi surface corresponding to a Hamiltonian and Fermi energy if the Fermi surface 
!  contains parts from list of bands. The electronic structure is calculated in a cube at ndel points 
!  per spatial dimension. kbasis is the basis of the reciprocal lattice. clusterdata contains the 
!  data for the lattice construction. The path used in electronic structure calculations can be given by BZpath.
! INPUT
!  ham         - parametrized Hamiltonian
!  fermi       - Fermi energy
!  nbands      - list bands where parts of the Fermi surface are expected
!  ndel        - number of points per spatial direction
!  kbasis      - basis vectors of the reciprocal lattice
!  cluster     - {cut,smin,smax} defines vectors for BZ construction
!  kpath       - Path in BZ
! OPTIONS
!  * GOVerbose        - additional information  (False)
!  
! OUTPUT
!  Plot of the Fermi surface with ListPlot3D
!  
! ERROR MESSAGES
!  -
! GTPack MODULES
!  GTVoronoiCell
! GTPack NOTEBOOKS 
!  Symmetry_Analysis-Coefficients_TB
! DESCRIPTION
!  
! LITERATURE
! 
! TODO
!  more tests
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTTbWannier90Hamiltonian[data_ , bas_, OptionsPattern[]] :=  Module[
	{eps,cut,verb,bands,vecs,fc,ham,dist,ind,set,i,tab,rv,dd,distl,dmax,h0,a,
	 m,n,fac,ik,k},
	 eps = 10^(-5);
  (*--- options ---*)
  	cut   = OptionValue[GOCutOff];
  	verb  = OptionValue[GOVerbose];
  (*--- extract information from data ---*)
    bands = data[[1]];
    vecs  = data[[2]];
    fc    = data[[3]];
    tab   = data[[4]];
    ham   = Table[0, {bands}, {bands}];
    dist  = {};
    Do[ 
       ind = bands^2 (i - 1);
       set = tab[[ind + 1]];
       rv  = set[[1]]*bas[[1]] + set[[2]]*bas[[2]] + set[[3]]*bas[[3]];
       dd  = Norm[rv] // FullSimplify; dist = Append[dist, dd];
    , {i, 1, vecs}];
    distl = Tally[dist];
    distl = Sort[Transpose[distl][[1]], #1 < #2 &];    
    If[cut == 0,
       dmax = distl[[Length[distl]]],
       dmax = distl[[cut + 1]]
    ];
    If[verb,
    	Print["dmax      : ",dmax];
        Print["Distances : ", distl];	
       None	
    ];
    Do[
       h0  = Table[0, {bands}, {bands}];  
       ind = bands^2 (i - 1);   
       set = tab[[ind + 1]];
       rv  = set[[1]]*bas[[1]] + set[[2]]*bas[[2]] + set[[3]]*bas[[3]];        
       fac = Exp[2 \[Pi]  I*rv.{\[Xi], \[Eta], \[Zeta]}];
       If[Norm[rv] <= dmax,
          Do[
             ik  = ind + k;
             set = tab[[ik]];        
             m   = set[[4]]; 
             n   = set[[5]];
             If[Abs[set[[7]]] < eps, 
             	a = 0, 
              	a = set[[7]]
             ];
            h0[[m, n]] = set[[6]] + I*a // Chop   
          , {k, 1, bands^2}];
          ham = ham + fac*h0/fc[[i]],
          None
       ]
   , {i, 1, vecs}]; 
   Return[ham]
]


    
(*
***)

(*-------------------------- Attributes ------------------------------*)

Attributes[GTTbSpinMatrix]={Protected, ReadProtected}
Attributes[GTTbMatrixElement]={Protected, ReadProtected}
Attributes[GTTbSimplify]={Protected, ReadProtected}
Attributes[GTTbSymbol2C]={Protected, ReadProtected}
Attributes[GTTbHamiltonianElement]={Protected, ReadProtected}
Attributes[GTTbHamiltonian]={Protected, ReadProtected}
Attributes[GTTbSpinOrbit]={Protected, ReadProtected}
Attributes[GTSymmetryBasisFunctions]={Protected, ReadProtected}
Attributes[GTTbAngularPart]={Protected, ReadProtected}
Attributes[GTIrepInfo]={Protected, ReadProtected}
Attributes[GTTbNumberOfIntegrals]={Protected, ReadProtected}
Attributes[GTTbRealSpaceMatrix]={Protected, ReadProtected}
Attributes[GTTbHamiltonianRS]={Protected, ReadProtected}
Attributes[GTTbParameterSet]={Protected, ReadProtected}
Attributes[GTTbParameterNames]={Protected, ReadProtected}
Attributes[GTTbGetParameter]={Protected, ReadProtected}
Attributes[GTTbParmToRule]={Protected, ReadProtected}
Attributes[GTTbDatabaseInfo]={Protected, ReadProtected}
Attributes[GTTbDatabaseUpdate]={Protected, ReadProtected}
Attributes[GTTbDatabaseRetrieve]={Protected, ReadProtected}
Attributes[GTTbPrintParmSet]={Protected, ReadProtected}
Attributes[GTTbAtomicWaveFunction]={Protected, ReadProtected}
Attributes[GTTbBlochFunction]={Protected, ReadProtected}
Attributes[GTTbWaveFunction]={Protected, ReadProtected}
Attributes[GTTbParmToFortran]={Protected, ReadProtected}
Attributes[GTTbToFortran]={Protected, ReadProtected}
Attributes[GTFindStateNumbers]={Protected, ReadProtected}
Attributes[GTPlotStateWeights]={Protected, ReadProtected}
Attributes[GTHamiltonianPlot]={Protected, ReadProtected}


End[] (* End Private Context *)




(**/---------------------- End GTPack2.0 --------------------------------------*)


EndPackage[]
