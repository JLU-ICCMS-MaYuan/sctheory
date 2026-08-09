(* ::Package:: *)

(****m* /Auxiliary.m
!
! NAME
!  Auxiliary.m
! AUTHOR
!  W. Hergert, M. Geilhufe
! MODIFICATION HISTORY
!  * 24.02.2012	- initial documentation 
!  * 14.05.2013	- modified usage message
!  * 28.12.2017	- check headers of entire package 
!  * 19.05.2018	- check doku, headers and package header
! USAGE
!  
! GTPack MODULES
!
! --- Quaternion operations ---
!
!  * GTQAbs            - gives the absolute value of a quaternion
!  * GTQConjugate      - gives the conjugate quaternion of a quaternion
!  * GTQInverse        - gives the inverse of a quaternion
!  * GTQMultiplication - multiplies two quaternions
!  * GTQPolar          - gives the polar angle of a quaternion
!
! --- Logical operators ---
!
!  * GTEulerAnglesQ    - gives True if the input is a list of Euler angles, and gives False otherwise 
!  * GTQuaternionQ     - gives True if the input is a quaternion, and gives False otherwise 
!  * GTSymbolQ gives   - True if if the input is a symbol, and gives False otherwise
!
! --- Spherical Harmonics ---
!
!  * GTCartesianSphericalHarmonicY  - gives the spherical harmonics in cartesian form 
!  * GTCartesianTesseralHarmonicY   - gives the tesseral harmonics in cartesian form 
!  * GTTesseralHarmonicY            - gives the tesseral harmonic (real spherical harmonic) 
!  * GTGauntCoefficient             - calculates Integral over three spherical harmonics
!
! --- Graphic representation of point group relationships ---
! 
!  * GTPointGroups     - gives a graph according to the subgroup relationships of the 32 point groups
!  * GTGroupHierarchy  - gives a graph containing all subgroups of a group. Also all groups which contain the given group as a subgroup are included in the graph     
!  * GTGroupConnection - gives the connection of the groups in a list of groups as a graph                 
!  * GTCaleyGraph      - plots the CaleyGraph of a point group
!
! --- Input and output operations ---
! 
!  * GTReadFromFile     - reads an object from a file
!  * GTWriteToFile      - writes an object to a file 
!  * GTFermiSurfaceXSF  - exports data to XCrxysDen to plot Fermi surfaces
!  * GTExtractDatasets  - extracts various sample data sets used in the documentation
!  * GTSetMessage       - Set messages to Notebook or Message window
!
! --- ab initio data ---
!
!  * GTVASPBandsPlot  - reads and plots a VASP band structure
!  * GTVASPBands      - reads EIGENVAL and transforms for band structure plot
!
! --- Handling of cluster data ---
! 
!  * GTNeighborPlot   - is a tool to present the information about the input data to construct adjacency matrices in visual form
!  * GTClusterFilter  - removes certain types of atoms from a given cluster
!  * GTCompactStore   - allows to store sparse matrices in a compact form
!
! ---Special matrices ---
! 
!  * GTSU2Matrix    - gives the SU(2) matrix for a counterclockwise rotation around a 3D vector
!  * GTDiracMatrix - gives the Dirac matrices
! 
! --- Presentation of calculated data ---
!
!  * GTSetTableColors   - is used to change the color scheme of tables in GTPack
!  * GTBlueRed          - Color function from blue to red
!  * GTCharacterTableP  - gives pretty print of a character table
!
! --- Internal modules but not private 
!
!  * GTNotationToSFL      - is used to switch notation to SFL in string format 
!  * GTOverScript         - Construction of Hermann-Mauguin symbols (Overline)
!  * GTSubScript          - Construction of Hermann-Mauguin symbols (Index)
!  * GTSF                 - Construction of Schoenfliess space group notations
!  * GTConnectionPattern  - is used to construct the graphs
!  * GTDoubleNames        - searches for the double point group names switch to GTPack standard
!  * GTNotationText       - constructs Labels for the graphs
!  * GTSimplify           - simplifys expressions
!  * GTSF                 - construction of Schoenfliess space group notations
!
! --- iternal functions ---
!
!  * GTnlmcoeff            - coefficients for spherical tesseral harmonics  
!  * GTNotationToSFL       - used to switch notation to SFL in string format
!  * GTOverScript          - Construction of Hermann-Mauguin symbols (Overline)
!  * GTSubScript           - Construction of Hermann-Mauguin symbols (Index).
!  * GTSF                  - Construction op Schoenfliess space group notations
!  * GTSimplify            - simplifies a given expression.
!  * GTPointInSquareQ      - gives True if point lies inside the square defined by vertices
!  * GTConnectionPattern   - is used to construct the graphs
!  * GTDoubleNames         - Search for the double point group names, switch to GTPack standard
!  * GTNotationText        - Construct Labels for the graphs 
!
! DESCRIPTION
!  Auxiliary.m contains all additional functions that are needed by the GroupTheory package. 
!  These comprise of e.g. handling of quaternions, spherical harmonics or group relations. 
!
! LITERATURE
!  Wolfram Hergert, Matthias Geilhufe
!  Group Theory in Solid State Physics and Photonics - Problem solving with Mathematica
!
! TODO
!  
! PROBLEMS
!
!--------------------------------------------------------------------------------
!
***)

BeginPackage["GroupTheory`Auxiliary`",{"GroupTheory`Symbols`","GroupTheory`Basic`","GroupTheory`Install`"}]

(*-------------------------- Quaternion operations -----------------*)
 GTQAbs                         ::usage = "GTQAbs[\*StyleBox[\"q\", \"TI\"]] gives the absolute value of the quaternion \*StyleBox[\"q\", \"TI\"]."
 GTQConjugate                   ::usage = "GTQConjugate[\*StyleBox[\"q\", \"TI\"]] gives the conjugate quaternion \*SuperscriptBox[StyleBox[\"q\", \"TI\"], \"*\"] of the quaternion \*StyleBox[\"q\", \"TI\"]."
 GTQInverse                     ::usage = "GTQInverse[\*StyleBox[\"q\", \"TI\"]] gives the inverse of the quaternion \*StyleBox[\"q\", \"TI\"]."
 GTQMultiplication			  ::usage = "GTQMultiplication[\*StyleBox[\"q1, q2\", \"TI\"]] multiplies the quaternions \*StyleBox[\"q1\", \"TI\"] and \*StyleBox[\"q2\", \"TI\"]."
 GTQPolar                       ::usage = "GTQPolar[\*StyleBox[\"q\", \"TI\"]] gives the polar angle of the quaternion \*StyleBox[\"q\", \"TI\"]."
 
(*-------------------------- Logical operators ---------------------*)
 GTEulerAnglesQ                 ::usage = "GTEulerAnglesQ[\*StyleBox[\"A\", \"TI\"]] gives \*ButtonBox[\"True\", BaseStyle->\"Link\", ButtonData->\"paclet:ref/True\"] if \*StyleBox[\"A\", \"TI\"] is a list of Euler angles, and gives \*ButtonBox[\"False\", BaseStyle->\"Link\", ButtonData->\"paclet:ref/False\"] otherwise."
 GTQuaternionQ                  ::usage = "GTQuaternionQ[\*StyleBox[\"q\", \"TI\"]] gives \*ButtonBox[\"True\", BaseStyle->\"Link\", ButtonData->\"paclet:ref/True\"] if \*StyleBox[\"q\", \"TI\"] is a Quaternion, and gives \*ButtonBox[\"False\", BaseStyle->\"Link\", ButtonData->\"paclet:ref/False\"] otherwise."
 GTSymbolQ                      ::usage = "GTSymbolQ[\*StyleBox[\"A\", \"TI\"]] gives \*ButtonBox[\"True\", BaseStyle->\"Link\", ButtonData->\"paclet:ref/True\"] if \*StyleBox[\"A\", \"TI\"] is a symbol, and gives \*ButtonBox[\"False\", BaseStyle->\"Link\", ButtonData->\"paclet:ref/False\"] otherwise."

(*-------------------------- Spherical Harmonics -------------------*)
 GTCartesianSphericalHarmonicY	::usage = "GTCartesianSphericalHarmonicY[\*StyleBox[\"l,m,x,y,z\", \"TI\"]] gives the spherical harmonic \*Cell[BoxData[FormBox[RowBox[{SubsuperscriptBox[\"Y\", \"l\", \"m\"], \"(\", RowBox[{\"x\",\",\",\"y\",\",\",\"z\"}],\")\"}],TraditionalForm]],\"InlineMath\"] in cartesian form."
 GTCartesianTesseralHarmonicY	::usage = "GTCartesianTesseralHarmonicY[\*StyleBox[\"l,m,x,y,z\", \"TI\"]] gives the tesseral harmonic \*Cell[BoxData[FormBox[RowBox[{SubsuperscriptBox[\"S\", \"l\", \"m\"], \"(\", RowBox[{\"x\",\",\",\"y\",\",\",\"z\"}],\")\"}],TraditionalForm]],\"InlineMath\"] in cartesian form."
 GTGauntCoefficient             ::usage = "GTGauntCoefficient[\*StyleBox[\"l1,m1,l2,m2,l3,m3\", \"TI\"]] calculates the integral over three spherical harmonics \*Cell[BoxData[FormBox[RowBox[{\"\[Integral]\",RowBox[{SubsuperscriptBox[\"Y\", \"l1\", \"m1\"],SubsuperscriptBox[\"Y\", \"l2\", \"m2\"],SubsuperscriptBox[\"Y\", \"l3\", \"m3\"],RowBox[{\"\[DifferentialD]\", \"\[CapitalOmega]\"}]}]}],TraditionalForm]], \"InlineMath\"]."
 GTTesseralHarmonicY	        ::usage = "GTTesseralHarmonicY[\*StyleBox[\"l,m,\[Theta],\[Phi]\", \"TI\"]] gives the tesseral harmonic (real spherical harmonic)\*Cell[BoxData[FormBox[RowBox[{SubsuperscriptBox[\"S\", \"l\", \"m\"], \"(\", RowBox[{\"\[Theta]\",\",\",\"\[Phi]\"}],\")\"}],TraditionalForm]],\"InlineMath\"]."
 
(*----- graphic representation of point group relationships --------*)
 GTGroupConnection              ::usage = "GTGroupConnection[\*StyleBox[\"list of groups\", \"TI\"]] gives the connection of the groups in \*StyleBox[\"list of groups\", \"TI\"] as a graph."
 GTGroupHierarchy               ::usage = "GTGroupHierarchy[\*StyleBox[\"group\", \"TI\"]] gives a graph containing all subgroups of a \*StyleBox[\"group\", \"TI\"]. Also all groups which contain \*StyleBox[\"group\", \"TI\"] as a subgroup are included in the graph."      
 GTPointGroups                  ::usage = "GTPointGroups[\*StyleBox[\"Mode\", \"TI\"]] gives a graph according to the subgroup relationships of the 32 point groups."  
 GTCayleyGraph					::usage = "GTCayleyGraph[\*StyleBox[\"group\", \"TI\"]] plots the Cayley graph of a point group."
                  
 
(*------------------------- Input and Output Operations -------------*)
 GTExtractDatasets		  	    ::usage = "GTExtractDatasets[] extracts various sample data used in the documentation.\nGTExtractDatasets[\*StyleBox[\"destination\", \"TI\"]] extracts the sample data to \*StyleBox[\"destination\", \"TI\"] directory."
 GTFermiSurfaceXSF              ::usage = "GTFermiSurfaceXSF[\*StyleBox[\"file,Hamiltonian,kbasis,ndel,list of bands,Fermi energy,information\", \"TI\"]] calculates the energy band structure for a given \*StyleBox[\"list of bands\", \"TI\"] from \*StyleBox[\"Hamiltonian\", \"TI\"] in the primitive unit cell in k-space, defined by \*StyleBox[\"kbasis\", \"TI\"] on a regular mesh, defined by \*StyleBox[\"ndel\", \"TI\"]. The results are stored in \*StyleBox[\"file\", \"TI\"].bxsf as data readable in XCrySDen. \*StyleBox[\"information\", \"TI\"] is added to the data."
 GTReadFromFile                 ::usage = "GTReadFromFile[\*StyleBox[\"name\", \"TI\"]] reads an object from the file \*StyleBox[\"name\", \"TI\"]."
 GTWriteToFile                  ::usage = "GTWriteToFile[\*StyleBox[\"object,file\", \"TI\"]] writes \*StyleBox[\"object\", \"TI\"] to the file \*StyleBox[\"file\", \"TI\"]." 
 GTSetMessage                   ::usage = "GTSetMessage[\*StyleBox[\"set\", \"TI\"]] sets output of messages to \*StyleBox[\"set\", \"TI\"] (Notebook,Console)."
(*--------------------------Ab initio band structures ---------------*)
 GTVASPBandsPlot                ::usage = "GTVASPBandsPlot[\*StyleBox[\"file,number kpoints, number bands, kpoints\", \"TI\"]] reads and plots  a VASP bandstructure."
 GTVASPBands                    ::usage = "GTVASPBandsPlot[\*StyleBox[\"file,kpoints\", \"TI\"]] reads EIGENVAL and transforms data for band structure plot." 
(*-------------------------- Handling of Cluster data -------------------------*)
 GTClusterFilter                ::usage = "GTClusterFilter[\*StyleBox[\"cluster,types\", \"TI\"]] removes certain \*StyleBox[\"types\", \"TI\"] of atoms from a given \*StyleBox[\"cluster\", \"TI\"]." 
 GTCompactStore                 ::usage = "GTCompactStore[\*StyleBox[\"matrix\", \"TI\"]] allows to store sparse matrices in a compact form."
 GTNeighborPlot                 ::usage = "GTNeighborPlot[\*StyleBox[\"neighbors\", \"TI\"]] is a tool to present the information about the input data to construct adjacency matrices in visual form."

(*-------------------------- Special matrices --------------*)
 GTDiracMatrix                  ::usage = "GTDiracMatrix[\*StyleBox[\"i\", \"TI\"]] gives the \*Cell[BoxData[FormBox[SuperscriptBox[\"i\", \"th\"], TraditionalForm]], \"InlineMath\"] Dirac matrix."
 GTSU2Matrix                    ::usage = "GTSU2Matrix[\*StyleBox[\"\[Theta],w\", \"TI\"]] gives the SU(2) matrix for a counterclockwise rotation around the 3D vector \*StyleBox[\"w\", \"TI\"]." 
 
(*-------------------------- Presentation of calculated data    ----------------------*)
 GTBlueRed                      ::usage = "GTBlueRed[\*StyleBox[\"x\", \"TI\"]] defines a ColorFunction from blue (\*StyleBox[\"x\", \"TI\"]=0) to red (\*StyleBox[\"x\", \"TI\"]=1)."
 GTSetTableColors               ::usage = "GTSetTableColors[\*StyleBox[\"colors\", \"TI\"]] defines the color scheme for tables in GTPack." 
 GTCharacterTableP              ::usage = "GTCharacterTableP[\*StyleBox[\"charcter table\", \"TI\"]] prints a character table if elements only im matrix form given."
 
 GTReorderMatrix                ::usage = "GTReorderMatrix[dim_,list_] gives a matrix to interchange rows and columms in a matrix."
 GTGetIrepMatrix                ::usage = "GTGetIrepMatrix[\*StyleBox[\"group,matrices,selection\", \"TI\"]] prints the representation \*StyleBox[\"matrices\", \"TI\"] to the \*StyleBox[\"selection\", \"TI\"] of symmetry elements of an irreducible representation of \*StyleBox[\"group\", \"TI\"]."
 (*-------------------------- internal functions not private ------------*)
 
 GTNotationToSFL                ::usage = "GTNotationToSFL[symbol] is used to switch notation to SFL in string format." 
 GTOverScript                   ::usage = "Construction of Hermann-Mauguin symbols (Overline)."
 GTSubScript                    ::usage = "Construction of Hermann-Mauguin symbols (Index)."
 GTSF					        ::usage = "Construction op Schoenfliess space group notations."
 GTSimplify						::usage = "GTSimplify[expression] simplifies a given expression."
 GTPointInSquareQ               ::usage = "GTPointInSquareQ[vertices,point] gives True if point lies inside the square defined by vertices."
 GTConnectionPattern            ::usage = "GTConnectionPattern is used to construct the graphs."
 GTDoubleNames                  ::usage = "Search for the double point group names, switch to GTPack standard."
 GTNotationText                 ::usage = "Construct Labels for the graphs" 
 GTGetUnitaryTrafoMatrix        ::usage = "GTGetUnitaryTrafoMatrix determines a unitary matrix U with the property A = U*B*U^-1 and U^n = C for three matrices A, B, and C."

 
 

(*--------------------------- Options -----------------------------*) 
Options[GTCompactStore]      = {GOCompact -> False, GOMatrixType -> "General"}
Options[GTConnectionPattern] = {GONotation->"SFL",GODimension->3}
Options[GTDiracMatrix]       = {GORepresentation -> "Dirac"};
Options[GTDoubleNames]       = {GOVerbose -> False}
Options[GTGroupConnection]   = {GOSelectCase->"Restricted",GOPlotStyle->"Normal",GONotation->"SFL",GraphStyle -> "DiagramGold"}
Options[GTExtractDatasets]   = {GOOverWrite -> False}
Options[GTGroupHierarchy]    = {GOPlotStyle->"Normal",GONotation->"SFL",GraphStyle -> "DiagramGold"}
Options[GTNotationText]      = {GONotation -> "SFL", GODimension -> 3, GOVerbose -> False}
Options[GTPointGroups]       = {GODimension->3,GONotation->"SFL",GraphStyle -> "IndexLabeled"}
Options[GTSetTableColors]    = {GOVerbose -> True}
Options[GTVASPBandsPlot]     = {GOVerbose -> False, Joined -> True, FrameLabel -> {" ", "Energy (ev)"}, PlotRange -> All, PlotLabel -> "Band structure", PlotStyle -> Blue, GOShift -> 0, GOFermiEnergy -> 0}
Options[GTVASPBands]         = {GOVerbose -> True, GOStore -> 0}
Options[GTGetIrepMatrix]     = {GOVerbose -> True}
Options[GTCayleyGraph]      = {GOVerbose -> False, GOColorScheme -> {Blue, Red, Green},GOPlot -> True, VertexSize -> 0.8, VertexLabelStyle -> Directive[Red, 15], VertexShapeFunction -> "Circle"}

(*-------------------------- Modules ------------------------------*) 
 
 Begin["`Private`"] 
 (****h* /GTSetMessage
! NAME
!  GTSetMessage
! AUTHOR
!  W. Hergert
! PACKAGE
! Auxiliary.m
! MODIFICATION HISTORY
!  * 101.04.2023- 1st version 
! USAGE
!  GTSetMessage[name] sets the output of messages to Notebook or Console
! INPUT
!  * name    - output direction (Notebook, Console)
! OUTPUT
!   -
! GTPack OPTIONS
!   -
! DESCRIPTION
!  When we started to define our error messages in terms of the Mathematica message systems, the messages came to
!  the Notebook. Now the messages are going to the Message window. It seems to be better to have the methods directly 
!  in the Notebook. The command canbe used to switch between both options. 
! LITERATURE
!  Method found in the internet in MAthematica forum
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTSetMessage::console =" Messages are going to Message window.";
GTSetMessage::notebook=" Messages are going to Notebook.";

GTSetMessage[set_]:=Module[{current},
    If[Intersection[{set},{"Console","Notebook"}]=={},
       Print["Only Console or Notebook allowed for setting."];
       Return[],
       None
    ];
(*--- Notebook ---*)
    current=CurrentValue[$FrontEndSession,{MessageOptions,"KernelMessageAction"}];
    Print["current setting : ",current];
    If[current==="PrintToNotebook",
       If[set=="Notebook",
          Message[GTSetMessage::notebook],      
          current=CurrentValue[$FrontEndSession,{MessageOptions,"KernelMessageAction"}]="PrintToConsole";
          Message[GTSetMessage::console]
       ],
       None
    ];
(*--- Console --*);
    If[current==="PrintToConsole",
       If[set==="Console",
          Message[GTSetMessage::console] ,
          current=CurrentValue[$FrontEndSession,{MessageOptions,"KernelMessageAction"}]="PrintToNotebook";
          Message[GTSetMessage::notebook]
       ],
       None

    ]
]

(*
***) 

(****h* /GTCayleyGraph
! NAME
!  GTCayleyGraph
! AUTHOR
!  W. Hergert
! PACKAGE
! CrystalStructure.m
! MODIFICATION HISTORY
!  * 10.12.2021 - 1st version 
! USAGE
!  GTCayleyGraph[name] plots the Cayley graph of the corresponding point group
! INPUT
!  * name    - name of the point group
! OUTPUT
!   * cayley graph if GOPlot->True
!   * definition of graph for use with graph commands
! GTPack OPTIONS
!  * GOVerbose:
!
!     - False - no additional information
!     - True  - detailed information 
!  * GOPlot: 
!
!     - False - no plot, output of edges
!     - True  - plot of the Cayley graph
!  * GOColorScheme:
!
!    The colors for the edges that represent tehgenerators
!
! STANDARD OPTIONS
!  * VertexSize
!  * VertexLabelStyle
!  * VertexShapeFunction
! GTPack MODULES
!  GTInstallGroup
!  GTGenerators
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The existing implementations of the Cayley graph do not add the corresponding names to vertices and edges. Here we realize this.
!  All information for the graph is found by GTPack moduls. FFinally Mathematica is used to plot the graph.
! LITERATURE
!-
! TODO
!  A complete test of all point groups shold be done.
! RELEASE
!  1.0.0
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
! 
*)



GTCayleyGraph::name       = "`1` is not an allowed point group name.";
GTCayleyGraph::generators = "Number of indentified generators differs from true number.";
GTCayleyGraph::col        = "Not eneough colors for the generators.";



GTCayleyGraph[gruppe_, OptionsPattern[]] := Module[
	{verb, cols, s, ls, gruppe1, grp, nelm, gens, ngen, gdat, tal,edge, tab1, m, tab, vs, plt,names,ltab,graph}, 
     names = {"C1", "Ci", "C2", "Cs", "C2h", "D2", "C2v", "D2h", "C4", 
              "S4", "C4h", "D4", "C4v", "D2d", "D4h", "C3", "S6", "D3", "C3v", 
              "D3d", "C6", "C3h", "C6h", "D6", "C6v", "D3h", "D6h", "T", "Th", 
              "O", "Td", "Oh"};
  (*--- options ---*)
     verb = OptionValue[GOVerbose];
     cols = OptionValue[GOColorScheme];
     s    = OptionValue[VertexSize];
     ls   = OptionValue[VertexLabelStyle];
     vs   = OptionValue[VertexShapeFunction];
     plt  = OptionValue[GOPlot];
  (*--- check input ---*)
     gruppe1 = ToString[gruppe];
     If[Intersection[names, {gruppe1}] == {},
        Message[GTCayleyGraph::name, gruppe]; Return[],
        None
     ];
   (*  Print[$ContextPath];*)
  (*--- group and generators ---*)
     grp  = GTInstallGroup[gruppe1];
     Print[grp];
     nelm = Length[grp];
     gens = GTGenerators[grp];
     ngen = Length[gens];
     If[ngen > Length[cols],
        Message[GTCayleyGraph::col]; Return[],
        None
     ];
     Print[Grid[{{" ", gruppe}, {"elements", grp}, {"generators", gens}}, Frame -> All, 
           Background -> {1 -> GTBackGroudColor1, 1 -> GTCornerColor, {1, 1} -> White}]];
     tab = {};
     Do[ 
        tab = Append[tab, {gens[[i]], cols[[i]]}]
     , {i, 1, ngen}];
     Print[Grid[tab, Frame -> All,Background -> {1 -> GTBackGroudColor1, None}]];
  (*--- Table connection with generators ---*)
     ltab     = nelm + 1;
     tab      = Table[0, {ltab}, {ltab}];
     tab[[1]] = {0, grp} // Flatten;
     tab      = Transpose[tab];
     tab[[1]] = {gruppe, grp} // Flatten;
     Do[
        Do[
           Do[
              If[grp[[i]]\[SmallCircle]gens[[k]] == grp[[j]],
                 tab[[i + 1, j + 1]] = k,
                 None
              ]
          , {j, 1, nelm}]
       , {i, 1, nelm}]
     , {k, 1, ngen}];
     If[verb,
        Print[Grid[tab, Frame -> All, Background -> {1 -> GTBackGroudColor1,1 -> GTBackGroudColor1, {1, 1} -> GTCornerColor}]],
        None
     ];
  (*--- Generation of graph ---*)
     gdat = tab[[2 ;; ltab, 2 ;; ltab]] // Flatten; 
     tal  = Drop[Transpose[Tally[gdat]][[2]], 1];
     If[Length[tal] == ngen,
        None,
        Message[GTCaleyGraph1::generators]; Return[]
     ];
     edge = {};
     Do[
        tab1 = Table[0, {tal[[i]]}];
        edge = Append[edge, tab1], {i, 1, Length[tal]}];
        Do[m = 1;
           Do[
              Do[
                 If[tab[[i + 1, j + 1]] == k,
                    edge[[k, m]] = DirectedEdge[i, j]; m = m + 1,
                    None
                 ];
              , {i, 1, nelm}]
           , {j, 1, nelm}]
        , {k, 1, ngen}];
     Do[
        edge[[k]] = Style[#, cols[[k]]] & /@ edge[[k]]
     , {k, 1, ngen}];
     graph = edge // Flatten;
  (*--- Result ---*)
     panelLabel[lbl_] := 
     Panel[lbl, FrameMargins -> None, Background -> None];
     If[plt,
        Show[Graph[graph, VertexLabels -> Table[i -> Placed[grp[[i]] // StandardForm, 
        	       Center, panelLabel], {i, nelm}], VertexShapeFunction -> vs, VertexSize -> s, 
                   VertexLabelStyle -> ls]
            ],
            Return[graph]
     ]
  ]


(*
***) 
(****b* /GTCharacterTableP
! NAME
!  GTCharacterTableP
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  * 11.01.2017 : first version
!  * 28.12.2017 : check header and documentation
!  * 18.05.2018 : element from class in header instead of classe name, new docu  
! USAGE
!  GTCharacterTableP[character  table] gives a pretty print of a charaxcter table
! INPUT
!  character table
! OUTPUT
!  print of character table
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  
! GTPack MODULES
!  
! DESCRIPTION
!  Some times it is useful to get a pretty print of a character table which was 
!  already calculated by means of GTCharacterTable 
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
GTCharacterTableP[chartab_] := Module[
  {cls, ctab, names, ncl, elmcl, clname},
    cls    = chartab[[1]];
    ctab   = chartab[[2]];
    names  = chartab[[3]];
    ncl    = Length[cls]; 
    elmcl  = Length[#] & /@ cls; 
    clname = Table[ToString[elmcl[[i]]] cls[[i,1]], {i, 1, ncl}]; 
    clname = Prepend[clname, " "]; 
    ctab   = Prepend[Prepend[ctab // Transpose, names] // Transpose,  clname]; 
    Print[Grid[ctab, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
               Background -> {{1 -> GTBackGroundColor1}, {1 -> 
               GTBackGroundColor1}, {1, 1} -> GTCornerColor}
              ]
         ]
  ]
  
(*
***)

(****b* /GTReorderMatrix
! NAME
!  GTReorderMatrix
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  11.01.2017 : first version
!  28.12.2017 : check header and documentation
! USAGE
!  GTReorderMatrix[dim_,list_] creates matrix to interchange rows and columms in a matrix.
! INPUT
!   dim  - dimension
!   list - pairs of numbers definig the reordering
! OUTPUT
!  extracted data
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  
! GTPack MODULES
!  
! DESCRIPTION
!  Reordering of a matrix can be performed by means of matrix multiplication. The modul creates a 
!  a matrix T to reorder rows and columns simultaneously: A' = TAT
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

GTReorderMatrix[ndim_, list_] := 
 Module[{cmat, nl, test, rm, i1, i2, i, k},
  cmat = IdentityMatrix[ndim];
  nl = Length[list];
  If[list == {},
        Return[cmat],
   test = Length[#] & /@ list // Union;
   If[Length[test] == 1 && test[[1]] == 2,
         None,
         Print["Error: List has to be a list of pairs!"]; Abort
    ];
   (*--- remove diagonals --*)
   rm = Flatten[list];
   Do[i1 = list[[i, 1]]; i2 = list[[i, 2]];
         cmat[[i1, i2]] = cmat[[i2, i1]] = 1
    , {i, 1, nl}];
   Do[k = rm[[i]]; cmat[[k, k]] = 0
    , {i, 1, 2*nl}];
   Return[cmat]
   ]
  ]


(*
***)

(****b* /GTExtractDatasets
! NAME
!  GTExtractDatasets
! AUTHOR
!  S. Schenk
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  19/09/2016 : new version
!  18.05.2018 : check docu, GOOverWirte not documented
! USAGE
!  GTExtractDatasets[dest_] extracts the zip name to dir.
! INPUT
!  destination of extraction
! OUTPUT
!  extracted data
! GTPack OPTIONS
!  GOOVerWrite - regulates overwrite of old data 
! Standard OPTIONS
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
GTExtractDatasets[des_:False, OptionsPattern[]] := Module[{src, dest, ow}, 
	If[StringQ[des],
		dest = des, 
		If[FileExistsQ[$BaseDirectory <> "/Applications/GroupTheory/_datasets.zip"], 
			dest = $BaseDirectory <> "/Applications/GroupTheory/Documentation/English/ReferencePages/Symbols/datasets"];
		If[FileExistsQ[$UserBaseDirectory <> "/Applications/GroupTheory/_datasets.zip"], 
			dest = $UserBaseDirectory <> "/Applications/GroupTheory/Documentation/English/ReferencePages/Symbols/datasets"]
	];
	If[FileExistsQ[$BaseDirectory <> "/Applications/GroupTheory/_datasets.zip"], 
		src = $BaseDirectory <> "/Applications/GroupTheory/_datasets.zip"];
	If[FileExistsQ[$UserBaseDirectory <> "/Applications/GroupTheory/_datasets.zip"], 
		src = $UserBaseDirectory <> "/Applications/GroupTheory/_datasets.zip"];
	If[Head[des] == Rule, If[GOOverWrite /. des, ow = True]];
	If[DirectoryQ[dest] && TrueQ[OptionValue[GOOverWrite] || ow], 
		DeleteDirectory[dest, DeleteContents -> True]];
	If[! DirectoryQ[dest], CreateDirectory[dest];
		ExtractArchive[src, dest]];
]
(*
***)

(****b* /GTCartesianSphericalHarmonicY
! NAME
!  GTCartesianSphericalHarmonicY
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  12/17/2013 : new version
!  01/09/2014 : changes
!  18.05.2015 : check docu
! USAGE
!  GTCartesianSphericalHarmonicY[l_,m_,x_,y_,z_] gives the spherical harmonic Y_l^m(x,y,z) in cartesian form.
! INPUT
!  quantum numbers l and m
!  cartesian coordinates x,y,z
! OUTPUT
!  function of x, y and z
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  GTnlmcoeff
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
GTCartesianSphericalHarmonicY[l_, m_, x_, y_, z_] := Simplify[PowerExpand[GTnlmcoeff[l, m]*LegendreP[l, m,z/Sqrt[x^2 + y^2 + z^2]]*Sum[Binomial[Abs[m],k]*(x/Sqrt[x^2+y^2])^k*(y/Sqrt[x^2+y^2])^(Abs[m]-k)*(Cos[\[Pi]*(Abs[m]-k)/2]+I*Sign[m]*Sin[\[Pi]*(Abs[m]-k)/2]),{k,0,Abs[m]}]]]

(*
***)

(****b* /GTCartesianTesseralHarmonicY
! NAME
!  GTCartesianTesseralHarmonicY
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  12/17/2013 : new version
!  01/09/2014 : new version according to Vogel, Podolski (PhysRevB 2004)
!  12.05.2018 : check docu
! USAGE
!  GTCartesianTesseralHarmonicY[l_,m_,x_,y_,z_] gives the tesseral harmonic S_l^m (x,y,z) in cartesian form.
! INPUT
!  quantum numbers l and m
!  cartesian coordinates x,y,z
! OUTPUT
!  function of x, y and z
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  GTnlmCoeff
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
GTCartesianTesseralHarmonicY[l_, m_, x_, y_, z_] := Module[{k},
Simplify[
	PowerExpand[Which[	m > 0,	Sqrt[2]*((-1)^m)*GTnlmcoeff[l, m]*LegendreP[l, m, z/Sqrt[x^2 + y^2 + z^2]]*Sum[Binomial[m,k]*(x/Sqrt[x^2+y^2])^k*(y/Sqrt[x^2+y^2])^(m-k)*Cos[\[Pi]*(m-k)/2],{k,0,m}], 
  		m < 0,	Sqrt[2]*((-1)^m)*GTnlmcoeff[l, Abs[m]]*LegendreP[l, Abs[m], z/Sqrt[x^2 + y^2 + z^2]]*Sum[Binomial[Abs[m],k]*(x/Sqrt[x^2+y^2])^k*(y/Sqrt[x^2+y^2])^(Abs[m]-k)*Sin[\[Pi]*(Abs[m]-k)/2],{k,0,Abs[m]}],
  					m == 0, GTCartesianSphericalHarmonicY[l, 0, x, y, z]]]]]

(*
***)

(****b* /GTDiracMatrix
! NAME
!  GTDiracMatrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  28.12.2017 : Check header and documentation
! USAGE
!  GTDiracMatrix[i] gives the ith Dirac matrix.
! INPUT
!  index
! OUTPUT
!  Dirac matrix
! GTPack OPTIONS
!  GORepresentation : sets representation for Dirac Matrices
! Standard OPTIONS
!  
! DESCRIPTION
!  Dirac matrices can be generated in different representations
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

GTDiracMatrix[io_, OptionsPattern[]] := Module[
  {indlist, prefac, i},
  i = io;
  If[io == 0, i = 4];
  Which[OptionValue[GORepresentation] == "Dirac",
   indlist = {{2, 1}, {2, 2}, {2, 3}, {3, 0}, {1, 0}};
   prefac = {I, I, I, 1, 1},
   OptionValue[GORepresentation] == "Majorana",
   indlist = {{1, 3}, {3, 0}, {1, 1}, {1, 2}, {2, 1}};
   prefac = {I, I, -I, -1, -1},
   OptionValue[GORepresentation] == "Weyl",
   indlist = {{2, 1}, {2, 2}, {2, 3}, {1, 0}, {3, 0}};
   prefac = {I, I, I, 1, -1},
   OptionValue[GORepresentation] == "Liu",
   indlist = {{1, 1}, {2, 1}, {3, 1}, {0, 2}, {0, 3}};
   prefac = {1, 1, 1, 1, 1}];
  prefac[[i]]*
   KroneckerProduct[PauliMatrix[indlist[[i, 1]]], 
    PauliMatrix[indlist[[i, 2]]]]
  ]
(*
***)

(****b* /GTnlmcoeff
! NAME
!  GTTnlmcoeff
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  
! USAGE
!  GTTnlmcoeff[l_,m_] gives coefficients for the construction of the Tesseral Spherical Harmonics
! INPUT
!  quantum numbers l and m 
! OUTPUT
!  coefficient
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
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
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTnlmcoeff[l_, m_] := Sqrt[(2*l + 1)*(l - m)!/(4*\[Pi]*(l + m)!)]


(*
***)


(****b* /GTGetUnitaryTrafoMatrix
! NAME
!  GTGetUnitaryTrafoMatrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  RepresentationTheorySG.m
! MODIFICATION HISTORY
!   12/05/2018: first version
! USAGE
!   GTGetUnitaryTrafoMatrix determines a unitary matrix U with the property A = U*B*U^-1 and U^n = C for three matrices A, B, and C.
! INPUT
!   Three matrices A, B, C, matrix power n
! OUTPUT
!   a matrix U with the property A = U*B*U^-1 and U^n = C
! GTPack OPTIONS
!
! STANDARD OPTIONS
!
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  needed by GTSGGetIreps, GTGetIreps for the induction of two representations from an irreducible representationwith orbit 1
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
GTGetUnitaryTrafoMatrix[Amat_, Bmat_, Cmat_, n_] := 
 Module[{mdim, uu, U, inf1, inf2, eq1, eq2, eq4, sol},
  mdim = Length[Amat];
 (*Define a dummy matrix U*)
  uu = Table[U[i, j], {i, 1, mdim}, {j, 1, mdim}];
 (*Impose conditions U.A-B.U=0*)
  inf1 = uu . Amat - Bmat . uu // Simplify;
  eq1 = Table[inf1[[i, j]] == 0, {i, 1, mdim}, {j, 1, mdim}];
  (*Impose conditions U.U=C*)
  inf2 = MatrixPower[uu, n] // Simplify;
  eq2 = Table[inf2[[i, j]] == Cmat[[i, j]], {i, 1, mdim}, {j, 1, mdim}];
  (*inf3=uu.ConjugateTranspose[uu]-ConjugateTranspose[uu].uu//Simplify;
    eq3=Table[inf3[[i,j]]\[Equal]0,{i,1,mdim},{j,1,mdim}];*)
    eq4 = Det[uu] == 1;
  (*****************************************)
  (*First attempt to find a solution*)
  (*****************************************)
  (*Solve the equation system*)
  sol = If[mdim > 1, FindInstance[Flatten[{eq1, eq2, eq4}], Flatten[uu]], FindInstance[Flatten[{eq1, eq2}], Flatten[uu]]];
  Return[uu /. Flatten[sol]]]


GTSGGetUnitaryTrafoMatrixV2[Amat_, Bmat_,Cmat_,n_] := Module[{p,sol,eqs,ll,eva, evb,Um,Umsq},
  Print["Failed"];
  eva = FullSimplify[Orthogonalize[Eigenvectors[Amat]]];
  evb = FullSimplify[Orthogonalize[Eigenvectors[Bmat]]];
  ll=Length[eva];
  Um=Simplify[Sum[Exp[I p[i]] KroneckerProduct[eva[[i]],Conjugate[evb[[i]]]],{i,1,ll}]];
  Umsq=Simplify[MatrixPower[Um,n]];
  eqs=Flatten[{Table[Umsq[[i,j]]==Cmat[[i,j]],{i,1,ll},{j,1,ll}],Table[0<=p[i]<2 \[Pi],{i,1,ll}]}];
  sol=Flatten[FullSimplify[Normal[Quiet[Solve[eqs,Table[p[i],{i,1,ll}]]]]/.C[m_]->1]];
  Um=Simplify[Um/.sol/.C[m_]->1];
  If[Length[sol]==0,Print["Error: GTSGGetUnitaryTrafoMatrix: Couldn't determine unique U to induce representation matrices."];Print[Amat,Bmat,Cmat,n];Abort[],
  Return[Um]]
  ]
(*
***)


(****b* /GTTesseralHarmonicY
! NAME
!  GTTesseralHarmonicY
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  01/09/2014 : new version according to Vogel, Podolski (PhysRevB 2004)
!  1805.2018  : check docu
! USAGE
!  GTTesseralHarmonicY[l_,m_,t_,p_] gives the tesseral harmonic 
!  (real spherical harmonic) S^m_l (\theta,\phi)
! INPUT
!  quantum numbers l and m
!  Angles t, p
! OUTPUT
!  function of t,p
! GTPack OPTIONS
!
! Standard OPTIONS
!  
! GTPack MODULES
!  GTnlmcoeff
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

GTTesseralHarmonicY[l_, m_, t_, p_] := 
 Which[m > 0, 
  Sqrt[2]*((-1)^m)*GTnlmcoeff[l, m]*LegendreP[l, m, Cos[t]]*Cos[m*p], m == 0, 
  SphericalHarmonicY[l, m, t, p], m < 0, 
  Sqrt[2]*((-1)^m)*GTnlmcoeff[l, Abs[m]]*LegendreP[l, Abs[m], Cos[t]]*Sin[Abs[m]*p]]

(*
***)

(****b* /GTEulerAnglesQ
! NAME
!  GTEulerAnglesQ
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  5/14/13    : modified usage message
!  18.05.2018 : documentation checked
! USAGE
!  GTEulerAnglesQ[A] gives True if A is a list of Euler angles, and gives
!  False otherwise.
! INPUT
!  Object A
! OUTPUT
!  logical
! GTPack OOPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  -
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

GTEulerAnglesQ[a_] := If[(Length[a]==2)&&(Length[a[[1]]]==3)&&((a[[2]]==1)||(a[[2]]==-1)), True, False]
   
(*
***) 

(****b* /GTOverScript
! NAME
!  GTOverScript
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m
! MODIFICATION HISTORY
!  05/03/2013 : first version
! USAGE
!  internal function which is used for the construction of Hermann-Mauguin symbols
! INPUT
!  a number (type number)
! OUTPUT
!  the number with overline (type symbol)
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack Modules
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

GTOverScript[x_]:="\*OverscriptBox[\""<>ToString[x]<>"\",\"_\"]"    


(*
***)


(****b* /GTSubScript
! NAME
!  GTSubScript
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m
! MODIFICATION HISTORY
!  05/03/2013 : first version
! USAGE
!  internal function which is used for the construction of Hermann-Mauguin symbols
! INPUT
!  number and index
! OUTPUT
!  number with index in correct output form
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
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
GTSubScript[x_, y_] :=  "\*SubscriptBox[\"" <> ToString[x] <> "\",\"" <> ToString[y] <> "\"]"


(*
***)


(****b* /GTSF
! NAME
!  GTSF
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m
! MODIFICATION HISTORY
!  06/07/2013 : first version
! USAGE
!  internal function which is used for the construction of Schoenfliess space group symbol 
! INPUT
!  the space group information 
! OUTPUT
!  the space group symbol (type symbol)
! GTPack OPTIONS
!  - 
! Standard OPTIONS
!  -
! GTPack MODULES
!  -
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


GTSF[s_, low_, up_] := 
 "\*SubsuperscriptBox[\"" <> ToString[s] <> "\",\"" <> ToString[low] <>
   "\",\"" <> ToString[up] <> "\"]"
   
   
(*
***)  
   
(****b* /GTConnectionPattern
! NAME
!  GTConnectionPattern
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m
! MODIFICATION HISTORY
!  05/03/2013 : first version
!  January 2016 : fitted to the new scheme of pointgroup names
! USAGE
!  it is used to set up the graph of groups
! INPUT
!  list of symbols
! OUTPUT
!  symbol for the path in the graph 
! GTPack OPTIONS
!   GONotion    : "HM" or "SFL" gives  the notation to which the symbol has to be transformed
!	GODimension : 2D or 3D point groups
! STANDARD OPTIONS
! 
! GTPack MODULES
!  GTNotationText
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

GTConnectionPattern[x_, OptionsPattern[]] := 
  Module[{nt, t, x1, x2, dim}, 
  nt = OptionValue[GONotation];
  dim = OptionValue[GODimension];
  x1 = GTNotationText[x[[1]], GONotation -> nt, GODimension -> dim];
  If[NumberQ[x[[2]]], 
  	x2 = x[[2]], 
    x2 = GTNotationText[x[[2]], GONotation -> nt, GODimension -> dim]
  ];
  t = x1 -> x2;
  Return[t]
]


(*
***)


(****b* /GTNotationToSFL
! NAME
!  GTNotationToSFL
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  29/05/2015 : introduced in package
!  28.12.2017 : check header
! USAGE
!  The standard notation for point groups in GTPack is the SFL notation. 
!  GTNotationToSFL translates all possible input:
!  - HM notation in BacketingBars
!  - symbol SFL notation
!  - string SFL notation
!  into a string in SFL notation.          
!  internal function          
! INPUT
!  point group name
! OUTPUT
!  point grouup name in SFL notation
! GTPack Options
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTInstallGroup
!  GTGroupConnection
!  GTGroupHierarchy
!  
! DESCRIPTION
!  The BracketingBar construct allows the construction of HM symbols. All previously 
!  defind commands do not work with this construct. This command allows to keep the 
!  other commands unchanged.
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


GTNotationToSFL[symb0_, OptionsPattern[]] := 
 Module[{names, symbn, i,symb1},
  names = {{"Oh", "O", "Th", "Td", "D6h", "D4h", "T", "C6h", "C6v", 
     "D6", "D3h", "D3d", "C4h", "C4v", "D4", "D2h", "D2d", "C6", 
     "C3i", "C3h", "C3v", "D3", "C4", "S4", "C2h", "C2v", "D2", "Cs", 
     "Ci", "C2", "C1", 
     "C3"}, {\[LeftBracketingBar]m-3m\[RightBracketingBar], \
\[LeftBracketingBar]432\[RightBracketingBar], \[LeftBracketingBar]m-3\
\[RightBracketingBar], \[LeftBracketingBar]-43 \
m\[RightBracketingBar], \[LeftBracketingBar]6/
       mmm\[RightBracketingBar], \[LeftBracketingBar]4/
       mmm\[RightBracketingBar], \[LeftBracketingBar]23\
\[RightBracketingBar], \[LeftBracketingBar]6/
       m\[RightBracketingBar], \[LeftBracketingBar]6 mm\
\[RightBracketingBar], \[LeftBracketingBar]622\[RightBracketingBar], \
\[LeftBracketingBar]-6 m2\[RightBracketingBar], \
\[LeftBracketingBar]-3 m\[RightBracketingBar], \[LeftBracketingBar]4/
       m\[RightBracketingBar], \[LeftBracketingBar]4 mm\
\[RightBracketingBar], \[LeftBracketingBar]422\[RightBracketingBar], \
\[LeftBracketingBar]mmm\[RightBracketingBar], \[LeftBracketingBar]-42 \
m\[RightBracketingBar], \[LeftBracketingBar]6\[RightBracketingBar], \
\[LeftBracketingBar]-3\[RightBracketingBar], \[LeftBracketingBar]-6\
\[RightBracketingBar], \[LeftBracketingBar]3 m\[RightBracketingBar], \
\[LeftBracketingBar]32\[RightBracketingBar], \[LeftBracketingBar]4\
\[RightBracketingBar], \[LeftBracketingBar]-4\[RightBracketingBar], \
\[LeftBracketingBar]2/
       m\[RightBracketingBar], \[LeftBracketingBar]mm2\
\[RightBracketingBar], \[LeftBracketingBar]222\[RightBracketingBar], \
\[LeftBracketingBar]m\[RightBracketingBar], \[LeftBracketingBar]-1\
\[RightBracketingBar], \[LeftBracketingBar]2\[RightBracketingBar], \
\[LeftBracketingBar]1\[RightBracketingBar], \[LeftBracketingBar]3\
\[RightBracketingBar]}};
  If[Head[symb0] === BracketingBar,
       Do[
            If[symb0 == names[[2, i]],
                 symbn = names[[1, i]]
              , None
          ]
    , {i, 1, 32}]; Return[symbn]
   , None];
   If[Head[symb0] === BracketingBar,
       None,
              If[ StringQ[symb0],
                        symb1 = symb0,
                       symb1 = ToString[symb0]
                 ]; Return[symb1]]
  ]

(*
***)

(****b* /GTGauntCoefficient
! NAME
!  GTGauntCoefficient
! AUTHOR
!  M.Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  07/03/2016 : introduced in package
!  28.12.2017 : check header
! USAGE
!  GTGauntCoefficient[l1_, m1_, l2_, m2_, l3_, m3_] calculates the integral over 
!  three spherical harmonics.                      
! INPUT
!  l1,m1,l2,m2,l3,m3
! OUTPUT
!  Gaunt coefficient
! GTPack Options
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! DESCRIPTION
!  The Gaunt coefficients are represented by ThreeJSymbol implemented in Mathematica
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

GTGauntCoefficient[l1_, m1_, l2_, m2_, l3_, m3_] := Module[{},Quiet[Sqrt[(2 l1 + 1)*(2 l2 + 1)*(2 l3 + 1)/(4 \[Pi])] ThreeJSymbol[{l1, 0}, {l2, 0}, {l3,0}] ThreeJSymbol[{l1, m1}, {l2, m2}, {l3, m3}]]]

(*
***)

(****b* /GTGroupConnection
! NAME
!  GTGroupConnection
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  01/03/2013   : introduced in package
!  05/03/2013   : options introduced
!  02/01/2014   : VertexLabels in "Embedding"-Mode corrected
!  January 2015 : correct labels in the graph
!  28.12.2017   : check header and description
!  18.05.2018   : check of docu, loading package in options removed
!  06.04.2020   : some prpblem appeared with the labeling of the graph, perhaps caused by the new Mathematica version
!                 problem solved by using the stanard option GraphStyle
!  07.04.2020   : documentation page changed
! USAGE
!  GTGroupConnection[grps] gives the connection of the groups in list of groups as a graph.                      
! INPUT
!  grps  - list of point groups
! OUTPUT
!  plot of the graph
! GTPack Options
!  GOPlotStyle  : "Normal"     -> no embedding
!               : "Embedding"  -> Embedding in the whole connection graph
!  GOSelectCase : "All"        -> also connections to further subgroups
!               : "Restricted" -> Only connections between elements in list grps
!  GONotation   : "SFL"        -> Schoenfliess
!               : "HM"         -> Herman-Mauguin
! STandard OPTIONS
!  GraphStyle   : "DiagramGold"
! GTPack MODULES
!  GTConnectionPattern, GTNotationText
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
GTGroupConnection[grp0_, OptionsPattern[]] := 
 Module[{nt,names3D,names,grp1,lab,nn,con3D,pgg,sg,l,xx,i,str1,str2,l1,l2,j,k,m}, 
  nt = OptionValue[GONotation];
  (* rfunc=OptionValue[VertexRenderingFunction];*)
  names3D = {{"Oh", 1}, {"O", 2}, {"Th", 3}, {"Td", 4}, {"D6h", 
     5}, {"D4h", 6}, {"T", 7}, {"C6h", 8}, {"C6v", 9}, {"D6", 
     10}, {"D3h", 11}, {"D3d", 12}, {"C4h", 13}, {"C4v", 14}, {"D4", 
     15}, {"D2h", 16}, {"D2d", 17}, {"C6", 18}, {"C3i", 19}, {"C3h", 
     20}, {"C3v", 21}, {"D3", 22}, {"C4", 23}, {"S4", 24}, {"C2h", 
     25}, {"C2v", 26}, {"D2", 27}, {"Cs", 28}, {"Ci", 29}, {"C2", 
     30}, {"C1", 31}, {"C3", 32}
  };
  names = Map[GTConnectionPattern[#, GONotation -> nt] &, names3D];
  grp1 = Map[GTNotationText[#, GONotation -> nt] &, grp0];
  lab = {};
  Do[
  	 nn = GTNotationText[names3D[[i, 1]], GONotation -> nt];
     lab = Append[lab, nn -> {VertexLabels -> nn}]
  , {i, 1, 32}];
  con3D = {{"Oh", "Th"}, {"Oh", "Td"}, {"Oh", "O"}, {"Oh", 
     "D4h"}, {"D4h", "C4h"}, {"D4h", "C4v"}, {"D4h", "D4"}, {"D4h", 
     "D2h"}, {"D4h", "D2d"}, {"Th", "T"}, {"Th", "S6"}, {"Th", 
     "D2h"}, {"Td", "D2d"}, {"Td", "T"}, {"Td", "C3v"}, {"O", 
     "D4"}, {"O", "T"}, {"O", "D3"}, {"D6h", "D2h"}, {"D6h", 
     "C6h"}, {"D6h", "C6v"}, {"D6h", "D6"}, {"D6h", "D3h"}, {"D6h", 
     "D3d"}, {"C4h", "C4"}, {"C4h", "S4"}, {"C4h", "C2h"}, {"C4v", 
     "C4"}, {"C4v", "C2v"}, {"D4", "C4"}, {"D4", "D2"}, {"D2h", 
     "C2h"}, {"D2h", "C2v"}, {"D2h", "D2"}, {"C4", "C2"}, {"S4", 
     "C2"}, {"C2h", "Cs"}, {"C2h", "Ci"}, {"C2h", "C2"}, {"C2v", 
     "Cs"}, {"C2v", "C2"}, {"D2", "C2"}, {"Cs", "C1"}, {"Ci", 
     "C1"}, {"C2", "C1"}, {"C6h", "C2h"}, {"C6h", "C6"}, {"C6h", 
     "S6"}, {"C6h", "C3h"}, {"C6v", "C2v"}, {"C6v", "C6"}, {"C6v", 
     "C3v"}, {"D6", "D2"}, {"D6", "C6"}, {"D6", "D3"}, {"D3h", 
     "C2h"}, {"D3h", "C3h"}, {"D3h", "C3v"}, {"D3h", "D3"}, {"D3d", 
     "C2h"}, {"D3h", "C3i"}, {"D3d", "C3v"}, {"D3d", "D3"}, {"C6", 
     "C2"}, {"C6", "C3"}, {"C3", "C1"}, {"C3i", "Ci"}, {"C3i", 
     "C3"}, {"C3h", "Cs"}, {"C3h", "C3"}, {"C3v", "Cs"}, {"C3v", 
     "C3"}, {"D3", "C2"}, {"D3", "C3"}
  };
  pgg = Map[GTConnectionPattern[#, GONotation -> nt] &, con3D];
  sg = {};
  l = Length[pgg];
  xx = Map[ToString[#] &, pgg /. names];
  If[OptionValue[GOSelectCase] == "All",  
   Do[
     Do[
          str1 = ToString[grp1[[k]] /. names] <> " -> " <> "**";  
          str2 = "**" <> " -> " <> ToString[grp1[[m]] /. names];
          l1 = StringMatchQ[xx, str1];
          l2 = StringMatchQ[xx, str2];
          Do[
               Do[If[l1[[i]] && l2[[j]],
                   sg = Append[sg, pgg[[i]]],
                   None],
       {i, 1, l}],
      {j, 1, l}],
     {k, 1, Length[grp1]}]
    , {m, 1, Length[grp1]}];
   sg = Union[sg],
   Do[
        Do[
              
     str1 = ToString[grp1[[k]] /. names] <> " -> " <> 
       ToString[grp1[[m]] /. names];    
             l1 = StringMatchQ[xx, str1];
             Do[If[l1[[i]],
                 sg = Append[sg, pgg[[i]]],
                 None]
      , {i, 1, l}]
     , {k, 1, Length[grp1]}]
    , {m, 1, Length[grp1]}];
   sg = Union[sg]
   ];
  If[OptionValue[GOPlotStyle] == "Embedding", 
      HighlightGraph[ 
       	     Graph[pgg], sg, Properties -> lab,  
              VertexLabelStyle -> Directive[Black, 15]
          ],   
      LayeredGraphPlot[sg, GraphStyle -> OptionValue[GraphStyle]]
   ]
  ]
  

(*
***)

(****b* /GTGroupHierarchy
! NAME
!  GTGroupHierarchy
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  01/03/2013   : introduced in package
!  05/03/2013   : options introduced
!  02/01/2014   : VertexLabels in "Embedding"-Mode corrected
!  30/05/2015   : String input of point group names possible
!  January 2016 : correct labels in graph
!  28.12.2017   : check header
!  18.05.2018   : check docu
!  06.04.2020   : some prpblem appeared with the labeling of the graph, perhaps caused by the new Mathematica version
!                 problem solved by using the stanard option GraphStyle
!  07.04.2020   : documentation page changed
! USAGE
!  GTGroupHierarchy[group] gives a graph containing all subgroups of a group. Also all groups 
!  which contain group as a subgroup are included in the graph.
! INPUT
!  group   - point group
! OUTPUT
!  plot of the graph
! GTPack Options
!  GOPlotStyle  : "Normal"    -> no embedding
!               : "Embedding" -> Embedding in the whole connection graph
!  GONotation   : "SFL"       -> Schoenfliess
!               : "HM"        -> Herman-Mauguin    
! STandard OPTIONS
  GraphStyle    : "DiagramGold"
! GTPack MODULES
!  GTConnectionPattern, GTNotationText
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

GTGroupHierarchy[grp0_, OptionsPattern[]] := Module[{nt,names3D,names,grp1,lab,nn,i,con3D,str1,str2,pgg,l,l1,l2,sg,xx}, 
  nt = OptionValue[GONotation];
  names3D = {{"Oh", 1}, {"O", 2}, {"Th", 3}, {"Td", 4}, {"D6h", 
     5}, {"D4h", 6}, {"T", 7}, {"C6h", 8}, {"C6v", 9}, {"D6", 
     10}, {"D3h", 11}, {"D3d", 12}, {"C4h", 13}, {"C4v", 14}, {"D4", 
     15}, {"D2h", 16}, {"D2d", 17}, {"C6", 18}, {"S6", 19}, {"C3h", 
     20}, {"C3v", 21}, {"D3", 22}, {"C4", 23}, {"S4", 24}, {"C2h", 
     25}, {"C2v", 26}, {"D2", 27}, {"Cs", 28}, {"Ci", 29}, {"C2", 
     30}, {"C1", 31}, {"C3", 32}};
  names = Map[GTConnectionPattern[#, GONotation -> nt] &, names3D];
  grp1 = GTNotationText[grp0, GONotation -> nt];
  lab = {};
  Do[nn = GTNotationText[names3D[[i, 1]], GONotation -> nt];
   lab = Append[lab, nn -> {VertexLabels -> nn}], {i, 1, 32}];
  con3D = {{"Oh", "Th"}, {"Oh", "Td"}, {"Oh", "O"}, {"Oh", 
     "D4h"}, {"D4h", "C4h"}, {"D4h", "C4v"}, {"D4h", "D4"}, {"D4h", 
     "D2h"}, {"D4h", "D2d"}, {"Th", "T"}, {"Th", "S6"}, {"Th", 
     "D2h"}, {"Td", "D2d"}, {"Td", "T"}, {"Td", "C3v"}, {"O", 
     "D4"}, {"O", "T"}, {"O", "D3"}, {"D6h", "D2h"}, {"D6h", 
     "C6h"}, {"D6h", "C6v"}, {"D6h", "D6"}, {"D6h", "D3h"}, {"D6h", 
     "D3d"}, {"C4h", "C4"}, {"C4h", "S4"}, {"C4h", "C2h"}, {"C4v", 
     "C4"}, {"C4v", "C2v"}, {"D4", "C4"}, {"D4", "D2"}, {"D2h", 
     "C2h"}, {"D2h", "C2v"}, {"D2h", "D2"}, {"C4", "C2"}, {"S4", 
     "C2"}, {"C2h", "Cs"}, {"C2h", "Ci"}, {"C2h", "C2"}, {"C2v", 
     "Cs"}, {"C2v", "C2"}, {"D2", "C2"}, {"Cs", "C1"}, {"Ci", 
     "C1"}, {"C2", "C1"}, {"C6h", "C2h"}, {"C6h", "C6"}, {"C6h", 
     "S6"}, {"C6h", "C3h"}, {"C6v", "C2v"}, {"C6v", "C6"}, {"C6v", 
     "C3v"}, {"D6", "D2"}, {"D6", "C6"}, {"D6", "D3"}, {"D3h", 
     "C2h"}, {"D3h", "C3h"}, {"D3h", "C3v"}, {"D3h", "D3"}, {"D3d", 
     "C2h"}, {"D3h", "S6"}, {"D3d", "C3v"}, {"D3d", "D3"}, {"C6", 
     "C2"}, {"C6", "C3"}, {"C3", "C1"}, {"S6", "Ci"}, {"S6", 
     "C3"}, {"C3h", "Cs"}, {"C3h", "C3"}, {"C3v", "Cs"}, {"C3v", 
     "C3"}, {"D3", "C2"}, {"D3", "C3"}};
  str1 = ToString[grp1 /. names] <> " -> " <> "**";
  str2 = "**" <> " -> " <> ToString[grp1 /. names];
  pgg = Map[GTConnectionPattern[#, GONotation -> nt] &, con3D];
  l = Length[pgg]; xx = Map[ToString[#] &, pgg /. names];
  l1 = StringMatchQ[xx, str1]; l2 = StringMatchQ[xx, str2]; sg = {};
  Do[If[l1[[i]], sg = Append[sg, pgg[[i]]], Null], {i, 1, l}];
  Do[If[l2[[i]], sg = Append[sg, pgg[[i]]], Null], {i, 1, l}];
  If[OptionValue[GOPlotStyle] == "Embedding", 
      HighlightGraph[Graph[pgg], sg, Properties -> lab, 
      VertexLabelStyle -> Directive[Black, 10]], 
   LayeredGraphPlot[sg, GraphStyle-> OptionValue[GraphStyle]]
   ]
 ]

(*
***)

(****b* /GTPointGroups
! NAME
!  GTPointGroups
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  01/03/2013 : introduced in package
!  05/03/2013 : options introduced
!  28.12.2017 : check header and description
!  18.05.2018 : check docu
!  07.04.2020 : problems with graph representation, new option inserted, documentation changed
! USAGE
!  GTPointGroups[mode] gives a graph according to the subgroup relationships of 
!  the 32 point groups.
! INPUT
!  mode="Graph" gives a standard Graph representation,
!  mode="LayeredGraph" tries to order the nodes of the graph,
!  mode="Graph3D" shows a 3D plot of the graph
! OUTPUT
!  plot of the graph  
! GTPack OPTIONS
!  GODimension : space dimension 
!  GONotation  : Schoenfliess or Herman-Mauguin notation
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTConnectionPattern
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

GTPointGroups[mode_:"Graph",OptionsPattern[]]:=
	Module[{dim,nt,gt3d,gt2d,gtp,gtp1,gst},dim=OptionValue[GODimension];nt=OptionValue[GONotation];
		gst = OptionValue[GraphStyle];
		gt3d={{"Oh","Th"},{"Oh","Td"},{"Oh","O"},{"Oh","D4h"},{"D4h","C4h"},{"D4h","C4v"},{"D4h","D4"},
			{"D4h","D2h"},{"D4h","D2d"},{"Th","T"},{"Th","S6"},{"Th","D2h"},{"Td","D2d"},{"Td","T"},
			{"Td","C3v"},{"O","D4"},{"O","T"},{"O","D3"},{"D6h","D2h"},{"D6h","C6h"},{"D6h","C6v"},
			{"D6h","D6"},{"D6h","D3h"},{"D6h","D3d"},{"C4h","C4"},{"C4h","S4"},{"C4h","C2h"},{"C4v","C4"},
			{"C4v","C2v"},{"D4","C4"},{"D4","D2"},{"D2h","C2h"},{"D2h","C2v"},{"D2h","D2"},{"C4","C2"},
			{"S4","C2"},{"C2h","Cs"},{"C2h","Ci"},{"C2h","C2"},{"C2v","Cs"},{"C2v","C2"},{"D2","C2"},
			{"Cs","C1"},{"Ci","C1"},{"C2","C1"},{"C6h","C2h"},{"C6h","C6"},{"C6h","S6"},{"C6h","C3h"},
			{"C6v","C2v"},{"C6v","C6"},{"C6v","C3v"},{"D6","D2"},{"D6","C6"},{"D6","D3"},{"D3h","C2h"},
			{"D3h","C3h"},{"D3h","C3v"},{"D3h","D3"},{"D3d","C2h"},{"D3h","S6"},{"D3d","C3v"},{"D3d","D3"},
			{"C6","C2"},{"C6","C3"},{"C3","C1"},{"S6","Ci"},{"S6","C3"},{"C3h","Cs"},{"C3h","C3"},{"C3v","Cs"},
			{"C3v","C3"},{"D3","C2"},{"D3","C3"}};
		gt2d={{"C4v","C4"},{"C4v","C2v"},{"C4","C2"},{"C2v","Cs"},{"C2v","C2"},{"Cs","C1"},{"C2","C1"},
			{"C6v","C2v"},{"C6v","C6"},{"C6v","C3v"},{"C6","C2"},{"C6","C3"},{"C3v","Cs"},{"C3v","C3"},{"C3","C1"}};
		If[dim==3,gtp=gt3d,gtp=gt2d];gtp1=Map[GTConnectionPattern[#,GONotation->nt,GODimension->dim]&,gtp];
		Which[
			mode=="Graph",GraphPlot[gtp1,GraphStyle -> gst],
			mode=="LayeredGraph",LayeredGraphPlot[gtp1,GraphStyle -> gst],
			mode=="Graph3D",GraphPlot3D[gtp1,GraphStyle -> gst]
		]
]


(*
***)

(****b* /GTQuaternionQ
! NAME
!  GTQuaternionQ
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  5/14/13   : modified usage message
! 18.05.2018 : documentation checked
! USAGE
!  GTQuaternionQ[A] gives True if A is a quaternion, and False otherwise.
! INPUT
!  Object A
! OUTPUT
!  logical
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
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
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTQuaternionQ[a_] := If[(Length[a]==2)&&(NumberQ[N[a[[1]]]])&&(Length[a[[2]]] == 3)&&VectorQ[a[[2]]], True, False]
   
(*
***)   
   
   
   
(*-------------------------- Quaternion operations -----------------*)
   
(****b* /GTQAbs
! NAME
!  GTQAbs
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  * 14.05.2013 : modified usage message
!  * 28.12.2017 : check header
! USAGE
!  GTQAbs[q] gives the absolute value of the quaternion q.
! INPUT
!  Quaternion q in the form q = {a,{i,j,k}}
! OUTPUT
!  Absolute value (real number) of a quaternion.
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  
! DESCRIPTION
!  The absolute value of a quaternion is calculated by
! |latex \input{doc_tex/src/GTQAbs.tex}
! LITERATURE
!   Altman
! TODO
!  
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

 GTQAbs[a_] := Sqrt[Flatten[a] . Flatten[a]]

(*
***)

(****b* /GTQConjugate
! NAME
!  GTQConjugate
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  * 14.05.2013 : modified usage message
!  * 28.12.2017 : check header and documentation
! USAGE
!  GTQConjugate[q] gives the conjugate quaternion q* of the quaternion q.
! INPUT
!  Quaternion q in the form q = {a,{i,j,k}}
! OUTPUT
!  Conjugate Cuaternion in the form q* = {a,{-i,-j,-k}}
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  -  
! DESCRIPTION
!  For a quaternion and its conjugate quaternion the following relation holds 
!  |latex \input{doc_tex/src/GTQConjugate.tex}
! LITERATURE
!  Altman
! TODO
!  
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

 GTQConjugate[a_] := {a[[1]], -a[[2]]}

(*
***)

(****b* /GTQInverse
! NAME
!  GTQInverse
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  * 14.05.2013 : modified usage message
!  * 28.12.2017 : check header and documentation
! USAGE
!  GTQInverse[q] gives the inverse of the quaternion q.
! INPUT
!  Quaternion q in the form q = {a,{i,j,k}}
! OUTPUT
!  Inverse quaternion
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  GTQConjugate, GTQAbs
! DESCRIPTION
!  The inverse quaternion can be estimated by
!  |latex \input{doc_tex/src/GTQInverse.tex}
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

 GTQInverse[a_] := GTQConjugate[a]/GTQAbs[a]^2

(*
***)

(****b* /GTQMultiplication
! NAME
!  GTQMultiplication
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  * 04.02.2015    : programmed
!  * 28.12.2017 : check header and documentation
! USAGE
!  GTQMultiplication[q1,q2] calculates the multiplication of two quaternions.
! INPUT
!  Quaternion q1 and q2 in the form qi = {a,{i,j,k}}
! OUTPUT
!  quaternion.
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  - 
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

 GTQMultiplication[a_,b_] := {a[[1]] b[[1]] - a[[2]] . b[[2]], 
   a[[1]] b[[2]] + b[[1]] a[[2]] + Cross[a[[2]], b[[2]]]}
   
(*--------------------------- Symolic Operators -------------------*)
(*Quaternion multiplication*) 
 a_\[Diamond]b_ := GTQMultiplication[a,b]
(*
***)

(****b* /GTQPolar
! NAME
!  GTQPolar
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  * 14.05.2013 : modified usage message
!  * 28.12.2017 : check header and documentation
! USAGE
!  GTQPolar[q] gives the polar angle of the quaternion q.
! INPUT
!  Quaternion q in the form q = {a,{i,j,k}}
! OUTPUT
!  Polar Angle (real)
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  GTQAbs
! DESCRIPTION
!  The polar angle of a quaternion can be estimated by
!  |latex \input{doc_tex/src/GTQPolar.tex}
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
 GTQPolar[a_] := ArcCos[a[[1]]/GTQAbs[a]]
(*
***)

(****b* /GTSimplify
! NAME
!  GTSimplify
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  5/14/13    : modified usage message
!  28.12.2017 : check header 
! USAGE
!  GTSimplify[input] simplifies the input in a globally defined manner. 
!  internal function
! INPUT
!  input
! OUTPUT
!  simplified output
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  -  
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

GTSimplify[in_] := Expand[Simplify[in]];
   
(*
***) 

(****b* /GTSymbolQ
! NAME
!  GTSymbolQ
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  5/14/13 : modified usage message
!  28.12.2017 : check header and documentation
! USAGE
!  GTSymbolQ[A] gives True if A is a symbol, and False otherwise.
! INPUT
!  Object A
! OUTPUT
!  logical
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  - 
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

GTSymbolQ[a_] := If[Head[a] === Symbol,True,False]
   
(*
***) 

(****b* /GTWriteToFile
! NAME
!  GTWriteToFile
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  9/9/13     : 1st definition
!  10.02.2017 : bug corrected
!  28.12.2017 : check header and documentation
! USAGE
!  GTWriteToFile[object,file] writes object to the file.
! INPUT
!  object
!  file   - file name
! OUTPUT
!  -
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
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
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTWriteToFile[hop_,file_]:=Module[{s,file1,hop1},
	If[Head[file]===String,
	  file1=file; hop1=hop,
	  file1=hop;hop1=file
	 ];	
	s=OpenWrite[file1];Write[s,hop1];Close[s];
]


(*
***)




(****b* /GTReadFromFile
! NAME
!  GTReadFromFile
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  9/9/13 : 1st definition
!  28.8.16: FileExistQ added
!  28.12.2017 : check header and documentation, FileExistsQ was commented and now removed
! USAGE
!  GTReadFromFile[name] reads an object from the file name.
! INPUT
!  name
! OUTPUT
!  object
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES 
!  -
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
GTReadFromFile[file_] := Module[{s,hop},
	s=OpenRead[file];
	hop=Read[s];
	Close[s];
	Return[hop]
]



(*
***) 



(****b* /GTNeighborPlot
! NAME
!  GTNeighborPlot
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  1st version June 2014
!  28.12.2017 : check header and documentation
! USAGE
!  GTNeighborPlot[neighbors] is a tool to present the information about the input data to 
!  construct adjacency matrices in visual form.
! INPUT
!  neighbors -  {neighbors,atoms} 
!  neighbors is a matrix which contains the distances bewtween the 
!  different atoms which should be considered
!  
!  atoms is the list of atom names in the correct order with
!  respect to neighbors
! OUTPUT
!  nice plot of the matrix
! GTPack OPTIONS
!  -
! Standard OPTONS
!  -
! GTPack MODULES
!  -
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

GTNeighborPlot[nmat_] := Module[{tab, atoms, neighbors}, 
	atoms = nmat[[2]]; neighbors = nmat[[1]]; 
    tab = Join[{Join[{" "}, atoms]}, 
    Transpose[Join[{atoms}, neighbors]]]; 
    Grid[tab, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
    Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} -> GTCornerColor}]
]


(*
***) 


(****b* /GTClusterFilter
! NAME
!  GTClusterFilter
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  1st version June 2014
!  28.12.2017 : check header and documentation
!  18.05.2018 : check docu, link to GTClusterManipulate added.
! USAGE
!  GTClusterFilter[cluster,types] removes certain types of 
!  atoms from a given cluster.
! INPUT
!  cluster - cluster of atoms 
!  types   - types to remove
! OUTPUT
!  cluster
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  -
! DESCRIPTION
!
! LITERATURE
!
! TODO
!  
! PROBLEMS
!  where it is used? Module is not mentioned in the book.
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTClusterFilter[cl_, cutlist_] := 
 Module[{nat,cl1,i,ats}, nat = Length[cl]; cl1 = {};
  Do[
      If[Intersection[cutlist, {cl[[i, 2]]}] == {},
         cl1 = Append[cl1, cl[[i]]],
         None
      ]
   , {i, 1, nat}];
   Print[Length[cl1], " atoms"];
   ats=Transpose[cl1][[2]] // Tally;
   Print["Atoms in cluster : ",ats];
  Return[cl1]
  ]
  
  
(*
***) 

 
(****b* /GTCompactStore
! NAME
!  GTCompactStore
! AUTHOR
!  W. Hergert
! PACKAGE
!  Lattice.m
! MODIFICATION HISTORY
!  first version June 2014
!  28.12.2017 : check header and documentation
! USAGE
!  GTCompactStore[matrix] allows to store sparse matrices in a compact form.
! INPUT
!  matrix - square matrix (full/compact)
! OUTPUT
!  square matrix in opposite format to input (compact/full)
! GTPack OPTIONS
!  * GOCompact    - FALSE (standard) The INPUT is a full matrix, output will be compact storage 
!                 - TRUE  the matrix in compact form will be recasted to a full matrix
!  * GOMatrixType - "General" (standard) a general matrix is considered
!                 - "Symmetric" the matrix is symmetric, i.e. only one triangle has to be 
!                    considered. If only a triangle is stored and "Symmetric" is not set, you will
!                    get only this triangle, otherwise the whol matrix will be reconstructed.  
! Standard OPTIONS
!  -
! GTPack MODULES
!  -
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
GTCompactStore[mat_, OptionsPattern[]] := Module[{comp,symm,ndim,matc,ind,mato,i,j,l,k,n,elem},
  comp = OptionValue[GOCompact];
  symm = OptionValue[GOMatrixType];
  If[comp,
  (*--- reconstruct full matrix ---*)
     ndim = mat[[1]]; matc = mat[[2]]; ind = mat[[3]];
     mato = Table[0, {ndim}, {ndim}];
     Do[
        i = IntegerPart[(ind[[l]] - 1)/ndim] + 1;
        j = ind[[l]] - (i - 1)*ndim;
        mato[[j, i]] = matc[[l]]
     , {l, 1, Length[ind]}];
     If[symm == "Symmetric", 
        Do[
        	Do[ 
        	   mato[[i, j]] = mato[[j, i]]
        	, {j, i + 1, ndim}]
        , {i, 1, ndim}]],
  (*--- store compact ---*)
     ndim = Length[mat]; ind = {}; matc = {};
     Do[
     	If[symm == "General", 
     	   n = 1, 
     	   n = j
     	];
        Do[
           k = (j - 1)*ndim + i;
           elem = mat[[i, j]]; 
           If[elem === 0, 
           	  None,
              matc = Append[matc, elem]; ind = Append[ind, k]
           ]
        , {i, n, ndim}]
     , {j, 1, ndim}]; 
     mato = {ndim, matc, ind}
   ];
   Return[mato]
]

(*
***)

(****b* /GTSU2Matrix
! NAME
!  GTSU2Matrix
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Basic.m 
! MODIFICATION HISTORY
!  07/04/2014 : first commit
!  28.12.2017 : check header and documentation
!  18.05.2018 :check docu
! USAGE
!  GTSU2Matrix[angle,w] gives the SU(2) rotation matrix for a counterclockwise rotation 
!  around the 3D vector w.
! INPUT
!  angle, axis
! OUTPUT
!  2 dimensional spin rotation matrix
! GTPackOPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  - 
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
 
GTSU2Matrix[angle_, axis_] := (Cos[angle/2]*IdentityMatrix[2] - I*Sin[angle/2]*axis . Table[PauliMatrix[i], {i, 1, 3}])

(*
***)



(****b* /GTNotationText
! NAME
!  GTNotationText
! AUTHOR
! W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  January 2016 : first commit
!  28.12.2017 : check header 
! USAGE
!  GTNotationText[grp] constructes the text for the labels in the graphs
!  internal function
! INPUT
!  pointgroup name
! OUTPUT
!  text for label
! GTPack OPTIONS
!  GODimension : 2 or 3 for 2D or 3D point groups 
!  GONotation  : SFL or HM
!  GOVerbose   : if TRUE information about changes in names (double names)
!                will be printed
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -  
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
 

GTNotationText[grp_, 
  OptionsPattern[]] := 
 Module[{dim, names2D, names3D, names, text, text3D, text2D, symb,not,posf,pos,verb},
  dim = OptionValue[GODimension];
  not = OptionValue[GONotation];
  verb = OptionValue[GOVerbose];
  names3D = {"Oh", "O", "Th", "Td", "D6h", "D4h", "T", "C6h", "C6v", 
    "D6", "D3h", "D3d", "C4h", "C4v", "D4", "D2h", "D2d", "C6", "C3i",
     "C3h", "C3v", "D3", "C4", "S4", "C2h", "C2v", "D2", "Cs", "Ci", 
    "C2", "C1", "C3"};
  names2D = {"C4v", "C4", "C2v", "C2", "Cs", "C1", "C6v", "C6", "C3v",
     "C3"};
  text3D = {{"\!\(\*SubscriptBox[\(O\), \(h\)]\)", "O", 
     "\!\(\*SubscriptBox[\(T\), \(h\)]\)", 
     "\!\(\*SubscriptBox[\(T\), \(d\)]\)", 
     "\!\(\*SubscriptBox[\(D\), \(6  h\)]\)", 
     "\!\(\*SubscriptBox[\(D\), \(4  h\)]\)", "T", 
     "\!\(\*SubscriptBox[\(C\), \(6  h\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(6  v\)]\)", 
     "\!\(\*SubscriptBox[\(D\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(D\), \(3  h\)]\)", 
     "\!\(\*SubscriptBox[\(D\), \(3  d\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(4  h\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(4  v\)]\)", 
     "\!\(\*SubscriptBox[\(D\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(D\), \(2  h\)]\)", 
     "\!\(\*SubscriptBox[\(D\), \(2  d\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(3  i\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(3  h\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(3  v\)]\)", 
     "\!\(\*SubscriptBox[\(D\), \(3\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(S\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(2  h\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(2  v\)]\)", 
     "\!\(\*SubscriptBox[\(D\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(s\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(i\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(1\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(3\)]\)"}, 
     {"m\!\(\*OverscriptBox[\(3\), \(_\)]\)m", "432", "m3", 
     "\!\(\*OverscriptBox[\(4\), \(_\)]\)3m", "6/mmm", "4/mmm", "23", 
     "6/m", "6mm", "622", "\!\(\*OverscriptBox[\(6\), \(_\)]\)m2" , 
     "\!\(\*OverscriptBox[\(3\), \(_\)]\)m", "4/m", "4mm", "422", 
     "mmm", "42m", "6", "\!\(\*OverscriptBox[\(3\), \(_\)]\)", 
     "\!\(\*OverscriptBox[\(6\), \(_\)]\)", "3m", "32", "4", 
     "\!\(\*OverscriptBox[\(4\), \(_\)]\)", "2/m", "mm2", "222", "m", 
     "\!\(\*OverscriptBox[\(1\), \(_\)]\)", "2", "1", "3"}};
  text2D = {{"\!\(\*SubscriptBox[\(C\), \(4  v\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(4\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(2  v\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(s\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(1\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(6  v\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(6\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(3  v\)]\)", 
     "\!\(\*SubscriptBox[\(C\), \(3\)]\)"},
                   {"4mm", "4", "mm2", "2", "m", "1", "6mm", "6", 
     "3m", "3"}};
  (*--- transform into SFL string format ---*)
  
  symb = GTNotationToSFL[grp];
  (*--- check for non-standard names in SFL notation ---*)
  
  symb = GTDoubleNames[symb, GOVerbose -> verb];
  (*--- output of dcript form ---*)
  If[dim == 3,
        names = names3D; text = text3D,
        names = names2D; text = text2D
   ];
  posf = Flatten[Position[names, symb]];
  If[posf == {},
       Print["Error: GTNotationText: Symbol not found: ", symb]; 
   Abort[],
       pos = posf[[1]]
   ];
  If[not == "SFL",
       symb = text[[1, pos]],
       symb = text[[2, pos]]
   ];
  Return[symb]
  ]

(*
***)

(****b* /GTDoubleNames
! NAME
!  GTDoubleNames
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  January 2016 : first commit
!  28.12.2017   : check header 
! USAGE
!  GTDOobleNames[grp] checks if grp is a synonymous name to the GTPack standard. The standard
!  will be returned.
!  internal function
! INPUT
!  pointgroup name
! OUTPUT
!  pointgroup name  
! GTPack OPTIONS
!  GOVerbose: if TRUE information about changes in names (double names)
!             will be printed
! Standard OPTIONS
!   -
! GTPack MODULES
!   -
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
 

GTDoubleNames[grp_, OptionsPattern[]] := Module[{namesdbl,verb,pt,grp1,pos,new},
  namesdbl = {{"C1h", "S2", "S6", "S1", "V", "Vd", "Vh"}, {"Cs", "Ci",
      "C3i", "Cs", "D2", "D2d", "D2h"}};
  verb = OptionValue[GOVerbose];
  grp1 = GTNotationToSFL[grp];
  pos = Position[namesdbl, grp1];
  If[pos == {} ,
       new = grp1,
       If[pos[[1, 1]] == 1,
          pt = pos[[1, 2]];
          new = namesdbl[[2, pt]];
          If[verb,
               Print["GTDoubleNames: ", grp1, " is ", new],
              None
        ],
    new = grp1
     ]
   ];
  Return[new]
  ]

(*
***)


(****b* /GTSetTableColors
! NAME
!  GTSetTableColors
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  January 2016 : first commit
!  28.12.2017   : check header and documentation
! USAGE
!  GTSetTableColors[colors] csets a list of colors for tables. if the list is empty
!  the standard is set again.
! INPUT
!  list of colors
! OUTPUT
!  protected variables are set 
! GTPack OPTIONS
!  GTVerbose
! STandard OPTIONS
!  -
! GTPack MODULES
!  - 
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
*)

GTSetTableColors[colors_, OptionsPattern[]] := 
 Module[{m},
  Unprotect[GTBackGroundColor1,GTBackGroundColor2, GTCornerColor, GTDividerColor1, 
   GTDividerColor2];
   If[colors=={},
   (*--- reset to standard ---*)
     GTBackGroundColor1= Pink;
     GTBackGroundColor2= Gray;
     GTCornerColor     = Yellow;
     GTDividerColor1   = Red;
     GTDividerColor2   = Black,
   (*--- set new colors ---*)  	
     GTBackGroundColor1= colors[[1]];
     GTBackGroundColor2= colors[[2]];
     GTCornerColor     = colors[[3]];
     GTDividerColor1   = colors[[4]];
     GTDividerColor2   = colors[[5]]
  ];
  Protect[GTBackGroundColor1, GTBackGroundColor2, GTCornerColor, GTDividerColor1, 
   GTDividerColor2];
  If[OptionValue[GOVerbose],
  	If[colors=={},
  	   Print["Colors reset to standard."],
       Print["Colors in Tables changed to:"];
       m = {{"Background 1", colors[[1]]},
            {"Background 2", colors[[2]]},   
            {"Corner",       colors[[3]]}, 
            {"Dividers 1",   colors[[4]]}, 
            {"Dividers 2",   colors[[5]]}  
           };
       Print[TableForm[m]]
  	],
  	None
  ]
]

(*
***)


(****b* /GTPointInSquareQ
! NAME
!  GTPointInSquareQ
! AUTHOR
!  W. Hergert
! PACKAGE
!   Auxiliary.m
! MODIFICATION HISTORY
!   Oktober  2015  : first version
!   November 2015  : changed to logical function, point now relative to first vertex point
!   April    2016  : from Test.m to Auxiliary.m
!   28.12.2017     : check header
! USAGE
!  It is used to cut out of the honeycomb lattice those points which are lying
!  in the region to be rolled up to the nanotube.
!  (internal function)
! INPUT
!  o vert  - the vertices of the region
!  o point - point to be checked if it is inside
! OUTPUT
!  logical value: True if the point is inside 
! GTPack OPTIONS
!   -
! Standard OPTIONS
!   - 
! GTPack MODULES
!  -
! DESCRIPTION
!  It should be possible to solve this problem with region functions in a more
!  elegant way.
! LITERATURE
!  
! TODO
!  
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)


GTPointInSquareQ[vert_, point_] := Module[
  {xp, yp, v1, v2, glg, \[Alpha]1, \[Beta]1, \[Alpha], \[Beta], lsg, log1, log2, log},
   xp = point[[1]]-vert[[1,1]]; yp = point[[2]]-vert[[1,2]];
  (*--- first triangle ---*)
   v1 = vert[[2]] - vert[[1]]; 
   v2 = vert[[3]] - vert[[1]];
   glg = {xp == \[Alpha]1*v1[[1]] + \[Beta]1*v2[[1]], 
          yp == \[Alpha]1*v1[[2]] + \[Beta]1*v2[[2]]
         };
   lsg = Solve[glg, {\[Alpha]1, \[Beta]1}][[1]];
   \[Alpha] = \[Alpha]1 /. lsg; \[Beta] = \[Beta]1 /. lsg;
   If[0 <= \[Alpha] && \[Alpha] <= 1 && 
   	  0 <= \[Beta]  && \[Beta] <= 1  &&
      0 <= \[Alpha] + \[Beta] && \[Alpha] + \[Beta] <= 1, 
      log1 = True, 
      log1 = False
   ];
   (*--- second triangle  ---*)
   v1 = vert[[3]] - vert[[1]]; 
   v2 = vert[[4]] - vert[[1]];
   glg = {xp == \[Alpha]1*v1[[1]] + \[Beta]1*v2[[1]], 
          yp == \[Alpha]1*v1[[2]] + \[Beta]1*v2[[2]]
         };
   lsg = Solve[glg, {\[Alpha]1, \[Beta]1}][[1]];
   \[Alpha] = \[Alpha]1 /. lsg; \[Beta] = \[Beta]1 /. lsg;
   If[0 <= \[Alpha] && \[Alpha] <= 1 && 
   	  0 <= \[Beta]  && \[Beta] <= 1  &&
      0 <= \[Alpha] + \[Beta] && \[Alpha] + \[Beta] <= 1, 
      log2 = True, 
      log2 = False
   ]; 
   (*--- final check  ---*)
   If[log1 || log2, 
   	  log = True, 
   	  log = False
   ]; 
   Return[log]
]

(*
***)



(****b* /GTBlueRed
! NAME
!  GTBlueRed
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  April 2015 : implemented in Auxiliary.m
!  28.12.2017   : check header and documentation
! USAGE
!  GTBlueRed[x] defines a color function Blue - White - Red
! INPUT
!  x - mixing value
! OUTPUT
!  color
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  -
! DESCRIPTION
!  If x=0 the color iss blue, 1/2 results in white and 1 gives red finally
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

GTBlueRed[x_] :=Module[{}, 
    If[x < 1/2, 
       CMYKColor[(1 - 2 x), (1 - 2 x), 0, 0], 
       CMYKColor[0, 2 x - 1, 2 x - 1, 0]
    ]
]


(*
***)





(****b* /GTFermiSurfaceXSF
! NAME
!  GTFermiSurfaceXSF
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  May 2016     : first version
!  28.12.2017   : check header and documentation
!  18.05.2018   : check docu
! USAGE
!  GTFermiSurfaceXSF generates data for a Fermi surface plot by means of XCrysDEN
! INPUT
!  file   - output file name without extension
!  ham    - the hamiltonian, usually a TB Hamiltonian, but it shoud work also with pseudopotentials
!  kbasis - the basis vectors in k-space
!  ndel   - number of steps per direction in k-space
!  nbands - list of band numbers, which play a role 
!  fermi  - Fermi energy, this can be changed in XCrySDen
! OUTPUT
!  dat in file.bxsf
! GTPack OPTONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  -
! DESCRIPTION
!  
! LITERATURE
!  XCrySDen description of data formats (web page)
! TODO
!  
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTFermiSurfaceXSF[file_, ham_, kbasis_, ndel_, nbnds_, fermi_,  info_] := 
  Module[{nb,fname,head,tail,dela,delb,delc,nbas,ewlist,hdim,sel,i,j,k,kx,ky,kz,hc,ew,ewl,tt,bandspart,band,no,ni,nu,all},
  nb = Length[nbnds];
  fname = file <> ".bxsf"; 
  Print["Output in file: ", fname];
  (*--- heads and tails of files --*)
  head = {{"BEGIN_INFO"}}; 
  head = Flatten[Append[head, info], 1];
  head = Flatten[Append[head, 
  	     {{"   Fermi Energy:", fermi},
           "END_INFO",
           " ", 
           "BEGIN_BLOCK_BANDGRID_3D", 
           "Bandstructure", 
           "BANDGRID_3D_BANDS", 
            nb, 
           {ndel, ndel, ndel}, 
           {0.0, 0.0, 0.0},
           kbasis[[1]], kbasis[[2]], kbasis[[3]]
          }
          ], 1];
  tail = {"END_BANDGRID_3D", "END_BLOCK_BANDGRID_3D"};
  (*--- calculate the bandstructure ---*)
  {dela, delb, delc} = Map[Norm[#] &, kbasis]/(ndel - 1); 
  nbas = kbasis/Map[Norm[#] &, kbasis]; 
  ewlist = {};
  hdim = Length[ham]; 
  sel = Table[0, {hdim}]; 
  Do[
  	 sel[[nbnds[[i]]]] = 1
  , {i, 1, nb}];
  Do[
     Do[
        Do[
           {kx, ky, kz} =   (i - 1)*dela*nbas[[1]] 
                          + (j - 1)*delb*nbas[[2]] 
                          + (k - 1)*delc*nbas[[3]];
            hc     = ham /. {\[Xi] -> kx, \[Eta] -> ky, \[Zeta] -> kz};
            ew     = Sort[Eigenvalues[hc]// Chop,Less]; 
            ewl    = Pick[ew, sel, 1] ; 
            ewlist = Append[ewlist, ewl] ;
        , {k, 1, ndel}]
     , {j, 1, ndel}]
  , {i, 1, ndel}];
  (*--- output for xcrysden bxsf file ---*)
  tt = Transpose[ewlist]; 
  bandspart = {}; 
  no = Length[tt[[1]]];
  Do[
  	 band = Partition[tt[[i]], 6]; 
  	 ni   = Length[band];
  	 nu   = ni*6 + 1;
     bandspart = Append[bandspart, {"Band:", nbnds[[i]]}];
     Do[
        bandspart = Append[bandspart, band[[k]]]
     , {k, 1, ni}];
     If[no > nu,
        bandspart = Append[bandspart, Take[tt[[i]], {nu, no}]], None]
  , {i, 1, nb}];
  all = Flatten[{head, bandspart, tail}, 1];
  Export[fname, all, "Table"]
]

(*
***)

(****b* /GTFermiSurfaceXSF
! NAME
!  GTFermiSurfaceXSF
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  May 2016     : first version
!  28.12.2017   : check header and documentation
!  18.05.2018   : check docu
! USAGE
!  GTFermiSurfaceXSF generates data for a Fermi surface plot by means of XCrysDEN
! INPUT
!  file   - output file name without extension
!  ham    - the hamiltonian, usually a TB Hamiltonian, but it shoud work also with pseudopotentials
!  kbasis - the basis vectors in k-space
!  ndel   - number of steps per direction in k-space
!  nbands - list of band numbers, which play a role 
!  fermi  - Fermi energy, this can be changed in XCrySDen
! OUTPUT
!  dat in file.bxsf
! GTPack OPTONS
!  -
! Standard OPTIONS
!  -
! GTPack MODULES
!  -
! DESCRIPTION
!  
! LITERATURE
!  XCrySDen description of data formats (web page)
! TODO
!  
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTGetIrepMatrix[grp_, irep_, sym_, OptionsPattern[]] := 
 Module[{pos, lg, ld, imat, tab, tb, verb},
  (*---options---*)
  verb = OptionValue[GOVerbose];
  If[Head[sym] === List,
       lg = Length[sym];
       ld = IntegerPart[lg/10];
       imat = {};
       Do[
            pos = Flatten[Position[grp, sym[[i]]]][[1]];
           imat = Append[imat, irep[[pos]]]
    , {i, 1, lg}];
   tab = Join[{sym}, {MatrixForm[#] & /@ imat}] // Transpose;
   If[verb, Do[tb = Take[tab, {(i - 1)*10 + 1, i*10}] // Transpose;
     Print[
      Grid[tb, Frame -> All, 
       Dividers -> {Black, {2 -> GTDividerColor1}}, 
       Background -> {None, {1 -> GTBackGroundColor1}}]], {i, 1, 
      ld}];
    tb = Take[tab, {ld*10 + 1, lg}] // Transpose;
    Print[
     Grid[tb, Frame -> All, 
      Dividers -> {Black, {2 -> GTDividerColor1}}, 
      Background -> {None, {1 -> GTBackGroundColor1}}]], None];
   Return[Join[{sym}, {imat}] // Transpose], 
   pos = Flatten[Position[grp, sym]][[1]];
   tab = {{sym}, {MatrixForm[irep[[pos]]]}};
   If[verb, 
    Print[Grid[tab, Frame -> All, 
      Dividers -> {Black, {2 -> GTDividerColor1}}, 
      Background -> {None, {1 -> GTBackGroundColor1}}]], None];
   Return[{{sym, irep[[pos]]}}]];]





(*
***)

(****b* /GTVASPBandsPlot
! NAME
!  GTVASPBandsPlot
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  March 20182  : first version
!
! USAGE
!  GTVASPBandsPlot reads a VASP band structure file and plots the band structure.
!
! INPUT
!  file    - name of the VASP file
!  nkp     - total number of k-points used for the band structure plot
!  nbands  - number of bands
!  kpoints - contains the position and name of the symmetry points.
! OUTPUT
!  plot
! GTPack OPTONS
!  GOShift       - 0
!  GOFermiEnergy - 0
!  GOVerbose     - False 
! Standard OPTIONS
!  Joined     - True
!  PlotLabel  - "Band structure"
!  PlotRange  - All
!  PlotStyle  - Blue
!  FrameLabel - {" ","Energy (eV)"}
! GTPack MODULES
!  -
! DESCRIPTION
!  The file is not an EIGENVAL file, but a file which contains all bands one after each other.
!  The head contains information on the path in the BZ. This information has to be recasted
!  into the list kpoints in the form: {{k1,name1},{k2,name2},....} the names are the symbols for the 
!  symmetry points. The ki are the positions along the path.
! LITERATURE
!  -
! TODO
!  
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTVASPBandsPlot[file_, nkp_, nbands_, kpoints_, OptionsPattern[]] :=
 Module[{verb,join,plr,frame,pltl,plts,shft,s,ident,band,bs,eval,bt,xmax,emin,emax,bar,gs,gb,efe,efl,i,k},
  (*--- options ---*)
  verb  = OptionValue[GOVerbose];
  join  = OptionValue[Joined];
  plr   = OptionValue[PlotRange];
  frame = OptionValue[FrameLabel];
  pltl  = OptionValue[PlotLabel];
  plts  = OptionValue[PlotStyle];
  shft  = OptionValue[GOShift];
  efe  = OptionValue[GOFermiEnergy];
  (*--- read the data ---*)
  s = OpenRead[file]; bs = {};
  Do[ident = Read[file, String];
     If[verb,
        Print[ident],
        None
     ]; 
     band = {};
     Do[
        eval = Read[file, {Number, Number}];
        eval[[2]] = eval[[2]] - shft;
        band = Append[band, eval]
     , {k, 1, nkp}];
     bs = Append[bs, band]
  , {i, 1, nbands}];
  bt   = bs[[1]] // Transpose;
  xmax = Max[bt[[1]]]; 
  emin = Min[bt[[2]]];
  bt   = bs[[nbands]] // Transpose; 
  emax = Max[bt[[2]]];
  Print["Maximum Abscissa = ",xmax];
  bar  = Table[{{kpoints[[i, 1]], emin - 20}, {kpoints[[i, 1]], emax + 20}}, {i, 1, Length[kpoints]}];
  gb   = ListPlot[bs, Joined -> join, Frame -> True, PlotStyle -> plts, PlotLabel -> pltl,  PlotRange -> plr,
                      FrameLabel -> frame, FrameTicks -> {{Automatic, None}, {kpoints, kpoints}}, PlotLabel -> pltl
                 ];
  gs   = ListLinePlot[bar, PlotStyle -> Black, Ticks -> {None, Automatic}, PlotRange -> {{0, xmax}, {emin - 10, emax + 10}}];
  If[Abs[efe] > 0, 
     efl = Graphics[{Red, Line[{{0, efe}, {xmax, efe}}]}];
     Show[gb, gs, efl], 
     Show[gb, gs]
  ]
  ]

(*
***)


(****b* /GTVASPBands
! NAME
!  GTVASPBands
! AUTHOR
!  W. Hergert
! PACKAGE
!  Auxiliary.m 
! MODIFICATION HISTORY
!  March 20182  : first version
!
! USAGE
!  GTVASPBands reads an EIGENVAL file and transforms the band structure data
!  for a plot with GTBandsPlot
!
! INPUT
!  file    - name of the EIGENVAL-file
!  kpoints - contains the the symmetry points
! OUTPUT
!  data in GTPack form
! GTPack OPTONS
!  GOVerbose - True
!  GOSTore   - 0
! Standard OPTIONS
!  -
! GTPack MODULES
!  GTWriteToFile
! DESCRIPTION
!  The file is an EIGENVAL file.The list kpoints in the form: {{k1,name1},{k2,name2},....}
!  The names are the symbols for the  symmetry points. The ki are the corresponding
!  coordinates of the k-point.
! LITERATURE
!  -
! TODO
!  
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTVASPBands[file_, points_, OptionsPattern[]] :=Module[
	{verb,store,s,sys,kvecs,evals,nt,nkp,nbands,w,evec,k,j,i,bands,symp,del,tab},
 (*--- options ---*)
   verb  = OptionValue[GOVerbose];
   store = OptionValue[GOStore]; 
   s     = OpenRead[file];
   Read[file, {Number, Number, Number, Number}];
   Read[file, {Number, Number, Number, Number, Number}];
   Read[file, {Number}]; 
   Read[file, String];
   sys               = StringTrim[Read[file, String]];
   {nt, nkp, nbands} = Read[file, {Number, Number, Number}];
   If[verb,
      tab = {{"System ", " Number of k-points ",  "Number of bands  "}, {sys, nkp, nbands}} // Transpose;
      Print[Grid[tab, Frame -> All, Background -> {1 -> LightRed, None}]],
      None
   ];   
   kvecs = {}; 
   evals = {};
   Do[
      Read[file, String];
      {kx, ky, kz, w} = Read[file, {Number, Number, Number, Number}];
      kvecs           = Append[kvecs, {kx, ky, kz}]; 
      evec            = Table[0, {nbands}];
      Do[
         {k, evec[[j]], w} = Read[file, {Number, Number, Number}]
      , {j, 1, nbands}];
      evals = Append[evals, evec]
   , {i, 1, nkp}];
(*--- prepare data for GTPack ---*)
   bands = {};
   symp  = {};
   Do[
   	  If[i == 1,
         del = 0,
         del = del + Norm[kvecs[[i]] - kvecs[[i - 1]]]
      ];
      bands = Append[bands, {i, del, kvecs[[i]], evals[[i]]}]
   , {i, 1, nkp}];
   Do[
   	  Do[
   	  	 del = Norm[kvecs[[i]] - points[[j, 1]]];
         If[del < 10^(-6), 
         	symp = Append[symp, {i, points[[j, 2]]}], 
            None
         ]
      , {j, 1, Length[points]}]
   , {i, 1, nkp}];
   Return[{bands, symp}];
(*--- store ---*)
  If[Head[store] === String,
     GTWriteToFile[store, {bands, symp}],
     None
   ]   
]


(*
***)



(*-------------------------- Attributes ------------------------------*)

Attributes[GTQAbs]={Protected, ReadProtected}
Attributes[GTQConjugate]={Protected, ReadProtected}
Attributes[GTDiracMatrix]={Protected, ReadProtected}
Attributes[GTQInverse]={Protected, ReadProtected}
Attributes[GTQPolar]={Protected, ReadProtected}
Attributes[GTEulerAnglesQ]={Protected, ReadProtected}
Attributes[GTQuaternionQ]={Protected, ReadProtected}
Attributes[GTSymbolQ]={Protected, ReadProtected}
Attributes[GTCartesianSphericalHarmonicY]={Protected, ReadProtected}
(* Attributes[GTPointGroups]={Protected, ReadProtected}
Attributes[GTGroupHierarchy]={Protected, ReadProtected}
Attributes[GTGroupConnection]={Protected, ReadProtected}*)
Attributes[GTSU2Matrix]={Protected, ReadProtected}
Attributes[GTOverScript]={Protected, ReadProtected}
Attributes[GTSF]={Protected, ReadProtected}
Attributes[GTCayleyGraph1]={Protected, ReadProtected}

End[] (* End Private Context *)

EndPackage[]
