(*
!
! List of Changes:
!  13.06.2021   GTCluster    error message included if no atom in cluster
!  17.06.2021                stop if structure information has the wrong format
!  22.05.2021   GTShells     in comarision of distances with del an Abs[] was missing -> corrected
!
*)

(****m* /Lattice.m
!
! NAME
!  Lattice.m
! AUTHOR
!  W. Hergert, M. Geilhufe
! MODIFICATION HISTORY
!   * 15.08.2013 - initially created and documented
!   * 23.05.2018 - check of headers and documentation pages
! USAGE
!  This package contains modules for handling the direct lattice or the reciprocal lattice. 
!  
! GTPack MODULES
!
! --- Cluster construction and analysis  ---
!  
!  * GTCluster           - cluster construction
!  * GTShells            - reorder the cluster in shells
!  * GTAdjacencyMatrix   - constructs the adjacency matrix for a given cluster
!
! --- Lattices in real space and reciprocal space ---
!
!  * GTLatCluster        - spherical cluster of lattice points
!  * GTLatShells	         - reorder the lattice point cluster in shells
!
! --- Symmetry analysis of clusters in real space ---
!
!  * GTShellVectorsQlp   - minimal set of vectors Q_l^p representing the shell
!  * GTGroupGlp          - subgroups of the point group and their generators leaving Q_l^p invariant
!  * GTTransformToQlp    - symmetry operations which transform the shell vectors to the Q_l^p
!
! --- Reciprocal Lattice ---
!
!  * GTReciprocalBasis	- calculates reciprocal basis vectors from the basis of the direct lattice
!
! --- Brillouin zones ---
!
!  * GTBZPointMesh       - point mesh in irreducible part of the BZ
!  * GTBZPath            - standard path in the BZ for a given structure
!  * GTBZLines	         - k-points along a set of lines in the BZ
!  * GTGroupOfK          - group of the wavevector
!
! --- Brillouin zones and Wigner Seitz zones ---
!
!  * GTVoronoiCell       - construction of WS cell or BZ according to the basis
!
! DESCRIPTION
!  Lattice.m contains commands to handle clusters and lattices in direct or reciprocal space. The initial data for
!  the construction can be found in data bases. The constructed clusters or lattices are used mainly in the application
!  part: electronic structure, photonics, phonons.
!
!  Some special direct lattice commands are used in Egorov's tight binding method.
!
!  The reciprocal lattice commands include commands used for the calculation of band structures and 
!  densities of states, i.e. paths in the BZ or meshs in the BZ.
! LITERATURE
! 
***)


BeginPackage["GroupTheory`Lattice`",{"GroupTheory`Basic`","GroupTheory`Symbols`","GroupTheory`Auxiliary`","GroupTheory`Install`","GroupTheory`RepresentationTheory`"}]

(*-------------------------- Cluster construction and analysis  -----------------*)  
 GTAdjacencyMatrix         ::usage = "GTAdjacencyMatrix[\*StyleBox[\"cluster,neighbours\",\"TI\"]] constructs an adjacency matrix for \*StyleBox[\"cluster\",\"TI\"] with the information about interactions to include stored in \*StyleBox[\"neighbours\",\"TI\"]."
 GTCluster                 ::usage = "GTCluster[\*StyleBox[\"space group,radius\",\"TI\"]] constructs a spherical cluster with \*StyleBox[\"radius\",\"TI\"] for a given \*StyleBox[\"space group\",\"TI\"]."
 GTShells                  ::usage = "GTShells[\*StyleBox[\"cluster,basis,shells\",\"TI\"]] performs a reordering in \*StyleBox[\"shells\",\"TI\"] around the atoms in the \*StyleBox[\"basis\",\"TI\"]."

(*------------------------- Lattices in Real Space or reciprocal Space ----------*) 
 GTLatCluster              ::usage = "GTLatCluster[\*StyleBox[\"basis vectors,radius\",\"TI\"]] constructs a spherical cluster of \*StyleBox[\"radius\",\"TI\"] with lattice points definded by \*StyleBox[\"basis vectors\",\"TI\"]."
 GTLatShells               ::usage = "GTLatShells[\*StyleBox[\"cluster\",\"TI\"]] reorders a \*StyleBox[\"cluster\",\"TI\"] of lattice points in shells."

(*-------------------------- Symmetry Analysis of Real Structure Clusters --------------------------*)
 GTGroupGlp                ::usage = "GTGroupGlp[\*StyleBox[\"point group,vectors glp\",\"TI\"]] gives the subgroups \*Cell[BoxData[FormBox[SubsuperscriptBox[\"\[ScriptCapitalG]\",\"p\",\"l\"],TraditionalForm]],\"InlineMath\"] of a given point group and their generators, leaving the vectors \*Cell[BoxData[FormBox[SubsuperscriptBox[\"Q\",\"p\",\"l\"],TraditionalForm]],\"InlineMath\"] invariant."
 GTShellVectorsQlp         ::usage = "GTShellVectorsQlp[\*StyleBox[\"point group,shell vectors\",\"TI\"]] calculates the minimal set of vectors \*Cell[BoxData[FormBox[SubsuperscriptBox[\"Q\",\"p\",\"l\"],TraditionalForm]],\"InlineMath\"] representing the coordination spheres."
 GTTransformToQlp          ::usage = "GTTransformToQlp[\*StyleBox[\"point group,shell vectors,vectors qlp\",\"TI\"]] gives the symmetry operations of a \*StyleBox[\"point group\",\"TI\"] which transform the \*StyleBox[\"shell vectors\",\"TI\"] to the \*StyleBox[\"vectors \",\"TI\"] \*Cell[BoxData[FormBox[SubsuperscriptBox[\"Q\",\"p\",\"l\"],TraditionalForm]],\"InlineMath\"]."
  
(* -----------------------  Reciprocal Lattice ----------------------------------- *)
 GTReciprocalBasis         ::usage = "GTReciprocalBasis[\*StyleBox[\"lattice vectors\", \"TI\"]] calculates the reciprocal lattice vectors from a given set of \*StyleBox[\"lattice vectors\", \"TI\"]."

(*------------------------- Brillouin Zones --------------------------*)  
 GTBZLines                 ::usage = "GTBZLines[\*StyleBox[\"kpath,kplist\",\"TI\"]] generates \*StyleBox[\"k\",FontWeight->\"Bold\"]-points along a set of lines defined in the Brillouin zone. \*StyleBox[\"kplist\",\"TI\"] defines the number of points used for each part of the path."
 GTBZPath                  ::usage = "GTBZPath[\*StyleBox[\"structure\",\"TI\"]] generates a standard path in the Brillouin zone for a given \*StyleBox[\"structure\",\"TI\"]."
 GTBZPointMesh             ::usage = "GTBZPointMesh[\*StyleBox[\"n,a,structure\",\"TI\"]] gives a \*StyleBox[\"k\",FontWeight->\"Bold\"]-point mesh in the irreducible part of the Brillouin zone with a mesh refinement factor \*StyleBox[\"n\",\"TI\"] and a scaling factor \*StyleBox[\"a\",\"TI\"] for a \*StyleBox[\"structure\",\"TI\"]."
 GTBZMPBPointMesh          ::usage = "GTBZMPBPointMesh[file, kmesh,basis] exports a set of k-points in MPB data format. \*StyleBox[\"basis\",\"TI\"] is the reciprocal basis used in the MPB calculation."
 GTGroupOfK			       ::usage = "GTGroupOfK[\*StyleBox[\"group,\",\"TI\"]\*StyleBox[\"k\",\"TI\",FontWeight->\"Bold\"]\*StyleBox[\"vector,reciprocal lattice vectors\",\"TI\"]] gives the \*StyleBox[\"group\",\"TI\"] of the wave vector \*StyleBox[\"k\",\"TI\",FontWeight->\"Bold\"]."
  
 (*------------------------ Brilloin zones and Wigner-Seitz zones -----*)
 GTVoronoiCell             ::usage = "GTVoronoiCell[\*StyleBox[\"basis,cluster,BZpath\",\"TI\"]] constructs a Voronoi cell to a lattice with \*StyleBox[\"basis\",\"TI\"]. \*StyleBox[\"cluster\",\"TI\"] contains the data for the lattice construction. If a Brillouin zone is considered the path used in electronic structure calculations can be given by \*StyleBox[\"BZpath\",\"TI\"]."

(*--------------------------- Options ----------------------------*)
Options[GTAdjacencyMatrix] = {GOPlot -> False,GOSort->False}
Options[GTCluster]	       = {GOLattice -> {}, GOTbEquivalence -> True,GOPosition->0}
Options[GTGroupOfK]        = {GOTolerance->0.001, GOFast->GOFastValue}
Options[GTLatCluster]      = {GOVerbose->True}
Options[GTLatShells]       = {GOSort->False}
Options[GTShells]	       = {GOVerbose -> False, GOTbLattice ->{}, GOPosition -> "Absolute"}
Options[GTVoronoiCell]     = {GOVerbose -> False, GOBZPath -> True, GOOutput -> "Voronoi",VertexLabelStyle->Directive[Black, 20,Background -> Yellow]}
Options[GTBZMPBPointMesh]  = {GOVerbose -> True}
Options[GTGroupGlp]        = {GOVerbose -> True,GOCharTabs->True,GOIrepNotation->"Bouckaert"}

Begin["`Private`"] (* Begin Private Context *) 


(****h* /GTAdjacencyMatrix
! NAME
!  GTAdjacencyMatrix
! AUTHOR
!  W. Hergert
! PACKAGE
!  Lattice.m 
! MODIFICATION HISTORY
!  * 15.06.2014 - 1st version 
!  * 15.07.2014 - Introduction of GOSort
!  * 28.09.2017 - new Routine GTReplace
!  * 22.05.2018 - check header, documentation
! USAGE
!  GTAdjacencyMatrix[cluster,neighbours] constructs an adjacency matrix for cluster with the information about interactions to include stored in neighbors.
! INPUT
!  * cluster    - cluster (generated by GTCluster) 
!  * neighbours - list contains information about the neighbors to take into account 
!                 the first part of the list is square matrix, containing at position (i,j) 
!                 the distances of the two sorts which have to be taken into account. The second
!                 part is a list of atom names in the correct order with respect to the distance matrix.
!
! OUTPUT
!  adjacency matrix in compact storage form or plot of the matrix
!
! GTPack OPTIONS
!  * GOPlot:
!
!     - False : output of the adjacency matrix (standard)
!     - True  : MatrixPlot of the adjacency matrix
!  * GOSort:
!
!     - False : Names are not sorted (standard)
!     - True  : interaction "B,A" renamed to "A,B"
!
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTCompactStore, internal, without own header: GTReplace, GTSortNames

! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  * GTSortNames - during the construction of the adjacency matrix names like "A,B" and "B,A" occur.
!                  This is sometimes not nice for the setup of the Hamiltonian. The module rearranges the names 
!                  in canonical order, i.e. only "A,B" occurs if GOSort==True
!  * GTReplace   - the distances are replaced by the shell numbers.
!
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  * 27.09.2014 : The old version of GTReplace works with Position to find the distances. This doesn't work well
!                 if the distances are numerical values. The new version works with a tolerance.
! SOURCE
!--------------------------------------------------------------------------------
! 
*)
GTAdjacencyMatrix[cl_, nmat_, OptionsPattern[]] := 
         Module[{optp,nbas,mdist,dm,nat,ndi,l,name,i,j,ndim,mm,d,m0,mc,nm0,m1,opts}, 
         optp = OptionValue[GOPlot];
         opts = OptionValue[GOSort]; 
         (* set up the interactions for use *)
         nbas = Length[nmat[[2]]]; mdist = Max[nmat[[1]]]; dm = nmat[[1]]; 
         nat = Table[0, {nbas^2}]; ndi = Table[0, {nbas^2}]; l = 0;
         Do[
            Do[l = l + 1; 
               name = nmat[[2, i]] <> "," <> nmat[[2, j]];
               nat[[l]] = name;
               ndi[[l]] = {dm[[i, j]]}
            , {i, 1, nbas}]
         , {j, 1, nbas}];
         (*--- calculate distances ---*)
         ndim = Length[cl]; mm = Table[0, {ndim}, {ndim}];
         Do[
         	Do[
         	   d = Norm[cl[[i, 1]] - cl[[j, 1]]]; mm[[i, j]];
         (*--- cut out if distance larger than maximum distance ---*)
               If[d > mdist, None, 
                 (* mm[[i, j]] = {d, cl[[i, 2]] <> "," <> cl[[j, 2]]}*)
                 mm[[i, j]] = {d, cl[[j, 2]] <> "," <> cl[[i, 2]]}
               ]
            , {i, j + 1, ndim}]
         , {j, 1, ndim}];
         (*--- compact storage, replace distance by shell number ---*)
         mc = GTCompactStore[mm];
         m0 = mc[[2]]; nm0 = Length[m0];
         m1 = Map[GTReplace[#, nat, ndi] &, m0];
         (*--- canonical order for interaction names ---*)
         If[opts==True,
            mc[[2]] = Map[GTSortNames[#] &, m1],
            mc[[2]] =m1
         ];
         (*--- matrix plot or output of result in compact storage mode ---*)
         If[optp, 
         	m1 = Transpose[mc[[2]]][[1]]; mc[[2]] = m1;
            m0 = GTCompactStore[mc, GOCompact -> True, GOMatrixType -> "Symmetric"],
            None
         ];
         If[optp, 
         	MatrixPlot[m0], 
         	Return[mc]
         ]
]
GTReplace[matelm_, names_, dist_] := 
   Module[{mrp, posn, np, nd,eps,dist1,dd}, 
   	  mrp = matelm; eps = 10.0^(-5);
      posn = Position[names, mrp[[2]]] // Flatten;
      If[posn == {},
         Print["GTReplace : Name conflict",names," ",mrp[[2]]]; Abort[],
         np = posn[[1]];
         dist1 = dist[[np]] // Flatten; nd = Length[dist1];
         Do[
         	dd = Abs[dist1[[i]] - mrp[[1]]];
            If[dd <= eps,
               mrp[[1]] = i,
               None
            ]
         , {i, 1, nd}]
      ];
      Return[mrp]
]

GTSortNames[matelm_] :=Module[{me, ml, ne, nn,re}, 
	 me = matelm[[2]]; ne = matelm[[1]];
     If[ne == 0,
        re = {0, 0},
        ml = Sort[StringSplit[me, ","]]; 
        nn = ml[[1]] <> "," <> ml[[2]]; re = {ne, nn}
      ];
      Return[re]
]

(*
 ***)
 
 
(****h* /GTBZLines
! NAME
!  GTBZLines
! AUTHOR
!  W. Hergert
! PACKAGE
! Lattice.m 
! MODIFICATION HISTORY
!  * 15.11.2013 - 1st version 
!  * 23.05.2018 - check header and documentation
! USAGE
!  GTBZLines[kpath,points] generates k-points along a set of lines defined in the Brillouin zone. points defines the number of points used for each part of the path.
!
! INPUT
!  * kpath  -  Path in BZ. Can  be generated with GTBZPath.
!               Consists of a list of k-points, followed by a list of the 
!               corresponding labels.     
!          
!  * points - List of number of points used along the different lines.
!               Attention! Symmetry points at the beginning and end are
!               included.      
! OUTPUT
!  List contains all the kpoints for the band structure calculation and
!  the information about the labels.
!
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  k-points have to given in units (2 Pi)/a
!  Not only the k-points are generated. The second part of the output contains 
!  a list which contais the labels, but also the number of the kpoint, which 
!  is labeled accordingly : {{1,"Gamma"},{15,"X"},...}
!
!  It is developed from and old version of GTBandStructure, which was split into parts.
!
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  In the moment all routines containing the lattice are fixed more or less to cubic structures. We measure all lengths in absolute 
!  units in direct or reciprocal space. It would be helpful to introduce a measure in terms of lattice constants. 
! SOURCE
!--------------------------------------------------------------------------------
! 
*)


GTBZLines[kpath_,points_]:=Module[{kp,npath,kvec,npts,del,step,knorm,symb,kv,kpt,i,k},
	              kp=kpath[[1]];npath=Length[kp]-1;
                  If[npath==0,
                       npath=1,
                       None
                  ];
                  kvec={};kv=0;kpt=0;knorm=0;symb={};
                  Do[
                     npts=points[[i]];
                     del=(kp[[i+1]]-kp[[i]])/(npts-1);
                     step=Sqrt[del.del];kpt=kp[[i]];kv=kv+1;
                     kvec=Append[kvec,{kv,knorm,kpt}];
                     symb=Append[symb,{kv,kpath[[2,i]]}];
                     Do[
                        kv=kv+1;knorm=knorm+step;kpt=kpt+del;
                        kvec=Append[kvec,{kv,knorm,kpt}];
                     ,{k,2,npts-1}];
                     knorm=knorm+step;
                 ,{i,1,npath}];
                 symb=Append[symb,{kv+1,kpath[[2,npath+1]]}];
                 kvec=Append[kvec,{kv+1,knorm,kp[[npath+1]]}];
                 kvec=Append[{kvec},symb];
                 Return[kvec]
]


(*
***) 
 
(****h* /GTBZPath
! NAME
!  GTBZPath
! AUTHOR
!  W. Hergert
! PACKAGE
!  Lattice.m 
! MODIFICATION HISTORY
!  * 15.09.2013 - 1st version
!  * 17.02.2014 - modification, correction of units
!  * 23.05.2018 - check header and documentation
!  * 19.09.2018 - new message system, Switch instead a lot of If's
! USAGE
!  GTBZPath[structure] generates a standard path in the Brillouin zone for a given structure.
!
! INPUT
!  * strc  - structure as a string, the following structures are implemented:  fcc, bcc,sc,Honeycomb, SquareLattice
!     
! OUTPUT
!  path contains the high symmetry points of the path and their names
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  All k-points are defined with respect to units (2 Pi)/aa.  This is necessary to avoid discrepancies between the TB, Pseudopotential
!  and Photonics parts.
! LITERATURE
!  -
! TODO
!  implementation of more structures 
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

 GTBZPath::badarg = 
  "Argument `1` in GTBZPath[`1`]  is not one of the the following structures:";

GTBZPath[struc_:"Help"] := 
 Module[{path,plist},
  plist={"fcc","bcc","sc","Honeycomb","SquareLattice"};
  (*--- check input ---*)
  If[struc=="Help",
  	 Print["The following paths are implemented: ",plist];Return[],
  	 If[Intersection[plist, {struc}] == {},
        Message[GTBZPath::badarg, struc]; Print[plist]; Return[], 
        None
     ];
  ];	 
  Switch[struc,
  	"fcc", 
   path = {{{0, 0, 0}, {0, 0, 1}, {0, 1/2, 1}, {1/2, 1/2, 1/2}, {0, 0, 0},{3/4,3/4,0}, {1,
        1, 0}}, {"\[CapitalGamma]", "X", "W", "L", "\[CapitalGamma]","K", 
      "X"}},
 "bcc", 
   path = {{{0, 0, 0}, {0, 0, 1}, {0, 1/2, 1/2}, {0, 0, 0}, {1/2, 1/2, 1/2}, {0,
        1/2, 1/2}}, {"\[CapitalGamma]", "H", "N", "\[CapitalGamma]", "P", 
      "N"}},
   "sc", 
   path = {{{0, 0, 0}, {0, 0, 1/2}, {0, 1/2, 1/2}, {1/2, 1/2, 1/2}, {0, 0, 0}, {0,
        1/2, 1/2}, {1/2, 1/2, 1/2}}, {"\[CapitalGamma]", "X", "M", "R", 
      "\[CapitalGamma]", "M", "R"}}, 
 "Honeycomb", 
   path = {{{-2/3/Sqrt[3], 0, 0}, {0, 0, 0}, {1/2/Sqrt[3], -1/6, 
       0}, {2/3/Sqrt[3], 0, 0}}, {"K'", "\[CapitalGamma]", "M", "K"}}, 
   "SquareLattice", 
   path = {{{0, 0, 0}, {1/2, 0, 0}, {1/2, 1/2, 0}, {0, 0, 
    0}}, {"\[CapitalGamma]", "X", "M", "\[CapitalGamma]"}},
  _,
  Message[GTBZPath::badarg, struc]; Print[plist]; Return[]
  ];
  Return[path]
 ]
   
(*
***) 

(****h* /GTBZPointMesh
! NAME
!  GTBZPointMesh
! AUTHOR
!  W. Hergert
! PACKAGE
! Lattice.m 
! MODIFICATION HISTORY
!  * 15.09.2013 - 1st version 
!  * 23.05.2018 - check header and documentation
!  * 19.09.2018 - Message system changed, Switch instead of all the If's
!  * 08.02.2023 - nice table for help
! USAGE
!  GTBZPointMesh[n,a,structure] gives a k-point mesh in the irreducible part of the Brillouin zone with a mesh refinement factor n and a scaling factor a for a structure.
!
! INPUT
!  * n          - index defines 
!  * a          - List of the vectors in the different shells (output of GTShells)
!  * structure  - the options defines the structure:
!
!               - Help - no calculation, info about implemented cases (standard)
!               - line
!               - square
!               - square_full
!               - sc  (simple cubic)
!               - fcc (face centered cubic)
!               - bcc (body centered cubic)              
! OUTPUT
!  information about the number of calculated points, gives the k-point list back
!
! GTPack Options
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The  module is based on the k-meshes developed in the past from the Dresden group
!  for the DOS calculation with the tetrahedron method 
! LITERATURE
!  G. Lehmann, M. Taut, On the Numerical Calculation of the Density of States and Related Properties,
!  pss (b) 54, 469-477 (1972)
!
! TODO
!  implementation of hexagonal structure   
! RELEASE
!  1.0.1
! PROBLEMS
!  A more general method could be implemented starting from the LMTO codes. This should be connected with the implementation
!  of the thetrahedron method.
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTBZPointMesh::badstruc =  "Argument `1` in GTBZPointMesh not allowed, has to be one of the following:";


GTBZPointMesh[n_,a_,info_:"Help"]:=Module[{q,dk,qs,kx,ky,kz,n1,kxm,kym,kzm,n2,nkp,plist,test,dd,dx1,dx2,dy,i,j},	
	 plist = {"line", "square", "square_full","hexagonal", "sc", "fcc", "bcc"};
     test   = Intersection[{info}, plist];
     If[info=="Help",
     	None,
     	If[test == {},
           Message[GTBZPointMesh::badstruc, info]; Print[plist];Return[],
           None
        ]
     ];
    If[info=="Help",
       Print["k-point meshs implemented: line,square,square_full,hexagonal,sc,fcc,bcc (set third parameter accordingly)"];
       tab = {{"structure", "line ", "square"  , "square_full", "hexagonal", 
               "sc", "fcc", "bcc"}, {"k-points", "nkp = 2 n-1", "nkp = n(1+n)/2", 
               "nkp = n*n ", "nkp = n(1+n)/2 ", "nkp = n(1+n)(2+n)/6 ", 
               "nkp = 9 n^2+1+2*n*(8*n^2+7)/3", "nkp = (n+1)*(n+2)*(2n+3)/6 "}};
       Print[Grid[tab // Transpose, Alignment -> Left, Frame -> All, Background -> {{1 -> GTBackGroundColor1}, {1 -> 
             GTBackGroundColor2}, {1, 1} -> GTCornerColor}]];Return[], 
  (*     Print["line         : nkp = 2 n-1 "];
       Print["square       : nkp = n(1+n)/2 "];
       Print["square_full  : nkp = n*n "];
       Print["hexagonal    : nkp = n(1+n)/2 "];
       Print["sc           : nkp = n(1+n)(2+n)/6 "];
       Print["fcc          : nkp = 9 n^2+1+2*n*(8*n^2+7)/3"];
       Print["bcc          : nkp = (n+1)*(n+2)*(2n+3)/6 "];Return[],  *)
       None
    ];
    q={};
    Switch[info,
    (*--- line ---*)
    "line",
       nkp = 2*n - 1;
       dk=a/(n-1);
       Do[
          qs = {(kx - 1)*dk - a, 0, 0}; 
          q  = Append[q, qs]
       , {kx, 1, 2*n - 1}],      
    (*--- square ---*)
    "square",
       nkp = n*(n + 1)/2;
       dk=a/(n-1);
       Do[
          Do[
          	 qs={(kx-1)*dk,(ky-1)*dk,0};q=Append[q,qs]
          ,{ky,1,kx}]
        ,{kx,1,n}],
        (*--- square full ---*)
    "square_full",
       nkp = n*n;
       dk=a/(n-1);
       Do[
          Do[
          	 qs={(kx-1)*dk,(ky-1)*dk,0};q=Append[q,qs]
          ,{ky,1,n}]
        ,{kx,1,n}],
    (*--- sc ---*)
    "sc",
       nkp = n (1 + n) (2 + n)/6;
       dk=a/2/(n-1);
       Do[
          Do[
             Do[
                qs={(kx-1)*dk,(ky-1)*dk,(kz-1)*dk};
                q=Append[q,qs]
             ,{kz,1,kx}]
          ,{ky,kx,n}]
       ,{kx,1,n}],
    (*--- bcc ---*)
    "bcc",
       nkp = (n + 1)*(n + 2)*(2 n + 3)/6;
       n1=2*n;kzm=n+1;
       Do[kym=n-kz+2;
          Do[kxm=2*(n-kz-ky)+5;
             Do[
                qs={(kx+ky+kz-3)/n1*a,(ky+kz-2)/n1*a,(kz-1)/n1*a};
                q=Append[q,qs]
             ,{kx,1,kxm}]
          ,{ky,1,kym}]
       ,{kz,1,kzm}],
    (*--- fcc ---*)
     "fcc", 
       nkp = 9 n^2 + 1 + 2*n*(8*n^2 + 7)/3;
       n1=4*n;kzm=2*n+1;
       Do[kym=3*n+2-kz-IntegerPart[kz/2];
          Do[kxm=n1-kz-ky+3;n2=6*n+6-3*kz-2*ky;
             If[n2<kxm,kxm=n2,None];
             Do[
                qs={(ky+kz-2)/n1*a,(kx+ky+kz-3)/n1*a,(kz-1)/n1*a};
                q=Append[q,qs]
             ,{kx,1,kxm}]
          ,{ky,1,kym}]
       ,{kz,1,kzm}],
           (*--- hexagonal ---*)
     "hexagonal",
       nkp = n*(n+1)/2;
       dx1 = a/(n - 1); 
       dx2 = a*0.25/(n - 1); 
       dy  = a*Sqrt[3.]/4./(n - 1);
       q  = {};
       Do[
          dd = (i - 1)*dx1;
          Do[
             kx = dd - (j - 1)*dx2;
             ky = dy*(j - 1);
             q = Append[q, {kx,ky,0}]
         , {j, 1, i}]
       , {i, 1, n}],
       _,
       Message[GTBZPointMesh::badstruc, info]; Return[]
    ];
    Print[nkp, " k-points calculated"]; Return[q]
]

(*
***) 



 

(****h* /GTCluster
! NAME
!  GTCluster
! AUTHOR
!  W. Hergert
! PACKAGE
! Lattice.m 
! MODIFICATION HISTORY
!  * 15.08.2013 - 1st version 
!  * 01.02.2014 - Extension to 2D case 
!  * 15.06.2014 - Handling auf equivalent atoms in the basis
!  * 22.12.2016 - FullSimplify in the evaluation of vec added, avoids problems with accuracy
!  * 23.05.2018 - check of header and documentation
!  * 13.06.2021 - error message if no atoms in cluster
!  * 17.06.2021 - stop if structure information has the wrong format
! USAGE
!  GTCluster[space group,radius] constructs a spherical cluster with radius for a given space group.
! INPUT
!  * space group  - information about a structure in the format of CrystalStructure.m
!  * radius       - cutoff distance
! OUTPUT
!  cluster of atoms
!
! GTPack OPTIONS
!  * GOLattice:
!
!       -  {} (standard) -> original lattice constant information is used
!       -  [a->1, b->2, c->3} , list of substiution rules for the lattice parameters
!  * GOTBEquivalence: 
!
!       - True  : no special investigation with respect to equivalent atoms int the basis (Standard)
!       - False : rename atoms in the basis to make all atoms inequivalent      
! * GOPosition -  shifts the whole cluster, standard is 0, i.e. no shift
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!-
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

intChangeCoord::coord =  "Type of coordinates is `1`. It has to be R or C."
                         

intChangeCoord[struc_, OptionsPattern[{GOVerbose->False}]] := Module[
	{ln, st, struc1, pos, gv, nb, apos, npos, apos1, verb}, 
     verb   = OptionValue[GOVerbose];
     struc1 = struc;
  (*---check type of data set---*)ln = Length[struc];
     If[ln == 9,
        st = struc1[[9]];
        If[Intersection[{st}, {"R", "C"}] == {},
           Message[intChangeCoord::coord, st]; Abort[],
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
GTCluster::error = "There are no atoms in the cluster. Increase radius."
GTCluster::struc = "Wrong format of structure file. 3d-lattice vectors necessary."

GTCluster[spcgrp_, rc_, OptionsPattern[]] :=  Module[
	{lat, nb, bas, rule, nc, fac, vec, ib, i, j, k, rv, nc1,tbas,mult,tb,nat,sub,t,shift,ats,spcgr1}, 
	 (*--- Transformation to Cartesioan  basis vectors ---*)  
     spcgr1= intChangeCoord[spcgrp];
     lat   = spcgr1[[6]];
      If[Length[lat] == 3,
       None,   
       Message[GTCluster::struc];Return[]
      ];

     bas = spcgr1[[7]]; nb = Length[bas];
     (*--- check for nonequivalent basis atoms ---*)
     shift=OptionValue[GOPosition];
     If[OptionValue[GOTbEquivalence] == True, 
     	None,
        tbas = Transpose[bas];
        mult = Tally[tbas[[2]]]; nat = Length[mult];
        tb = {}; sub = {}; 
        Do[
           t  = Table[mult[[i, 1]] <> ToString[j], {j, 1, mult[[i, 2]]}]; 
           tb = Append[tb, t] // Flatten;
        , {i, 1, nat}]; 
        tbas[[2]] = tb; bas = Transpose[tbas];
     ];
     (*--- cluster construction ---*)
     If[OptionValue[GOLattice] == {}, 
     	rule = spcgrp[[8]], 
        rule = OptionValue[GOLattice]
     ];
     lat = lat /. rule;
     bas = bas /. rule;
     fac = Min[Select[Map[Norm[#] &, lat], # > 0 &]];
     nc = IntegerPart[rc/fac] + 3; rv = {};
     If[Min[Map[Norm[#] &, lat]] == 0, 
        nc1 = 0, 
        nc1 = nc
     ];
     Do[
        Do[
       	   Do[
              Do[
                 vec = i*lat[[1]] + j*lat[[2]] + k*lat[[3]] + bas[[ib, 1]]-shift//FullSimplify;
                 If[Norm[vec] <= rc, 
                    rv = Append[rv, {vec, bas[[ib, 2]]}], 
                    None
                 ]
              , {ib, 1, nb}]
           , {i, -nc, nc}]
        , {j, -nc, nc}]
     , {k, -nc1, nc1}];
     If[Length[rv]==0,
     	Message[GTCluster::error];Return[],
     	None
     ];
     Print[Length[rv], " atoms"];
     ats=Transpose[rv][[2]] // Tally;
     Print["Atoms in cluster: ",ats];
     Return[rv]
]
 
(*
***) 
		



(****h* /GTGroupGlp
! NAME
!  GTGroupGlp
! AUTHOR
!  W. Hergert
! PACKAGE
! Lattice.m 
! MODIFICATION HISTORY
!  * 13.08.2013 - 1st version
!  * 01.02.2017 - extended to several vectors per shell
!  * 23.05.2018 - check of header and documentation
! USAGE
!  GTGroupGlp[point group,vectors glp] gives the subgroups G_l^p of a given point group and their generators, leaving the vectors Q_l^p invariant.
!  All vectors of a shell can be represented by a minimum number of vectors Q_l^p. We find the groups G_l^p
!  consisting of symmetry tranformations of the point group that leave Q_l^p invariant.
! INPUT
!  * pg    - Operations of the point group
!  * qvec  - set of vectors Q_l^p (output of GTShellVectorsQlp)          
! OUTPUT
!  The point groups G_l^p and a list of their generators
! GTPack OPTIONS
!  * GOVerbose       - level of verbosity (option for GTCharacterTable)
!  * GOIrepNotation  - notation for IREPs (option for GTCharacterTable)
!  * GTCharTabs      - print of character tables or not
! STANDARD OPTIONS
!
! GTPack MODULES
!   * GTGetMatrix
!   * GTCharacterTable
!   * GTGenerators
!   * GTGetSymbol
!   * GTGroupQ
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  R.F. Egorov et. al. , phys. stat. sol. 26, 391 (1968), p. 396
! TODO
!  The options are not organized very clever.  There is not a big difference between GOVerbose and 
!  GOCharTabs. Perhaps a simplifications is possible. The structure of the output seems to be also a bit weired.
!  Can the lists be simplified?
! RELEASE
!  1.0.0
! PROBLEMS
!  might be more vectors per shell / extended to include several vectors per shell 1.2.2017
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

 GTGroupGlp[pg_,qvec_,OptionsPattern[]]:=Module[
 	{qv,gen,grps,op,gpt,n,nop,gp1,is,ns,verb,sgrps,sgens,k,nvec,ls,ctabs,notat,grp,tab,
 	 tabs,stabs},
  (*--- options ---*)
 	 verb  = OptionValue[GOVerbose];
 	 ctabs = OptionValue[GOCharTabs];
     notat = OptionValue[GOIrepNotation];
  (*---    ---*)    
	 nop = Length[pg];
	 ns  = Length[qvec];
	 gen = {};
	 grps= {};
	 tabs= {};
  (*--- loop over shells ---*)	 
     Do[
     	sgrps = {};
     	sgens = {};
     	stabs = {};
  (*--- loop over vectors per shell ---*)
        nvec=Length[qvec[[is]]];
        Do[	 
           gpt   = {}; 
  (*---  loop over point group ---*)          	
           Do[
      	      op = GTGetMatrix[pg[[n]]];
      	      qv = qvec[[is,1]][[1]];
      	      If[op.qv==qv,
      	 	     gpt=Append[gpt,op],
      	 	     None
      	      ];
           ,{n,1,nop}];
           grp = Map[GTGetSymbol[#] &, gpt];
           gp1 = GTGenerators[grp];
           sgens  = Append[sgens,gp1];
           sgrps  = Append[sgrps,grp];
           If[ctabs,
              tab = GTCharacterTable[grp, GOIrepNotation -> notat, GOVerbose -> verb];
              stabs = Append[stabs, tab],
              None
           ]
        ,{k,1,nvec}];
        grps = Append[grps,sgrps];
        gen  = Append[gen,sgens];
        tabs = Append[tabs, stabs]
     ,{is,1,ns}];
     If[verb,
        ls = {};
        Do[
           ls = Append[ls, Map[GTGroupQ[#] &, grps[[is]]]]//Flatten
        , {is, 1, ns}];
        Print["Check if the sets are groups: ", ls],
        None
     ];
     If[ctabs,
        Return[{grps, gen, tabs}],
        Return[{grps, gen}]
     ]
 ]

(*
***) 

(****h* /GTGroupOfK
! NAME
!  GTGroupOfK
! AUTHOR
!  M. Geilhufe
! PACKAGE
!   Lattice.m
! MODIFICATION HISTORY
!  * 12.03.2013 - first version
!  * 26.03.2013 - allow 2-dimensional lattices
!  * 23.08.2016 - better compatibility with space group elements, concentrate on rotational part only, use GTGetRotationMatrix instead of switching standard representation
!  * 23.05.2018 - check of header and documentation
! USAGE
!  GTGroupOfK[group,kvector,reciprocal lattice vectors] gives the group of the wave vector k.
! INPUT
!  * full point group
!  * k-vector
!  * reciprocal lattice-vectors
! OUTPUT
!  group of k
! GTPack OPTIONS
!  * GOFast - toggle time consuming checks in input checking
!  * GOTolerance - maximal deviation to consider numbers as equal
! STANDARD OPTIONS
!  - 
! GTPack MODULES
!  * GTGroupQ 
!  * GTGetRotationMatrix  
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  Dresselhaus
! TODO
!  Write the documentation notebooks. Seems to work, but check the module!
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! ERROR MESSAGES
! * $ErrorNoGroup
! * $ErrorDimensionKvec
!--------------------------------------------------------------------------------
! SOURCE
*)
GTGroupOfK[grp0_, vec0_, lattice0_, OptionsPattern[]] := Module[{grp, go, trvec, eq, sol, add, grpout, c, dp,lattice, vec},
   
   If[Not[Length[vec0]==Length[lattice0]],Print[$ErrorDimensionKvec];Abort[]];
  
  (*--- Transform the input to O(3) matrices ---*)
  (* concentrate on rotational part only for space group elements*)
  go = Length[grp0];
  grp = Table[GTGetRotationMatrix[grp0[[i]]],{i,1,go}];
  
  (*--- this is the part of GOFast ---*) 
  If[OptionValue[GOFast], None, If[GTGroupQ[grp], None, Print[$ErrorNoGroup]; Return[]]];

  (*--- check if the lattice is 2-dimensional ---*)
  (* make it 3-diemnsional *)
  dp = Length[lattice0];
  If[dp == 2,
   lattice = Table[Flatten[{lattice0[[i, 1 ;; 2]], 0}], {i, 1, dp}];
   vec = Flatten[{vec0[[1 ;; 2]], 0}],
   lattice = lattice0;
   vec = vec0];
   
  (*--- Main algortihm ---*)
  trvec = Table[Inverse[grp[[i]]].vec, {i, 1, go}];
  eq = Table[trvec[[i]] == vec + Sum[c[j]*lattice[[j]], {j, 1, dp}], {i, 1,go}];
  sol = Rationalize[Round[Table[Table[c[j], {j, 1, dp}] /. Flatten[Solve[eq[[i]], Table[c[j], {j, 1, dp}]]], {i, 1, go}],OptionValue[GOTolerance]]];
  If[dp == 3,
   add = Table[
     IntegerQ[sol[[i, 1]]] && IntegerQ[sol[[i, 2]]] && 
      IntegerQ[sol[[i, 3]]], {i, 1, go}],
   add = Table[
      IntegerQ[sol[[i, 1]]] && IntegerQ[sol[[i, 2]]], {i, 1, go}];
   ];
  grpout = {};
  Do[If[add[[i]], grpout = Append[grpout, grp0[[i]]],None], {i, 1, go}];
  
  
  Return[grpout]
 ]

(*
***)


(****h* /GTLatCluster 
! NAME
!  GTLatCluster
! AUTHOR
!  W. Hergert
! PACKAGE
!  Lattice.m 
! MODIFICATION HISTORY
!  * 23.10.2013 - 1st version 
!  * 15.11.2013 - output changed to 3D vectors
!  * 15.02.2016 - GOVerbose 
!  * 23.05.2018 - check of header and documentation
! USAGE
!  GTLatCluster[basis vectors,radius] constructs a spherical cluster of radius with lattice points definded by basis vectors.   
! INPUT
!  * basis vectors  - lattice vectors
!  * radius         - cutoff distance
!
! OUTPUT
!  list of lattice points in a sphere of radius rc
! GTPack OPTIONS
!  * GOVerbose - verbosity level
! STANDARD OPTIONS
!
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The basis vectors may be provide as 3 vectors for a 3D problem containing x,y,z components.
!  For a 2D problem the vectors are given as two vectors with x an y components only.
!  The constructed vectors ccontain always x,y,z components.
! LITERATURE
!  -
! TODO
!  The error handling can be definitely improved. Perhaps this should be done if the entire message list is
!  rewritten.
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTLatCluster[lat_,rc_,OptionsPattern[]]:=Module[{fac,nc,rv,vec,i,j,k,err,nl},err=False;nl=Length[lat];
	                     If[nl==3,
	                     	If[Length[lat[[1]]]==3,
	                     	   fac=Min[Norm[lat[[1]]],Norm[lat[[2]]],Norm[lat[[3]]]];lat1=lat,
	                     	   err=True
	                     	],
	                     	None
	                      ];  	
	                      If[nl==2,
	                        If[Length[lat[[1]]]==2,
	                     	   fac=Min[Norm[lat[[1]]],Norm[lat[[2]]]];lat1=Append[lat,{0,0}],
	                     	   err=True
	                     	],
	                     	None
	                      ];  	
                         nc=IntegerPart[rc/fac]+6;
                         rv={};
                         If[err,
                         	Print["Error : Wrong Lattice vectors!"],
                            Do[
                              Do[
                                Do[
                                   vec=i*lat1[[1]]+j*lat1[[2]]+k*lat1[[3]];
                                   If[Norm[vec]<=rc,
                                      If[nl==3,
                                 	     rv=Append[rv,vec],
                                 	     rv=Append[rv,Flatten[{vec,0}]]
                                 	  ],	
                                	None
                                  ]
                               ,{k,-nc,nc}]
                             ,{j,-nc,nc}]
                           ,{i,-nc,nc}];
                          rv=DeleteDuplicates[rv];
                          If[OptionValue[GOVerbose],
                             Print[Length[rv]," lattice vectors"],
                             None
                          ];
                          Return[rv]
                         ]  
]

(*
 ***)
 
  
(****h* /GTLatShells
! NAME
!  GTLatShells
! AUTHOR
!  W. Hergert
! PACKAGE
!  Lattice.m 
! MODIFICATION HISTORY
!  * 15.10.2013 - 1st version 
!  * 23.05.2018 - check of header and documentation
! USAGE
!  GTLatShells[cluster] reorders a cluster of lattice points in shells.
! INPUT
!  cluster - cluster of lattices points (generated by GTLatCluster) 
! OUTPUT
!  cluster reordered in shells, or in canonical order
! GTPack OPTIONS
!  * GOSort:
!
!     - False - order in shells (standard)
!     - True  - canonical order only 
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
! SOURCE
!--------------------------------------------------------------------------------
! 
*)
GTLatShells[clat_,OptionsPattern[]]:=Module[{gv1,gv2,sh,ts,k,i},
                   If[OptionValue[GOSort],
                   	  sh=Sort[clat,Norm[#1]<Norm[#2]&],
                   	 gv1=Sort[clat,Norm[#1]<Norm[#2]&];
                     gv2=Map[Norm[#]&,gv1];sh={};gv2=Sort[Union[gv2],#1<#2&];
                     Do[ts={};
                        Do[
                           If[Norm[clat[[k]]]==gv2[[i]],
                              ts=Append[ts,clat[[k]]],
                              None
                           ]
                        ,{k,1,Length[clat]}];
                        sh=Append[sh,ts]
                     ,{i,1,Length[gv2]}]
                   ]; 
                  Return[sh]
]

(*
 ***)

(****h* /GTReciprocalBasis
! NAME
!  GTReciprocalBasis
! AUTHOR
!  M. Geilhufe
! PACKAGE
!  Lattice.m 
! MODIFICATION HISTORY
!  * 23.05.2018 : check of header and documentation
! USAGE
!  GTReciprocalBasis[lattice vectors] calculates the reciprocal lattice vectors from a given set of lattice vectors.
! INPUT
!  set of lattice vectors
! OUTPUT
!  set of reciprocal lattice vectors
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
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTReciprocalBasis[latticeinp_] := 
 Module[{norm, lattice, latdim, rec, recout},
  latdim = Length[latticeinp];
  If[latdim==2,lattice={{latticeinp[[1,1]],latticeinp[[1,2]],0},{latticeinp[[2,1]],latticeinp[[2,2]],0},{0,0,1}},
  	If[latdim==3,lattice=latticeinp,Print["Invalid Input!"]]
  ];
  norm = (lattice[[1]]\[Cross]lattice[[2]]).lattice[[3]]; 
  rec = 2 \[Pi]/norm {lattice[[2]]\[Cross]lattice[[3]], 
    lattice[[3]]\[Cross]lattice[[1]], 
    lattice[[1]]\[Cross]lattice[[2]]};
  If[latdim==2,recout={{rec[[1,1]],rec[[1,2]]},{rec[[2,1]],rec[[2,2]]}},
  	If[latdim==3,recout=rec,None]
  ];
  Return[recout]
  ]
    
(*
 ***)
 
 
  
 
 
(****h* /GTShells
! NAME
!  GTShells
! AUTHOR
!  W. Hergert
! PACKAGE
! Lattice.m 
! MODIFICATION HISTORY
!   * 05.08.2013 - 1st version
!   * 10.08.2013 - modified to create also output for construction of Hamiltonians
!   * 28.02.2014 - totally new version to include more than two basis atoms
!   * 05.06.2014 - handling of equivalent atoms in the basis is included
!   * 22.12.2016 - changes connected to test of hexagonal structures, small change in finding distances to get simpler exressions
!   * 23.05.2018 - check header and documentation
!   * 22.05.2021 - in comarision of distances with del an Abs[] was missing -> corrected
! USAGE
!  GTShells[cluster,basis,shells] performs a reordering in shells around the atoms in the basis.
! INPUT
!  cluster     - cluster of atoms (output of GTCluster)
!  basis       - the basis list like in the spcgrp files 
!  shells      - number of neighbor shells around the basis atoms to contruct 
!
! OPTIONS
!  *GOVerbose:
!      - True  - Output of information in table form (standard)
!      - False - no tables 
!  * GOPosition  
!      - "Absolute" - vectors from cluster
!      - "Relative" - vectors relative to the basis atoms
!  * GOTbLattice: 
!      - None - Cluster reordered only into shells with respect to the basis atoms
!      - {{"A,B",{1,3}},{"A,C",{2}},...}:
!                Is used, if the information for the Hamiltonian construction is necessary. "A,B" defines
!                the interaction between two basis atoms. The list like {1,3} tells that the interaction should be taken from 
!                the first and third neighbor sphere of "A".
!                Usually we run GTShells with the standard options first. From the table the GOTbLattice input can be easily selected.               
! OUTPUT
!  list of shells
! ERROR MESSAGES
!  Error message amd Abort, if basis and shells do not have the same length.
!  Error message and Abort, if atom names in cluster and basis are in conflict.
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  The module becomes realtively involved due to handling f equivalent atoms in teh basis. The origin of the problems lies in the fact, that
!  for the construction of the TB Hamiltonians have to distinguished, but we do not want to blow up the parameter sets.
! LITERATURE
!  -
! TODO
!   may some additional checks will be necessary to avoid a crash if lists are empty
!   may be this can be done in a simpler way
! RELEASE
!  1.0.0
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTShells::badstruc =  "Lengths of basis list and shell list are different.";
GTShells::badnames =  "Names of atom in cluster and modified basis are different: `1` ";
GTShells[cl_, basis_, shells_, OptionsPattern[]] := 
    Module[{gvb, gol, gpos, cmd, nb, nat, shlist, prt, balist, center, type, shmax, distlist, prtemp, clt, tpos, k, l, tapp,  
            ns, divl, re, spw, tab, pos, txt, type1, type2, nsh, ts, ts1, ts2,basis1,tbas,mult,nat1,lbas,tb,sub,t,i,j,r,test,
            atl,mul,nt,nl,sym,sym1,nc,a,b,p,multa,multb,shl1,shl1t,s,s1,delta},
   (*--- Options ---*)
   delta=10^(-5);
   gvb = OptionValue[GOVerbose];
   gol = OptionValue[GOTbLattice];
   If[OptionValue[GOPosition] == "Absolute", gpos = False, gpos = True];
   (*--- messages ---*)
   cmd = "GTShells : ";
   (*--- handle equal atoms in the basis ---*)
   basis1 = basis;
   If[gol == {}, 
  	  None,
  (*--- new basis names ---*)	 
      tbas = Transpose[basis1];
      mult = Tally[tbas[[2]]]; nat = Length[mult];
      nat1 = Plus @@ Transpose[mult][[2]]; 
      If[nat1 == nat, 
       	lbas = False, 
     	lbas = True;tb = {}; sub = {}; 
     	Do[
           t  = Table[mult[[i, 1]] <> ToString[j], {j, 1, mult[[i, 2]]}]; 
           tb = Append[tb, t] // Flatten;
           r  = Table[t[[j]] -> mult[[i, 1]], {j, 1, mult[[i, 2]]}]; 
           sub = Append[sub, r]
        , {i, 1, nat}]; 
        sub = sub // Flatten; 
        tbas[[2]] = tb; basis1 = Transpose[tbas];
        Print[cmd, "modified basis : ", basis1];
        test = Complement[Transpose[cl][[2]] // Union, 
        Transpose[basis1][[2]] // Union]; 
        If[test == {}, 
          None, 
          Message[GTShells::badnames,test]; Abort[]
        ];
  (*--- new lattice transformation ---*)
     atl = Transpose[mult][[1]];
     mul = Transpose[mult][[2]];
     nt = Length[gol]; test = gol; nl = {}; 
     Do[
     	sym = test[[i, 1]]; nc = test[[i, 2, 1]];
        sym1 = StringSplit[sym, ","]; a = sym1[[1]]; b = sym1[[2]];
        If[a == b;
           p = Flatten[Position[atl, a]][[1]]; mult = mul[[p]];
           If[nc == 0,
              Do[
              	 nl = Append[nl, {StringJoin[a <> ToString[j] <> "," <> a <> ToString[j]], {0}}]
              , {j, 1, mult}],
           None
        ];
        If[nc > 0, 
           Do[
           	  Do[
           	  	 nl = Append[nl, {StringJoin[a <> ToString[i] <> "," <> a <> ToString[j]], {nc}}] 
           	  , {i, j + 1, mult}]
           , {j, 1, mult - 1}]],
           None
        ];
        If[a == b, 
           None,
           p = Flatten[Position[atl, a]][[1]]; multa = mul[[p]]; 
           p = Flatten[Position[atl, b]][[1]]; multb = mul[[p]];
           Do[
           	  Do[ 
                 nl = Append[nl, {StringJoin[a <> ToString[i] <> "," <> b <> ToString[j]], {nc}}] 
              , {i,1, multa}]
           , {j, 1, multb}]]
     , {i, 1, nt}]; 
  (*--- new substitution list and rules for back transformation are generated ---*)   
     gol = nl]
   ];
  (*--- consistency check ---*)
  nb = Length[basis1];
  nat = Length[cl];
  If[nb == Length[shells], 
  	 None, 
  	 Message[GTShells::badstruc]; Abort[]
  ];
  (*--- shell construction ---*)
  shlist = {}; prt = {};
  Do[
  	 balist = {};
     center = basis1[[i, 1]]; type = basis1[[i, 2]];
     shmax = shells[[i]] + 1;
  (*--- find distances with respect to an atom of the basis ---*)
     distlist = Take[Sort[Map[Norm[#] &, Map[center - # &, Transpose[cl][[1]]]] // Union, Less], {2, shmax}];
     prtemp = {};
     Do[
     	clt = {};
        Do[
           tpos = cl[[l, 1]] - center;
 (*          If[Norm[tpos] == distlist[[k]], *)
              If[Abs[Norm[tpos]-distlist[[k]]]< delta,
              If[gpos, 
           	     tapp = {tpos, cl[[l, 2]]}, 
           	     tapp = cl[[l]]
              ]; 
              clt = Append[clt, tapp], 
              None
            ]
        , {l, 1, nat}];
        balist = Append[balist, clt];
        prtemp = Append[prtemp, distlist[[k]]];
        prtemp = Append[prtemp, Transpose[clt][[2]] // Tally]
     , {k, 1, shells[[i]]}];
     shlist = Append[shlist, {type, balist}];
     prt = Append[prt, Flatten[{type, prtemp}, 1]]
  , {i, 1, nb}];
  (*--- Print table ---*)
  If[gvb === True, 
  	 ns = Max[shells]; divl = {{2 -> GTDividerColor1}};
     Do[
     	divl = Append[divl, {2 l + 1 -> Black, 2 l + 2 -> GTDividerColor2}]
     , {l, 1,ns}];
     divl = divl // Flatten; re = {"Basis"};
     Do[
     	re = Append[re, {"dist", "Atoms"}]
     , {i, 1, ns}];
     re = Flatten[re]; spw = Join[{re}, prt];
     tab = Grid[spw, Dividers -> {divl, {2 -> GTDividerColor1}}, Frame -> True, Background -> {{GTBackGroundColor1, None}, {GTBackGroundColor1, None},{1,1}->GTCornerColor}]; 
     Print[tab], 
     None
  ];
  (*--- info for Hamiltonian construction ---*)
  If[gol == {}, 
  	 Return[shlist], 
  	 ns = Length[gol]; shl1 = {}; prt = {};
     Do[
     	pos = {}; txt = StringSplit[gol[[i, 1]], ","]; 
     	type1 = txt[[1]];
        type2 = txt[[2]]; nsh = gol[[i, 2]];
        If[Length[nsh] == 1 && nsh[[1]] == 0, 
           shl1 = Append[shl1, {gol[[i, 1]], 0, {{{0, 0, 0}}}}], 
           ts = Select[shlist, #[[1]] == type1 &][[1, 2]];
           Do[
           	  ts1 = ts[[nsh[[k]]]];
              ts2 = Select[ts1, #[[2]] == type2 &];
              pos = Append[pos, Transpose[ts2][[1]]]
           , {k, 1, Length[nsh]}];
           shl1 = Append[shl1, {gol[[i, 1]], Length[nsh], pos}]
         ]
     , {i, 1, ns}];
   (*--- Replace names --*)
     If[lbas == True, 
     	shl1t = Transpose[shl1];
        s  = Map[StringSplit[#, ","] &, shl1t[[1]]] /. sub;
        s1 = Map[StringJoin[#[[1]] <> "," <> #[[2]]] &, s]; 
        shl1t[[1]] = s1; shl1 = Transpose[shl1t];
        Return[{shl1, Transpose[gol][[1]]}],
        Return[shl1]
     ];
       
  ]
]

(*
***) 		
	

(****h* /GTShellVectorsQlp
! NAME
!  GTShellVectorsQlp
! AUTHOR
!  W. Hergert
! PACKAGE
! Lattice.m 
! MODIFICATION HISTORY
!  * 13.08.2013 - first version
!  * 14.02.2015 - GTShellVectorQ renamed to GTShellVectorsQlp
!  * 01.02.2017 - output from GTShells becomes directly input of GTShellVectorQlp
!                 necessary data will be selected within the command.
!  * 24.05.2018 - check header and documentation
!
! USAGE
!  GTShellVectorsQlp[point group,shell vectors] calculates the minimal set of vectors SQ_l^p representing the coordination spheres.
!  All vectors of a shell can be created by a minimum number of vectors Q_l^p. All vectors
!  of the shell can be constructed from those vectors by means of the point group operations.
!  The module finds those set of vectors.
! INPUT
!  * pg     - Operations of the point group
!  * shell  - the cluster sorted into shells (output og GTShells)
!                
! OUTPUT
!  The vectors Q_l^p corresponding to the different shells.
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTGetMatrix
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  R.F. Egorov et. al. , phys. stat. sol. 26, 391 (1968), p. 396
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
	
	
 GTShellVectorsQlp[grp_,shell_]:=Module[
	{grp1,tmp,pos,sort,atom,qv,qvt,qvs,is,ns},
  (*--- take the right part from shells ---*)
	shell1 = shell[[1,2]];
	ns     = Length[shell1];
	qvt    = {};
	grp1   = GTGetMatrix[grp];
    Do[
       tmp = Transpose[shell1[[is]]];
       pos = tmp[[1]]//FullSimplify;
       sort= tmp[[2]];
       If[Length[Union[sort]]>1,
          Print["More than one sort in the shell!"],      
          Print["shell ",is," with ",Length[sort], " atoms of type ",sort[[1]]]
       ];
       qv   = pos[[1]];
       atom = sort[[1]];
       qvs  = {{qv,atom}};
       While[pos!={},
            pos = Complement[Union[grp1.qv]//FullSimplify,pos];
            If[pos!={},
               qv = pos[[1]];
               qvs=Append[qvs,{qv,atom}],
               None
            ]
       ];
       qvt = Append[qvt,qvs]
    ,{is,1,ns}];
    Return[qvt]
 ]
(*
***)    



(****h* /GTTransformToQlp
! NAME
!  GTTransformToQlp
! AUTHOR
!  W. Hergert
! PACKAGE
! Lattice.m 
! MODIFICATION HISTORY
!  * 13.08.2013 - 1st version 
!  * 14.02.2015 - GTTransformToQ renamed in GTTransformToQlp
!  * 01.02.2017 - output of GTShells will be directly used
!  * 24.05.2018 - check of header and documentation
!
! USAGE
!  GTTransformToQlp[point group,shell vectors,vectors qlp] gives the symmetry operations of a point group which transform the shell vectors to the vectors Q_l^p.
!  A set of vectors in a shell is represented by a vector Q_l^p. We find such pointgroup
!  operations, that transform the equivalent vectors of the shell to Q_l^p
!
! INPUT
!  * pg         - Operations of the point group
!  * shell      - List of the vectors in the different shells (output of GTShells)
!  * qvec       - set of vectors Q_l^p in the shells (output of GTShellVectorQ)               
! OUTPUT
!  List of symmetry transformations for the different Q_l^p and shells
! GTPAck OPTIONS
!  - 
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  * GTGetMatrix
!  * GTGetSymbol
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  -
! LITERATURE
!  R.F. Egorov et. al. , phys. stat. sol. 26, 391 (1968), p. 396
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

 GTTransformToQlp[pg_,shell1_,qvec_]:=Module[
	{nop,ns,opsets,nqlp,opsh,shv,nvec,opqv,qs,k,ind,op,m,l,shell},
     opsets={};
  (*--- take the coorect part from shells1 ---*)
     shell = shell1[[1,2]];
  (*--- nop - number of operations in point group, ns - shells ---*)
	 nop = Length[pg];
	 ns  = Length[qvec];
	 If[ns==Length[shell],
	    None,
	    Print["Error: List of shells and Subsuperscript[\[ScriptCapitalG], l, p] have to have the same length"];Abort[]
	 ];
  (*--- loop over shells ---*)
     Do[
     	nqlp = Length[qvec[[i]]];
     	opsh = {};
        shv  = Transpose[shell[[i]]][[1]];
        nvec = Length[shv];
        opqv = {};
        Do[
           qs=qvec[[i,k]][[1]];
           Do[
           	  ind=True;
              Do[
              	 op = GTGetMatrix[pg[[l]]];
                 If[op.shv[[m]]==qs&&ind==True,
            	    opqv = Append[opqv,GTGetSymbol[op]];
            	    ind  = False,
            	    None
                 ]
              ,{l,1,nop}]
           ,{m,1,nvec}];
           opsh = Append[opsh,opqv]
        ,{k,1,nqlp}];
        opsets = Append[opsets,opsh]
     ,{i,1,ns}];
     Return[opsets]
 ]

(*
***) 

(****h* /GTVoronoiCell
! NAME
!  GTVoronoiCell
! AUTHOR
!  W. Hergert
! PACKAGE
! Lattice.m 
! MODIFICATION HISTORY
!  * 23.03.2016 - first implementation
!  * 12.08.2016 - change of the design in 2D. Use lattice vectors instead of x,y (Matthias)
!  * 24.05.2018 - check header and documentation 
!  * 08.02.2023 : labelling of k-points changed according to new Mathematica version 
! USAGE
!  GTVoronoiCell[basvec,cldat,kpath] calculates the WS cell or the Brillouin zone
! INPUT
!  * basvec - basis vectors of the lattice, either for 2D or 3D case, for
!               direct or reciprocal lattice.
!  * cldat:
!
!           - cldat={cut,shmin,shmax},   
!           - cut: is the cut radius for the cluster construction,
!           - shmin,shmax: minimal and maximal shell number for the latice vector to consider
!  * kpath  - path in the BZ which can be shown if a Brillouin zone is considered
! OUTPUT
!WS cell or BZ
! GTPack OPTIONS
!  * GOBZPath
!  * GOOutPut
!  * GOVerbose
! STANDARD OPTIONS
!  * VertexLabelStyle - Directive[Black, 20,Background -> Yellow]  style for the symmetry point labels.
! GTPack MODULES
!  * GTGetMatrix
!  * GTGetSymbol
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  WS celles and BZs are based on the same construction principle, one in the direct and the other in the reciprocsal
!  lattice. 
! LITERATURE
!  -
! TODO
!   08.02.2023 : In GTFermiSurface and GTVoronoiCell Graphs are used for the k-paths. The representation of graphs
!               has change, thus problems appeared. A solution is found, but it is not optimal yet.
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
! 
*)

GTVoronoiCell[basvec_, cldat_, kpath_, OptionsPattern[]] :=
  Module[ {verb, bz, vec, dmax, lpt0, lpt1, test, dim,cor, path, ledge, x1, x2, t, tt, i, k, trple, ltr,
           p1, p2, p3, v1, v2, v3, n1, n2, n3, mat, lst, x, y, z, glg, res, nsol,modes,plt,clat,pts0,pts1,kpath1,
           sol,sol1,mesh,kos,gg,length,lstyle},
 (*--- Options ---*)
   modes = {"Voronoi", "Mesh", "Output"};
   verb  = OptionValue[GOVerbose];
   bz    = OptionValue[GOBZPath];
   plt   = OptionValue[GOOutput];
   lstyle= OptionValue[VertexLabelStyle];
   If[Intersection[{plt}, modes] == {},
      Print["Error: Mode ", plt, " not allowed!"]; Abort[],
      None
   ];
 (*--- Data for cluster construction ---*)
  {dmax, v1, v2} = cldat;
   clat = GTLatCluster[basvec, dmax, GOVerbose -> verb];
   vec  = GTLatShells[clat];
   pts0 = Flatten[Take[vec, {v1, v2}], 1];     lpt0 = Length[pts0];
   pts1 = Flatten[Take[vec, {v1, v2 + 2}], 1]; lpt1 = Length[pts1];
  (*--- Check of dimensionality ---*)  
  test = Transpose[pts0][[3]] // Union;
  If[Length[test] == 1,
     dim = 2;
     pts0 = Take[Transpose[pts0], {1, 2}] // Transpose;
     pts1 = Take[Transpose[pts1], {1, 2}] // Transpose,
     dim = 3;
     Print[Green, " x-axis ", Blue, " y-axis ", Red, " z-axis"];
  ];
  If[verb,
  	 Print["Dimension: ", dim];
  	 Print[lpt0, " vectors in search for Voronoi cell."];
     Print[lpt1, " vectors in check for correctness."];
     Print["Selected vectors from the cluster:"];
     Print[pts0],       
     None
  ];
  lst = Table[k, {k, 1, lpt0}]; trple = Subsets[lst, {dim}]; 
  ltr = Length[trple];
  If[verb,
     Print[ltr, " equations."], 
     None
  ];
  (*--- Calculate vertices for 2 or 3 dimensions ---*)
  sol = {}; sol1 = {};
  Switch[dim,
         2,
           Do[
              p1 = trple[[k, 1]]; p2 = trple[[k, 2]];
              v1 = pts0[[p1]]; v2 = pts0[[p2]];
              n1 = v1/Norm[v1]; n2 = v2/Norm[v2]; 
              mat = {n1, n2};
              If[Det[mat] == 0,
                 None,
                 glg = {n1.{x, y} == v1.n1/2, 
                 	    n2.{x, y} == v2.n2/2
                 	   };
                 res = Solve[glg, {x, y}];
                 nsol= {x, y} /. res // Flatten;
                 sol = Append[sol, nsol]
               ]
           , {k, 1, ltr}],
         3,
           Do[
              p1 = trple[[k, 1]]; p2 = trple[[k, 2]]; p3 = trple[[k, 3]];
              v1 = pts0[[p1]]; v2 = pts0[[p2]]; v3 = pts0[[p3]];
              n1 = v1/Norm[v1]; n2 = v2/Norm[v2];n3 = v3/Norm[v3];
              mat = {n1, n2, n3};
              If[Det[mat] == 0,
                 None,
                 glg = {n1.{x, y, z} == v1.n1/2, 
                 	    n2.{x, y, z} == v2.n2/2, 
                        n3.{x, y, z} == v3.n3/2
                       };
                 res = Solve[glg, {x, y, z}];
                 nsol= {x, y, z} /. res // Flatten;
                 sol = Append[sol, nsol]
              ]
           , {k, 1, ltr}],
          _,
            Print["Error: dimension not implemented!"]
  ];
 (*--- Check against other points 2D/3D ---*)
  sol = Union[sol];
  Do[
     t  = Norm[sol[[k]]];
     tt = Min[Map[Norm[sol[[k]] - #] &, pts1]];
     If[t <= tt,
        sol1 = Append[sol1, sol[[k]]],
        None
     ]
  , {k, 1, Length[sol]}];
  (*--- the Voronoi cell as a convex hull ---*)
  mesh = ConvexHullMesh[sol1];      
  (*--- Coordinate System ---*)
  Switch[dim,
         2,
           length=0.55*Max@Map[Norm, basvec];
           kos = Graphics[{{Gray, Thick, Arrow[{{0, 0}, basvec[[1]](*{length, 0}*)}]}, 
           	               {Gray, Thick, Arrow[{{0, 0}, basvec[[2]](*{0, length}*)}]}
           	              }],
         3, 
           kos = Graphics3D[{{Green, Arrowheads[.04], Arrow[Tube[{{0, 0, 0}, {1.75, 0, 0}}]]}, 
           	                 {Blue , Arrowheads[.04], Arrow[Tube[{{0, 0, 0}, {0, 1.75, 0}}]]}, 
           	                 {Red  , Arrowheads[.04], Arrow[Tube[{{0, 0, 0}, {0, 0, 1.75}}]]}
           	                }],
        _, 
          Print["Error: Arrows not implemented!"]
  ];
 (*--- k-path ---*)
 If[bz,
    kpath1=kpath;
    If[dim==2&&Length[kpath[[1,1]]]==3,
  	   kpath1[[1]] = Take[Transpose[kpath1[[1]]], {1, 2}] // Transpose,
  	   None
     ],
     None
  ];
  If[bz,
     cor = {}; path = {};
     ledge = Length[kpath[[1]]];
     Do[
        x1  = kpath1[[1, i]]; 
        x2  = kpath1[[2, i]];
        t   = x2 -> x1;
        cor = Append[cor, t];
        If[i < ledge,
           x1   = kpath1[[2, i]]; 
           x2   = kpath1[[2, i + 1]];
           t    = x1 -> x2;
           path = Append[path, t],
           None
        ]
     , {i, 1, ledge}];
     cor = Union[cor];
     If[dim == 2, 
        gg = GraphPlot[path, VertexCoordinateRules -> cor, VertexLabeling -> True, PlotStyle -> {Thick, Red},PlotTheme->"LargeGraph",VertexLabelStyle -> lstyle],
        gg = GraphPlot3D[path, VertexCoordinateRules -> cor, VertexLabeling -> True, PlotStyle -> {Thick, Red},PlotTheme->"LargeGraph",VertexLabelStyle -> lstyle]
     ],
     None
  ];
(*--- Display Voronoi cell ---*)
 If[plt === "Mesh" || plt === "Voronoi",
   If[plt === "Voronoi",
      If[bz,
      	 If[dim==2,
            Show[HighlightMesh[mesh, {Style[1,None], Style[2, Opacity[.5, RGBColor[70/256,111/256,242/256]]]}], kos, gg],
            Show[HighlightMesh[mesh, {Style[1, {Thick, Yellow}], Style[2, Opacity[.3, Pink]]}],kos,gg]
      	 ],
         If[dim==2,
            Show[HighlightMesh[mesh, {Style[1, None], Style[2, Opacity[.5, RGBColor[70/256,111/256,242/256]]]}], kos],
            Show[HighlightMesh[mesh, {Style[1, {Thick, Yellow}], Style[2, Opacity[.3, Pink]]}],kos]
      	 ]
      ],
       {ptt, lines} = {MeshCoordinates[mesh], MeshCells[mesh, 1]};
      If[bz,
         Show[MeshRegion[ptt, lines], kos, gg],
         Show[MeshRegion[ptt, lines], kos]
      ]   
   ],  
   If[bz,
      Return[{sol1, kos, gg}],
      Return[{sol1, kos}]
   ]      
 ]
]
  
 (*
 ***) 

(****h* /GTBZMPBPointMesh
! NAME
!  GTBZMPBPointMesh
! AUTHOR
!  W. Hergert
! PACKAGE
!  Lattice.m 
! MODIFICATION HISTORY
!  * 27.08.2016 - first version
!  * 24.05.2018 - check header and documentation 
!  * 03.10.2018 - recalculation of k-vectors implemented, new message system
! 
! USAGE
!  GTBZMeshExports[file, kmesh] exports a kmesh in MPB data format.
!
! INPUT
!   * file      - output file name
!   * kpoints   - list of kpoints
!  
! OUTPUT
!  data stored in file
! GTPack OPTIONS
!   GOVerbose
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  MPB_Symmetry_Analysis.nb
! DESCRIPTION
!   -
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)


GTBZMPBPointMesh::badbasis = "Error : Check your basis.";
GTBZMPBPointMesh::badvecs  = "Error : Basis vectors of different Length.";


GTBZMPBPointMesh[file_, kpts_, basis_, OptionsPattern[]] := Module[
	{nkp, data, i, c, vector, verb, file1, kpts1,chck,lchck,basis1,rec,kvec}, 
     verb = OptionValue[GOVerbose];
     If[Head[file] === String, 
     	file1 = file; 
     	kpts1 = kpts, 
     	file1 = kpts; 
        kpts1 = file
     ];
     nkp  = Length[kpts1];
  (*--- check basis ---*)
     chck  = Length[#] & /@ basis // Union; 
     lchck = Length[chck];
     If[lchck > 1,
         Message[GTBZMPBPointMesh::badvecs]; Return[],
         None
     ];
     chck = Length[#] & /@ basis;
     If[Length[basis] == 3 && Length[basis[[1]]] == 3,
        basis1 = basis,
        If[Length[basis] == 2 && Length[basis[[1]]] == 2,
           basis1 = {
           	{basis[[1, 1]], basis[[1, 2]], 0},
            {basis[[2, 1]], basis[[2, 2]], 0}, 
            {0, 0, 1}},
            Message[GTBZMPBPointMeshT::badbasis]; Return[]
       ]
     ];
   (*--- recalculate the k-points ---*) 
     rec  = Inverse[basis1] // Transpose; 
     data = {"(set! k-points (list"};
     Do[
     	kvec   = rec.kpts1[[i]];
        c      = Map[ToString[#] &, kvec // N];
        vector = 
        "(vector3 " <> c[[1]] <> " " <> c[[2]] <> " " <> c[[3]] <> ")";
        data   = Append[data, vector]
     , {i, 1, nkp}];
     data = Append[data, "))"];
     If[verb, 
        Print["Use ", file1, " in control file of MPB"];
        Print["Output in file    : ", file1];
        Print["Number of k-points: ", nkp]
     ];
     Export[file1, data, "Table"];
 ]


(* 
GTBZMPBPointMesh[file_, kpts_, OptionsPattern[{GOVerbose -> True}]] := Module[{nkp,data,i,c,vector,verb,file1,kpts1},
  verb=OptionValue[GOVerbose];
  If[Head[file]===String,
  	      file1=file;kpts1=kpts,
  	      file1=kpts;kpts1=file
  ];	      
  nkp = Length[kpts1];
  data = {"(set! k-points (list"};
  Do[
  	 c      = Map[ToString[#] &, kpts1[[i]] // N];
     vector = "(vector3 " <> c[[1]] <> " " <> c[[2]] <> " " <> c[[3]] <> ")";
     data   = Append[data, vector];
  , {i, 1, nkp}]; 
  data = Append[data , "))"];
  If[verb,
  	 Print["Use ",file1," in control file of MPB"];
  	 Print["Output in file    : ", file1];
  	 Print["Number of k-points: ", nkp]
  ];
  Export[file1, data, "Table"];
  ]
  
 *)
	 
	 
	 
	 
(*
***)	  
  
  
  
 (*
 ***) 
 (*-------------------------- Attributes ------------------------------*)
 Attributes[GTGroupOfK]={Protected, ReadProtected}
 (* Attributes[GTAdjacencyMatrix]={Protected, ReadProtected}*)
 Attributes[GTLatCluster]={Protected, ReadProtected}
 Attributes[GTCluster]={Protected, ReadProtected}
 Attributes[GTShells]={Protected, ReadProtected}
 Attributes[GTLatShells]={Protected, ReadProtected}
 Attributes[GTShellVectorsQlp]={Protected, ReadProtected}
 Attributes[GTGroupGlp]={Protected, ReadProtected}
 Attributes[GTTransformToQlp]={Protected, ReadProtected}
 Attributes[GTBZPointMesh]={Protected, ReadProtected}
 Attributes[GTBZPath]={Protected, ReadProtected}
 Attributes[GTBZLines]={Protected, ReadProtected}
End[] (* End Private Context *)



  
EndPackage[]
(*
***)
