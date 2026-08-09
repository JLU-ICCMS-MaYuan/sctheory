(****m* /Photonics.m
!
! NAME
!  Photonics.m
! AUTHOR
!  W. Hergert
! MODIFICATION HISTORY
!  * There was a previous old package, which is the basis of this development.
!  *    09.2015 : Check of documentation of Photonics.m
!  * 24.02.2023 : Error messages are completely in Mathematica style now. Old versions of 
!                 commands have been removed  
! USAGE
!  Contains all modules corresponding to the investigation of photonic band structures
!
! GTPack MODULES
!  
! --- Structure factor Fourier Transform of inverse permittivity ---
! 
!  * GTPhCuboid           - structure factor of a cuboid
!  * GTPhDCObjects        - Fourier transform of a inverse permittivity of an ensemble of objects
!  * GTPhDCPixel          - Fourier transform of a inverse permittivity from a pixel map
!  * GTPhDielectric       - Fourier transform of inverse permittivity
!  * GTPhEllipticRod      - structure factor of an elliptic rod (CircularRod as special case)
!  * GTPermittivity       - gives eps^(-1)(G) for a range of reciprocal lattice vectors
!  * GTPermittivityMatrix - gives a complete permitivity matrix eps^(-1)(G,G')
!  * GTPhPixelSmooth      - gives a smoothed pixel map
!  * GTPhPrismaticRod     - structure factor for a prismatic rod (RectangularRod as special case)
!  * GTPhRodSmooth        - smoothe structure factor of prismatic rod
!  * GTPhSlab             - structure factor for a slab 
!  * GTPhSlabSmooth       - smoothed structure factor of a slab
!  * GTPhSphere           - structure factor for a sphere
!
! --- Visualization of distribution of permittivity ---
!
!  * GTPhPixelStructure   - allows to inspec and modify the permittivity distribution given a s a pixel map
!  * GTPhShowStructure    - shows the arrangement of the odielectric obk
!
! --- Master equation in matrix form (plane wave basis) ---
!
!  * GTPhMaster           - constructs master equation for a slab
!  * GTPhMasterEquation   - constructs the "Hamiltonian" corresponding to the master equation
!  * GTPhMasterObjects    - constructs master equation from a list of objects
!  * GTPhMasterPixel      - constructs master eqaution from a pixel map  
!  * GTPhBandsObjects     - master equation with numerical evaluation of inverse eps for each k-point
!
! --- Symmetry analysis ---
!
!  * GTPhFields	         - calculation of electromagnetic fields 
!  * GTPhSymmetryBands   - symmetry analysis photonic band structure
!  * GTPhSymmeryField    - symmetry analysis of single electromagnetic field
!  * GTPhSymmetryPoint   - symmetry analysis at a single k-point
!  * GTPhUncoupledBands  - calculation of uncoupled bands
!  
! 
! --- Connection to MPB ---
!
!  * GTPhMPBBands        - reads photonic band structure from MPB
!  * GTPhMPBFields       - reads electromagnetic field calculated by means of MPB
!
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  Usually the calculation of photzonic band structures need a high number of plane waves in the expansion. Mathematica
!  will be too slow to do calculations with high accuracy. This part of GTPack can be used for principle understanding and
!  discussion of symmetry properties. For more accurate calculations the codes of the Joannopoulos group 
!
!	* MPB - photonic bands ( http://ab-initio.mit.edu/mpb/ )
!	* MEEP - FDTD code for electromagnetic systems ( http://ab-initio.mit.edu/meep/ )
!
!  should be used.
!
!  All modules start with GTPh#.
! LITERATURE
! 
!
***)

BeginPackage["GroupTheory`Photonics`",{"GroupTheory`Symbols`","GroupTheory`Lattice`","GroupTheory`Auxiliary`",
	                                   "GroupTheory`RepresentationTheory`","GroupTheory`Basic`","GroupTheory`ElectronicStructure`","GroupTheory`CrystalStructure`"}]
(*----------------------- Structure factors and Fourier transform of inverse permittivity -------*)
 GTPhCuboid               ::usage = "GTPhCuboid[\*StyleBox[\"reciprokal lattice vector,structure\",\"TI\"]] gives the structure factor for a cuboid."
 GTPhDCObjects            ::usage = "GTPhCuboid[\*StyleBox[\"reciprocal lattice vector,permittivity background, structure\",\"TI\"]] gives the Fourier transform of the inverse permittivity of an ensemble of geometric objects."
 GTPhDCPixel              ::usage = "GTPhDCPixel[\*StyleBox[\"reciprocal lattice vector,permittivity map, unit cell dimensions\",\"TI\"]] gives the Fourier transform of the inverse permittivity if the permittivity is defined by a pixel map."
 GTPhDielectric           ::usage = "GTPhDielectric[\*StyleBox[\"reciprocal lattice vector,permittivity structure, permittivity background, geometry\", \"TI\"]] gives the Fourier transform of \*Cell[BoxData[FormBox[RowBox[{SuperscriptBox[\"\[Epsilon]\",RowBox[{\"-\",\"1\"}]],\"(\",StyleBox[\"r\",FontWeight->\"Bold\"],\")\"}],TraditionalForm]],\"InlineMath\"]."
 GTPhEllipticRod          ::usage = "GTPhEllipticRod[\*StyleBox[\"reciprocal lattice vector,structure\",\"TI\"]] gives the structure factor for an elliptic rod."
 GTPhPermittivity         ::usage = "GTPhPermittivity[\*StyleBox[\"reciprocal lattice basis,cutoff,objects,permittivity background\",\"TI\"]] gives permittivity matrix or inverse permittivity matrix.."
 GTPhPermittivityMatrix   ::usage = "GTPhPermittivityMatrix[\*StyleBox[\"reciprocal lattice basis,cutoff,objects,permittivity background\",\"TI\"]] gives \[Epsilon]^-1 (Overscript[G, \[RightVector]],Overscript[G, \[RightVector]]')."
 GTPhPixelSmooth          ::usage = "GTPhPixelSmooth[\*StyleBox[\"pixel map\",\"TI\"]] gives a smoothed permittivity distribution defined by means of a pixel map."
 GTPhPrismaticRod         ::usage = "GTPhPrismaticRod[\*StyleBox[\"reciprocal lattice vector,structure\",\"TI\"]] gives the structure factor for a prismatic rod."
 GTPhRodSmooth            ::usage = "GTPhRodSmooth[\*StyleBox[\"reciprocal lattice vector,structure\",\"TI\"]] gives the structure factor for a smoothed rectangular rod."
 GTPhSlab                 ::usage = "GTPhSlab[\*StyleBox[\"reciprocal lattice vector,structure\",\"TI\"]] gives the structure factor for a slab."
 GTPhSlabSmooth           ::usage = "GTPhSlabSmooth[\*StyleBox[\"reciprocal lattice vector,structure\",\"TI\"]] gives the structure factor for a smoothed slab."
 GTPhSphere               ::usage = "GTPhSphere[\*StyleBox[\"reciprocal lattice vector,structure\",\"TI\"]] gives the structure factor for a sphere."

(*---------------------- Visualization of distribution of permittivity ---------------------------------*)
 GTPhPixelStructure       ::usage = "GTPhPixelStructure[\*StyleBox[\"permittivity map,modification,structure\",\"TI\"]] allows to modify the permittivity distribution given in pixel values."
 GTPhShowStructure        ::usage = "GTPhPixelStructure[\*StyleBox[\"structure,basis,number of unit cells\",\"TI\"]] generates an image of the arrangement of dielectric objects in the unit cell."

(* ---------------------- Master equation in matrix form (plane wave basis)  ------------------- *)
 GTPhMaster               ::usage = "GTPhMasterPixel[\*StyleBox[\"reciprocal lattice vector,permittivity map,unit cell dimensions\",\"TI\"]] constructs the master equation from precalculated Fourier transforms of 1/\[Epsilon]."
 GTPhMasterEquation       ::usage = "GTPhMasterEquation[\*StyleBox[\"structure,permittivity structure,permittivity background,geometry,reciproval lattice vectors\",\"TI\"]] constructs the matrix formulation of the master equation to calculate photonic band structures."
 GTPhMasterObjects        ::usage = "GTPhMasterObjects[\*StyleBox[\"list of objects,background permittivity,recirocal lattice vectors \",\"TI\"]] constructs the master equation from a list of geometric objects."    
 GTPhMasterPixel          ::usage = "GTPhMasterPixel[\*StyleBox[\"reciprocal lattice vector,permittivity map,unit cell dimensions\",\"TI\"]] constructs the master equation from a pixel map of the permittivity."    
 GTPhBandsObjects         ::usage = "GTPhBandsObjects[\*StyleBox[\"structure, permitivity background, reciprocal basis, cutoff, kpoints, number of bands\",\"TI\"]] calculates the photonic at bandstructure at kpoints for a structure defined by a list of objects."

(* ---------------------- Symmetry Analysis of fields ------------------------------------------ *)
 GTPhFields               ::usage = "GTPhFields[\*StyleBox[\"file name\", \"TI\"]\*StyleBox[\"k\",\"TI\",FontWeight->\"Bold\"]\*StyleBox[\"point,band number,\", \"TI\"]\*StyleBox[\"G\",\"TI\",FontWeight->\"Bold\"]\*StyleBox[\"vectors,geometry\", \"TI\"]] calculates the electro magnetic field in a photonic crystal."
 GTPhSymmetryBands        ::usage = "GTPhSymmetryBands[\*StyleBox[\"fileb,filef,klist,minb,maxb,objects,nmesh,group,basis\",\"TI\"]]  performs the symmetry analysis at a series of point s of a bandstructure. The results are presented in a table. Also a plot of the band structure with the assigned Ireps is possible."
 GTPhSymmetryField        ::usage = "GTPhSymmetryField[\*StyleBox[\"field,character table,number of points\",\"TI\"]] performs the symmetry analysis of field with respect to a point group with character table."
 GTPhSymmetryPoint        ::usage = "GTPhSymmetryPoint[\*StyleBox[\"filen,minb,maxb,objects,nmesh,group,basis\",\"TI\"]] performs the symmetry analysis at a certain k-point for band numbers min <= n <= maxb. objects defines the fields taken from the h5 files. nmesh sets the mesh size, group is the point group of the photonic crystal. basis defines the basis of the reciprocal lattice."
 GTPhUncoupledBands       ::usage = "GTPhUncoupledBands[\*StyleBox[\"reciprocal lattice vectors, band structure, matrix\",\"TI\"]] finds uncoupled bands from \*StyleBox[\"bands\", \"TI\"] if the bands are calculated with the reciprocal lattice vectors \*StyleBox[\"gvecs\",\"TI\"]. \*StyleBox[\"matrix\",\"TI\"] is a transformation matrix used to express the symmetry to find the uncoupled bands."
(*----------------------- Connection to MPB ----------------------------------------------------- *)
 GTPhMPBBands             ::usage = "GTPhMPBBands[\*StyleBox[\"band structure,symmetry points,basis\", \"TI\"]] reads a photonic band structure calculated by MPB for band structure plots or DOS calculations. \*StyleBox[\"basis\", \"TI\"] is the reciprocal basis used in the MPB calculation."
 GTPhMPBFields            ::usage = "GTPhMPBFields[\*StyleBox[\"file,objects\", \"TI\"]] reads a field calculated by MPB. "
 


(*----------------------- Internal Moduls ------------------------------------------------------- *)
(* GTToString              ::usage = "GTToString[n] number to String for file name construction"             *)

(*--------------------------- Options ----------------------------*)         
 Options[GTPhDCObjects]          = {GOVerbose -> False}
 Options[GTPhDCPixel]            = {GOVerbose -> False}
 Options[GTPhMaster]             = {GOPhPol -> "Automatic", GOVerbose -> True,GODCMethod->"MatrixInverse"}
 Options[GTPhMasterEquation]     = {GOPhPol->"TM"}
 Options[GTPhMasterObjects]      = {GOPhPol -> "Automatic", GOVerbose -> True}
 Options[GTPhMasterPixel]        = {GOPhPol -> "Automatic", GOVerbose -> True}
 Options[GTPhMPBBands]           = {GOVerbose -> True, GOPlotBands -> True}
 Options[GTPhMPBFields]          = {GOPlot -> False, GOTransformation->False,GOStore -> 0}
 Options[GTPhPermittivity]       = {GOVerbose -> False, GOPixel -> True}
 Options[GTPhPixelSmooth]        = {GOSmooth -> "5Points"}
 Options[GTPhPixelStructure]     = {GOPlot->True,GOColorScheme->"PearlColors"}  
 Options[GTPhShowStructure]      = {GOVerbose -> False,GOColorScheme->"PearlColors"}
 Options[GTPhSymmetryBands]      = {GOVerbose -> False, GOIrepNotation -> "Mulliken", PlotStyle -> {{Thin, Black}}, Joined -> True, PlotRange -> All,GOLabelShift->{0.05,0.02},
 	                                GOPlot->True,PlotLabel->"Photonic band structure",FrameLabel->{" ","\[Omega] a/( 2 \[Pi]c)"},GOPrecision->4,GOLabelStyle->{}}
 Options[GTPhSymmetryField]      = {GOPlot -> False, GOVerbose -> False, GOPlotStyle->"Re",ColorFunction->""}
 Options[GTPhSymmetryPoint]      = {GOVerbose -> False, GOIrepNotation -> "Mulliken", GOPlot->False, GOPlotStyle->"Re"}
 Options[GTPhBandsObjects]       = {GOVerbose -> True, GOPhPol -> "Automatic", GOPlotBands -> True, GOStore -> 0,GODCMethod->"Direct",GODecimals->4}
 Options[GTPhPermittivityMatrix] = {GOVerbose -> False, GOPixel -> False,GODCMethod->"Direct"}
 Options[GTPhUncoupledBands]     = {GOPlotBands -> True,GOVerbose->False}
    
Begin["`Private`"] (* Begin Private Context *) 


(****k* /GTPhSlab
! NAME
!  GTPhSlab
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 21.02.2016 : implementation
!  * 28.06.2018 : check header and documentation
! USAGE
!  structure factor of a dielectric slab of thickness a. The lengths are
!  mesured with respect to the lattice constant, i.e. the only important geometric 
!  quantity is the filling factor in this case.
! INPUT
!  gvec - wave vector
!
!  geom - structural data (geom[[1]] - filling factor)       
!
! OUTPUT
!  structure factor
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

GTPhSlab[gvec_, geom_] := Module[{arg, sfac}, 
  arg = geom[[1]]*\[Pi]*Norm[gvec];
  sfac = SphericalBesselJ[0, arg]; 
  Return[sfac]
]
  
(*
***)  


(****k* /GTPhSlabSmooth
! NAME
!  GTPhSlabSmooth
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 12.03.2016 : implementation
!  * 01.05.2016 : special cases g=0 and half width=0 implemented.
!  * 28.06.2018 : check header and documentation
! USAGE
!  structure factor of a dielectric slab of thickness a. 
!
!  The lengths are measured with respect to the lattice constant, i.e. the only important geometric 
!  quantity is the filling factor in this case.
!  
!  The profil is smoothed by convolution with a Gaussian
! INPUT
!  gvec - wave vector
!  
!  geom - structural data
 
!       o geom[[1]] - filling factor
!       o geom[[2]] -  half width of Guassian       
!
! OUTPUT
!  structure factor
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
! Please note, that in this case  analytically gamma S(G) was evaluated, i.e.
! we divide by gamma here,  to leave all other things in the GTPhDCObjects
! the same.
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! SO
URCE
!--------------------------------------------------------------------------------
!
*)

GTPhSlabSmooth[gvec_, geom_] := 
   Module[{s, gv, arg1, gam, sig, ef, ef1, ef2},
   sig = geom[[2]]; 
   gam = geom[[1]]; 
   If[sig==0,
   	  s=GTPhSlab[gvec,geom],
      arg1= (2 Sqrt[2] sig);
      gv  = Norm[gvec]; 
      If[gv==0,
      	 s=1.0,
         ef  = Erf[(-1 + gam)/arg1] + Erf[(1 + gam)arg1];
         ef  = ef*(-1 + Exp[2 I gv Pi])*Exp[gv Pi (I (-1 + gam) + 2 gv Pi sig^2)];
         ef1 = Erf[(-1 + gam - 4 I gv Pi sig^2)/arg1] - 
               Erf[(1 + gam - 4 I gv Pi sig^2)/arg1];
         ef2 =-Erf[(-1 + gam + 4 I gv Pi sig^2)/arg1] + 
               Erf[(1 + gam + 4 I gv Pi sig^2)/arg1];
         ef2 = ef2*Exp[2 I gv Pi gam];
         s   = ef + ef1 + ef2;
         s   = -s*Exp[-gv Pi (I gam + 2 gv Pi sig^2)]*I/(4 gv Pi); 
         s = s/gam
     ]
   ];
   Return[s]
]

(*
***)  

(****k* /GTPhRodSmooth
! NAME
!  GTPhRodSmooth
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 03.12.2016 : implementation
!  * 05.01.2016 : special case introduced, if half width =0
!  * 28.06.2018 : check header and documentation
! USAGE
!  structure factor of a dielectric rectangular rod with dimensions a,b. 
!  The lengths are measured with respect to the lattice constant.
!  The profile is smoothed by convolution with a Gaussian.
! INPUT
!  gvec - wave vector
!
!  geom - structural data 
!       geom[[1]]   - filling factor
!       geom[[2]]   - a
!       geom[[3]]]  - b
!       geom[[4]]   - half width of Gaussian       
!
! OUTPUT
!  structure factor
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
! Please note, that int this case  analytically gamma S(G) was evaluated, i.e.
! we divide by gamma here,  to leave all other things in the GTPhDCObjects
! the same
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

GTPhRodSmooth[gvec_, geom_] := 
 Module[{s, gam, sig, sx, sy, geo, arg,geom1,i},
  sig = geom[[4]];
  If[sig == 0,
  	 geom1 = Table[0, {5}];
     Do[
     	geom1[[i]] = geom[[i]]
     , {i, 1, 3}]; 
     geom1[[4]] = \[Pi]/2;
     s = GTPhPrismaticRod[gvec, geom1],
     arg = 1/(2 Sqrt[2] sig);
     (*--- x-direction ---*)
     gam = geom[[2]]; geo = {gam, sig};
     If[gvec[[1]] == 0,
        sx = 1/2 (Erf[(-1 + gam)*arg] + Erf[(1 + gam)*arg]);
        sx = sx + 1/2 gam (-Erf[(-1 +gam)*arg] +  Erf[(1 + gam)*arg]);
        sx = sx - Exp[-(1 + gam)^2/8/sig^2] (-1 + Exp[gam/(2 sig^2)]) Sqrt[2/Pi] sig,
        sx = GTPhSlabSmooth[gvec[[1]], geo]*gam
     ];
     (*--- y-direction ---*)
     gam = geom[[3]]; geo = {gam, sig};
     If[gvec[[2]] == 0,
        sy = 1/2 (Erf[(-1 + gam)*arg] + Erf[(1 + gam)*arg]);
        sy = sy + 1/2 gam (-Erf[(-1 +gam)*arg] +  Erf[(1 + gam)*arg]);
        sy = sy - Exp[-(1 + gam)^2/8/sig^2] (-1 + Exp[gam/(2 sig^2)]) Sqrt[2/Pi] sig,
        sy = GTPhSlabSmooth[gvec[[2]], geo]*gam
     ];
     s = sx*sy/geom[[1]]
  ];
  Return[s]
]


(*
***)  



(****k* /GTPhSphere
! NAME
!  GTPhSphere
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 21.02.2016 : implementation
!  * 28.06.2018 : check documentation and header
! USAGE
!  structure factor of a dielectric sphere of radius R. 
!  The radius is measured with respect to the lattice constant.
!
! INPUT
!  gvec - wave vector
!
!  geom - structural data  
!   o geom[[1]] - filling factor
!   o geom[[2]] - radius       
!
! OUTPUT
!  structure factor
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

GTPhSphere[gvec_, geom_] := Module[{arg, sfac}, 
  arg = 2* \[Pi] *Norm[gvec]*geom[[2]]; 
  If[arg<1.*10^(-6),
  	 sfac=1,
     sfac = 3 SphericalBesselJ[1, arg]/arg
  ]; 
  Return[sfac]
]
  
(*
***)  

(****k* /GTPhPrismaticRod
! NAME
!  GTPhPrismaticRod
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 21.02.2016 : implementation
!  * 28.06.2018 : check header and documentation
! USAGE
!  structure factor of a prismatic rod defined by the side lengths ax,ay
!  and the angle alpha.
! 
!  RectangularRod is used as a token in the package. This spacial case is
!  calculated with this modul.
!  
! Lengths is measured with respect to the lattice constant.
!
! INPUT
!  gvec - wave vector
!
!  geom - structural data 
!         geom[[1]] - filling factor, 
!         geom[[2]] - ax,
!         geom[[3]] - ay,
!         geom[[4]] - alpha.                 
!
! OUTPUT
!  structure factor 
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
GTPhPrismaticRod[gvec_, geom_] := Module[{a, b, alpha, sfac, gx, gy},
  a = geom[[2]]; 
  b = geom[[3]]; 
  alpha = geom[[4]]; 
  gx = gvec[[1]]; 
  gy = gvec[[2]]; 
  sfac = SphericalBesselJ[0, \[Pi]*gx*a]*SphericalBesselJ[0, b*\[Pi]*(gx Cos[alpha] + gy Sin[alpha])];
  Return[sfac]
]

(*
***)  


(****k* /GTPhEllipticRod
! NAME
!  GTPhEllipticRod
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 21.02.2016 : implementation
!  * 28.02.2018 : cheack header and documentation
! USAGE
!  structure factor of an elliptic rod defined by the major and minor half axis
!  CircularRod is used as a token in the package. This special case is
!  calculated with this modul.
!
!  The lengths are measured with respect to the lattice constant.
!
! INPUT
!  gvec - wave vector
!  geom - structural data 
!         geom[[1]] - filling factor, 
!         geom[[2]] - minor axis,
!         geom[[3]] - major axis.               
! OUTPUT
!  structure factor 
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

GTPhEllipticRod[gvec_, geom_] := Module[{minor, major, sfac, gy, gx, arg},
  minor = geom[[2]]; 
  major = geom[[3]]; 
  gx = gvec[[1]]; 
  gy = gvec[[2]]; 
  arg = 2 Pi minor*Sqrt[gx^2 + (major/minor) gy^2];
  If[Abs[arg] < 10^(-6), 
  	 sfac = 1, 
  	 sfac = 2 BesselJ[1, arg]/arg
  ]; 
  Return[sfac]
]


(*
***)  

(****k* /GTPhCuboid
! NAME
!  GTPhCuboid
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 21.02.2016 : implementation
!  * 28.06.2018 : check header and documentation
!  * 04.10.2018 : check header and documentation
! USAGE
!  Structure factor of a cuboid defined by the side lengths.
!  The lengths are measured with respect to the lattice constant.
!
! INPUT
!  gvec - wave vector
!  
!  geom - structural data 
!         geom[[1]] - side length x, 
!         geom[[2]] - side length y,
!         geom[[3]] - side length z.               
!
! OUTPUT
!  structure factor
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

GTPhCuboid[gvec_, geom_] := Module[{a, b, c, sfac, gx, gy, gz}, 
  a  = geom[[2]]; 
  b  = geom[[3]]; 
  c  = geom[[4]];
  gx = gvec[[1]]; 
  gy = gvec[[2]]; 
  gz = gvec[[3]];
  sfac = SphericalBesselJ[0, \[Pi] gx a] *
         SphericalBesselJ[0, \[Pi] gy b] *
         SphericalBesselJ[0, \[Pi] gz c]; 
  Return[sfac]
]

(*
***)
(****k* /GTPhDCObjects
! NAME
!  GTPhDCObjects
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 25.02.2016 : first version
!  * 12.03.2016 : implementation of smoothed structures
!  * 29.06.2018 : check header and documentation
!  * 04.10.2018 : check header and documentation
!  * 24.02.2023 : Mathematica Message system
!
! USAGE
!  GTPhDCObjects[file, objects] gives the Fourier transform of the inverse permittivity for an arrangement
!  of simple geometric forms in the unt cell. The Fourier transform is calculated analytically.
!  
!  The list struct describes  all objects: struct={object1,object2,object3,..}
!  
!  The following objects are defined: ElliptcRod, CircularRod, RectangularRod, PrismaticRod, Sphere, Cuboid, SlabSmooth, RodSmooth
!  
!  For each object the description has the form: object={identifier,permittivity,geometry,shift}
!
!   o shift is a vector to shift the object in the unit cell (in terms of the basis of the direct lattice)
!   o geometry  is a list that defines the size and shape of the object. The first entry ist the filling factor. 
!     The last one if meaningful, a rotation angle.
! INPUT
!   o gvec      - reciprocal lattice vector
!   o eback     - permittivity of the background
!   o struct    - definition of the structure
!
OUTPUT
!  Fourier transform of the arrangements of objects
! GTPack OPTIONS
!  o GOVerbose
!
!     - True   - additional information
!     - False  - no information
!
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTReadFromFile
! GTPack NOTEBOOKS 
!  GTPhMPBBands.nb 
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
! There was a problem with the definition of the structures, I have added a fith element in the geometry (rotation angle) some time ago, the full list is only necessary if the rotation will be really used 
!
!  29.06.2018 : There is still a problem with the rotations. 
! RELEASE
!  1.0.0
! PROBLEMS
!  04.10.2018 : May be that this has to be checked again. Seems that example on documentation page works correct.
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTPhDCObjects::name       = "`1` is not implemented.";

GTPhDCObjects [gvec_, eback_, struct_, OptionsPattern[]] := 
    Module[{verb,nst,info,eps,str,shift,kappam,ieps,sel,geom, i,sft,ftf, gv, part},
    (*---  Options ---*)
    verb = OptionValue[GOVerbose];
    (*--- analyse data ---*)
    nst   = Length[struct];
    info  = Transpose[struct][[1]];
    eps   = Transpose[struct][[2]];
    str   = Transpose[struct][[3]];
    shift = Transpose[struct][[4]];
    If[verb,
       Print["structures in elementary cell : " , info],
       None
    ];
    If[Norm[gvec]*1. < 10^(-6),
    (*--- mean value of kappa ---*)
       kappam = 1/eback;
       Do[
          kappam = kappam + str[[i, 1]] (1/eps[[i]] - 1/eback)
       , {i, 1, nst}];
       ieps = kappam;
       If[verb,
          Print["Mean value of 1/\[Epsilon]: ", kappam, " = ", kappam*1.0],
          None
       ],
    (*--- sum up the structures ---*)
       ieps=0;
       Do[
       	  sel  = info[[i]];
          geom = str[[i]];
          ftf  = geom[[1]]*(1/eps[[i]] - 1/eback);
          If[info[[i]]=="CircularRod"||info[[i]]=="EllipticRod"||
          	 info[[i]]=="PrismaticRod"||info[[i]]=="RectangularRod"||
          	 info[[i]]=="RodSmooth",
          	 If[Length[gvec] == 3, 
          	    gv = Take[gvec, {1, 2}], 
          	    gv = gvec
             ];
             If[Length[shift[[i]]] == 3, 
          	    sft = Take[shift[[i]], {1, 2}], 
                sft = shift[[i]]
             ];
             gv = RotationMatrix[geom[[5]]].gv,
             gv = gvec; sft = shift[[i]]
          ];      
          Switch[sel,
                "Slab",
                    part = ftf*Exp[-I 2 Pi gv[[1]] sft[[1]]] GTPhSlab[gv, geom];
                    ieps = ieps + part,
                "SlabSmooth",
                    part = ftf*Exp[-I 2 Pi gv[[1]] sft[[1]]] GTPhSlabSmooth[gv, geom];
                    ieps = ieps + part,
                "CircularRod",        
                    part = ftf*Exp[-I 2 Pi  gv.sft]*GTPhEllipticRod[gv, geom];
                    ieps = ieps + part,
                "EllipticRod",        
                    part = ftf*Exp[-I  2 Pi gv.sft]*GTPhEllipticRod[gv, geom];
                    ieps = ieps + part,
                "PrismaticRod",
                    part = ftf*Exp[-I 2 Pi gv.sft] GTPhPrismaticRod[gv, geom];
                    ieps = ieps+part,        
                "RectangularRod",
                    part = ftf*Exp[-I 2 Pi gv.sft] GTPhPrismaticRod[gv, geom];
                    ieps = ieps+part,           
                "RodSmooth",
                    part = ftf*Exp[-I 2 Pi gv.sft] GTPhRodSmooth[gv, geom];
                    ieps = ieps + part,
                "Sphere",
                    part = ftf*Exp[-I 2 Pi gvec.shift[[i]]] GTPhSphere[gvec, geom];
                    ieps = ieps + part,
                "Cuboid",
                    part = ftf*Exp[-I 2 Pi gvec.shift[[i]]] GTPhCuboid[gvec, geom];
                    ieps = ieps + part,
                 _, 
                    Message[GTPhDCObjects::name, info[[i]]]; Abort[]
          ]
       , {i, 1, nst}]
     ];
     Return[Chop[ieps]]
 ]
 
(*
***)
  
(****k* /GTPhPixelStructure
! NAME
!  GTPhPixelstructure
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!   * 08.02.2016 : first version
!   * 29.06.2018 : check header and documentation
! USAGE
!  GTPhPixelstructure allows to change the permittivity values of some pixels, if the permittivity is defined
!  in Form of a pixel map.
! INPUT
!   o permittivity map -  list of the permittivity values in the pixels
!   o modification     - list contains the changes {{i1,j1,eps1},i2,j2,eps2},...}
!   o geom             - defines the structure of the unit cell {a,b,phi}                
!
! OUTPUT
!  The permittivity distribution is plotted, or the list containing the distribution is in the output.
! GTPack OPTIONS
!   o ColorScheme  - defines the color scheme for the plot
!   o GOPlot
!
!       True  - Plot of the modified structure
!       False - output of permittivity distribution
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

GTPhPixelStructure[map_, modify_, geom_, OptionsPattern[]] := 
      Module[{a, b, \[Phi], mapi, maxc, i, j, nx, ny, nch, lx, ly, plt, points, pg,g,pixel,col}, 
      plt = OptionValue[GOPlot];
      col = OptionValue[GOColorScheme];
      mapi = map;
      (*--- define rhombus ---*)
      {a, b, \[Phi]} = geom;
      points = {{0, 0}, {a, 0}, {a + b Cos[\[Phi]], 
      b Sin[\[Phi]]}, {b Cos[\[Phi]], b Sin[\[Phi]]}}; 
      pg = Polygon[points];
      (*--- modification ---*)
      nch = Length[modify];
      If[nch > 0,
         Do[
         	lx = modify[[i, 1]]; 
         	ly = modify[[i, 2]];
            mapi[[lx, ly]] = modify[[i, 3]],
         {i, 1, nch}],
         None
      ]; 
      maxc=Max[mapi];
      g = {}; 
      ny = Length[mapi];
      Do[
      	 nx = Length[mapi[[j]]]; 
         Do[
         	pixel = {EdgeForm[Black], ColorData[col][(mapi[[i, j]] - 1)/(maxc - 1)], 
                    Translate[pg, {(i - 1)*a + (j - 1)*b*
                    Cos[\[Phi]], (j - 1) b Sin[\[Phi]]}]};
            g = Append[g, pixel]
         , {i, 1, nx}]
      , {j, 1, ny}];
      (*--- plot the new matrix ---*)
      If[plt,
         Grid[{{Graphics[g, 
              ImageSize -> 500]}, {BarLegend[{col, {1, maxc}}, 
              LegendLabel -> "Permittivity \[Epsilon](\!\(\*StyleBox[\"r\",\nFontWeight->\"Plain\"]\))", 
              LegendLayout -> "Row", 
              LegendMarkerSize -> 200]
              }}],
              Return[mapi]
       ]
]

(*
***)

(****k* /GTPhShowStructure
! NAME
!  GTPhShowStructure
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 21.02.2016 : first version
!  * 23.05.2016 : correct representation of the background color implemented.
!  * 29.06.2018 : check header and documentation 
!  * 24.02.2023 : Mathematica message system
!
! USAGE
!  GTPhShowStructure shows the arrangement of the geometrical objects defining the 
!  distribution of the permittivity in the unit cell
! INPUT
!  structure - list difinig the structures  
!  basis     - basis vectors of th direct lattice
!  nez       - number of unit cells (nez=0 - 1 unit cell)               
!
OUTPUT
!  Plot of the distribution of the permittivity.
! GTPack OPTIONS
!   o GOColorScheme  - Color scheme for the plot
!   o GOVerbose
!        
!        True    - additional information 
!        False   - no additional information (standard)!
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
GTPhShowStructure::name       = "`1` is not implemented.";

GTPhShowStructure[struct_, eback_, basis_, nez_, OptionsPattern[]] := 
   Module[{verb,bas,b1,b2,nst,info,eps,maxc,str,shift,g,go,sft,i,a,b,alpha,v1,v2,gt,gl,uz,vec,j,col,got},
  (*---  Options ---*)
  verb = OptionValue[GOVerbose];
  col = OptionValue[GOColorScheme];
  (*--- analyse data ---*)  
  bas = Take[Transpose[basis], {1, 2}] // Transpose;
  b1 = bas[[1]]; b2 = bas[[2]];
  nst    = Length[struct];
  info   = Transpose[struct][[1]];
  eps    = Transpose[struct][[2]];
  maxc   = Max[{eback,eps}];
  str    = Transpose[struct][[3]];
  shift = Transpose[struct][[4]];
  If[verb,
     Print["structures in elementary cell : " , info];
     Print["maximum value of permittivity : ",maxc],
     None
  ];
  (*--- Construction of the graphical objects ---*)
  g = {};
  Do[
     sft = Take[shift[[i]], {1, 2}];
     sft = sft[[1]]*b1 + sft[[2]]*b2;       
     Switch[info[[i]],
        "CircularRod",
           go = Disk[sft, str[[i, 2]]],
        "EllipticRod",
           got = Disk[sft, {str[[i, 2]], str[[i, 3]]}];
           If[Length[str[[i]]] == 5, 
           	  go = Rotate[got, str[[i, 5]]], 
           	  go=got
           ],
        "RectangularRod",    
           a = str[[i, 2]];
           b = str[[i, 3]];
           v1= {a , b}/2;
           v2= {-a, b }/2;
           got = Translate[Polygon[{v1, v2, -v1, -v2}], sft];
           If[Length[str[[i]]] == 5, 
           	  go = Rotate[got, str[[i, 5]]], 
           	  go = got
           ],
        "PrismaticRod",
           a = str[[i, 2]]; b = str[[i, 3]]; alpha = str[[i, 4]];
           v1 = {a + b Cos[alpha], b Sin[alpha]}/2; 
           v2 = {b Cos[alpha] - a, b Sin[alpha]}/2;
           got = Translate[Polygon[{v1, v2, -v1, -v2}], sft];
           If[Length[str[[i]]] == 5, 
           	  go = Rotate[got, str[[i, 5]]], 
         	  go = got
           ],
        _,
           Message[GTPhShowStructure::name, info[[i]]]; Abort[]           
     ];
     g = Append[g, {ColorData[col][(eps[[i]]-1)/(maxc - 1)], go}]
  , {i, 1, nst}];
  gt = {}; gl = {};
  (*--- construction of the unit cells *)
  uz = {EdgeForm[Directive[Thick, Black]], ColorData[col][(eback - 1)/(maxc - 1)], 
        Polygon[{{0, 0}, b1, b1 + b2, b2}]
       };
  Do[
  	 Do[
        vec = i*b1 + j*b2;
        gt = Append[gt, Translate[g, vec]];
        gl = Append[gl, Translate[uz, vec]]
    , {i, 0, nez}]
  , {j, 0, nez}];
  Grid[{{Graphics[{gl,gt}, 
         ImageSize -> 500]}, {BarLegend[{col, {1, maxc}}, 
         LegendLabel -> "Permittivity \[Epsilon](r)", LegendLayout -> "Row", 
         LegendMarkerSize -> 200]
       }}
  ]
]
 
 
 
(*
***)


(****k* /GTPhDCPixel
! NAME
!  GTPhDCPixel
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 03.03.2016 : first version
!  * 28.06.2018 : check header and documentation
!  * 04.10.2018 : check header and documentation
! USAGE
!  GTPhDCPixel gives the Fourier transform of the inverse permittivity if the permittivity
!  is defined by a pixel map.
!
! INPUT
!   o gvec     - reciprocal lattice vector
!   o pixelmap - list of the permittivity values in the pixels
!   o basis    - dimension of the unit cell (a,b,alpha - side lengths and enclosed angle)             
!
! OUTPUT
!  eps^(-1)(G)
! GTPack OPTIONS
!  GOVerbose
!    o True  - additional information 
!    o False - no sadditional information (standard)
! STANDARD OPTIONS
!   -
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

GTPhDCPixel[gvec_, pixelmap_, basis_, OptionsPattern[]] :=
    Module[{verb, npixa, npixb, a, b, pa, pb, geom, kappam, ieps, sfac, i, j,tau,s},
    (*--- Options ---*)
    verb = OptionValue[GOVerbose];
    (*--- Analysis of the data ---*)
    npixa = Length[pixelmap[[1]]]; 
    npixb = Length[pixelmap];
    a = Norm[basis[[1]]]; b = Norm[basis[[2]]];
    pa = basis[[1]]/a; pb = basis[[2]]/b; \[Alpha] = ArcCos[pa.pb];
    pa = a*pa/npixa; pb = pb*b/npixb;
    geom = {0, N[a/npixa], N[b/npixb], \[Alpha]};
    If[verb,
       Print[npixa, " divisions of basis vector ", basis[[1]]];
       Print[npixb, " divisions of basis vector ", basis[[2]]];
       Print["pixel data ", geom],
       None
     ];
     (*--- mean value of kappa ---*)
     kappam = N[Plus @@ Flatten[1/pixelmap]/npixa/npixb];
     If[verb,
        Print["Mean value of 1/\[Epsilon]: ", kappam, " = ", kappam],
        None
     ];
     If[Norm[gvec]*1. < 10^(-6),
       ieps = kappam,
       (*--- structure independent factor ---*)
       sfac = GTPhPrismaticRod[gvec, geom];
       If[verb,
          Print["structure independent factor ", sfac],
          None
       ];
       (*--- evaluate structure factor ---*)
       s = 0;
       Do[
          Do[
             tau = (i - 1)*pa + (j - 1)*pb + (pa + pb)/2;
             s = s + Exp[-I N[2 Pi gvec.tau]]/pixelmap[[i, j]];
          , {i, 1, npixa}]
       , {j, 1, npixb}];
       ieps = s*sfac/npixa/npixb
     ];
     Return[ieps // Chop]
 ]

(*
***) 

(****k* /GTPhMasterPixel
! NAME
!  GTPMhasterCPixel
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 03.03.2016 : first version
!  * 29.06.2018 : check header and documentation
! USAGE
!  GTPhMasterPixel constructs the master equation from a pixel map of the permittivity.  
!  is defined by a pixel map.
!
! INPUT
!  o gvec     - reciprocal lattice vector
!  o pixelmap - list of the permittivity values in the pixels
!  o basis    - dimension of the unit cell (a,b,alpha - side lengths and enclosed angle)             
! 
! OUTPUT
!   photonic Hamiltonian
! GTPack OPTIONS
!   o GOPhPol
!
!     - "Automatic"  - no special polarization (1D or 3D)
!     - "E" or "TM"  - TM polarization
!     - "H" or "TE"  - TE polarization
!
!   o GOVerbose
!
!     - True  - additional information (standard)
!     - False - no additional information
! GTPack MODULES
!  GTPhDCPixel
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
! SOURCE
!--------------------------------------------------------------------------------
!
*)  

GTPhMasterPixel[gvec_, pixelmap_, basis_, OptionsPattern[]] := 
  Module[{k, poltypes,opol,verb,ndim,ham,gi,gj,i,j,gv,fteps},
  (*--- Definitions ---*)
  k = xez;
  poltypes = {"E", "H", "TM", "TE", "Automatic"};
  (*--- Options ---*)
  opol = OptionValue[GOPhPol];
  verb = OptionValue[GOVerbose];
  (*--- check polarization ---*)
  If[Position[poltypes, opol] == {},
     Print["Error : Polarization ", opol, " wrong!"]; Abort[],
     None
  ];
  If[opol == "E",
     opol = "TM";
     If[verb,
        Print["Polarization E corresponds to TM"],
        None
     ],
     None
  ];
  If[opol == "H",
     opol = "TE";
     If[verb,
        Print["Polarization H corresponds to TE"],
        None
     ],
     None
  ];
  ndim = Length[gvec];
  ham = Table[0, {ndim}, {ndim}];
  Do[
  	 gi = gvec[[i]];
     Do[
     	gj = gvec[[j]]; gv = gi - gj;
        fteps = GTPhDCPixel[gv, pixelmap, basis];
        If[opol == "TM",
           ham[[i, j]] = fteps*Norm[k + gj]^2,
           None
        ];
        If[opol == "TE",
           ham[[i, j]] = fteps*(k + gi).(k + gj),
           None
        ]
     , {j, 1, ndim}]
  , {i, 1, ndim}];
  Return[ham]
]



(*
***)

(****k* /GTPhPermittivity
! NAME
!  GTPhPermittivity
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 03.03.2016 : first version
!  * 28.06.2018 : check header and documentation
!  * 24.02.2023 : Mathematica messages
! USAGE
!  GTPhPermittivity  gives eps^(-1)(G) for all potential reciprocal lattice vectors in the master equation
!
! INPUT
!  basis    - basis vectors reciprocal lattice
!
!  ncut     - value of |G| for cut of reciprocal lattice vectors
!
!  objects  
!      - "Pixel"   : {pixelmap,unitcell}
!      - "Objects" : list of objects         
!
!  eback    - background value of permittivity
!
! OUTPUT
!   eps^(-1)(G)
! GTPack OPTIONS
!  GOVerbose
!    o True  - additional information 
!    o False - no sadditional information (standard)
!  
!  GOPixel
!    o True   - definition as a pixel map (standard)
!    o False - definition by an analytical expression
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

GTPhPermittivity::option = "GOPixel has wrong value."

GTPhPermittivity[basis_, ncut_, objects_, eback_, OptionsPattern[]] := 
    Module[{verb,meth,gvec,gd,gv,i,j,lg,ld,pmap,basvec,ft},
    (*--- Options ---*) 	
     verb = OptionValue[GOVerbose];
     meth = OptionValue[GOPixel];
     If[meth === True || meth === False,
        None,
        Message[GTPhPermittivity::option]; Abort[]
     ];
     (*--- reciprocal lattice vectors ---*)
     gvec = GTLatCluster[basis, ncut, GOVerbose -> verb];
     lg=Length[gvec];gd = {};
     Do[
        Do[
           gv = gvec[[i]] - gvec[[j]];
           If[Position[gd, gv] == {}, 
           	  gd = Append[gd, gv], 
           	  None
           ]
         , {i, 1, lg}]
     , {j, 1, lg}]; ld = Length[gd];
     If[verb,
        Print["Number of differences of recip. lattice vectors: ", ld],
        None
     ];
     If[meth,
          pmap   = objects[[1]]; 
          basvec = objects[[2]];
          If[Length[basvec[[1]]]==2,
          	 basvec=Append[Transpose[basvec], {0, 0}] // Transpose,
          	 None
          ];
          ft = Map[GTPhDCPixel[#, pmap, basvec] &, gd]
       ,
          ft = Map[GTPhDCObjects[#, eback, objects] &, gd]
     ];
     Return[{ft, gd}]
 
 ]
(*
***) 



(****k* /GTPhPermittivityMatrix
! NAME
!  GTPhPermittivityMatrix
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 24.11.2016 : first version
!  * 28.06.2018 : check header and documentation
!  * 24.02.2023 : Mathematica messages
! USAGE
!  GTPhPermittivityMatrix calculates the the Inverse of the  Fourier transforms of the permittivity,i.e.
!  eps^(-1)(G,G')
!
! INPUT
!  basis    - basis vectors reciprocal lattice
!
!  ncut     - value of |G| for cut of reciprocal lattice vectors
!
!  objects
!    - "Pixel"   : {pixelmap,unitcell}
!    - "Objects" : list of objects         
! 
!  eback    - background value of permittivity
!
! OUTPUT
!  eps^(-1)(G,G')
! GTPack OPTIONS
!  GOVerbose
!    o True  - additional information 
!    o False - no sadditional information (standard)
!  
!  GOPixel
!    o True   - definition as a pixel map 
!    o False - definition by an analytical expression (standard)
!
!  GODCMethod
!    o "Direct"         - direct calculation
!    o "MatrixInverse"  - Ho's Method   
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  In the usual method the Fourier transform of the inverse permittivity is calculated. In Ho's method
!  the matrix eps(G,G') is calculated first and then inverted. The elements of the matrix are used in the setup
!  of the Hamiltonian. This improves the convergence.
!  If Ho's method is used, the inverse of the permittivity values have to bee calulated first to get the correct
!  values in GTPhDCObjects.
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

GTPhPermittivityMatrix::option = "GODCMethod describes wrong method."
GTPhPermittivityMatrix::pixel  = "GOPixel has wrong value."


GTPhPermittivityMatrix[basis_, ncut_, objects_, eback_, OptionsPattern[]] := Module[
	{verb,meth,gvec,lg,gd,gd1,i,j,pmap,basvec,ebacki,st,ft,mho},
  (*---Options---*)
    verb = OptionValue[GOVerbose];
    meth = OptionValue[GOPixel];
    mho  = OptionValue[GODCMethod];
    If[mho == "MatrixInverse" || mho == "Direct",
      None,
      Message[GTPhPermittivityMatrix::option]; Abort[]
    ];
    If[meth === True || meth === False,
       None,
       Message[GTPhPermittivityMatrix::pixel]; Abort[]
    ];
  (*---reciprocal lattice vectors---*)
    gvec = GTLatCluster[basis, ncut, GOVerbose -> verb];
    lg   = Length[gvec]; 
    gd   = Table[gvec[[i]] - gvec[[j]], {i,1,lg}, {j,1,lg}];
    gd1  = Flatten[gd, 1];
  (*--- calculate matrix, change structure data accordinglyy ---*)  
  (*--- pixel map ---*)
    If[meth,
       pmap   = objects[[1]]; 
       If[mho=="MatrixInverse",
       	  pmap=1./pmap,
       	  None
       ];
       basvec = objects[[2]];
       If[Length[basvec[[1]]] == 2, 
          basvec = Append[Transpose[basvec], {0, 0}] // Transpose, 
          None
       ];
       ft = Map[GTPhDCPixel[#, pmap, basvec] &, gd1], 
   (*--- predefined objects ---*)  
       If[mho=="MatrixInverse",
       ebacki= 1./eback  ;
       st=objects//Transpose;st[[2]]=1./st[[2]];st=st//Transpose,
       ebacki=eback;st=objects
       ];
       ft = Map[GTPhDCObjects[#, ebacki, st] &, gd1]
    ]; 
    gd = Partition[ft, lg];
    If[mho=="MatrixInverse",
       gd = Inverse[gd // N],
       gd=gd//N
    ];
    Return[{gvec, gd}]
]

(*
***) 

(****k* /GTPhMaster
! NAME
!  GTPhMaster
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 15.02.2016 : first version
!  * 26.11.2016 : reformulation for Ho's method
!  * 24.02.2023 : Mathematica messages , old version deleted
! USAGE
!  GTPhMaster constructs the master equation from precalculated Fourier transforms of the inverse permittivity.
!
! INPUT
!  * eback    - background permittivity
!  * gvec     - reciprocal lattice vectors
!  * ftperm   - {gdif, ftt} 
!               gdif is the list of reciprocal lattice vectors that appear in the construction of the master equation
!               ftt contains the corresponding Fourier transfroms 
!
! OUTPUT
!  master equation - photonic "Hamiltonian"
! GTPack OPTIONS
!  o GOVerbose
!
!     - True  - additional information 
!     -False - no sadditional information (standard)
! 
!  o GOPhPol
!
!     - "Automatic"  - no special polarization (1D or 3D)
!     - "E" or "TM"  - TM polarization
!     - "H" or "TE"  - TE polarization
! 
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
GTPhMaster::pol    = "Wrong polarisation for two-dimensional case."
GTPhMaster::method = "Method is not implemented."

GTPhMaster[eback_, gvec_, epsmat_, OptionsPattern[]] := Module[
	{k,c,poltypes,opol,verb,dim,test,ndim,ham,gi,gj,fteps,i,j,gv1,ngv1,pv1,npv,pv,pol2,
	 ki,kj,p1,p2},
  (*---Definitions---*)
    k = xez; 
    c = {1, 1.7, 2.051}; c = c/Norm[c];
    poltypes = {"E", "H", "TM", "TE", "Automatic"};
 (*---Options---*)
    opol = OptionValue[GOPhPol];
    verb = OptionValue[GOVerbose];
 (*---check polarization---*)
    If[Position[poltypes, opol] == {}, 
      Print["Error : Polarization ", opol, " wrong!"]; Abort[], 
      None
    ];
    If[opol == "E", 
       opol = "TM";
       If[verb, 
       	  Print["Polarization E corresponds to TM"], 
       	  None
       ], 
       None
    ];
    If[opol == "H", 
       opol = "TE";
       If[verb, 
       	  Print["Polarization H corresponds to TE"], 
   	      None
   	   ], 
       None
    ];
 (*---check dimension of the problem---*)
    dim = 3;
    test = Map[Union[#] &, Transpose[gvec]];
    If[test[[3]] == {0}, 
       dim = dim - 1, 
       None
    ];
    If[test[[2]] == {0}, 
       dim = dim - 1, 
       None
    ];
    If[verb, 
       Print["Dimensionality: ", dim], 
       None
    ];
    ndim = Length[gvec];
    Switch[dim,
   2, 
     If[(dim == 2 && opol == "TM") || (dim == 2 && opol == "TE"), 
        None, 
        Message[GTPhMaster::pol]; Abort[]
     ];
     ham = Table[0, {ndim}, {ndim}];
     Do[
     	gi = gvec[[i]];
        Do[
           gj = gvec[[j]];
           fteps = epsmat[[i, j]];
           If[opol == "TM", 
           	ham[[i, j]] = fteps*Norm[k + gj]^2, 
           	None
           ];
           If[opol == "TE", 
           	 ham[[i, j]] = fteps*(k + gi).(k + gj), 
           	 None
           ]
       , {j, 1, ndim}]
     , {i, 1, ndim}],
   3, 
     gv1  = Map[k + # &, gvec];
     ngv1 = Map[Norm[#] &, gv1];
     pv1  = Map[Cross[c, #] &, gv1];
     npv  = Map[Norm[#] &, pv1];
     pv1  = pv1/npv; 
     pv   = {};
     Do[
     	pol2 = Cross[pv1[[i]], gv1[[i]]];
        npv  = Norm[pol2];
        pv   = Append[pv, {pv1[[i]], pol2/npv}]
     , {i, 1, ndim}];
     ham = Table[0, {2*ndim}, {2*ndim}];
     Do[
        Do[
           Do[
              Do[
                 fteps = epsmat[[i, j]];
                 ki    = i + (p1 - 1)*ndim;
                 kj    = j + (p2 - 1)*ndim;
                 ham[[ki, kj]] = 
                   fteps*ngv1[[i]]*ngv1[[j]]*pv[[i, p1]].pv[[j, p2]]*(-1)^(p1 + p2)
              , {j, 1, ndim}]
           , {i, 1,ndim}]
        , {p1, 1, 2}]
     , {p2, 1, 2}], 
   _, 
     Message[GTPhMaster::method]; Abort[]
   ];
   Return[ham]
]



(*
***)


(****k* /GTPhPMasterObjects
! NAME
!  GTPhMasterObjects
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 14.03.2016 : first version
!  * 09.04.2016 : only half of the matrices are calculated for 2D and 3D problems, use of internal cross product
!  * 29.06.2018 : Check header and documentation
! USAGE
!  GTPhMasterObjects constructs the master equation from analytical expressions of the Fourier coefficients of the inverse
!  permittivity
!
! INPUT
!   o eback    - background permittivity
!   o gvec     - reciprocal lattice vectors
!   o struct   - list of structural elements in the unit cell: {object1, object2,..}
!                object={identifier,permittivity,geometry,shift}
!
!             - shift is a vector to shift the object in the unit cell (interms of the basis of the direct lattice)
!             - geometry  is a list that defines the size and shape of the object. The first entry ist the filling factor. 
!               The last one if meaningful, a rotation angle.
!             
!
! OUTPUT
!  master equation - photonic "Hamiltonian"
! GTPack OPTIONS
!   o GOPhPol
!
!     - "Automatic"  - no special polarization (1D or 3D)
!     - "E" or "TM"  - TM polarization
!     - "H" or "TE"  - TE polarization
!
!   o GOVerbose
!
!     - True  - additional information (standard)
!     - False - no additional information
! STANDARD OPTIONS
!   -
! GTPack MODULES
!   -
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

GTPhMasterObjects::struc    = "Structure `1`is not implemented."
GTPhMasterObjects::object   = "Objects `1`  not implemented."
GTPhMasterObjects::pol      = "Polarization `1`  is wrong."
GTPhMasterObjects::pol2D    = "Wrong polarization for two-dimensional case."
GTPhMasterObjects::problem  = "Problem not defined."
 
GTPhMasterObjects[struct_, eback_, gvec_, OptionsPattern[]] := 
  Module[{k,poltypes,types,opol,verb,ndim,ham,gkv,i,nst,info,test,dim,gi,gj,j,gv,fteps,gv1,ngv1,pv,pv1,npv,pol2,
  	      ki,kj,p1,p2},
  (*--- Definitions ---*)
   k = xez;
   poltypes = {"E", "H", "TM", "TE", "Automatic"};
   types = {"CircularRod", "EllipticRod", "RectangularRod", "PrismaticRod", "Sphere", "Cuboid", 
   	        "Slab", "EmptyLattice","SlabSmooth","RodSmooth"};
  (*--- Options ---*)
   opol = OptionValue[GOPhPol];
   verb = OptionValue[GOVerbose];
  (*--- empty lattice ---*)
   ndim = Length[gvec];
   ham  = Table[0, {ndim}, {ndim}];
   If[Head[struct] === String,
     If[struct == "EmptyLattice",
       If[verb,
          Print["Empty lattice band structure"],     
          None
       ];
       gkv = Map[Sqrt[(xez + #).(xez+#)] &, gvec];
       Do[
          ham[[i, i]] = gkv[[i]]
       , {i, 1, ndim}];
       Return[ham],
       Message[GTPhMasterObjects::struc, struct]; Abort[]
    ],
    None
   ];
  (*--- all other cases ---*)
   If[Head[struct] == List,
      nst  = Length[struct];
      info = Transpose[struct][[1]];
      If[verb,
        Print["structures in elementary cell: " , info],
        None
   ];
  (*--- check type of objects ---*)
   test = Complement[info, types];
   If[test == {},
      None,
      Message[GTPhMasterObjects::object, test];
      Abort[]
   ];
  (*--- check dimension of the problem ---*)
   dim = 3;
   test = Map[Union[#] &, Transpose[gvec]];
   If[test[[3]] == {0},
      dim = dim - 1,
      None
   ];
   If[test[[2]] == {0},
      dim = dim - 1,
      None 
   ];
   If[verb,
      Print["Dimensionality: ", dim], 
      None
   ];
  (*--- check polarization ---*)
   If[Position[poltypes, opol] == {},
   	Message[GTPhMasterObjects::pol, opol]; Abort[],
      None
   ];
   If[opol == "E",
      opol = "TM";
      If[verb,
         Print["Polarization E corresponds to TM"],
         None
      ],
      None
   ];
   If[opol == "H",
      opol = "TE";
      If[verb,
         Print["Polarization H corresponds to TE"],
         None
      ],
      None
   ];
  (*--- Setup Hamiltonian ---*)
   Switch[dim,
  (*--- 1D problem ---*)
       1,
         Do[
         	gi=gvec[[i]];
         	Do[
         	   gj=gvec[[j]];gv = gi - gj;
               fteps = GTPhDCObjects[gv, eback, struct, GOVerbose -> False];        
               ham[[i,j]]=fteps*(k+gi).(k+gj);             
            ,{j,1,ndim}]
         ,{i,1,ndim}],    	
  (*--- 2D-problem ---*)
       2,
         If[(dim == 2 && opol == "TM") || (dim == 2 && opol == "TE"),
   	        None,
            Message[GTPhMasterObjects::pol2D];Abort[]
         ];
         Do[
            gi = gvec[[i]];
            Do[
               gj = gvec[[j]]; gv = gi - gj;
               fteps = GTPhDCObjects[gv, eback, struct, GOVerbose -> False];
               If[opol == "TM",
                  ham[[i, j]] = fteps*Norm[k + gj]*Norm[k + gi],
                  None
               ];
               If[opol == "TE",
                  ham[[i, j]] = fteps*(k + gi).(k + gj),
                  None
               ];
            , {j, 1, ndim}]
         , {i, 1, ndim}],
   (*--- 3D-Problem ---*)       
       3,
         gv1  = Map[k + # &, gvec];
         ngv1 = Map[Norm[#] &, gv1];
         pv1  = Map[Cross[{0, 0, 1}, #] &, gv1];
         npv  = Map[Norm[#] &, pv1];
         pv1  = pv1/npv; pv = {};
         Do[
         	pol2 = Cross[pv1[[i]], gv1[[i]]];
            npv = Norm[pol2];
            pv = Append[pv, {pv1[[i]], pol2/npv}]
         , {i, 1, ndim}];
         ham = Table[0, {2*ndim}, {2*ndim}];
         Do[
            Do[
               Do[
                  Do[
                  	 gv = gvec[[i]] - gvec[[j]];
                     fteps = GTPhDCObjects[gv, eback, struct, GOVerbose -> False];         
                     ki = i + (p1 - 1)*ndim;
                     kj = j + (p2 - 1)*ndim;
                     ham[[ki, kj]] = fteps*ngv1[[i]]*ngv1[[j]]*pv[[i, p1]].pv[[j, p2]]*(-1)^(p1 + p2);
                  , {j, 1, ndim}]
               , {i, 1, ndim}]
           , {p1, 1, 2}]
        , {p2, 1, 2}],
       _,
         Message[GTPhMasterObjects::problem]
     ]
    ];
    Return[ham//N]
 ]
     
                  
(*
***)          

(****k* /GTPhPixelSmooth
! NAME
!  GTPhPixelSmooth
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 25.02.2016 : first version
!  * 28.06.2018 : check header and documentation
!  * 24.02.2023 : Mathematica messages
! USAGE
!  GTPhPixelSmooth
! INPUT
!  permittivity map -  list of the permittivity values in the pixels
!
!  modification     - list contains the changes {{i1,j1,eps1},i2,j2,eps2},...}
!
!  geom             - defines the structure of the unit cell {a,b,phi}                
!
! OUTPUT
!  Usually the permittivity distribution is plotted.
! GTPack OPTIONS
!
! STANDARD OPTIONS
!  
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
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTPhPixelSmooth::option       = "Method `1` is notimplemented.";

GTPhPixelSmooth[map_, OptionsPattern[]] := Module[{sm,nx,ny,ms1,i,j,ms},
  (*--- Options ---*)
  sm = OptionValue[GOSmooth];
  If[sm=="5Points",
  	 None,
  	 Message[GTPhPixelSmooth::option, sm];Abort[]
  ];	  
  (*--- Prepare data ---*)
  nx = Length[map];
  ny = Length[map[[1]]];
  ms1 = Table[0, {nx}, {ny}]; 
  ms = Table[0, {nx + 2}, {ny + 2}];
  (*--- Prepare table for cyclic interpolation ---*) 
  Do[
     Do[
        ms[[i + 1, j + 1]] = map[[i, j]]
     , {i, 1, nx}]
  , {j, 1, ny}];
  Do[
  	 ms[[i + 1, 1]] = map[[i, ny]];
     ms[[i + 1, ny + 2]] = map[[i, 1]]
  , {i, 1, nx}];
  Do[
  	 ms[[1, i + 1]] = map[[nx, i]];
     ms[[nx + 2, i + 1]] = map[[1, i]]
  , {i, 1, ny}];
  Switch[sm,
    "5Points",          
         Do[
         	Do[
         	   ms1[[i - 1, j - 1]] = 
         	   (ms[[i, j]] + ms[[i - 1, j]] + ms[[i + 1, j]] + ms[[i, j - 1]] + ms[[i, j + 1]])/5
            ,{i, 2, nx + 1}]
         , {j, 2, ny + 1}],
     _, 
         Message[GTPhPixelSmooth::option, sm]
  ];
  Return[ms1]
]

(*
***)

(****k* /GTPhDielectric
! NAME
!  GTPhDielectric
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!   * 28.06.2018 : check header and documentation
!   * 04.10.2018 : check header and documentation
!   * 24.02.2023 : Mathematica messages
! USAGE
!  GTPhDielectric[reciprocal lattice vector,structure,permittivity structure,permittivity background,geometry] gives the
!  Fourier transform of \[Epsilon]^-1(r).
! INPUT
!  gabs         -  reciprocal lattice vector or norm of the reciprocal lattice vector
!
!  type         -  type of photonic "atoms"
!
!                  RectangularRod - dielectric rectangular rods in air, or alternatively holes
!                  CircularRod    - dielectric circular rod in air, or alternatively holes  
!                  EllipticalRod  - elliptical rod
!                  Sphere
!                  Cuboid      
!
!  estruc       -  permittivity of the structure
!
!  eback        -  permittivity of the background
!
!  geometry     -  list containing paramters of the structure 
!                  Rod : {radius of cylinder}
!
! OUTPUT
!  Fourier transform of inverse  epsilon to the given reciprocal lattice vector
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
GTPhDielectric::struc = "Structure `1` is not implemented.";

GTPhDielectric[gab_,type_,estruc_,eback_,geometry_]:=Module[
	{ftf,ff,gabs,types,rad,argx,argy,argz,ae,be,arg,phi,gx,gy,gf,tab,geom,types1},
	       types={"CircularRod","RectangularRod","EllipticalRod","Sphere","Cuboid","One-dimensional"};
           types1={"CircularRod    ",
           	       "RectangularRod ",
           	       "EllipticalRod  ",
           	       "Sphere         ",
           	       "Cuboid         ",
           	       "One-dimensional"};
           geom={"filling factor, radius",
           	     "filling factor, ax, ay",
           	     "filling factor, a, b, \[Phi]",
           	     "filling factor, radius",
           	     "filling factor, ax,ay,az",
           	     "filling factor"
           	     };
           If[type=="Help",
           	  tab = {types1, geom} // Transpose; 
           	  tab = Prepend[tab, {"Structure", "Geometry"}]; 
           	  Print[Grid[tab, Alignment -> Left, Frame -> All, Background -> {{1 -> LightBlue}, 
           	  	    {1 -> GTBackGroundColor1}, {1, 1} ->GTBackGroundColor1}]]; Return[],
           	  None
           ];	  
           gabs=Norm[gab];
           If[type=="CircularRod",           	 	  
           	  ff=geometry[[1]];	    	   
           	  rad=geometry[[2]];           	  
              If[gabs==0,
              	ftf=ff/estruc+(1-ff)/eback,
                arg=2*Pi*gabs*rad;
              	ftf=2*ff*(1/estruc-1/eback)*BesselJ[1,arg]/arg
              ],
              None
            ];       
            If[type=="RectangularRod",         
               ff=geometry[[1]];
               argx=gab[[1]]*geometry[[2]];argy=gab[[2]]*geometry[[3]]; 
               If[gabs==0,
              	  ftf=ff/estruc+(ff-1)/eback,
              	  ftf=ff*(1/estruc-1/eback)*Sin[argx]*Sin[argy]/argx/argy
               ],
               None
            ];  
            If[type=="EllipticalRod",         
               ff=geometry[[1]];
               ae=geometry[[2]];be=geometry[[3]];phi=geometry[[4]];   
               gx=gab[[1]];gy=gab[[2]];
               gf=gx^2(Cos[phi]^2+(be/ae)^2 Sin[phi]^2)+gy^2(Sin[phi]^2+(be/ae)Cos[phi]^2) 
                                        +2 gx gy Cos[phi] Sin[phi](1-(be/ae)); 
               If[gabs==0,
              	  ftf=ff/estruc+(ff-1)/eback,
              	  arg=ae*Sqrt[gf];
              	  ftf=2*ff*(1/estruc-1/eback)*BesselJ[1,arg]/arg
               ],
               None
            ];     
            If[type=="Sphere",	  
           	  ff=geometry[[1]];	
           	  rad=geometry[[2]];           	  
              If[gabs==0,
              	 ftf=ff/estruc+(1-ff)/eback,
                 arg=2*Pi*gabs*rad;
              	 ftf=3*ff*(1/estruc-1/eback)*SphericalBesselJ[1,arg]/arg
               ]
              ,None
            ];          
            If[type=="Cuboid",	  
           	  ff=geometry[[1]];	
           	  argx=gab[[1]]*geometry[[2]];argy=gab[[2]]=geometry[[3]]; 
           	  argz=gab[[3]]*geometry[[3]];
              If[gabs==0,
              	 ftf=ff/estruc+(1-ff)/eback,
              	 ftf=ff*(1/estruc-1/eback)*Sin[argx]*Sin[argy]+Sin[argy]/argx/argy/argz
                 ]
              ,None
           ];  
           If[type == "One-dimensional",
             ff = geometry[[1]];  
             arg = ff*\[Pi]*gab[[1]];
             If[gabs == 0,
                ftf = ff/estruc + (1 - ff)/eback,
                ftf = ff*(1/estruc - 1/eback)*Sin[arg]/arg
             ],
           None
           ];
           If[Complement[{type},types]=={},
              Return[ftf],
              Message[GTPhDielectric::struc, type]
     	   ]
]


(*
***) 


(****k* /GTPhMasterEquation
! NAME
!  GTPhMasterEquation
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 11.11.2013 : first version2013
!  * 03.03.2015 : case "Empty lattice" added
!  * 11.01.2016 : one-dimensional case added
!  * 29.06.2018 : check header and documentation
!  * 24.02.2023 : Mathematica messages, old version removed
! USAGE
!  GTPhMasterEquation[structure,permittivity structure, permittivity background,geometry,reciproval lattice vectors] constructs the matrix formulation 
!   of the master equation to calculate photonic band structures.
! INPUT
!  * type         -  type of photonic "atoms" (see GTPhDielectric)                  
!  * estruc       -  permittivity of the structure
!  * eback        -  permittivity of the background
!  * geo          -  list containing paramters of the structure 
!  * gvec         -  reciprocal lattice vectors             
!
OUTPUT
!  Fourier transform of epsilon to the given reciprocal lattice vector
! GTPack OPTIONS
!   o GOPhPol
!
! GTPack MODULES
!   GTPhDielectric, 
!   PhCrossP (internal only)
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!   The formulation     
! LITERATURE
!   K. Sakoda, Optical Properties of Photonic Crystals
!   J.D. Joannopoulos, Photonic Crystals - Molding the flow of light
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  o Problems with som k-vectors /0 appears -> one component 0.0001
!    nort elegant, but it works
!  o it is very slow
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTPhMasterEquation::struc = "Type `1` is not implemented."

GTPhMasterEquation[type_, estruc_, eback_, geo_, gvec_, OptionsPattern[]] :=
     Module[{k, types, ps, pol, ndim, ham, gkv, i,j, gi, gj, gv, fteps,gv1,ngv1,pv1,npv,pv,
     	     pol2,p1,p2,ki,kj},
     k = xez; ndim = Length[gvec];
    (*--- implemented structures ---*)
     types = {"CircularRod", "EllipticalRod", "RectangularRod", "Sphere", "Cuboid", "EmptyLattice","One-dimensional"};
     If[Position[types, type] == {},
     	Message[GTPhMasterEquation::struc,type]; Abort[],      
        ps = Flatten[Position[types, type]][[1]]
     ];
     (*--- Polarisation ---*)
     pol = OptionValue[GOPhPol];
     If[pol == "E",
        pol = "TM",
        None
     ];
     If[pol == "H",
        pol = "TE",
       None
     ];
     (*--- Empty Lattice ---*)
     If[type === "EmptyLattice",
        ham = Table[0, {ndim}, {ndim}];
        gkv = Map[(Norm[xez + #]^2) &, gvec];
        Do[
           ham[[i, i]] = gkv[[i]]
        , {i, 1, ndim}],
        None
     ];
     (*--- 1D case ---*)
     If[type === "One-dimensional", 
     	ham = Table[0, {ndim}, {ndim}];
        Do[gi = gvec[[i]];
           Do[gj = gvec[[j]]; gv = gi - gj;
              fteps = GTPhDielectric[gv, type, estruc, eback, geo];
              ham[[i, j]] = fteps*(k + gi).(k + gj)
           , {j, 1, ndim}]
        , {i, 1, ndim}],
        None
     ];
    (*--- 2D case ---*)
    If[ps <= 3,
       ham = Table[0, {ndim}, {ndim}];
       Do[gi = gvec[[i]];
          Do[gj = gvec[[j]]; gv = gi - gj;
             fteps = GTPhDielectric[gv, type, estruc, eback, geo];
             If[pol == "TM",
                ham[[i, j]] = fteps*Norm[k + gj]^2,
                None
             ];
             If[pol == "TE",
                ham[[i, j]] = fteps*(k + gi).(k + gj),
                None
             ]
          , {j, 1, ndim}]
       , {i, 1, ndim}],
       None
    ];
    (*--- 3D case ---*)  
    If[ps == 4 || ps == 5, 
       gv1 = Map[k + # &, gvec];
       ngv1 = Map[Norm[#] &, gv1];
       pv1 = Map[PhCrossP[{0, 0, 1}, #] &, gv1];
       npv = Map[Norm[#] &, pv1];
       pv1 = pv1/npv; pv = {};
       Do[
       	  pol2 = PhCrossP[pv1[[i]], gv1[[i]]];
          npv = Norm[pol2];
          pv = Append[pv, {pv1[[i]], pol2/npv}]
       , {i, 1, ndim}];
       ham = Table[0, {2*ndim}, {2*ndim}];
       Do[
   	      Do[
   	      	 Do[
   	      	 	Do[
   	      	 	   gv = gvec[[i]] - gvec[[j]];
                   fteps = GTPhDielectric[gv, type, estruc, eback, geo];
                   ki = i + (p1 - 1)*ndim;
                   kj = j + (p2 - 1)*ndim;
                   ham[[ki, kj]] = fteps*ngv1[[i]]*ngv1[[j]]*
                                   pv[[i, p1]].pv[[j, p2]]*(-1)^(p1 + p2)
                , {j, 1, ndim}]
             , {i, 1,ndim}]
          , {p1, 1, 2}]
       , {p2, 1, 2}],
       None
    ];
    Return[ham]
]





(*
***) 


(****k* /PhCrossP
! NAME
!  PhCrossP
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  1st version November 2013
! USAGE
!  Crossproduct of two vectors
! INPUT
!  vectors a,b
!
OUTPUT
!  c= a x b
! ERROR MESSAGES
!  
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  
! DESCRIPTION
!  There was a problem  to use the crossproduct of the VectorAnalysis package.
!  If this proble is resolved this Module is obsolete. 
!
! LITERATURE
!  
! TODO
!  check problem with VectorAnalysis package
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

PhCrossP[a_, b_] := Module[{ax, ay, az, bx, by, bz, cp},
              ax = a[[1]]; ay = a[[2]]; az = a[[3]];
              bx = b[[1]]; by = b[[2]]; bz = b[[3]];
              cp = {ay*bz - az*by, az*bx - ax*bz, ax*by - ay*bx};
              Return[cp]
]
		
(*
***) 


(****k* /GTPhFields
! NAME
!  GTPhFields
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 11.11.2013 : first version
!  * 29.06.2018 : check header and documentation
! USAGE
!  GTPhFields[file name,kpoint,band number,Gvectors,geometry] calculates the electromagnetic field in a photonic crystal
!  from the eigenvectors as solution of the master equation
! INPUT
!   o file  - Name of file which contains the eigenvectors of the problem at a certain k-point 
!   o kpt   - the corresponding k-point
!   o bdnr  - band number
!   o gvec  - reciprocal lattice vectors
!   o geom  - geometry information of the system
!
! OUTPUT
!  values of the field 
! GTPack Modules
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTReadFromFile
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
! SOURCE
!--------------------------------------------------------------------------------
!
*)


GTPhFields[file_,kpt_,bdnr_,gvec_,geo_]:=Module[{wav,nkpt,ng,i,cvecs,nx,ny,xmax,xmin,dx,ymax,ymin,dy,field,ka,fld,rvec,j,s,k },
                         wav=GTReadFromFile[file];nkpt=Length[wav];ng=Length[gvec];
                         Do[
                            If[wav[[i,1]]==kpt,
                               cvecs=wav[[i,2]],
                               None
                              ] 
                         ,{i,1,nkpt}];
                         nx=geo[[3]];ny=geo[[6]];
                         xmax=geo[[1]];xmin=geo[[2]];dx=(xmax-xmin)/(nx-1);
                         ymax=geo[[4]];ymin=geo[[5]];dy=(ymax-ymin)/(ny-1);   
                         field={};ka=kpt;
                         Do[fld={};
                            Do[
                               rvec={xmin+(i-1)*dx,ymin+(j-1)*dy,0}; 
                               s=0;         
                               Do[
                               	  s=s+cvecs[[bdnr,k]]*Exp[I 2 \[Pi] (ka+gvec[[k]]).rvec]
                               ,{k,1,ng}];
                               fld=Append[fld,s]
                             ,{i,1,nx}];
                             field=Append[field,fld]
                         ,{j,1,ny}];
                         Return[field]
]


(*
***) 


(****k* /GTPhMPBBands
! NAME
!  GTPhMPBBands
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 17.02.2016 : first version
!  * 09.05.2016 : now also possible to read data for DOS calculation
!  * 28.06.2018 : check header and documentation
!  * 24.08.2018 : basis as argument included for correct handling of nonrectangular lattices
!  * 04.10.2018 : problem with detection of symmetry points resolved.
!  * 24.02.2023 : Mathematica messages, old version removed
!
! USAGE
!  GTPhMPBBands[file, symmetry points,basis] is used to read a photonic band structure from file  calculated by MPB. 
!  The band structure is calculated along lines defined by the symmetry points. Also frequencies for DOS Calculations can
!  be read. k-vectors are recalculated in Cartesian basis.
!
!  The scripts GSMPBBands2D/GSMPBBands3D is used to construct a data set that can be used by the GTBandsPlot command. 
!  Output of the script are the files mbnds_te and mbnds_tm containing the data in the appropriate form.
! INPUT
!   file             - output file from GSMPBBands2D /GSMPBBands3D

!   symmetry points  - symmetry points used  by MPB {{point_1,name_1},{point_2,name_2}...} 
!
!   basis            - reciprocal basis used in MPB calculation
!
!  
! OUTPUT
!  a list that can be used in GTBandsPlot or in DOS calculations
! GTPack OPTIONS
!  GOVerbose
!     o True  - additional information (standard)
!     o False - no additional information
!
!  GOPlotBands
!     o True  - prepare for band structure plot (True)
!     o False - prepare for DOS calculation 
! GTPack MODULES
!  GTReadFromFile
! GTPack NOTEBOOKS 
!  GTPhMPBBands.nb 
! DESCRIPTION
!  -
! LITERATURE
!  - 
! TODO
!  -
! PROBLEMS
!  24.08.2018
!  Up to now all worked fine, because we have tested only a square lattice. In this case the coordinates 
!  related to the reciprocal basis are the same like the coordinates in cartesian basis. For hexagonal lattice as 
!  an example MPB gives the k-vectors in reciprocal basis. All commands in GTPack work in Cartesian basis. This causes
!  problems if one compares results with GTPack calculations. Therefore a recalculation to Cartesian basis is included. 
!  Therefore the basis used in the MPB calculation has to be provided. If MPB bands are compared with GTPack bands the x-axis
!  is now automatically the same.
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTPhMPBBands::badfile =  "Error: File `1` not in this directory.";

GTPhMPBBands[file_, points_,basis_, OptionsPattern[]] := Module[
	{verb, fn, isl, bands, i, nbands, nkpt, absz, gpt, bnds, pt1, pt2, j, kpt, del,nplt,g}, 
  (*--- options ---*)	
	verb = OptionValue[GOVerbose];
    nplt = OptionValue[GOPlotBands];
  (*--- check file names ---*)  
    fn   = FileNames[]; 
    isl  = Intersection[fn, {file}];
    If[isl == {}, 
       Message[GTPhMPBBands::badfile,file]; Return[], 
       bands = GTReadFromFile[file]
    ];
    isl = {};
    Do[
       If[IntegerQ[bands[[i]]], 
       	  isl = Append[isl, bands[[i]]], 
          None
       ]
    , {i, 1, Length[bands]}];
    nbands = Length[bands]/Max[isl] - 5;
    g      = Partition[bands, Length[bands]/Max[isl]];
    nkpt   = Length[g];
    If[verb, 
       Print["Number of bands   : ", nbands];
       Print["Number of k-points: ", nkpt], 
       None
    ];
    If[nplt,
 (*--- prepare bands ---*)
    absz = 0; gpt = {};
    bnds = {{g[[1, 1]], 0, Take[g[[1]], {2, 4}], Take[g[[1]], {6, 5 + nbands}]}};
    Do[
       pt1  = Take[g[[i - 1]], {2, 4}]; 
       pt2  = Take[g[[i]], {2, 4}];
  (*--- recalculate k-vectors: reciprocal basis ->Cartesian basis ---*)
       If[Length[pt1]==3&&Length[basis[[1]]]==2,
       	  pt1=Take[pt1,{1,2}];
       	  pt2=Take[pt2,{1,2}],
       	 None
       ];
       pt1=Total[pt1*basis];
       pt2=Total[pt2*basis];
       absz = absz + Norm[pt2 - pt1];
       kpt  = {g[[i, 1]], absz, Take[g[[i]], {2, 4}], Take[g[[i]], {6, 5 + nbands}]};
       bnds = Append[bnds, kpt]
    , {i, 2, nkpt}];
    If[verb, 
       Print["Maximum abszissa  : ", absz], 
       None
    ];
(*    If[IntegerQ[points[[1, 1]]], 
       gpt = points, *)
       gpt={};
       Do[
       	  Do[
       	  	 del = Norm[Take[g[[i]], {2, 4}] - points[[j, 1]]];
             If[del < 10^(-4), 
             	gpt = Append[gpt, {i, points[[j, 2]]}], 
                None
             ]
         , {j, 1, Length[points]}]
      , {i, 1, nkpt}];
      (* ];  *)
      Return[{bnds, gpt}],
 (*--- prepare DOS ---*)
      bnds = {};
      Do[
      	 kpt = {Take[g[[i]], {2, 4}], Take[g[[i]], {6, 5 + nbands}]};
         bnds = Append[bnds, kpt]
      ,{i, 1, nkpt}];
      Return[bnds]
    ]
  ]
  
	
(*
***) 


(****k* /GTPhMPBFields
! NAME
!  GTPhMPBBands
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 24.02.2016 : first version
!  * 07.03.2016 : improvents with respect to field analysis
!  * 18.10.2016 : inplementation of 3D case
!  * 02.07.2018 : check header and documentation
!  * 24.02.2023 : Mathematica messages
! USAGE
!  GTPhMPBFields[file, objects] is used to read a field calculated by means of MPB.
!  
!  The fields have to be treated by means of mpb-data befor the symmetry analysis.
!  The module imports directly h5 files 

! INPUT
!   file  - h5 files of MPB after treatment with mpb-data
!   
! OUTPUT
!  data in readable form for further analysis or a plot of the fields
! GTPack OPTIONS
!  o GOVerbose  - some additional output
!
!  o GOPlot     
!      - False  - no plot (standard)
!      - True   - a contour plot of the field
!  
!  o GOStore    
!      - 0      - no output (standard)
!      - name   - output into the file name
!
!  o GOTransformation   
!       - False - import of data to readable form
!       - True  - data are prepared
! 
! GTPack MODULES
!  GTReadFromFile, GTWriteToFile
! GTPack NOTEBOOKS 
!  GTPhMPBBands.nb ,MPB_Symmetry_Analysis.nb
! DESCRIPTION
!  The extension to three dimensions is realized internally. The function should realize ..?
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! BUG FIXES 
!  o 02.07.2018 
!
!    When the Datataset was only "Bloch wavevector" the out put has to be the k-vector. Thi was not the case
!    anymore. The line datstore=impdat solved the problem.
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTPhMPBFields::file = "File `1`is not in this directory."
GTPhMPBFields::plot = "This plot mode is not implemented."

GTPhMPBFields[file_, object_, OptionsPattern[]] := Module[{geom,store,dim,plt,fl,lst,impdat,kvec,data,datstore,resx,resy,resz,i,j,k,np,no,pl1,pl2,pl3,pl4,
	                                                       pld,go,datax,datay,dataz,dens}, 
  (*--- options ---*)	
	geom  = OptionValue[GOTransformation];
    store = OptionValue[GOStore];
    plt   = OptionValue[GOPlot];
  (*---check if file name is valid---*)
    fl = Intersection[FileNames[], {file}];
    If[fl == {},
       Message[GTPhMPBFields::file, file]; Return[],
       None
    ];
  (*---information about objects in data file---*)  
    If[object == "Datasets",
       lst = Import[file, {"Datasets"}];
       Print[lst]; Return[],
       None
    ];
  (*---import data---*) 
    impdat = Import[file, {"Datasets", object}];
    datstore=impdat;
    fl = Intersection[object, {"Bloch wavevector"}];
    If[geom,
       If[fl == {},
          Print["Error: Transformation needs object 'Bloch wavevector' "]; Abort[],
          If[Length[impdat] == 7, 
             dim = 3, 
             dim = 2
          ]
       ];
  (*---correct boundaries with phase factors for 2d and 3D ---*)  
     If[dim == 2,
        kvec = impdat[[1]];
        data = impdat[[2]] + I impdat[[3]];
        AppendTo[data, Exp[I 2 \[Pi] kvec[[1]]] data[[1]]];
        data = Map[Append[#, Exp[I 2 \[Pi] kvec[[2]]] #[[1]]] &, data];
        datstore = {kvec, data},
        None
     ];
     If[dim == 3,
        kvec = impdat[[1]];
        resx = Length[impdat[[2]]];
        resy = Length[impdat[[2, 1]]];
        resz = Length[impdat[[2, 1, 1]]];
        data = Table[
                 {impdat[[2, i, j, k]] + I impdat[[3, i, j, k]],                          
                  impdat[[4, i, j, k]] + I impdat[[5, i, j, k]],                        
                  impdat[[6, i, j, k]] + I impdat[[7, i, j, k]]}
               , {i, 1, resx}, {j, 1, resy}, {k, 1, resz}];
        Do[
           AppendTo[data[[k]], Exp[I 2 \[Pi] kvec[[1]]] data[[k, 1]]];
           data[[k]] = Map[Append[#, Exp[I 2 \[Pi] kvec[[2]]] #[[1]]] &, data[[k]]];
        , {k, 1, resz}];
        data = AppendTo[data, Exp[I 2 \[Pi] kvec[[3]]] data[[1]]];
        datstore = {kvec, data},
        None
     ],
     None
  ];
  If[store == 0,
     None,
     Print["Written to file: ", store];
     GTWriteToFile[datstore, store]
  ];
  (*--- Plot of modes or components of modes for 2D and 3D case ---*)
  If[plt,
     If[fl == {},
        np = 0,
        np = 1
     ];
     no = Length[object] - np;
     Switch[no,
  (*--- real or imaginary part of scalar field --*)
            1,
              ListDensityPlot[impdat[[np + 1]], ColorFunction -> GTBlueRed, ColorFunctionScaling -> True,PlotLabel -> object[[np + 1]]] ,
            2,
  (*--- real and imaginary part of one scalar field or component of \vector field---*)
                 
              pl1 = ListDensityPlot[impdat[[np + 1]], ColorFunction -> GTBlueRed, ColorFunctionScaling -> True, PlotLabel -> object[[np + 1]]] ;
              pl2 = ListDensityPlot[impdat[[np + 2]], ColorFunction -> GTBlueRed, ColorFunctionScaling -> True, PlotLabel -> object[[np + 2]]] ;
              pld = Abs[impdat[[np + 1]] + I impdat[[np + 2]]];
              pl3 = ListDensityPlot[pld, ColorFunction -> GTBlueRed, ColorFunctionScaling -> True, PlotLabel -> "Abs"] ;
              go  = Grid[{{pl1, pl2, pl3}}]; Return[go],
  (*--- all components of the vector field ---*)
            6, 
              datax = impdat[[np + 1]] + I impdat[[np + 2]]; 
              datay = impdat[[np + 3]] + I impdat[[np + 4]];   
              dataz = impdat[[np + 5]] + I impdat[[np + 6]];
              pl1  = ListContourPlot3D[Abs[datax], Mesh -> None, Contours -> 1, ContourStyle -> Directive[ Opacity[0.3], LightBlue], 
                                       PlotLabel -> "Abs x-component"];
              pl2  = ListContourPlot3D[Abs[datay], Mesh -> None, Contours -> 1, ContourStyle -> Directive[ Opacity[0.3], LightBlue], 
                                       PlotLabel -> "Abs y-component"];   
              pl3  = ListContourPlot3D[Abs[dataz], Mesh -> None, Contours -> 1, ContourStyle -> Directive[ Opacity[0.3], LightBlue], 
                                       PlotLabel -> "Abs z-component"];
              dens = datax*Conjugate[datax] + datay*Conjugate[datay] + dataz*Conjugate[dataz] // Chop;
              dens = Sqrt[dens];
              pl4 = ListContourPlot3D[dens, Mesh -> None, Contours -> 1, ContourStyle -> Directive[ Opacity[0.3], LightBlue], 
                                      PlotLabel -> "Density Field"];
               go = Grid[{{pl1, pl2, pl3, pl4}}]; Return[go],
            _,Message[GTPhMPBFields::plot]; Return[]
     ],
     Return[datstore]
  ]
]





(*
***)	
	

(****k* /GTPhSymmetryField
! NAME
!  GTPhSymmetryField
! AUTHOR
!  W. Hergert, S. Thomas, M. Daene
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 03.07.2016 : first version
!  * 14.12.2016 : implementation of 3D case
!  * 02.04.2017 : correction dimensions od problem
!  * 24.02.2023 : MAthematica messages, old version removed
!
! USAGE
!  GTPhSymmetryField[field, character table, number of points] performs the symmetry analysis with respect tp
!  a pointgroup with character table. 
!  
! INPUT
!   o field    - field to analyse
!   o chartab  - character table of the corresponding point group
!   o mpoints  - number of interpolation points in mesh
OUTPUT
!  the name of the Irep or table and picture of the transformed fields
! 
! GTPack OPTIONS
!  o GOVerbose   
!
!      - True   -  some additional output
!      - False  - suppress additional output
!
!  o GOPlot
!     
!      - False - no plot and table (standard)
!      - True  - a table of the deviations and a plot of the fields
!
!  o GOPlotStyle  
!
!      - "Re"   - real part (standard)
!      - "Im"   - imaginary part         
!      - "Abs"  - absolute value
!
! GTPack MODULES
!  
! GTPack NOTEBOOKS 
!  MPB_Symmetry_Analysis.nb
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  14.12.2016 The definition of the interpolation in 3D case can be done in a more clever way probably.
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTPhSymmetryField::field = "Only Re, Im, Abs is allowed." 

GTPhSymmetryField[field_, chartab_, mpoints_, OptionsPattern[]] := Module[
	{dim,verb,plt,style,colf,colfq,colfi,classes,chars,names,symbs,nreps,go,rmats,rmatsi,chfactors,
	 xval,yval,zval,rad,mesh,x,ml,f,f1,f2,f3,pic0,xn,yn,zn,ftrans,temp,ii,temp2,pic1,rep,i,final,
	 norm0,ndevs,normtab,pos,nfigs,scale,frame,label,k,y,pics,ptf},
  (*--- dimension of the problem ---*)	
     If[Length[Dimensions[field]]>2,
     	dim=3,
     	dim=2
     ];
  (*--- options ---*)
     verb   = OptionValue[GOVerbose];
     plt    = OptionValue[GOPlot];
     style  = OptionValue[GOPlotStyle];
     colf   = OptionValue[ColorFunction];
     If[colf=="GrayScale",
     	colfq=True,
     	colfq=False;
     	If[colf=="",
     	   colfi=GTBlueRed,
     	   colfi=colf
     	]
     ]; 	   
  (*--- switch off messages ---*) 
     Off[Power::infy, Infinity::indet, General::stop];
     If[verb,
        Print["Messages: Power::infy, Infinity::indet, General::stop switched off"],
        None
     ];
  (*--- Analyse symmetry input ---*)
     classes = chartab[[1]];
     chars   = chartab[[2]];
     names   = chartab[[3]];
     symbs   = Prepend[names, " "];
     nreps   = Length[classes];
     go      = Length[Flatten[classes]];
  (*--- inverse of the transformation matrices ---*)
     If[dim == 2,
        rmats  = Take[#, {1, 2}, {1, 2}] & /@ GTGetMatrix[Flatten[classes, 1]];
        rmatsi = Inverse[#] & /@ rmats,
        rmats  = GTGetMatrix[Flatten[classes, 1]]; 
        rmatsi = Inverse[#] & /@ rmats
     ];
  (*--- character corresponding to each group element ---*)
     chfactors = Flatten[0 classes + #] & /@ chars;
  (*---generate the mesh---*)
     If[dim == 2,
        xval = Length[field[[1]]];
        yval = Length[field];
        rad  = Min[{xval, yval}]/Sqrt[2.],
        xval = Length[field[[1]]];
        yval = Length[field[[1, 1]]];
        zval = Length[field[[1, 1, 1]]];
        rad  = Min[{xval, yval, zval}]/Sqrt[2.]
     ];
     mesh = Table[x, {x, -rad, rad, 2. rad/(mpoints - 1)}];
     ml   = Length[mesh];
  (*--- interpolating function of the field---*)
  (* warum Transpose ?? *)
     If[dim == 2,
  (*      f  = ListInterpolation[Transpose[field], {{-xval, xval}, {-yval, yval}}],*)
         f  = ListInterpolation[field, {{-xval, xval}, {-yval, yval}}],
  (*--- interpolation for all 3 components of the field , define pure function ---*)
        f1 = ListInterpolation[field[[All, All, All, 1]], {{-xval, xval}, {-yval, yval}, {-zval, zval}}];
        f2 = ListInterpolation[field[[All, All, All, 2]], {{-xval, xval}, {-yval, yval}, {-zval, zval}}];
        f3 = ListInterpolation[field[[All, All, All, 3]], {{-xval, xval}, {-yval, yval}, {-zval, zval}}];
        f   = {f1[#1, #2, #3], f2[#1, #2, #3], f3[#1, #2, #3]} &;
     ];
  (*--- initial field with interpolating function ---*)
     If[dim == 2, 
        pic0 = Table[f @@ {mesh[[xn]], mesh[[yn]]}, {yn, 1, ml}, {xn, 1, ml}],
        pic0 = Table[f @@ {mesh[[xn]], mesh[[yn]], mesh[[zn]]}, {xn, 1, Length[mesh]}, {yn, 1, Length[mesh]}, {zn, 1, Length[mesh]}]
     ];
  (*--- calculate a list of transformed fields ---*)
     If[dim == 2,
        ftrans = Table[
        	              temp = N[#.{mesh[[xn]], mesh[[yn]]} & /@ rmatsi]; f @@ # & /@ temp
                 , {yn, 1, ml}, {xn, 1, ml}],
        ftrans = Table[                    
                       temp  = N[#.{mesh[[xn]], mesh[[yn]], mesh[[zn]]} & /@ rmatsi]; 
                       temp2 = f @@ # & /@ temp;      
                       Table[rmats[[ii]].temp2[[ii]], {ii, 1, Length[temp2]}]
                 , {xn, 1, Length[mesh]}, {yn, 1,Length[mesh]}, {zn, 1, Length[mesh]}]
     ];
  (*---calculate the projections of the fields---*)
     If[dim == 2,
        pic1 = Table[
        	         Table[
        	         	   chars[[rep, 1]]/go chfactors[[rep]].# & /@ ftrans[[i]]
        	         , {i, 1, Length[ftrans]}]
        	   , {rep, 1, nreps}],
        pic1 = Table[
        	         chars[[rep, 1]]/Length[rmatsi] Sum[chfactors[[rep, i]] ftrans[[All, All, All, i]]
        	   , {i, 1, Length[chfactors[[rep]]]}], {rep, 1, Length[chfactors]}]
     ];
  (*--- the complete list of field plots ---*) 
     final = Flatten[{{pic0}, pic1}, 1];
  (*---  calculate Frobenuius Norm ---*)
     norm0 = Flatten[Abs[final[[1]]]]; norm0 = Sqrt[norm0.norm0];
  (*--- relative norm of the deviations- --*)    
     ndevs   = Table[temp = Flatten[Abs[final[[1]] - final[[i]]]]; Sqrt[temp.temp], {i, 2, nreps + 1}]/norm0;
     normtab = {symbs, Prepend[100 (1 - ndevs)/(Plus @@ (1 - ndevs)), "Norm"]};
  (*--- position of the best fit ---*)  
     pos = Flatten[Position[ndevs, Sort[ndevs][[1]]]][[1]];
     If[verb,
        Print["Best fit: ", names[[pos]]],
        None
     ];
     If[plt,
        Print[
              Grid[normtab, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
                            Background -> {{1 -> GTBackGroundColor1}, {1 ->  GTBackGroundColor1}, {1, 1} -> GTCornerColor}
                  ]
             ],
        None
     ];
   (*--- frame for graphical output only for 2D case ---*)
     If[plt && dim == 2,
        nfigs = Length[final];
        scale = Length[final[[1, 1]]];
        frame = Flatten[
        	            {Thickness[0.003], Table[Line[{{i scale + 1/2, 1}, {i scale + 1/2, scale}}], {i, 1,nfigs - 1}], 
        	             Line[{{1, 1}, {scale nfigs, 1}}], 
                         Line[{{1, scale}, {scale nfigs, scale}}], 
                         Line[{{1, 1}, {1, scale}}], 
                         Line[{{scale nfigs, 1}, {scale nfigs, scale}}]
                        }
                       ];
   (*---construct final graphics---*)
       label = {Table[{(i - 1/2) scale, symbs[[i]], 0, RGBColor[1, 1, 0]}, {i, 1, nfigs}], None, None, None};
       pics = Table[
       	            Flatten[Table[final[[k, y]], {k, 1, nfigs}]]
       	      , {y, 1, scale}];
       Switch[style,
             "Re",
                 ptf = Re[pics],
             "Im",
                ptf = Im[pics],
            "Abs",
                ptf = Abs[pics],
               _,
                Message[GTPhSymmetryField::field]; Return[]
       ]; 
       If[colfq,
          ColorConvert[
             ListDensityPlot[ptf, ColorFunction ->GTBlueRed, 
                            AspectRatio -> Automatic, FrameStyle -> Black, FrameTicks -> label,
                            ImageSize -> 500, Epilog -> frame, PlotRange -> Full],"GrayScale"],                  
              ListDensityPlot[ptf, ColorFunction ->colfi, 
                            AspectRatio -> Automatic, FrameStyle -> Black, FrameTicks -> label,
                            ImageSize -> 500, Epilog -> frame, PlotRange -> Full]
         ],
         names[[pos]]
     ]
 ]
   


(*
***)	
	

(****k* /GTToString
! NAME
!  GTToString
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  03/07/2016 : first version
!
! USAGE
!  GTToString[number]converts an inter to a string in such aform, thati it can be used
!  for the construction of MPB file names
! INPUT
!   number - integer
OUTPUT
!  string
! ERROR MESSAGES
!  -
! GTPack OPTIONS
!  -
! GTPack MODULES
! -
! GTPack NOTEBOOKS 
! -
! DESCRIPTION
!  In MPB file names the number start with 01,02,...,09,10,... 
!  The module constructes strings in this scheme from intergers 1,2,3,4,...
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


GTToString[num_] := Module[{str},
  If[num >= 10,
     str = ToString[num],
     str = "0" <> ToString[num]
  ]; 
  Return[str]
]

(*
***)

(****k* /GTPhSymmetryPoint
! NAME
!  GTPhSymmetryPoint
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 03.07.2016 : first version
!  * 22.052016 : output slightly changed, changes caused by changes inGTPhMPBFields (renamed option) implemented.
!  * 25.08.2016 : check of basis ->basis1 to handle 2D case correctly
!  * 29.06.2018 : check header and documentation
!  * 24.08.2018 : recalculation of k-vector to Cartesian basis added
! USAGE
!  GTPhSymmetryPoint[filen, minb, maxb, objects, nmesh, group, basis] performs the symmetry analysis at a certain k-point
!  for band numbers min <= n <= maxb. objects defines the fields taken from the h5 files. nmesh sets the mesh size.
!  group is the point group of the photonic crystal. basis defines the basis of the reciprocal lattice.
!
! INPUT
!   o filen      - the filename in the MPB format
!   o minb,maxb  - minimum and maximum band number
!   o objects    - objects in h5 file to use, i.e. {"Bloch wavevector","z.r-new","z.i-new"}
!   o nmesh      - mesh size for interpolation
!   o group      - point group of the photonic crystal
!   o basis      - basis of the reciprocal lattice
OUTPUT
!  list of Ireps according to the analysed bands
! GTPack OPTIONS
!   o GOVerbose   
!
!      - True   -  some additional output
!      - False  - suppress additional output
!
!   o GOPlot
!     
!      - False - no plot and table (standard)
!      - True  - a table of the deviations and a plot of the fields
!
!   o GOPlotStyle  
!
!      - "Re"   - real part (standard)
!      - "Im"   - imaginary part         
!      - "Abs"  - absolute value
!
!   o GOIrepNotation - notation od Ireps to use (standard is "Mulliken")
!
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTPhMPBFields, GTPhSymmetryField, GTToString
! GTPack NOTEBOOKS 
!  MPB_Symmetry_Analysis.nb
! DESCRIPTION
!   The module reads directly h5 files, i.e. it is designed to be used to analyse MPB results.
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  The k-vector is taken from the MPB file to select the correct little group. A recalculation to Cartesian coordinates 
!  is implemented to get the correct little group.
!
! SOURCE
!--------------------------------------------------------------------------------
!
*)


GTPhSymmetryPoint[filen_, minb_, maxb_, objects_, nmesh_, group_, basis_, OptionsPattern[]] :=
    Module[{verb, chars, pos, part1, part2, repname, fu, fo, objlist, kvec, field, gofk, ctab, ireplist, irep, fname, i,plot, style,
    	    bds,nms,kvec1},
    (*--- options ---*)
       verb    = OptionValue[GOVerbose];
       repname = OptionValue[GOIrepNotation];
       plot    = OptionValue[GOPlot];
       style   = OptionValue[GOPlotStyle];
    (*--- analyse file names ---*)
       chars = Characters[filen]; 
       pos = Flatten[Position[chars, "b"]][[1]];
       part1 = StringJoin[Take[chars, {1, pos}]]; 
       part2 = StringJoin[Take[chars, {pos + 3, Length[chars]}]];
       If[verb,
          fu = part1 <> GTToString[minb] <> part2; 
          fo = part1 <> GTToString[maxb] <> part2;
          Print["Files from ", fu, " to ", fo],
          None
       ];
       objlist = Prepend[objects, "Bloch wavevector"];
       {kvec, field} = GTPhMPBFields[filen, objlist, GOTransformation -> True];
       If[Length[kvec]==3&&Length[basis[[1]]]==2,
       	 kvec=Take[kvec,{1,2}],
       	 None
       ];
     (*--- recalculate kvec: relative to basis vectors -> Cartesian ---*)
        kvec1 = kvec; 
        kvec=Total[kvec*basis];
       If[verb,
          Print["k-vector reciprocal basis  ", kvec1];
          Print["K-vector in Cartesian basis", kvec]
       ];
    (*--- construct group of wavevector and character table ---*)
       gofk = GTGroupOfK[group, kvec, basis];
       If[verb,
       	  Print["Group of k : ",gofk],
       	  None
       ];
       ctab = GTCharacterTable[gofk, GOVerbose -> verb, GOIrepNotation -> repname];
       ireplist = {};
       Do[
          fname = part1 <> GTToString[i] <> part2;
          field = GTPhMPBFields[fname, objlist, GOTransformation -> True];
          irep  = GTPhSymmetryField[field[[2]], ctab, nmesh, GOVerbose -> verb,GOPlot->plot,GOPlotStyle->style];
          AppendTo[ireplist, irep]
       , {i, minb, maxb}];
       If[plot,
          bds = Table[i, {i,minb,maxb}]; nms = Table["band : ", {i,minb,maxb}];
          Grid[{nms, bds, ireplist} // Transpose],
          Return[ireplist]
       ]
  ]

	

(*
***)	



(****k* /GTPhSymmetryBands
! NAME
!  GTPhSymmetryBands
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  03/07/2016 : first version
!  05/22/2016 : Correction of minor mistakes in output
!  02.04.2017 : option precision in table included, call of GTPhMPBBands modified.
!  04.04.2020 : plot didn't contain lines at symmetry points, klist add to bandstr before plot.
! USAGE
!  GTPhSymmetryBands[fileb, filef, klist, minb, maxb, objects, nmesh, group, basis] performs the symmetry analysis  at
!  a series of point s of a bandstructure. The results are presented in a table. Also a plot of the band structure with the assigned Ireps
!  is possible.
!
! INPUT
!   fileb      - band structure from MPB processed for plot mit GTBanndsPlot
!   filef      - the filename  for the field in  MPB format
!   klist      - list of k-points for field analysis, for example: {{1, "\[CapitalGamma]"}, {52, "X"}, {103, "M"}, {154, "\[CapitalGamma]"}}
!   minb,maxb  - minimum and maximum band number
!   objects    - objects in h5 file to use, i.e. {"Bloch wavevector","z.r-new","z.i-new"}
!   nmesh      - mesh size for interpolation
!   group      - point group of the photonic crystal
!   basis      - basis of the reciprocal lattice
OUTPUT
!  table with reslts of analysis, band strucutre with Ireps
! ERROR MESSAGES
!   -
! GTPack OPTIONS
!  GOVerbose      - some additional output
!  GOIrepNotation - notation od Ireps to use (standard is "Mulliken")
!  GOPlot         - False - no plot and table (standard)
!                 - True  - a table of the deviations and a plot of the fields
!  GOLabelShift   - shift the positions of the Irep labels
!  PlotStyle
!  Joined
!  PlotRange 
! GTPack MODULES
!  GTPhMPBBands, GTPhSymmetryPoint, GTBandsPlot
! GTPack NOTEBOOKS 
!  MPB_Symmetry_Analysis.nb
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


GTPhSymmetryBands[fileb_, filef_, klist_, minb_, maxb_, objects_, nmesh_, group_, basis_, OptionsPattern[]] :=
    
 Module[{verb,repname,style,range,join,del1,del2,nkp,bandstr,bands,bndstab,nb,k,pos,chars,part1,part2,flist,
 	     tab,ireptab,ireplist,i,elist,blist,plot,clist,plist,bplot,text,j,tab1,kk1,nbnds,plabel,flabel,nround,lstyle},
  (*--- options ---*)
    verb        = OptionValue[GOVerbose];
    repname     = OptionValue[GOIrepNotation];
    style       = OptionValue[PlotStyle];
    range       = OptionValue[PlotRange];
    join        = OptionValue[Joined];
    plot        = OptionValue[GOPlot];
    plabel      = OptionValue[PlotLabel];
    nround      = OptionValue[GOPrecision];
    flabel      = OptionValue[FrameLabel];
    lstyle      = OptionValue[GOLabelStyle];
    {del1,del2} = OptionValue[GOLabelShift];
    nkp         = Length[klist];
  (*--- input bandstructure ---*)
    bandstr = GTPhMPBBands[fileb, klist,basis, GOVerbose -> False];
    bands = bandstr[[1]]; bndstab = {}; nbnds=Length[bands[[1,4]]];
   
  (*--- select energies of points to analyse ---*)
    Do[
       nb = klist[[k, 1]];
       pos = Flatten[Position[Transpose[bands][[1]], nb]][[1]];
       ew     = Round[#*10.^nround]/10.^nround& /@bands[[pos, 4]];
       bndstab = Append[bndstab, ew]
(*       bndstab = Append[bndstab, bands[[pos, 4]]] *)
    , {k, 1, nkp}];
  (*--- names of fields ---*)
       chars = Characters[filef]; 
       pos   = Flatten[Position[chars, "k"]][[1]];
       part1 = StringJoin[Take[chars, {1, pos}]]; 
       part2 = StringJoin[Take[chars, {pos + 3, Length[chars]}]];
       flist = {};
       Do[
          flist = Append[flist, part1 <> GTToString[klist[[k, 1]]] <> part2]
       , {k, 1, nkp}];
       If[verb,
          Print["Fields to analyse: ", flist],
          None
       ];
  (*--- symmetry analysis ---*)
       tab1 = {}; ireptab = {};
       Do[
          ireplist = GTPhSymmetryPoint[flist[[k]], minb, maxb, objects, nmesh, group, basis, GOVerbose -> verb, GOIrepNotation -> repname];
          ireptab  = Append[ireptab, ireplist];
          tab1      = Append[tab1, {ireplist, Take[bndstab[[k]],{minb,maxb}]}]
       , {k, 1, nkp}];
  (*--- construct a table with results ---*)
       elist = Table["E(k)", {2*nkp}]; kk1 = Flatten[klist];
       Do[
          elist[[2 i-1]] = kk1[[2 i]]
       , {i, 1, nkp}];
       elist = Prepend[elist, " "];
       blist = Table[i, {i, minb, maxb}]; 
       tab = Prepend[Transpose[Prepend[Flatten[tab1, 1], blist]], elist] // Transpose;
           Print[Grid[tab, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, 
     	       {2 -> GTDividerColor1}},Background -> {{1 -> GTBackGroundColor1}, 
               {GTBackGroundColor1, {GTBackGroundColor2,White}},{1, 1} -> GTCornerColor}
              ]
  (*     Print[Grid[tab, Frame -> All, Dividers -> {{2 -> GTDividerColor1}, {2 -> GTDividerColor1}}, 
             Background -> {{1 -> GTBackGroundColor1}, {1 -> GTBackGroundColor1}, {1, 1} -> GTCornerColor}]*)
       ];
  (*--- prepare the plot ---*)
        clist = Transpose[bands][[2]]; 
        plist = Map[clist[[#]] + del1 &, Transpose[klist][[1]]];
        text = {}; 
        bandstr[[2]] = klist;
          Do[
        	 Do[
       	  	  (*  tt   = Text[ireptab[[i, j]], {plist[[i]], bndstab[[i, j+minb-1]] + del2}]; 
                text = Append[text, tt]*)
                If[lstyle=={},
                   text  = Append[text, Text[ireptab[[i, j]], {plist[[i]], bndstab[[i, j+minb-1]] + del2}]],
                   text  = Append[text, Text[Style[ireptab[[i, j]],lstyle], {plist[[i]], bndstab[[i, j+minb-1]] + del2}]]
              
           ];
             , {j, 1,maxb-minb+1}]
          , {i, 1, nkp}];
       If[plot,
         
          bplot = GTBandsPlot[bandstr, nbnds, PlotStyle -> style, Joined -> join, PlotRange -> range,FrameLabel->flabel,PlotLabel->plabel];
          Show[bplot, Graphics[text]],
          Return[{ireptab,bndstab,text}]
      ]
  ]





(*
***)	



(****k* /GTPhBandsObjects
! NAME
!  GTPhBandsObjects
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 03.07.2016 : first version
!  * 22.05.2016 : Correction of minor mistakes in output
!  * 11.08.2016 : output of table restricted to 15 bands
!  * 29.06.2018 : check header and documenttion
!  * 24.02.2023 : Mathematica messages
! USAGE
!  GTPhBandsObjects[structure, permitivity background, reciprocal basis, cutoff, kpoints, number of bands] calculates the photonic
!  bandstructure at kpoints for a structure defined by a list of objects."
!
! INPUT
!   o struct     - defines the structure as a list of objects
!   o eback      - permittivity of background
!   o basvec     - basis of reciprocal lattice 
!   o cut        - cutoff of G-vector cluster
!   o kpoints    - list of k-points
!   o nev        - number of bands to calculate
!  
! OUTPUT
! photonic bandstructure, ready to use in bandstructure plots or calulation of DOS
!
! GTPack OPTIONS
!   o GOPhPol
!
!     - "Automatic"  - no special polarization (1D or 3D)
!     - "E" or "TM"  - TM polarization
!     - "H" or "TE"  - TE polarization
!
!   o GOVerbose
!
!     - True  - additional information (standard)
!     - False - no additional information
!
!  o GOPlotBands   
!      -  True  - calculation for bandstructure 
!      -  False - calculation for DOS
!  
!  o GOStore    - stores to filename, if filename is given 
!
! GTPack MODULES
!  GTLatCluster, GTLatShells, GTPhDCObjects
! GTPack NOTEBOOKS 
!  MPB_Symmetry_Analysis.nb
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
GTPhBandsObjects::method   = "GODCMethod describes wrong method."
GTPhBandsObjects::struc    = "Structure `1 is not implemented."
GTPhBandsObjects::objects  = "Objects `1  not implemented."
GTPhBandsObjects::pol      = "Polarization `1  is wrong."
GTPhBandsObjects::pol2D    = "Wrong polarization for two-dimensional case."
GTPhBandsObjects::problem  = "Problem not defined."


GTPhBandsObjects[struct_, eback_, basvec_, cut_, kpoints_, nev_, OptionsPattern[]] := Module[
	{eps,c,poltypes,bands,types,opol,nplt,ndt,mho,gvec,epsmat,kvec,lx,klx,kt,verb,ndim,nklx,nkp,kv,k,
	 ev,info,test,dim,ham,gi,gj,i,j,gv,fteps,gv1,ngv1,pv,npv,pv1,pol2,ki,kj,p1,p2,nbds,npt,tab,ind,
	 nrd,tmp},
  (*---set constants and definitions---*)
	eps      = 10^(-7);
    c        = {1.1, 3.7, 2.01}; c = c/Norm[c];
    poltypes = {"E", "H", "TM", "TE", "Automatic"};
    types    = {"CircularRod", "EllipticalRod", "RectangularRod", 
                "PrismaticRod", "Sphere", "Cuboid", "Slab", "EmptyLattice", 
                "SlabSmooth", "RodSmooth"};
    bands = {}; 
  (*---Options---*)
    verb = OptionValue[GOVerbose];
    opol = OptionValue[GOPhPol];
    nplt = OptionValue[GOPlotBands];
    ndt  = OptionValue[GOStore];
    mho  = OptionValue[GODCMethod];
    nrd  = OptionValue[GODecimals];
    If[mho == "MatrixInverse" || mho == "Direct",
      None,
      Message[GTPhBandsObjects::method]; Abort[]
    ];
  (*--- calculation of G-vectors and eps matrix ---*)  
    {gvec,epsmat} = GTPhPermittivityMatrix[basvec, cut, struct, eback,GODCMethod->mho];
    ndim = Length[gvec];
    If[verb, 
       Print["number of G-vectors : ", ndim], 
       None
    ];
  (*---prepare the kpoints----*)
    If[nplt, 
       kvec = Transpose[kpoints[[1]]][[3]]; 
       lx   = kpoints[[2]];
       klx  = lx[[1, 1]]; kt = 1, 
       kvec = kpoints
    ];
    nklx = Length[lx];
    nkp = Length[kvec];
    If[verb, 
       Print["number of k-vectors : ", nkp], 
       None
    ];
  (*---empty lattice---*)
    If[Head[struct] === String, 
       If[struct == "EmptyLattice", 
          If[verb, 
          	 Print["Empty lattice band structure"], 
          	 None
          ];
          Do[
          	 kv    = kvec[[k]];
             ev    = Map[Norm[kv + #] &, gvec];
             ev    = Take[Sort[ev, Less], {1, nev}];
             bands = Append[bands, {kpoints[[1, k, 1]], kpoints[[1, k, 2]], kv, ev}]
          , {k, 1, nkp}], 
          Message[GTPhBandsObjects::struc]; Abort[]
       ], 
       None
    ];
  (*---objets in elementary cell---*)
    If[Head[struct] == List, 
       info = Transpose[struct][[1]];
       If[verb, 
       	  Print["structures in elementary cell: ", info], 
       	  None
       ];
  (*---check type of objects---*)
       test = Complement[info, types];
       If[test == {}, 
          None, 
          Message[GTPhBandsObjects::objects,test]; Abort[]
       ], 
       None
    ];
  (*---check dimension of the problem---*)dim = 3;
    test = Map[Union[#] &, Transpose[gvec]];
    If[test[[3]] == {0}, 
       dim = dim - 1, 
       None
    ];
    If[test[[2]] == {0}, 
       dim = dim - 1, 
       None
    ];
    If[verb, 
       Print["Dimensionality: ", dim], 
       None
    ];
  (*---check polarization---*)
    If[Position[poltypes, opol] == {}, 
    	Message[GTPhBandsObjects::pol,opol]; Abort[], 
       None
    ];
    If[opol == "E", opol = "TM";
      If[verb, 
      	Print["Polarization E corresponds to TM"], 
      	None
      ], 
      None
    ];
    If[opol == "H", 
       opol = "TE";
       If[verb, 
       	  Print["Polarization H corresponds to TE"], 
       	  None
       ], 
       None
    ];
  (*---1D,2D,3D problems---*)
    Switch[dim,
  (*---1D problem---*)
  1, 
    Do[
       kv = kvec[[k]];
       ham = Table[0, {ndim}, {ndim}];
       Do[
          gi = gvec[[i]];
          Do[
       	     gj    = gvec[[j]]; 
       	     gv    = gi - gj;
             fteps = epsmat[[i,j]];
             ham[[i, j]] = fteps*(kv + gi).(kv + gj);
             ham[[j, i]] = Conjugate[ham[[i, j]]]
          , {j, 1, i}]
       , {i, 1, ndim}];
       ev = Eigenvalues[N[ham], -nev] // Chop;
       ev = Sort[Map[Sqrt[#] &, ev], Less];
       If[nplt, 
          bands = Append[bands, {kpoints[[1, k, 1]], kpoints[[1, k, 2]], kvec[[k]], ev}], 
          bands = Append[bands, {kvec[[k]], ev}]
       ]
    , {k, 1, nkp}],
  (*---2D-problem---*)
  2, 
    If[(dim == 2 && opol == "TM") || (dim == 2 && opol == "TE"), 
       None, 
       Message[GTPhBandsObjects::pol2D]; Abort[]
    ];
    Do[
      If[(verb && nplt) && k == klx, 
         Print["k-point ", k, " symmetry point ", lx[[kt, 2]]];
         kt = kt + 1;
         If[kt <= nklx, 
      	    klx = lx[[kt, 1]], 
      	    None
         ], 
         None
      ];
      If[verb && ! nplt, 
         If[k == IntegerPart[nkp/2.], Print["50% calculated"], 
         	None
         ], 
         None
      ];
      kv  = kvec[[k]];
      ham = Table[0, {ndim}, {ndim}];
  (*--- avoid numerical problems ---*)
      If[Norm[kv] < eps, 
      	 kv = kv + {eps, 0, 0}, 
      	 None
      ];
      If[opol == "TM", 
      	 Do[
      	 	gi = gvec[[i]];
            Do[
               gj    = gvec[[j]]; 
               gv    = gi - gj;
               fteps = epsmat[[i,j]];
               ham[[i, j]] = fteps*Norm[kv + gj]*Norm[kv + gi];
               ham[[j, i]] = Conjugate[ham[[i, j]]]
            , {j, 1, i}]
         , {i, 1, ndim}], 
         None
      ];
      If[opol == "TE", 
      	 Do[
      	 	gi = gvec[[i]];
            Do[
               gj    = gvec[[j]]; 
               gv    = gi - gj;
               fteps = epsmat[[i,j]];
               ham[[i, j]] = fteps*(kv + gi).(kv + gj);
               ham[[j, i]] = Conjugate[ham[[i, j]]]
            , {j, 1, i}]
         , {i, 1, ndim}], 
         None
      ];
      ev = Eigenvalues[N[ham], -nev] // Chop;
      ev = Sort[Map[Sqrt[#] &, ev], Less];
      If[nplt, 
         bands = Append[bands, {kpoints[[1, k, 1]], kpoints[[1, k, 2]], kvec[[k]], ev}], 
         bands = Append[bands, {kvec[[k]], ev}]
      ]
    , {k, 1, nkp}],
  (*---3D-Problem---*)
   3, 
     Do[
     	If[(verb && nplt) && k == klx, 
           Print["k-point ", k, " symmetry point ", lx[[kt, 2]]];
           kt = kt + 1;
           If[kt <= nklx,
           	  klx = lx[[kt, 1]], 
           	  None
           ], 
           None
        ];
        If[verb && ! nplt, 
           If[k == IntegerPart[nkp/2.], 
           	  Print["50% calculated"], 
           	  None
           ], 
           None
        ];
        kv = kvec[[k]];
  (*--- avoid numerical problems ---*)
        If[Norm[kv] < eps, 
        	kv = kv + {eps, 0, 0}, 
        	None
        ];
        gv1  = Map[kv + # &, gvec];
        ngv1 = Map[Norm[#] &, gv1];
        pv1  = Map[Cross[c, #] &, gv1];
        npv  = Map[Norm[#] &, pv1];
        pv1  = pv1/npv; pv = {};
        Do[
           pol2 = Cross[pv1[[i]], gv1[[i]]];
           npv  = Norm[pol2];
           pv   = Append[pv, {pv1[[i]], pol2/npv}], {i, 1, ndim}];
           ham  = Table[0, {2*ndim}, {2*ndim}];
           Do[
           	  Do[
           	     Do[
           	   	    Do[
                       fteps = epsmat[[i, j]];
                       ki = i + (p1 - 1)*ndim;
                       kj = j + (p2 - 1)*ndim;
                       ham[[ki, kj]] = fteps*ngv1[[i]]*ngv1[[j]]*pv[[i, p1]].pv[[j, p2]]*(-1)^(p1 + p2);
                       ham[[kj, ki]] = Conjugate[ham[[ki, kj]]]
                   , {j, 1, i}]
                , {i, 1, ndim}]
              , {p1, 1, 2}]
           , {p2, 1, 2}];
           ev = Eigenvalues[N[ham], -nev] // Chop;
           ev = Sort[Map[Sqrt[#] &, ev], Less];
           If[nplt, 
              bands = Append[bands, {kpoints[[1, k, 1]], kpoints[[1, k, 2]], kv, ev}], 
              bands = Append[bands, {kvec[[k]], ev}]
           ]
     , {k, 1, nkp}],
  (*---default case---*)
  _, 
    Message[GTPhBandsObjects::problem]; Abort[]
    ];
  (*---Print list of eigenvalues---*)
   If[verb, 
      If[nev > 15, 
         nbds = 15; 
         Print["Warning: Output restricted to 12 bands!"],
         nbds = nev
      ];
      If[nplt,
  (*---symmetry points only---*)
        npt = Length[kpoints[[2]]];
        tab = Table[Table[0, {nbds + 2}], {npt}];
        Do[
           ind         = lx[[i, 1]];
           tab[[i, 1]] = lx[[i, 2]];
           tab[[i, 2]] = bands[[ind, 3]];
        Do[
           tmp             = bands[[ind, 4, k]];
           tab[[i, k + 2]] = ToString[PaddedForm[N[tmp], {nrd, nrd}]];
        , {k, 1, nbds}]
     , {i,1, npt}];
     Print[
           Grid[Join[{{Text["Label"], Text["k-point"], Text["Eigenvalues"]}}, tab], 
           	    Frame -> True, Background -> {{GTBackGroundColor1}, {GTBackGroundColor1}}, 
                Dividers -> {{2 -> GTDividerColor1, 3 -> GTDividerColor1}, {All,2 -> GTDividerColor1}}
               ]
          ],
  (*---all k-points---*)
     npt = Length[kvec];
     tab = Table[Table[0, {nbds + 1}], {npt}];
     Do[
     	tab[[i, 1]] = bands[[i, 1]];
        Do[
           tmp              = bands[[i, 2, k]];
           tab[[i, k + 1]] = ToString[PaddedForm[N[tmp], {nrd, nrd}]]
        , {k, 1, nbds}]
     , {i, 1, npt}];
     Print[
           Grid[Join[{{Text["k-point"], Text["Eigenvalues"]}}, tab], 
          Frame -> True, Background -> {{GTBackGroundColor1}, {GTBackGroundColor1}}, 
          Dividers -> {{2 -> GTDividerColor1}, {All, 2 -> GTDividerColor1}}]]], None];
  (*---store to files---*)
     If[Head[ndt] === String, 
     	Print["Write to file: ", ndt];
        If[nplt, 
           GTWriteToFile[{bands, lx}, ndt], 
           GTWriteToFile[bands, ndt]
        ], 
        None
     ];
  (*---Return data---*)
     If[nplt, 
        Return[{bands, lx}], 
        Return[bands]
     ]
 ]

(*
***)	


(****k* /GTPhUnCoupledBands
! NAME
!  GTPhUnCoupledBands
! AUTHOR
!  W. Hergert
! PACKAGE
!  Photonics.m 
! MODIFICATION HISTORY
!  * 15.05.2016 : first implementation
!  * 01.07.2018 : check header and implement documentation
!                 a bug is fixed, see code
!  * 11.09.2018 : bug fix
! USAGE
!  GTPhBandsObjects[gvec, bands,rules] cuts the uncoupled bands out of bands. the reciprocal lattice vectors gvec are used to construct
!  the Hamiltonian. rules are used to find the uncopled bands.
!
! INPUT
!   * gvec     - defines the structure as a list of objects
!   * bands    - permittivity of background
!   * mat      - transformation matrix {x,y,z}->{x',y',z'} to express the symmetry
!  
OUTPUT
!  band structure without uncoupled bands
! GTPack OPTIONS
!  * GOVerbose     
!  
!      - True   - additional output
!      - False  - no additional output
!
!  GOPlotBands    - calculation for bandstructure or DOS
!
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
! 
! PROBLEMS
!  11.09.2018: I could not solve the problems with the rule and the internal functions. Now a transformation matrix in the input 
!              is used to express the transformation rule
! SOURCE
!--------------------------------------------------------------------------------
!
*)


GTPhUncoupledBands[gvec_, bands_, mat_, OptionsPattern[{}]] := Module[
 {f, g,pltb,ebands,pbands,kpts,evecs,nkp,min,nbnds,ng,verb,tab,bandsn,bd,vec,kp,vec1,nb,tv,lb,extr,head,res,res1,k,x,y,z},
  pltb = OptionValue[GOPlotBands];
  verb = OptionValue[GOVerbose];
  (*--- extract bands and eigenvectors ---*)
  If[pltb,
       ebands = bands[[1, 1]];
       pbands = Transpose[ebands][[4]];  
       kpts   = Transpose[bands[[2]]][[1]];
       evecs  = Transpose[bands[[2]]][[2]],
       pbands = Transpose[bands[[1]]][[2]];
       kpts   = Transpose[bands[[1]]][[1]];
       evecs  = Transpose[bands[[2]]][[2]]
   ];
  nkp = Length[kpts]; 
  min = nbnds = Length[evecs[[1]]]; 
  ng = Length[evecs[[1, 1]]];
  bandsn = Table[0, {nkp - 2}];
  Do[
  	 bd   = {};
  	 vec  = Exp[I 2 \[Pi] (kpts[[kp]] + #).{x, y, z}] & /@ gvec;
     vec1 = Exp[I 2 \[Pi] (kpts[[kp]] + #).(mat.{x, y, z})] & /@ gvec;
     Do[
        f  = evecs[[kp, nb]].vec;
        g  = evecs[[kp, nb]].vec1;
        tv = f - g // Chop;
        If[
        	   Head[tv] === Integer || Abs[tv] < 10.^(-7),
           bd = Append[bd, pbands[[kp, nb]]],
           None
        ]
     , {nb, 1, nbnds}];
     bandsn[[kp - 1]] = bd
  , {kp, 2, nkp - 1}];
  lb   = Map[Length[#] &, bandsn] // Min; 
  If[verb, 
  	 tab = {{"Number of k-points", nkp}, {"number of bands", nbnds}, {"number of uncoupled bands", nbnds - lb}, {"recip. lattice vectors", ng}};
     Print[Grid[tab, Frame -> All, Background -> {1 -> LightRed}]], 
     None
  ];
  extr = Map[Take[#, {1, lb}] &, bandsn];
  If[pltb,
     head = Take[Transpose[ebands], {1, 3}] // Transpose;
     res  = Table[{head[[k, 1]], head[[k, 2]], head[[k, 3]], extr[[k - 1]]}, {k, 2, nkp - 1}];
     res1 = Prepend[res, {head[[1, 1]], head[[1, 2]], head[[1, 3]], extr[[2]]}];
     res  = Append[res1, {head[[nkp, 1]], head[[nkp, 2]], head[[nkp, 3]], extr[[nkp - 2]]}];
     res  = {res, bands[[1, 2]]},
     res  = Table[{kpts[[k]], extr[[k - 1]]}, {k, 2, nkp - 1}]
  ];
  Return[res]
  ]

(*
***)

   
End[] (* End Private Context *)


EndPackage[]
(*
***)
