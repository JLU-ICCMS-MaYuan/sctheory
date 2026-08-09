(****m* /ElectronicStructure.m
!
! NAME
!  ElectronicStructure.m
! AUTHOR
!  W. Hergert, M. Geilhufe
! MODIFICATION HISTORY
!  * 24.09.2013 - initially created and documented  
!  * 16.08.2015 - check of documentation 
!  * 28.12.2017 - check header and documentation
!  * 22.06.2018 - check headers and documentation
!  * 24.06.2018 - check completed
! USAGE
!  Contains modules for electronic structure calculations.
!
! GTPack MODULES
!  
! --- Calculation of Band Structures ---
!
!  * GTBands               - band energies for a given Hamiltonian at a set of k-points
!  * GTBandStructure       - calculates and plots a bandstructure
!  * GTBandsPlot           - takes data from file and plots a bandstructure
!  * GTFatBandsPlot        - plot of fat bands
!  * GTBandsPlotImprove    - improves a plot from symmetry analysis
!
! --- Calculation of Density Of States ---
!
!  * GTDensityOfStates	   - calculates and plots the DOS, Hamiltonian is given
!  * GTPartialDOS		   - calcualtes the partial density of states
!  * GTDensityOfStatesRS   - calulates DOS from real space Hamiltonian
!  * GTDensityOfStatesPlot - calculates DOS from precalculated energy bands 	
!  * GTFermiSurface        - calculates the Fermi surface	
!  * GTFermiSurfaceCut	   - calculates a cut through the Fermi surface
!  * GTBandsDOSPlot        - creates combined plots of bands, DOS, IDOS and partial DOS
!
! ---  Group Theory and Band Structures ---
!
!  * GTCompatibility       - constructs compatibility relations
!
! ---  internal functions ---
!
!  * GTSortEigensystem     - sorts output of eigensystem calculation
!  * BZGauss               - sum of Gaussians for DOS calculation 
!
! DESCRIPTION
!  The commands are not method specific, i.e. different Hamiltonians can be used 
!  to calculate a band structure or the density of states, for instance. Therefore it works for 
!  tight-binding or pseudopotential models and photonic bandstructures as well. Also phonons
!  can be considered.
!
! LITERATURE
!  W. Hergert, M. Geilhufe 
!  Group Theory in Solid State Physics and Photonics. - Problem Solving
!  with Mathematica 
!
***)


BeginPackage["GroupTheory`ElectronicStructure`",{"GroupTheory`Basic`","GroupTheory`Symbols`","GroupTheory`Install`","GroupTheory`Lattice`","GroupTheory`TightBinding`","GroupTheory`RepresentationTheory`","GroupTheory`Auxiliary`","GroupTheory`ThirdParty`"}]

(*---------------- Calculations of Band Structures ---------------*) 
 GTBands                 ::usage = "GTBands[\*StyleBox[\"Hamiltionian,kpoints,nev\",\"TI\"]] calculates the lowest \*StyleBox[\"nev\",\"TI\"] band energies for a given \*StyleBox[\"Hamiltionian\",\"TI\"] at certain \*StyleBox[\"kpoints\",\"TI\"]."
 GTBandsPlot             ::usage = "GTBandsPlot[\*StyleBox[\"file,nband\",\"TI\"]] plots \*StyleBox[\"nband\",\"TI\"] bands from the data in \*StyleBox[\"file\",\"TI\"]."
 GTBandsPlotImprove      ::usage = "GTBandsPlotImpove[\*StyleBox[\"bands,text,mod,minb,maxb\",\"TI\"]] helps to improve the positions of the Irep labels in a bandstructure plot."
 GTBandStructure         ::usage = "GTBandStructure[\*StyleBox[\"Hamiltonian,kpath,npoints,nband\",\"TI\"]] calculates and plots the band structure for a given \*StyleBox[\"Hamiltionian\",\"TI\"]."
 GTFatBandsPlot          ::usage = "GTFatBandsPlot[\*StyleBox[\"bands,nband,scale\",\"TI\"]] plots \*StyleBox[\"nband\",\"TI\"] fat bands from the data in \*StyleBox[\"bands\",\"TI\"]. \*StyleBox[\"scale\",\"TI\"] scales thickness of the fat bands."
(*---------------- Calculations of Density of States and Fermi surface ------------*)
 GTDensityOfStates       ::usage = "GTDensityOfStates[\*StyleBox[\"Hamiltonian,structure,parameters\",\"TI\"]] calculates and plots the density of states for a given \*StyleBox[\"Hamiltonian\",\"TI\"]."
 GTPartialDOS            ::usage = "GTPartialDOS[\*StyleBox[\"Hamiltonian,structure,parameters\",\"TI\"]] calculates and plots partial density of states for a given \*StyleBox[\"Hamiltonian\",\"TI\"]."
 GTDensityOfStatesRS     ::usage = "GTDensityOfStatesRS[\*StyleBox[\"Hamiltonian,parameters\",\"TI\"]] is used to calculate the density of states of a real space \*StyleBox[\"Hamiltonian\",\"TI\"]."
 GTDensityOfStatesPlot   ::usage = "GTDensityOfStatesPlot[\*StyleBox[\"band structure,parameters\",\"TI\"]] calculates the density of states from a precalculated \*StyleBox[\"band structure\",\"TI\"]."
 GTFermiSurface          ::usage = "GTFermiSurface[Hamiltonian, Fermi energy,  list of bands, ndel,  kbasis, clusterdata, kpath ] calculates the Fermi surface corresponding to a Hamiltonian and Fermi energy if the Fermi surface 
                                    contains parts from list of bands. The electronic structure is calculated in a cube at ndel points per spatial dimension. kbasis is the basis of the reciprocal lattice. clusterdata contains the data for the lattice construction. The path used in electronic structure calculations can be given by BZpath."
 GTFermiSurfaceCut       ::usage = "GTFermiSurfaceCut[\*StyleBox[\"Hamiltonian, Fermi energy,  list of bands, ndel,  transformation,sizze of region\",\"TI\"]] calculates a cut through the Fermi surface on a given cutting plane."
 GTBandsDOSPlot          ::usage = "GTBandsDOSPlot[\*StyleBox[\"bands, dos, nbands\",\"TI\"]] creates a combined plot of bandstructure and density of states."
(*---------------- Group Theory and Band Structures -------------*)
 GTCompatibility         ::usage = "GTCompatibility[\*StyleBox[\"group1,group2\",\"TI\"]] calculates the compatibility relations of two point groups."

(*--------------------------- Options --------------------------------*)
 Options[GTBands]                    = {GOEigenvectors -> False, GOStore -> 0, GOVerbose -> False, GOTbOrthogonal -> True, GOPlotBands -> True, GOPhotonic -> False, GOShift -> 0,GODecimals->4}
 Options[GTBandsPlot]                = {PlotRange->All,Joined->False,FrameLabel->{" ","Energy"},PlotLabel->"Band structure",GOShift->0,PlotStyle->Blue,GOFermiEnergy->0}
 Options[GTBandsPlotImprove]         = {PlotRange -> All, PlotStyle -> {{Thin, Black}}, Joined -> True, GOLabelShift -> {0.0, 0.0}, PlotLabel -> "Band structure", FrameLabel -> {" ", "a.u."}}
 Options[GTBandStructure]            = {GOStore->0,GOTbOrthogonal->True,GOPhotonic->False,Joined->False,PlotRange->All,FrameLabel->{" ","Energy"},PlotLabel->"Band structure",GOShift->0,LabelStyle->{},PlotStyle->Blue,GOIrepTextStyle->{20,Darker[Red]},GOFermiEnergy->-100000}
 Options[GTCompatibility]            = {GOVerbose->False,GOIrepNotation->"Bethe",GOFast->GOFastValue}
 Options[GTDensityOfStates]          = {GOBands->"ALL", GOPhotonic -> False, GOStore -> 0,PlotRange -> All,FrameLabel -> {"Energy","DOS"},PlotLabel -> "Density of States", GOPlotDos -> "DOS",PlotStyle->Blue,GOFermiEnergy->-100000}
 Options[GTDensityOfStatesPlot]      = {GOBands->"ALL",GOStore->0,PlotRange->All,FrameLabel->{"Energy","DOS"},PlotLabel->"Density of States",GOFermiEnergy->-100000,PlotStyle->Blue,GOPlotDos -> "DOS"}
 Options[GTDensityOfStatesRS]        = {GOStore -> 0, PlotRange -> All, FrameLabel -> {"Energy", "DOS"}, GOFermiEnergy -> -100000, PlotLabel -> "Density of States", PlotStyle -> Red, GOPlotDos -> "DOS"}
 Options[GTFermiSurface]             = {GOVerbose -> False, GOBZPath -> False, GOBZ -> True, GORegionFunction -> False,VertexLabelStyle->Directive[Black, 20,Background -> Yellow]}
 Options[GTFermiSurfaceCut]          = {GOVerbose -> False,GOPlot->True,ContourStyle->{}}
 Options[GTFatBandsPlot]             = {PlotRange -> All, Joined -> False, FrameLabel -> {" ", "Energy"}, PlotLabel -> "Band structure", PlotStyle -> Blue, GOFermiEnergy -> -100000, GOFatBands -> 0}
 Options[GTPartialDOS]               = {GOStore -> 0, PlotRange -> All, FrameLabel -> {"Energy", "DOS"}, PlotLabel -> "Density of States", GOPlotDos -> "DOS", PlotStyle -> Blue, GOFermiEnergy -> -100000}
 Options[GTBandsDOSPlot]             = {PlotRange -> All, Filling->5,FrameLabel -> {{" ", "Energy"}}, PlotLabel -> {"Bands", "DOS / IDOS", "partial DOS"}, GOShift -> 0, PlotStyle -> Blue, GOFermiEnergy -> 0, GOVerbose -> False, GOPlotDos -> "DOS", GOPlotGrid -> {800, 400, 40}}


Begin["`Private`"] (* Begin Private Context *) 

(****f* /GTBandsDOSPlot
! NAME
!  GTBandsDOSPlot
! AUTHOR
!  W. Hergert 
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 30.05.2020 : first complete version, also documentation page finished
!  * 21.06.2020 : final check of implementation
! USAGE
!  GTBandsDOSPlotlDOS[bands,dos,nbands] creates a combined plot of bandstructure and density of states
! INPUT
!  * bands  - band strucuture
!  * dos    - density of states 
!  * nbands - number of bands to be plotted from bands
!
! OUTPUT
!  a combined plot prepared by means of plotGrid
!   
! GTPack OPTIONS
!
!   * GOPlotDos:
!
!           "DOS"   - plot DOS (Standard) 
!           "IDOS"  - plt IDOS
!           "PDOS"  - plt PDOS
!           "ALL"   - plot DOS/IDOS and PDOS
!   * GOFermiEnergy:
!
!           -10000  - Fermi energy not calculated (standard) 
!          {nel,es} - nel - number of electrons, es - energy to start the search
!               nel - number of electrons, input for FindRoot set automatically
!   * GOVerbose:
!
!			"True"  - additional information
!			"False" - no information (standard)
!
!   * GOPlotGrid: this an option containing parameters for plotGrid (ThirdParty.m) (width, height, padding)=(800,400,40)
!             
! STANDARD OPTIONS
!  * PlotRange     - standard : All
!  * FrameLabel    - standard: {"Energy","DOS"}
!  * PlotLabel     - standard: "Density of states"
!  * PlotStyle     - standard: Blue
! 
! GTPack MODULES
!   plotGrid, GTBandsPlot
!  
! GTPack NOTEBOOKS
!  -  
! DESCRIPTION
!  I thought it might be a good idea to have the opportunity to create plots of a band structure and the corresponding
!  DOS in manner like it can be found in literatur, i.e. the band structure and the DOS aside, turned by 90 degree. This 
!  might be better in some cases instead having two different plots.
!  Finally it is implemented in a quite general way, i.e. if GTPlotDos is "DOS" or "IDOS" or "PDOS" we have the bandstructure and 
!  one DOOS plot. In case of "ALL" we have the band structure plot and two DOS plots.
! LITERATURE
!  -
! TODO
!  check the module. seems to be space for improvements
! RELEASE
!  1.2.
! PROBLEMS
!  The data have to be "rotated" for ListPlot. This means that the axis have to be constructed by hand. The corresponding information 
!  is provide by PlotLabel.
!  For "IDOS" and "ALL" the integrated DOS has to be plotted together with the DOS. The inttegrated DOS is scaled and a second axis has 
!  to be constructed. 
!  It might be tht those constructions are not stable enough.
!
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTBandsDOSPlot::plot   =  "Wrong Pattern in PlotRange";
GTBandsDOSPlot::bnddat =  "Band bandstructure data. Perhaps no k-path.";
GTBandsDOSPlot::dosp   =  "Selector for kind of plot is `1`. It has to be DOS, IDOS, PDOS or ALL.";
GTBandsDOSPlot::dos    =  "Error in Dpreparation of plots.";

GTBandsDOSPlot[bands_, dos_, nbdplt_, OptionsPattern[]] := Module[
  {frame,pltl,plts,shft,efl,verb,width,height,padd,dosp,pltr,kp,energ,pdos,tick,plband,pldos,dosm,dim,tmp1,tmp2,
   dosn,elmax,daxis,i,gbands,dos1,dplot,plist,dos2,npart,dpplot,fill},
(*--- options ---*)
  frame = OptionValue[FrameLabel];
  pltl  = OptionValue[PlotLabel];
  plts  = OptionValue[PlotStyle];
  fill   = OptionValue[Filling];
  shft  = OptionValue[GOShift];
  efl   = OptionValue[GOFermiEnergy];
  verb  = OptionValue[GOVerbose];
  {width, height, padd} = OptionValue[GOPlotGrid];
(*--- evaluate options ---*)
  dosp  = OptionValue[GOPlotDos];
  If[Intersection[{dosp}, {"DOS", "IDOS", "PDOS", "ALL"}]=={},
  	 Message[GTBandsDOSPlot::dosp,dosp]; Return[],
  	 None
  ];  
  pltr  = OptionValue[PlotRange];
  If[Length[pltr] == 4,
   {kp, energ, pdos, tick} = pltr;
   plband = {kp, energ};
   pldos  = {pdos, energ};
   dosm   = pdos[[2]],
     Message[GTBandsDOSPlot::plot]; Return[]
   ];
   If[Length[bands] != 2,
     Message[GTBandsDOSPlot::bnddat]; Return[],
     None
   ];  
(*--- band structure plot, use standard Module GTBandsPlot ---*)
   gbands = GTBandsPlot[bands, nbdplt, PlotRange -> plband, GOFermiEnergy -> efl, Joined -> True, GOShift -> shft, 
            PlotLabel -> pltl[[1]], FrameLabel -> frame[[1]]];
(*--- prepare data for DOS plots ---*)
  dim = Dimensions[dos];
  If[Length[dim] > 2,
(*--- DOS and IDOS from GTDensityOfStates  adapt data set---*)   
     tmp1 = dos[[1]] // Transpose; 
     tmp2 = dos[[2]] // Transpose;
     dosn = Append[tmp1, tmp2[[2]]] // Transpose,
(*--- PDOS data set no adaption necessary ---*)
     dosn = dos
  ];
(*--- scale IDOS if necessary ---*)
  If[dosp == "IDOS"||dosp=="ALL",
     tmp1 = dosn // Transpose; 
     elmax = Max[tmp1[[3]]];
     If[verb,
        Print["Maximum Number of electrons= ", elmax],
        None
     ];
     tmp1[[3]] = tmp1[[3]]*dosm/elmax/2.  ;
     dosn      = tmp1 // Transpose;
   (*--- axis electron numbers for DOS plot ---*)
     daxis     = Table[{dosm/elmax/2. (i - 1), i-1}, {i, 1, 2 IntegerPart[elmax]}],
     None
  ];
(*  Print[dosn];*)
(*--- data for plot ---*)
  Switch[dosp,
  "DOS",
     dos1  = Map[{#[[2]], #[[1]]} &, dosn];
     dplot = ListPlot[dos1, Joined -> True, PlotRange -> pldos, Frame -> True, PlotLabel -> pltl[[2]],Filling -> fill, FrameLabel -> frame[[2]], 
                      FrameTicks -> {{True, Automatic}, {{#, ToString[Abs[#]]} & /@ FindDivisions[{0, dosm}, tick], False}}]; 
     plist = {{gbands, dplot}},
  "IDOS",
     dos1  = {Map[{#[[2]], #[[1]]} &, dosn], Map[{#[[3]], #[[1]]} &, dosn]}; 
     dplot = ListPlot[dos1, Joined -> True, PlotRange -> pldos, PlotLabel -> pltl[[2]],Frame -> True, Filling -> fill, FrameLabel -> frame[[2]], 
     	              FrameTicks -> {{True, False}, {{#, ToString[Abs[#]]} & /@ FindDivisions[{0, dosm}, tick], daxis}}]; 
     plist = {{gbands, dplot}},
  "PDOS",
      dos2  = {}; 
      npart = Length[dosn[[1]]];
      Do[
      	 dos2 = Append[dos2, Map[{-#[[i]], #[[1]]} &, dosn]]
      , {i, 4, npart}];
      dplot = ListPlot[dos2, Joined -> True, PlotRange -> {{-dosm, 0}, energ}, PlotLabel -> pltl[[2]],Frame -> True, Filling -> fill, FrameLabel -> frame[[2]], 
                       FrameTicks -> {{True, False}, {{#, ToString[Abs[#]]} & /@ FindDivisions[{-dosm, 0}, tick], False}}]; 
      plist = {{dplot, gbands}},
  "ALL",
      dos1  = {Map[{#[[2]], #[[1]]} &, dosn], Map[{#[[3]], #[[1]]} &, dosn]};
      dos2  = {}; 
      npart = Length[dosn[[1]]];
      Do[
      	 dos2 = Append[dos2, Map[{-#[[i]], #[[1]]} &, dosn]]
      , {i, 4, npart}];
      dpplot = ListPlot[dos2, Joined -> True, PlotRange -> {{-dosm, 0}, energ},PlotLabel -> pltl[[3]], Frame -> True, Filling -> fill, 
                        FrameLabel -> frame[[2]], FrameTicks -> {{True, False}, {{#, ToString[Abs[#]]} & /@ FindDivisions[{-dosm, 0}, tick], False}}];
    (*  dplot  = ListPlot[dos1, Joined -> True, PlotRange -> pldos, PlotLabel -> pltl[[2]],Frame -> True, Filling -> 5, FrameLabel -> frame[[2]], 
                        FrameTicks -> {{True, False}, {{#, ToString[Abs[#]]} & /@ FindDivisions[{0, dosm}, tick], daxis}}]; *)
                        dplot = ListPlot[dos1, Joined -> True, PlotRange -> pldos, PlotLabel -> pltl[[2]],Frame -> True, Filling -> fill, FrameLabel -> frame[[2]], 
     	              FrameTicks -> {{True, False}, {{#, ToString[Abs[#]]} & /@ FindDivisions[{0, dosm}, tick], daxis}}]; 
      plist = {{dpplot, gbands, dplot}},
  _,
      Message[GTBandsDOSPlot::dos]; Return[]
   ];
  (*--- plot with plotGrid ---*)
  
  plotGrid[plist, width, height, ImagePadding -> padd]
  ]



(*
***)


(****f* /GTPartialDOS
! NAME
!  GTPartialDOS
! AUTHOR
!  W. Hergert 
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 20.09.2018 : first complete version
! USAGE
!  GTPartialDOS[hop,struc,parm] calculates density of states in TB models
! INPUT
!  * hop   - tight-binding Hamiltonian
!  * struc - structure (fcc,bcc,...)
!  * parm  - list of parameters for DOS calculation:
!
!            1 -  n    - index k-mesh
!            2 -  a    - scaling factor k-mesh
!            3 - emin  - minimum energy
!            4 - emax  - maximum enery
!            5 - ne    - number of energy points  
!            6 - sigma - width of the Gaussians (if = o then equal to energy steps)
!            7 - scale - scaling factor for the DOS
!            8 - pDOS  - control of partial DOS calculation list like {{"p",{2,3,4}}
! OUTPUT
!  graph, or multiple graphs
!   
! GTPack OPTIONS
!
!   * GOPlotDos:
!
!           "DOS"   - plot DOS (Standard) 
!           "IDOS"  - plt IDOS
!           "PDOS"  - plt PDOS
!           "ALL"   - plot DOS/IDOS and PDOS
!   * GOFermiEnergy:
!
!           -10000  - Fermi energy not calculated (standard) 
!          {nel,es} - nel - number of electrons, es - energy to start the search
!   * GOStore       - DOS will be stored, if a filename is given   
! STANDARD OPTIONS
!  * PlotRange     - standard : All
!  * FrameLabel    - standard: {"Energy","DOS"}
!  * PlotLabel     - standard: "Density of states"
!  * PlotStyle     - standard: Blue
! 
! GTPack MODULES
!   GTBZPointMesh, BZGauss, GTWriteToFile, TwoAxisListPlot
!  
! GTPack NOTEBOOKS
!  -  
! DESCRIPTION
!  The priciple scheme of the 
! LITERATURE
!  -
! TODO
!  implementation of the true partial DOS
! RELEASE
!  1.0.1
! PROBLEMS
!  GTPartialDOS has the same philosophy like GTDensityOfStates. All is calculated from the Hamiltonian.
!  In pinciple, one could also add something like GTPartiaDOSPlot that reads the bands calculated by 
!  GTBands to calculate only the DOS. Seems to be not so urgent.
!
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTDensityOfStates::dosp    =  "Selector for kind of plot is `1`. It has to be DOS, IDOS, PDOS or ALL.";
GTDensityOfStates::plots   =  "Wrong Pattern in Plotstyle";
GTDensityOfStates::idos    =  "Use IDOS or ALL in GoPlotDos to calculate Fermi energy.";
GTDensityOfStates::bnds    =  "Wrong descriptor in GOBands.";

GTPartialDOS::badarg = "Selector for kind of plot is `1`. It has to be DOS, IDOS, PDOS or ALL.";
GTPartialDOS::badsel = "Selector for partial density of states not defined";
GTPartialDOS::badcol = "Wrong description of colors in PlotStyle";
GTPartialDOS::badef  = "Use IDOS or ALL in GOPlotDos to calculate Fermi energy.";
GTPartialDOS::frame  = "Wrong FrameLabel for ALL in GOPlotDos.";
GTPartialDOS::plotl  = "Wrong PlotLabel for ALL in GOPlotDos.";

GTPartialDOS[hop_, struc_, parm_, OptionsPattern[]] := 
   Module[{data,plr,plt,frame,frame1,pltl,plt1,plt2,style,fen,test,names,ranges,cdos,cidos,flog,n,a,qp,nkp,de,emin,emax,ne,scale,axis,ndos,npdos,
   	       dos,nbands,ham,k,ew,ev,fill,tmp,evec,fac,nr,j,i,m,bs,func,x,xs,nel,id,rule,xf,gf,dosdat,pldos,idosdat,pdosdat,plidos,ppdos,pall},
  	flog  = False;  	      
  (*---options and consistency tests---*)
    data  = OptionValue[GOStore];
    plr   = OptionValue[PlotRange];
    plt   = OptionValue[GOPlotDos];
    test  = Intersection[{plt}, {"DOS", "IDOS", "ALL", "PDOS"}];
    If[test == {},
       Message[GTPartialDOS::badarg, plt]; Return[],
       None
    ];
    frame = OptionValue[FrameLabel];
    pltl  = OptionValue[PlotLabel];
    frame1 = frame;
    plt1 = plt2 = pltl;
    If[plt == "ALL",
       If[Length[Dimensions[frame]]==2,
          frame1 = {frame[[2, 1]], frame[[1, 1]]},
          Message[GTPartialDOS::frame, plt]; Return[]
       ];
       If[Head[pltl] === List,
          plt1 = pltl[[1]];
          plt2 = pltl[[2]],
          Message[GTPartialDOS::plotl, plt]; Return[]
       ];
    ];
    style = OptionValue[PlotStyle];
    If[Head[style] === RGBColor,
       cdos = cidos = style,
       If[Length[style] == 2,
            cdos = style[[1]];   
            cidos = style[[2]],
            Message[GTPartialDOS::badcol]; Return[]
       ]
    ];
    If[plt == "PDOS" ||plt=="ALL",
    	  If[Length[parm] == 7,
         Message[GTPartialDOS::badsel]; Return[],
         None
    	  ];
      {names, ranges} = Transpose[parm[[8]]];
       names          = Prepend[names, "DOS"]
    ];  
    fen   = OptionValue[GOFermiEnergy];    
    If[fen > 0||Head[fen]===List,
    	  If[{}==Intersection[{plt},{"IDOS","ALL"}],
    	  	 Message[GTPartialDOS::badef];Return[],
    	  	 None
    	  ]
    ];   	  	
  (*---k-mesh---*)
    n   = parm[[1]];
    a   = parm[[2]];
    qp  = GTBZPointMesh[n, a, struc];
    nkp = Length[qp];
  (*---prepare for DOS---*)
    If[parm[[6]] == 0,
        \[Sigma] = de,
        \[Sigma] = parm[[6]]
    ];
    emin  = parm[[3]];
    emax  = parm[[4]];
    ne    = parm[[5]];
    scale = parm[[7]];
    de    = (emax - emin)/(ne - 1.);
    axis  = Table[(i - 1)*de + emin, {i, 1, ne}];
    If[plt == "PDOS" || plt == "ALL",
      ndos  = Length[parm[[8]]] + 3;
      npdos = Length[parm[[8]]];
      fill={1->None};
    	  Do[
    	     fill = Append[fill, {i + 1 -> Axis}]
    	  , {i, 1, npdos}];
    	  fill=fill//Flatten,
      ndos = 2
    ];
    dos    = Table[0, {ndos}, {ne}];
    nbands = Length[hop];
  (*--- loop over k-values ---*)
    Do[
       ham = hop /. {\[Xi] -> qp[[k, 1]], \[Eta] -> qp[[k, 2]], \[Zeta] -> qp[[k, 3]]};
  (*---Calculate eigenvalues and eigenvectors---*)
       {ew, ev} = Eigensystem[ham] // Sort;
        ew      = ew // Chop;
       Do[
          tmp      = scale*Map[BZGauss[#, \[Sigma], ew[[i]]] &, axis];
          dos[[1]] = dos[[1]] + tmp
       , {i, 1, nbands}] ;
   (*--- partial DOS ---*)
       If[plt == "PDOS" || plt == "ALL",
          Do[
    	         evec = ev[[j]];
             Do[
             	fac=0;
             	Do[
             	   nr=ranges[[i,m]];	 
             	   fac=fac+evec[[nr]]*Conjugate[	evec[[nr]]]
             	,{m,1,Length[ranges[[i]]]}];
             	fac=fac//Chop;
                tmp = scale*fac*Map[BZGauss[#, \[Sigma], ew[[j]]] &, axis];
                dos[[i + 3]] = dos[[i + 3]] + tmp/nkp;
         dos[[3]]     = dos[[3]] + tmp/nkp
                
                (* nur part bands aufsummiert! *)
             , {i, 1 npdos}]
          , {j, 1, nbands}],
          None
       ]
    , {k, 1, nkp}];
    dos[[1]] = dos[[1]]/nkp;
    bs       = Transpose[{axis, dos[[1]]}];
  (*---Calculate IDOS---*)
    If[plt == "IDOS" || plt == "ALL",
       func = Interpolation[bs];
       Do[
          dos[[2, i]] = Integrate[func[x], {x, emin, axis[[i]]}]
       , {i, 1, ne}];
     Print["integrated DOS at E= ",axis[[ne]]," is ", dos[[2,ne]]];
  (*---calculation of Fermi energy---*)
    If[Head[fen] === List,
       flog = True;  
       xs   = fen[[2]];  
       nel  = fen[[1]],
       If[fen > 0,
          flog = True; 
          xs = (axis[[ne]] + axis[[1]])/2.; 
          nel = fen,
          flog = False
        ]
     ];
     If[flog, 
     	id=Transpose[{axis, dos[[2]]}];
     	func =Interpolation[id];
        rule = FindRoot[func[x] == nel, {x, xs}];
        xf = x /. rule;
        Print["Fermi Energy \!\(\*SubscriptBox[\(E\), \(F\)]\)= ", xf, 
             " (units see band structure)"
            ];
        gf = Graphics[{Red, Thick,Line[{{xf, -1}, {xf, Max[dos[[1]]] + 10}}]}],
       None
       ], 
       None
    ];
  (*---  store data ---*) 
    If[data == 0,
       None,
        tmp = Prepend[dos, axis] // Transpose;
       GTWriteToFile[tmp, data]
   ];
 (*--- plot the data ---*)
   dosdat = Transpose[{axis, dos[[1]]}];
  Switch[plt,
   "DOS",  
     pldos = ListPlot[dosdat, Frame -> True, PlotStyle -> cdos,FrameLabel -> frame1, Joined -> True, 
     	     PlotRange -> plr,PlotLabel -> plt1];,
   "IDOS",
     idosdat = Transpose[{axis, dos[[2]]}];
     plidos  = ListPlot[idosdat, Frame -> True, PlotStyle -> cdos, FrameLabel -> frame1, Joined -> True, 
     	       PlotRange -> plr, PlotLabel -> plt1];,
   "PDOS",
     pdosdat = {dosdat};
     Do[
     	pdosdat = Append[pdosdat, Transpose[{axis, dos[[i + 3]]}]]
     , {i, 1, npdos}];    
     ppdos = ListPlot[pdosdat, Frame -> True, FrameLabel -> frame1, Joined -> True, PlotRange -> plr, 
     	     PlotLabel -> plt2, Filling->fill,PlotLegends -> names];,
   "ALL",
     pdosdat = {dosdat};
     Do[
     	pdosdat = Append[pdosdat, Transpose[{axis, dos[[i + 3]]}]]
     , {i, 1, npdos}];
     ppdos   = ListPlot[pdosdat, Frame -> True, FrameLabel -> frame1,Joined -> True, PlotRange -> plr, 
     	       PlotLabel -> plt2, PlotLegends -> names,Filling->fill];
     idosdat = Transpose[{axis, dos[[2]]}];
     pall    = TwoAxisListPlot[{dosdat, idosdat}, Frame -> True,PlotStyle -> {cdos, cidos}, FrameLabel -> frame, 
               PlotLabel -> plt1, Joined -> True];,
   _,
   Abort[]
   ];
  Switch[plt,
   "DOS",
    If[flog,
   	   Return[Show[pldos,gf]],
   	   Return[Show[pldos]]
    ],
   "IDOS",
    If[flog,
   	   Return[Show[plidos,gf]],
   	   Return[Show[plidos]]
    ],
   "PDOS",
    If[flog,
   	   Return[Show[ppdos,gf]],
       Return[Show[ppdos]]
    ],
   "ALL",
   If[flog,
      Return[Grid[{{Show[pall,gf],Show[ppdos,gf]}}]],
      Return[Grid[{{Show[pall],Show[ppdos]}}]]
    ],   
   _,
   Abort[]
   ];
  ]

(*
***)

(****f* /GTFatbandsPlot
! NAME
!  GTFatBandsPlot
! AUTHOR
!  W. Hergert 
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!   * 20.09.2018 : first complete version
! USAGE
!  GTFatBanndsPlot[bands,bandnumb,scale] gives a bands structure plot with indication of the angular momentum character (fsat bands).
! INPUT
!  * bands    - band structure calculated by means of GTBands
!  * bandnumb - number of bands used in the plot
!  * scale    - scaling factor to control the width of the fat bands
! OUTPUT
!  plot of band structure with fat bands
! GTPack OPTIONS
!    * GOFermiEnergy:  plot of
!
!            0      - no Fermi energy givem
!           .NE. 0  - plot Fermi energy
!          
!   * GOFatBands:
!   
!           0             - fat bands not defined, normal band structure plot
!           {"p",{2,3,4}} - angular momentum and list of orbitals in the Hamiltonian
! Standard OPTIONS
!  * PlotRange     - standard : All
!  * FrameLabel    - standard: {"Energy","DOS"}
!  * PlotLabel     - standard: "Band structure"
!  * PlotStyle     - standard: Blue (Plotstyle for fat bands, others are fixed)
!  * Joined        - standard: False
!  PlotRange
! GTPack MODULES
!  -
! DESCRIPTION
!   At the end of the day it was simple to realize. From the eigenvector at the corresponding k-point
!   the weight factor for the fat band is calculated. Than half of the weight is subtracted and added to 
!   the norma band structure, giving lower and upper limit of the fat band. scale is finally used to 
!   scale the fat band width. The normal band structure is shown together with the fat band, a filled region 
!   between the upper and lower limit of the fat band.
! LITERATURE
!  CrystallographicPointGroups-author.nb from Wolfram Demonstration project
! TODO
!  Find out more details about the group theory implementation in Mathemaatica. 
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTFatBandsPlot::badbnds =  "Warning: No eigenvectors are given, normal band structure plot.";
GTFatBandsPlot::badfat  =  "Warning: Fat bands not defined, normal band structure plot.";
GTFatBandsPlot::badnumb =  "Warning: Not enough bands `1` bands used.";

GTFatBandsPlot[bands_, nbdplt_, scale_, OptionsPattern[]] := 
 Module[{bds,lx,ev,fb,npt,nl,join,plr,frame,pltl,plts,ef,fat,axis,evs,i,nbnd,mi,ma,bs,xmax,bar,lab,band,
 	     k,gb,gs,bandu,bando,fac,evec,n,gf,efl,eu,eo,tmp,nr,nfat},
   (*--- check input of bands ---*)
     If[Length[bands[[1]]] == 2,
        bds = bands[[1, 1]] // N;
        lx  = bands[[1, 2]];
        ev  = bands[[2]];
        fb  = True,
        bds = bands[[1]] // N;
        lx  = bands[[2]];
        fb  = False;
        Message[GTFatBandsPlot::badbnds]
     ];     
     npt   = Length[bds];
     nl    = Length[lx];
  (*--- options ---*)
     join  = OptionValue[Joined];
     plr   = OptionValue[PlotRange];
     frame = OptionValue[FrameLabel];
     pltl  = OptionValue[PlotLabel];
     plts  = OptionValue[PlotStyle];
     ef    = OptionValue[GOFermiEnergy];
     fat   = OptionValue[GOFatBands];
     If[fb == True && Head[fat] === List,
        nfat=Length[fat[[2]]],
        Message[GTFatBandsPlot::badfat];
        fb=False
     ];
     axis = Table[0, {i, 1, npt}]; 
     evs  = Table[0, {i, 1, npt}];
     nbnd = Length[bds[[1, 4]]];
     If[nbnd >= nbdplt,
        nbnd = nbdplt,
        Message[GTFatBandsPlot::badnumb, ToString[nbnd]]
     ];
     (*--- normal band structure, scaling plot, labels ---*)
     Do[
        axis[[i]] = bds[[i, 2]];
        evs[[i]] = Sort[bds[[i, 4]], Less]
     , {i, 1, npt}];
     mi   = Min[evs]; 
     ma   = Max[evs]; 
     bs   = {};
     xmax = Max[axis] // N;
     Print["Maximum Abscissa = ", xmax];
     bar  = Table[{{axis[[lx[[i, 1]]]], mi - 10}, {axis[[lx[[i, 1]]]], ma + 10}}, {i, 1, nl}];
     lab  = Table[{axis[[lx[[i, 1]]]], lx[[i, 2]]}, {i, 1, nl}];
     Do[
     	band = {};
        Do[
           band = Append[band, {axis[[i]], evs[[i, k]]}]
        , {i, 1, npt}];
        bs = Append[bs, band]
     , {k, 1, nbnd}];
 (*--- plots normal band structure, bars and labels ---*) 
     gb = ListPlot[bs, Frame -> True, FrameTicks -> {{Automatic, None}, {lab, lab}}, PlotRange -> plr, 
          Joined -> join, FrameLabel -> frame, PlotLabel -> pltl, PlotStyle -> {{Opacity[0.3], Red}}];
     gs = ListLinePlot[bar, PlotStyle -> Black, Ticks -> {None, Automatic}, PlotRange -> plr];
 (*--- construction of fat bands ---*)    
     If[fb == True && Head[fat] === List,
        bs = {};
        Do[
           bandu = {}; 
           bando = {};
           Do[
              evec  = ev[[i, 2, k]];
              fac=0;
              Do[
              	 nr=fat[[2,i]];
              	 fac=fac+evec[[nr]]*Conjugate[evec[[nr]]]
              ,{i,1,nfat}];
              fac=fac//Chop;	
            (*  fac   = evec[[evu ;; evo]].Conjugate[evec[[evu ;; evo]]] // Chop; *)
              eu    = evs[[i, k]] - scale*fac/2;
              eo    = evs[[i, k]] + scale*fac/2;
              If[eu<eo,
              	 None,
              	 tmp=eu;eu=eo;eo=tmp
              ];
              bandu = Append[bandu, {axis[[i]],eu}];
              bando = Append[bando, {axis[[i]],eo}]
           , {i, 1, npt}];
           bs = Append[bs, bandu];
           bs = Append[bs, bando];
        , {k, 1, nbnd}];
        gf = ListPlot[bs, Frame -> True, FrameTicks -> {{Automatic, None}, {lab, lab}}, PlotRange -> plr, 
        Joined -> join, FrameLabel -> frame, PlotLabel -> pltl, 
        Filling -> Table[n -> {n + 1}, {n, 1, 2*nbnd - 2, 2}], 
        PlotStyle -> {{Opacity[0.5],plts}}];
     ];
    If[Abs[ef] > 0, 
    	   efl = Graphics[{Green, Line[{{0, ef}, {xmax, ef}}]}];
       If[fb, 
          Show[gb,gf,gs,efl],
          Show[gb,gs,efl]
       ],
       If[fb, 
          Show[gb,gf,gs],
          Show[gb,gs]
       ]
    ]              
  ]

(*
***)


(****f* /GTSortEigensystem
! NAME
!  GTSortEigensystem
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 08.02.2017 : first implemenetation
!  * 28.12.2017 : check header and documentation
!  * 23.06.2017 : check header and documentation
! USAGE
!  GTSortEigensystem[eigensystem_] sorts the eigensystem in canonical order for real eigenvalues.
! INPUT
!  eigensystem  - output of Eigensystem
! OUTPUT
!  eigensystem sorted with respect to the eigenvalues
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  -
! GTPack Modules
!   -
! GTPack NOTEBOOKS
!  - 
! DESCRIPTION 
!  Eigensystem gives the eigenvalues sorted with respect to the ABSOLUTE value. 
!  We need the eigenvalues in canonical order. This module is doing the job.  
!
! LITERATURE
!   from Mathematica StackExchange
! TODO
!  
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

  GTSortEigensystem[es_] := Module[
  	{eigs, vecs, l1, l2},
     {eigs, vecs} = es;
      l1 = Partition[Riffle[eigs, vecs], 2];
      l2 = Sort[l1, #1[[1]] < #2[[1]] &] // Transpose;
  Return[l2]
  ]

(*
***)

(****f* /GTBands
! NAME
!  GTBands
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 15.11.2013 : first version
!  * 07.02.2015 : check options, GOplt renamed to GOPlotBands
!  * 23.05.2016 : implementation of nonorthogonal problems, new output if 
!				 GOVerbose=True, Print with Grid as already used in other commands 
!  * 08.11.2016 : output in table restricted to 15 bands   
!  * 30.12.2016 : sort the  energies in the table
!  * 28.12.2017 : check header and documentation
!  * 20.12.2018 : check of header and docu for version 1.0
!  * 23.06.2018 : check of header and documentation
!  * 14.06.2020 :new message system, output extented to k-path
! USAGE
!  GTBands[Hamiltionian,kpoints,nev] calculates the lowest nev band energies for a given Hamiltonian at certain kpoints.
! INPUT
!  * Hamiltionian - "Hamiltonian", for a non-orthogonal calculation a list of two matrices
!  * kpoints	  - complete list of k-points, as generated with GTBZLines for band structures
!                   or GTBZmesh for DOS calcualtions     
!  * nev          - number of bands
! 
! OUTPUT
!  Output of the calculated band structure in the form controlled by the options.
! 
! GTPack OPTIONS
!   * GOEigenvectors - calculation eigenvectors  (standard - no eigenvectors)
!   * GOTbOrthogonal - non-orthogonal problem (standard - orthogonal)
!   * GOPhotonic     - photonic Hamiltonian ( standard -electronic structure)                      
!   * GOStore        - output of the calculation to a file ( standard -  no output)  
!   * GOPlotBands    - gives information for the band structure plot, or the bare eigenvalues (standard - band structure plot) 
!   * GOShift        - shift of all energies by a certain amount
!   * GODecimals     - decimals for output in table (standard is 4)
!   * GOVerbose:
!
!	   - False : no additional output of information  (standard)   
!	   - True  : energies at high symmetry points in a table
!
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTWriteToFile
! GTPack NOTEBOOKS
!  -
! DESCRIPTION 
!  Especially for 3D calculations of photonic crystals the combination with the plot is not
!  very efficient because of the long time to calculate the band structure.
!
!  In principle an orthogonal and a non-orthogonal scheme can be used. In the non-orthogonal 
!  scheme a list consisting of Hamiltonian and overlap matrix has to be provided.
!  
!  In the photonic case the square root of the eigenvalues corresponds to the eigenfrequencies.
!  GOPhotonic is used to distinguish the cases. If phonons are considered also this option has to be used.
!
!  If band structure and eigenvectors should be stored a list of two file names is expected.
!  
! LITERATURE
! 
! TODO
!  The nonorthogonal scheme is not fully tested  yet. One has
!  to test if this has some influence on the calculation time. Is it necessary to
!  distinguish between calculation with and without the calculation of eigenvectors?
!
! RELEASE
!  1.0.0
! PROBLEMS
!  Test nonorthogonal problems before the firstofficial release.
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTBands::badarg  = "If GOTbOrthogonal->`1`  {Hamiltonian, Overlap} necessary";
GTBands::badkval = "GOTPlotBands not in agreement with k-lines/k-mesh";
GTBands::bands   = "Warning: Output restricted to 15 bands";

GTBands[hopin_, kpoints_, nev_, OptionsPattern[]] := 
  Module[{nop,nwv,ndt,nvb,phot,nplt,shft,hop,ov,kpts,npt,bands,wav,kpa,i,hop1,ov1,bd,tab,ind,lx,k,nbds,tmp,nround},
  (*---Options---*)
  nop    = OptionValue[GOTbOrthogonal];
  nwv    = OptionValue[GOEigenvectors];
  ndt    = OptionValue[GOStore];
  nvb    = OptionValue[GOVerbose];
  nround = OptionValue[GODecimals];
  phot   = OptionValue[GOPhotonic];
  nplt   = OptionValue[GOPlotBands];
  shft   = OptionValue[GOShift];
  (*---check form of Hamiltonian---*)
  If[nop,
     hop = hopin,
     If[Length[hopin] != 2,
        Message[GTBands::badarg, nop]; Return[],
        hop = hopin[[1]]; ov = hopin[[2]]
     ]
  ];
  (*--- prepare k-points ---*)
  If[nplt,
  	If[Length[kpoints] == 2,
      None,
      Message[GTBands::badkval]; Return[]
    ];
  (*--- kpts is the bare list of k-points,  lx the list of symmetry points ---*)
    kpts = Transpose[kpoints[[1]]][[3]]; lx = kpoints[[2]],
  (*--- mesh data for DOS calculaton ---*)
    kpts = kpoints
  ];
  npt = Length[kpts]; bands = {}; wav = {};
  If[nwv,
  (*--- Calculation of eigenvalues and eigenvectors ---*)          
     Do[
        kpa  = kpts[[i]];
        hop1 = hop /. {\[Xi] -> kpa[[1]], \[Eta] -> kpa[[2]], \[Zeta] -> kpa[[3]]} // Chop;
        If[nop,  
           bd  = Eigensystem[hop1] // Chop;
           bd  = GTSortEigensystem[bd],
           ov1 = ov /. {\[Xi] -> kpa[[1]], \[Eta] -> kpa[[2]], \[Zeta] -> kpa[[3]]} // Chop;
           bd = Eigensystem[{hop1, ov1}] // Chop;
           bd  = GTSortEigensystem[bd]
        ];
        If[phot,
           bd[[1]] = Sqrt[bd[[1]]],
           None
        ];
        bd[[1]] = bd[[1]] - shft // Chop;
        If[nplt,          
           bands = Append[bands, {kpoints[[1, i, 1]], kpoints[[1, i, 2]], kpa,bd[[1,1;;nev]]}],
           bands = Append[bands, {kpa, bd[[1,1;;nev]]}]
        ];
        wav = Append[wav, {kpa, bd[[2,1;;nev]]}]
     , {i, 1, npt}],
  (*--- Calculation of eigenvalues only ---*)
     Do[
     	kpa = kpts[[i]];
        hop1 = hop /. {\[Xi] -> kpa[[1]], \[Eta] -> kpa[[2]], \[Zeta] -> kpa[[3]]} // Chop;
        If[nop, 
           bd  = Eigenvalues[hop1] // Chop;
           bd  = Sort[bd,#1<#2&];
           bd  = bd[[1;;nev]],  
           ov1 = ov /. {\[Xi] -> kpa[[1]], \[Eta] -> kpa[[2]], \[Zeta] -> kpa[[3]]} // Chop;       
           bd  = Eigenvalues[{hop1, ov1}] // Chop;
           bd  = Sort[bd,#1<#2&];
           bd  = bd[[1;;nev]]
        ];
        If[phot,
           bd = Sqrt[bd],
           None
        ];   
        bd = bd - shft // Chop;
        If[nplt,
           bands = Append[bands, {kpoints[[1, i, 1]], kpoints[[1, i, 2]], kpa, bd}],
           bands = Append[bands, {kpa, bd}]]
      , {i, 1, npt}]
  ];     
  (*--- Print list of eigenvalues ---*)
  If[nvb,
  	If[nev>15,
  	   nbds=15;Message[GTBands::bands],
  	   nbds=nev
  	 ];
     If[nplt,
        npt = Length[kpoints[[2]]];
        tab = Table[Table[0, {nbds + 2}], {npt}];
        Do[
           ind         = lx[[i, 1]];
           tab[[i, 1]] = lx[[i, 2]];
           tab[[i, 2]] = bands[[ind, 3]];
           Do[
           	  tmp=bands[[ind, 4, k]];
           	  If[nround==0,
                 tab[[i, k + 2]] = tmp,
                 tab[[i, k + 2]] = ToString[PaddedForm[N[tmp], {nround, nround}]]
              ];
           , {k, 1, nbds}]
        , {i, 1, npt}];
        Print[
        	  Grid[Join[{{Text["Label"], Text["k-point"], Text["Eigenvalues"]}}, tab], 
        	  	   Frame      -> True, 
                   Background -> {{GTBackGroundColor1}, {GTBackGroundColor1}}, 
                   Dividers   -> {{2 -> GTDividerColor1, 3 -> GTDividerColor1}, {All, 2 -> GTDividerColor1}}
              ]
        ],
        npt = Length[kpts];
        tab = Table[Table[0, {nbds + 1}], {npt}];
        Do[
           tab[[i, 1]] = bands[[i, 1]];
           Do[
           	  tmp= bands[[i, 2, k]];
             (* tab[[i, k + 1]] = N[bands[[i, 2, k]],4]*)
             tab[[i, k + 1]] = ToString[PaddedForm[N[tmp], {nround, nround}]]
           , {k, 1, nbds}]
        , {i, 1, npt}];
        Print[
              Grid[Join[{{Text["k-point"], Text["Eigenvalues"]}}, tab], 
                   Frame      -> True, 
                   Background -> {{GTBackGroundColor1}, {GTBackGroundColor1}}, 
                   Dividers   -> {{2 -> GTDividerColor1}, {All, 2 -> GTDividerColor1} }
              ]
        ],
        None
     ]
  ];
  (*--- store to files ---*)
  If[Head[ndt] == List, 
  	If[nplt,
  	   GTWriteToFile[{bands,lx}, ndt[[1]]]; 
       Print["Write bands along symmetry lines to file: ", ndt[[1]]];
       If[nwv, 
     	  GTWriteToFile[wav, ndt[[2]]]; 
          Print["Write eigenvectors along symmetry lines to file: ", ndt[[2]]], 
          None
       ],
       GTWriteToFile[bands, ndt[[1]]]; 
       Print["Write bands for DOS to file: ", ndt[[1]]];
       If[nwv, 
     	  GTWriteToFile[wav, ndt[[2]]]; 
          Print["Write eigenvectors for DOS to file: ", ndt[[2]]], 
          None
       ]
  	]   
  ];
  If[Head[ndt] == String, 
  	 GTWriteToFile[{bands,lx}, ndt]; 
     Print["Write bands to file: ", ndt], 
     None
  ];
  (*---Return data---*)
  If[nplt, 
  	 If[nwv, 
  	 	Return[{{bands, lx}, wav}], 
  	 	Return[{bands, lx}]
  	 ], 
     If[nwv, 
     	Return[{bands, wav}], 
     	Return[bands]
     ]
  ]
]

(*
***) 


(****f* /GTBandStructure
! NAME
!  GTBandStructure
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 15.08.2013 : first version
!  * 01.11.2013 : first revision
!  * 20.2.2014  : Correction Non-Orthogonal calculation, tested for graphene
!  * 08.12.2014 : PlotStyle, LabelStyle added to Options
!  * 15.03.2015 : possibility for adding names of irreducible representations and option GOIrepTextStyle
!  * 19.05.2016 : PlotStyle was set as option for the vertical lines at X,M etc, but not for the bandstructure itself, corrected.
!  * 27.12.2016 : lines at symmetry points long enough, GOFermiEnergy introduced
!  * 28.12.2017 : check header and documentation
!  * 23.06.2018 : check header and documentation
!  * 14.06.2020 : In documentation it was said, that one can provide the structure instad a list of k-points. This didn' work.
!                 It is corrected. Output not changed.
!  
! USAGE
!  GTBandStructure[Hamiltonian,kpath,npoints,nband] calculates and plots the band structure for a given Hamiltonian.
!  
! INPUT
!  * Hamiltonian    -- Hamiltonian or list of Hamiltoninan and overlap matrix if nonorthogonal calculation should be performed
!  * kpath          -- list of symmetry points, defining the path. The path is given in a form as it comes from GTBZPath         
!  * npoints        -- number of k-points along a line, the same number for all lines
!  * nband          -- number of bands to plot
! OUTPUT
!  The output is the plot of the corresponding band structure.
!
!  The x-coordinate of the final point in the bandstructure plot is printed. This allows to optimize the 
!  plot in a second step. The band structure plot is given in the usual form.
! 
! GTPack OPTIONS
!  * GOStore         - band structure will be stored if a filename is provided with the option (standard  - no output)
!  * GOTbOrthogonal  - orthogonal/nonorthogonal calculation (standard - orthogonal)
!  * GOShift         - shift of all energies
!  * GOIrepTextStyle - defines style of labels
!  * GOFermiEnergy   - Fermi energy
!  * GOPhotonic
! 
!	   - True   :   photonic bands or phonons (square root of eigenvalues)
!	   - False  :   standard (eigenvalues are energies) 
! 
!  * GOVerbose
!
!      - True   : additional information 
!      - False  : no additional information
!
! STANDARD OPTIONS
!  * Joined       - option of ListPlot
!  * PlotRange    - option of ListPlot
!  * FrameLabel   - option of ListPlot
!  * PlotLabel    - option of ListPlot
!  * LabelStyle   - option of ListPlot
!  * PlotStyle    - option of ListPlot
! GTPack MODULES  
!  GTWriteToFile
! GTPack NOTEBOOKS
!  - 
! DESCRIPTION
!  Calculation and plot of the bandstructure along given lines in the BZ is performed.
!  The Hamiltonian is has to be constructed before. The Hamiltoninan is given as an analytical 
!  form in dependence on the k-vector. The parameters have to substitued by the 
!  actual parameter values befor the calculation can start.
!
!  The Fermi energy has to be provided and will be indicated in the plot. THe calculation of the 
!  Fermi energy has to be done via a DOS calculation.
! LITERATURE
!   -
! TODO
!   -
! RELEASE
!  1.0.0
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)
GTBandStructure::badarg = "If GOTbOrthogonal->`1`  {Hamiltonian, Overlap} necessary";
GTBandStructure::bands   = "Warning : too much bands asked for plot, reset to `1` bands.";

GTBandStructure[hopin_,struc0_,points_,nbdplt_,OptionsPattern[]] :=
                  Module[{hop1,bands,bs,axis,bar,kp,lx,lab,npath,step,del,kpa,x,k,nu,npmax,mi,ma,i,j,band,nbnd,gb,gs,plr,kpath,nop,bd,
                  	      plotl,hop,ov,ov1,join,phot,pev,ndt,nplot,frame,shft,labelstyle,plotstyle,struc,nrnames,kpaxesval,irindex,
                  	      plotnames,Ireptextstyle,ef,efl,xmax,tab}, 
                  Ireptextstyle = OptionValue[GOIrepTextStyle];
                  join=OptionValue[Joined];
                  phot=OptionValue[GOPhotonic];
                  nop =OptionValue[GOTbOrthogonal];
                  plr =OptionValue[PlotRange];
                  ndt =OptionValue[GOStore];
                  frame=OptionValue[FrameLabel];
                  shft      =OptionValue[GOShift];
                  plotl     =OptionValue[PlotLabel];
                  labelstyle=OptionValue[LabelStyle];
                  plotstyle=OptionValue[PlotStyle];
                  ef        =OptionValue[GOFermiEnergy];
                  (*--- Check if names for the irreducible representations are present, prepare k-path ---*)
                  If[StringQ[struc0],
                    kpath = GTBZPath[struc0];
                    tab   = Table[Table["", {j, 1, nbdplt}], {i, 1, Length[kpath[[1]]]}];
                    kpath = Append[kpath, tab];
                    struc = kpath,
                    struc = {None,None,None};
                    struc[[1]] = struc0[[1]];
                    struc[[2]] = struc0[[2]];
                    struc[[3]] = Table[Table["",{j,1,nbdplt}],{i,1,Length[struc[[1]]]}];
                    If[Length[struc0]==3,
                  	   nrnames = Min[Length[struc0[[3,1]]],nbdplt];
                  	   Do[
                  	   	  struc[[3,i,1;;nrnames]] = struc0[[3,i,1;;nrnames]]
                  	   ,{i,1,Length[struc[[1]]]}];
                    ];
                    kpath=struc;
                  ];  
                  kpath=struc;         
                  (*--- Check  if nono-othogonal calculation should be done ---*)                  
                  If[nop,
                  	 hop=hopin,
                  	  If[Length[hopin] != 2,
                         Message[GTBandStructure::badarg, nop]; Return[],
                         hop = hopin[[1]]; ov = hopin[[2]]
                      ]  
                  ];
                  (*--- prepare calculation and plot ---*)	 
	              nbnd=Length[hop];bands={};bs={};axis={};bar={};
	              If[nbdplt>nbnd,
	              	 Message[GTBandStructure::bands, nbnd] ;nplot=nbnd,
	              	 nplot=nbdplt
	              ];            
	              If[StringQ[struc],
	                 kpath=GTBZPath[struc],
	                 kpath=struc
	              ];   
	              kp=kpath[[1]];lx=kpath[[2]];
	              npath=Length[kp]-1;If[npath==0,npath=1,None];
	              kpaxesval = Table[0, {i, 1, npath + 1}];
                  Do[
                  	 del=(kp[[i+1]]-kp[[i]])/(points-1);
                  	 step=Sqrt[del.del];
                  	 If[i==1,
                  	 	x=-step;nu=1,
                  	 	nu=2
                  	  ];
                      Do[
                      	 kpa=kp[[i]]+(k-1)*del;
                      	 x=x+step;axis=Append[axis,x];
                      	 (*--- orthogonal or non-orthogonal calculation ---*)
                      	 hop1=hop /. {\[Xi]->kpa[[1]],\[Eta]->kpa[[2]],\[Zeta]->kpa[[3]]}//Chop;
                      	 If[nop,  
                      	 	If[phot,
                      	 	   pev=Eigenvalues[hop1]//Sort//Chop;bands=Append[bands,Take[Sqrt[pev],{1,nplot}]]
                              ,bd=Eigenvalues[hop1]//Sort//Chop;bd=Take[bd,{1,nplot}];bd=bd-shft;bands=Append[bands,bd]
                            ],
                            ov1=ov /. {\[Xi]->kpa[[1]],\[Eta]->kpa[[2]],\[Zeta]->kpa[[3]]}//Chop; 
                            bd=Eigenvalues[{hop1,ov1}]//Sort//Chop; bd=Take[bd,{1,nplot}]; bd=bd-shft;    	  
                            bands=Append[bands,bd]
                      	 ]               	   
                      ,{k,nu,points}];
                      kpaxesval[[i + 1]] = x;
                   ,{i,1,npath}];
                   (*--- prepare the plot ---*)
                   npmax=Length[axis];mi=Min[bands];ma=Max[bands];
                   bar=Table[{{axis[[i*points-(i-1)]],mi-10},{axis[[i*points-(i-1)]],ma+10}},{i,1,npath}];
                   lab=Table[{axis[[(i-1)*points-(i-1)+1]],lx[[i]]},{i,1,npath+1}];
				   plotnames = {};
                   Do[
                   	  irindex = 1;
                   	  band={};
                      Do[
                      	If[axis[[i]] == kpaxesval[[irindex]], 
                      		plotnames = Append[plotnames, {Text[Style[struc[[3,irindex, k]], {Ireptextstyle}],{kpaxesval[[irindex]], bands[[i, k]]}]}];
  							irindex = irindex + 1];
                      	 band=Append[band,{axis[[i]],bands[[i,k]]}]
                      ,{i,1,npmax}];
                      bs=Append[bs,band]
                   ,{k,1,nplot}];
                   xmax=Max[axis]//N;
                   Print["Maximum Abscissa = ",xmax];	                       
                   gb=ListPlot[bs,Frame->True,PlotStyle->plotstyle,FrameTicks->{lab,Automatic},PlotRange->plr,Joined->join,FrameLabel->frame,PlotLabel->plotl,LabelStyle->labelstyle,Epilog -> plotnames];
                   gs=ListLinePlot[bar,PlotStyle->Black,Ticks->{None,Automatic},PlotRange->plr];  
                   (*--- store to file, output ---*)                                                       
                   If[Head[ndt]===String,
                      GTWriteToFile[bands,ndt];Print["Write to file: ",ndt],
                      None
                   ];
                   If[Abs[ef]>0,
                   	  efl=Graphics[{Red,Line[{{0,ef},{xmax,ef}}]}];
                   	  Show[gb,gs,efl],
                      Show[gb,gs]   
                   ]
                                                             
]
(*
***) 



(****f* /GTBandsPlot
! NAME
!  GTBandsPlot
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!   * 01.11.2013 : 1st version 
!   * 17.02.2016 : PlotStyle implemented
!   * 10.06.2016 : FrameTicks corrected, works now also for one line correctly
!   * 27.12.2016 : lines at symmetry points long enough, GOFermiEnergy introduced
!   * 23.06.2018 : check header and documentation
!   * 14.06.2020 : Message system changed
! USAGE
!  GTBandsPlot[file,nband] plots nband bands from the data in file.
! INPUT
!  * file   - contains the band structure information (k-points , energies/frequencies) and
!             all the information about the labels along the lines. Tis is the output of GTBands.
!  * nband  - number of bands to plot
!
! OUTPUT
!  Plot of the energy bands which are calculated before by means of GTBands.
! 
! GTPack OPTIONS 
!  * GOShift         - shift of the whole band structure
!  * GOFermiEnergy   - Fermi energy is plotted
! STANDARD OPTIONS
!  * PlotRange   - option of ListPlot
!  * Joined      - option of ListPlot 
!  * FrameLabel  - option of ListPlot
!  * PlotLabel   - option of ListPlot
!  * PlotStyle   - option of ListPlot
!
! GTPack MODULES
!  -
! GTPack NOTEBOOKS
!  - 
! DESCRIPTION
!  GTBandsPlot is useful if band structures are calculated by other programs or GTBands and have to be plotted
!  in the GTPack format.
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

GTBandsPlot::bands   = "Warning : not enough bands `1` bands will be used.";

GTBandsPlot[bands_,nbdplt_,OptionsPattern[]]:=Module[
	                   {bds,lx,npt,nl,join,plr,axis,i,evs,nbnd,mi,ma,bs,frame,pltl,bar,lab,band,k,gb,gs,shft,plts,
	                   ef,efl,xmax},
	                   bds=bands[[1]]//N;lx=bands[[2]];npt=Length[bds];nl=Length[lx];
	                  (*--- options ---*) 
                       join  = OptionValue[Joined];
                       plr   = OptionValue[PlotRange];
                       frame = OptionValue[FrameLabel];
                       pltl  = OptionValue[PlotLabel];
                       plts  = OptionValue[PlotStyle];
                       shft  = OptionValue[GOShift];
                       ef    = OptionValue[GOFermiEnergy];
                       axis=Table[0,{i,1,npt}];evs=Table[0,{i,1,npt}];
                       nbnd=Length[bds[[1,4]]];
                       If[nbnd>=nbdplt,
                       	  nbnd=nbdplt,
                       	  Message[GTBandsPlot::bands, nbnd]                    
                       ];
                       Do[
                         axis[[i]]=bds[[i,2]];
                         evs[[i]]=Sort[bds[[i,4]],Less]-shft
                        ,{i,1,npt}];
                        mi=Min[evs];ma=Max[evs];bs={};
                        xmax=Max[axis]//N;
                       Print["Maximum Abscissa = ",xmax];
                       bar=Table[{{axis[[lx[[i,1]]]],mi-10},{axis[[lx[[i,1]]]],ma+10}},{i,1,nl}];
                       lab=Table[{axis[[lx[[i,1]]]],lx[[i,2]]},{i,1,nl}];
                       Do[
                   	      band={};
                          Do[
                      	     band=Append[band,{axis[[i]],evs[[i,k]]}]
                          ,{i,1,npt}];
                          bs=Append[bs,band]
                       ,{k,1,nbnd}];                   
                      gb=ListPlot[bs,Frame->True,PlotStyle->plts,
                      	          FrameTicks->{{Automatic,None},{lab,lab}},PlotRange->plr,Joined->join,FrameLabel->frame,PlotLabel->pltl];
                      gs=ListLinePlot[bar,PlotStyle->Black,Ticks->{None,Automatic},PlotRange->plr]; 
                      If[Abs[ef]>0,
                      	 efl=Graphics[{Red,Line[{{0,ef},{xmax,ef}}]}];
                      	 Show[gb,gs,efl],                                                
                         Show[gb,gs]  
                      ]                                               
]

(*
***) 

(****f* /BZGauss
! NAME
!  BZGauss
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 01.09.2013 : 1st version 
!  * 28.12.2017 : check header and documentation
!  * 23.06.2018 : check header and documentation
! USAGE
!  Calculation of DOS
! INPUT
!   * e     - energy value on the energy axis
!   * sigma - width at half maximum of the Gaussian 
!   * e0    - energetic position of the correponding band   
! OUTPUT
!   contribution of one point of the band structure to the DOS
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
! -
! GTPack NOTEBOOKS
!  - 
! DESCRIPTION
!  The DOS is calculated via a sum over Gaussian peaks located at the positions of the bands of the correponding k-value
! LITERATURE
! 
! TODO
!
! RELEASE
!  1.0.0
!  
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

BZGauss[e_,\[Sigma]_,e0_]:=Module[{},If[e<e0-6*\[Sigma]||e>e0+6*\[Sigma],0,Exp[-(e-e0)^2/2./\[Sigma]^2]/Sqrt[2\[Pi] \[Sigma]^2]]]

(*
***)

(****f* /GTCompatibility
! NAME
!  GTCompatibility
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!   * 1st version October 2013
!   * 06.08.2016 : groups installed inside the command. Now also input as lists is possible.
!   * 07.10.2016 : perhaps problem mit Listable of GTGetSymbol in calculation of class1, class2, changed to Map
!   * 28.12.2017 : check header and documentation
!   * 23.06.2018 : check header and documentation
! USAGE
!  GTCompatibility[group1,group2] calculates the compatibility relations of two point groups.
! INPUT
!  * group1 - group of K at a certain k-point
!  * group2 - a subgroup of group1, otherwise groups are reordered.
! OUTPUT
!   table of compatibility relations
! GTPack OPTIONS
!   * GOIrepNotation - notation of the irreducible representations (standard is "Bethe")
!   * GOFast         - is set to $GOFastValue, can be used to avoid checking and to fasten calculations
!   * GOVerbose:
!
!	   - False :   no information (standard)
!	   - True  :   additional information about the installation of the groups
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  GTGroupQ, GTGetSymbol, GTInstallGroup, GTSubGroupQ, GTCharacterTable
! GTPack NOTEBOOKS
!  - 
! DESCRIPTION
!  At first it is checked which group is the group with lower order. This group 
!  becomes group group2. It is checked, if the group group2 is a subgroupp of group1. If
!  this is not the case the procedure is stopped.
!
!  To get the compatibility relations it is checked how often the Ireps of the group 
!  with lower order, a subgroup of the other group, appear.
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

GTCompatibility[gp1_, gp2_, OptionsPattern[]] := 
  Module[{go, mgp1, mgp2, nt1, nt2, mt, ord1, ord2, ct1, ct2, class1, class2, chars1, chars2, names1, names2, mlt1, mlt2, pos, rep1, 
          rep2, rel, i, k, comp, np, not, notation,fast,grd},
  (*--- options ---*)
  go    = OptionValue[GOVerbose];
  not  = OptionValue[GOIrepNotation];
  fast = OptionValue[GOFast];
  (*--- check Irep notation ---*)
  notation = {"Bethe", "Mulliken", "Bouckaert"};
  If[Intersection[{not}, notation] == {},
     Print["Error: Notation of Ireps not known"]; Abort[],
     None
  ];
  (*--- check and interprete input ---*)
  If[Head[gp1] === List,
     If[fast,
        If[GTGroupQ[gp1],
           None,
           Print["First argument in input is not a group"] ; Abort[]
        ]
     ];
     mgp1 = GTGetSymbol[gp1],
     mgp1 = GTGetSymbol[GTInstallGroup[gp1, GOVerbose -> go]]
  ];
  ord1 = Length[mgp1]; nt1 = gp1;
  If[Head[gp2] === List,
     If[fast,
        If[GTGroupQ[gp2],
           None,
           Print["Second argument in input is not a group"] ; Abort[]
        ]
     ];
     mgp2 = GTGetSymbol[gp2],
     mgp2 = GTGetSymbol[GTInstallGroup[gp2, GOVerbose -> go]]
  ];
  ord2 = Length[mgp2]; nt2 = gp2;
  If[ord1 == ord2,
     Print["Error: the groups are of same order"],
     If[ord1 > ord2,
  (*--- reorder groups ---*)
        None,
        mt = mgp1; mgp1 = mgp2; mgp2 = mt;
        mt = nt1; nt1 = nt2; nt2 = mt;
        ord1 = Length[mgp1]; ord2 = Length[mgp2]
     ];   
     If[GTSubGroupQ[mgp2, mgp1],
        ct1    = GTCharacterTable[mgp1, GOVerbose -> go, GOIrepNotation -> not];
        ct2    = GTCharacterTable[mgp2, GOVerbose -> go, GOIrepNotation -> not];
        class1 = GTGetSymbol[#]& /@ ct1[[1]];
        class2 = GTGetSymbol[#]& /@ ct2[[1]];
        chars1 = ct1[[2]]; chars2 = ct2[[2]]; 
        names1 = ct1[[3]]; names2 = ct2[[3]];
        mlt1   = Length[ct1[[1]]]; mlt2 = Length[ct2[[1]]];
        rel    = {};
    (*--- loop over all reps of the large group---*)
        Do[
           rep1 = {};
           Do[
              pos  = Position[class1, mgp2[[k]]][[1, 1]];
              rep1 = Append[rep1, chars1[[i, pos]]]
           , {k, 1, ord2}];
     (*---loop over all reps of the small group---*)            
           comp = {}; 
           Do[
              rep2 = {};
              Do[
                 pos = Position[class2, mgp2[[k]]][[1, 1]];
                 rep2 = Append[rep2, chars2[[j, pos]]]
              , {k, 1, ord2}];
              np = rep1.rep2/ord2;
              If[np == 0,
                 None,
                 comp = Append[comp, names2[[j]]]
              ]
            , {j, 1, mlt2}];
            rel = Append[rel, {names1[[i]], comp}]
        , {i, 1, mlt1}];
        grd=Grid[rel // Transpose, Frame -> All,Dividers -> {Black, {2 -> GTDividerColor1}},
                                   Background -> {None, {1 -> GTBackGroundColor1}}
                ];Print[grd],
        Print["Error: no subgroup relationship"]
     ]
  ]  
]



(*
***)



(****f* /GTDensityOfStates
! NAME
!  GTDensityOfStates
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!   * 05.09.2013 : 1st version 
!   * 15.11.2013 : revision 
!   * 01.09.2014 : check and revision during PP implementation
!   * 27.09.2014 : scaling factor for DOS and option PlotStyle introduced
!   * 10.02.2015 : GOPlot -> GOPlotDos
!   * 17.06.2016 : option GOFermiEnergy inserted
!   * 30.11.2016 : unification of all DOS calculations, use of TwoAxisListPlot
!   * 28.12.2017 : check header and documentation 
!   * 24.06.2018 : check header and documentation  
!   * 29.05.2020 : new mesage system implemented, better error handling, improvement of documentation
!
! USAGE
!  GTDensityOfStates[Hamiltonian,structure,parameters] calculates and plots the density of states for a given Hamiltonian.
! INPUT
!  * ham    - Hamiltonian or bands 
!  * struc  - structure of the system    
!  * parm   - list of parameters for DOS calculation:
!
!            1 -  n    - index k-mesh
!            2 -  a    - scaling factor k-mesh
!            3 - emin  - minimum energy
!            4 - emax  - maximum enery
!            5 - ne    - number of energy points  
!            6 - sigma - width of the Gaussians (if = o then equal to energy steps)
!            7 - scale - scaling factor for the DOS
!  
! OUTPUT
!  Plot of the DOS and/or the integrated DOS
!  
! GTPack OPTIONS
!   * GOBands:
!
!          "All"  - all bands 
!           n     - the first n bands
!           {..}  - the bands in the list will be used
!   * GOPhotonic:
!   
!           True  - switch on photonic case  
!           False - electronic structure (standard)
!   * GOPlotDos:
!
!           "DOS"   - plot DOS (Standard) 
!           "IDOS"  - plt IDOS
!           "ALL"   - plot DOS and IDOS 
!   * GOFermiEnergy:
!
!           0       - Fermi energy not calculated (standard) 
!          {nel,es} - nel - number of electrons, es - energy to start the search
!   * GOStore      - DOS will be stored, if a filename is given   
! STANDARD OPTIONS
!  * PlotRange     - Mathematica option for Plot commands (standard : All)
!  * FrameLabel    - standard: {"Energy","DOS"}
!  * PlotLabel     - standard: "Density of states"
!  * PlotStyle     - standard: Blue
! 
! GTPack MODULES
!   GTBZPointMesh, BZGauss, GTWriteToFile, TwoAxisListPlot
!  
! GTPack NOTEBOOKS
!  - 
! DESCRIPTION
!  The DOS is calculated via a sum of over Gaussian peaks.
!  TwoAxisListPlot is an external package
! LITERATURE
!  -
! TODO
!  implementation of the true partial DOS
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTDensityOfStates::dosp    =  "Selector for kind of plot is `1`. It has to be DOS, IDOS, PDOS or ALL.";
GTDensityOfStates::plots   =  "Wrong Pattern in Plotstyle";
GTDensityOfStates::idos    =  "Use IDOS or ALL in GoPlotDos to calculate Fermi energy.";
GTDensityOfStates::bnds    =  "Wrong descriptor in GOBands.";

 
(* ::Subsubsection:: *)
(* Region Title *)
GTDensityOfStates[hop_, struc_, parm_, OptionsPattern[]] := Module[
 	{data,bds,phot,plr,plt,frame,pltl,style,fen,test,cdos,cidos,n,a,qp,nkp,nbands,bind,i,ii,
 	 ew,ham,ev0,ev,bands,emin,emax,ne,de,scale,axis,dos,k,bs,func,idos,id,xs,nel,flog,rule,
 	 xf,gf,plot,x},
  (*--- options- --*)
    data  = OptionValue[GOStore];
    bds   = OptionValue[GOBands];
    If[bds=="ALL"||Head[bds]===List||Head[bds]===Integer,
       None,
       Message[GTDensityOfStates::bnds]; Return[]
    ];   	
    phot  = OptionValue[GOPhotonic];
    plr   = OptionValue[PlotRange];
    plt   = OptionValue[GOPlotDos];
    frame = OptionValue[FrameLabel];
    pltl  = OptionValue[PlotLabel];
    style = OptionValue[PlotStyle];
    fen   = OptionValue[GOFermiEnergy];
    If[(Head[fen] === List||fen>0),
    	If[plt=="IDOS"||plt=="ALL",
           None,
           Message[GTDensityOfStates::idos]; Return[]
    	],
    	None   
    ];
    test  = Intersection[{plt}, {"DOS", "IDOS", "ALL"}];
    If[test == {},
        Message[GTDensityOfStates::dosp,plt]; Return[],
       None
    ];
    If[Head[style] === RGBColor,
       cdos = cidos = style,
       If[Length[style] == 2,
          cdos = style[[1]]; cidos = style[[2]],
          Message[GTDensityOfStates::plots]; Return[]
       ]
    ];
  (*--- k-mesh ---*)
    n   = parm[[1]];
    a   = parm[[2]];
    qp  = GTBZPointMesh[n, a, struc];
    nkp = Length[qp];
  (*--- bands to use in calculation ---*)
    nbands = Length[hop];
    If[bds === "ALL", 
       nbands = Length[hop]; 
       bind = Table[i, {i, 1, nbands}], 
       If[Head[bds] === List, 
       	  bind   = bds; 
       	  nbands = Length[bind], 
          nbands = bds; 
          bind = Table[i, {i, 1, nbands}]
       ]
    ];
  (*---Calculaton of the bands---*)
    ew = {};
    Do[
       ham = hop /. {\[Xi] -> qp[[i, 1]], \[Eta] -> qp[[i, 2]], \[Zeta] -> qp[[i, 3]]};
       ev0 = Eigenvalues[ham] // Sort;
       ev  = Table[0, {nbands}];
       Do[
       	  ev[[ii]] = ev0[[bind[[ii]]]]
       , {ii, 1, nbands}];
       If[phot,
       	  ew = Append[ew, Sqrt[ev] // Chop], 
          ew = Append[ew, ev] // Chop
       ]
    , {i, 1, nkp}];
    If[nbands <= 1, 
       bands = {Flatten[ew]}; nbands = 1, 
       bands = Transpose[ew]
    ];
  (*---Prepataton of calculation of the DOS---*)
    emin  = parm[[3]]; 
    emax  = parm[[4]]; 
    ne    = parm[[5]];
    scale = parm[[7]];
    de    = (emax - emin)/(ne - 1.);
    If[parm[[6]] == 0, 
       \[Sigma] = de, 
       \[Sigma] = parm[[6]]
    ];
    axis = Table[(i - 1)*de + emin, {i, 1, ne}]; 
    dos  = Table[0, {ne}];
  (*---Calculation of the DOS---*)
    Do[
       Do[
       	  dos = dos + scale*
          Map[BZGauss[#, \[Sigma], bands[[k, i]]] &, axis]
       , {i, 1, nkp}]
   , {k, 1, nbands}];
   bs = Transpose[{axis, dos/nkp}];
  (*---Calculation of integrated DOS---*)  
   If[plt == "IDOS" || plt == "ALL",
      func = Interpolation[bs]; idos = Table[0, {ne}];
      Do[
         idos[[i]] = Integrate[func[x], {x, emin, axis[[i]]}]
     , {i, 1, ne}];
     id = Transpose[{axis, idos}];
     Print["integrated DOS = ", idos[[ne]]];
  (*---calculation of Fermi energy---*)
     If[Head[fen] === List,
        flog = True;  
        xs   = fen[[2]];  
        nel  = fen[[1]],
        If[fen > 0,
           flog = True; xs = (axis[[ne]] + axis[[1]])/2.; nel = fen,
           flog = False
        ]
     ];
     If[flog,
        func = Interpolation[id];          
        rule = FindRoot[func[x] == nel, {x, xs}];
        xf = x /. rule;
       Print["Fermi Energy \!\(\*SubscriptBox[\(E\), \(F\)]\)= ", xf, 
             " (units see band structure)"
            ];
       gf = Graphics[{Red,Thick, Line[{{xf, -1}, {xf, Max[dos/nkp] + 10}}]}],
       None
     ], 
     None 
    ];
  (*---Output to File---*)
    If[Head[data] === String,
       If[plt == "DOS",
          GTWriteToFile[bs, data],
          GTWriteToFile[{bs, id}, data]
       ];
       Print["Write to file: ", data],
             None
       ];
  (*--- Plots ---*)   
       If[plt == "DOS",
          plot = ListPlot[bs, Frame -> True, PlotStyle -> cdos, FrameLabel -> frame, Joined -> True, 
          	              PlotRange -> plr, PlotLabel -> pltl
          	             ],
          None
       ];
       If[plt == "IDOS",        
          plot = ListPlot[id, Frame -> True, PlotStyle -> cidos, FrameLabel -> frame, 
          	              PlotLabel -> pltl, PlotRange -> plr, Joined -> True
          	             ],     
                                
         None
       ];
       If[plt == "ALL",
          plot = TwoAxisListPlot[{bs, id}, Frame -> True, PlotStyle -> {cdos, cidos}, 
          	     FrameLabel -> frame,PlotLabel -> pltl,Joined -> True],
          None
       ];
       If[flog, 
       	  plot = {plot, gf}, 
       	  None
      ];
      Show[plot]
  ]
  

(*
***) 

(****f* /GTDensityOfStatesPlot
! NAME
!  GTDensityOfStatesPlot
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 1st versionNovember 2013
!  * 29.09.2014 : instead of single parameters plist of parameters introduced
!  * 10.02.2015 : correction GOBands
!  * 26.02.2016 : GTDosPlot -> GTDensityOfStatesPlot
!  * 17.06.2016 : new option GOFermiEnergy, better output of data
!  * 26.08.2016 : new option GOPlotDos to select what is plotted, The interpretation of the values of Plotstyle is changed.
!                 if it is a list , the colors are used for DOS and IDOS , otherwise both curves have the same color.
!  * 30.11.2016 : TwoAxisListPlot introduced for ALL
!  * 28.12.2017 : check header and documentation  
!  * 24.06.2017 : check header and documentation   
!  * 29.05.2020 : new mesage system implemented, better error handling, improvement of documentation
!       
! USAGE
!  GTDensityOfStatesPlot[band structure,parameters] calculates the density of states from a precalculated band structure.
! INPUT
!
!  * bnds   - band structure as calculated with GTBands
!  * parm={emin,emax,ne,s0,scale}:
!
!     - emin   - minimum energy
!     - emax   - maximum enery
!     - ne     - number of energy points  
!     - s0     - width of the Gaussians (if = o then equal to energy steps)
!     - scale  - scaling factor for the DOS
! 
! OUTPUT
!  Plot of the DOS and/or the integrated DOS. DOS can be stored in file.
!
! GTPack OPTIONS
!   * GOBands:
!
!          "All"  - all bands 
!           n     - the first n bands
!           {..}  - the bands in the list will be used
!   * GOPlotDos:
!
!           "DOS"   - plot DOS (Standard) 
!           "IDOS"  - plt IDOS
!           "ALL"   - plot DOS and IDOS 
!   * GOFermiEnergy:
!
!           0       - Fermi energy not calculated (standard) 
!          {nel,es} - nel - number of electrons, es - energy to start the search
!   * GOStore      - DOS will be stored, if a filename is given   
!
! STANDARD OPTIONS
!  PlotRange     - option of ListPlot
!  FrameLabel    - option of ListPlot
!  PlotLabel     - option of ListPlot
!  PlotStyle     - option of ListPlot
! GTPack MODULES
!   BZGauss, GTWriteToFile, TwoAxisListPlot 
! 
! GTPack NOTEBOOKS
!  -  
! DESCRIPTION
!  The DOS is calculated via a sum of over Gaussian peaks.
!  TwoAxisListPlot is an external package
! LITERATURE
!  -
! TODO
!  implementation of the true partial DOS
! RELEASE
!  1.0.0
! PROBLEMS
!  in priciple one could arrange only one ListPlot instead of two. Then the interpretation of PlotStyle woud be directly
!  that of ListPlot. Seems to be a general problem: If Mathematica commands are used, then options should be kept alive in 
!  their original meaning. (something for GTPack2.0) -> What does this mean?
! SOURCE
!
!--------------------------------------------------------------------------------
!
*)

GTDensityOfStatesPlot::dosp    =  "Selector for kind of plot is `1`. It has to be DOS, IDOS, PDOS or ALL.";
GTDensityOfStatesPlot::plots   =  "Wrong Pattern in Plotstyle";
GTDensityOfStatesPlot::idos    =  "Use IDOS or ALL in GoPlotDos to calculate Fermi energy.";
GTDensityOfStatesPlot::bnds    =  "Wrong descriptor in GOBands.";

GTDensityOfStatesPlot[bnds_, parm_, OptionsPattern[]] := 
 Module[{data, plr, sbnds, nbands, bind, i, de, axis, dos, bands, nkp,
         k, bs, x, id, frame, pltl, emin, emax, ne, s0, scale, 
         fen, idos, func,flog,nel,xs,xf,rule,gf,pstyle,pldos,test,plot,cdos,cidos},
  (*--- parameters ---*)
  emin  = parm[[1]]; 
  emax  = parm[[2]]; 
  ne    = parm[[3]]; 
  s0    = parm[[4]]; 
  scale = parm[[5]];
  (*--- check options ---*)
  data  = OptionValue[GOStore];
  sbnds = OptionValue[GOBands];
  If[sbds=="ALL"||Head[sbds]===List||Head[sbds]===Integer,
     None,
     Message[GTDensityOfStatesPlot::bnds]; Return[]
  ];  
  fen   = OptionValue[GOFermiEnergy];
  plr   = OptionValue[PlotRange];
  frame = OptionValue[FrameLabel];
  pltl  = OptionValue[PlotLabel];
  pstyle= OptionValue[PlotStyle];
  pldos = OptionValue[GOPlotDos];
  test  = Intersection[{pldos},{"DOS","IDOS","ALL"}];
  If[test=={},
  	 Message[GTDensityOfStatesPlot::dosp,pldos]; Return[],
  	 None
  ];
  If[(Head[fen] === List||fen>0),
      If[pldos=="IDOS"||pldos=="ALL",
         None,
         Message[GTDensityOfStatesPlot::idos]; Return[]
      ],
      None   
  ];
  If[Head[pstyle] === RGBColor,
     cdos = cidos = pstyle,
     If[Length[pstyle] == 2 ,
        cdos = pstyle[[1]]; cidos = pstyle[[2]],
         Message[GTDensityOfStatesPlot::plots]; Return[]
   ]  
  ];
  If[Head[sbnds] === String,
     nbands = Length[bnds[[1, 2]]]; 
     bind = Table[i, {i, 1, nbands}],
     If[Head[sbnds] === List,
        bind = sbnds;  nbands = Length[bind],
        nbands = sbnds; bind = Table[i, {i, 1, nbands}]
     ]
  ];
  (*--- prepare Data for plot ---*)
  de = (emax - emin)/(ne - 1.);
  axis = Table[(i - 1)*de + emin, {i, 1, ne}]; dos = Table[0, {ne}];
  If[s0 == 0,
     \[Sigma] = de,
     \[Sigma] = s0
  ];
  (*--- DOS calculation ---*)
  bands = Transpose[Transpose[bnds][[2]]];
  nkp = Length[bands[[1]]];
  Do[
     Do[
        dos = dos + 
        scale*Map[BZGauss[#, \[Sigma], bands[[bind[[k]], i]]] &, axis] 
     , {i, 1, nkp}]
  , {k, 1, nbands}];
  bs = Transpose[{axis, dos/nkp}];
  (*--- integrated DOS ---*)
  If[Head[fen] === List||fen>0,
  	 flog=True,
     flog=False
  ];
  If[pldos=="IDOS"||pldos=="ALL"||flog,
     idos = Table[0, {ne}];
     func = Interpolation[bs]; idos = Table[0, {ne}];
     Do[
        idos[[i]] = Integrate[func[x], {x, emin, axis[[i]]}]
     , {i, 1, ne}];
     id = Transpose[{axis, idos}];
     Print["integrated DOS = ", idos[[ne]]];
  (*--- calculation of Fermi energy ---*)
    If[Head[fen] === List,
       flog = True; xs = fen[[2]]; nel = fen[[1]],
       If[fen > 0,
          flog = True; xs = (axis[[ne]] + axis[[1]])/2.; nel = fen, 
          flog = False
       ]
    ];
    If[flog,
       func = Interpolation[id];
       rule = FindRoot[func[x] == nel, {x, xs}];
       xf = x /. rule;
       Print["Fermi Energy \!\(\*SubscriptBox[\(E\), \(F\)]\)= ", xf, " (units see band structure)"];
       gf = Graphics[{Red,Thick, Line[{{xf, -1}, {xf, Max[dos/nkp] + 10}}]}],
       None
    ],
    None
  ];
  (*--- store data---*)  
  If[Head[data] === String,
  	If[pldos=="DOS", 
  	   GTWriteToFile[bs, data], 
  	   GTWriteToFile[{bs, id}, data]
  	]; 
    Print["Write to file: ", data], 
    None
  ];
   (*--- Plots ---*)   
       If[pldos == "DOS",
          plot = ListPlot[bs, Frame -> True, PlotStyle -> cdos, FrameLabel -> frame, Joined -> True, 
          	              PlotRange -> plr, PlotLabel -> pltl
          	             ],
          None
       ];
       If[pldos == "IDOS",        
          plot = ListPlot[id, Frame -> True, PlotStyle -> cidos, FrameLabel -> frame, 
          	              PlotLabel -> pltl, PlotRange -> plr, Joined -> True
          	             ],     
                                
         None
       ];
       If[pldos == "ALL",
          plot = TwoAxisListPlot[{bs, id}, Frame -> True, PlotStyle -> {cdos, cidos}, 
          	     FrameLabel -> frame,PlotLabel -> pltl,Joined -> True],
          None
       ];
       If[flog, 
       	  plot = {plot, gf}, 
       	  None
      ];
      Show[plot]
  ]
 

(*
***) 

(****f* /GTDensityOfStatesRS
! NAME
!  GTDensityOfStatesRS
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!   * 1st version September 2014
!   * 28.09.2014 : scaling factor introduced
!   * 29.09.2014 : GTDensityOfStatesRealSpace -> GTDensityOfStatesRS
!   * 10.02.2015 : GOPlot ->GOPlotDos
!   * 25.10.2016 : calculation of Fermi energy implemented, The style of plotting and the choice of colors is 
!                  arranged like in GTDensityOfStatesPlot 
!   * 30,11,2016 : TwoAxisListPlot introduced
!   * 28.12.2017 : check header and documentation   
!   * 06.05.2018 : possibility added to read externally calculated eigenvalues
!   * 24.06.2018 : check header and documentation 
!   * 29.05.2020 : new mesage system implemented, better error handling, improvement of documentation  
!
! USAGE
!  GTDensityOfStatesRS[Hamiltionian,parameters] is used to calculate the density of states of a real space Hamiltionian.
! INPUT
!  * ham    - real space Hamiltonian or list of eigenvalues 
!  * dosp   - list of parameters for DOS calculation:
!
!            1 - nb    - number of atoms in the basis
!            2 - natom - number of atoms in the cluster
!            3 - emin  - minimum energy
!            4 - emax  - maximum enery
!            5 - ne    - number of energy points  
!            6 - sigma - width of the Gaussians (if = o then equal to energy steps)
!            7 - scale - scaling factor for the DOS
! OUTPUT
!  Plot of the DOS and/or integrated DOS

! GTPack OPTIONS
!   * GOPlotDos:
!
!           "DOS"   - plot DOS (Standard) 
!           "IDOS"  - plt IDOS
!           "ALL"   - plot DOS and IDOS 
!   * GOFermiEnergy:
!
!           0       - Fermi energy not calculated (standard) 
!          {nel,es} - nel - number of electrons, es - energy to start the search
!   * GOStore      - DOS will be stored, if a filename is given   
!
! STANDARD OPTIONS
!  * PlotRange  - Mathematica option for Plot commands (standard : All)
!  * FrameLabel - standard: {"Energy","DOS"}
!  * PlotLabel  - standard: "Density of states"
!  * PlotStyle  - standard: "Red"
!  
! GTPack MODULES
!   BZGauss, GTWriteToFile, TwoAxisListPlot
! GTPack NOTEBOOKS 
!  -
! GTPack NOTEBOOKS
!  - 
! DESCRIPTION
!  The DOS is calculated via a sum of over Gaussian peaks.
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

GTDensityOfStatesRS::dosp    =  "Selector for kind of plot is `1`. It has to be DOS, IDOS, PDOS or ALL.";
GTDensityOfStatesRS::plots   =  "Wrong Pattern in Plotstyle";
GTDensityOfStatesRS::idos    =  "Use IDOS or ALL in GoPlotDos to calculate Fermi energy.";

GTDensityOfStatesRS[hop_, dosp_, OptionsPattern[]] := Module[
 {nbas,natom,emin,emax,ne,de,axis,energ,data,plr,plt,pls,frame,pltl,ev,dos,i,bs,func,idos,x,scale,test,
  \[Sigma],flog,id,fen,xs,nel,cdos,cidos,rule,xf,gf,plot},
  (*--- interpretation of options ---*)  
     data  = OptionValue[GOStore];
     plr   = OptionValue[PlotRange];
     fen   = OptionValue[GOFermiEnergy];
     plt   = OptionValue[GOPlotDos];
     pls   = OptionValue[PlotStyle];
     frame = OptionValue[FrameLabel];
     pltl  = OptionValue[PlotLabel];
     test  = Intersection[{plt}, {"DOS", "IDOS", "ALL"}];
     If[test == {},
        Message[GTDensityOfStatesRS::dosp,plt]; Return[],
        None
     ];
     If[(Head[fen] === List||fen>0),
      If[plt=="IDOS"||plt=="ALL",
         None,
         Message[GTDensityOfStatesRS::idos]; Return[]
      ],
      None   
  ];
     If[Head[pls] === RGBColor,
        cdos = cidos = pls,
        If[Length[pls] == 2,
           cdos = pls[[1]]; cidos = pls[[2]],
           Message[GTDensityOfStatesRS::plots]; Return[]
        ]
     ];
  (*--- parameters and energy axis ---*)   
     nbas  = dosp[[1]]; natom = dosp[[2]]; 
     emin  = dosp[[3]]; emax  = dosp[[4]]; 
     ne    = dosp[[5]]; scale = dosp[[7]];
     de    = (emax - emin)/(ne - 1);
     axis  = Table[energ, {energ, emin, emax, de}];
  (*--- calculate the DOS from Hamiltonian or energies ---*)
     If[Length[Dimensions[hop]] == 2,
        ev = Eigenvalues[hop],
        ev = hop
     ];
     dos   = Table[0, {Length[axis]}];
     If[dosp[[6]] == 0,
        \[Sigma] = de,
        \[Sigma] = dosp[[6]]
     ];
     Do[
        dos = dos + scale*Map[BZGauss[#, \[Sigma], ev[[i]]] &, axis]
     , {i, 1, Length[ev]}];
     bs = Transpose[{axis, dos*nbas/natom}];
  (*--- integrated DOS ---*)
     flog = False;
     If[plt == "IDOS" || plt == "ALL",
        func = Interpolation[bs]; idos = Table[0, {ne}];
        Do[
           idos[[i]] = Integrate[func[x], {x, emin, axis[[i]]}]
        , {i, 1, ne}];
        id = Transpose[{axis, idos}];
        Print["integrated DOS = ", idos[[ne]]];
  (*--- calculation of Fermi energy ---*)
        If[Head[fen] === List,
           flog = True;  
           xs   = fen[[2]];  
           nel  = fen[[1]],
           If[fen > 0,
              flog = True; xs = (axis[[ne]] + axis[[1]])/2.; nel = fen,
              flog = False
           ]
        ];
        If[flog,
          func = Interpolation[id];          
          rule = FindRoot[func[x] == nel, {x, xs}];
          xf = x /. rule;
          Print["Fermi Energy \!\(\*SubscriptBox[\(E\), \(F\)]\)= ", xf, " (units see band structure)"];
          gf = Graphics[{Green, Line[{{xf, -1}, {xf, Max[dos*nbas/natom] + 10}}]}],
          None
        ], 
        None
     ];
  (*--- Output to File ---*)
     If[Head[data] === String,
        If[plt == "DOS",
           GTWriteToFile[bs, data],
           GTWriteToFile[{bs, id}, data]
        ];
        Print["Write to file: ", data],
        None
     ];
  (*--- Plots ---*)   
     If[plt == "DOS",
        plot = ListPlot[bs, Frame -> True, PlotStyle -> cdos, FrameLabel -> frame, Joined -> True, 
        	            PlotRange -> plr, PlotLabel -> pltl
        	           ],
                       None
       ];
       If[plt == "IDOS",
          plot = ListPlot[id, Frame -> True, PlotStyle -> cidos, FrameLabel -> frame, PlotLabel -> pltl,      
                          PlotRange -> plr, Joined -> True
                         ],
          None
      ];
      If[plt == "ALL",           
         plot = TwoAxisListPlot[{bs, id}, Frame -> True, PlotStyle -> {cdos, cidos}, FrameLabel -> frame, 
                                PlotLabel -> pltl,Joined -> True
                               ],
         None
     ];
     If[flog, 
     	plot = {plot, gf}, 
     	None
     ];
     Show[plot]
  ]
 


(*
***)


(****f* /GTBandsPlotImprove
! NAME
!  GTBandsPlotImprove
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 03/20/2016 : first implementation
!  * 28.12.2017 : check header and documentation
!  * 24.06.2018 : check header and documentation    
!
! USAGE
!  GTBandsPlotImprove[bands,text,mod,minb,maxb] is used to improve the positions of the Irep labels in a band structure plot.
! INPUT
!   * bands         - band structure a a result of the band structure calculation 
!   * text          - Irep labels {{Irep,{x,y}},{Irep,{x,y}},...}
!   * mod           - list of modifications: {{kp,bnd,shft},{kp,bnd,shft},...}
!   * [minb,maxb]   - bands to be considered
!  
! OUTPUT
!  A table with the original and the changed Irep positions, the new plot of the band structure including the modified labels
!
! GTPack OPTIONS
!   * GOLabelShift  - rigid shift of all labels  {0,0} is standard
! 
! STANDARD OPTIONS
!  * Joined        - True
!  * PlotRange     - Mathematica option for Plot commands (standard : All)
!  * FrameLabel    - standard: {"","Energy (eV)"}
!  * PlotLabel     - standard: "BandStrucutre"
!  * PlotStyle     - standard: {{Thin,Black}}
!
! GTPack MODULES
!   GTBandsPlot
! GTPack NOTEBOOKS
!  - 
! DESCRIPTION
!  The tool can be used in an interactive way to optimize the positions of the labels. It is possible to delete 
!  labels. If shift="D" the corresponding label will be deleted.
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  This task can be done also by means of the graphics palette.
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTBandsPlotImprove[bands_, text_, mod_, minb_, maxb_, OptionsPattern[]] := Module[
	{style,range,join,delx,dely,nmod,posx,posy,reps,pos1,pos2,reps1,nl,kp,bnd,shft,table,i,headline,txt,frame,plab,bplot},
  (*--- options ---*)
   style = OptionValue[PlotStyle];
   range = OptionValue[PlotRange];
   frame = OptionValue[FrameLabel];
   plab  = OptionValue[PlotLabel];
   join  = OptionValue[Joined];
   delx  = OptionValue[GOLabelShift][[1]];
   dely  = OptionValue[GOLabelShift][[2]];
  (*--- rigid shift of all labels ---*)
   nmod = Length[mod]; posx = {}; posy = {}; reps = {};
   Do[
      posx = Append[posx, text[[i, 2, 1]] + delx];
      posy = Append[posy, text[[i, 2, 2]] + dely];
      reps = Append[reps, text[[i, 1]]]
   , {i, 1, Length[text]}];
  (*--- change selected labels ---*)
   nl = maxb - minb + 1;
   pos1 = Partition[posy, nl]; pos2 = pos1; reps1 = Partition[reps, nl];
  (*--- modifications ---*)
   Do[
   	  kp   = mod[[i, 1]]; 
   	  bnd  = mod[[i, 2]]-minb+1; 
   	  shft = mod[[i, 3]];
      If[shft == "D",
         pos2[[kp, bnd]] = "D",
         pos2[[kp, bnd]] = pos1[[kp, bnd]] + shft
      ];
   , {i, 1, nmod}];
   table = Table[{reps1[[i]], pos1[[i]], pos2[[i]]}, {i, 1, Length[reps1]}];
   headline = Table[i, {i, minb, maxb}]; 
   table = Flatten[Prepend[table, {headline}], 1];
   Print[Grid[table, Frame -> All, Background -> {None, {GTBackGroundColor1, {GTBackGroundColor2, White, White}}}]];
  (*--- reconstruct text ---*)
  txt  = {}; 
  posy = Flatten[pos2, 1];
  Do[
  	 If[posy[[i]] == "D",
        None,
        txt = Append[txt, Text[reps[[i]], {posx[[i]], posy[[i]]}]]
     ]
  , {i, 1, Length[text]}];
  bplot = GTBandsPlot[bands, maxb, PlotStyle -> style, Joined -> join, PlotRange -> range,FrameLabel->frame,PlotLabel->plab];
  Show[bplot,Graphics[txt]]
]
  

(*
***)



(****f* /GTFermiSurfaceCut
! NAME
!  GTFermiSurfaceCut
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 15.05.2016 : first implementation
!  * 28.12.2017 : check header and documentation   
!  * 24.06.2018 : check header, documentation page created
!
! USAGE
!  GTFermiSurfaceCut[ham,ef,nbands,npt,trans,mima] calculates a cut through the Fermi surface. 
! INPUT
!  * ham    - parametrized Hamiltonian
!  * ef     - Fermi energy
!  * nbands - list of bands tat are expected to cut the Fermi energy, i.e. contribute to the Femi surface
!  * npt    - number of points per spatial direction to define k-mesh
!  * trans  - {nv,shft}:
!
!       - nv   - normal vector of the plane
!       - shft - shift vector with respect to the Gamma point
!  * mima   - {min,max} defines the quadratic region of the cutting plane
!
! OUTPUT
!  Plot of the Fermi surface cut  with ContourPlot
!  
! GTPack OPTIONS
!  * GOVerbose:
!
!          - True  -  additional information
!          - False -  no additional information (standard)
!
!  * GOPlot:
!
!          - True  -   plot the cuts (Standard)  
!          - False -   list of graphics objects 
! STANDARD OPTIONS
!  * ContourStyle - Option of ContourPlot
! GTPack MODULES
!  -
! GTPack NOTEBOOKS
!  - 
! DESCRIPTION
!  The Fermi surface corresponding to a Hamiltonian and Fermi energy  is calculated. A plane is defined that cuts the 
!  surce a square int e plane is defined to get a quadratic part of the cut .
! LITERATURE
! 
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

GTFermiSurfaceCut[ham_, ef_, nbands_, npt_, trans_, mima_,  OptionsPattern[]] := Module[
	{verb,plot,xmin,xmax,del,kpt,kmesh,x,y,i,j,nkp,rv,rotplane,rotmat,kptt,test,sel,cst,
	 tab,nb ,ib,hc,ew,ewl,gr,fl,k,fu,xg,yg},
  (*---options---*)
     verb  = OptionValue[GOVerbose];
     plot  = OptionValue[GOPlot];
     cst   = OptionValue[ContourStyle];
  (*---k-mesh---*)
     xmin  = mima[[1]]; 
     xmax  = mima[[2]];
     del   = (xmax - xmin)/(npt - 1);
     kpt   = {}; 
     kmesh = {};
     Do[
        Do[
           x     = (i - 1)*del + xmin;
           y     = (j - 1)*del + xmin;
           kpt   = Append[kpt, {x, y, 0}];
           kmesh = Append[kmesh, {x, y}]
        , {i, 1, npt}]
     , {j, 1, npt}];
     nkp = Length[kpt];
  (*---transform k-mesh---*)
  (*--- rotation in plane ---*) 
     rv = {trans[[1, 1]], trans[[1, 2]], 0};
     If[Norm[rv] > 0,
        rotplane = RotationMatrix[{{1, 0, 0}, rv}];
        kpt      = rotplane.# & /@ kpt // N,
        None
     ];
  (*--- final rotation ---*)
     rotmat = RotationMatrix[{{0, 0, 1}, trans[[1]]}];
     kptt   = rotmat.# + trans[[2]] & /@ kpt // N;
     If[verb,
  	    Print[nkp, " k-points"];
        Print["rotation of the mesh corners"];
        Print["start with:"];
        test = {{xmin, xmin, 0}, {xmin, xmax, 0}, {xmax, xmax, 0}, {xmax, xmin, 0}};
        Print[test];
        If[Norm[rv] > 0,
           Print["rotation in plane:"];
           test = rotplane.# & /@ test // FullSimplify;
           Print[test],
           None
        ];
        Print["final rotation:"];
        test = rotmat.# & /@ test // N;
        Print[test],
        None
     ];
  (*---select bands---*)
     sel = Table[0, {Length[ham]}]; tab = {};
     nb  = Length[nbands];
     Do[
     	sel[[nbands[[ib]]]] = 1
     , {ib, 1, nb}];
  (*---calculate band structure---*)tab = {};
     Do[
     	Clear[\[Xi], \[Eta], \[Zeta]];
        hc  = ham /. {\[Xi] -> kptt[[i, 1]], \[Eta] -> kptt[[i, 2]], \[Zeta] -> kptt[[i, 3]]};
        ew  = Sort[Eigenvalues[hc]];
        ewl = Pick[ew, sel, 1];
        tab = Append[tab, ewl]
     , {i, 1, Length[kptt]}];
  (*---generate plots--*)
     gr = Table[0, {nb}];
     fl = Table[0, {nb}];
     Do[
     	fu = {};
        Do[
           fu      = Append[fu, {kmesh[[k]], tab[[k, i]]}]
        , {k, 1, nkp}];
        fl[[i]] = Interpolation[fu];
        gr[[i]] = ContourPlot[fl[[i]][xg, yg] == ef, {xg, xmin, xmax}, {yg, xmin, xmax}, 
                                 AspectRatio -> 1,ContourStyle->cst]
     , {i, 1, nb}];
     If[plot,
        GraphicsGrid[{gr}],
        Return[gr]
     ]
  ]

(*
***)


(****f* /GTFermiSurface
! NAME
!  GTFermiSurface
! AUTHOR
!  W. Hergert
! PACKAGE
!  ElectronicStructure.m 
! MODIFICATION HISTORY
!  * 21.12.2016 : first implementation
!  * 28.12.2017 : check header and documentation  
!  * 24.06.2018 : check header and documentation  
!  * 08.02.2023 : labelling of k-points changed according to new Mathematica version 
!
! USAGE
!  GTFermiSurface[ham,fermi,nbands,ndel,kbasis,cluster,kpath]  calculates the Fermi surface.
! INPUT
!  * ham			- parametrized Hamiltonian
!  * fermi		- Fermi energy
!  * nbands		- list bands where parts of the Fermi surface are expected
!  * ndel		- number of points per spatial direction
!  * kbasis		- basis vectors of the reciprocal lattice
!  * cluster		- {cut,smin,smax} defines vectors for BZ construction
!  * kpath		- Path in BZ
!  
! OUTPUT
!  Plot of the Fermi surface with ListPlot3D
!
! GTPack OPTIONS
!  * GORegionFunction - can be used to cut parts of the Fermi surface that are outside of the BZ.
!  * GOBZ             - plot BZ
!  * GOBZPath         - plot BZPath
!  * GOVerbose:
!
!     - True  - additional information  
!     - False - no additional information (standard) 

! STANDARD OPTIONS
!  * VertexLabelStyle - Directive[Black, 20,Background -> Yellow]  style for the symmetry point labels.
! GTPack MODULES
!  GTVoronoiCell
! GTPack NOTEBOOKS
!  - 
! DESCRIPTION
!  Calculates the Fermi surface corresponding to a Hamiltonian and Fermi energy if the Fermi surface 
!  contains parts from list of bands. The electronic structure is calculated in a cube at ndel points 
!  per spatial dimension. kbasis is the basis of the reciprocal lattice. clusterdata contains the 
!  data for the lattice construction. The path used in electronic structure calculations can be given by BZpath.
! LITERATURE
!  -
! TODO
!  08.02.2023 : In GTFermiSurface and GTVoronoiCell Graphs are used for the k-paths. The representation of graphs
!               has change, thus problems appeared. A solution is found, but it is not optimal yet.
! RELEASE
!  1.0.0
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)


GTFermiSurface[ham_, fermi_,  nbands_, ndel_, kbasis_,cluster_, kpath_, OptionsPattern[]] := Module[
	{nb,del,verb,bzpath,bz,region,zone,sel,tab,ib,hc,ew,ewl,plot,bnr,band,fs,plt1,lstyle},
  nb = Length[nbands];del = 2/(ndel - 1);
  (*--- interpretation of options ---*)
  verb    = OptionValue[GOVerbose];
  bzpath  = OptionValue[GOBZPath];
  bz      = OptionValue[GOBZ];
  region  = OptionValue[GORegionFunction];
  lstyle  = OptionValue[VertexLabelStyle];
  (*--- calculate Brillouin zone ---*)
  If[bz || region,
     zone = GTVoronoiCell[kbasis, cluster, kpath, GOVerbose -> verb, GOBZPath -> bzpath, GOOutput -> "Output",VertexLabelStyle->lstyle],
     zone = {}
  ];
  (*--- calculate band structure data for plot ---*)
  If[verb,
     Print["Calculation of eigenvalues"],
     None
  ];
  sel = Table[0, {Length[ham]}]; tab = {}; 
  Do[
  	 sel[[nbands[[ib]]]] = 1
  ,{ib, 1, nb}];
  Do[
     Do[
        Do[
           Clear[\[Xi], \[Eta], \[Zeta]];
           hc = ham /. {\[Xi] -> x, \[Eta] -> y, \[Zeta] -> z}; 
           ew = Sort[Eigenvalues[hc]];
           ewl = Pick[ew, sel, 1]; tab = Append[tab, ewl];
        , {x, -1, 1, del}]
    , {y, -1, 1, del}]
  , {z, -1, 1, del}];
  (*--- construct the plots ---*)
  If[verb,
     If[region,
     	Print["Construction of the plots: cutting away parts outside BZ takes a while"],
     	Print["Construction of the plots"]
     ],	
     None
  ];
  plot = {}; 
  bnr = Table["Band number :" <> ToString[nbands[[ib]]], {ib, 1, nb}];
  Do[
     band = Partition[#, ndel] & /@ Partition[Transpose[tab][[ib]], ndel*ndel];
     If[region,
        fs = ListContourPlot3D[band, 
        	                   Contours -> {fermi}, Mesh -> None, 
                               DataRange -> {{-1, 1}, {-1, 1}, {-1, 1}},
                               RegionFunction -> 
                               Function[{x, y, z}, {x, y, z} \[Element] ConvexHullMesh[zone[[1]]]]
                              ],
        fs = ListContourPlot3D[band, 
     	                       Contours -> {fermi}, Mesh -> None, 
                               DataRange -> {{-1, 1}, {-1, 1}, {-1, 1}}
                              ]
     ];
     If[bz,    
        If[bzpath, 
           fs = {bnr[[ib]], Show[HighlightMesh[ConvexHullMesh[zone[[1]]], {Style[1, {Thick, Yellow}], 
                 Style[2, Opacity[.3, Pink]]}], zone[[2]], zone[[3]], fs]
                },
           fs = {bnr[[ib]], Show[HighlightMesh[ConvexHullMesh[zone[[1]]], {Style[1, {Thick, Yellow}], 
                 Style[2, Opacity[.3, Pink]]}], zone[[2]], fs]
                }
        ],
        fs = {bnr[[ib]], Show[fs]}
     ];
     plot = Append[plot, fs]
  , {ib, 1, nb}]; 
  plt1 = plot // Transpose;
  Return[Grid[plt1]]
]
  
(*
***)




(*-------------------------- Attributes ------------------------------*)
 Attributes[GTBands]={Protected, ReadProtected}
(* Attributes[GTBandStructure]={Protected, ReadProtected}*)
 Attributes[GTBandsPlot]={Protected, ReadProtected}
 Attributes[GTCompatibility]={Protected, ReadProtected}
(* Attributes[GTDensityOfStates]={Protected, ReadProtected}*)
 Attributes[GTDensityOfStatesRS]={Protected, ReadProtected}
 Attributes[GTDensityOfStatesPlot]={Protected, ReadProtected}




End[] 

EndPackage[]
(*
***)
