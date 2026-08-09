(****m* /SimPack.m
!
! NAME
!  SimPack.m
! AUTHOR
!  W. Hergert, M. Geilhufe
! MODIFICATION HISTORY
! * 26.10.2012 : initially created and documented 
! * 20.05.2018 : RoboDoc checked 
! * 29.06.2022 : package more or less completed and headers checked
!
! USAGE
!  Modules used to connect GTPack to external FORTRAN codes
!
! GTPack MODULES
!
! --- import and export of parameter sets
!
!   * GTTbParmExport  - exports a parameter set to be used in SimPack
!   * GTTbParmImport  - reads a parameter set and gives it in the usual form of a parameter set
!
! --- import results of calculations
!
!   * GTTbReadBands    - Read bands from Fortran calculation
!   * GTTbReadDOS      - Read DOS from Fortran calculation
!
! --- other things ---
!
!   * GTTbFitExport    - exports TB band structure data to be used in fitting procedures tests
!   * GTBZTbPointMesh  - export k-mesh to be used in Fortran
!   
! 
! DESCRIPTION
!  In some cases the solution of problems with TB-Hamiltonians seems to be complicated in the
!  framework of GTPack. Thus, a connection to external Fortran codes is established by the modules
!  in this package. External tasks can be:
!   - band structure calculations
!   - calculation of DOS by means of smearing and tetrahedron method
!   - fitting of TB-models to ab initio band structures
!
! LITERATURE
! 
! TODO
!
! PROBLEMS
!
!--------------------------------------------------------------------------------
!
***)

BeginPackage["GroupTheory`SimPack`",{"GroupTheory`Basic`","GroupTheory`Lattice`","GroupTheory`Auxiliary`","GroupTheory`ElectronicStructure`","GroupTheory`Symbols`"}]

  GTBZTbPointMesh             ::usage=  "GTBZTbPointMesh[\*StyleBox[\"file, kpts, w, dig\", \"TI\"]] exports a k-mesh \*StyleBox[\"kpts\",\"TI\"] to \*StyleBox[\"file\",\"TI\"] with the format of reals \*StyleBox[\"Ew.d\",\"TI\"]."
  GTTbParmExport              ::usage=  "GTTbParmExport[\*StyleBox[\"parmset, file, info\", \"TI\"]]  exports a parameter set \*StyleBox[\"parmset\",\"TI\"] to \*StyleBox[\"file\",\"TI\"] with additional information \*StyleBox[\"info\",\"TI\"]."
  GTTbFitExport               ::usage=  "GTTbFitExport[\*StyleBox[\"file, bands\", \"TI\"]] exports TB band structure \*StyleBox[\"bands\",\"TI\"] to \*StyleBox[\"file\",\"TI\"] to be used to test fitting procedures."
  GTTbParmImport              ::usage=  "GTTbParmImport[\*StyleBox[\"file, parmset\", \"TI\"]] reads a parameter set \*StyleBox[\"parmset\",\"TI\"]  from \*StyleBox[\"file\",\"TI\"]."
  GTTbReadBands               ::usage=  "GTTbReadBands[\*StyleBox[\"file,kmesh\", \"TI\"]] reads band structure calculated by Fortran code from \*StyleBox[\"file\",\"TI\"] that corresponds to \*StyleBox[\"kmesh\",\"TI\"]."
  GTTbReadDOS                 ::usage=  "GTTbReadBands[\*StyleBox[\"file,names\", \"TI\"]] reads DOS calculated by Fortran code from \*StyleBox[\"file\",\"TI\"]."

  Options[GTBZTbPointMesh]  = {GOVerbose -> True,GOOutput->"Bands"}
  Options[GTTbFitExport]    = {GOVerbose -> True, GORestrict -> {}}
  Options[GTTbParmImport]   = {GOVerbose -> False, GODecimals -> 5}
  Options[GTTbReadBands]    = {GOVerbose -> True}
  Options[GTTbReadDOS]      = {GOVerbose -> True, GOPlotDos -> "DOS",PlotRange -> All, FrameLabel -> {"Energy", "DOS"}, 
                               PlotLabel -> "Density of States", PlotStyle -> {Blue, Orange},Filling -> None}


Begin["`Private`"] 




(****o* /GTTbFitExport
! NAME
!  GTTbFitExport
! AUTHOR
!  W. Hergert
! PACKAGE
!  SimPack.m
! MODIFICATION HISTORY
!   * January 2015  : first version
!   * 29.06.2022    : check header
!   * 27.01.2023    : Option GOWeight renamed to GORestrict, because it is here a restriction 
!                     of the number of bands. 
!                     Additional weight factors can be defined in the fit itself.
! USAGE
!  GTFitExport exports bands structures to be used in tests of fitting procedurs
! INPUT
!  * file    - name of the file for export
!  * bands   - band structure calculate by means of GTPack (GTBands)
! OPTIONS
!
!  * GOVerbose 
!
!			"True"  - additional information (standard)
!			"False" - no information
!
!  * GORestricts   - select bands for the fit 
!           {}      - no selection (standard)
!
! OUTPUT
!  a bands structure file that can be used to test fit procedures in SimPack
! 
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  In the Fortran code different methods are implemented (LSQR, simmplex according to Nelder-Mead, genetic algorithm)
!  Sometimes it will be necessary to test the fit with synthetic data.
! LITERATURE
!  -
! TODO
!  The export and use of weights has to be checked carefully!
! PROBLEMS
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTbFitExport[file_, bands_, OptionsPattern[]] := Module[
	{verb,weight,nbands, nkp,wfac,exp,temp,i},
  (*--- options ---*)
     verb   = OptionValue[GOVerbose];
     weight = OptionValue[GORestrict];
  (*--- prepare the data ---*)
     nbands = bands[[1, 2]] // Length;
     nkp = bands // Length;
     If[verb,
        Print["number of k-points = ", nkp];
        Print["number of bands    = ", nbands],
        None
     ];
     If[weight == {},
       wfac = Table[nbands, {nkp}],
       wfac = weight
     ];
     exp = {nkp, nbands};
     Do[
        temp = {{wfac[[i]]}, bands[[i]]*1.} // Flatten; 
        exp  = Append[exp, temp]
     , {i, 1, nkp}];
     Export[file, exp, "Table"]
]

(*
***) 
  

(****o* /GTTbParmImport
! NAME
!  GTTbParmImport
! AUTHOR
!  W. Hergert
! PACKAGE
!  SimPack.m
! MODIFICATION HISTORY
!   January 2015  : first version
!   29.07.2022    : check of header
!   27.01.2023    : new version with output of a table 
! USAGE
!  GTTbParmImport is used to import parameter sets created as a result of the fitting procedures
! INPUT
!  * file    - name of the file contaiming the fitted parameters
!  * parmset - a GTPack parameter set having the same structure as the fitted one
!
! OPTIONS
!
!  * GOVerbose 
!
!			"True"  - additional information (standard)
!			"False" - no information
!
! OUTPUT
!  parameterset in GTPack form
! 
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
! -
! DESCRIPTION
!  The command reads parameter sets which result from fits to ab initio or some test bandstructures.
!  Getting the parameters back to GTPack allows further calculations.
!
! LITERATURE
!  -
! TODO
!  -
! PROBLEMS
! -
!--------------------------------------------------------------------------------
! SOURCE
*)

(*
GTTbParmImport[file_, parms_, OptionsPattern[]] :=Module[
	{nparm, names, verb,pnew,nparm1},
  (*--- options ---*)
     verb  = OptionValue[GOVerbose];
     names = Transpose[parms][[1]];
     nparm = Length[parms];
     pnew  = ReadList[file, Number];
     If[verb,
        Print[pnew],
        None
     ];
     nparm1 = Length[pnew];
     If[nparm1 > nparm,
        Print["Error Length paramter vectors!"]; 
        Abort[],
        None
     ];
     pnew = {pnew, Table[0, {nparm - nparm1}]} // Flatten;
     pnew = {names, pnew} // Transpose;
     Return[pnew]
]
*)

GTTbParmImport::parm = "Number of the parameters in the two sets disagree."

GTTbParmImport[file_, parms_, OptionsPattern[]] := Module[
	{verb,nround,names,nnames,pnew,nparm,i,no,n1,plt},
  (*---options---*)
  verb   = OptionValue[GOVerbose];
  nround = OptionValue[GODecimals];
  (* dataset for names *)
  names = Transpose[parms][[1]];
  nnames = Length[names];
  (* new data *)
  pnew = ReadList[file, Number];
  nparm = Length[pnew];
  (* check *)
  If[nparm == nnames,
  	 None,
     Message[GTTbParmImport::parm]; Return[]    
  ];
  (* number of decimals *)
  If[nround == 0,
   None,
   Do[
   	  pnew[[i]] =ToString[PaddedForm[N[pnew[[i]]], {nround, nround}]]
     , {i, 1, nparm}]
   ];
  (* print  nice table *)
  If[verb,
     no = IntegerPart[nnames/10]; 
     n1 = nnames - no*10;
     Do[
   	    plt = {Take[names, {(i - 1)*10 + 1, i*10}], 
        Take[pnew, {(i - 1)*10 + 1, i*10}]};
        Print[Grid[plt, Frame -> All, 
              Dividers -> {Black, {2 -> GTDividerColor1}}, 
              Background -> {None, {1 -> GTBackGroundColor1}}]]
     , {i, 1,no}];
     plt = {Take[names, {no*10 + 1, nnames}], Take[pnew, {no*10 + 1, nnames}]};
     Print[Grid[plt, Frame -> All, 
           Dividers -> {Black, {2 -> GTDividerColor1}}, 
           Background -> {None, {1 -> GTBackGroundColor1}}]],
     None
  ];
  (* output of dataset *)
  pnew = {names, pnew} // Transpose;
  Return[pnew]
]

(*
***) 



(****o* /GTTbParmExport
! NAME
!  GTTbParmExport
! AUTHOR
!  W. Hergert
! PACKAGE
!  Simpack.m
! MODIFICATION HISTORY
!  * 15.01.2015 : first version
!  * 29.06.2022 : implemented finally as GTTbParmExport (priviously also GTTbParmToFortran)
! USAGE
!  GTTbParmExport exports a tight-inding parameter set in a form useful for other
!  external packages
! INPUT
!  * parmset - a parameter set from a database (use GTTbDatabaseRetrieve)
!  * file    - name of the file for export
!  * info    - a comment line is written at the beginning of the dataset
!  * ham     - Hamiltonian (optional argument) 
! OUTPUT
!  * paramter sets readeable from FORTRAN code in SimPack
! GTPack OPTIONS
!  -
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  Mathematica allows to construct names for the tight-binding parameters,
!  which are close to the nomenclature used in the literature. Such nameing
!  conventions are not very useful in numeric programs written in FORTRAN.
!  The names are transformed in an appropriate form.
!  In some cases you might have a parameter set, but the actual Hamiltonian
!  does not contain all parameters. If the Hamiltonian is provided only the necessary 
!  parameters will be exported.
!
! LITERATURE
!  -
! TODO
!  seems to work properly
! PROBLEMS
! nothing to do in the moment
!
!--------------------------------------------------------------------------------
! SOURCE
*)

GTTbParmExport[pset_, file_, info_, ham_: {}] := Module[
	{parms0, vals0, np, parms, vals, rule, i, pos, vp, pp, s, parm, txt, n1, no, plt},
  (*--- select parameters ---*)
    parms0 = Transpose[pset][[1]];
    vals0  = Transpose[pset][[2]];   
  (*--- ham not given ---*)  
    If[ham == {},
       parms = parms0;
       vals  = vals0;
       np    = Length[parms],
  (*--- ham is given ---*)     
       rule  = Select[parms0, Not[FreeQ[ham, #]] &];
       parms = {}; 
       vals  = {};
       np    = Length[rule];
       Do[
          pos = Position[parms0, rule[[i]]] // Flatten;
          If[pos == {},
             vp = 0,
             vp = vals0[[pos[[1]]]];
             pp = parms0[[pos[[1]]]]
          ];
          parms = Append[parms, pp];
          vals  = Append[vals, vp]
       , {i, 1, np}]
    ];
  (*--- Export the paramter set ---*)
    parms = parms // Flatten;
    s = OpenWrite[file, FormatType -> StandardForm];
    Write[s, info];
    Do[
    	   parm = parms[[i]]; 
    	   txt  = parm;
       If[Head[parm] == Subscript, 
          txt = StringJoin[parm[[1]], "_", ToString[parm[[2]]]], 
          None
       ];
       If[Head[parm] == Subsuperscript, 
          txt = StringJoin[parm[[1]], "_", ToString[parm[[2]]], "^(", parm[[3]], ")"], 
          None
       ];
       If[Head[parm] == Superscript, 
          txt = StringJoin[parm[[1]], "^(", parm[[2]], ")"], 
          None
       ];
       txt = Characters[txt] /. {"\[Sigma]" -> "s", "\[Pi]" -> "p", "\[Delta]" -> "d", "," -> ";"};
       txt = "'" <> txt <> "'";
       Write[s, txt, "    ", vals[[i]]]
    , {i, 1, np}];
    Close[s];
  (*--- Protocol ---*)
    no = IntegerPart[np/10]; n1 = np - no*10;
    Do[
    	   plt = {Table[k, {k, (i - 1)*10 + 1, i*10}], Take[parms, {(i - 1)*10 + 1, i*10}], Take[vals, {(i - 1)*10 + 1, i*10}]};
       Print[Grid[plt, Frame -> All, Dividers -> {Black, {2 -> GTDividerColor1}}, Background -> {None, {1 -> GTBackGroundColor1, 
                  2 -> LightYellow}}
                 ]
            ]
    , {i, 1, no}];
     plt = {Table[k, {k, no*10 + 1, np}], Take[parms, {no*10 + 1, np}], Take[vals, {no*10 + 1, np}]};
     Print[Grid[plt, Frame -> All, Dividers -> {Black, {2 -> GTDividerColor1}}, Background -> {None, {1 -> GTBackGroundColor1, 
     	             2 -> LightYellow}}
     	       ]
     	  ];
     Return[{parms,vals}//Transpose]	  
  ]
(*
***)

(****o* /GTBZTbPointMesh
! NAME
!  GTBZTbPointMesh
! AUTHOR
!  W. Hergert
! PACKAGE
!  SimPack.m 
! MODIFICATION HISTORY
!  * 01.04.2022 - first version
!  * 28.06.2022 - introduced in SimPack 
! 
! USAGE
!  GTBZTbPointMesh[file,kpts,w,dig] exports a kmesh for calculations with the FORTRAN codes.
!
! INPUT
!   * file      - output file name
!   * kpoints   - list of kpoints
!   * w         - format of the REALs Ew.d
!   * d         - format of the REALs Ew.d
!  
! OUTPUT
!  data stored in file
! GTPack OPTIONS
!
!   GOVerbose
!     o True    - additional information (standard)
!     o False   - no additional information
!   GOOutput  
!     o "Bands" - band structure calculation (standard)
!     o "DOS"   - DOS calculation 
! STANDARD OPTIONS
!  -
! GTPack MODULES
!  - f77EForm - change REALs in format Ew.d  (internal module)
!    
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!   The modules GTBZLines and GTBZPointMesh are used to generate k-vectors along a line ore in the part of the
!   BZ. To be used in the FORTRAN codes the module exports the k-vectors in the appropriate format.
! LITERATURE
!  see for f77EForm : https://mathematica.stackexchange.com/questions/19387/how-to-export-data-files-using-specific-number-format
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

GTBZTbPointMesh::option = " `1` in not an allowed value for GOOutput." 


GTBZTbPointMesh[file_, kpts_,fw_Integer, ndig_Integer, OptionsPattern[]] := Module[
	{nkp,s,verb,opt,v}, 
     verb = OptionValue[GOVerbose];
     opt  = OptionValue[GOOutput];
     If[Intersection[{opt},{"Bands","DOS"}]=={},
     	Message[GTBZTbPointMesh::option];Return[],
     	None
     ];
     If[verb, 
        Print["Use k-points in FORTRAN-ETBM calculation for ",opt];
        Print["Output in file    : ", file]
     ];  
(* output for band structure *)     
     s = OpenWrite[file, FormatType -> StandardForm];
     If[opt=="Bands",
        nkp=Length[kpts[[1]]];
        If[verb, 
           Print["Number of k-points: ", nkp]
        ];       
        Write[s,nkp];
        Do[
           v=kpts[[1,i]]//Flatten;
     	   Write[s, i,"  ",f77Eform[v[[3]], fw,ndig], "  ",
     	   	               f77Eform[v[[4]], fw,ndig], "  ", f77Eform[v[[5]],fw,ndig]]
     , {i, 1, nkp}],
        None
     ];
     If[opt=="DOS",
        nkp=Length[kpts];
        If[verb, 
           Print["Number of k-points: ", nkp]
        ];       
        Write[s,nkp];
        Do[
           v=kpts[[i]];
     	   Write[s, i,"  ",f77Eform[v[[1]], fw,ndig], "  ",
     	   	               f77Eform[v[[2]], fw,ndig], "  ", f77Eform[v[[3]],fw,ndig]]
     , {i, 1, nkp}],
        None
     ];
     Close[s]
 ]

f77Eform[x_?NumericQ, fw_Integer, ndig_Integer] := 
 Module[{sig, s, p, ps}, {s, p} = MantissaExponent[x];
  {sig, ps} = {ToString[Round[10^ndig Abs[s]]], ToString[Abs[p]]};
  StringJoin @@ 
   Join[Table[" ", {fw - ndig - 7}], 
    If[x < 0, "-", " "], {"0."}, {sig}, 
    Table["0", {ndig - StringLength[sig]}], {"E"}, 
    If[p < 0, {"-"}, {"+"}], 
    Table["0", {2 - StringLength[ps]}], {ps}]]

(*
***)

(****k* /GTTbReadBands
! NAME
!  GTTbReadBands
! AUTHOR
!  W. Hergert
! PACKAGE
!  Simpack.m 
! MODIFICATION HISTORY
!  * 01.04.2022 : first version
!  * 29.06.2022 : check header
!
! USAGE
!  GTTbReadBands[file,kpts] is used to read band structures calculate with the Fortran code
!
! INPUT
!   *  file  - filename of the data
!   *  kpts  - k-vectors generated with GTBZLines or GTBZPointMesh
!  
! OUTPUT
!  a band structure to be directly used in GTBandsPlot
!
! GTPack OPTIONS
!
!  GOVerbose
!     o True  - additional information (standard)
!     o False - no additional information
!
! GTPack MODULES
!  -
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  GTPack can be used to plot the results of the external band structure calculation by  means of the Fortran codes.
!  Thus, this module reads the data in the appropriate form.
!  kpts has to be the kmesh data set for which the external calculation has been performed. The information in this
!  file is used to complete the plot.
!  A similar module is GTTbReadDOS.
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
GTTbReadBands::badfile =  "File `1` not in this directory.";
GTTbReadBands::kpoints =  "Disagreement in k-point number.";	


GTTbReadBands[file_,kpts_,OptionsPattern[ ]] := Module[
	{verb,fn,isl,s,nkp,nkpt1,ndim,tab,bands,v,i}, 
  (*--- options ---*)	
	verb = OptionValue[GOVerbose];
  (*--- check file names ---*)  
    fn   = FileNames[]; 
    isl  = Intersection[fn, {file}];
    If[isl == {}, 
       Message[GTTbReadBands::badfile,file]; Return[], 
       s=OpenRead[file]
    ];    
  (*--- read and prepare the data ---*)
      nkp=kpts[[1]]//Length; 	
      If[verb,
      	 Print["Number of k-points : ",nkp],
      	 None
      ];
      {nkpt1,ndim}=Read[s,{Number, Number}];
      tab=Table[Number,ndim+1];
      If[nkp==nkpt1,
      	 None,
      	 Message[GTTbReadBands::kpoints]; Return[]
      ];
      bands={};
      Do[
      	 v=Drop[Read[s,tab],1];
      	 bands=Append[bands,{i,kpts[[1,i,2]],kpts[[1,i,3]],v}]
      ,{i,1,nkp}];
      bands=Append[{bands},kpts[[2]]];  	  
      Return[bands]
    ]
    
    
(*
***)



(****k* /GTTbReadDOS
! NAME
!  GTTbReadDOS
! AUTHOR
!  W. Hergert
! PACKAGE
!  Simpack.m 
! MODIFICATION HISTORY
!  * 01.04.2022 : first version
!  * 10.05.2022 : revised version after completion of FORTRAN code
!  * 29.06.2022 : check of the header
!  * 25.01.2023 : calculation of Ef removed because it comes from SimPack already
!  * 04.02.2023 : revised version with better structure
!
! USAGE
!  GTTbReadDOS[file,sum] is used to read density of states calculated with the FORTRAN codes.
!  The data can be directy plotted.
!
! INPUT
!   *  file  - filename of the data
!   *  names - names for the partial densities of states
!               
!  
! OUTPUT
!   either the Dos data read in or a plot of the DOS
!
! GTPack OPTIONS
!
!  GOVerbose
!     o True  - additional information (standard)
!     o False - no additional information
!
!  GOPlotDOS
!     o DOS, IDOS,ALL,PDOS,DOS+IDOS - plots different selection of DOS
!     o NONE                        - no plots
!
! STANDARD OPTIONS
!   PlotRange  - All
!   PlotLabel  - "Density of states"
!   FrameLabel - {"Energy","DOS"}
!   PlotStyle  - {Blue,Orange}
!   Filling   -  {}
!     
! GTPack MODULES
!   dr - module to read a DOS (internal module)
! GTPack NOTEBOOKS 
!  -
! DESCRIPTION
!  GTPack will be used to read and plot the results of a DOS calculation  by  means of the FORTRAN codes.
!  First, the energies for the DOS calculation are calculated by means of the Hamiltonian generated in GTPack 
!  and exported to FORTRAN with GTTbParmExport. Next the DOS is calculated with FORTRAN implementation
!  of Gaussian smearing or the tetrahedron method in SimPack.
!
!  GTTbReadDOS reads the results in the appropriate form.
!
!  The DOS file is strucured as following: 
!  3 lines of control data: 
!  - total DOS
!  - integrated DOS
!  - total DOS
!  - partial DOSes
!
!  A similar module is GTTbReadBands.

! LITERATURE
!  see Documentation of SimPack
! TODO
!  -
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)


GTTbReadDOS::badfile = "File `1` not in this directory.";
GTTbReadDOS::plot    = "Number of densities `1` : not enough data to plot PDOS.";
GTTbReadDOS::badoption = " `1` for GOPlotDos is not correct. It has to be either DOS, IDOS, ALL, or DOS.";

GTTbReadDOS[file_, names_, OptionsPattern[]] :=
 Module[{verb, pdos, frl, plt, pll, plst, fill, fn, isl, s, plo,dd,dd1,emin, de, emax, ef, nel, np, nenin, energ, tab, nl, 
  dos, temp,i,nam,func,eln,nam1,efl,plg,plg1,g1,g2,lef},
   plo = {"DOS", "IDOS", "ALL", "PDOS", "NONE", "DOS+IDOS"};
   lef=True;
 (*--- options ---*)
   verb = OptionValue[GOVerbose];
   pdos = OptionValue[GOPlotDos];
   frl  = OptionValue[FrameLabel];
   plt  = OptionValue[PlotRange];
   pll  = OptionValue[PlotLabel];
   plst = OptionValue[PlotStyle];
   fill = OptionValue[Filling];
 (*--- check file name ---*)
   fn   = FileNames[];
   isl  = Intersection[fn, {file}];
   If[isl == {},
      Message[GTTbReadDOS::badfile, file]; Return[],
      s = OpenRead[file]
   ];
 (* check plot option *)
   isl = Intersection[plo, {pdos}];
   If[isl == {},
      Message[GTTbReadDOS::badoption, pdos]; Return[],
      s = OpenRead[file]
   ];
 (*--- Read DOS ---*)
   {emin, de, emax, ef, nel} = Read[s, {Number, Number, Number, Number, Number}];
   If[ef == emax || ef < 0,
      ef = "not calc.";
      lef=False,
      
      None
   ];
 (*--- number of densities, numbr of energy values ---*)
   {np, nenin} = Read[s, {Number, Number}];
   If[verb,
     tab = {{"emin", "de ", "emax"," Fermi energy \!\(\*SubscriptBox[\(E\), \(F\)]\)", "electrons"}, 
     	     {emin, de, emax, ef, nel}};
     Print[Grid[tab, Frame -> All, Background -> {None, {1 -> GTBackGroundColor1}}]];
     Print["Number of densities ", np, "  energy points = ", nenin],
     None
   ];  
 (*--- check consistency ---*) 
   If[np == 2 && (pdos == "PDOS"|| pdos=="ALL"),
      Message[GTTbReadDOS::plot,np];
      pdos = "DOS",
      None
   ];
 (*--- generate energy axis ---*)
   energ = Table[emin + (i - 1)*de, {i, 1, nenin}];
 (*--- read all DOSes ---*)
   nl  = IntegerPart[nenin/8.];
   tab = Table[Number, {8}];
   dos = {};
   Do[
      temp = dosin[nl, nenin, tab, s];
      temp = {energ, temp} // Transpose;
      dos  = Append[dos, temp]
   , {i, 1, np}];
   Close[s];
 (*--- total electron number ---*)
   nam = {"\!\(\*SubscriptBox[\(E\), \(F\)]\)", "DOS", "\!\(\*SubscriptBox[\(N\), \(el\)]\)"};
   If[np == 2&&lef,
      func = Interpolation[dos[[2]]];
      eln  = func[ef];
      func = Interpolation[dos[[1]]];
      dd   = func[ef];
      tab  = {nam, {"total", eln, dd}} // Transpose;
      Print[Grid[tab, Frame -> All, Background -> {None, {1 -> GTBackGroundColor1}, {1, 1} -> Yellow}]],
      None
   ];
 (*---partial DOSes---*)
   If[np > 3&&lef,
      dd   = {"DOS"};
      eln  = {"\!\(\*SubscriptBox[\(N\), \(el\)]\)"};
      func = Interpolation[dos[[2]]];
      eln  = Append[eln, func[ef]];
      func = Interpolation[dos[[1]]];
      dd   = Append[dd, func[ef]];
      Do[
         func = Interpolation[dos[[i + 3]]];
         dd   = Append[dd, func[ef]];
         eln  = Append[eln, Integrate[func[x], {x, emin, ef}]]
      , {i, 1, np - 3}];
      nam1 = {"\!\(\*SubscriptBox[\(E\), \(f\)]\)", "total", names} // Flatten;
      tab  = {nam1, dd, eln} ;
      Print[Grid[tab, Frame -> All, Background -> {None, {1 -> GTBackGroundColor1}, {1, 1} -> Yellow}]],           
      None
   ];
 (* --- Plots ---*)
   If[lef,
      efl = Graphics[{Thick, Red, Line[{{ef, -1}, {ef, 1000}}]}],
      efl = Graphics[{Black, Line[{{emax, -1}, {emax, 1000}}]}]
   ];
   Switch[pdos,
     "NONE",
       Return[dos],
     "DOS",
       dd = dos[[1]]; plg = "DOS",
     "IDOS",
       dd = dos[[2]]; plg = "IDOS",
     "DOS+IDOS",
       dd = dos[[1 ;; 2]]; plg = {"DOS", "IDOS"},
     "PDOS",
       dd = dos[[3 ;; np]]; plg = Prepend[names, "DOS"] // Flatten,
     "ALL",
       dd = dos[[3 ;; np]]; plg = Prepend[names, "DOS"] // Flatten;
      dd1 = dos[[1 ;; 2]]; plg1 = {"DOS", "IDOS"},
     _,
      Message[GTTbReadDOS::badoption, pdos]; Return[]
   ];
   If[pdos == "ALL",
      g1 = {ListPlot[dd1, Joined -> True, Frame -> True, PlotLegends -> plg1, FrameLabel -> frl, PlotLabel -> pll, 
            Filling -> fill, PlotStyle -> plst, PlotRange -> plt], efl};
      g2 = {ListPlot[dd, Joined -> True, Frame -> True, PlotLegends -> plg, FrameLabel -> frl, PlotLabel -> pll, 
            Filling -> fill, PlotStyle -> plst, PlotRange -> plt], efl};
      Grid[{{Show[g1], Show[g2]}}],
      Show[{ListPlot[dd, Joined -> True, Frame -> True,PlotLegends -> plg, FrameLabel -> frl, PlotLabel -> pll, 
            Filling -> fill, PlotStyle -> plst, PlotRange -> plt], efl}
          ]
   ]
]


dosin[nl_, nenin_, tab_, s_] := Module[{v, dd, k, tab1}, dd = {};
  Do[v = Read[s, tab];
   dd = Append[dd, v], {k, 1, nl}];
  If[nenin > nl*8, tab1 = Table[Number, {nenin - nl*8}];
   v = Read[s, tab1];
   dd = Append[dd, v], None];
  dd = dd // Flatten;
  Return[dd]]


End[]

EndPackage[]