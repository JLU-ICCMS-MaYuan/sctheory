(****m* /TwoAxisListPlot.m
!
! NAME
!   TwoAxisLiistPlot.m
! AUTHOR
!  M. Honeychurh
! MODIFICATION HISTORY
!  * 27.06.2018 : cleaned up the header 
!
! USAGE
!  TwoAxisListPlot is used in DOS plots
!
! MODULES
!
! --- two axis plots ---
!
!    o TwoAxisListPlot        - plots data sets givenby points 
!    o TwoAxisListLinePlot    - plots a line through the point
!
! DESCRIPTION
!   This is an external package from the internet  (see: http://library.wolfram.com/infocenter/MathSource/7507)
! LITERATURE
!  -
! TODO
!  -
! RELEASE
!  1.0.0
! PROBLEMS
!  -
!--------------------------------------------------------------------------------
!
***)


BeginPackage["GroupTheory`ThirdParty`"]
(* Exported symbols added here with SymbolName::usage *) 
TwoAxisListPlot::usage = 
  "TwoAxisListPlot[{data 1,data 2},opts] rescales the \
second data set and plots both data sets on the same plot with the \
ticks for the first data set on the left and those for the second \
data set on the right. Ticks are automatically colored to match the \
data set they represent. The function takes any option that can be \
used with ListPlot.";

TwoAxisListLinePlot::usage = 
  "TwoAxisListLinePlot[{data 1,data 2},opts] rescales the \
second data set and plots both data sets on the same plot with the \
ticks for the first data set on the left and those for the second \
data set on the right. Ticks are automatically colored to match the \
data set they represent. The function takes any option that can be \
used with ListLinePlot.";

plotGrid::usage = 
   "PlotGrid[l_list, width,height] helps to plot combined figures \
l_list - list of the plots, width  - width of the individual plot,  \
height - height of the individual plot";


Begin["`Private`"]
(* Implementation of the package *)

Unprotect[FindDivisions,TwoAxisListPlot,TwoAxisListLinePlot,plotGrid];

FindDivisions = 
 If[$VersionNumber < 6.1, Developer`FindDivisions, 
  FindDivisions]; 
  
  setcol[styl : _RGBColor | _Hue | _CMYKColor | _GrayLevel] := {styl}

setcol[styl_List] := 
 With[{col = 
    Select[styl, 
     MemberQ[{RGBColor, Hue, CMYKColor, GrayLevel}, Head[#]] &]}, 
  If[col =!= {}, {col}, col]]

setcol[_] := {}

rools2 = {{a : _RGBColor | _Hue | _CMYKColor | _GrayLevel List, {}} \
:> {a, ColorData[1][2]}, {{}, 
     b : _RGBColor | _Hue | _CMYKColor | _GrayLevel} :> {ColorData[1][
      1], b}, {{}, {}} :> {ColorData[1][1], ColorData[1][2]}};

frameproc[opts__] := Module[{styls},
  styls = Select[{opts}, #[[1]] === PlotStyle &][[1, 2]];
  If[Length[styls] == 1, styls = Join[styls, styls], Take[styls, 2]];
  Map[setcol, styls] /. rools2
  ]
  
  Options[TwoAxisListPlot] = {AlignmentPoint -> Center, 
   AspectRatio -> 1/GoldenRatio, Axes -> True, AxesLabel -> None, 
   AxesOrigin -> Automatic, AxesStyle -> {}, Background -> None, 
   BaselinePosition -> Automatic, BaseStyle -> {}, 
   ClippingStyle -> None, ColorFunction -> Automatic, 
   ColorFunctionScaling -> True, ColorOutput -> Automatic, 
   ContentSelectable -> Automatic, 
   CoordinatesToolOptions -> Automatic, DataRange -> Automatic, 
   DisplayFunction :> $DisplayFunction, Epilog -> {}, Filling -> None,
    FillingStyle -> Automatic, FormatType :> TraditionalForm, 
   Frame -> True, FrameLabel -> None, 
   FrameStyle -> {{ColorData[1][1], ColorData[1][2]}, Automatic}, 
   FrameTicks -> Automatic, FrameTicksStyle -> {}, GridLines -> None, 
   GridLinesStyle -> {}, ImageMargins -> 0.`, ImagePadding -> All, 
   ImageSize -> Automatic, InterpolationOrder -> None, 
   Joined -> False, LabelStyle -> {}, MaxPlotPoints -> \[Infinity], 
   Mesh -> None, MeshFunctions -> {#1 &}, MeshShading -> None, 
   MeshStyle -> Automatic, Method -> Automatic, 
   PerformanceGoal :> $PerformanceGoal, PlotLabel -> None, 
   PlotMarkers -> None, PlotRange -> Automatic, 
   PlotRangeClipping -> True, PlotRangePadding -> Automatic, 
   PlotRegion -> Automatic, PlotStyle -> Automatic, 
   PreserveImageOptions -> Automatic, Prolog -> {}, 
   RotateLabel -> True};

Options[TwoAxisListLinePlot] = {AlignmentPoint -> Center, 
   AspectRatio -> 1/GoldenRatio, Axes -> True, AxesLabel -> None, 
   AxesOrigin -> Automatic, AxesStyle -> {}, Background -> None, 
   BaselinePosition -> Automatic, BaseStyle -> {}, 
   ClippingStyle -> None, ColorFunction -> Automatic, 
   ColorFunctionScaling -> True, ColorOutput -> Automatic, 
   ContentSelectable -> Automatic, 
   CoordinatesToolOptions -> Automatic, DataRange -> Automatic, 
   DisplayFunction :> $DisplayFunction, Epilog -> {}, Filling -> None,
    FillingStyle -> Automatic, FormatType :> TraditionalForm, 
   Frame -> True, FrameLabel -> None, 
   FrameStyle -> {{ColorData[1][1], ColorData[1][2]}, Automatic}, 
   FrameTicks -> Automatic, FrameTicksStyle -> {}, GridLines -> None, 
   GridLinesStyle -> {}, ImageMargins -> 0.`, ImagePadding -> All, 
   ImageSize -> Automatic, InterpolationOrder -> None, Joined -> True,
    LabelStyle -> {}, MaxPlotPoints -> \[Infinity], Mesh -> None, 
   MeshFunctions -> {#1 &}, MeshShading -> None, 
   MeshStyle -> Automatic, Method -> Automatic, 
   PerformanceGoal :> $PerformanceGoal, PlotLabel -> None, 
   PlotMarkers -> None, PlotRange -> Automatic, 
   PlotRangeClipping -> True, PlotRangePadding -> Automatic, 
   PlotRegion -> Automatic, PlotStyle -> Automatic, 
   PreserveImageOptions -> Automatic, Prolog -> {}, 
   RotateLabel -> True};

(****w* /TwoAxisListPlot
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

TwoAxisListPlot[{data1_List, data2_List}, opts : OptionsPattern[]] := 
 Module[{pltopts, d2vals, d1mn, d1mx, d2mn, d2mx, newdata, rescl, 
   tcklocsm, tcklocsM, d2tcks},
  pltopts = 
   If[Select[{opts}, #[[1]] === PlotStyle &] === {}, 
    Join[{opts}, {FrameStyle -> {{ColorData[1][1], ColorData[1][2]}, 
        Automatic}}], 
    Join[{opts}, {FrameStyle -> {frameproc[opts], Automatic}}]];
  d2vals = If[MatrixQ[data2], data2[[All, 2]], data2];
  {d1mn, d1mx} = 
   If[MatrixQ[data1], {Min[data1[[All, 2]]], 
     Max[data1[[All, 2]]]}, {Min[data1], Max[data1]}];
  {d2mn, d2mx} = {Min[d2vals], Max[d2vals]};
  rescl[num_] := Rescale[num, {d2mn, d2mx}, {d1mn, d1mx}];
  newdata = 
   If[MatrixQ[data2], 
    Transpose[{data2[[All, 1]], rescl[#] & /@ d2vals}], 
    rescl[#] & /@ d2vals];
  {tcklocsM, tcklocsm} = 
   FindDivisions[{d2mn, d2mx} + ((d2mx - d2mx) .05) {-1, 1}, {10, 
     5}];
  d2tcks = 
   Join[MapThread[{rescl[#1], 
       NumberForm[N[#2], {5, 2}], {0.00675`, 
        0.`}, {AbsoluteThickness[0.125`]}} &, {tcklocsM, 
      tcklocsM}], {rescl[#], 
       "", {0.00375`, 
        0.`}, {AbsoluteThickness[0.125`]}} & /@ (Flatten[
       tcklocsm[[All, 2 ;; -2]]])];
  ListPlot[{data1, newdata}, 
   FrameTicks -> {{Automatic, d2tcks}, {Automatic, Automatic}}, 
   Evaluate@pltopts, Options[TwoAxisListPlot]]
  
  ]
  
(*
***) 
 
(****w* /TwoAxisListLinePlot
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
  
  TwoAxisListLinePlot[{data1_List, data2_List}, 
  opts : OptionsPattern[]] := 
 Module[{pltopts, d2vals, d1mn, d1mx, d2mn, d2mx, newdata, rescl, 
   tcklocsm, tcklocsM, d2tcks},
  pltopts = 
   If[Select[{opts}, #[[1]] === PlotStyle &] === {}, 
    Join[{opts}, {FrameStyle -> {{ColorData[1][1], ColorData[1][2]}, 
        Automatic}}], 
    Join[{opts}, {FrameStyle -> {frameproc[opts], Automatic}}]];
  d2vals = If[MatrixQ[data2], data2[[All, 2]], data2];
  {d1mn, d1mx} = 
   If[MatrixQ[data1], {Min[data1[[All, 2]]], 
     Max[data1[[All, 2]]]}, {Min[data1], Max[data1]}];
  {d2mn, d2mx} = {Min[d2vals], Max[d2vals]};
  rescl[num_] := Rescale[num, {d2mn, d2mx}, {d1mn, d1mx}];
  newdata = 
   If[MatrixQ[data2], 
    Transpose[{data2[[All, 1]], rescl[#] & /@ d2vals}], 
    rescl[#] & /@ d2vals];
  {tcklocsM, tcklocsm} = 
   FindDivisions[{d2mn, d2mx} + ((d2mx - d2mx) .05) {-1, 1}, {10, 
     5}];
  d2tcks = 
   Join[MapThread[{rescl[#1], 
       NumberForm[N[#2], {5, 2}], {0.00675`, 0.`}, {GrayLevel[0.`], 
        AbsoluteThickness[0.125`]}} &, {tcklocsM, 
      tcklocsM}], {rescl[#], 
       "", {0.00375`, 0.`}, {GrayLevel[0.`], 
        AbsoluteThickness[0.125`]}} & /@ (Flatten[
       tcklocsm[[All, 2 ;; -2]]])];
  ListLinePlot[{data1, newdata}, 
   FrameTicks -> {{Automatic, d2tcks}, {Automatic, Automatic}}, 
   Evaluate@pltopts, Options[TwoAxisListLinePlot]] 
  ]

(*
***) 

(****h* /plotGrid
! NAME
!  plotGrid
! AUTHOR
!  Jens
! PACKAGE
!  ThirdParty.m
! MODIFICATION HISTORY
!  * 22.05.2020 : first test version in test package
! USAGE
!  plotGrid helps to plot combined figures
! INPUT
!	* l_list - list of the plots
!	* width  - width of the individual plot
!	* height - height of the individual plot	
! OUTPUT
!  list with weighting factors and energy derivatives
! GTPack OPTIONS
!  -
! Standard OPTIONS
!  *
! GTPack MODULES
!  -
! DESCRIPTION
!   
! LITERATURE
!  Source is Mathematica stackexchange
!   https://mathematica.stackexchange.com/a/6882/53872
! TODO
!  
! PROBLEMS
!  
! SOURCE
!--------------------------------------------------------------------------------
!
*)

Options[plotGrid] = {ImagePadding -> 40};
plotGrid[l_List, w_, h_, opts : OptionsPattern[]] := 
 Module[{nx, ny, sidePadding = OptionValue[plotGrid, ImagePadding], 
   topPadding = 0, widths, heights, positions, i,j,
   frameOptions = 
    FilterRules[{opts}, 
     FilterRules[Options[Graphics], 
      Except[{ImagePadding, Frame, FrameTicks}]]]}, {ny, nx} = 
   Dimensions[l];
  widths = (w - 2 sidePadding)/nx Table[1, {nx}];
  widths[[1]] = widths[[1]] + sidePadding;
  widths[[-1]] = widths[[-1]] + sidePadding;
  heights = (h - 2 sidePadding)/ny Table[1, {ny}];
  heights[[1]] = heights[[1]] + sidePadding;
  heights[[-1]] = heights[[-1]] + sidePadding;
  positions = 
   Transpose@
    Partition[
     Tuples[Prepend[Accumulate[Most[#]], 0] & /@ {widths, heights}], 
     ny];
  Graphics[
   Table[Inset[
     Show[l[[ny - j + 1, i]], 
      ImagePadding -> {{If[i == 1, sidePadding, 0], 
         If[i == nx, sidePadding, 0]}, {If[j == 1, sidePadding, 0], 
         If[j == ny, sidePadding, topPadding]}}, AspectRatio -> Full],
      positions[[j, i]], {Left, Bottom}, {widths[[i]], 
      heights[[j]]}], {i, 1, nx}, {j, 1, ny}], 
   PlotRange -> {{0, w}, {0, h}}, ImageSize -> {w, h}, 
   Evaluate@Apply[Sequence, frameOptions]]]

(*
**)

End[]

Protect[FindDivisions,TwoAxisListPlot,TwoAxisListLinePlot,plotGrid];

EndPackage[]

