(****m* /Additional.m
!
! NAME
!  Additional.m
! AUTHOR
!  W. Hergert
! MODIFICATION HISTORY
!  * 14.09.2012 : initially created and documented  
!  * 27.06.2018 : check header
! USAGE
!  Helpful things for example notebooks
!
! GTPack MODULES
!
! --- Rotations ---
!
!  o GTActiveRotation   - shows active rotation definition
!  o GTPassiveRotation  - shows passive rotation definition
!
! DESCRIPTION
!  Additional.m contains all additional functions that are NOT needed directly by the GroupTheory package. 
!  We come across such definitions in the formulation of example notebooks. As an example: A simulation of the active and 
!  passive definition or rotations might be useful for understanding, but is not necessary in GTPack context.
!
! LITERATURE
! 
! TODO
!  -
! RELEASE
!  1.0.0 
! PROBLEMS
!  - 
!--------------------------------------------------------------------------------
!
***)

BeginPackage["GroupTheory`Additional`",{"GroupTheory`Symbols`"}]

 GTActiveRotation  ::usage = "GTActiveRotation[] simulates the active definition of a rotation."
 GTPassiveRotation ::usage = "GTPassiveRotation[] simulates the active definition of a rotation."
Begin["`Private`"] 

(****a* /GTActiveRotation
! NAME
!  GTActiveRotation
! AUTHOR
!  W. Hergert
! PACKAGE
!  Additional.m 
! MODIFICATION HISTORY
!  * 17.03.2015: cosmetics
!  * 27.06.2018 : check header
! USAGE
!  GTActiveRotation simulates the active definition of a rotation
! INPUT
!  use sliders for input
! OUTPUT
!  interactive manipulation
! GTPackOPTIONS
!  -
! STANDRD OPTIONS
!  -
! GTPack MODULES
!  - 
! GTPack NOTEBOOKS 
!  GTRotation.nb
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

GTActiveRotation[]:= Module[{x0,y0},
Manipulate[
	       Show[
	            ListPlot[{{x0,y0}},PlotRange->{{-8,8},{-8,8}},AspectRatio->1,PlotStyle->PointSize[Large]]
               ,Graphics[{{Red,Arrow[{{0,0},{x0,y0}}]},
           	              {Dashed,Line[{{x0,y0},{0,y0}}],Line[{{x0,y0},{x0,0}}]},
           	              {Blue,Arrow[{{0,0},{Cos[\[Theta]]x0-Sin[\[Theta]]y0,Sin[\[Theta]]x0+Cos[\[Theta]]y0}}]},
           	              {Dashed,Line[{{Cos[\[Theta]]x0-Sin[\[Theta]]y0,Sin[\[Theta]]x0+
           	             	             Cos[\[Theta]]y0},{0,Sin[\[Theta]]x0+Cos[\[Theta]]y0}}],
           	            	      Line[{{Cos[\[Theta]]x0-Sin[\[Theta]]y0,Sin[\[Theta]]x0+
           	          	        	     Cos[\[Theta]]y0},{Cos[\[Theta]]x0-Sin[\[Theta]]y0,0}}]},
           	          	        	     {Pink,Circle[{0,0},Sqrt[x0^2+y0^2]]}}
           	            ]
               ]
          ,{{x0,2,"x"},-8,8}
          ,{{y0,2,"y"},-8,8}
          ,{{\[Theta],0,"\[Theta]"},0,2\[Pi]}
]
]

(*
***) 

(****a* /GTPassiveRotation
! NAME
!  GTPassiveRotation
! AUTHOR
!  W. Hergert
! PACKAGE
!  Additional.m 
! MODIFICATION HISTORY
!  * 16.03.2015 : cosmetics
!  * 27.06.2018 : check header
! USAGE
!  GTPassiveRotation simulates the passive definition of a rotation
! INPUT
!  use sliders for input
! OUTPUT
!  interactive manipulation
! GTPackOPTIONS
!  -
! STANDRD OPTIONS
!  -
! GTPack MODULES
!  - 
! GTPack NOTEBOOKS 
!  GTRotation.nb
! DESCRIPTION
!  -
! LITERATURE
!  -
! TODO
!  -
RELEASE
!  1.0.0
! PROBLEMS
!  -
! SOURCE
!--------------------------------------------------------------------------------
!
*)

GTPassiveRotation[] := Module[{x0,y0},
Manipulate[
	       Show[
	       	    ListPlot[{{x0,y0}},PlotRange->{{-8,8},{-8,8}},AspectRatio->1],
                Graphics[{{Red,Arrow[{{0,0},{x0,y0}}]},
                	      {Dashed,Line[{{x0,y0},{0,y0}}],Line[{{x0,y0},{x0,0}}]},
                	      {Blue,Arrow[{{-2Sqrt[x0^2+y0^2]Cos[\[Theta]],-2Sqrt[x0^2+y0^2]Sin[\[Theta]]},
                	      	           {2Sqrt[x0^2+y0^2]Cos[\[Theta]],2Sqrt[x0^2+y0^2]Sin[\[Theta]]}}]},
                	      {Blue,Arrow[{{+2Sqrt[x0^2+y0^2]Sin[\[Theta]],-2Sqrt[x0^2+y0^2]Cos[\[Theta]]},
                	      	           {-2Sqrt[x0^2+y0^2]Sin[\[Theta]],2Sqrt[x0^2+y0^2]Cos[\[Theta]]}}]},
                	      {Dashed,Magenta,Line[{{x0,y0},{(Cos[\[Theta]]x0+Sin[\[Theta]]y0)Cos[\[Theta]],
                	      	                             (Cos[\[Theta]]x0+Sin[\[Theta]]y0)Sin[\[Theta]]}}],
                	                      Line[{{x0,y0},{-(-Sin[\[Theta]]x0+Cos[\[Theta]]y0)Sin[\[Theta]],
                	                      	              (-Sin[\[Theta]]x0+Cos[\[Theta]]y0)Cos[\[Theta]]}}]}
                	       ,{Pink,Circle[{0,0},Sqrt[x0^2+y0^2]]}}
                	     ]
                ]
             ,{{x0,2,"x"},-8,8},
              {{y0,3,"y"},-8,8},
              {{\[Theta],0,"\[Theta]"},0,2\[Pi]}
]
]

(*
***) 

End[] (* End Private Context *)

EndPackage[]
