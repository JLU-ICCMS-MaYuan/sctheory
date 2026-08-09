(****m* /LandauTheory.m
!
! NAME
!  LandauTheory.m
! AUTHOR
!  M. Geilhufe
! MODIFICATION HISTORY
!  * Oct/24/2022 - initially created and documented  
! USAGE
! Contains all molecules related to Landau theory, i.e., the phenomenological theory of order based on free energy expansion
!
!  GTPack MODULES
! 
!   * LandauExpansion         - gives a symmetry-adapted free energy expansion of generalized order parameters.
!
! DESCRIPTION
!   -
! LITERATURE
!   -
! TODO
!  -
! PROBLEMS
! -
!--------------------------------------------------------------------------------
!
***)

BeginPackage["GroupTheory`LandauTheory`",{"GroupTheory`RepresentationTheory`","GroupTheory`Basic`","GroupTheory`Symbols`","GroupTheory`Auxiliary`"}]
 

(*---------- Basic Landau expansion ----------*)
 GTLandauExpansion   ::usage=   "GTLandauExpansion[\*StyleBox[\"irreducible representations, basis, order\",\"TI\"]] gives a symmetry-adapted free energy expansion of generalized order parameters."
 
  
(*---------- Options ---------------------------------------------*) 	
Options[GTLandauExpansion] = {GOSuperconductor -> {}}
 
Begin["`Private`"]


GTLandauExpansion[irepso_, basiso_, order_, OptionsPattern[]] := 
 Module[{truepos, ireps, basis, rep, repo, baso, gl, coeff, phase, 
   rule, evals, evecs},
  (*group order*)
  gl = Length[irepso[[1]]];
  (*check which terms correspond to superconducting order*)
  truepos = 
   Flatten[Position[
     OptionValue[GOSuperconductor][[1 ;; Min[Length[irepso],Length[OptionValue[GOSuperconductor]]]]], True]];
  (*form U(1) invariant product for superconductors*)
  ireps = irepso;
  basis = basiso;
  Do[
   AppendTo[basis, Conjugate[basis[[tr]]]];
   AppendTo[ireps, Conjugate[ireps[[tr]]]];
   , {tr, truepos}];
  rule = Flatten[
    Table[
     phase = If[MemberQ[truepos,b], I RandomReal[{1, 3}], 0];
      Table[bi -> RandomReal[{1, 3}] + phase, {bi, basis[[b]]}], {b, 1, Length[basiso]}]];
  rep = ireps[[1]];
  Do[rep = GTClebschGordanSum[rep, ir, GOFast -> True], {ir, 
    ireps[[2 ;;]]}];
  repo = rep;
  basis = Flatten[basis];
  baso = basis;
  Sum[
   repo = GTDirectProductRep[repo, rep, GOFast -> True];
   baso = FullSimplify[Flatten[KroneckerProduct[baso, basis]]];

   (*{evals, evecs} = Eigensystem[Total[repo]/gl];
   coeff = Pick[evecs, (# == 1) & /@ evals];*)
   coeff = Select[Orthogonalize[Total[repo]/gl], Total[Abs[#]] > 0 &];
   If[Length[coeff] > 0,
    coeff = Select[Union[FullSimplify[coeff . baso]], Abs[Im[# /. rule]] < 0.001 &];
    Sum[C[o, i] coeff[[i]], {i, 1, Length[coeff]}],
    0]
   , {o, 2, order}]
  ]
  
  
(*
***)


End[]
	 
EndPackage[]
  