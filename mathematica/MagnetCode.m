(* ::Package:: *)

(* ::Section:: *)
(*magcode.m*)


(* ::Text:: *)
(*A collection of functions for calculating the force between magnets and coil/magnet systems.*)
(**)
(*Please report problems and suggestions at the GitHub issue tracker:  *)
(*	http://github.com/wspr/magcode/issues*)
(**)
(*Copyright 2010*)
(*Will Robertson*)


(* ::Section:: *)
(*Licence*)


(* ::Text:: *)
(*This package consists of the file MagnetCode.m. It may be freely distributed and modified under the terms & conditions of the Apache License, v2.0:*)
(*	http://www.apache.org/licenses/LICENSE-2.0*)


(* ::Section:: *)
(*Preamble*)


BeginPackage["MagnetCode`"];


MagnetCoilForce::usage = "";


Begin["`Private`"]


(* ::Section:: *)
(*Package*)


Options[MagnetCoilForce]={
  MagnetRadius->0,
  MagnetLength->0,
  CoilRadii->{0,0},
  CoilLength->0,
  Displacement->0,
  Magnetisation->0,
  CoilTurns->0,
  Current->0,
  IntegrationPrecision->2
};


MagnetCoilForce[OptionsPattern[]] := Module[
  {
   expr,
   coilarea,force
  },

   prec=OptionValue[IntegrationPrecision];
   magr=OptionValue[MagnetRadius];
   coilr=OptionValue[CoilRadii][[1]];
   coilR=OptionValue[CoilRadii][[2]];
   magl=OptionValue[MagnetLength];
   coill=OptionValue[CoilLength];
   displ=OptionValue[Displacement];
   current=OptionValue[Current];
   magn=OptionValue[Magnetisation];
   turns=OptionValue[CoilTurns];

  coilarea=coill (coilR-coilr);

  force = 2 current turns magn / coilarea NIntegrate[
      Sum[(-1)^(a+b) MagnetCoilForceKernel[(-1)^a magl, displ+(-1)^b  coill]
          ,{a,0,1},{b,0,1}],
    {r,0,magr},{R,coilr,coilR},
    PrecisionGoal->prec];

  force
]
MagnetCoilForceKernel[l_,L_]:=
  (
    (-l+L) r R EllipticPi[-((4 r R)/(r-R)^2),-((4 r R)/((-l+L)^2+(r-R)^2))]
  ) /
  ( Sqrt[(-l+L)^2+(r-R)^2] (r-R) )


(* ::Section:: *)
(*End*)


End[];
EndPackage[];
