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
CalculateCoilParams::usage = "";

MagnetCoilForceEccentricKernel::usage="";


Begin["`Private`"]


(* ::Section:: *)
(*Package*)


Options[MagnetCoilForce]={

  MagnetRadius->0,
  MagnetLength->0,
  Magnetisation->1,

  CoilRadii->{0,0},
  CoilLength->0,
  CoilTurns->0,
  Current->0,

  Displacement->0.0,
  Eccentricity->0.0,

  IntegrationPrecision->2
};


MagnetCoilForce[OptionsPattern[]] := Module[
  {
   coilarea,
   force,
   magr  = OptionValue[MagnetRadius],
   magl  = OptionValue[MagnetLength],
   magn  = OptionValue[Magnetisation],
   coilr = OptionValue[CoilRadii][[1]],
   coilR = OptionValue[CoilRadii][[2]],
   coill = OptionValue[CoilLength],
   turns = OptionValue[CoilTurns],
   curr  = OptionValue[Current],
   displ = OptionValue[Displacement],
   eccen = OptionValue[Eccentricity],
   prec  = OptionValue[IntegrationPrecision]
  },
  
  coilarea=coill (coilR-coilr);
  
  If[ displ==0.0 ,
    force = 0 , 
    If[ eccen==0.0 ,
      force = 2 curr turns magn / coilarea NIntegrate[
        Sum[e1 e2 MagnetCoilForceKernel[ e1 magl, displ + e2 coill ]
            ,{e1,{1,-1}},{e2,{1,-1}}],
        {r,0,magr},{R,coilr,coilR},
        PrecisionGoal->prec
      ]
    ,
      force = curr turns magn / ( 4 \[Pi] coilarea ) NIntegrate[
        Sum[e1 e2 MagnetCoilForceEccentricKernel[e2 magl/2,displ+e1 coill/2,eccen],
          {e1,{1,-1}},{e2,{1,-1}}]
        ,
        {r2,0,magr}, {\[Phi]2,0,2\[Pi]}, {R,coilr,coilR},
        PrecisionGoal->prec
       ]
    ]
  ];

  force
]
MagnetCoilForceKernel[l_,L_]:=
  ( (-l+L) r R EllipticPi[-((4 r R)/(r-R)^2),-((4 r R)/((-l+L)^2+(r-R)^2))] ) /
  ( Sqrt[(-l+L)^2+(r-R)^2] (r-R) )


MagnetCoilForceEccentricKernel[l_,L_,e_]=
  ( 2 R r2 (l-L) *
    Sqrt[
         1-(2 r R ( Cos[\[Phi]]-1 ) ) /
           ( (l-L)^2+(R-r)^2 )
        ] 
  ) /
  (
    (R-r) Sqrt[(l-L)^2+r^2+R^2-2 r R Cos[\[Phi]]]
  ) *
  ( 
     EllipticPi[-(4 r R)/(R-r)^2,\[Phi]/2  ,-(4 r R)/( (l-L)^2+(R-r)^2 )]
   - EllipticPi[-(4 r R)/(R-r)^2,\[Phi]/2-\[Pi],-(4 r R)/( (l-L)^2+(R-r)^2 )]
  )//.
          { r -> Sqrt[x^2+y^2] , \[Phi] -> ArcTan[y,x] } //.
          { x -> x2 + e , y -> y2 } //.
          { x2 -> r2 Cos[\[Phi]2] , y2 -> r2 Sin[\[Phi]2] }//.{(e+r2 Cos[\[Phi]2])^2+r2^2 Sin[\[Phi]2]^2->FullSimplify[ExpandAll[(e+r2 Cos[\[Phi]2])^2+r2^2 Sin[\[Phi]2]^2]]};


CalculateCoilParams[parameters__] := 
  ReleaseHold[{
               Hold@Symbol["CoilRadii"]->{RadiusInner,RadiusOuter},
               Hold@Symbol["CoilTurns"]->TotalTurns,
               Hold@Symbol["CoilLength"]->CoilLength,
               Hold@Symbol["Current"]->Current
              } //. 
  Flatten@{{
    RadiusOuter -> RadiusInner + ( WireLength (2WireRadius+2WireCoating)^2 ) /
                                ( 2 \[Pi] CoilLength (RadiusInner+WireRadius+WireCoating) ) ,
    WireLength -> CoilResistance WireArea / WireResistivity ,
    WireResistivity -> 1.7*10^-8, (* copper *)
    WireArea -> \[Pi] (WireRadius+WireCoating)^2,
    TotalTurns -> Round[TurnsZ TurnsR],
    TurnsZ -> CoilLength / (2WireRadius+2WireCoating),
    TurnsR -> (RadiusOuter-RadiusInner) / (2WireRadius+2WireCoating)
    },parameters}//.{cm->0.01,mm->0.001,volts->1,ohms->1}]


(* ::Section:: *)
(*End*)


End[];
EndPackage[];
