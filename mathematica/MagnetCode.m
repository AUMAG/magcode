(* ::Package:: *)

(* ::Section:: *)
(*MagnetCode.m*)


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


Begin["`Private`"]


(* ::Section:: *)
(*Package*)


Options[MagnetCoilForce]={

  MagnetRadius->0,
  MagnetLength->0,
  Magnetisation->1,

  CoilRadius->0,
  CoilThickness->0,
  CoilLength->0,
  CoilTurns->0,
  Current->0,

  Displacement->0.0,
  Eccentricity->0,

  IntegrationPrecision->2
};


MagnetCoilForce[OptionsPattern[]] := Module[
  {
   coilarea,
   force,
   magr  = OptionValue[MagnetRadius],
   magl  = OptionValue[MagnetLength],
   magn  = OptionValue[Magnetisation],
   coilr = OptionValue[CoilRadius],
   coilR = OptionValue[CoilRadius]+OptionValue[CoilThickness],
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
      If[ OptionValue[CoilThickness] == 0.0 ,
        force =  magn turns curr / ( 2 coill ) MagnetThinCoilForceKernel[
          coilr,magr,-coill/2,coill/2,displ-magl/2,displ+magl/2]
      ,
        force = 2 curr turns magn / coilarea NIntegrate[
          Sum[ e1 e2 MagnetCoilForceKernel[ e1 magl/2, displ + e2 coill/2 ]
              ,{e1,{1,-1}},{e2,{1,-1}}],
          {r,0,magr},{R,coilr,coilR},
          PrecisionGoal->prec
        ]
      ]
    ,
      force = curr turns magn / ( 4 \[Pi] coilarea ) NIntegrate[
        Sum[
          e1 e2 MagnetCoilForceEccentricKernel[
            e1 magl/2, displ+e2 coill/2, eccen
          ],
          {e1,{1,-1}},{e2,{1,-1}}
        ]
        ,
        {r2,0,magr}, {\[Phi]2,0,2\[Pi]}, {R,coilr,coilR},
        PrecisionGoal->prec
       ]
    ]
  ];

  force
]


MagnetThinCoilForceKernel[r1_,r2_,z1_,z2_,z3_,z4_]:=
  fff2[r1,r2,z2,z3,z4]-fff2[r1,r2,z1,z3,z4]
fff2[r1_,r2_,zt_,z3_,z4_]:=
  fff3[r1,r2,zt,z4]-fff3[r1,r2,zt,z3]
fff3[r1_,r2_,zt_,ztt_]:=
  fff4[
    zt-ztt,r1-r2,r1+r2,
    (r1-r2)^2+(zt-ztt)^2,Sqrt[(r1+r2)^2+(zt-ztt)^2],
    (4 r1 r2)/((r1+r2)^2+(zt-ztt)^2)
  ]
fff4[c1_,c2_,c3_,c4_,c5_,c6_]:=
  If[c1==0,0,
    (c4 c5 )/c1 EllipticK[c6]-
    c1 c5 EllipticE[c6]-
    If[c2==0,0,(c3^2 c4)/(c1 c5) EllipticPi[-((c1^2 c6)/c2^2),c6]]
  ]


MagnetCoilForceKernel[l_,L_]:=
  EllipticPi[
    -4 r R / (r-R)^2,
    -4 r R / ( (L-l)^2 + (r-R)^2 )
  ] * (L-l) r R /
  ( (r-R) Sqrt[(L-l)^2+(r-R)^2] )


MagnetCoilForceEccentricKernel[Z_,L_,e_] =
  2 r R (L-Z) / ( (r-R) Sqrt[ (r-R)^2+(L-Z)^2 ] ) *
  ( EllipticPi[ (-4 r R)/(r-R)^2,\[Phi]/2  , -4 r R /( (r-R)^2+(L-Z)^2 ) ] -
    EllipticPi[ (-4 r R)/(r-R)^2,\[Phi]/2-\[Pi], -4 r R /( (r-R)^2+(L-Z)^2 ) ] ) //.
    {r-> Sqrt[x1^2+y1^2], \[Phi]->ArcTan[y1,x1]} //.
    {x1->x2+e, y1->y2} //.
    {x2->r2 Cos[\[Phi]2], y2->r2 Sin[\[Phi]2]};


(* ::Section:: *)
(*End*)


End[];
EndPackage[];
