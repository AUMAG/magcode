(* ::Package:: *)

(* ::Title:: *)
(*MagnetCode.m*)


(* ::Text:: *)
(*A collection of functions for calculating the force between magnets and coil/magnet systems.*)
(**)
(*Please report problems and suggestions at the GitHub issue tracker:  *)
(*	http://github.com/wspr/magcode/issues*)
(**)
(*Copyright 2010-2011*)
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
CoilCoilForce::usage = "";
MagnetCode::unknownMethod="The method `1` is unknown.";


Begin["`Private`"]


(* ::Section:: *)
(*Options and main function*)


Options[MagnetCoilForce]={

  MagnetRadius->0,
  MagnetLength->0,
  Magnetisation->1,

  CoilRadius->0,
  CoilThickness->0,
  CoilLength->0,
  CoilTurns->0,
  CoilTurnsR->0,
  CoilTurnsZ->0,
  Current->0,

  Displacement->0.0,
  Eccentricity->0,

  IntegrationPrecision->2,
  Babic->True,
  Method->Automatic,
  RadialForce->True,
  Verbose->False
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
   turnsR = OptionValue[CoilTurnsR],
   turnsZ = OptionValue[CoilTurnsZ],
   curr  = OptionValue[Current],
   displ = OptionValue[Displacement],
   eccen = OptionValue[Eccentricity],
   method  = OptionValue[Method],
   prec  = OptionValue[IntegrationPrecision],
   radialForceBool = OptionValue[RadialForce]
  },

  If [OptionValue[Verbose], Diagnostic[x_] := Print[x], Diagnostic[x_] := Null];
  
  If[ turns == 0 , turns = turnsR*turnsZ; ];

  If[ turns == 0 , force = 0 ,  

  If[ eccen == 0.0 ,
  If[ displ == 0.0 ,
    force = 0 , 
      If[ OptionValue[CoilThickness] == 0.0,
        Switch[method,
        Automatic,
          Diagnostic["Coaxial thin coil:"];
          force = MagnetThinCoilForce[curr turns,magn,magr,magl,coilr,coill,displ];
        ,
        "Filament",
          Diagnostic["Coaxial thin coil: (filament)"];
          force = ThinCoilForceFilaments[curr,turns,coilr,coill,magn,magr,magl,displ];
        ,
        _, Message[MagnetCode::unknownMethod,method]
        ]
      ,

        Switch[method,
          Automatic,
          Diagnostic["Coaxial thick coil (plain integral):"];
          force = MagnetThickCoilForce[
            curr turns,magn,magr,magl,coilr,coilR,coill,displ,PrecisionGoal->prec];
          ,
          "Babic",
          Diagnostic["Coaxial thick coil (Babic):"];
          force = ThinThickCoilAxialForce[
            curr turns,magn*magl/(4 \[Pi] 10^-7),magr,coilr,coilR,coill,magl,displ,PrecisionGoal->prec];
          ,
          "Filament",
          Diagnostic["Coaxial thick coil (filament):"];
          force = MagnetThickCoilForceFilament[curr,turns,turnsZ,turnsR,magn,magr,magl,coilr,coilR,coill,displ]
          ,
          "Shell",
          Diagnostic["Coaxial thick coil (shell):"];
          force = MagnetThickCoilForceShell[curr,turnsZ,turnsR,magn,magr,magl,coilr,coilR,coill,displ]
          ,
          _, Message[MagnetCode::unknownMethod,method]
          ]
          
        ]
      ]
    ,
      Switch[method,
      Automatic,
        Diagnostic["Eccentric thick coil:"];
        force = MagnetCoilEccentricAxialForce[
          curr * turns, magn, magr, magl, coilr, coilR, coill, displ, eccen, PrecisionGoal->prec];
        ,
      "Filament",
        Diagnostic["Eccentric thick coil: (filament)"];
        force = MagnetThickCoilForceFilament[curr,turns,turnsZ,turnsR,magn,magr,magl,coilr,coilR,coill,displ,eccen,radialForceBool,PrecisionGoal->prec]
        ,
        _, Message[MagnetCode::unknownMethod,method]
      ]
  ];
  ];

  force
]


(* ::Section:: *)
(*Filament models*)


CoilCoilForce[I1_,I2_,r_,R_,z_] = With[
  { m = 4 r R/((r+R)^2+z^2) },
  I1 I2 4 \[Pi] 10^-7 Sqrt[m] z
  ( 
    2 EllipticK[m]-
    (2-m)/(1-m) EllipticE[m]
  ) / (4 Sqrt[r R])
];

CoilCoilForce[I1_,I2_,r_,R_,z_,e_,rtf_,param__] := Module[
{xs,ys,V,VV,m,m1,m1sq,KK,EE,fz,fe},

  If[ e == 0.0 ,
    (* use zero eccentricity solution (zero radial term): *)
    If[rtf,
      force = {CoilCoilForce[I1,I2,r,R,z],0.0};
    ,
      force = CoilCoilForce[I1,I2,r,R,z];
    ]
  ,

  (* else, with non-zero eccentricity: *)
  xs = e - R Sin[t];
  ys =   + R Cos[t];

  VV = xs^2+ys^2;
  V = Sqrt[VV];

  m1 = 4 / ( (r+V)^2 + z^2 );
  m1sq = Sqrt[m1];
  m = r V m1;
  KK = EllipticK[m];
  EE = EllipticE[m];

  fz = (* axial term: *)
   z NIntegrate[
    m1sq / VV ( ys Cos[t] - xs Sin[t] )( 2 KK - EE (2-m)/(1-m) )
    ,{t,0,2\[Pi]},param];

  If[rtf,
    fe = (* radial term: *)
     NIntegrate[
      m1sq Sin[t] ( -2 KK + EE (2 - m( r/V + 1 ) )/(1-m) )
     ,{t,0,2\[Pi]},param];

    force = 0.5 10^-7 R I1 I2 {fz,fe};
  ,
    force = 0.5 10^-7 R I1 I2 fz;
  ]
  ];

  force
]


ThinCoilForceFilaments[Current_,CoilTurns_,CoilRadius_,CoilLength_,Magnetisation_,MagnetRadius_,MagnetLength_,Displacement_] := Module[
{f1,f2},

  f1=Table[-CoilLength/2+(n-1) CoilLength/(CoilTurns-1),{n,1,CoilTurns}];
  f2=Table[-MagnetLength/2+(n-1) MagnetLength/(CoilTurns-1),{n,1,CoilTurns}];

  Sum[
    CoilCoilForce[Current, Magnetisation MagnetLength/(4\[Pi] 10^(-7) CoilTurns),CoilRadius,MagnetRadius,Displacement+z1+z2]
    ,{z1,f1},{z2,f2}
  ]
]


(* ::Section:: *)
(*Thin coil models*)


MagnetThinCoilForce[currturns_,magn_,magr_,magl_,coilr_,coill_,displ_] :=
  magn currturns / ( 2 coill ) MagnetThinCoilForceKernel[
    coilr,magr,-coill/2,coill/2,displ-magl/2,displ+magl/2
   ]

MagnetThinCoilForceKernel[r1_,r2_,z1_,z2_,z3_,z4_] :=
  fff2[r1,r2,z2,z3,z4]-fff2[r1,r2,z1,z3,z4]

fff2[r1_,r2_,zt_,z3_,z4_] :=
  fff3[r1,r2,zt-z4]-fff3[r1,r2,zt-z3]

fff3[r1_,r2_,z_] :=
  If[z==0,0,
    fff4[z,(r1-r2)^2/z^2+1,
      Sqrt[(r1+r2)^2+z^2],
      (4 r1 r2)/((r1+r2)^2+z^2)
    ]
  ]

fff4[m1_,m2_,m3_,m4_]:=
  m1 m3 ( m2 EllipticK[m4] - EllipticE[m4] +
    If[m2==1,0,
      m2 ( (m1/m3)^2 - 1 ) EllipticPi[m4/(1-m2),m4]
    ]
  )


(* ::Section:: *)
(*Thick coil models*)


(* ::Subsection:: *)
(*Direct integral solution*)


MagnetThickCoilForce[currturns_,magn_,magr_,magl_,coilr_,coilR_,coill_,displ_,param_] :=
  currturns magn / ( coill (coilR-coilr) ) NIntegrate[
    MagnetThickCoilForceKernel[ magr, magl, R, Z, displ ],
    {R,coilr,coilR}, {Z,-coill/2,coill/2}, param
   ]

MagnetThickCoilForceKernel[r_,maglength_,R_,Z_,d_]:=
   MagnetThickCoilForceKernel2[r,d+maglength/2,R,Z]-
    MagnetThickCoilForceKernel2[r,d-maglength/2,R,Z]

MagnetThickCoilForceKernel2[r_,z_,R_,Z_]:=
  With[ { m1 = 4 r R / ((r+R)^2+(z-Z)^2) },
    Sqrt[(r+R)^2+(z-Z)^2]
    ( (1-m1/2)EllipticK[m1]-EllipticE[m1] )
 ]


(* ::Subsection:: *)
(*Babic's solution*)


(* curr turns,magn*magl/(4 \[Pi] 10^-7),magr,coilr,coilR,coill,magl,displ,... *)
ThinThickCoilAxialForce[NI1_,NI2_,R_,R1_,R2_,Z1_,Z2_,D_,opt_] :=
  (2 \[Pi] 10^-7 NI1 NI2 R^3)/(3 Z1 Z2(R2-R1))*
  Sum[ ii jj kk ThinThickCoilAxialForceKernel[ii,jj,kk,R1,R2,R,Z1,Z2,D,opt] ,
      {ii,{1,-1}}, {jj,{1,-1}}, {kk,{1,-1}}]

ThinThickCoilAxialForceKernel[ii_,jj_,kk_,R1_,R2_,R_,Z1_,Z2_,D_,opt_] :=
  Block[{\[Rho],t,m,m2,result},
  \[Rho] = 1/(2R) (R1+R2+kk (R2-R1));
  t = (D+ii Z1/2+jj Z2/2)/R;

  If[t==0,result=0,
    m = (4 \[Rho])/((\[Rho]+1)^2+t^2);
    m2 = Sqrt[t^2+1];
    result = Sqrt[m \[Rho]] (-EllipticK[m] ( (m2+2)/(m2+1) (t^2-2)+(\[Rho]^2+\[Rho]+2)-2/(\[Rho]+1))+EllipticE[m] (4\[Rho])/m)+
      (\[Pi]/2)/Abs[t] (\[Rho] (\[Rho]^2-3)Sign[\[Rho]-1]( 1 - HeumanLambda[Abs[ArcSin[(\[Rho]-1)/(\[Rho]+1) Sqrt[1/(1-m)]]],m] ) +
        (t^2-2)Sqrt[t^2+1] ( 1-HeumanLambda[Abs[ArcSin[t/(1+m2)]],m] +
          Sign[\[Rho]-Sqrt[t^2+1]](1-HeumanLambda[Abs[ArcSin[t/(1+m2) Sqrt[1/(1-m)]]],m]) ) )+
      -6 NIntegrate[ArcSinh[(\[Rho]+Cos[2\[CurlyPhi]])/Sqrt[Sin[2\[CurlyPhi]]^2+t^2]],{\[CurlyPhi],0,\[Pi]/2},opt];
  ];
  - t * result
]

HeumanLambda[\[Phi]_,m_] := EllipticF[\[Phi],1-m]/EllipticK[1-m] + 2/\[Pi] EllipticK[m]JacobiZeta[\[Phi],1-m]


(* ::Subsection:: *)
(*Filament method*)


MagnetThickCoilForceFilament[curr_,turns_,turnsZ_,turnsR_,magn_,magr_,magl_,coilr_,coilR_,coill_,displ_] :=
  MagnetThickCoilForceFilament[curr,turns,turnsZ,turnsR,magn,magr,magl,coilr,coilR,coill,displ,0,False,Null]

MagnetThickCoilForceFilament[curr_,turns_,turnsZ_,turnsR_,magn_,magr_,magl_,coilr_,coilR_,coill_,displ_,eccen_,rtf_,param__] :=
  Module[{fr,f1,f2,n},

  fr=coilr+Table[(n-1) (coilR-coilr)/(turnsR-1),{n,1,turnsR}];
  f1=Table[-coill/2+(n-1) coill/(turnsZ-1),{n,1,turnsZ}];
  f2=Table[-magl/2+(n-1) magl/(turns-1),{n,1,turns }];

  Sum[
    CoilCoilForce[curr, magn*magl/(4\[Pi] 10^(-7) turns) ,rr,magr,displ+z1+z2,eccen,rtf,param]
    ,{rr,fr},{z1,f1},{z2,f2}
  ]
  ]


(* ::Subsection:: *)
(*Shell method*)


MagnetThickCoilForceShell[curr_,turnsZ_,turnsR_,magn_,magr_,magl_,coilr_,coilR_,coill_,displ_] :=
Module[{fr},
  fr=coilr + Table[(n-1) (coilR-coilr)/(turnsR-1),{n,1,turnsR}];
  1/turnsR Sum[
    MagnetThinCoilForce[curr turnsR turnsZ,magn,magr,magl,rr,coill,displ]
    ,{rr,fr}
  ]
]
MagnetThickCoilForceShell[curr_,turnsZ_,1,magn_,magr_,magl_,coilr_,coilR_,coill_,displ_] :=
  MagnetThinCoilForce[curr turnsZ,magn,magr,magl,(coilR+coilr)/2,coill,displ]


(* ::Subsection:: *)
(*Eccentric integral solution*)


MagnetCoilEccentricAxialForce[currturns_, magn_, magr_, magl_, coilr_, coilR_, coill_, displ_, eccen_, param_] :=
  force = - currturns magn / ( \[Pi] coill (coilR-coilr) ) *
    NIntegrate[
      MagnetCoilForceEccentricKernel[displ+magl/2, eccen, magr] -
      MagnetCoilForceEccentricKernel[displ-magl/2, eccen, magr]
    ,
      {\[Phi]2,0,2\[Pi]},{R,coilr,coilR},{Z,-coill/2,coill/2},
      param
    ]

MagnetCoilForceEccentricKernel[z_, e_, r2_] :=
  Module[{x1,y1,r,\[Phi],m1},
    x1 = e+r2 Cos[\[Phi]2];
    y1 = r2 Sin[\[Phi]2];
    r = Sqrt[x1^2 + y1^2];
    \[Phi] = ArcTan[y1, x1];
    m1 = 4 r R / ((r+R)^2+(z-Z)^2);
    Sqrt[(r R)/m1] (EllipticE[m1]-(1-m1/2) EllipticK[m1])
  ]


(* ::Section:: *)
(*End*)


End[];
EndPackage[];
