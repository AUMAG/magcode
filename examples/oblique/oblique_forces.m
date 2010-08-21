function [ displ, forces ] = oblique_forces( varargin )
%OBLIQUE_FORCES
%
% INPUTS:         DEFAULT  DESCRIPTION
%
%  'magn'         (1T)     magnet magnetisation
%  'unitlength'   (0.01)   cube root of the magnet volume, metres
%  'dispratio'    (1)      displacement range in terms of 'unitlength'
%  'points'       (50)     number of displacements to calculate over
%  'standoff'     (1e-6)   minimum distance between "touching" magnets
%  'magangle'   * (45)     inclination angle(s) of the magnets, degrees
%  'magratio'   * (0.5)    ratio(s) between magnet height/weight
%  'gapratio'   * (0)      horizontal offset(s) at zero displacement
%                          in terms of 'unitlength' 
%
% Starred inputs may take vector inputs.
%
% OUTPUTS:
%
%  DISPL    max size:  [  <magangles> <gaps> <mag ratios> <points>]
%  FORCES   max size:  [3 <magangles> <gaps> <mag ratios> <points>]
%

p = inputParser;
p.addParamValue('magn',1);
p.addParamValue('unitlength',0.01);
p.addParamValue('magratio',0.5);
p.addParamValue('dispratio',1);
p.addParamValue('magangle',45);
p.addParamValue('gapratio',0);
p.addParamValue('points',50);
p.addParamValue('standoff',0.0000001);

p.parse(varargin{:});

unitlength = p.Results.unitlength;
magratio = p.Results.magratio;
magn     = p.Results.magn;
dispratio= p.Results.dispratio;
magangle = p.Results.magangle;
gapratio = p.Results.gapratio;
points   = p.Results.points;
standoff = p.Results.standoff;

%%

volume = unitlength^3;

%% Rotation matrices

R21 = @(T) [cosd(T) sind(T) 0;-sind(T) cosd(T) 0; 0 0 1];
R12 = @(T) R21(-T);
R31 = @(T) R21(180-T);
R13 = @(T) R21(T-180);

%% Theory

% Displacements between magnets in the magnets' coordinate system

D1 = @(dx,dy,T,a,d) ...
     [ a; 0; 0 ] + R21(T)*[dx+d; dy; 0];

D2 = @(dx,dy,T,a,d) ...
     [ a; 0; 0 ] + R31(T)*[dx-d; dy; 0];

% Magnet setup

mag = @(a,b,magdir) ...
  struct('dim',[a b b],'magn',magn,'magdir',[magdir 0 0]);

% Forces

f1 = @(dx,dy,T,a,b,d) ...
   magnetforces(mag(a,b,1),mag(a,b,-1),D1(dx,dy,T,a,d));

f2 = @(dx,dy,T,a,b,d) ...
   magnetforces(mag(a,b,1),mag(a,b,-1),D2(dx,dy,T,a,d));

F = @(dx,dy,T,a,b,d) ...
  R12(T)*f1(dx,dy,T,a,b,d) + R13(T)*f2(dx,dy,T,a,b,d);



%% Setup

Ng = length(gapratio);
Nm = length(magratio);
NT = length(magangle);


% initialise variables
forces = repmat(NaN,[3 NT Ng Nm points]);
displ  = repmat(NaN,[  NT Ng Nm points]);

%% Forces calc

for mm = 1:Nm
    
  a = (volume*magratio(mm)^2)^(1/3);
  b = (volume/magratio(mm))^(1/3);
  
  for gg = 1:Ng
    
    d = gapratio(gg)*unitlength;
    
    for tt = 1:NT
    
      displ(tt,gg,mm,:) = linspace(standoff,standoff+dispratio*unitlength,points);
      
      for yy = 1:points
        
        forces(:,tt,gg,mm,yy) = F(0,displ(tt,gg,mm,yy),magangle(tt),a,b,d);
        
      end
    end
  end
end

displ = squeeze(displ);
forces = squeeze(forces);

end

