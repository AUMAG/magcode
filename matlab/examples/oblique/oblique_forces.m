function [ displ, forces ] = oblique_forces( varargin )
%OBLIQUE_FORCES
%
% This function calculates the forces generated on a magnetic spring
% consisting of two pairs of inclined permanent magnets.
%
% INPUTS:         DEFAULT  DESCRIPTION
%
%  'magn'         1T      magnet magnetisation
%  'unitlength' * 0.01    cube root of the magnet volume, metres
%  'dispratio'    1       displacement range in terms of 'unitlength'
%  'points'       50      number of displacements to calculate over
%  'dispoffset'   [0;0;0] offset added to each displacement
%  'magangle'   * 45      inclination angle(s) of the magnets, degrees
%  'magratio'   * 0.4     ratio(s) between magnet height/weight
%  'gapratio'   * 0       horizontal offset(s) at zero displacement
%                           in terms of 'unitlength' 
%
% Starred inputs may take vector inputs.
%
% OUTPUTS:
%
%  DISPL   max size: [  <lengths> <magangles> <gaps> <mag ratios> <points>]
%  FORCES  max size: [3 <lengths> <magangles> <gaps> <mag ratios> <points>]
%
% Outputs are SQUEEZEd, so non-vector input will reduce the number
% of dimensions in the output arrays.
%
% Coordinate system:
%    X - "horizontal" in the plane of the magnet rotations
%    Y - "vertical"
%    Z - "out of plane"
%

%% Inputs and initial setup

p = inputParser;
p.addParamValue('magn',1);
p.addParamValue('unitlength',0.01);
p.addParamValue('magratio',0.4);
p.addParamValue('dispratio',1);
p.addParamValue('magangle',45);
p.addParamValue('gapratio',0);
p.addParamValue('points',50);
p.addParamValue('dispoffset',[0; 0; 0]);

p.parse(varargin{:});

unitlength = p.Results.unitlength;
magratio   = p.Results.magratio;
magn       = p.Results.magn;
dispratio  = p.Results.dispratio;
magangle   = p.Results.magangle;
gapratio   = p.Results.gapratio;
points     = p.Results.points;
dispoffset = p.Results.dispoffset;

%% Theory

% Rotation matrices
%
% R21 rotates a vector in frame 1 into frame 2.
% (Consider it a CCW rotation.) E.g., if frame 2 is +45° from frame 1,
% a vector at 60° in frame 1 is at 15° in frame 2.
%
% R12 is obviously the opposite, but it can also be considered to be simply
% a regular rotation matrix; a vector at 15° applied to R12(45°) will then
% be at 60°.

R21 = @(T) [cosd(T) sind(T) 0;-sind(T) cosd(T) 0; 0 0 1];
R12 = @(T) R21(-T);
R31 = @(T) R21(180-T);
R13 = @(T) R21(T-180);

% Magnet setup

mag = @(a,b,magdir) ...
  struct('dim',[a b b],'magn',magn,'magdir',[magdir 0 0]);

% Displacements between magnets in the magnets' coordinate system

D1 = @(dy,T,a,d)  [ a; 0; 0 ] + R21(T)*([ d; dy; 0]+dispoffset);
D2 = @(dy,T,a,d)  [ a; 0; 0 ] + R31(T)*([-d; dy; 0]+dispoffset);

% Forces

f1 = @(dy,T,a,b,d)  magnetforces(mag(a,b,1),mag(a,b,-1),D1(dy,T,a,d));
f2 = @(dy,T,a,b,d)  magnetforces(mag(a,b,1),mag(a,b,-1),D2(dy,T,a,d));
F  = @(dy,T,a,b,d)  R12(T)*f1(dy,T,a,b,d) + R13(T)*f2(dy,T,a,b,d);

%% Setup

Ns = length(unitlength);
Ng = length(gapratio);
Nm = length(magratio);
NT = length(magangle);

% initialise variables
forces = repmat(NaN,[3 Ns NT Ng Nm points]);
displ  = repmat(NaN,[  Ns NT Ng Nm points]);

%% Forces calc

for uu = 1:Ns
  volume = unitlength(uu)^3;
  for mm = 1:Nm
    a = (volume*magratio(mm)^2)^(1/3);
    b = (volume/magratio(mm))^(1/3);
    for gg = 1:Ng
      d = gapratio(gg)*unitlength(uu);
      for tt = 1:NT
        displ(uu,tt,gg,mm,:) = ...
          linspace(eps,dispratio*unitlength(uu),points);
        for yy = 1:points
          forces(:,uu,tt,gg,mm,yy) = ...
            F(displ(uu,tt,gg,mm,yy),magangle(tt),a,b,d);
        end
      end
    end
  end
end

displ = squeeze(displ);
forces = squeeze(forces);

end

