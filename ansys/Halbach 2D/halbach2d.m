function [] = halbach2d(N)
%% This script creates Ansys input code, runs Ansys in batch mode and reads in results by using readAnsysHalbach.m.

%clear all;close all;clc;
!del C:\magcode\ansys\halbach2d.txt
mu = 1.05;
mu0 = 4*pi*10^-7;
h1 = 1.2/(mu*mu0);                  % Coercivity inner magnet
h2 = 1.2/(mu*mu0);                  % Coercivity outer magnet
%pi = 4*atan(1);
%N = 8;                              % Number of magnets
theta = pi/(N/2);                   % Angle of 'cube' rotation
angle = theta/(2*pi)*360;           % Angle of 'cube' rotation in degrees
mesh = .003;                        % Mesh size
M = N;      % M is used as number of magnets in readAnsys.m

a = .3;                             % Airgap width
b = .3;                             % Airgap height
%c = .075;                            % Outer ring radius
d = .05;                            % Inner ring radius
%e = .015;                           % Outer magnet radial length
%f = .015;                           % Outer magnet tangential length
%g = .01;                           % Inner magnet radial length
%h = .01;                           % Inner magnet tangential length 
e=sqrt(.0054/N);
f=e;
g=sqrt(.0024/N);
h=g;
c=d+g/2+.023+e/2;   % 0.023 is fixed distance between magnets of inner and outer ring


% Run Ansys
runAnsysHalbach(a,b,N,h1,c,theta,angle,e,f,d,h2,g,h,mesh);

% Run results reader and make scatter and contour plot
readAnsysHalbach(c,d,e,f,g,h,theta,M);
%Arrows(c,e,f,d,g,h,theta,M);
end