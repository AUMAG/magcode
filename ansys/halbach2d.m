%% This script creates Ansys input code, runs Ansys in batch mode and reads in results by using readAnsysHalbach.m.
clear all;close all;clc;
!del C:\magcode\ansys\halbach2d.txt
mu = 1.05;
mu0 = 4*pi*10^-7;
h1 = 1.2/(mu*mu0);                  % Coercivity inner magnet
h2 = 1.2/(mu*mu0);                  % Coercivity outer magnet
pi = 4*atan(1);
N = 16;                              % Number of magnets
theta = pi/(N/2);                   % Angle of 'cube' rotation
angle = theta/(2*pi)*360;           % Angle of 'cube' rotation in degrees
mesh = .005;                        % Mesh size

a = .4;                             % Airgap width
b = .4;                             % Airgap heigth
c = .15;                            % Outer ring radius
d = .06;                            % Inner ring radius
e = .040;                           % Outer magnet radial length
f = .040;                           % Outer magnet tangential length
g = .015;                           % Inner magnet radial length
h = .015;                           % Inner magnet tangential length 

runAnsysHalbach(a,b,N,h1,c,theta,angle,e,f,d,h2,g,h,mesh);

M = N;      % M is used as number of magnets in readAnsys.m

% Run results reader and make scatter and contour plot
readAnsysHalbach(c,d,e,f,g,h,theta,M);