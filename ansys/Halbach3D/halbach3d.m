function [] = halbach3d(n)
%% This script creates Ansys input code, runs Ansys in batch mode and reads in results by using readAnsysHalbach.m.

% clear all;close all;clc;
!del C:\magcode\ansys\Halbach3D\halbach3d.txt
mu = 1.05;
mu0 = 4*pi*10^-7;
h1 = 1.2/(mu*mu0);                  % Coercivity inner magnet
h2 = 1.2/(mu*mu0);                  % Coercivity outer magnet

N = 20;                              % Number of magnets
theta = pi/(N/2);                   % Angle of 'cube' rotation
angle = theta/(2*pi)*360;           % Angle of 'cube' rotation in degrees
mesh = 2;                        % Mesh size
M = N;      % M is used as number of magnets in readAnsys.m

a = .3;                             % Airgap width
b = .3;                             % Airgap height
c = .2;                             % Airgap depth
%d = .1;                           % Outer ring radius
e = .05;                            % Inner ring radius
% f = .025;                           % Outer magnet radial length
f = sqrt(.0054/N);
% g = .025;                           % Outer magnet tangential length
g=f;
% h = .025;                           % Outer magnet depth
h=f;
% k = .02;                           % Inner magnet radial length
k = sqrt(.0024/N);
% l = .02;                           % Inner magnet tangential length 
l = k;
% m = .02;                           % Inner magnet depth
m=l;
% n = .01;                            % Inner magnet offset
d = e+k/2+.023+f/2;

% Run Ansys
runAnsysHalbach3d(a,b,k,l,m,n,N,h1,c,theta,angle,e,f,d,h2,g,h,mesh);

% Run results reader and make scatter and contour plot
%readAnsysHalbach3d(c,d,e,f,g,h,theta,M);
%Arrows(c,e,f,d,g,h,theta,M);
end