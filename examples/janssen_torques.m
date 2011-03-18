%% Setup
%
% In case you don't have the various bits'n'pieces that I use to create
% my Matlab graphics (probably likely).

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  simple_graph = true;
else
  simple_graph = false;
end

%% Parallel magnets example
%
% Reproducing the results of Janssen et al (2010)

a = 0.01/2;
b = 0.026/2;
c = 0.014/2;

alpha = 0;
beta  = -0.008;
gamma = 0.015;

magnet_fixed.dim = [2*a 2*b 2*c];
magnet_float.dim = [2*a 2*b 2*c];

magnet_fixed.magn = 1.23;
magnet_float.magn = 1.23;

magnet_fixed.magdir = [0 0  1]; %  z
magnet_float.magdir = [0 0 -1]; % -z

N = 51;
offset = repmat([alpha; beta; gamma],[1 N]);
lever  = repmat([0; 0; -0.047],[1 N]);
displ = linspace(0, 0.035, N);
displ_range = offset+[1; 0; 0]*displ;
lever_range = lever-displ_range;

magnet_float.lever = lever_range;

torque_xyz = magnetforces(magnet_fixed,magnet_float,displ_range,'torque');

willfig('janssen'); clf; hold on
plot(displ,torque_xyz(1,:))
plot(displ,torque_xyz(2,:))
plot(displ,torque_xyz(3,:))
xlabel('Horizontal $x$ displacement, m')
ylabel('Torques, Nm')

legend('x','y','z')

if ~simple_graph
  draworigin([0 0],'h','--')
  colourplot
  % matlabfrag('fig/akoun-repro','dpi',3200);
end

