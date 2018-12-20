%% Setup

clear all
close all

%% Parallel magnets example
%
% Reproducing the results of Janssen et al (2010).
% N.B. they don't quite match the graphs in the paper but
% have been verified with Janssen's own code

a = 0.01/2;
b = 0.026/2;
c = 0.014/2;

alpha = 0;
beta  = -0.008;
gamma = 0.015;

lever  = [0; 0; -0.047];

magnet_fixed.type = 'cuboid';
magnet_float.type = 'cuboid';

magnet_fixed.dim = [2*c 2*b 2*a];
magnet_float.dim = [2*a 2*b 2*c];

magnet_fixed.magn = 1.23;
magnet_float.magn = 1.23;

magnet_fixed.magdir = [0 0  1]; %  z
magnet_float.magdir = [0 0 -1]; % -z

N = 51;
offset = repmat([alpha; beta; gamma],[1 N]);
displ = linspace(0, 0.035, N);
displ_range = offset+[1; 0; 0]*displ;

magnet_float.lever = [0; 0; -0.047];

torque_xyz = magnetforces(magnet_fixed,magnet_float,displ_range,'torque');

figure(1); clf; hold on
plot(1000*displ,torque_xyz(1,:))
plot(1000*displ,torque_xyz(2,:))
plot(1000*displ,torque_xyz(3,:))
xlabel('Horizontal displacement, mm')
ylabel('Torques, Nm')
legend('x','y','z')
grid on



%% Orthogonal magnets example
%
% Reproducing the results of Janssen et al (2011).

clear all

a1 = 0.010/2;
b1 = 0.026/2;
c1 = 0.014/2;

a2 = 0.014/2;
b2 = 0.026/2;
c2 = 0.010/2;

alpha = 0;
beta  = -0.008;
gamma = 0.015;

magnet_fixed.type = 'cuboid';
magnet_float.type = 'cuboid';

magnet_fixed.dim = [2*a1 2*b1 2*c1];
magnet_float.dim = [2*a2 2*b2 2*c2];

magnet_fixed.magn = 1.23;
magnet_float.magn = 1.23;

magnet_fixed.magdir = [0 0 1];
magnet_float.magdir = [0 1 0];

N = 51;
offset = repmat([alpha; beta; gamma],[1 N]);
lever  = repmat([0; 0; -0.047],[1 N]);
displ = linspace(0, 0.05, N);
displ_range = offset+[1; 0; 0]*displ;
lever_range = lever-displ_range;

magnet_float.lever = lever_range;

torque_xyz = magnetforces(magnet_fixed,magnet_float,displ_range,'torque');

figure(2); clf; hold on
plot(displ,torque_xyz(1,:))
plot(displ,torque_xyz(2,:))
plot(displ,torque_xyz(3,:))
xlabel('Horizontal $x$ displacement, m')
ylabel('Torques, Nm')

legend('x','y','z')
grid on

