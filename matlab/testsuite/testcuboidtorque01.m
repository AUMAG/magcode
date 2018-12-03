
clear all
disp('=================')
fprintf('TEST cuboid torques: ')

magnet_fixed.type = 'cuboid';
magnet_float.type = 'cuboid';

magnet_fixed.dim = [0.04 0.04 0.02];
magnet_float.dim = magnet_fixed.dim;

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.magdir  = [0 0 1];
magnet_float.magdir  = [0 0 1]; % must be (anti-)parallel

T = magnetforces(magnet_fixed,magnet_float,[0.02; 0; 0.03],'torque');

assert( all( round(1e6*T) == [0; 422829; 0] ), 'incorrect reference torques between parallel magnets' )


magnet_fixed.dim = [0.04 0.02 0.04];
magnet_float.dim = magnet_fixed.dim;

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.magdir  = [0 1 0];
magnet_float.magdir  = [0 1 0]; % must be (anti-)parallel

T2 = magnetforces(magnet_fixed,magnet_float,[0; 0.03; 0.02],'torque');

assert( all( round(1e6*T2) == [422829; 0; 0] ), 'incorrect reference torques between parallel magnets' )


magnet_fixed.dim = [0.02 0.04 0.04];
magnet_float.dim = magnet_fixed.dim;

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.magdir  = [1 0 0];
magnet_float.magdir  = [1 0 0]; % must be (anti-)parallel

T3 = magnetforces(magnet_fixed,magnet_float,[0.03; 0.02; 0],'torque');

assert( all( round(1e6*T3) == [0; 0; 422829] ), 'incorrect reference torques between parallel magnets' )

fprintf('passed\n')
disp('=================')



