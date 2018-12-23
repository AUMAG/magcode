
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
magnet_float.magdir  = [0 1 0];

T1 = magnetforces(magnet_fixed,magnet_float,[0.02; 0; 0.03],'torque');

magnet_fixed.dim = [0.04 0.02 0.04];
magnet_float.dim = magnet_fixed.dim;

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.magdir  = [0 1 0];
magnet_float.magdir  = [0 0 1];

T2 = magnetforces(magnet_fixed,magnet_float,[0.02; 0.03; 0],'torque');

assert( all( round(1e6*T1) ==  [746760; -1371; -1017634] ), 'incorrect reference torques between parallel magnets' )
assert( all( round(1e6*T2) == [-746760; 1017634; -1371] ), 'incorrect reference torques between parallel magnets' )

fprintf('passed\n')
disp('=================')



