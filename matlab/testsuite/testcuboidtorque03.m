
clear all
disp('=================')
fprintf('TEST cuboid torques singularities: ')

magnet_fixed.type = 'cuboid';
magnet_float.type = 'cuboid';

magnet_fixed.dim = [0.04 0.04 0.04];
magnet_float.dim = magnet_fixed.dim;

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.magdir  = [0 0 1];
magnet_float.magdir  = [0 0 1]; % must be (anti-)parallel

T1 = magnetforces(magnet_fixed,magnet_float,[0.04 0; 0.05 0; 0 0.05],'torque');

assert( all(~isnan(T1(:))) , 'no nans' )

magnet_fixed.magdir  = [0 1 0];
magnet_float.magdir  = [0 0 1];

T2 = magnetforces(magnet_fixed,magnet_float,[0.04; 0.05; 0],'torque');
check2 = [ -1582256; 1474241; -3679 ];

assert( all( round(1e6*T2) == check2 ), 'incorrect reference torques between parallel magnets' )


fprintf('passed\n')
disp('=================')



