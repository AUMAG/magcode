
clear all

magnet_fixed.type = 'cuboid';
magnet_float.type = 'cuboid';

magnet_fixed.dim = [0.03 0.04 0.05];
magnet_float.dim = [0.02 0.03 0.04];

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.magdir  = [0 0 1];
magnet_float.magdir  = [0 0 1]; % must be (anti-)parallel

T1 = magnetforces(magnet_fixed,magnet_float,[0.03; 0.04; 0.05],'torque');

assert( all(isreal(T1)) ,'T1 not real')

check1 = [-548582; 459351; -34579];

assert( all( round(1e6*T1) == check1 ), 'incorrect reference torques between parallel magnets' )




