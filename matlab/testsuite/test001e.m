
clear all

magnet_fixed.dim = [0.03 0.04 0.05];
magnet_float.dim = [0.055 0.045 0.035];

magnet_fixed.type = 'cuboid';
magnet_float.type = 'cuboid';

magnet_fixed.magn = 1;
magnet_float.magn = 1;

magnet_fixed.magdir = [30  50];
magnet_float.magdir = [60  45];

displ = [0.1 0.09 0.11];

f_all = magnetforces(magnet_fixed,magnet_float,displ);
k_all = magnetforces(magnet_fixed,magnet_float,displ,'stiffness');

assert( all( round(f_all*1e6) == [-1391969;    -1140254;    -1042102]) , ...
  'Forces components appear incorrect.')

assert( all( round(k_all*1e6) == [ -5424370 ;    2953150 ;    2471220]) , ...
  'Stiffness components appear incorrect.')



