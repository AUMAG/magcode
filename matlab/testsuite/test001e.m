
clear all
disp('=================')
fprintf('TEST 001e: ')

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
f_x = magnetforces(magnet_fixed,magnet_float,displ,'x');
f_y = magnetforces(magnet_fixed,magnet_float,displ,'y');
f_z = magnetforces(magnet_fixed,magnet_float,displ,'z');

assert( all(f_all==[f_x(1); f_y(2); f_z(3)]) , ...
  'Forces components calculated separately shouldn''t change.')

k_all = magnetforces(magnet_fixed,magnet_float,displ,'stiffness');
k_x = magnetforces(magnet_fixed,magnet_float,displ,'stiffness','x');
k_y = magnetforces(magnet_fixed,magnet_float,displ,'stiffness','y');
k_z = magnetforces(magnet_fixed,magnet_float,displ,'stiffness','z');

assert( all(k_all==[k_x(1); k_y(2); k_z(3)]) , ...
  'Stiffness components calculated separately shouldn''t change.')

fprintf('passed\n')
disp('=================')



