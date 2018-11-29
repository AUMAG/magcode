
clear all
disp('=================')
fprintf('TEST cylinder eccentric force consistency: ')

magnet_fixed.dim = [0.01 0.02];
magnet_float.dim = [0.01 0.02];

magnet_fixed.magn = 1;
magnet_float.magn = 1;

magnet_fixed.type = 'cylinder';
magnet_float.type = 'cylinder';

magnet_fixed.dir  = [0; 0; 1];
magnet_float.dir  = [0; 0; 1]; % must be same

magnet_fixed.magdir  = [0; 0;  1];
magnet_float.magdir  = [0; 0; -1];

zdispl = 0.03;
F1 = magnetforces(magnet_fixed,magnet_float,[0;   0;   zdispl]);
F2 = magnetforces(magnet_fixed,magnet_float,[eps; eps; zdispl]);

check = [round(1000*F1(3)) round(1000*F2(3))];
assert( check(1) == check(2) , 'Inconsistent force with negligible offset (coaxial %i ~= eps offset %i)',check(1),check(2));

fprintf('passed\n')
disp('=================')


