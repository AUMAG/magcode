
clear all
disp('=================')
fprintf('TEST cylinder forces: ')

magnet_fixed.dim = [0.02 0.04];
magnet_float.dim = [0.02 0.04];

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.type = 'cylinder';
magnet_float.type = 'cylinder';

magnet_fixed.dir  = [0 0 1];
magnet_float.dir  = [0 0 1]; % must be same

magnet_fixed.magdir  = [0 0  1];
magnet_float.magdir  = [0 0 -1];

F = magnetforces(magnet_fixed,magnet_float,[0 0 0.05]);

assert( round(1000*F(3)) == 265537 , 'Incorrect reference force between cylindrical magnets' );

fprintf('passed\n')
disp('=================')


