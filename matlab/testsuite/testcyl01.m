
clear all

magnet_fixed.type = 'cylinder';
magnet_float.type = 'cylinder';

magnet_fixed.radius = 0.02;
magnet_float.radius = 0.02;

magnet_fixed.height = 0.04;
magnet_float.height = 0.04;

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.dir  = [0; 0; 1];
magnet_float.dir  = [0; 0; 1]; % must be same

magnet_fixed.magdir  = [0; 0;  1];
magnet_float.magdir  = [0; 0; -1];


%% Equal radii & height

F = magnetforces(magnet_fixed,magnet_float,[0; 0; 0.05]);

check = [round(1000*F(3)) 265537];
assert( check(1) == check(2) , 'Incorrect reference force between cylindrical magnets (calc %i ~= ref %i)',check(1),check(2));


%% Different radii & height

magnet_float.radius = 0.015;
magnet_float.height = 0.03;

F = magnetforces(magnet_fixed,magnet_float,[0; 0; 0.05]);

check = [round(1000*F(3)) 111490];
assert( check(1) == check(2) , 'Incorrect reference force between cylindrical magnets (calc %i ~= ref %i)',check(1),check(2));


%% Different height

magnet_fixed.radius = 0.02;
magnet_float.radius = 0.02;

magnet_fixed.height = 0.02;
magnet_float.height = 0.01;

F = magnetforces(magnet_fixed,magnet_float,[0; 0; 0.05]);

check = [round(1000*F(3)) 18953];
assert( check(1) == check(2) , 'Incorrect reference force between cylindrical magnets (calc %i ~= ref %i)',check(1),check(2));



