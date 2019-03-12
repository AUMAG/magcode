
clear all

magnet_fixed.dim = [0.04 0.04 0.04];
magnet_float.dim =  magnet_fixed.dim;

magnet_fixed.type = 'cuboid';
magnet_float.type = 'cuboid';

% Fixed parameters:
magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;
magnet_fixed.magdir = [0 0 1];
displ = 0.12*[1 1 1];


magnet_float.magdir = [1 1 0];
f1 = magnetforces(magnet_fixed,magnet_float,displ);

% Components:
magnet_float.magdir = [1 0 0];
fc1 = magnetforces(magnet_fixed,magnet_float,displ);

magnet_float.magdir = [0 1 0];
fc2 = magnetforces(magnet_fixed,magnet_float,displ);

f2 = (fc1+fc2)/sqrt(2);

assert ( ...
             isequal ( round( f1 , 4 ) , round( f2 , 4 ) ) , ...
             'Components should sum due to superposition (1)' ...
           )


magnet_float.magdir = [0 1 1];
f1 = magnetforces(magnet_fixed,magnet_float,displ);

% Components:
magnet_float.magdir = [0 1 0];
fc1 = magnetforces(magnet_fixed,magnet_float,displ);

magnet_float.magdir = [0 0 1];
fc2 = magnetforces(magnet_fixed,magnet_float,displ);

f2 = (fc1+fc2)/sqrt(2);

assert ( ...
             isequal ( round( f1 , 4 ) , round( f2 , 4 ) ) , ...
             'Components should sum due to superposition (2)' ...
           )


magnet_float.magdir = [1 1 1];
f1 = magnetforces(magnet_fixed,magnet_float,displ);

% Components:
magnet_float.magdir = [1 0 0];
fc1 = magnetforces(magnet_fixed,magnet_float,displ);

magnet_float.magdir = [0 1 0];
fc2 = magnetforces(magnet_fixed,magnet_float,displ);

magnet_float.magdir = [0 0 1];
fc3 = magnetforces(magnet_fixed,magnet_float,displ);

f2 = (fc1+fc2+fc3)/sqrt(3);

assert ( ...
             isequal ( round( f1 , 4 ) , round( f2 , 4 ) ) , ...
             'Components should sum due to superposition (3)' ...
           )


