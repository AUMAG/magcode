
clear all
disp('=================')
fprintf('TEST ''grade'' specification: ')

displ = [0.03 0.05 0.07];
magnet_fixed.type = 'cuboid';
magnet_fixed.dim = [0.01 0.02 0.03];
magnet_fixed.magn = 2*sqrt(42/100); % = 'N42'
magnet_fixed.magdir  = [1 0 0];
magnet_float = magnet_fixed;

magnet_fixed2.type = 'cuboid';
magnet_fixed2.dim = [0.01 0.02 0.03];
magnet_fixed2.grade = 'N42';
magnet_fixed2.magdir  = [1 0 0];
magnet_float2 = magnet_fixed2;

magnet_fixed3.type = 'cuboid';
magnet_fixed3.dim = [0.01 0.02 0.03];
magnet_fixed3.grade = '42';
magnet_fixed3.magdir  = [1 0 0];
magnet_float3 = magnet_fixed3;

magnet_fixed4.type = 'cuboid';
magnet_fixed4.dim = [0.01 0.02 0.03];
magnet_fixed4.grade = '42';
magnet_fixed4.magdir  = [1 0 0];
magnet_float4 = magnet_fixed4;

F1 = magnetforces(magnet_fixed, magnet_float, displ);
F2 = magnetforces(magnet_fixed2,magnet_float2,displ);
F3 = magnetforces(magnet_fixed3,magnet_float3,displ);
F4 = magnetforces(magnet_fixed4,magnet_float4,displ);

assert( all( round(1e6*F1) == round(1e6*F2) ), 'grade spec should be consistent' )
assert( all( round(1e6*F1) == round(1e6*F3) ), 'grade spec should be consistent' )
assert( all( round(1e6*F1) == round(1e6*F4) ), 'grade spec should be consistent' )

fprintf('passed\n')
disp('=================')






