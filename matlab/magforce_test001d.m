  
disp('=================') 
fprintf('TEST 001d: ') 
 
magnet_fixed.dim  =  [0.04 0.04 0.04]; 
magnet_float.dim  =   magnet_fixed.dim; 
 
% Fixed parameters: 
magnet_fixed.magn  =  1.3; 
magnet_float.magn  =  1.3; 
magnet_fixed.magdir  =  [0  90];  % $z$ 
displ  =  0.12 * [1 1 1]; 
 
  
magnet_float.magdir  =  [45  0];  % $\vec e_x+\vec e_y$ 
f1  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
% Components: 
magnet_float.magdir  =  [0  0];  % $\vec e_x$ 
fc1  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
magnet_float.magdir  =  [90  0];  % $\vec e_y$ 
fc2  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
f2  =  (fc1+fc2)/sqrt(2); 
 
 
 assert (  ... 
      isequal ( chop( f1 , 4 ) , chop ( f2 , 4 ) ) ,  ... 
      'Components should sum due to superposition'  ... 
    ) 
 
 
 
 
 
  
magnet_float.magdir  =  [0  45];  % $\vec e_y+\vec e_z$ 
f1  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
% Components: 
magnet_float.magdir  =  [0  0];  % $\vec e_x$ 
fc1  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
magnet_float.magdir  =  [0  90];  % $\vec e_z$ 
fc2  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
f2  =  (fc1+fc2)/sqrt(2); 
 
 
 assert (  ... 
      isequal ( chop( f1 , 4 ) , chop ( f2 , 4 ) ) ,  ... 
      'Components should sum due to superposition'  ... 
    ) 
 
 
 
 
 
  
[t p r]  =  cart2sph(1/sqrt(3),1/sqrt(3),1/sqrt(3)); 
magnet_float.magdir  =  [t p] * 180/pi;  % $\vec e_y+\vec e_z+\vec e_z$ 
f1  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
% Components: 
magnet_float.magdir  =  [0  0];  % $\vec e_x$ 
fc1  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
magnet_float.magdir  =  [90  0];  % $\vec e_y$ 
fc2  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
magnet_float.magdir  =  [0  90];  % $\vec e_z$ 
fc3  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
f2  =  (fc1+fc2+fc3)/sqrt(3); 
 
 
 assert (  ... 
      isequal ( chop( f1 , 4 ) , chop ( f2 , 4 ) ) ,  ... 
      'Components should sum due to superposition'  ... 
    ) 
 
 
 
 
 
 
fprintf('passed\n') 
disp('=================') 
 

