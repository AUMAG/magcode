 
 
clc; 
 
f  =  []; 
 
magnet_fixed.dim  =  [0.04 0.04 0.04]; 
magnet_float.dim  =   magnet_fixed.dim; 
 
magnet_fixed.magn  =  1.3; 
magnet_float.magn  =  1.3; 
 
magnet_fixed.magdir  =  [0  90];  % $z$ 
for ii  =  [1, -1] 
  for jj  =  [1, -1] 
 
    magnet_float.magdir  =  ii * [90 0];  % $\pm y$ 
    displ  =  jj * [1e-12 1e-12 0.1];  % $\pm z$ 
    f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
    pause 
 
  end 
end 
 
 
% chop: 
f( abs(f)<1e-10 )  =  0; 
