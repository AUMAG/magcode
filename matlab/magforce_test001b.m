 
 
 
f  =  []; 
 
magnet_fixed.dim  =  [0.04 0.04 0.04]; 
magnet_float.dim  =   magnet_fixed.dim; 
 
magnet_fixed.magn  =  1.3; 
magnet_float.magn  =  1.3; 
magnet_fixed.magdir  =  [0  0]; % x 
 
for ii  =  [1 -1] 
 
magnet_float.magdir  =  ii * [0 90]; % z 
displ  =  ii * [0 0 0.1]; 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,displ+eps); 
 
magnet_float.magdir  =  ii * [90 0]; % y 
displ  =  ii * [0 0.1 0]; 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,displ+eps); 
 
end 
 
 
