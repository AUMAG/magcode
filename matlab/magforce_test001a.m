 
 
f  =  []; 
 
magnet_fixed.dim  =  [0.04 0.04 0.04]; 
magnet_float.dim  =     magnet_fixed.dim; 
 
magnet_fixed.magn  =  1.3; 
magnet_float.magn  =    magnet_fixed.magn; 
magnet_fixed.magdir  =  [0 90]; % vertical 
magnet_float.magdir  =  magnet_fixed.magdir; 
 
displ  =  [0 0 0.1]; 
 
magnet_float.magn  =  1.3; 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,displ+eps); 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,-displ+eps); 
 
magnet_float.magn  =  -1.3; 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,displ+eps); 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,-displ+eps); 
 
magnet_fixed.magdir  =  [0 0]; % x 
magnet_float.magdir  =  magnet_fixed.magdir; 
displ  =  [0.1 0 0]; 
 
magnet_float.magn  =  1.3; 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,displ+eps); 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,-displ+eps); 
 
magnet_float.magn  =  -1.3; 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,displ+eps); 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,-displ+eps); 
 
magnet_fixed.magdir  =  [90 0]; % y 
magnet_float.magdir  =  magnet_fixed.magdir; 
displ  =  [0 0.1 0]; 
 
magnet_float.magn  =  1.3; 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,displ+eps); 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,-displ+eps); 
 
magnet_float.magn  =  -1.3; 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,displ+eps); 
f(:,end+1)  =  magnetforces(magnet_fixed,magnet_float,-displ+eps); 
 
assert( chop( f(3,1) , 6 ) == - chop ( f(3,2) , 6) ); 
assert( chop( f(3,2) , 6 ) ==   chop ( f(3,3) , 6) ); 
assert( chop( f(3,3) , 6 ) == - chop ( f(3,4) , 6) ); 
 
disp('Tests passed'); 
 

