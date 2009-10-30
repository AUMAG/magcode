 
 
function [forces_out]  =  magnetforces(magnet_fixed, magnet_float, displ) 
 
 
 
%% MAGNETFORCES  Calculate forces between two cuboid magnets 
% 
% Finish this off later. 
% 
 
 
 
 
 
 
 
 
 
a1  =  0.5 * magnet_fixed.dim(1); 
b1  =  0.5 * magnet_fixed.dim(2); 
c1  =  0.5 * magnet_fixed.dim(3); 
size1  =  [a1 b1 c1]; 
a2  =  0.5 * magnet_float.dim(1); 
b2  =  0.5 * magnet_float.dim(2); 
c2  =  0.5 * magnet_float.dim(3); 
size2  =  [a2 b2 c2]; 
 
J1r  =  magnet_fixed.magn; 
J2r  =  magnet_float.magn; 
J1t  =  magnet_fixed.magdir(1); 
J2t  =  magnet_float.magdir(1); 
J1p  =  magnet_fixed.magdir(2); 
J2p  =  magnet_float.magdir(2); 
 
 
 
 
J1  =  [ J1r  *  cosd(J1p)  *  cosd(J1t) ,  ... 
       J1r  *  cosd(J1p)  *  sind(J1t) ,  ... 
       J1r  *  sind(J1p) ]; 
 
J2  =  [ J2r  *  cosd(J2p)  *  cosd(J2t) ,  ... 
       J2r  *  cosd(J2p)  *  sind(J2t) ,  ... 
       J2r  *  sind(J2p) ]; 
 
 
 
 
force_components  =  repmat(NaN,[3 3 3]); 
 
 
 
Rx  =  @(theta) [1 0 0; 0 cosd(theta) -sind(theta); 0 sind(theta) cosd(theta)] ; 
Ry  =  @(theta) [cosd(theta) 0 sind(theta); 0 1 0; -sind(theta) 0 cosd(theta)] ; 
Rz  =  @(theta) [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1] ; 
 
rotate_z_to_x  =  @(vec)  Ry( 90) * vec' ; 
rotate_x_to_z  =  @(vec)  Ry(-90) * vec' ; 
flip_x_z  =  @(vec) abs(rotate_x_to_z(vec)) ; 
 
rotate_z_to_y  =  @(vec)  Rx(-90) * vec' ; 
rotate_y_to_z  =  @(vec)  Rx( 90) * vec' ; 
flip_y_z  =  @(vec) abs(rotate_y_to_z(vec)) ; 
 
 
 
disp('x-x'); 
 
 
size1rot  =  flip_x_z(size1); 
size2rot  =  flip_x_z(size2); 
drot   =  rotate_x_to_z(displ); 
J1rot  =  rotate_x_to_z(J1); 
J2rot  =  rotate_x_to_z(J2); 
forces_xyz  =  forces_parallel(size1rot,size2rot,drot,J1rot,J2rot); 
force_components(1,1,:)  =  rotate_z_to_x(forces_xyz); 
 
 
 
disp('x-y'); 
 
 
size1rot  =  flip_x_z(size1); 
size2rot  =  flip_x_z(size2); 
drot   =  rotate_x_to_z( displ ); 
J1rot  =  rotate_x_to_z(J1); 
J2rot  =  rotate_x_to_z(J2); 
forces_xyz  =  forces_orthogonal_y(size1rot,size2rot,drot,J1rot,J2rot); 
force_components(1,2,:)  =  rotate_z_to_x(forces_xyz); 
 
 
 
disp('x-z'); 
 
 
size1rot  =  flip_x_z(size1); 
size2rot  =  flip_x_z(size2); 
drot   =  rotate_x_to_z( displ ); 
J1rot  =  rotate_x_to_z(J1); 
J2rot  =  rotate_x_to_z(J2); 
forces_xyz  =  forces_orthogonal_y(size1rot,size2rot,drot,J1rot,J2rot); 
force_components(1,3,:)  =  rotate_z_to_x(forces_xyz); 
 
 
 
 
disp('y-x'); 
 
 
size1rot  =  flip_x_z(size1); 
size2rot  =  flip_x_z(size2); 
drot      =  rotate_y_to_z( displ ); 
J1rot     =  rotate_y_to_z( J1    ); 
J2rot     =  rotate_y_to_z( J2    ); 
forces_xyz  =  forces_orthogonal_x(size1rot,size2rot,drot,J1rot,J2rot); 
force_components(2,1,:)  =  rotate_z_to_y(forces_xyz); 
 
 
 
disp('y-y'); 
 
 
size1rot  =  flip_y_z(size1); 
size2rot  =  flip_y_z(size2); 
drot   =  rotate_y_to_z( displ ); 
J1rot  =  rotate_y_to_z(J1); 
J2rot  =  rotate_y_to_z(J2); 
forces_xyz  =  forces_parallel(size1rot,size2rot,drot,J1rot,J2rot); 
force_components(2,2,:)  =  rotate_z_to_y(forces_xyz); 
 
 
 
disp('y-z'); 
 
 
size1rot  =  flip_y_z(size1); 
size2rot  =  flip_y_z(size2); 
drot   =  rotate_y_to_z( displ ); 
J1rot  =  rotate_y_to_z(J1); 
J2rot  =  rotate_y_to_z(J2); 
forces_xyz  =  forces_orthogonal_y(size1rot,size2rot,drot,J1rot,J2rot); 
force_components(2,3,:)  =  rotate_z_to_y(forces_xyz); 
 
 
 
disp('z-x'); 
 
 
forces_xyz  =  forces_orthogonal_x( size1,size2,displ,J1,J2 ); 
force_components(3,1,:)  =  forces_xyz; 
 
 
disp('z-y'); 
 
 
forces_xyz  =  forces_orthogonal_y( size1,size2,displ,J1,J2 ); 
force_components(3,2,:)  =  forces_xyz; 
 
 
disp('z-z'); 
 
 
forces_xyz  =  forces_parallel(size1,size2,displ,J1,J2); 
force_components(3,3,:)  =  forces_xyz; 
 
 
 
forces_out  =  squeeze(sum(sum(force_components,1),2)); 
 
 
 
end 
 
 
 
 
 
function forces_xyz  =  forces_parallel(size1,size2,offset,J1,J2) 
% You probably want to call 
%   warning off MATLAB:divideByZero 
%   warning off MATLAB:log:logOfZero 
 
if length(J1) == 3 
  J1  =  J1(3); 
end 
if length(J2) == 3 
  J2  =  J2(3); 
end 
 
if (J1==0 || J2==0) 
  disp('Zero magnetisation (parallel)') 
  forces_xyz  =  [0 0 0]; 
  return; 
end 
 
dx  =  offset(1); 
dy  =  offset(2); 
dz  =  offset(3); 
 
a  =  size1(1); 
b  =  size1(2); 
c  =  size1(3); 
A  =  size2(1); 
B  =  size2(2); 
C  =  size2(3); 
 
[index_h, index_j, index_k, index_l, index_p, index_q]  =  ndgrid([0 1]); 
index_sum  =  (-1).^(index_h+index_j+index_k+index_l+index_p+index_q); 
 
% (Using this method is actually LESS efficient than using six for 
% loops for h..q over [0 1], but it looks a bit nicer, huh?) 
 
u  =  dx + A * (-1).^index_j - a * (-1).^index_h; 
v  =  dy + B * (-1).^index_l - b * (-1).^index_k; 
w  =  dz + C * (-1).^index_q - c * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
f_x  =   ... 
  + 0.5 * (v.^2-w.^2).*log(r-u)  ... 
  + u.*v.*log(r-v)  ... 
  + v.*w.*atan(u.*v./r./w)  ... 
  + 0.5 * r.*u; 
 
f_y  =   ... 
  + 0.5 * (u.^2-w.^2).*log(r-v)  ... 
  + u.*v.*log(r-u)  ... 
  + u.*w.*atan(u.*v./r./w) ... 
  + 0.5 * r.*v; 
 
f_z  =   ... 
  - u.*w.*log(r-u)  ... 
  - v.*w.*log(r-v)  ... 
  + u.*v.*atan(u.*v./r./w)  ... 
  - r.*w; 
 
fx  =  index_sum.*f_x; 
fy  =  index_sum.*f_y; 
fz  =  index_sum.*f_z; 
 
magconst  =  J1 * J2/(4 * pi * (4 * pi * 1e-7)); 
forces_xyz  =  magconst.*[ sum(fx(:)) sum(fy(:)) sum(fz(:)) ] 
 
end 
 
 
 
 
function forces_xyz  =  forces_orthogonal_x(size1,size2,offset,J1,J2) 
 
swap_x_y  =  @(vec) [vec(2) vec(1) vec(3)]; 
 
forces_xyz  =  forces_orthogonal_y( ... 
  swap_x_y(size1), swap_x_y(size2), swap_x_y(offset), ... 
  J1, swap_x_y(J2) ); 
 
forces_xyz  =  swap_x_y( forces_xyz ); 
 
end 
 
function forces_xyz  =  forces_orthogonal_y(size1,size2,offset,J1,J2) 
 
if length(J1) == 3 
  J1  =  J1(3); 
end 
if length(J2) == 3 
  J2  =  J2(2); 
end 
 
if (J1==0 || J2==0) 
  disp('Zero magnetisation (orth)') 
  forces_xyz  =  [0 0 0]; 
  return; 
end 
 
dx  =  offset(1); 
dy  =  offset(2); 
dz  =  offset(3); 
 
a  =  size1(1); 
b  =  size1(2); 
c  =  size1(3); 
A  =  size2(1); 
B  =  size2(2); 
C  =  size2(3); 
 
[index_h, index_j, index_k, index_l, index_p, index_q]  =  ndgrid([0 1]); 
index_sum  =  (-1).^(index_h+index_j+index_k+index_l+index_p+index_q); 
 
% (Using this method is actually LESS efficient than using six for 
% loops for h..q over [0 1], but it looks a bit nicer, huh?) 
 
u  =  dx + A * (-1).^index_j - a * (-1).^index_h; 
v  =  dy + B * (-1).^index_l - b * (-1).^index_k; 
w  =  dz + C * (-1).^index_q - c * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
f_x  =   ... 
  - v .* w .* log( r-u )  ... 
  + v .* u .* log( r+w )  ... 
  + w .* u .* log( r+v )  ... 
  - 0.5  *  u.^2 .* atan( v .* w ./ ( u .* r) )  ... 
  - 0.5  *  v.^2 .* atan( u .* w ./ ( v .* r) )  ... 
  - 0.5  *  w.^2 .* atan( u .* v ./ ( w .* r) ); 
 
f_y  =   ... 
  0.5  *  ( u.^2 - v.^2 ) .* log( r+w )  ... 
  - u .* w .* log ( r-u )  ... 
  - u .* v .* atan( u .* w ./ ( v .* r) )  ... 
  - 0.5  *  w .* r; 
 
f_z  =   ... 
  0.5  *  ( u.^2 - w.^2 ) .* log( r+v )  ... 
  - u .* v .* log ( r-u )  ... 
  - u .* w .* atan( u .* v ./ ( w .* r) )  ... 
  - 0.5  *  v .* r; 
 
fx  =  index_sum.*f_x; 
fy  =  index_sum.*f_y; 
fz  =  index_sum.*f_z; 
 
magconst  =  J1 * J2/(4 * pi * (4 * pi * 1e-7)); 
forces_xyz  =  magconst.*[ sum(fx(:)) sum(fy(:)) sum(fz(:)) ] 
 
 
end 
 
 
 
 
 
 
 

