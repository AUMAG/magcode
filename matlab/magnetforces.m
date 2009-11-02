 
 
function [forces_out]  =  magnetforces(magnet_fixed, magnet_float, displ) 
 
 
 
%% MAGNETFORCES  Calculate forces between two cuboid magnets 
% 
% Finish this off later. 
% 
 
 
 
 
 
 
 
 
 
a1  =  0.5 * magnet_fixed.dim(1); 
b1  =  0.5 * magnet_fixed.dim(2); 
c1  =  0.5 * magnet_fixed.dim(3); 
size1  =  [a1; b1; c1]; 
a2  =  0.5 * magnet_float.dim(1); 
b2  =  0.5 * magnet_float.dim(2); 
c2  =  0.5 * magnet_float.dim(3); 
size2  =  [a2; b2; c2]; 
 
J1r  =  magnet_fixed.magn; 
J2r  =  magnet_float.magn; 
J1t  =  magnet_fixed.magdir(1); 
J2t  =  magnet_float.magdir(1); 
J1p  =  magnet_fixed.magdir(2); 
J2p  =  magnet_float.magdir(2); 
 
if (J1r<0 || J2r<0) 
  error('By convention, magnetisation must be positive; change the angle to reverse direction.') 
end 
 
 
 
 
displ  =  reshape(displ,[3 1]); % column vector 
 
J1  =  [ J1r  *  cosd(J1p)  *  cosd(J1t)  ;  ... 
       J1r  *  cosd(J1p)  *  sind(J1t)  ;  ... 
       J1r  *  sind(J1p) ]; 
 
J2  =  [ J2r  *  cosd(J2p)  *  cosd(J2t)  ;  ... 
       J2r  *  cosd(J2p)  *  sind(J2t)  ;  ... 
       J2r  *  sind(J2p) ]; 
 
 
 
 
force_components  =  repmat(NaN,[9 3]); 
 
 
 
Rx  =  @(theta) [1 0 0; 0 cosd(theta) -sind(theta); 0 sind(theta) cosd(theta)] ; 
Ry  =  @(theta) [cosd(theta) 0 sind(theta); 0 1 0; -sind(theta) 0 cosd(theta)] ; 
Rz  =  @(theta) [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1] ; 
 
rotate_z_to_x  =  @(vec)  Ry( 90) * vec ; 
rotate_x_to_z  =  @(vec)  Ry(-90) * vec ; 
flip_x_z  =  @(vec) abs(rotate_x_to_z(vec)) ; 
 
rotate_z_to_y  =  @(vec)  Rx(-90) * vec ; 
rotate_y_to_z  =  @(vec)  Rx( 90) * vec ; 
flip_y_z  =  @(vec) abs(rotate_y_to_z(vec)) ; 
 
reverse_z  =  @(vec) Ry(180) * vec ; 
reverse_y  =  @(vec) Rz(180) * vec ; 
 
reverse_none  =  @(vec) vec ; 
 
 
 
 
 
size1_rot  =  flip_x_z(size1); 
size2_rot  =  flip_x_z(size2); 
d_rot   =  rotate_x_to_z(displ); 
J1_rot  =  rotate_x_to_z(J1); 
J2_rot  =  rotate_x_to_z(J2); 
 
if J2_rot(2) < 0 
  coord_transform  =  reverse_y 
else 
  coord_transform  =  reverse_none 
end 
 
   d_trans  =  coord_transform( d_rot); 
  J1_trans  =  coord_transform(J1_rot); 
  J2_trans  =  coord_transform(J2_rot); 
 
forces_x_x  =  forces_z_z(size1_rot,size2_rot,d_trans,J1_trans,J2_trans); 
force_components(1,:)  =  rotate_z_to_x( coord_transform(forces_x_x) ); 
 
forces_x_y  =  forces_z_y(size1_rot,size2_rot,d_trans,J1_trans,J2_trans); 
force_components(2,:)  =  rotate_z_to_x( coord_transform(forces_x_y) ); 
 
forces_x_z  =  forces_z_y(size1_rot,size2_rot,d_trans,J1_trans,J2_trans); 
force_components(3,:)  =  rotate_z_to_x( coord_transform(forces_x_z) ); 
 
 
 
 
size1_rot  =  flip_y_z(size1); 
size2_rot  =  flip_y_z(size2); 
d_rot      =  rotate_y_to_z( displ ); 
J1_rot     =  rotate_y_to_z( J1    ); 
J2_rot     =  rotate_y_to_z( J2    ); 
 
if J2_rot(2) < 0 
  coord_transform  =  reverse_y ; 
else 
  coord_transform  =  reverse_none ; 
end 
 
   d_trans  =  coord_transform( d_rot); 
  J1_trans  =  coord_transform(J1_rot); 
  J2_trans  =  coord_transform(J2_rot); 
 
forces_y_x  =  forces_z_x(size1_rot,size2_rot,d_trans,J1_trans,J2_trans); 
force_components(4,:)  =  rotate_z_to_y( coord_transform(forces_y_x) ); 
 
forces_y_y  =  forces_z_z(size1_rot,size2_rot,d_trans,J1_trans,J2_trans); 
force_components(5,:)  =  rotate_z_to_y( coord_transform(forces_y_y) ); 
 
forces_y_z  =  forces_z_y(size1_rot,size2_rot,d_trans,J1_trans,J2_trans); 
force_components(6,:)  =  rotate_z_to_y( coord_transform(forces_y_z) ); 
 
 
 
 
 
if J1(3) < 0 
  coord_transform  =  reverse_z 
else 
  coord_transform  =  reverse_none 
end 
 
   d_trans  =  coord_transform(displ); 
  J1_trans  =  coord_transform(J1); 
  J2_trans  =  coord_transform(J2); 
 
forces_z_z  =  forces_z_z( size1,size2,d_trans,J1_trans,J2_trans ); 
force_components(7,:)  =  coord_transform(forces_z_z); 
 
forces_z_y  =  forces_z_y( size1,size2,d_trans,J1_trans,J2_trans ); 
force_components(8,:)  =  coord_transform(forces_z_y); 
 
forces_z_x  =  forces_z_x( size1,size2,d_trans,J1_trans,J2_trans ); 
force_components(9,:)  =  coord_transform(forces_z_x); 
 
 
 
 
 
disp('  ') 
disp('CALCULATING FORCES') 
disp('==================') 
disp('Displacement:') 
disp(displ) 
disp('Magnetisations:') 
disp(J1) 
disp(J2) 
 
 
 
disp('Forces z-x:') 
disp(forces_z_x) 
disp('Forces z-y:') 
disp(forces_z_y) 
disp('Forces z-z:') 
disp(forces_z_z) 
 
 
 
disp('Forces x-x:') 
disp(forces_x_x) 
disp('Forces x-y:') 
disp(forces_x_y) 
disp('Forces x-z:') 
disp(forces_x_z) 
 
 
 
disp('Forces y-x:') 
disp(forces_y_x) 
disp('Forces y-y:') 
disp(forces_y_y) 
disp('Forces y-z:') 
disp(forces_y_z) 
 
 
 
 
forces_out  =  sum(force_components); 
 
 
 
end 
 
 
 
 
 
function forces_xyz  =  forces_z_z(size1,size2,offset,J1,J2) 
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
  forces_xyz  =  [0; 0; 0]; 
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
 
% (Using this vectorised method is actually less efficient than using six |for| 
% loops over |[0, 1]|. To be addressed.) 
 
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
forces_xyz  =  magconst.*[ sum(fx(:)) sum(fy(:)) sum(fz(:)) ] ; 
 
end 
 
 
 
 
function forces_xyz  =  forces_z_x(size1,size2,offset,J1,J2) 
 
swap_x_y  =  @(vec) [vec(2); vec(1); vec(3)]; 
 
forces_xyz  =  forces_z_y( ... 
  swap_x_y(size1), swap_x_y(size2), swap_x_y(offset), ... 
  J1, swap_x_y(J2) ); 
 
forces_xyz  =  swap_x_y( forces_xyz ); 
 
end 
 
function forces_xyz  =  forces_z_y(size1,size2,offset,J1,J2) 
 
if length(J1) == 3 
  J1  =  J1(3); 
end 
if length(J2) == 3 
  J2  =  J2(2); 
end 
 
 
 
if (J1==0 || J2==0) 
  forces_xyz  =  [0; 0; 0]; 
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
 
% (Using this vectorised method is actually less efficient than using six |for| 
% loops over |[0, 1]|. To be addressed.) 
 
u  =  dx + A * (-1).^index_j - a * (-1).^index_h; 
v  =  dy + B * (-1).^index_l - b * (-1).^index_k; 
w  =  dz + C * (-1).^index_q - c * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
 
 
 
 
 
if (J1<0 || J2<0) 
  error('Positive magnetisations only!') 
end 
 
f_x  =   ... 
  - v .* w .* log( r-u )  ... 
  + v .* u .* log( r+w )  ... 
  + u .* w .* log( r+v )  ... 
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
forces_xyz  =  magconst.*[ sum(fx(:)) sum(fy(:)) sum(fz(:)) ] ; 
 
 
end 
 
 
 
 
 

