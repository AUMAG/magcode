 
 
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
  error(['By convention, magnetisation must be positive; ',  ... 
         'change the angle to reverse direction.']) 
end 
 
 
 
 
swap_x_y  =  @(vec) vec([2 1 3]); 
swap_x_z  =  @(vec) vec([3 2 1]); 
swap_y_z  =  @(vec) vec([1 3 2]); 
 
Rx  =  @(theta) [1 0 0; 0 cosd(theta) -sind(theta); 0 sind(theta) cosd(theta)] ; 
Ry  =  @(theta) [cosd(theta) 0 sind(theta); 0 1 0; -sind(theta) 0 cosd(theta)] ; 
Rz  =  @(theta) [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1] ; 
 
Rx_180  =  Rx(180); 
Rx_090  =  Rx( 90); 
Rx_270  =  Rx(-90); 
Ry_180  =  Ry(180); 
Ry_090  =  Ry( 90); 
Ry_270  =  Ry(-90); 
Rz_180  =  Rz(180); 
 
identity_function  =  @(inp) inp; 
 
rotate_round_x  =  @(vec) Rx_180 * vec ; 
rotate_round_y  =  @(vec) Ry_180 * vec ; 
rotate_round_z  =  @(vec) Rz_180 * vec ; 
rotate_none  =  identity_function ; 
 
rotate_z_to_x  =  @(vec)  Ry_090 * vec ; 
rotate_x_to_z  =  @(vec)  Ry_270 * vec ; 
 
rotate_z_to_y  =  @(vec)  Rx_090 * vec ; 
rotate_y_to_z  =  @(vec)  Rx_270 * vec ; 
 
 
 
 
displ  =  reshape(displ,[3 1]); % column vector 
 
J1  =  [ J1r  *  cosd(J1p)  *  cosd(J1t)  ;  ... 
       J1r  *  cosd(J1p)  *  sind(J1t)  ;  ... 
       J1r  *  sind(J1p) ]; 
 
J2  =  [ J2r  *  cosd(J2p)  *  cosd(J2t)  ;  ... 
       J2r  *  cosd(J2p)  *  sind(J2t)  ;  ... 
       J2r  *  sind(J2p) ]; 
 
 
 
 
force_components  =  repmat(NaN,[9 3]); 
 
 
 
disp('  ') 
disp('CALCULATING FORCES') 
disp('==================') 
disp('Displacement:') 
disp(displ') 
disp('Magnetisations:') 
disp(J1') 
disp(J2') 
 
 
 
 
 
size1_rot  =  swap_x_z(size1); 
size2_rot  =  swap_x_z(size2); 
d_rot   =  rotate_x_to_z(displ); 
J1_rot  =  rotate_x_to_z(J1); 
J2_rot  =  rotate_x_to_z(J2); 
 
disp('Forces x-x:') 
forces_x_x  =  forces_calc_z_z(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(1,:)  =  rotate_z_to_x( forces_x_x ); 
 
disp('Forces x-y:') 
forces_x_y  =  forces_calc_z_y(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(2,:)  =  rotate_z_to_x( forces_x_y ); 
 
disp('Forces x-z:') 
forces_x_z  =  forces_calc_z_y(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(3,:)  =  rotate_z_to_x( forces_x_z ); 
 
 
 
 
 
size1_rot  =  swap_y_z(size1); 
size2_rot  =  swap_y_z(size2); 
d_rot      =  rotate_y_to_z( displ ); 
J1_rot     =  rotate_y_to_z( J1    ); 
J2_rot     =  rotate_y_to_z( J2    ); 
 
disp('Forces y-x:') 
forces_y_x  =  forces_calc_z_x(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(4,:)  =  rotate_z_to_y( forces_y_x ); 
 
disp('Forces y-y:') 
forces_y_y  =  forces_calc_z_z(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(5,:)  =  rotate_z_to_y( forces_y_y ); 
 
disp('Forces y-z:') 
forces_y_z  =  forces_calc_z_y(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(6,:)  =  rotate_z_to_y( forces_y_z ); 
 
 
 
 
 
 
disp('Forces z-z:') 
forces_z_z  =  forces_calc_z_z( size1,size2,displ,J1,J2 ); 
force_components(7,:)  =  forces_z_z; 
 
disp('Forces z-y:') 
forces_z_y  =  forces_calc_z_y( size1,size2,displ,J1,J2 ); 
force_components(8,:)  =  forces_z_y; 
 
disp('Forces z-x:') 
forces_z_x  =  forces_calc_z_x( size1,size2,displ,J1,J2 ); 
force_components(9,:)  =  forces_z_x; 
 
 
 
 
forces_out  =  sum(force_components); 
 
 
 
 
 
 
function forces_xyz  =  forces_calc_z_z(size1,size2,offset,J1,J2) 
% You probably want to call 
%   warning off MATLAB:divideByZero 
%   warning off MATLAB:log:logOfZero 
 
J1  =  J1(3); 
J2  =  J2(3); 
 
if (J1==0 || J2==0) 
  disp('Zero magnetisation.') 
  forces_xyz   =   [0; 0; 0]; 
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
 
% (Using this vectorised method is less efficient than using six |for| 
% loops over |[0, 1]|. To be addressed.) 
 
u  =  dx + A * (-1).^index_j - a * (-1).^index_h; 
v  =  dy + B * (-1).^index_l - b * (-1).^index_k; 
w  =  dz + C * (-1).^index_q - c * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
 
 
 
f_x  =   ... 
  + 0.5 * (v.^2-w.^2).*log(r-u)  ... 
  + u.*v.*log(r-v)  ... 
  + v.*w.*atan2(u.*v,r.*w)  ... 
  + 0.5 * r.*u; 
 
f_y  =   ... 
  + 0.5 * (u.^2-w.^2).*log(r-v)  ... 
  + u.*v.*log(r-u)  ... 
  + u.*w.*atan2(u.*v,r.*w) ... 
  + 0.5 * r.*v; 
 
f_z  =   ... 
  - u.*w.*log(r-u)  ... 
  - v.*w.*log(r-v)  ... 
  + u.*v.*atan2(u.*v,r.*w)  ... 
  - r.*w; 
 
fx  =  index_sum.*f_x; 
fy  =  index_sum.*f_y; 
fz  =  index_sum.*f_z; 
 
magconst  =  J1 * J2/(4 * pi * (4 * pi * 1e-7)); 
forces_xyz  =  magconst.*[ sum(fx(:)) ; sum(fy(:)) ; sum(fz(:)) ] ; 
 
disp(forces_xyz') 
 
end 
 
 
 
 
function forces_xyz  =  forces_calc_z_y(size1,size2,offset,J1,J2) 
 
J1m  =  J1(3); 
J2m  =  J2(2); 
 
if (J1m==0 || J2m==0) 
  disp('Zero magnetisation.') 
  forces_xyz  =  [0; 0; 0]; 
  return; 
end 
 
if     ( J1m>0 && J2m>0 ) 
  rotate_transform  =  rotate_none; 
elseif ( J1m<0 && J2m>0 ) 
  rotate_transform  =  rotate_round_y; 
elseif ( J1m>0 && J2m<0 ) 
  rotate_transform  =  rotate_round_z; 
elseif ( J1m<0 && J2m<0 ) 
  rotate_transform  =  rotate_round_x; 
end 
 
forces_tmp  =  forces_calc_z_y_plusplus(  ... 
            size1,size2,  ... 
         rotate_transform(offset),  ... 
         rotate_transform(J1),      ... 
         rotate_transform(J2)       ... 
          ); 
forces_xyz  =  rotate_transform( forces_tmp); 
disp(forces_xyz') 
 
end 
 
 
 
function forces_xyz  =  forces_calc_z_x(size1,size2,offset,J1,J2) 
 
forces_xyz  =  forces_calc_z_y( ... 
  swap_x_y(size1), swap_x_y(size2), swap_x_y(offset), ... 
  J1, swap_x_y(J2) ); 
 
forces_xyz  =  swap_x_y( forces_xyz ); 
 
end 
 
 
 
 
 
function forces_xyz  =  forces_calc_z_y_plusplus(size1,size2,offset,J1,J2) 
 
J1  =  J1(3); 
J2  =  J2(2); 
 
if (J1<0 || J2<0) 
  error('Positive magnetisations only!') 
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
 
% (Using this vectorised method is less efficient than using six |for| 
% loops over |[0, 1]|. To be addressed.) 
 
u  =  dx + A * (-1).^index_j - a * (-1).^index_h; 
v  =  dy + B * (-1).^index_l - b * (-1).^index_k; 
w  =  dz + C * (-1).^index_q - c * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
 
 
 
f_x  =   ... 
  - multiply_x_log_y ( v .* w , r-u )  ... 
  + multiply_x_log_y ( v .* u , r+w )  ... 
  + multiply_x_log_y ( u .* w , r+v )  ... 
  - 0.5  *  u.^2 .* atan1( v .* w , u .* r )  ... 
  - 0.5  *  v.^2 .* atan1( u .* w , v .* r )  ... 
  - 0.5  *  w.^2 .* atan1( u .* v , w .* r ); 
 
f_y  =   ... 
  0.5  *  multiply_x_log_y( u.^2 - v.^2 , r+w )  ... 
  - multiply_x_log_y( u .* w , r-u )  ... 
  - u .* v .* atan1( u .* w , v .* r )  ... 
  - 0.5  *  w .* r; 
 
f_z  =   ... 
  0.5  *  ( u.^2 - w.^2 ) .* log( r+v )  ... 
  - multiply_x_log_y( u .* v , r-u )  ... 
  - u .* w .* atan1( u .* v , w .* r )  ... 
  - 0.5  *  v .* r; 
 
f_x  =  index_sum.*f_x; 
f_y  =  index_sum.*f_y; 
f_z  =  index_sum.*f_z; 
 
forces_xyz  =  J1 * J2/(4 * pi * (4 * pi * 1e-7)) .*  ... 
  [ sum(f_x(:)) ; sum(f_y(:)) ; sum(f_z(:)) ] ; 
 
end 
 
 
 
 
function out  =  multiply_x_log_y(x,y) 
  out  =  x.*log(y); 
  out(isnan(out)) = 0; 
end 
 
 
 
function out  =  atan1(x,y) 
  out  =  zeros(size(x)); 
  ind  =  x~=0 & y~=0; 
  out(ind)  =  atan(x(ind)./y(ind)); 
end 
 
 
 
 
 
end 
 

