 
 
[forces torques]  =  function magnetforces(magnet_fixed, magnet_float, magnet_disp); 
 
 
 
%% MAGNETFORCES  Calculate forces between two cuboid magnets 
% 
% Finish this off later. 
% 
 
 
 
 
 
 
 
 
 
  a1  =  0.5 * magnet_fixed.dim(1); 
  b1  =  0.5 * magnet_fixed.dim(2); 
  c1  =  0.5 * magnet_fixed.dim(3); 
  a2  =  0.5 * magnet_float.dim(1); 
  b2  =  0.5 * magnet_float.dim(2); 
  c2  =  0.5 * magnet_float.dim(3); 
 
  J1r  =  magnet_fixed.magn; 
  J2r  =  magnet_float.magn; 
  J1t  =  magnet_fixed.magdir(1); 
  J2t  =  magnet_float.magdir(1); 
  J1p  =  magnet_fixed.magdir(2); 
  J2p  =  magnet_float.magdir(2); 
 
  dx  =  magnet_disp(1); 
  dy  =  magnet_disp(2); 
  dz  =  magnet_disp(3); 
 
 
 
 
[J1x J1y J1z]  =  sph2cart(J1t,J1p,J1r); 
[J2x J2y J2z]  =  sph2cart(J2t,J2p,J2r); 
 
 
 
 
 
 
 
Rx  =  @(theta) [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)] ; 
Ry  =  @(theta) [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)] ; 
Rz  =  @(theta) [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1] ; 
 
rotate_z_to_x  =  @(vec) Ry(pi/2) * vec' 
rotate_x_to_z  =  @(vec) Ry(-pi/2) * vec' 
 
rotate_z_to_y  =  @(vec) Rx(-pi/2) * vec' 
rotate_y_to_z  =  @(vec) Rx(pi/2) * vec' 
 
 
 
 
force_components_x  =  zeros(3); 
force_components_y  =  zeros(3); 
force_components_z  =  zeros(3); 
 
 
 
if ~( J1x == 0 || J2x == 0 ) 
   
 
[dxr dyr dzr]  =  rotate_x_to_z([dx dy dz]); 
[Fx Fy Fz]  =  forces_parallel(c1,b1,a1,c2,b2,a2,dxr,dyr,dzr,J1x,J2x); 
force_components(1,1,:)  =  rotate_z_to_x([Fx Fy Fz]); 
 
 
end 
if ~( J1y == 0 || J2y == 0 ) 
   
 
[dxr dyr dzr]  =  rotate_y_to_z([dx dy dz]); 
[Fx Fy Fz]  =  forces_parallel(a1,c1,b1,a2,c2,b2,dxr,dyr,dzr,J1y,J2y); 
force_components(2,2,:)  =  rotate_z_to_y([Fx Fy Fz]); 
 
 
end 
if ~( J1z == 0 || J2z == 0 ) 
   
 
[Fx Fy Fz]  =  forces_parallel(a1,b1,c1,a2,b2,c2,dx,dy,dz,J1z,J2z); 
force_components(3,3,:)  =  [Fx Fy Fz]; 
 
 
end 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
end 
 
 
 
 
 
function [Fx Fy Fz]  =  forces_parallel(a,b,c,A,B,C,dx,dy,dz,J,J2) 
% You probably want to call 
%   warning off MATLAB:divideByZero 
%   warning off MATLAB:log:logOfZero 
 
if nargin < 11 
  J2 = J; 
elseif nargin < 10 
  error('Wrong number of input arguments.') 
end 
 
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
 
magconst  =  J * J2/(4 * pi * (4 * pi * 1e-7)); 
Fx  =  magconst * sum(fx(:)); 
Fy  =  magconst * sum(fy(:)); 
Fz  =  magconst * sum(fz(:)); 
 
end 
 
 
 
 
function [Kx Ky Kz]  =  stiffness_parallel(a,b,c,A,B,C,dx,dy,dz,J,J2) 
% You probably want to call 
%   warning off MATLAB:divideByZero 
%   warning off MATLAB:log:logOfZero 
 
if nargin < 11 
  J2 = J; 
elseif nargin < 10 
  error('Wrong number of input arguments.') 
end 
 
[index_h, index_j, index_k, index_l, index_p, index_q]  =  ndgrid([0 1]); 
index_sum  =  (-1).^(index_h+index_j+index_k+index_l+index_p+index_q); 
 
% Using this method is actually less efficient than using six for 
% loops for h..q over [0 1]. To be addressed. 
 
u  =  dx + A * (-1).^index_j - a * (-1).^index_h; 
v  =  dy + B * (-1).^index_l - b * (-1).^index_k; 
w  =  dz + C * (-1).^index_q - c * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
 
k_x  =   ... 
  - r  ... 
  - (u.^2. * v)./(u.^2+w.^2)  ... 
  - v.*log(r-v) ; 
k_y  =   ... 
  - r  ... 
  - (v.^2. * u)./(v.^2+w.^2)  ... 
  - u.*log(r-u) ; 
 
k_z  =  -k_x-k_y; 
 
kx  =  index_sum.*k_x; 
ky  =  index_sum.*k_y; 
kz  =  index_sum.*k_z; 
 
magconst  =  J * J2/(4 * pi * (4 * pi * 1e-7)); 
Kx  =  magconst * sum(kx(:)); 
Ky  =  magconst * sum(ky(:)); 
Kz  =  magconst * sum(kz(:)); 
 
end 
 
 
 
 
function [Fx Fy Fz]  =  forces_orthogonal(a,b,c,A,B,C,dx,dy,dz,J,J2) 
 
if nargin < 11 
  J2 = J; 
elseif nargin < 10 
  error('Wrong number of input arguments.') 
end 
 
[index_h, index_j, index_k, index_l, index_p, index_q]  =  ndgrid([0 1]); 
index_sum  =  (-1).^(index_h+index_j+index_k+index_l+index_p+index_q); 
 
% (Using this method is actually LESS efficient than using six for 
% loops for h..q over [0 1], but it looks a bit nicer, huh?) 
 
u  =  dx + A * (-1).^index_j - a * (-1).^index_h; 
v  =  dy + B * (-1).^index_l - b * (-1).^index_k; 
w  =  dz + C * (-1).^index_q - c * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
f_x  =   ... 
  - v .* w .* ln( r-u )  ... 
  + v .* u .* ln( r+w )  ... 
  + w .* u .* ln( r+v )  ... 
  - 0.5  *  u.^2 .* arctan( v .* w ./ ( u .* r) )  ... 
  - 0.5  *  v.^2 .* arctan( u .* w ./ ( v .* r) )  ... 
  - 0.5  *  w.^2 .* arctan( u .* v ./ ( w .* r) ); 
 
fy  =   ... 
  0.5  *  ( u.^2 - v.^2 ) .* ln( r+w )  ... 
  - u .* w .* ln ( r-u )  ... 
  - u .* v .* arctan( u .* w ./ ( v .* r) )  ... 
  - 0.5  *  w .* r; 
 
fz  =   ... 
  0.5  *  ( u.^2 - w.^2 ) .* ln( r+v )  ... 
  - u .* v .* ln ( r-u )  ... 
  - u .* w .* arctan( u .* v ./ ( w .* r) )  ... 
  - 0.5  *  v .* r; 
 
fx  =  index_sum.*f_x; 
fy  =  index_sum.*f_y; 
fz  =  index_sum.*f_z; 
 
magconst  =  J * J2/(4 * pi * (4 * pi * 1e-7)); 
Fx  =  magconst * sum(fx(:)); 
Fy  =  magconst * sum(fy(:)); 
Fz  =  magconst * sum(fz(:)); 
 
 
 
 
% not yet calculated 
 
 
 
 
 
 

