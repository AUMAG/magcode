 
  [forces torques]  =  function magnetforces(magnet_fixed, magnet_float, magnet_disp); 
 
   
 
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
 
 
   
 
J1x  =  J1r * cos(J1p) * sin(J1t); 
J1y  =  J1r * sin(J1p) * sin(J1t); 
J1z  =  J1r * cos(J1t); 
 
J2x  =  J2r * cos(J2p) * sin(J2t); 
J2y  =  J2r * sin(J2p) * sin(J2t); 
J2z  =  J2r * cos(J2t); 
 
 
   
 
 
   
 
 
   
 
 
 
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
 
 
 
 
 
  
 
 
 

