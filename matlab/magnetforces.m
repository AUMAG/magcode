 
 
function [varargout]  =  magnetforces(magnet_fixed, magnet_float, displ, varargin) 
 
  
%% MAGNETFORCES  Calculate forces between two cuboid magnets 
% 
% Finish this off later. 
% 
 
 
 
 
 
  
if size(displ,1) == 3 
  % all good 
elseif size(displ,2) == 3 
  displ  =  transpose(displ); 
else 
  error('Displacements matrix should be of size (3, D) where D is the number of displacements.') 
end 
 
Ndispl  =  size(displ,2); 
 
Nvargin  =  length(varargin); 
debug_disp  =  @(str) disp([]); 
calc_force_bool  =  false; 
calc_stiffness_bool  =  false; 
 
for ii  =  1:Nvargin 
  switch varargin{ii} 
    case 'debug' 
      debug_disp  =  @(str) disp(str); 
    case 'force' 
      calc_force_bool  =  true; 
    case 'stiffness' 
      calc_stiffness_bool  =  true; 
    otherwise 
      error(['Unknown calculation option ''',varargin{ii},'''']) 
  end 
end 
 
if ~calc_force_bool && ~calc_stiffness_bool 
  calc_force_bool  =  true; 
end 
 
if calc_force_bool 
  forces_out  =  repmat(NaN,[3 Ndispl]); 
end 
if calc_stiffness_bool 
  stiffnesses_out  =  repmat(NaN,[3 Ndispl]); 
end 
 
 
  
size1  =  reshape(magnet_fixed.dim/2,[3 1]); 
size2  =  reshape(magnet_float.dim/2,[3 1]); 
 
if length(magnet_fixed.magdir)==2 
  J1r  =  magnet_fixed.magn; 
  J1t  =  magnet_fixed.magdir(1); 
  J1p  =  magnet_fixed.magdir(2); 
  J1   =  [ J1r  *  cosd(J1p)  *  cosd(J1t)  ;  ... 
          J1r  *  cosd(J1p)  *  sind(J1t)  ;  ... 
          J1r  *  sind(J1p) ]; 
else 
  if all(magnet_fixed.magdir == [0 0 0]) 
    J1  =  [0; 0; 0]; 
  else 
    J1  =  magnet_fixed.magn * magnet_fixed.magdir/norm(magnet_fixed.magdir); 
    J1  =  reshape(J1,[3 1]); 
  end 
end 
 
if length(magnet_float.magdir)==2 
  J2r  =  magnet_float.magn; 
  J2t  =  magnet_float.magdir(1); 
  J2p  =  magnet_float.magdir(2); 
  J2   =  [ J2r  *  cosd(J2p)  *  cosd(J2t)  ;  ... 
          J2r  *  cosd(J2p)  *  sind(J2t)  ;  ... 
          J2r  *  sind(J2p) ]; 
else 
  if all(magnet_float.magdir == [0 0 0]) 
    J2  =  [0; 0; 0]; 
  else 
    J2  =  magnet_float.magn * magnet_float.magdir/norm(magnet_float.magdir); 
    J2  =  reshape(J2,[3 1]); 
  end 
end 
 
 
  
magconst  =  1/(4 * pi * (4 * pi * 1e-7)); 
 
[index_i, index_j, index_k, index_l, index_p, index_q]  =  ndgrid([0 1]); 
 
index_sum  =  (-1).^(index_i+index_j+index_k+index_l+index_p+index_q); 
 
 
 
  
swap_x_z  =  @(vec) vec([3 2 1]); 
swap_y_z  =  @(vec) vec([1 3 2]); 
 
rotate_z_to_x  =  @(vec)  [0 0  1; 0 1 0; -1 0 0] * vec ; % Ry( 90) 
rotate_x_to_z  =  @(vec)  [0 0 -1; 0 1 0;  1 0 0] * vec ; % Ry(-90) 
 
rotate_y_to_z  =  @(vec)  [1 0 0; 0 0 -1; 0  1 0] * vec ; % Rx( 90) 
rotate_z_to_y  =  @(vec)  [1 0 0; 0 0  1; 0 -1 0] * vec ; % Rx(-90) 
 
rotate_x_to_y  =  @(vec)  [0 -1 0;  1 0 0; 0 0 1] * vec ; % Rz( 90) 
rotate_y_to_x  =  @(vec)  [0  1 0; -1 0 0; 0 0 1] * vec ; % Rz(-90) 
 
 
 
 
  
if calc_force_bool 
  for ii  =  1:Ndispl 
    forces_out(:,ii)   =   single_magnet_force(displ(:,ii)); 
  end 
end 
 
if calc_stiffness_bool 
  for ii  =  1:Ndispl 
    stiffnesses_out(:,ii)   =   single_magnet_stiffness(displ(:,ii)); 
  end 
end 
 
 
  
varargout{1}  =  forces_out; 
for ii  =  1:Nvargin 
  switch varargin{ii} 
    case 'force' 
      varargout{ii}  =  forces_out; 
    case 'stiffness' 
      varargout{ii}  =  stiffnesses_out; 
  end 
end 
 
 
 
 
  
function force_out  =  single_magnet_force(displ) 
 
force_components  =  repmat(NaN,[9 3]); 
 
  
debug_disp('  ') 
debug_disp('CALCULATING THINGS') 
debug_disp('==================') 
debug_disp('Displacement:') 
debug_disp(displ') 
debug_disp('Magnetisations:') 
debug_disp(J1') 
debug_disp(J2') 
 
 
  
size1_rot  =  swap_x_z(size1); 
size2_rot  =  swap_x_z(size2); 
d_rot   =  rotate_x_to_z(displ); 
J1_rot  =  rotate_x_to_z(J1); 
J2_rot  =  rotate_x_to_z(J2); 
 
debug_disp('Forces x-x:') 
forces_x_x  =  forces_calc_z_z(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(1,:)  =  rotate_z_to_x( forces_x_x ); 
 
debug_disp('Forces x-y:') 
forces_x_y  =  forces_calc_z_y(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(2,:)  =  rotate_z_to_x( forces_x_y ); 
 
debug_disp('Forces x-z:') 
forces_x_z  =  forces_calc_z_x(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(3,:)  =  rotate_z_to_x( forces_x_z ); 
 
 
  
size1_rot  =  swap_y_z(size1); 
size2_rot  =  swap_y_z(size2); 
d_rot      =  rotate_y_to_z( displ ); 
J1_rot     =  rotate_y_to_z( J1    ); 
J2_rot     =  rotate_y_to_z( J2    ); 
 
debug_disp('Forces y-x:') 
forces_y_x  =  forces_calc_z_x(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(4,:)  =  rotate_z_to_y( forces_y_x ); 
 
debug_disp('Forces y-y:') 
forces_y_y  =  forces_calc_z_z(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(5,:)  =  rotate_z_to_y( forces_y_y ); 
 
debug_disp('Forces y-z:') 
forces_y_z  =  forces_calc_z_y(size1_rot,size2_rot,d_rot,J1_rot,J2_rot); 
force_components(6,:)  =  rotate_z_to_y( forces_y_z ); 
 
 
  
debug_disp('z-z force:') 
force_components(9,:)  =  forces_calc_z_z( size1,size2,displ,J1,J2 ); 
 
debug_disp('z-y force:') 
force_components(8,:)  =  forces_calc_z_y( size1,size2,displ,J1,J2 ); 
 
debug_disp('z-x force:') 
force_components(7,:)  =  forces_calc_z_x( size1,size2,displ,J1,J2 ); 
 
 
 
force_out  =  sum(force_components); 
end 
 
 
  
function stiffness_out  =  single_magnet_stiffness(displ) 
 
stiffness_components  =  repmat(NaN,[9 3]); 
 
  
debug_disp('  ') 
debug_disp('CALCULATING THINGS') 
debug_disp('==================') 
debug_disp('Displacement:') 
debug_disp(displ') 
debug_disp('Magnetisations:') 
debug_disp(J1') 
debug_disp(J2') 
 
 
  
size1_rot  =  swap_x_z(size1); 
size2_rot  =  swap_x_z(size2); 
d_rot   =  rotate_x_to_z(displ); 
J1_rot  =  rotate_x_to_z(J1); 
J2_rot  =  rotate_x_to_z(J2); 
 
debug_disp('x-z stiffness:') 
stiffness_components(3,:)  =  rotate_z_to_x( stiffnesses_calc_z_x( size1_rot,size2_rot,d_rot,J1_rot,J2_rot ) ); 
 
debug_disp('x-y stiffness:') 
stiffness_components(2,:)  =  rotate_z_to_x( stiffnesses_calc_z_y( size1_rot,size2_rot,d_rot,J1_rot,J2_rot ) ); 
 
debug_disp('x-x stiffness:') 
stiffness_components(1,:)  =  rotate_z_to_x( stiffnesses_calc_z_z( size1_rot,size2_rot,d_rot,J1_rot,J2_rot ) ); 
 
 
  
size1_rot  =  swap_y_z(size1); 
size2_rot  =  swap_y_z(size2); 
d_rot      =  rotate_y_to_z( displ ); 
J1_rot     =  rotate_y_to_z( J1    ); 
J2_rot     =  rotate_y_to_z( J2    ); 
 
debug_disp('y-z stiffness:') 
stiffness_components(6,:)  =  rotate_z_to_y( stiffnesses_calc_z_y( size1_rot,size2_rot,d_rot,J1_rot,J2_rot ) ); 
 
debug_disp('y-y stiffness:') 
stiffness_components(5,:)  =  rotate_z_to_y( stiffnesses_calc_z_z( size1_rot,size2_rot,d_rot,J1_rot,J2_rot ) ); 
 
debug_disp('y-x stiffness:') 
stiffness_components(4,:)  =  rotate_z_to_y( stiffnesses_calc_z_x( size1_rot,size2_rot,d_rot,J1_rot,J2_rot ) ); 
 
 
 
  
debug_disp('z-z stiffness:') 
stiffness_components(9,:)  =  stiffnesses_calc_z_z( size1,size2,displ,J1,J2 ); 
 
debug_disp('z-y stiffness:') 
stiffness_components(8,:)  =  stiffnesses_calc_z_y( size1,size2,displ,J1,J2 ); 
 
debug_disp('z-x stiffness:') 
stiffness_components(7,:)  =  stiffnesses_calc_z_x( size1,size2,displ,J1,J2 ); 
 
 
 
stiffness_out  =  sum(stiffness_components); 
end 
 
 
  
  
function calc_out  =  forces_calc_z_z(size1,size2,offset,J1,J2) 
 
J1  =  J1(3); 
J2  =  J2(3); 
 
  
if (J1==0 || J2==0) 
  debug_disp('Zero magnetisation.') 
  calc_out   =   [0; 0; 0]; 
  return; 
end 
 
u  =  offset(1) + size2(1) * (-1).^index_j - size1(1) * (-1).^index_i; 
v  =  offset(2) + size2(2) * (-1).^index_l - size1(2) * (-1).^index_k; 
w  =  offset(3) + size2(3) * (-1).^index_q - size1(3) * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
 
 
component_x  =   ... 
  + multiply_x_log_y( 0.5 * (v.^2-w.^2), r-u )  ... 
  + multiply_x_log_y( u.*v, r-v )  ... 
  + v.*w.*atan1(u.*v,r.*w)  ... 
  + 0.5 * r.*u; 
 
component_y  =   ... 
  + multiply_x_log_y( 0.5 * (u.^2-w.^2), r-v )  ... 
  + multiply_x_log_y( u.*v, r-u )  ... 
  + u.*w.*atan1(u.*v,r.*w) ... 
  + 0.5 * r.*v; 
 
component_z  =   ... 
  - multiply_x_log_y( u.*w, r-u )  ... 
  - multiply_x_log_y( v.*w, r-v )  ... 
  + u.*v.*atan1(u.*v,r.*w)  ... 
  - r.*w; 
 
  
component_x  =  index_sum.*component_x; 
component_y  =  index_sum.*component_y; 
component_z  =  index_sum.*component_z; 
 
calc_out  =  J1 * J2 * magconst .*  ... 
  [ sum(component_x(:)) ; 
    sum(component_y(:)) ; 
    sum(component_z(:)) ] ; 
 
debug_disp(calc_out') 
 
end 
 
 
 
 
 
 
 
 
  
function calc_out  =  forces_calc_z_y(size1,size2,offset,J1,J2) 
 
J1  =  J1(3); 
J2  =  J2(2); 
 
  
if (J1==0 || J2==0) 
  debug_disp('Zero magnetisation.') 
  calc_out   =   [0; 0; 0]; 
  return; 
end 
 
u  =  offset(1) + size2(1) * (-1).^index_j - size1(1) * (-1).^index_i; 
v  =  offset(2) + size2(2) * (-1).^index_l - size1(2) * (-1).^index_k; 
w  =  offset(3) + size2(3) * (-1).^index_q - size1(3) * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
 
 
component_x  =   ... 
  - multiply_x_log_y ( v .* w , r-u )  ... 
  + multiply_x_log_y ( v .* u , r+w )  ... 
  + multiply_x_log_y ( u .* w , r+v )  ... 
  - 0.5  *  u.^2 .* atan1( v .* w , u .* r )  ... 
  - 0.5  *  v.^2 .* atan1( u .* w , v .* r )  ... 
  - 0.5  *  w.^2 .* atan1( u .* v , w .* r ); 
 
component_y  =   ... 
  0.5  *  multiply_x_log_y( u.^2 - v.^2 , r+w )  ... 
  - multiply_x_log_y( u .* w , r-u )  ... 
  - u .* v .* atan1( u .* w , v .* r )  ... 
  - 0.5  *  w .* r; 
 
component_z  =   ... 
  0.5  *  multiply_x_log_y( u.^2 - w.^2 , r+v )  ... 
  - multiply_x_log_y( u .* v , r-u )  ... 
  - u .* w .* atan1( u .* v , w .* r )  ... 
  - 0.5  *  v .* r; 
 
allag_correction  =  -1; 
component_x  =  allag_correction * component_x; 
component_y  =  allag_correction * component_y; 
component_z  =  allag_correction * component_z; 
 
if 0 
  
S = u; 
T = v; 
U = w; 
R = r; 
 
component_x_ii  =   ... 
    ( 0.5 * atan1(U,S)+0.5 * atan1(T.*U,S.*R) ).*S.^2  ... 
    + T.*S - 3/2 * U.*S - multiply_x_log_y( S.*T , U+R )-T.^2 .* atan1(S,T)  ... 
    + U.* ( U.* (  ... 
       0.5 * atan1(S,U)+0.5 * atan1(S.*T,U.*R)  ... 
     )  ... 
   - multiply_x_log_y( T , S+R )+multiply_x_log_y(S,R-T)  ... 
   )  ... 
 + 0.5 * T.^2 .* atan1(S.*U,T.*R) ... 
; 
 
component_y_ii  =   ... 
 0.5 * U.*(R-2 * S)+ ... 
 multiply_x_log_y( 0.5 * (T.^2-S.^2) , U+R )+ ... 
 S.*T.*( atan1(U,T)+atan1(S.*U,T.*R) )+ ... 
 multiply_x_log_y( S.*U , R-S ) ... 
; 
 
component_z_ii  =   ... 
 0.5 * T.*(R-2 * S)+ ... 
 multiply_x_log_y( 0.5 * (U.^2-S.^2), T+R )+ ... 
 S.*U.*( atan1(T,U)+atan1(S.*T,U.*R) )+ ... 
 multiply_x_log_y( S.*T , R-S ) ... 
; 
 
if 1 
xx  =  index_sum.*component_x; 
xx_ii  =  index_sum.*component_x_ii; 
assert( abs(sum(xx(:)) - sum(xx_ii(:))) < 1e-8 ) 
end 
 
if 1 
yy  =  index_sum.*component_y; 
yy_ii  =  index_sum.*component_y_ii; 
assert( abs(sum(yy(:)) - sum(yy_ii(:))) < 1e-8 ) 
end 
 
if 1 
zz  =  index_sum.*component_z; 
zz_ii  =  index_sum.*component_z_ii; 
assert( abs(sum(zz(:)) - sum(zz_ii(:))) < 1e-8 ) 
end 
 
if 1 
component_x  =  component_x_ii; 
component_y  =  component_y_ii; 
component_z  =  component_z_ii; 
end 
 
 
end 
 
  
component_x  =  index_sum.*component_x; 
component_y  =  index_sum.*component_y; 
component_z  =  index_sum.*component_z; 
 
calc_out  =  J1 * J2 * magconst .*  ... 
  [ sum(component_x(:)) ; 
    sum(component_y(:)) ; 
    sum(component_z(:)) ] ; 
 
debug_disp(calc_out') 
 
end 
 
 
 
 
  
function calc_out  =  forces_calc_z_x(size1,size2,offset,J1,J2) 
 
forces_xyz  =  forces_calc_z_y( ... 
  abs(rotate_x_to_y(size1)), abs(rotate_x_to_y(size2)), rotate_x_to_y(offset), ... 
  J1, rotate_x_to_y(J2) ); 
 
calc_out  =  rotate_y_to_x( forces_xyz ); 
 
end 
 
 
 
  
function calc_out  =  stiffnesses_calc_z_z(size1,size2,offset,J1,J2) 
 
J1  =  J1(3); 
J2  =  J2(3); 
 
  
if (J1==0 || J2==0) 
  debug_disp('Zero magnetisation.') 
  calc_out   =   [0; 0; 0]; 
  return; 
end 
 
u  =  offset(1) + size2(1) * (-1).^index_j - size1(1) * (-1).^index_i; 
v  =  offset(2) + size2(2) * (-1).^index_l - size1(2) * (-1).^index_k; 
w  =  offset(3) + size2(3) * (-1).^index_q - size1(3) * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
 
 
component_x  =   ... 
  - r  ... 
  - (u.^2 .*v)./(u.^2+w.^2)  ... 
  - v.*log(r-v) ; 
 
component_y  =   ... 
  - r  ... 
  - (v.^2 .*u)./(v.^2+w.^2)  ... 
  - u.*log(r-u) ; 
 
component_z  =  - component_x - component_y; 
 
  
component_x  =  index_sum.*component_x; 
component_y  =  index_sum.*component_y; 
component_z  =  index_sum.*component_z; 
 
calc_out  =  J1 * J2 * magconst .*  ... 
  [ sum(component_x(:)) ; 
    sum(component_y(:)) ; 
    sum(component_z(:)) ] ; 
 
debug_disp(calc_out') 
 
end 
 
 
 
 
 
 
  
function calc_out  =  stiffnesses_calc_z_y(size1,size2,offset,J1,J2) 
 
J1  =  J1(3); 
J2  =  J2(2); 
 
  
if (J1==0 || J2==0) 
  debug_disp('Zero magnetisation.') 
  calc_out   =   [0; 0; 0]; 
  return; 
end 
 
u  =  offset(1) + size2(1) * (-1).^index_j - size1(1) * (-1).^index_i; 
v  =  offset(2) + size2(2) * (-1).^index_l - size1(2) * (-1).^index_k; 
w  =  offset(3) + size2(3) * (-1).^index_q - size1(3) * (-1).^index_p; 
r  =  sqrt(u.^2+v.^2+w.^2); 
 
 
 
component_x  =   -((u.^2 .*v)./(u.^2 + v.^2)) - (u.^2 .*w)./(u.^2 + w.^2)  ... 
     + u.*atan1(v.*w,r.*u) - multiply_x_log_y( w , r + v ) +  ... 
     - multiply_x_log_y( v , r + w ); 
 
component_y  =   v/2 - (u.^2 .*v)./(u.^2 + v.^2) + (u.*v.*w)./(v.^2 + w.^2)  ... 
     +  u.*atan1(u.*w,r.*v) + multiply_x_log_y( v , r + w ); 
 
component_z  =  - component_x - component_y; 
 
allag_correction  =  -1; 
component_x  =  allag_correction * component_x; 
component_y  =  allag_correction * component_y; 
component_z  =  allag_correction * component_z; 
 
  
component_x  =  index_sum.*component_x; 
component_y  =  index_sum.*component_y; 
component_z  =  index_sum.*component_z; 
 
calc_out  =  J1 * J2 * magconst .*  ... 
  [ sum(component_x(:)) ; 
    sum(component_y(:)) ; 
    sum(component_z(:)) ] ; 
 
debug_disp(calc_out') 
 
end 
 
 
 
 
 
  
function calc_out  =  stiffnesses_calc_z_x(size1,size2,offset,J1,J2) 
 
stiffnesses_xyz  =  stiffnesses_calc_z_y( ... 
  abs(rotate_x_to_y(size1)), abs(rotate_x_to_y(size2)), rotate_x_to_y(offset), ... 
  J1, rotate_x_to_y(J2) ); 
 
calc_out  =  rotate_y_to_x(stiffnesses_xyz); 
 
end 
 
 
 
 
  
function out  =  multiply_x_log_y(x,y) 
  out  =  x.*log(y); 
  out(~isfinite(out)) = 0; 
end 
 
  
function out  =  atan1(x,y) 
  out  =  zeros(size(x)); 
  ind  =  x~=0 & y~=0; 
  out(ind)  =  atan(x(ind)./y(ind)); 
end 
 
 
 
 
 
 
 
end 
 

