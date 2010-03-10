 
 
function [varargout]  =  magnetforces(magnet_fixed, magnet_float, displ, varargin) 
 
  
%% MAGNETFORCES  Calculate forces between two cuboid magnets 
% 
% Finish this off later. Please read the PDF documentation instead for now. 
% 
 
 
 
 
 
  
debug_disp  =  @(str) disp([]); 
calc_force_bool  =  false; 
calc_stiffness_bool  =  false; 
 
% Undefined calculation flags for the three directions: 
calc_xyz  =  [-1 -1 -1]; 
 
for ii  =  1:length(varargin) 
  switch varargin{ii} 
    case 'debug',      debug_disp  =  @(str) disp(str); 
    case 'force',      calc_force_bool      =  true; 
    case 'stiffness',  calc_stiffness_bool  =  true; 
    case 'x',  calc_xyz(1)  =  1; 
    case 'y',  calc_xyz(2)  =  1; 
    case 'z',  calc_xyz(3)  =  1; 
    otherwise 
      error(['Unknown calculation option ''',varargin{ii},'''']) 
  end 
end 
 
% If none of |'x'|, |'y'|, |'z'| are specified, calculate all. 
if all( calc_xyz == -1 ) 
  calc_xyz  =  [1 1 1]; 
end 
 
calc_xyz( calc_xyz == -1 )  =  0; 
 
if ~calc_force_bool && ~calc_stiffness_bool 
  calc_force_bool  =  true; 
end 
 
 
  
if size(displ,1) == 3 
  % all good 
elseif size(displ,2) == 3 
  displ  =  transpose(displ); 
else 
  error(['Displacements matrix should be of size (3, D)', ... 
         'where D is the number of displacements.']) 
end 
 
Ndispl  =  size(displ,2); 
 
if calc_force_bool 
  forces_out  =  repmat(NaN,[3 Ndispl]); 
end 
 
if calc_stiffness_bool 
  stiffnesses_out  =  repmat(NaN,[3 Ndispl]); 
end 
 
 
  
size1  =  reshape(magnet_fixed.dim/2,[3 1]); 
size2  =  reshape(magnet_float.dim/2,[3 1]); 
 
J1  =  resolve_magnetisations(magnet_fixed.magn,magnet_fixed.magdir); 
J2  =  resolve_magnetisations(magnet_float.magn,magnet_float.magdir); 
 
  
magconst  =  1/(4 * pi * (4 * pi * 1e-7)); 
 
[index_i, index_j, index_k, index_l, index_p, index_q]  =  ndgrid([0 1]); 
 
index_sum  =  (-1).^(index_i+index_j+index_k+index_l+index_p+index_q); 
 
 
  
swap_x_y  =  @(vec) vec([2 1 3]); 
swap_x_z  =  @(vec) vec([3 2 1]); 
swap_y_z  =  @(vec) vec([1 3 2]); 
 
rotate_z_to_x  =  @(vec) [  vec(3);  vec(2); -vec(1) ] ; % Ry( 90) 
rotate_x_to_z  =  @(vec) [ -vec(3);  vec(2);  vec(1) ] ; % Ry(-90) 
 
rotate_y_to_z  =  @(vec) [  vec(1); -vec(3);  vec(2) ] ; % Rx( 90) 
rotate_z_to_y  =  @(vec) [  vec(1);  vec(3); -vec(2) ] ; % Rx(-90) 
 
rotate_x_to_y  =  @(vec) [ -vec(2);  vec(1);  vec(3) ] ; % Rz( 90) 
rotate_y_to_x  =  @(vec) [  vec(2); -vec(1);  vec(3) ] ; % Rz(-90) 
 
size1_x  =  swap_x_z(size1); 
size2_x  =  swap_x_z(size2); 
J1_x     =  rotate_x_to_z(J1); 
J2_x     =  rotate_x_to_z(J2); 
 
size1_y  =  swap_y_z(size1); 
size2_y  =  swap_y_z(size2); 
J1_y     =  rotate_y_to_z(J1); 
J2_y     =  rotate_y_to_z(J2); 
 
 
 
  
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
for ii  =  1:length(varargin) 
  switch varargin{ii} 
    case 'force' 
      varargout{ii}  =  forces_out; 
    case 'stiffness' 
      varargout{ii}  =  stiffnesses_out; 
  end 
end 
 
 
 
 
  
function J  =  resolve_magnetisations(magn,magdir) 
 
if length(magdir)==2 
  J_r  =  magn; 
  J_t  =  magdir(1); 
  J_p  =  magdir(2); 
  J    =  [ J_r  *  cosd(J_p)  *  cosd(J_t)  ;  ... 
          J_r  *  cosd(J_p)  *  sind(J_t)  ;  ... 
          J_r  *  sind(J_p) ]; 
else 
  if all(magdir == zeros(size(magdir)) ) 
    J  =  [0; 0; 0]; 
  else 
    J  =  magn * magdir/norm(magdir); 
    J  =  reshape(J,[3 1]); 
  end 
end 
 
end 
 
 
  
function force_out  =  single_magnet_force(displ) 
 
force_components  =  repmat(NaN,[9 3]); 
 
  
d_x   =  rotate_x_to_z(displ); 
d_y   =  rotate_y_to_z(displ); 
 
 
  
debug_disp('  ') 
debug_disp('CALCULATING THINGS') 
debug_disp('==================') 
debug_disp('Displacement:') 
debug_disp(displ') 
debug_disp('Magnetisations:') 
debug_disp(J1') 
debug_disp(J2') 
 
 
  
calc_xyz  =  swap_x_z(calc_xyz); 
 
debug_disp('Forces x-x:') 
force_components(1,:)  =   ... 
  rotate_z_to_x( forces_calc_z_z(size1_x,size2_x,d_x,J1_x,J2_x) ); 
 
debug_disp('Forces x-y:') 
force_components(2,:)  =   ... 
  rotate_z_to_x( forces_calc_z_y(size1_x,size2_x,d_x,J1_x,J2_x) ); 
 
debug_disp('Forces x-z:') 
force_components(3,:)  =   ... 
  rotate_z_to_x( forces_calc_z_x(size1_x,size2_x,d_x,J1_x,J2_x) ); 
 
calc_xyz  =  swap_x_z(calc_xyz); 
 
 
  
calc_xyz  =  swap_y_z(calc_xyz); 
 
debug_disp('Forces y-x:') 
force_components(4,:)  =   ... 
  rotate_z_to_y( forces_calc_z_x(size1_y,size2_y,d_y,J1_y,J2_y) ); 
 
debug_disp('Forces y-y:') 
force_components(5,:)  =   ... 
  rotate_z_to_y( forces_calc_z_z(size1_y,size2_y,d_y,J1_y,J2_y) ); 
 
debug_disp('Forces y-z:') 
force_components(6,:)  =   ... 
  rotate_z_to_y( forces_calc_z_y(size1_y,size2_y,d_y,J1_y,J2_y) ); 
 
calc_xyz  =  swap_y_z(calc_xyz); 
 
 
  
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
 
  
d_x   =  rotate_x_to_z(displ); 
d_y   =  rotate_y_to_z(displ); 
 
 
  
debug_disp('  ') 
debug_disp('CALCULATING THINGS') 
debug_disp('==================') 
debug_disp('Displacement:') 
debug_disp(displ') 
debug_disp('Magnetisations:') 
debug_disp(J1') 
debug_disp(J2') 
 
 
  
debug_disp('z-x stiffness:') 
stiffness_components(7,:)  =   ... 
  stiffnesses_calc_z_x( size1,size2,displ,J1,J2 ); 
 
debug_disp('z-y stiffness:') 
stiffness_components(8,:)  =   ... 
  stiffnesses_calc_z_y( size1,size2,displ,J1,J2 ); 
 
debug_disp('z-z stiffness:') 
stiffness_components(9,:)  =   ... 
  stiffnesses_calc_z_z( size1,size2,displ,J1,J2 ); 
 
calc_xyz  =  swap_x_z(calc_xyz); 
 
debug_disp('x-x stiffness:') 
stiffness_components(1,:)  =   ... 
  swap_x_z( stiffnesses_calc_z_z( size1_x,size2_x,d_x,J1_x,J2_x ) ); 
 
debug_disp('x-y stiffness:') 
stiffness_components(2,:)  =   ... 
  swap_x_z( stiffnesses_calc_z_y( size1_x,size2_x,d_x,J1_x,J2_x ) ); 
 
debug_disp('x-z stiffness:') 
stiffness_components(3,:)  =   ... 
  swap_x_z( stiffnesses_calc_z_x( size1_x,size2_x,d_x,J1_x,J2_x ) ); 
 
calc_xyz  =  swap_x_z(calc_xyz); 
 
calc_xyz  =  swap_y_z(calc_xyz); 
 
debug_disp('y-x stiffness:') 
stiffness_components(4,:)  =   ... 
  swap_y_z( stiffnesses_calc_z_x( size1_y,size2_y,d_y,J1_y,J2_y ) ); 
 
debug_disp('y-y stiffness:') 
stiffness_components(5,:)  =   ... 
  swap_y_z( stiffnesses_calc_z_z( size1_y,size2_y,d_y,J1_y,J2_y ) ); 
 
debug_disp('y-z stiffness:') 
stiffness_components(6,:)  =   ... 
  swap_y_z( stiffnesses_calc_z_y( size1_y,size2_y,d_y,J1_y,J2_y ) ); 
 
calc_xyz  =  swap_y_z(calc_xyz); 
 
 
 
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
 
 
 
if calc_xyz(1) 
  component_x  =   ... 
    + multiply_x_log_y( 0.5 * (v.^2-w.^2), r-u )  ... 
    + multiply_x_log_y( u.*v, r-v )  ... 
    + v.*w.*atan1(u.*v,r.*w)  ... 
    + 0.5 * r.*u; 
end 
 
if calc_xyz(2) 
  component_y  =   ... 
    + multiply_x_log_y( 0.5 * (u.^2-w.^2), r-v )  ... 
    + multiply_x_log_y( u.*v, r-u )  ... 
    + u.*w.*atan1(u.*v,r.*w)  ... 
    + 0.5 * r.*v; 
end 
 
if calc_xyz(3) 
  component_z  =   ... 
    - multiply_x_log_y( u.*w, r-u )  ... 
    - multiply_x_log_y( v.*w, r-v )  ... 
    + u.*v.*atan1(u.*v,r.*w)  ... 
    - r.*w; 
end 
 
  
if calc_xyz(1) 
  component_x  =  index_sum.*component_x; 
else 
  component_x  =  0; 
end 
 
if calc_xyz(2) 
  component_y  =  index_sum.*component_y; 
else 
  component_y  =  0; 
end 
 
if calc_xyz(3) 
  component_z  =  index_sum.*component_z; 
else 
  component_z  =  0; 
end 
 
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
 
 
 
allag_correction  =  -1; 
 
if calc_xyz(1) 
  component_x  =   ... 
    - multiply_x_log_y ( v .* w , r-u )  ... 
    + multiply_x_log_y ( v .* u , r+w )  ... 
    + multiply_x_log_y ( u .* w , r+v )  ... 
    - 0.5  *  u.^2 .* atan1( v .* w , u .* r )  ... 
    - 0.5  *  v.^2 .* atan1( u .* w , v .* r )  ... 
    - 0.5  *  w.^2 .* atan1( u .* v , w .* r ); 
  component_x  =  allag_correction * component_x; 
end 
 
if calc_xyz(2) 
  component_y  =   ... 
    0.5  *  multiply_x_log_y( u.^2 - v.^2 , r+w )  ... 
    - multiply_x_log_y( u .* w , r-u )  ... 
    - u .* v .* atan1( u .* w , v .* r )  ... 
    - 0.5  *  w .* r; 
  component_y  =  allag_correction * component_y; 
end 
 
if calc_xyz(3) 
  component_z  =   ... 
    0.5  *  multiply_x_log_y( u.^2 - w.^2 , r+v )  ... 
    - multiply_x_log_y( u .* v , r-u )  ... 
    - u .* w .* atan1( u .* v , w .* r )  ... 
    - 0.5  *  v .* r; 
  component_z  =  allag_correction * component_z; 
end 
 
  
if calc_xyz(1) 
  component_x  =  index_sum.*component_x; 
else 
  component_x  =  0; 
end 
 
if calc_xyz(2) 
  component_y  =  index_sum.*component_y; 
else 
  component_y  =  0; 
end 
 
if calc_xyz(3) 
  component_z  =  index_sum.*component_z; 
else 
  component_z  =  0; 
end 
 
calc_out  =  J1 * J2 * magconst .*  ... 
  [ sum(component_x(:)) ; 
    sum(component_y(:)) ; 
    sum(component_z(:)) ] ; 
 
debug_disp(calc_out') 
 
end 
 
 
 
 
  
function calc_out  =  forces_calc_z_x(size1,size2,offset,J1,J2) 
 
calc_xyz  =  swap_x_y(calc_xyz); 
 
forces_xyz  =  forces_calc_z_y( ... 
  swap_x_y(size1), swap_x_y(size2), rotate_x_to_y(offset), ... 
  J1, rotate_x_to_y(J2) ); 
 
calc_xyz  =  swap_x_y(calc_xyz); 
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
 
 
 
if calc_xyz(1) || calc_xyz(3) 
  component_x  =  - r - (u.^2 .*v)./(u.^2+w.^2) - v.*log(r-v) ; 
end 
 
if calc_xyz(2) || calc_xyz(3) 
  component_y  =  - r - (v.^2 .*u)./(v.^2+w.^2) - u.*log(r-u) ; 
end 
 
if calc_xyz(3) 
  component_z  =  - component_x - component_y; 
end 
 
  
if calc_xyz(1) 
  component_x  =  index_sum.*component_x; 
else 
  component_x  =  0; 
end 
 
if calc_xyz(2) 
  component_y  =  index_sum.*component_y; 
else 
  component_y  =  0; 
end 
 
if calc_xyz(3) 
  component_z  =  index_sum.*component_z; 
else 
  component_z  =  0; 
end 
 
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
 
 
 
if calc_xyz(1) || calc_xyz(3) 
  component_x  =   ((u.^2 .*v)./(u.^2 + v.^2)) + (u.^2 .*w)./(u.^2 + w.^2)  ... 
       - u.*atan1(v.*w,r.*u) + multiply_x_log_y( w , r + v ) +  ... 
       + multiply_x_log_y( v , r + w ); 
end 
 
if calc_xyz(2) || calc_xyz(3) 
  component_y  =  - v/2 + (u.^2 .*v)./(u.^2 + v.^2) - (u.*v.*w)./(v.^2 + w.^2)  ... 
       -  u.*atan1(u.*w,r.*v) - multiply_x_log_y( v , r + w ); 
end 
 
if calc_xyz(3) 
  component_z  =  - component_x - component_y; 
end 
 
  
if calc_xyz(1) 
  component_x  =  index_sum.*component_x; 
else 
  component_x  =  0; 
end 
 
if calc_xyz(2) 
  component_y  =  index_sum.*component_y; 
else 
  component_y  =  0; 
end 
 
if calc_xyz(3) 
  component_z  =  index_sum.*component_z; 
else 
  component_z  =  0; 
end 
 
calc_out  =  J1 * J2 * magconst .*  ... 
  [ sum(component_x(:)) ; 
    sum(component_y(:)) ; 
    sum(component_z(:)) ] ; 
 
debug_disp(calc_out') 
 
end 
 
 
 
 
 
  
function calc_out  =  stiffnesses_calc_z_x(size1,size2,offset,J1,J2) 
 
calc_xyz  =  swap_x_y(calc_xyz); 
 
stiffnesses_xyz  =  stiffnesses_calc_z_y( ... 
  swap_x_y(size1), swap_x_y(size2), rotate_x_to_y(offset), ... 
  J1, rotate_x_to_y(J2) ); 
 
calc_xyz  =  swap_x_y(calc_xyz); 
calc_out  =  swap_x_y(stiffnesses_xyz); 
 
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
 

