
function [varargout] = magnetforces(magnet_fixed, magnet_float, displ, varargin)


%% MAGNETFORCES  Calculate forces between two cuboid magnets
%
% Finish this off later. Please read the PDF documentation instead for now.
%





debug_disp = @(str) disp([]);
calc_force_bool = false;
calc_stiffness_bool = false;
calc_torque_bool = false;

% Undefined calculation flags for the three directions:
calc_xyz = [-1 -1 -1];

for ii = 1:length(varargin)
  switch varargin{ii}
    case 'debug',      debug_disp = @(str) disp(str);
    case 'force',      calc_force_bool     = true;
    case 'stiffness',  calc_stiffness_bool = true;
    case 'torque',     calc_torque_bool    = true;
    case 'x',  calc_xyz(1) = 1;
    case 'y',  calc_xyz(2) = 1;
    case 'z',  calc_xyz(3) = 1;
    otherwise
      error(['Unknown calculation option ''',varargin{ii},''''])
  end
end

% If none of |'x'|, |'y'|, |'z'| are specified, calculate all.
if all( calc_xyz == -1 )
  calc_xyz = [1 1 1];
end

calc_xyz( calc_xyz == -1 ) = 0;

if ~calc_force_bool && ~calc_stiffness_bool && ~calc_torque_bool
  calc_force_bool = true;
end


if size(displ,1) == 3
  % all good
elseif size(displ,2) == 3
  displ = transpose(displ);
else
  error(['Displacements matrix should be of size (3, D) ',...
         'where D is the number of displacements.'])
end

Ndispl = size(displ,2);

if calc_force_bool
  forces_out = repmat(NaN,[3 Ndispl]);
end

if calc_stiffness_bool
  stiffnesses_out = repmat(NaN,[3 Ndispl]);
end

if calc_torque_bool
  torques_out = repmat(NaN,[3 Ndispl]);
end


if ~isfield(magnet_fixed,'type')
  if length(magnet_fixed.dim) == 2
    magnet_fixed.type = 'cylinder';
  else
    magnet_fixed.type = 'cuboid';
  end
end

if ~isfield(magnet_float,'type')
  if length(magnet_float.dim) == 2
    magnet_float.type = 'cylinder';
  else
    magnet_float.type = 'cuboid';
  end
end

if ~strcmp(magnet_fixed.type, magnet_float.type)
  error('Magnets must be of same type')
end

magtype = magnet_fixed.type;

if strcmp(magtype,'cuboid')

  size1 = reshape(magnet_fixed.dim/2,[3 1]);
  size2 = reshape(magnet_float.dim/2,[3 1]);

  J1 = resolve_magnetisations(magnet_fixed.magn,magnet_fixed.magdir);
  J2 = resolve_magnetisations(magnet_float.magn,magnet_float.magdir);

  if calc_torque_bool
    if ~isfield(magnet_float,'lever')
      magnet_float.lever = [0 0 0];
    end
  end

elseif strcmp(magtype,'cylinder')

  size1 = reshape(magnet_fixed.dim,[2 1]);
  size2 = reshape(magnet_float.dim,[2 1]);

  if ~isfield(magnet_fixed,'dir')
    magnet_fixed.dir = [0 0 1];
  end
  if ~isfield(magnet_float,'dir')
    magnet_float.dir = [0 0 1];
  end
  if abs(magnet_fixed.dir) ~= abs(magnet_float.dir)
    error('Cylindrical magnets must be oriented in the same direction')
  end

  if ~isfield(magnet_fixed,'magdir')
    magnet_fixed.magdir = [0 0 1];
  end
  if abs(magnet_fixed.dir) ~= abs(magnet_fixed.magdir)
    error('Cylindrical magnets must be magnetised in the same direction as their orientation')
  end

  if ~isfield(magnet_float,'magdir')
    magnet_float.magdir = [0 0 1];
  end
  if abs(magnet_float.dir) ~= abs(magnet_float.magdir)
    error('Cylindrical magnets must be magnetised in the same direction as their orientation')
  end

  cyldir = find(magnet_float.magdir ~= 0);
  cylnotdir = find(magnet_float.magdir == 0);
  if length(cyldir) ~= 1
    error('Cylindrical magnets must be aligned in one of the x, y or z directions')
  end

  magnet_float.magdir = magnet_float.magdir(:);
  magnet_fixed.magdir = magnet_fixed.magdir(:);
  magnet_float.dir = magnet_float.dir(:);
  magnet_fixed.dir = magnet_fixed.dir(:);

  if ~isfield(magnet_fixed,'magn')
    magnet_fixed.magn = 4*pi*1e-7*magnet_fixed.turns*magnet_fixed.current/magnet_fixed.dim(2);
  end
  if ~isfield(magnet_float,'magn')
    magnet_float.magn = 4*pi*1e-7*magnet_float.turns*magnet_float.current/magnet_float.dim(2);
  end

  J1 = magnet_fixed.magn*magnet_fixed.magdir;
  J2 = magnet_float.magn*magnet_float.magdir;

end


magconst = 1/(4*pi*(4*pi*1e-7));

[index_i, index_j, index_k, index_l, index_p, index_q] = ndgrid([0 1]);

index_sum = (-1).^(index_i+index_j+index_k+index_l+index_p+index_q);


if strcmp(magtype,'cuboid')

swap_x_y = @(vec) vec([2 1 3]);
swap_x_z = @(vec) vec([3 2 1]);
swap_y_z = @(vec) vec([1 3 2]);

rotate_z_to_x = @(vec) [  vec(3);  vec(2); -vec(1) ] ; % Ry( 90)
rotate_x_to_z = @(vec) [ -vec(3);  vec(2);  vec(1) ] ; % Ry(-90)

rotate_y_to_z = @(vec) [  vec(1); -vec(3);  vec(2) ] ; % Rx( 90)
rotate_z_to_y = @(vec) [  vec(1);  vec(3); -vec(2) ] ; % Rx(-90)

rotate_x_to_y = @(vec) [ -vec(2);  vec(1);  vec(3) ] ; % Rz( 90)
rotate_y_to_x = @(vec) [  vec(2); -vec(1);  vec(3) ] ; % Rz(-90)

size1_x = swap_x_z(size1);
size2_x = swap_x_z(size2);
J1_x    = rotate_x_to_z(J1);
J2_x    = rotate_x_to_z(J2);

size1_y = swap_y_z(size1);
size2_y = swap_y_z(size2);
J1_y    = rotate_y_to_z(J1);
J2_y    = rotate_y_to_z(J2);

end



if strcmp(magtype,'cuboid')

if calc_force_bool
  for ii = 1:Ndispl
    forces_out(:,ii)  =  single_magnet_force(displ(:,ii));
  end
end

if calc_stiffness_bool
  for ii = 1:Ndispl
    stiffnesses_out(:,ii)  =  single_magnet_stiffness(displ(:,ii));
  end
end

if calc_torque_bool
  torques_out  =  single_magnet_torque(displ,magnet_float.lever);
end

elseif strcmp(magtype,'cylinder')

if strcmp(magtype,'cylinder')
  if any(displ(cylnotdir,:)~=0)
    error(['Displacements for cylindrical magnets may only be axial. ',...
           'I.e., only in the direction of their alignment.'])
  end
end

if calc_force_bool
  forces_out = magnet_fixed.dir*...
    forces_cyl_calc(size1, size2, squeeze(displ(cyldir,:)), J1(cyldir), J2(cyldir));
end

if calc_stiffness_bool
  error('Stiffness cannot be calculated for cylindrical magnets yet.')
end

if calc_torque_bool
  error('Torques cannot be calculated for cylindrical magnets yet.')
end


end


ii = 0;
if calc_force_bool
  ii = ii + 1;
  varargout{ii} = forces_out;
end

if calc_stiffness_bool
  ii = ii + 1;
  varargout{ii} = stiffnesses_out;
end

if calc_torque_bool
  ii = ii + 1;
  varargout{ii} = torques_out;
end




function J = resolve_magnetisations(magn,magdir)

if length(magdir)==2
  J_r = magn;
  J_t = magdir(1);
  J_p = magdir(2);
  J   = [ J_r * cosd(J_p) * cosd(J_t)  ; ...
          J_r * cosd(J_p) * sind(J_t)  ; ...
          J_r * sind(J_p) ];
else
  if all(magdir == zeros(size(magdir)) )
    J = [0; 0; 0];
  else
    J = magn*magdir/norm(magdir);
    J = reshape(J,[3 1]);
  end
end

end


function force_out = single_magnet_force(displ)

force_components = repmat(NaN,[9 3]);


d_x  = rotate_x_to_z(displ);
d_y  = rotate_y_to_z(displ);

if calc_torque_bool
  l_x = rotate_x_to_z(lever);
  l_y = rotate_y_to_z(lever);
end


debug_disp('  ')
debug_disp('CALCULATING THINGS')
debug_disp('==================')
debug_disp('Displacement:')
debug_disp(displ')
debug_disp('Magnetisations:')
debug_disp(J1')
debug_disp(J2')


calc_xyz = swap_x_z(calc_xyz);

debug_disp('Forces x-x:')
force_components(1,:) = ...
  rotate_z_to_x( forces_calc_z_z(size1_x,size2_x,d_x,J1_x,J2_x) );

debug_disp('Forces x-y:')
force_components(2,:) = ...
  rotate_z_to_x( forces_calc_z_y(size1_x,size2_x,d_x,J1_x,J2_x) );

debug_disp('Forces x-z:')
force_components(3,:) = ...
  rotate_z_to_x( forces_calc_z_x(size1_x,size2_x,d_x,J1_x,J2_x) );

calc_xyz = swap_x_z(calc_xyz);


calc_xyz = swap_y_z(calc_xyz);

debug_disp('Forces y-x:')
force_components(4,:) = ...
  rotate_z_to_y( forces_calc_z_x(size1_y,size2_y,d_y,J1_y,J2_y) );

debug_disp('Forces y-y:')
force_components(5,:) = ...
  rotate_z_to_y( forces_calc_z_z(size1_y,size2_y,d_y,J1_y,J2_y) );

debug_disp('Forces y-z:')
force_components(6,:) = ...
  rotate_z_to_y( forces_calc_z_y(size1_y,size2_y,d_y,J1_y,J2_y) );

calc_xyz = swap_y_z(calc_xyz);


debug_disp('z-z force:')
force_components(9,:) = forces_calc_z_z( size1,size2,displ,J1,J2 );

debug_disp('z-y force:')
force_components(8,:) = forces_calc_z_y( size1,size2,displ,J1,J2 );

debug_disp('z-x force:')
force_components(7,:) = forces_calc_z_x( size1,size2,displ,J1,J2 );


force_out = sum(force_components);
end


function torques_out = single_magnet_torque(displ,lever)

force_components = repmat(NaN,[9 3]);


d_x  = rotate_x_to_z(displ);
d_y  = rotate_y_to_z(displ);

if calc_torque_bool
  l_x = rotate_x_to_z(lever);
  l_y = rotate_y_to_z(lever);
end


debug_disp('  ')
debug_disp('CALCULATING THINGS')
debug_disp('==================')
debug_disp('Displacement:')
debug_disp(displ')
debug_disp('Magnetisations:')
debug_disp(J1')
debug_disp(J2')


debug_disp('Torque: z-z:')
torque_components(9,:) = torques_calc_z_z( size1,size2,displ,lever,J1,J2 );

debug_disp('Torque z-y:')
torque_components(8,:) = torques_calc_z_y( size1,size2,displ,lever,J1,J2 );

debug_disp('Torque z-x:')
torque_components(7,:) = torques_calc_z_x( size1,size2,displ,lever,J1,J2 );

calc_xyz = swap_x_z(calc_xyz);

debug_disp('Torques x-x:')
torque_components(1,:) = ...
  rotate_z_to_x( torques_calc_z_z(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );

debug_disp('Torques x-y:')
torque_components(2,:) = ...
  rotate_z_to_x( torques_calc_z_y(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );

debug_disp('Torques x-z:')
torque_components(3,:) = ...
  rotate_z_to_x( torques_calc_z_x(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );

calc_xyz = swap_x_z(calc_xyz);

calc_xyz = swap_y_z(calc_xyz);

debug_disp('Torques y-x:')
torque_components(4,:) = ...
  rotate_z_to_y( torques_calc_z_x(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );

debug_disp('Torques y-y:')
torque_components(5,:) = ...
  rotate_z_to_y( torques_calc_z_z(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );

debug_disp('Torques y-z:')
torque_components(6,:) = ...
  rotate_z_to_y( torques_calc_z_y(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );

calc_xyz = swap_y_z(calc_xyz);



torques_out = sum(torque_components);
end




function stiffness_out = single_magnet_stiffness(displ)

stiffness_components = repmat(NaN,[9 3]);


d_x  = rotate_x_to_z(displ);
d_y  = rotate_y_to_z(displ);

if calc_torque_bool
  l_x = rotate_x_to_z(lever);
  l_y = rotate_y_to_z(lever);
end


debug_disp('  ')
debug_disp('CALCULATING THINGS')
debug_disp('==================')
debug_disp('Displacement:')
debug_disp(displ')
debug_disp('Magnetisations:')
debug_disp(J1')
debug_disp(J2')


debug_disp('z-x stiffness:')
stiffness_components(7,:) = ...
  stiffnesses_calc_z_x( size1,size2,displ,J1,J2 );

debug_disp('z-y stiffness:')
stiffness_components(8,:) = ...
  stiffnesses_calc_z_y( size1,size2,displ,J1,J2 );

debug_disp('z-z stiffness:')
stiffness_components(9,:) = ...
  stiffnesses_calc_z_z( size1,size2,displ,J1,J2 );

calc_xyz = swap_x_z(calc_xyz);

debug_disp('x-x stiffness:')
stiffness_components(1,:) = ...
  swap_x_z( stiffnesses_calc_z_z( size1_x,size2_x,d_x,J1_x,J2_x ) );

debug_disp('x-y stiffness:')
stiffness_components(2,:) = ...
  swap_x_z( stiffnesses_calc_z_y( size1_x,size2_x,d_x,J1_x,J2_x ) );

debug_disp('x-z stiffness:')
stiffness_components(3,:) = ...
  swap_x_z( stiffnesses_calc_z_x( size1_x,size2_x,d_x,J1_x,J2_x ) );

calc_xyz = swap_x_z(calc_xyz);

calc_xyz = swap_y_z(calc_xyz);

debug_disp('y-x stiffness:')
stiffness_components(4,:) = ...
  swap_y_z( stiffnesses_calc_z_x( size1_y,size2_y,d_y,J1_y,J2_y ) );

debug_disp('y-y stiffness:')
stiffness_components(5,:) = ...
  swap_y_z( stiffnesses_calc_z_z( size1_y,size2_y,d_y,J1_y,J2_y ) );

debug_disp('y-z stiffness:')
stiffness_components(6,:) = ...
  swap_y_z( stiffnesses_calc_z_y( size1_y,size2_y,d_y,J1_y,J2_y ) );

calc_xyz = swap_y_z(calc_xyz);




stiffness_out = sum(stiffness_components);
end



function calc_out = forces_calc_z_z(size1,size2,offset,J1,J2)

J1 = J1(3);
J2 = J2(3);


if (J1==0 || J2==0)
  debug_disp('Zero magnetisation.')
  calc_out  =  [0; 0; 0];
  return;
end

u = offset(1) + size2(1)*(-1).^index_j - size1(1)*(-1).^index_i;
v = offset(2) + size2(2)*(-1).^index_l - size1(2)*(-1).^index_k;
w = offset(3) + size2(3)*(-1).^index_q - size1(3)*(-1).^index_p;
r = sqrt(u.^2+v.^2+w.^2);


if calc_xyz(1)
  component_x = ...
    + multiply_x_log_y( 0.5*(v.^2-w.^2), r-u ) ...
    + multiply_x_log_y( u.*v, r-v ) ...
    + v.*w.*atan1(u.*v,r.*w) ...
    + 0.5*r.*u;
end

if calc_xyz(2)
  component_y = ...
    + multiply_x_log_y( 0.5*(u.^2-w.^2), r-v ) ...
    + multiply_x_log_y( u.*v, r-u ) ...
    + u.*w.*atan1(u.*v,r.*w) ...
    + 0.5*r.*v;
end

if calc_xyz(3)
  component_z = ...
    - multiply_x_log_y( u.*w, r-u ) ...
    - multiply_x_log_y( v.*w, r-v ) ...
    + u.*v.*atan1(u.*v,r.*w) ...
    - r.*w;
end


if calc_xyz(1)
  component_x = index_sum.*component_x;
else
  component_x = 0;
end

if calc_xyz(2)
  component_y = index_sum.*component_y;
else
  component_y = 0;
end

if calc_xyz(3)
  component_z = index_sum.*component_z;
else
  component_z = 0;
end

calc_out = J1*J2*magconst .* ...
  [ sum(component_x(:)) ;
    sum(component_y(:)) ;
    sum(component_z(:)) ] ;

debug_disp(calc_out')

end







function calc_out = forces_calc_z_y(size1,size2,offset,J1,J2)

J1 = J1(3);
J2 = J2(2);


if (J1==0 || J2==0)
  debug_disp('Zero magnetisation.')
  calc_out  =  [0; 0; 0];
  return;
end

u = offset(1) + size2(1)*(-1).^index_j - size1(1)*(-1).^index_i;
v = offset(2) + size2(2)*(-1).^index_l - size1(2)*(-1).^index_k;
w = offset(3) + size2(3)*(-1).^index_q - size1(3)*(-1).^index_p;
r = sqrt(u.^2+v.^2+w.^2);


allag_correction = -1;

if calc_xyz(1)
  component_x = ...
    - multiply_x_log_y ( v .* w , r-u ) ...
    + multiply_x_log_y ( v .* u , r+w ) ...
    + multiply_x_log_y ( u .* w , r+v ) ...
    - 0.5 * u.^2 .* atan1( v .* w , u .* r ) ...
    - 0.5 * v.^2 .* atan1( u .* w , v .* r ) ...
    - 0.5 * w.^2 .* atan1( u .* v , w .* r );
  component_x = allag_correction*component_x;
end

if calc_xyz(2)
  component_y = ...
    0.5 * multiply_x_log_y( u.^2 - v.^2 , r+w ) ...
    - multiply_x_log_y( u .* w , r-u ) ...
    - u .* v .* atan1( u .* w , v .* r ) ...
    - 0.5 * w .* r;
  component_y = allag_correction*component_y;
end

if calc_xyz(3)
  component_z = ...
    0.5 * multiply_x_log_y( u.^2 - w.^2 , r+v ) ...
    - multiply_x_log_y( u .* v , r-u ) ...
    - u .* w .* atan1( u .* v , w .* r ) ...
    - 0.5 * v .* r;
  component_z = allag_correction*component_z;
end


if calc_xyz(1)
  component_x = index_sum.*component_x;
else
  component_x = 0;
end

if calc_xyz(2)
  component_y = index_sum.*component_y;
else
  component_y = 0;
end

if calc_xyz(3)
  component_z = index_sum.*component_z;
else
  component_z = 0;
end

calc_out = J1*J2*magconst .* ...
  [ sum(component_x(:)) ;
    sum(component_y(:)) ;
    sum(component_z(:)) ] ;

debug_disp(calc_out')

end




function calc_out = forces_calc_z_x(size1,size2,offset,J1,J2)

calc_xyz = swap_x_y(calc_xyz);

forces_xyz = forces_calc_z_y(...
  swap_x_y(size1), swap_x_y(size2), rotate_x_to_y(offset),...
  J1, rotate_x_to_y(J2) );

calc_xyz = swap_x_y(calc_xyz);
calc_out = rotate_y_to_x( forces_xyz );

end



function calc_out = stiffnesses_calc_z_z(size1,size2,offset,J1,J2)

J1 = J1(3);
J2 = J2(3);


if (J1==0 || J2==0)
  debug_disp('Zero magnetisation.')
  calc_out  =  [0; 0; 0];
  return;
end

u = offset(1) + size2(1)*(-1).^index_j - size1(1)*(-1).^index_i;
v = offset(2) + size2(2)*(-1).^index_l - size1(2)*(-1).^index_k;
w = offset(3) + size2(3)*(-1).^index_q - size1(3)*(-1).^index_p;
r = sqrt(u.^2+v.^2+w.^2);


if calc_xyz(1) || calc_xyz(3)
  component_x = - r - (u.^2 .*v)./(u.^2+w.^2) - v.*log(r-v) ;
end

if calc_xyz(2) || calc_xyz(3)
  component_y = - r - (v.^2 .*u)./(v.^2+w.^2) - u.*log(r-u) ;
end

if calc_xyz(3)
  component_z = - component_x - component_y;
end


if calc_xyz(1)
  component_x = index_sum.*component_x;
else
  component_x = 0;
end

if calc_xyz(2)
  component_y = index_sum.*component_y;
else
  component_y = 0;
end

if calc_xyz(3)
  component_z = index_sum.*component_z;
else
  component_z = 0;
end

calc_out = J1*J2*magconst .* ...
  [ sum(component_x(:)) ;
    sum(component_y(:)) ;
    sum(component_z(:)) ] ;

debug_disp(calc_out')

end





function calc_out = stiffnesses_calc_z_y(size1,size2,offset,J1,J2)

J1 = J1(3);
J2 = J2(2);


if (J1==0 || J2==0)
  debug_disp('Zero magnetisation.')
  calc_out  =  [0; 0; 0];
  return;
end

u = offset(1) + size2(1)*(-1).^index_j - size1(1)*(-1).^index_i;
v = offset(2) + size2(2)*(-1).^index_l - size1(2)*(-1).^index_k;
w = offset(3) + size2(3)*(-1).^index_q - size1(3)*(-1).^index_p;
r = sqrt(u.^2+v.^2+w.^2);


if calc_xyz(1) || calc_xyz(3)
  component_x =  ((u.^2 .*v)./(u.^2 + v.^2)) + (u.^2 .*w)./(u.^2 + w.^2) ...
       - u.*atan1(v.*w,r.*u) + multiply_x_log_y( w , r + v ) + ...
       + multiply_x_log_y( v , r + w );
end

if calc_xyz(2) || calc_xyz(3)
  component_y = - v/2 + (u.^2 .*v)./(u.^2 + v.^2) - (u.*v.*w)./(v.^2 + w.^2) ...
       -  u.*atan1(u.*w,r.*v) - multiply_x_log_y( v , r + w );
end

if calc_xyz(3)
  component_z = - component_x - component_y;
end


if calc_xyz(1)
  component_x = index_sum.*component_x;
else
  component_x = 0;
end

if calc_xyz(2)
  component_y = index_sum.*component_y;
else
  component_y = 0;
end

if calc_xyz(3)
  component_z = index_sum.*component_z;
else
  component_z = 0;
end

calc_out = J1*J2*magconst .* ...
  [ sum(component_x(:)) ;
    sum(component_y(:)) ;
    sum(component_z(:)) ] ;

debug_disp(calc_out')

end





function calc_out = stiffnesses_calc_z_x(size1,size2,offset,J1,J2)

calc_xyz = swap_x_y(calc_xyz);

stiffnesses_xyz = stiffnesses_calc_z_y(...
  swap_x_y(size1), swap_x_y(size2), rotate_x_to_y(offset),...
  J1, rotate_x_to_y(J2) );

calc_xyz = swap_x_y(calc_xyz);
calc_out = swap_x_y(stiffnesses_xyz);

end



function calc_out = torques_calc_z_z(size1,size2,offset,lever,J1,J2)

br1 = J1(3);
br2 = J2(3);

if br1==0 || br2==0
  debug_disp('Zero magnetisation')
  calc_out = 0*offset;
  return
end

a1 = size1(1)/2;
b1 = size1(2)/2;
c1 = size1(3)/2;

a2 = size2(1)/2;
b2 = size2(2)/2;
c2 = size2(3)/2;

a = offset(1);
b = offset(2);
c = offset(3);

d = a+lever(1);
e = b+lever(2);
f = c+lever(3);

Tx=0;
Ty=0;
Tz=0;

for ii=[0,1]
  for jj=[0,1]
    for kk=[0,1]
      for ll=[0,1]
        for mm=[0,1]
          for nn=[0,1]
                        
            Cw=(-1)^mm.*c1-f;
            Cv=(-1)^kk.*b1-e;
            Cu=(-1)^ii.*a1-d;
            
            w=c-(-1)^mm.*c1+(-1)^nn.*c2;
            v=b-(-1)^kk.*b1+(-1)^ll.*b2;
            u=a-(-1)^ii.*a1+(-1)^jj.*a2;
            
            s=sqrt(u.^2+v.^2+w.^2);

            Ex=(1/8).*(-2.*Cw.*(-4.*v.*u+s.^2+2.*v.*s)-w.*...
              (-8.*v.*u+s.^2+8.*Cv.*s+6.*v.*s)+4.*(2.*Cv.*u.*...
              w.*acoth(u./s)+w.*(v.^2+2.*Cv.*v-w.*(2.*Cw+w))...
              *acoth(v./s)-u.*(2.*w.*(Cw+w).*atan(v./w)+2*v.*...
              (Cw+w).*log(s-u)+(w.^2+2.*Cw.*w-v.*(2.*Cv+v)).*...
              atan(u.*v./(w.*s))))+2.*(2.*Cw+w).*(u.^2+w.^2).*log(v+s));

            Ey=1/8*((2.*Cw+w).*u.^2-8.*u.*v.*(Cw+w)+8.*u.*v.*(Cw+w).*log(s-v)...
              +4.*Cw.*u.*s+6.*w.*s.*u+(2.*Cw+w).*(v.^2+w.^2)+...
              4.*w.*(w.^2+2.*Cw.*w-u.*(2.*Cu+u)).*acoth(u./s)+...
              4.*v.*(-2.*Cu.*w.*acoth(v./s)+2.*w.*(Cw+w).*atan(u./w)...
              +(w.^2+2.*Cw.*w-u.*(2.*Cu+u)).*atan(u.*v./(w.*s)))...
              -2.*(2.*Cw+w).*(v.^2+w.^2).*log(u+s)+8.*Cu.*w.*s);

            Ez=(1/36).*(-u.^3-18.*v.*u.^2-6.*u.*(w.^2+6.*Cu...
              .*v-3.*v.*(2.*Cv+v)+3.*Cv.*s)+v.*(v.^2+6.*(w.^2+...
              3.*Cu.*s))+6.*w.*(w.^2-3.*v.*(2.*Cv+v)).*atan(u./w)...
              -6.*w.*(w.^2-3.*u.*(2.*Cu+u)).*atan(v./w)-9.*...
              (2.*(v.^2+2.*Cv.*v-u.*(2.*Cu+u)).*w.*atan(u.*v./(w.*s))...
              -2.*u.*(2.*Cu+u).*v.*log(s-u)-(2.*Cv+v).*(v.^2-w.^2)...
              .*log(u+s)+2.*u.*v.*(2.*Cv+v).*log(s-v)+(2.*Cu+...
              u).*(u.^2-w.^2).*log(v+s)));

            Tx=Tx+(-1)^(ii+jj+kk+ll+mm+nn)*Ex;
            Ty=Ty+(-1)^(ii+jj+kk+ll+mm+nn)*Ey;
            Tz=Tz+(-1)^(ii+jj+kk+ll+mm+nn)*Ez;
            
          end
        end
      end
    end
  end
end

calc_out = real([Tx; Ty; Tz].*br1*br2/(16*pi^2*1e-7));
                        
end

function calc_out = torques_calc_z_y(size1,size2,offset,lever,J1,J2)

if J1(3)~=0 && J2(2)~=0
  error('Torques cannot be calculated for orthogonal magnets yet.')
end

calc_out = 0*offset;

end

function calc_out = torques_calc_z_x(size1,size2,offset,lever,J1,J2)

if J1(3)~=0 && J2(1)~=0
  error('Torques cannot be calculated for orthogonal magnets yet.')
end

calc_out = 0*offset;

end



function calc_out = forces_cyl_calc(size1,size2,h_gap,J1,J2)

% inputs

r1 = size1(1);
r2 = size2(1);

% implicit

z = nan(4,length(h_gap));
z(1,:) = -size1(2)/2;
z(2,:) =  size1(2)/2;
z(3,:) = h_gap - size2(2)/2;
z(4,:) = h_gap + size2(2)/2;

C_d = zeros(size(h_gap));
  
for ii = [1 2]
  
  for jj = [3 4]

    a1 = z(ii,:) - z(jj,:);
    a2 = 1 + ( (r1-r2)./a1 ).^2;
    a3 = sqrt( (r1+r2).^2 + a1.^2 );
    a4 = 4*r1.*r2./( (r1+r2).^2 + a1.^2 );
    
    [K, E, PI] = ellipkepi( a4./(1-a2) , a4 );
    
    a2_ind = a2 == 1;
    if any(a2_ind)
      % singularity at a2=1 (i.e., equal radii)
      PI_term(a2_ind) = 0;
      PI_term(~a2_ind) = (1-a1.^2/a3.^2).*PI;
    else
      PI_term = (1-a1.^2./a3.^2).*PI;
    end
    
    f_z = a1.*a2.*a3.*( K - E./a2 - PI_term );

    f_z(abs(a1)<eps)=0; % singularity at a1=0 (i.e., coincident faces)
    
    C_d = C_d + (-1)^(ii+jj).*f_z;

  end

end

calc_out = J1*J2/(8*pi*1e-7)*C_d;

end


function [k,e,PI] = ellipkepi(a,m)

% Complete elliptic integrals calculated with the arithmetric-geometric mean
% algorithms contained here: http://dlmf.nist.gov/19.8
%
% Valid for a <= 1 and m <= 1

a0 = 1;
g0 = sqrt(1-m);
s0 = m;
nn = 0;

p0 = sqrt(1-a);
Q0 = 1;
Q1 = 1;
QQ = Q0;

w1 = ones(size(m));

while max(Q1(:)) > eps

  % for Elliptic I
  a1 = (a0+g0)/2;
  g1 = sqrt(a0.*g0);

  % for Elliptic II
  nn = nn + 1;
  c1 = (a0-g0)/2;
  w1 = 2^nn*c1.^2;
  s0 = s0 + w1;

  % for Elliptic III
  rr = p0.^2+a0.*g0;
  p1 = rr./(2.*p0);
  Q1 = 0.5*Q0.*(p0.^2-a0.*g0)./rr;
  QQ = QQ+Q1;

  a0 = a1;
  g0 = g1;
  Q0 = Q1;
  p0 = p1;

end

k = pi./(2*a1);
e = k.*(1-s0/2);
PI = pi./(4.*a1).*(2+a./(1-a).*QQ);

im = find(m == 1);
if ~isempty(im)
  k(im) = inf;
  e(im) = ones(length(im),1);
  PI(im) = inf;
end

end




function out = multiply_x_log_y(x,y)
  out = x.*log(y);
  out(~isfinite(out))=0;
end


function out = atan1(x,y)
  out = zeros(size(x));
  ind = x~=0 & y~=0;
  out(ind) = atan(x(ind)./y(ind));
end



end

