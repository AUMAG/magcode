%% MAGNETFORCES  Calculate forces between two cuboid magnets
%
% I haven't written proper inline help, sorry!
% Please read the PDF documentation instead for now.

% \START

% \section{The \texttt{magnetforces()} function}

function [varargout] = magnetforces(magnet_fixed, magnet_float, displ, varargin)

%% \subsection{Main function body}

magconst = 1/(4*pi*(4*pi*1e-7));

[index_i, index_j, index_k, index_l, index_p, index_q] = ndgrid([0 1]);

index_sum = (-1).^(index_i+index_j+index_k+index_l+index_p+index_q);


%% \subsubsection{Wrangling user input and output}
% We now have a choice of calculations to take based on the user input.
% This chunk and the next are used in both \texttt{magnetforces.m} and
% \texttt{multipoleforces.m}.

debug_disp = @(str) disp([]);
calc_force_bool     = false;
calc_stiffness_bool = false;
calc_torque_bool    = false;

for iii = 1:length(varargin)
  switch varargin{iii}
    case 'debug',      debug_disp = @(str) disp(str);
    case 'force',      calc_force_bool     = true;
    case 'stiffness',  calc_stiffness_bool = true;
    case 'torque',     calc_torque_bool    = true;
    case 'x',  warning("Options 'x','y','z' are no longer supported.");
    case 'y',  warning("Options 'x','y','z' are no longer supported.");
    case 'z',  warning("Options 'x','y','z' are no longer supported.");
    otherwise
      error(['Unknown calculation option ''',varargin{iii},''''])
  end
end

if ~calc_force_bool && ~calc_stiffness_bool && ~calc_torque_bool
  varargin{end+1} = 'force';
  calc_force_bool = true;
end

%% \subsubsection{Organise input displacements}
% Gotta check the displacement input for both functions.
% After sorting that out, we can initialise the output variables now we
% know how big they need to me.

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
  forces_out = nan([3 Ndispl]);
end

if calc_stiffness_bool
  stiffnesses_out = nan([3 Ndispl]);
end

if calc_torque_bool
  torques_out = nan([3 Ndispl]);
end

%% \subsubsection{Variables and data structures}
% First of all, address the data structures required for the input and output.
% Because displacement of a single magnet has three components, plus sizes of
% the faces another three, plus magnetisation strength and direction (two) makes
% nine in total, we use a structure to pass the information into
% the function. Otherwise we'd have an overwhelming number of input arguments.
%
% The input variables |magnet.dim| should be the entire side lengths of the
% magnets; these dimensions are halved when performing all of the calculations.
% (Because that's just how the maths is.)
%
% We use spherical coordinates to represent magnetisation angle, where |phi| is
% the angle from the horizontal plane ($-\pi/2\le \phi \le\pi/2$) and |theta|
% is the angle around the horizontal plane ($0\le\theta\le2\pi$). This follows
% Matlab's definition; other conventions are commonly used as well. Remember:
% \begin{quote}
% $(1,0,0)_{\text{cartesian}} \equiv (0,0,1)_{\text{spherical}}$\\
% $(0,1,0)_{\text{cartesian}} \equiv (\pi/2,0,1)_{\text{spherical}}$\\
% $(0,0,1)_{\text{cartesian}} \equiv (0,\pi/2,1)_{\text{spherical}}$
% \end{quote}
% Cartesian components can also be used as input as well, in which case they
% are made into a unit vector before multiplying it by the magnetisation
% magnitude.
% Either way (between spherical or cartesian input), |J1| and |J2| are made into
% the magnetisation vectors in cartesian coordindates.

if ~isfield(magnet_fixed,'fndefined')
  magnet_fixed = magnetdefine(magnet_fixed);
end
if ~isfield(magnet_float,'fndefined')
  magnet_float = magnetdefine(magnet_float);
end



if strcmp(magnet_fixed.type, 'coil')
  
  if ~strcmp(magnet_float.type, 'cylinder')
    error('Coil/magnet forces can only be calculated for cylindrical magnets.')
  end
  
  coil = magnet_fixed;
  magnet = magnet_float;
  magtype = 'coil';
  coil_sign = +1;
  
end

if strcmp(magnet_float.type, 'coil')
  
  if ~strcmp(magnet_fixed.type, 'cylinder')
    error('Coil/magnet forces can only be calculated for cylindrical magnets.')
  end
  
  coil = magnet_float;
  magnet = magnet_fixed;
  magtype = 'coil';
  coil_sign = -1;
  
end


if ~strcmp(magnet_fixed.type, magnet_float.type)
  error('Magnets must be of same type')
end
magtype = magnet_fixed.type;


if strcmp(magtype,'cuboid')
  
  size1 = magnet_fixed.dim(:)/2;
  size2 = magnet_float.dim(:)/2;
  
  J1 = magnet_fixed.magn*magnet_fixed.magdir;
  J2 = magnet_float.magn*magnet_float.magdir;

  swap_x_y = @(vec) vec([2 1 3],:);
  swap_x_z = @(vec) vec([3 2 1],:);
  swap_y_z = @(vec) vec([1 3 2],:);
  
  rotate_z_to_x = @(vec) [  vec(3,:);  vec(2,:); -vec(1,:) ] ; % Ry( 90)
  rotate_x_to_z = @(vec) [ -vec(3,:);  vec(2,:);  vec(1,:) ] ; % Ry(-90)
  
  rotate_y_to_z = @(vec) [  vec(1,:); -vec(3,:);  vec(2,:) ] ; % Rx( 90)
  rotate_z_to_y = @(vec) [  vec(1,:);  vec(3,:); -vec(2,:) ] ; % Rx(-90)
  
  rotate_x_to_y = @(vec) [ -vec(2,:);  vec(1,:);  vec(3,:) ] ; % Rz( 90)
  rotate_y_to_x = @(vec) [  vec(2,:); -vec(1,:);  vec(3,:) ] ; % Rz(-90)
  
  size1_x = swap_x_z(size1);
  size2_x = swap_x_z(size2);
  J1_x    = rotate_x_to_z(J1);
  J2_x    = rotate_x_to_z(J2);
  
  size1_y = swap_y_z(size1);
  size2_y = swap_y_z(size2);
  J1_y    = rotate_y_to_z(J1);
  J2_y    = rotate_y_to_z(J2);
  
elseif strcmp(magtype,'cylinder')
  
  size1 = magnet_fixed.dim(:);
  size2 = magnet_float.dim(:);
  
  if any(abs(magnet_fixed.dir) ~= abs(magnet_float.dir))
    error('Cylindrical magnets must be oriented in the same direction')
  end
  if any(abs(magnet_fixed.magdir) ~= abs(magnet_float.magdir))
    error('Cylindrical magnets must be oriented in the same direction')
  end
  if any(abs(magnet_fixed.dir) ~= abs(magnet_fixed.magdir))
    error('Cylindrical magnets must be magnetised in the same direction as their orientation')
  end
  if any(abs(magnet_float.dir) ~= abs(magnet_float.magdir))
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
  
  J1 = magnet_fixed.magn*magnet_fixed.magdir;
  J2 = magnet_float.magn*magnet_float.magdir;
  debug_disp('Magnetisation vectors:')
  debug_disp(J1)
  debug_disp(J2)
  
end


%% \subsubsection{Calculate for each displacement}
% The actual mechanics.
% The idea is that a multitude of displacements can be passed to the
% function and we iterate to generate a matrix of vector outputs.

if strcmp(magtype,'coil')
  
  for iii = 1:Ndispl
    forces_out(:,iii) = coil_sign*coil.dir*...
      forces_magcyl_shell_calc(...
        magnet.dim, ...
        coil.dim, ...
        squeeze(displ(cyldir,:)), ...
        J1(cyldir), ...
        coil.current, ...
        coil.turns);
  end
  
elseif strcmp(magtype,'cuboid')
  
  if calc_force_bool
    for iii = 1:Ndispl
      forces_out(:,iii) = single_magnet_force(displ(:,iii));
    end
  end
  
  if calc_stiffness_bool
    for iii = 1:Ndispl
      stiffnesses_out(:,iii) = single_magnet_stiffness(displ(:,iii));
    end
  end
  
  if calc_torque_bool
    torques_out = single_magnet_torque(displ,magnet_float.lever);
  end
  
elseif strcmp(magtype,'cylinder')
  
  if calc_force_bool
    for iii = 1:Ndispl
      forces_out(:,iii)  =  single_magnet_cyl_force(displ(:,iii));
    end
  end
  
  if calc_stiffness_bool
    error('Stiffness cannot be calculated for cylindrical magnets yet.')
  end
  
  if calc_torque_bool
    error('Torques cannot be calculated for cylindrical magnets yet.')
  end
  
end

%% \subsubsection{Return all results}
% After all of the calculations have occured, they're placed back into
% |varargout|. (This happens at the very end, obviously.)
% Outputs are ordered in the same order as the inputs are specified, which
% makes the code a bit uglier but is presumably a bit nicer for the user
% and/or just a bit more flexible.

argcount = 0;

for iii = 1:length(varargin)
  switch varargin{iii}
    case 'force',     argcount = argcount+1;
    case 'stiffness', argcount = argcount+1; 
    case 'torque',    argcount = argcount+1;
  end
end

varargout = cell(argcount,1);

argcount = 0;

for iii = 1:length(varargin)
  switch varargin{iii}
    case 'force',      argcount = argcount+1; varargout{argcount} = forces_out;
    case 'stiffness',  argcount = argcount+1; varargout{argcount} = stiffnesses_out;
    case 'torque',     argcount = argcount+1; varargout{argcount} = torques_out;
  end
end

% That is the end of the main function.


%% \subsection{Nested functions}

      
% \begin{mfunction}{single_magnet_cyl_force}
  function forces_out = single_magnet_cyl_force(displ)
    
    forces_out = nan(size(displ));
    
    ecc = sqrt(sum(displ(cylnotdir).^2));
    
    if ecc < eps
      debug_disp('Coaxial')
      magdir = [0;0;0];
      magdir(cyldir) = 1;
      forces_out = magdir*cylinder_force_coaxial(J1(cyldir), J2(cyldir), size1(1), size2(1), size1(2), size2(2), displ(cyldir)).';    
    else
      debug_disp('Non-coaxial')
      ecc_forces = cylinder_force_eccentric(size1, size2, displ(cyldir), ecc, J1(cyldir), J2(cyldir)).';  
      forces_out(cyldir) = ecc_forces(2);
      forces_out(cylnotdir(1)) = displ(cylnotdir(1))/ecc*ecc_forces(1);
      forces_out(cylnotdir(2)) = displ(cylnotdir(2))/ecc*ecc_forces(1);
      % Need to check this division into components is correct...
    end
    
  end
% \end{mfunction}

%\begin{mfunction}{single_magnet_force}

% The |x| and |y| forces require a rotation to get
% the magnetisations correctly aligned.
% In the case of the magnet sizes, the lengths are just flipped rather than
% rotated (in rotation, sign is important).
% After the forces are calculated, rotate them back to the original
% coordinate system.

  function force_out = single_magnet_force(displ)
    
    force_components = nan([9 3]);
    
    d_x  = rotate_x_to_z(displ);
    d_y  = rotate_y_to_z(displ);
    
    debug_disp('  ')
    debug_disp('CALCULATING THINGS')
    debug_disp('==================')
    debug_disp('Displacement:')
    debug_disp(displ')
    debug_disp('Magnetisations:')
    debug_disp(J1')
    debug_disp(J2')
        
    debug_disp('Forces x-x:')
    force_components(1,:) = ...
      rotate_z_to_x( cuboid_force_z_z(size1_x,size2_x,d_x,J1_x,J2_x) );
    
    debug_disp('Forces x-y:')
    force_components(2,:) = ...
      rotate_z_to_x( cuboid_force_z_y(size1_x,size2_x,d_x,J1_x,J2_x) );
    
    debug_disp('Forces x-z:')
    force_components(3,:) = ...
      rotate_z_to_x( cuboid_force_z_x(size1_x,size2_x,d_x,J1_x,J2_x) );
     
    debug_disp('Forces y-x:')
    force_components(4,:) = ...
      rotate_z_to_y( cuboid_force_z_x(size1_y,size2_y,d_y,J1_y,J2_y) );
    
    debug_disp('Forces y-y:')
    force_components(5,:) = ...
      rotate_z_to_y( cuboid_force_z_z(size1_y,size2_y,d_y,J1_y,J2_y) );
    
    debug_disp('Forces y-z:')
    force_components(6,:) = ...
      rotate_z_to_y( cuboid_force_z_y(size1_y,size2_y,d_y,J1_y,J2_y) );
        
    debug_disp('z-z force:')
    force_components(9,:) = cuboid_force_z_z( size1,size2,displ,J1,J2 );
    
    debug_disp('z-y force:')
    force_components(8,:) = cuboid_force_z_y( size1,size2,displ,J1,J2 );
    
    debug_disp('z-x force:')
    force_components(7,:) = cuboid_force_z_x( size1,size2,displ,J1,J2 );
    
    
    force_out = sum(force_components);
  end

% \end{mfunction}


% \begin{mfunction}{single_magnet_torque}

  function torques_out = single_magnet_torque(displ,lever)
    
    torque_components = nan([size(displ) 9]);
    
    
    d_x  = rotate_x_to_z(displ);
    d_y  = rotate_y_to_z(displ);
    
    l_x = rotate_x_to_z(lever);
    l_y = rotate_y_to_z(lever);
    
    
    debug_disp('  ')
    debug_disp('CALCULATING THINGS')
    debug_disp('==================')
    debug_disp('Displacement:')
    debug_disp(displ')
    debug_disp('Magnetisations:')
    debug_disp(J1')
    debug_disp(J2')
    
    
    debug_disp('Torque: z-z:')
    torque_components(:,:,9) = cuboid_torque_z_z( size1,size2,displ,lever,J1,J2 );
    
    debug_disp('Torque z-y:')
    torque_components(:,:,8) = torques_calc_z_y( size1,size2,displ,lever,J1,J2 );
    
    debug_disp('Torque z-x:')
    torque_components(:,:,7) = torques_calc_z_x( size1,size2,displ,lever,J1,J2 );
        
    debug_disp('Torques x-x:')
    torque_components(:,:,1) = ...
      rotate_z_to_x( cuboid_torque_z_z(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );
    
    debug_disp('Torques x-y:')
    torque_components(:,:,2) = ...
      rotate_z_to_x( torques_calc_z_y(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );
    
    debug_disp('Torques x-z:')
    torque_components(:,:,3) = ...
      rotate_z_to_x( torques_calc_z_x(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );
        
    debug_disp('Torques y-x:')
    torque_components(:,:,4) = ...
      rotate_z_to_y( torques_calc_z_x(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );
    
    debug_disp('Torques y-y:')
    torque_components(:,:,5) = ...
      rotate_z_to_y( cuboid_torque_z_z(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );
    
    debug_disp('Torques y-z:')
    torque_components(:,:,6) = ...
      rotate_z_to_y( torques_calc_z_y(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );    
    
    torques_out = sum(torque_components,3);
  end

%\end{mfunction}



  function stiffness_out = single_magnet_stiffness(displ)
    
    stiffness_components = nan([9 3]);
    
    
    d_x  = rotate_x_to_z(displ);
    d_y  = rotate_y_to_z(displ);
    
    
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
        
    debug_disp('x-x stiffness:')
    stiffness_components(1,:) = ...
      swap_x_z( stiffnesses_calc_z_z( size1_x,size2_x,d_x,J1_x,J2_x ) );
    
    debug_disp('x-y stiffness:')
    stiffness_components(2,:) = ...
      swap_x_z( stiffnesses_calc_z_y( size1_x,size2_x,d_x,J1_x,J2_x ) );
    
    debug_disp('x-z stiffness:')
    stiffness_components(3,:) = ...
      swap_x_z( stiffnesses_calc_z_x( size1_x,size2_x,d_x,J1_x,J2_x ) );
        
    debug_disp('y-x stiffness:')
    stiffness_components(4,:) = ...
      swap_y_z( stiffnesses_calc_z_x( size1_y,size2_y,d_y,J1_y,J2_y ) );
    
    debug_disp('y-y stiffness:')
    stiffness_components(5,:) = ...
      swap_y_z( stiffnesses_calc_z_z( size1_y,size2_y,d_y,J1_y,J2_y ) );
    
    debug_disp('y-z stiffness:')
    stiffness_components(6,:) = ...
      swap_y_z( stiffnesses_calc_z_y( size1_y,size2_y,d_y,J1_y,J2_y ) );
         
    stiffness_out = sum(stiffness_components);
  end






% \begin{mfunction}{stiffnesses_calc_z_z}

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
    
    component_x = - r - (u.^2 .*v)./(u.^2+w.^2) - v.*log(r-v) ;
    
    component_y = - r - (v.^2 .*u)./(v.^2+w.^2) - u.*log(r-u) ;
    
    component_z = - component_x - component_y;
    
    component_x = index_sum.*component_x;
    component_y = index_sum.*component_y;
    component_z = index_sum.*component_z;
    
    calc_out = J1*J2*magconst .* ...
      [ sum(component_x(:)) ;
      sum(component_y(:)) ;
      sum(component_z(:)) ] ;
    
    debug_disp(calc_out')
    
  end

% \end{mfunction}

% \begin{mfunction}{stiffnesses_calc_z_y}

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
    
    component_x =  ((u.^2 .*v)./(u.^2 + v.^2)) + (u.^2 .*w)./(u.^2 + w.^2) ...
      - u.*atan1(v.*w,r.*u) + multiply_x_log_y( w , r + v ) + ...
      + multiply_x_log_y( v , r + w );
    component_y = - v/2 + (u.^2 .*v)./(u.^2 + v.^2) - (u.*v.*w)./(v.^2 + w.^2) ...
      -  u.*atan1(u.*w,r.*v) - multiply_x_log_y( v , r + w );
    component_z = - component_x - component_y;

    component_x = index_sum.*component_x;
    component_y = index_sum.*component_y;
    component_z = index_sum.*component_z;
    
    calc_out = J1*J2*magconst .* ...
      [ sum(component_x(:)) ;
      sum(component_y(:)) ;
      sum(component_z(:)) ] ;
    
    debug_disp(calc_out')
    
    
    % \subsubsection{Helpers}
    %  The equations contain two singularities. Specifically, the equations
    %  contain terms of the form $x \log(y)$, which becomes |NaN| when both $x$
    %  and $y$ are zero since $\log(0)$ is negative infinity.
    %
    % \begin{mfunction}{multiply_x_log_y}
    %  This function computes $x \log(y)$, special-casing the singularity to output
    %  zero, instead. (This is indeed the value of the limit.)
    function out = multiply_x_log_y(x,y)
      out = x.*log(y);
      out(~isfinite(out))=0;
    end
    % \end{mfunction}
    
    % \begin{mfunction}{atan1}
    % We're using |atan| instead of |atan2| (otherwise the wrong results
    %	are calculated --- I guess I don't totally understand that), which becomes
    %	a problem when trying to compute |atan(0/0)| since |0/0| is |NaN|.
    function out = atan1(x,y)
      out = zeros(size(x));
      ind = x~=0 & y~=0;
      out(ind) = atan(x(ind)./y(ind));
    end
    % \end{mfunction}

  end

% \end{mfunction}

% \begin{mfunction}{stiffnesses_calc_z_x}

  function calc_out = stiffnesses_calc_z_x(size1,size2,offset,J1,J2)
        
    stiffnesses_xyz = stiffnesses_calc_z_y(...
      swap_x_y(size1), swap_x_y(size2), rotate_x_to_y(offset),...
      J1, rotate_x_to_y(J2) );
    
    calc_out = swap_x_y(stiffnesses_xyz);
    
  end

% \end{mfunction}

% \begin{mfunction}{torques_calc_z_y}
  function calc_out = torques_calc_z_y(size1,size2,offset,lever,J1,J2)
    
    if J1(3)~=0 && J2(2)~=0
      error('Torques cannot be calculated for orthogonal magnets yet.')
    end
    
    calc_out = 0*offset;
    
  end
% \end{mfunction}

% \begin{mfunction}{torques_calc_z_x}
  function calc_out = torques_calc_z_x(size1,size2,offset,lever,J1,J2)
    
    if J1(3)~=0 && J2(1)~=0
      error('Torques cannot be calculated for orthogonal magnets yet.')
    end
    
    calc_out = 0*offset;
    
  end
% \end{mfunction}


% \begin{mfunction}{forces_magcyl_shell_calc}
  function Fz = forces_magcyl_shell_calc(magsize,coilsize,displ,Jmag,Nrz,I)
    
    Jcoil = 4*pi*1e-7*Nrz(2)*I/coil.dim(3);
    
    shell_forces = nan([length(displ) Nrz(1)]);
    
    for rr = 1:Nrz(1)
      
      this_radius = coilsize(1)+(rr-1)/(Nrz(1)-1)*(coilsize(2)-coilsize(1));
      shell_size = [this_radius, coilsize(3)];
      
      shell_forces(:,rr) = cylinder_force_coaxial(magsize,shell_size,displ,Jmag,Jcoil);
      
    end
    
    Fz = sum(shell_forces,2);
    
  end
% \end{mfunction}


end
