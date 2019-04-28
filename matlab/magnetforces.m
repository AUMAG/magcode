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

calc_force_bool     = false;
calc_stiffness_bool = false;
calc_torque_bool    = false;
break_ind = numel(varargin);

for iii = 1:break_ind
  switch varargin{iii}
    case 'force',      calc_force_bool     = true;
    case 'stiffness',  calc_stiffness_bool = true;
    case 'torque',     calc_torque_bool    = true;
    case 'x',  warning("Options 'x','y','z' are no longer supported.");
    case 'y',  warning("Options 'x','y','z' are no longer supported.");
    case 'z',  warning("Options 'x','y','z' are no longer supported.");
    otherwise
      break_ind = iii;
      break
  end
end

if break_ind == numel(varargin)
  plain_opts = varargin;
  var_opts = {};
else
  plain_opts = varargin(1:(break_ind-1));
  var_opts   = varargin(break_ind:end);
end

all_methods = {'auto','dipole'};

ip = inputParser;
ip.addParameter('method','auto',@(x) any(validatestring(x,all_methods)));
ip.addParameter('figure',gcf);
ip.addParameter('draw',0);
ip.addParameter('drawN',2);
ip.addParameter('drawpath',false);
ip.addParameter('drawpathopt',{'--','color','black'});
ip.addParameter('markpath',false);
ip.addParameter('markpathN',10);
ip.addParameter('markpathopt',{'color','black','markerfacecolor','black'});
ip.parse(var_opts{:});

if ~calc_force_bool && ~calc_stiffness_bool && ~calc_torque_bool
  plain_opts = {'force'};
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

%% \subsubsection{Set up magnets}
%
% First of all, address the data structures required for the input and output.
% Because displacement of a single magnet has three components, plus sizes of
% the faces another three, plus magnetisation strength and direction (two) makes
% nine in total, we use a structure to pass the information into
% the function. Otherwise we'd have an overwhelming number of input arguments.
%
% The input variables |magnet.dim| should be the entire side lengths of the
% magnets; these dimensions are halved when performing all of the calculations.
% (Because that's just how the maths is.)

if ~isfield(magnet_fixed,'fndefined')
  magnet_fixed = magnetdefine(magnet_fixed);
end
if ~isfield(magnet_float,'fndefined')
  magnet_float = magnetdefine(magnet_float);
end

if strcmp(ip.Results.method, 'dipole')

  magtype = 'dipole';
  
else
  
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
    error('Magnets must be of same type (cuboid/cuboid or cylinder/cylinder)')
  end
  magtype = magnet_fixed.type;
  
end

%% \subsubsection{Magnet calculations}

switch magtype
  case 'dipole'
    
    [forces_out,torques_out] = dipole_forcetorque(magnet_fixed.dipolemoment,magnet_float.dipolemoment,displ);
    
  case 'cuboid'
    
    size1 = magnet_fixed.dim(:)/2;
    size2 = magnet_float.dim(:)/2;
    
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
    J1_x    = rotate_x_to_z(magnet_fixed.magM);
    J2_x    = rotate_x_to_z(magnet_float.magM);
    
    size1_y = swap_y_z(size1);
    size2_y = swap_y_z(size2);
    J1_y    = rotate_y_to_z(magnet_fixed.magM);
    J2_y    = rotate_y_to_z(magnet_float.magM);
    
    if calc_force_bool
      forces_out = cuboid_magnet_force(displ);
    end
    
    if calc_stiffness_bool
      for iii = 1:Ndispl
        stiffnesses_out(:,iii) = cuboid_magnet_stiffness(displ(:,iii));
      end
    end
    
    if calc_torque_bool
      torques_out = cuboid_magnet_torque(displ,magnet_float.lever);
    end
    
  case 'cylinder'
    
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
    
    cyldir    = find(magnet_float.magdir ~= 0);
    cylnotdir = find(magnet_float.magdir == 0);
    if length(cyldir) ~= 1
      error('Cylindrical magnets must be aligned in one of the x, y or z directions')
    end
    
    if calc_force_bool
      if magnet_fixed.isring && magnet_float.isring
        for iii = 1:Ndispl
          forces_out(:,iii) = ring_magnet_force(displ(:,iii));
        end
      else
        for iii = 1:Ndispl
          forces_out(:,iii) = cyl_magnet_force(displ(:,iii));
        end
      end
    end
    
    if calc_stiffness_bool
      error('Stiffness cannot be calculated for cylindrical magnets yet.')
    end
    
    if calc_torque_bool
      error('Torques cannot be calculated for cylindrical magnets yet.')
    end
    
    
  case 'coil'
    
    warning('Code for coils in Matlab has never been completed :( See the Mathematica code for more details!')
    for iii = 1:Ndispl
      forces_out(:,iii) = coil_sign*coil.dir*...
        forces_magcyl_shell_calc(...
        magnet.dim, ...
        coil.dim, ...
        squeeze(displ(cyldir,:)), ...
        magnet.magM(cyldir), ...
        coil.current, ...
        coil.turns);
    end
    
end


%% \subsubsection{Drawing}

markpath = ip.Results.markpath;
if ischar(ip.Results.markpath)
  markpath = {ip.Results.markpath,ip.Results.markpath,ip.Results.markpath};
end

if ip.Results.draw || ip.Results.drawpath || iscell(ip.Results.markpath)
  draw_everything(ip)
end


%% \subsubsection{Return all results}
%
% After all of the calculations have occured, they're placed back into
% |varargout|.
% Outputs are ordered in the same order as the inputs are specified, which
% makes the code a bit uglier but is presumably a bit nicer for the user
% and/or just a bit more flexible.
%
% TODO: with new inputParser this can be simplified

argcount = 0;

for iii = 1:length(plain_opts)
  switch plain_opts{iii}
    case 'force',     argcount = argcount+1;
    case 'stiffness', argcount = argcount+1;
    case 'torque',    argcount = argcount+1;
  end
end

varargout = cell(argcount,1);

argcount = 0;

for iii = 1:length(plain_opts)
  switch plain_opts{iii}
    case 'force',      argcount = argcount+1; varargout{argcount} = forces_out;
    case 'stiffness',  argcount = argcount+1; varargout{argcount} = stiffnesses_out;
    case 'torque',     argcount = argcount+1; varargout{argcount} = torques_out;
  end
end

% That is the end of the main function.


%% \subsection{Nested functions}

% \begin{mfunction}{draw_everything}

  function draw_everything(ip)
    
    figure(ip.Results.figure);
    fig_was_held_bool = ishold;
    if ~fig_was_held_bool, hold on; end
    
    % draw magnets
    if ip.Results.draw
      M = min(Ndispl,ip.Results.drawN);
      Mind = round(linspace(1,Ndispl,M));
      
      magnetdraw(magnet_fixed,[0;0;0])
      magnetdraw(magnet_float,displ(:,Mind))
      
    end
    
    % draw path
    if ip.Results.drawpath || iscell(ip.Results.markpath)
      plot3(displ(1,:),displ(2,:),displ(3,:),ip.Results.drawpathopt{:})
    end
    
    % draw markers
    if iscell(ip.Results.markpath)
      M = min(Ndispl,ip.Results.markpathN);
      Mind = round(linspace(1,Ndispl,M));
      
      plot3(displ(1,Mind(1)),      displ(2,Mind(1))      ,displ(3,Mind(1))      ,'linestyle','none','marker',markpath{1},ip.Results.markpathopt{:})
      plot3(displ(1,Mind(2:end-1)),displ(2,Mind(2:end-1)),displ(3,Mind(2:end-1)),'linestyle','none','marker',markpath{2},ip.Results.markpathopt{:})
      plot3(displ(1,Mind(end)),    displ(2,Mind(end))    ,displ(3,Mind(end))    ,'linestyle','none','marker',markpath{3},ip.Results.markpathopt{:})
      
    end
    
    if ~fig_was_held_bool, hold off; end
    
  end

% \end{mfunction}

% \begin{mfunction}{single_magnet_cyl_force}
  function forces_out = cyl_magnet_force(displ)

    forces_out = nan(size(displ));

    ecc = sqrt(sum(displ(cylnotdir).^2));

    if ecc < eps
      magdir = [0;0;0];
      magdir(cyldir) = 1;
      forces_out = magdir*cylinder_force_coaxial(magnet_fixed.magM(cyldir), magnet_float.magM(cyldir), magnet_fixed.dim(1), magnet_float.dim(1), magnet_fixed.dim(2), magnet_float.dim(2), displ(cyldir)).';
    else
      ecc_forces = cylinder_force_eccentric(magnet_fixed.dim, magnet_float.dim, displ(cyldir), ecc, magnet_fixed.magM(cyldir), magnet_float.magM(cyldir)).';
      forces_out(cyldir) = ecc_forces(2);
      forces_out(cylnotdir(1)) = displ(cylnotdir(1))/ecc*ecc_forces(1);
      forces_out(cylnotdir(2)) = displ(cylnotdir(2))/ecc*ecc_forces(1);
      % Need to check this division into components is correct...
    end

  end
% \end{mfunction}

% \begin{mfunction}{single_magnet_ring_force}
  function forces_out = ring_magnet_force(displ)

    forces_out = nan(size(displ));

    ecc = sqrt(sum(displ(cylnotdir).^2));

    if ecc < eps
      magdir = [0;0;0];
      magdir(cyldir) = 1;
      forces11 = magdir*cylinder_force_coaxial(-magnet_fixed.magM(cyldir), -magnet_float.magM(cyldir), magnet_fixed.dim(1), magnet_float.dim(1), magnet_fixed.dim(3), magnet_float.dim(3), displ(cyldir)).';
      forces12 = magdir*cylinder_force_coaxial(-magnet_fixed.magM(cyldir), +magnet_float.magM(cyldir), magnet_fixed.dim(1), magnet_float.dim(2), magnet_fixed.dim(3), magnet_float.dim(3), displ(cyldir)).';
      forces21 = magdir*cylinder_force_coaxial(+magnet_fixed.magM(cyldir), -magnet_float.magM(cyldir), magnet_fixed.dim(2), magnet_float.dim(1), magnet_fixed.dim(3), magnet_float.dim(3), displ(cyldir)).';
      forces22 = magdir*cylinder_force_coaxial(+magnet_fixed.magM(cyldir), +magnet_float.magM(cyldir), magnet_fixed.dim(2), magnet_float.dim(2), magnet_fixed.dim(3), magnet_float.dim(3), displ(cyldir)).';
      forces_out = forces11 + forces12 + forces21 + forces22;
    else
      ecc_forces = cylinder_force_eccentric(magnet_fixed.dim, magnet_float.dim, displ(cyldir), ecc, magnet_fixed.magM(cyldir), magnet_float.magM(cyldir)).';
      forces_out(cyldir) = ecc_forces(2);
      forces_out(cylnotdir(1)) = displ(cylnotdir(1))/ecc*ecc_forces(1);
      forces_out(cylnotdir(2)) = displ(cylnotdir(2))/ecc*ecc_forces(1);
      % Need to check this division into components is correct...
    end

  end
% \end{mfunction}

% \begin{mfunction}{single_magnet_force}

% The |x| and |y| forces require a rotation to get
% the magnetisations correctly aligned.
% In the case of the magnet sizes, the lengths are just flipped rather than
% rotated (in rotation, sign is important).
% After the forces are calculated, rotate them back to the original
% coordinate system.

  function force_out = cuboid_magnet_force(displ)

    force_out = zeros(size(displ));

    d_x  = rotate_x_to_z(displ);
    d_y  = rotate_y_to_z(displ);

    force_out = force_out + ...
      cuboid_force_x_x(size1,size2,displ,magnet_fixed.magM,magnet_float.magM);

    force_out = force_out + ...
      rotate_z_to_x( cuboid_force_z_y(size1_x,size2_x,d_x,J1_x,J2_x) );

    force_out = force_out + ...
      rotate_z_to_x( cuboid_force_z_x(size1_x,size2_x,d_x,J1_x,J2_x) );

    force_out = force_out + ...
      rotate_z_to_y( cuboid_force_z_x(size1_y,size2_y,d_y,J1_y,J2_y) );

    force_out = force_out + ...
      cuboid_force_y_y(size1,size2,displ,magnet_fixed.magM,magnet_float.magM);

    force_out = force_out + ...
      rotate_z_to_y( cuboid_force_z_y(size1_y,size2_y,d_y,J1_y,J2_y) );

    force_out = force_out + ...
      cuboid_force_z_z( size1,size2,displ,magnet_fixed.magM,magnet_float.magM );

    force_out = force_out + ...
      cuboid_force_z_y( size1,size2,displ,magnet_fixed.magM,magnet_float.magM );

    force_out = force_out + ...
      cuboid_force_z_x( size1,size2,displ,magnet_fixed.magM,magnet_float.magM );

  end

% \end{mfunction}


% \begin{mfunction}{single_magnet_torque}
%
% For the magnetforces code we always assume the first magnet is fixed.
% But the Janssen code assumes the torque is calculated on the first magnet
% and defines the lever arm for that first magnet. Therefore we need to
% flip the definitions a bit.

  function torques_out = cuboid_magnet_torque(displ,lever)

    torque_components = nan([size(displ) 9]);

    d_x  = rotate_x_to_z(displ);
    d_y  = rotate_y_to_z(displ);
    d_z  = displ;

    l_x = rotate_x_to_z(lever);
    l_y = rotate_y_to_z(lever);
    l_z = lever;

    torque_components(:,:,9) = cuboid_torque_z_z( size1,size2,d_z,l_z,magnet_fixed.magM,magnet_float.magM );

    torque_components(:,:,8) = cuboid_torque_z_y( size1,size2,d_z,l_z,magnet_fixed.magM,magnet_float.magM );

    torque_components(:,:,7) = torques_calc_z_x( size1,size2,d_z,l_z,magnet_fixed.magM,magnet_float.magM );

    torque_components(:,:,1) = ...
      rotate_z_to_x( cuboid_torque_z_z(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );

    torque_components(:,:,2) = ...
      rotate_z_to_x( cuboid_torque_z_y(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );

    torque_components(:,:,3) = ...
      rotate_z_to_x( torques_calc_z_x(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );

    torque_components(:,:,4) = ...
      rotate_z_to_y( torques_calc_z_x(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );

    torque_components(:,:,5) = ...
      rotate_z_to_y( cuboid_torque_z_z(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );

    torque_components(:,:,6) = ...
      rotate_z_to_y( cuboid_torque_z_y(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );

    torques_out = sum(torque_components,3);
  end

% \end{mfunction}


% \begin{mfunction}{single_magnet_stiffness}

  function stiffness_out = cuboid_magnet_stiffness(displ)

    stiffness_components = nan([9 3]);

    d_x  = rotate_x_to_z(displ);
    d_y  = rotate_y_to_z(displ);

    stiffness_components(7,:) = ...
      stiffnesses_calc_z_x( size1,size2,displ,magnet_fixed.magM,magnet_float.magM );

    stiffness_components(8,:) = ...
      stiffnesses_calc_z_y( size1,size2,displ,magnet_fixed.magM,magnet_float.magM );

    stiffness_components(9,:) = ...
      stiffnesses_calc_z_z( size1,size2,displ,magnet_fixed.magM,magnet_float.magM );

    stiffness_components(1,:) = ...
      swap_x_z( stiffnesses_calc_z_z( size1_x,size2_x,d_x,J1_x,J2_x ) );

    stiffness_components(2,:) = ...
      swap_x_z( stiffnesses_calc_z_y( size1_x,size2_x,d_x,J1_x,J2_x ) );

    stiffness_components(3,:) = ...
      swap_x_z( stiffnesses_calc_z_x( size1_x,size2_x,d_x,J1_x,J2_x ) );

    stiffness_components(4,:) = ...
      swap_y_z( stiffnesses_calc_z_x( size1_y,size2_y,d_y,J1_y,J2_y ) );

    stiffness_components(5,:) = ...
      swap_y_z( stiffnesses_calc_z_z( size1_y,size2_y,d_y,J1_y,J2_y ) );

    stiffness_components(6,:) = ...
      swap_y_z( stiffnesses_calc_z_y( size1_y,size2_y,d_y,J1_y,J2_y ) );

    stiffness_out = sum(stiffness_components);
  end

% \end{mfunction}





% \begin{mfunction}{stiffnesses_calc_z_z}

  function calc_out = stiffnesses_calc_z_z(size1,size2,offset,J1,J2)

    J1 = J1(3);
    J2 = J2(3);


    if (J1==0 || J2==0)
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

  end

% \end{mfunction}

% \begin{mfunction}{stiffnesses_calc_z_y}

  function calc_out = stiffnesses_calc_z_y(size1,size2,offset,J1,J2)

    J1 = J1(3);
    J2 = J2(2);


    if (J1==0 || J2==0)
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


% \subsubsection{Helpers}
%
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
