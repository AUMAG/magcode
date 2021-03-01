%% magnetdefine()
%
% Create a pseudo-object representing a magnet or coil.

% \START

% \begin{mfunction}{magnetdefine}

function [mag] = magnetdefine(varargin)


if nargin == 1
  mag = varargin{1};
else
  mag = struct(varargin{:});
end

if ~isfield(mag,'type')
  warning('Magnets should always define their "type". E.g., {''type'',''cuboid''} for a cuboid magnet.')
  if length(mag.dim) == 2
    mag.type = 'cylinder';
  else
    mag.type = 'cuboid';
  end
end

mag.iscoil = false;
if strcmp(mag.type,'coil')
  mag.iscoil = true;  
  mag.type = 'cylinder';
end

if isfield(mag,'grade')
  if isfield(mag,'magn')
    error('Cannot specify both ''magn'' and ''grade''.')
  else
    mag.magn = grade2magn(mag.grade);
  end
end

mag = attempt3Mvector(mag,'lever');
mag = make_unit_vector(mag,'magdir');
mag = make_unit_vector(mag,'dir');

% defaults

if ~isfield(mag,'lever')
  mag.lever = [0; 0; 0];
end

if ~isfield(mag,'position')
  mag.position = [0; 0; 0];
end

if ~isfield(mag,'euler_order')
  mag.euler_order = 'XYZ';
end

if isfield(mag,'rotation')
  if all(size(mag.rotation) == [3 3])
    % good
  elseif numel(mag.rotation) == 3
    mag.rotation = eulang2rotmat(mag.rotation,mag.euler_order);
  else
    error('Magnet rotation must either be 3 euler angles or 3x3 rotation matrix.')
  end
else
  mag.rotation = eye(3);
end


switch mag.type
  case 'cylinder', mag = definecylinder(mag);
  case 'cuboid',   mag = definecuboid(mag);
  otherwise
    error('Magnet type "%s" unknown',mag.type)
end


if isfield(mag,'magdir') && isfield(mag,'magn')
  mag.magM = mag.magdir*mag.magn;
  mag.dipolemoment = 1/(4*pi*1e-7)*mag.magM*mag.volume;
end

mag.fndefined = true;

end


function mag = definecuboid(mag)

    if ~isfield(mag,'magdir')
      warning('Magnet direction ("magdir") not specified; assuming +z.')
      mag.magdir = [0; 0; 1];
    else
      mag = make_unit_vector(mag,'magdir');
    end
    
    if isfield(mag,'dim')
      mag.volume = prod(mag.dim);
      mag.vertices = transpose([
          mag.dim(1)/2, mag.dim(2)/2, mag.dim(3)/2;
          mag.dim(1)/2,-mag.dim(2)/2, mag.dim(3)/2;
         -mag.dim(1)/2,-mag.dim(2)/2, mag.dim(3)/2;
         -mag.dim(1)/2, mag.dim(2)/2, mag.dim(3)/2;
          mag.dim(1)/2,-mag.dim(2)/2,-mag.dim(3)/2;
         -mag.dim(1)/2,-mag.dim(2)/2,-mag.dim(3)/2;
         -mag.dim(1)/2, mag.dim(2)/2,-mag.dim(3)/2;
          mag.dim(1)/2, mag.dim(2)/2,-mag.dim(3)/2
        ]);
    end

end


function mag = definecylinder(mag)
    
    % dimensions
    
    if isfield(mag,'radius') && isfield(mag,'height')
      mag.dim = [mag.radius(:); mag.height];
    elseif isfield(mag,'radiusInner') && isfield(mag,'radiusOuter') && isfield(mag,'height')
      mag.dim = [mag.radiusInner; mag.radiusOuter; mag.height];
    end
    
    if isfield(mag,'dim')
      if numel(mag.dim) == 2
        mag.isring = false;
        mag.radiusInner = 0;
        mag.radiusOuter = mag.dim(1);
        mag.radius      = mag.dim(1);
        mag.height      = mag.dim(2);
      elseif numel(mag.dim) == 3
        if ~mag.iscoil 
          mag.isring = true;
        else
          mag.isring = false;
        end
        mag.radiusInner = mag.dim(1);
        mag.radiusOuter = mag.dim(2);
        mag.height      = mag.dim(3);
      end
    else
      error('Dimension specifications appear inconsistent / not defined.')
    end
    
    if mag.radiusOuter <= mag.radiusInner
      error('Ring radii must be defined as [ri ro] with ro > ri.')
    end
    mag.volume = pi*(mag.radiusOuter^2-mag.radiusInner^2)*mag.height;

    
    % directions - default to +Z magnetisation
        
    if ~isfield(mag,'dir')
      if ~isfield(mag,'magdir')
        warning('Magnet direction and magnetisation direction ("dir" and "magdir") not specified; assuming +z for both.')
        mag.dir    = [0; 0; 1];
        mag.magdir = [0; 0; 1];
      else
        mag.dir = mag.magdir;
      end
    else
      if ~isfield(mag,'magdir')
        mag.magdir = mag.dir;
      end
    end
    
    % Convert from current/turns to equiv magnetisation. Note for a thick
    % coil this is the equiv magn for each layer of the shell!
    
    if mag.iscoil
      mag = coil_check_turns(mag);
    end
    
    if mag.iscoil
      mag.magn = 4*pi*1e-7*mag.turnsAxial*mag.current/mag.height;
    end
    

end



%\begin{mfunction}{attempt3Mvector}
% If a series of vectors (column arrays) are stacked they should create a
% [3 M] size array.
% If [M 3] is entered, transpose it.
% (If size [3 3], leave as is!)
function mag = attempt3Mvector(mag,vecname)

if isfield(mag,vecname)

  if (size(mag.(vecname),1)~=3) && (size(mag.(vecname),2)==3)
    mag.(vecname) = transpose(mag.(vecname)); % attempt [3 M] shape
  end
  
end

end
%\end{mfunction}


%\begin{mfunction}{grade2magn}
%  Magnet `strength' can be specified using either "magn" or "grade".
%  In the latter case, this should be a string such as "'N42'", from which
%  the |magn| is automatically calculated using the equation
%  \[
%    B_r = 2\sqrt{\mu_0 [BH]_{\mathrm{max}}}
%  \]
%  where $[BH]_{\mathrm{max}}$ is the numeric value given in the grade in \si{MG.Oe}.
%  I.e., an N42 magnet has $[BH]_{\mathrm{max}}=\SI{42}{MG.Oe}$.
%  Since $\SI{1}{MG.Oe}=100/(4\pi)\,\si{kJ/m^3}$, the calculation simplifies to
%  \[
%    B_r = 2\sqrt{ N/100 }
%  \]
%  where $N$ is the numeric grade in \si{MG.Oe}. Easy.

function magn = grade2magn(grade)

if isnumeric(grade)
  magn = 2*sqrt(grade/100);
else
  if strcmp(grade(1),'N')
    grade = grade(2:end);
  end
  magn = 2*sqrt(str2double(grade)/100);
end

end
%\end{mfunction}

%\begin{mfunction}{make_unit_vector}
function mag = make_unit_vector(mag,vecname)
% Magnetisation directions are specified in cartesian coordinates.
% Although they should be unit vectors, we don't assume they are.

if ~isfield(mag,vecname)
  return
end

vec_in = mag.(vecname);

if isnumeric(vec_in)
  
  if numel(vec_in) ~= 3
    error(['"',vecname,'" has wrong number of elements (should be 3x1 vector or string input like ''+x'').'])
  end
  norm_vec_in = norm(vec_in);
  if norm_vec_in < eps
    norm_vec_in = 1; % to avoid 0/0
  end
  vec = vec_in(:)/norm_vec_in;
  
elseif ischar(vec_in)
  
  switch vec_in
    case  'x'; vec = [1;0;0];
    case  'y'; vec = [0;1;0];
    case  'z'; vec = [0;0;1];
    case '+x'; vec = [1;0;0];
    case '+y'; vec = [0;1;0];
    case '+z'; vec = [0;0;1];
    case '-x'; vec = [-1; 0; 0];
    case '-y'; vec = [ 0;-1; 0];
    case '-z'; vec = [ 0; 0;-1];
    otherwise, error('Vector string %s not understood.',vec_in);
  end
  
else
  error('Strange input (should be 3x1 vector or string input like ''+x'').')
end

mag.(vecname) = vec;

end

% \end{mfunction}


%\begin{mfunction}{coil_check_turns}
% |'turns'| is UI. After this function all code should use |'turnsAxial'|
% and |'turnsRadial'| instead.
function mag = coil_check_turns(mag)

if isfield(mag,'magn')
  error('Do not specific "magn" for coils; use "turns" (or "turnsAxial" & "turnsRadial") instead.')
end
if ~isfield(mag,'current')
  error('"current" must be specified for coils.')
end
if ~isfield(mag,'turns') && ~isfield(mag,'turnsAxial')
  error('"turns" or "turnsAxial" must be specified for coils.')
end
if isfield(mag,'turns') && isfield(mag,'turnsAxial')
  error('Both "turns" or "turnsAxial" cannot be specified.')
end
if isfield(mag,'turns')
  mag.turnsAxial = mag.turns;
else
  if ~isfield(mag,'turnsRadial')
    mag.turnsRadial = 1;
  else
    if ~isfield(mag,'turnsAxial')
      error('"turnsAxial" must be defined if "turnsRadial" is.')
    end
  end
end

end
%\end{mfunction}


% \end{mfunction}



