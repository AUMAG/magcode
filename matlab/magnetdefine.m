%% magnetdefine()
%
% Create a pseudo-object representing a magnet or coil.

% \START

% \section{The \texttt{magnetdefine()} function}

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

if isfield(mag,'grade')
  if isfield(mag,'magn')
    error('Cannot specify both ''magn'' and ''grade''.')
  else
    mag.magn = grade2magn(mag.grade);
  end
end


if ~isfield(mag,'lever')
  mag.lever = [0; 0; 0];
else
  ss = size(mag.lever);
  if (ss(1)~=3) && (ss(2)==3)
    mag.lever = mag.lever.'; % attempt [3 M] shape
  end
end



if strcmp(mag.type,'cylinder')
  
  % default to +Z magnetisation
  if ~isfield(mag,'dir')
    if ~isfield(mag,'magdir')
      mag.dir    = [0 0 1];
      mag.magdir = [0 0 1];
    else
      mag.dir = mag.magdir;
    end
  else
    if ~isfield(mag,'magdir')
      mag.magdir = mag.dir;
    end
  end
  
  % convert from current/turns to equiv magnetisation:
  if ~isfield(mag,'magn')
    mag.magn = 4*pi*1e-7*mag.turns*mag.current/mag.dim(2);
  end
  
else
  
  if ~isfield(mag,'magdir')
    warning('Magnet direction ("magdir") not specified; assuming +z.')
    mag.magdir = [0; 0; 1];
  else
    mag.magdir = resolve_magdir(mag.magdir);
  end

end

mag.fndefined = true;

end

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

%\begin{mfunction}{resolve_magdir}

function magdir = resolve_magdir(magdir)

% Magnetisation directions are specified in either cartesian or spherical
% coordinates.
%
% We don't use Matlab's |sph2cart| here, because it doesn't calculate zero
% accurately (because it uses radians and |cos(pi/2)| can only be evaluated
% to machine precision of pi rather than symbolically).

  if numel(magdir) == 2
    theta = magdir(1);
    phi = magdir(2);
    magdir = [ cosd(phi)*cosd(theta); cosd(phi)*sind(theta); sind(phi) ];
  elseif numel(magdir) == 3
    if all(magdir == 0)
      magdir = [0; 0; 0]; % this looks redundant but ensures column vector
    else
      magdir = magdir(:)/norm(magdir);
    end
  elseif numel(magdir) == 1
    switch magdir
      case  'x'; magdir = [1;0;0];
      case  'y'; magdir = [0;1;0];
      case  'z'; magdir = [0;0;1];
      case '+x'; magdir = [1;0;0];
      case '+y'; magdir = [0;1;0];
      case '+z'; magdir = [0;0;1];
      case '-x'; magdir = [-1; 0; 0];
      case '-y'; magdir = [ 0;-1; 0];
      case '-z'; magdir = [ 0; 0;-1];
      otherwise, error('Magnetisation %s not understood.',magdir);
    end
  else
    error('magdir has wrong number of elements.')
  end


end

%\end{mfunction}
