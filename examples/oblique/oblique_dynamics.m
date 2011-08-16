function [ T, X, param ] = oblique_dynamics( varargin )
%
% Perform a dynamic simulation of the 3DOF oblique magnet system.
% Some informational messages are output. Sorry about that.
%
% INPUTS
% ======
%
% Takes many of the same options as OBLIQUE_FORCES3, as well as:
%
%   mass                 1
%   momentofinertia*
%   dampingratio         [0.1 0.1 0.1]
%   maxdispl             0.2
%   time                 [0 10]
%   perturb              [0.001 0.001 1] % units of [m m °]
%   RelTol & AbsTol      1e-4
%
% If the moment of inertia is not given, it will be automatically
% calculated with J = 1/3*m*L^2 where m is the mass and L is the lever arm.
%
% The 'maxdispl' parameter is used when constructing the curve for
% finding the equilibrium position. For very light masses you may wish to
% increase it.
%
% The time vector can either be of length two, in which case the time
% steps are calculated automatically, or can be an entire time vector
% over which the outputs are calculated.
%
% The actual dynamics that are solved depend on the initial perturbation!
% If the horiz or rotational state is not perturbed, it is assumed to be
% constrained. E.g., a perturbation of [0.001 0.001 0] will solve assuming
% vertical and horizontal translation *only*; no rotation is permitted.
% Thus, there are essentially four options:
%
%     [x y r] - no constraints
%     [x y 0] - constrain rotation
%     [0 y 0] - constrain rotation and horiz
%     [0 y r] - constrain horiz
% 
%
% OUTPUTS
% =======
%
%        T - time vector
%        X - state vector output
%    param - structure containing various plant parameters:
%
%         .y0              - equilibrium displacement (vertical)
%         .momentofinertia - if automatically calculated
%         .stiffness       - [1x3] stiffness in (x, y, r).
%         .damping         - [1x3] damping coefficients in (x, y, r).

p = inputParser;
p.addParamValue('mass',1);
p.addParamValue('momentofinertia',NaN);
p.addParamValue('dampingratio',[0.1 0.1 0.1]);
p.addParamValue('maxdispl',0.2);

p.addParamValue('time',[0 10]);
p.addParamValue('perturb',[1/1000 1/1000 1]);

p.addParamValue('RelTol',1e-4);
p.addParamValue('AbsTol',1e-4);

p.addParamValue('magn',1);
p.addParamValue('magangle',45);
p.addParamValue('unitlength',0.02);
p.addParamValue('gapratio',0.4);
p.addParamValue('magratio',0.4);
p.addParamValue('leverratio',2);

p.parse(varargin{:});

%%

fprintf('\n')

m = p.Results.mass;
time = p.Results.time;

J = p.Results.momentofinertia;
if(isnan(J))
  % assume the system is a uniformly distributed mass
  % pretty approximate but that's okay
  J = 1/3*m*(p.Results.leverratio*p.Results.unitlength)^2;
  fprintf('Calculated moment of inertia: %3.3f g m^2\n',1000*J)
end
param.momentofinertia = J;

zeta = p.Results.dampingratio;

displ = linspace(0.001, p.Results.maxdispl, 200);
magopt = {...
  'magn',      p.Results.magn,...
  'magangle',  p.Results.magangle,...
  'unitlength',p.Results.unitlength,...
  'gapratio',  p.Results.gapratio,...
  'magratio',  p.Results.magratio,...
  'leverratio',p.Results.leverratio...
};
f_m = oblique_forces3('displ',[0; 1; 0]*displ,'rotation',0,magopt{:});
y0 = interp1(f_m(2,:),displ,m*9.8);

param.y0 = y0;

x0 = 0;
r0 = 0; % degrees
X0 = [x0 y0 r0*pi/180];

% rotation input in degrees must be converted to radians for the dynamics!
dX = p.Results.perturb.*[1 1 pi/180];

incr = 1e-6;

%%
if(isnan(y0))
  error('Mass too great for this system. (Or maybe too light; check "maxdispl".)')
else
  fprintf('Equilibrium position is %2.2f mm\n',1000*y0)
end

dx = incr;
ky = -(interp1(displ,f_m(2,:),y0+dx)/9.8-m)/dx;

fprintf('Vertical stiffness is %2.2f N/m\n',ky)
if(ky<0)
  warning('oblique:unstablevert','Vertical forces are not stable.')
end

f_m2 = oblique_forces3('displ',[dx; y0; 0],'rotation',0,magopt{:});
k_m = -f_m2/dx;
kx = k_m(1);

fprintf('Horizontal stiffness is %2.2f N/m\n',kx)
if(kx<0)
  warning('oblique:unstablehoriz','Horizontal forces are not stable.')
end

dr = incr;
[~,t_m3] = oblique_forces3('displ',[0; y0; 0],'rotation',dr,magopt{:});
kr_m = -t_m3/dr;
kr = kr_m(3);

fprintf('Rotational stiffness is %2.2f Nm/rad\n',kr)
if(kr<0)
  warning('oblique:unstablerot','Torques are not stable.')
end

c = 2.*zeta.*sqrt([m*kx m*ky J*kr]);

param.stiffness = [kx ky kr];
param.damping = c;

options = odeset('RelTol',p.Results.RelTol,'AbsTol',p.Results.AbsTol);

tic
if dX(1)==0 && dX(3)==0
  disp('Constraining horiz. displacement and rotation')
  [T,X] = ode45(@oblique_solve_y,time,[y0+dX(2) 0],options);
elseif dX(3)==0  
  disp('Constraining rotation')
  [T,X] = ode45(@oblique_solve_xy,time,[x0+dX(1) 0 y0+dX(2) 0],options);
elseif dX(1)==0
  disp('Constraining horizontal displacement')
  [T,X] = ode45(@oblique_solve_yr,time,[y0+dX(2) 0 r0+dX(3) 0],options);
else
  disp('Unconstrained')
  [T,X] = ode45(@oblique_solve,time,[x0+dX(1) 0 y0+dX(2) 0 r0+dX(3) 0],options);
end
toc

%% Dynamic functions

  function dX = oblique_solve(~,X)
    % Our dynamic equation for displacement is:
    % $$ X1' = X2 ,
    %    X2' = -[0; g] + f_m(X1)/m - c/m X2 $$
    % For translation and something similar for rotation.
    %
    % So our states will be $[x; x'; y; y'; r; r']$
    
    [f_m t_m] = oblique_forces3('displ',[X(1); X(3); 0],'rotation',X(5)*180/pi,magopt{:});
    
    dX = zeros(6,1);
    dX(1) = X(2);
    dX(2) = f_m(1)/m - c(1)/m * X(2);
    dX(3) = X(4);
    dX(4) = -9.81 + f_m(2)/m - c(2)/m * X(4);
    dX(5) = X(6);
    dX(6) = t_m(3)/J - c(3)/J * X(6);
    
  end

  function dX = oblique_solve_y(~,X)
    % constrain rotation and horizontal displacement
    
    f_m = oblique_forces3('displ',[0; X(1); 0],'rotation',0,magopt{:});
    
    dX = zeros(2,1);
    dX(1) = X(2);
    dX(2) = -9.81 + f_m(2)/m - c(2)/m * X(2);

  end

  function dX = oblique_solve_xy(~,X)
    % constrain rotation
    
    f_m = oblique_forces3('displ',[X(1); X(3); 0],'rotation',0,magopt{:});
    
    dX = zeros(4,1);
    dX(1) = X(2);
    dX(2) = f_m(1)/m - c(1)/m * X(2);
    dX(3) = X(4);
    dX(4) = -9.81 + f_m(2)/m - c(2)/m * X(4);
    
  end

  function dX = oblique_solve_yr(~,X)
    % constrain horizontal displacement
        
    [f_m t_m] = oblique_forces3('displ',[0; X(1); 0],'rotation',X(3)*180/pi,magopt{:});
    
    dX = zeros(4,1);
    dX(1) = X(2);
    dX(2) = -9.81 + f_m(2)/m - c(2)/m * X(2);
    dX(3) = X(4);
    dX(4) = t_m(3)/J - c(3)/J * X(4);
    
  end

end

