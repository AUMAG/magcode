%% Stability in 3DOF for the oblique magnet spring
%
%

clc
close all

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  colourplot = @(varargin) disp('');
  draworigin = @(varargin) disp('');
  matlabfrag = @(varargin) disp('');
  moreticks  = @(varargin) disp('');
  axistight  = @(varargin) disp('');
end

%% Constant parameters for the dynamics

mass = 3;
zeta = [0.2 0.2 0.2];

%% Dynamic analysis (only vertical)
%
% Just to show the dynamic simulation works and gives reasonable results.

[T1, X1] = oblique_dynamics(...
  'mass',mass,...
  'dampingratio',zeta,...
  'time',[0 4],...
  'RelTol',1e-5,'AbsTol',1e-5,...
  'perturb',[0 0.001 0]...
  );

willfig('dyn-y'); clf; hold on
plot(T1,1000*X1(:,1))
xlabel('Time, s');
ylabel('Vertical displ., mm')


%% Dynamic analysis (vertical & horizontal)
%
% Specifically choose parameters such that the horizontal stiffness is
% negative.
%
% The results are unstable.

[T2, X2, param2] = oblique_dynamics(...
  'mass',0.5*mass,...
  'dampingratio',zeta,...
  'time',[0 0.2],...
  'RelTol',1e-5,'AbsTol',1e-5,...
  'perturb',[0.001 0.001 0]...
  );

figname = 'dyn-xy-unstable';
willfig(figname,'small'); clf; hold on
plot(T2,1000*X2(:,1),T2,1000*(X2(:,3)-param2.y0))
legend('$x$','$y$','location','northwest')
legend boxoff
legendshrink
xlabel('Time, s')
ylabel('Relative displacement, mm')
axis tight


%% Dynamic analysis (vertical & horizontal)
%
% Specifically choose parameters such that the horizontal stiffness is
% positive.
%
% This produces Figure 12.

[T2, X2, param2] = oblique_dynamics(...
  'mass',mass,...
  'dampingratio',zeta,...
  'time',[0 12],...
  'RelTol',1e-5,'AbsTol',1e-5,...
  'perturb',[0.0015 0.0015 0]...
  );

figname = 'dyn-xy';
willfig(figname,'small'); clf; hold on
plot(T2,1000*X2(:,1),T2,1000*(X2(:,3)-param2.y0))
legend('$x$','$y$')
legend boxoff
legendshrink
xlabel('Time, s')
ylabel('Relative displacement, mm')
axistight;
yl = ylim;
colourplot;
matlabfrag(['fig/mbq-',figname])

figname = 'dyn-xy-map';
willfig(figname,'small'); clf; hold on
cplot(1000*X2(:,1),1000*(X2(:,3)-param2.y0))
xlabel('Horizontal, mm')
ylabel('Vertical, mm')
axis equal
ylim(yl);
draworigin;
matlabfrag(['fig/mbq-',figname])


%% Dynamic analysis (vertical & horizontal)
%
% Specifically choose parameters such that the horizontal stiffness is
% positive, but the perturbation leads to instability.

[T2, X2, param2] = oblique_dynamics(...
  'mass',mass,...
  'dampingratio',zeta,...
  'time',[0 12],...
  'RelTol',1e-5,'AbsTol',1e-5,...
  'perturb',[0.00173 0.00173 0]...
  );

figname = 'dyn-xy';
willfig(figname,'small'); clf; hold on
plot(T2,1000*X2(:,1),T2,1000*(X2(:,3)-param2.y0))
legend('$x$','$y$')
legend boxoff
legendshrink
xlabel('Time, s')
ylabel('Relative displacement, mm')
axis tight
colourplot;

figname = 'dyn-xy-map';
willfig(figname,'small'); clf; hold on
cplot(1000*X2(:,1),1000*(X2(:,3)-param2.y0))
xlabel('Horizontal, mm')
ylabel('Vertical, mm')
axis equal


%% Vertical and rotational
%
% This produces Figure 13.

[T3, X3, param] = oblique_dynamics(...
  'mass',mass,...
  'dampingratio',zeta,...
  'time',linspace(0,5,200),...
  'RelTol',1e-5,'AbsTol',1e-5,...
  'perturb',[0 0.0015 3]...
  );

figname = 'dyn-yr-y';
willfig(figname,'small'); clf; hold on
plot(T3,1000*(X3(:,1)-param.y0))
moreticks;
axistight;
colourplot;
xlabel('Time, s')
ylabel('Relative displacement, mm')
matlabfrag(['fig/mbq-',figname])

figname = 'dyn-yr-r';
willfig(figname,'small'); clf; hold on
plot(T3,180/pi*X3(:,3))
moreticks;
axistight;
colourplot;
xlabel('Time, s')
ylabel('Rotation, deg.')
matlabfrag(['fig/mbq-',figname])


%% Vertical and rotational
%
% This time unstable.

[T3, X3, param] = oblique_dynamics(...
  'mass',mass,...
  'dampingratio',zeta,...
  'time',linspace(0,5,200),...
  'RelTol',1e-5,'AbsTol',1e-5,...
  'perturb',[0 0.0015 5]...
  );

figname = 'dyn-yr-y';
willfig(figname,'small'); clf; hold on
plot(T3,1000*(X3(:,1)-param.y0))
moreticks;
axis tight
colourplot;
xlabel('Time, s')
ylabel('Relative displacement, mm')
%matlabfrag(['fig/mbq-',figname])

figname = 'dyn-yr-r';
willfig(figname,'small'); clf; hold on
plot(T3,180/pi*X3(:,3))
moreticks;
axis tight
colourplot;
xlabel('Time, s')
ylabel('Rotation, deg.')
%matlabfrag(['fig/mbq-',figname])


%% All three: Unconstrained is unstable!
%
% This produces Figure 14.

fprintf('\n\n')
[T4, X4, param4] = oblique_dynamics(...
  'mass',mass,...
  'dampingratio',zeta,...
  'time',[0 0.6],...
  'RelTol',1e-5,'AbsTol',1e-5,...
  'perturb',[1e-9 -0.001 1e-9]...
  );

figname = 'dyn-xyr-x';
willfig(figname,'tiny'); clf; hold on
plot(T4,1000*X4(:,1))
xlabel('Time, s')
ylabel('Displacement, mm')
colourplot;
axistight([],[],'y','+x');
matlabfrag(['fig/mbq-',figname])

figname = 'dyn-xyr-y';
willfig(figname,'tiny'); clf; hold on
plot(T4,1000*(X4(:,3)-param4.y0))
xlabel('Time, s')
ylabel('Relative displacement, mm')
colourplot;
axistight([],[],'y','+x');
matlabfrag(['fig/mbq-',figname])

figname = 'dyn-xyr-r';
willfig(figname,'tiny'); clf; hold on
plot(T4,180/pi*X4(:,5))
xlabel('Time, s')
ylabel('Rotation, deg.')
colourplot;
axistight([],[],'y','+x');
matlabfrag(['fig/mbq-',figname])





