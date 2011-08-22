%% Rotations and torques for the 3DOF oblique magnetic spring
%
% For now a hodge-podge of examples.

clc
close all

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  colourplot = @(varargin) disp('');
  draworigin = @(varargin) disp('');
  matlabfrag = @(varargin) disp('');
  moreticks = @(varargin) disp('');
end

%% Schematic for paper
%
% Note here the option (used only here) "plotextras" that includes
% far more geometric information on the schematic.

willfig('angle-schem','huge'); clf

[f t] = oblique_forces3(...
  'displ',[0;1;0]*0.01,...
  'plot',true,'plotextras',true,...
  'plotvecscale',0.005,...
  'plottorquerratio',0.3,...
  'magangle',30,...
  'gapratio',0.5,...
  'magratio',0.4,...
  'rotation',15 ...
);

matlabfrag('fig/mbq-angle-schem')



%% Animation of forces under vertical displacement
%
% Just for visualisation purposes. And a bit of fun.

displ_range = linspace(0.005,0.02,50);

willfig('plotanim2');
figuresize(6,20,'centimeters')

[f t] = oblique_forces3(...
  'displ',[0;1;0]*displ_range,...
  'plot',true,...
  'plotvecscale',0.005,...
  'magangle',30,...
  'gapratio',0,...
  'rotation',0....
);

willfig('y displ; y force')
plot(1000*displ_range,squeeze(f(2,:)))
xlabel('Vertical displacement, mm');
ylabel('Vertical force, N')


%% Effect of magnet offset on torques
%
% The effect of the lever arm was relatively neglible and omitted from the
% results.

clc
rot = 0:0.5:10;
vert = 0.01;
lever = 2;
gaps = [0.25 0.5 1];

[f t ft] = oblique_forces3(...
  'displ',[0;1;0]*vert,...
  'magangle',30,...
  'leverratio',lever,...
  'gapratio',gaps,...
  'rotation',rot...
);

tz = 1000*squeeze(t(3,:,:));


%% Schematics of the three magnet offsets

ta = [pi/2 pi/4 pi/4];
for ii = 1:3
willfig(['show-rot-',num2str(ii)]); clf; hold on
[f t] = oblique_forces3(...
  'displ',[0;vert;0],...
  'leverratio',lever,...
  'plot',true,...
  'plotsize',0.1,...
  'plottorquearc',ta(ii),...
  'magangle',30,...
  'gapratio',gaps(ii),...
  'rotation',rot(round(end/2))...
);
matlabfrag(['fig/mbq-rot-ex-diag-',num2str(ii)])
end


%% Graphs of torque vs rotation

ltext = @(x,y,s) text(x,y,num2str(s),'horizontalalignment','c','userdata',...
  ['matlabfrag:\fboxsep=2pt\colorbox{white}{$',num2str(s),'$}']);

willfig('rotation; z torque','small'); clf; hold on
ind = 17*ones(size(gaps));
for gg = 2:length(gaps)
  plot(rot,squeeze(tz(gg,:)))
  ltext(rot(ind(gg)),tz(gg,ind(gg))-0.5,num2str(gaps(gg)))
end
for gg = 1
  plot(rot,squeeze(tz(gg,:)))
  ltext(rot(ind(gg)),tz(gg,ind(gg)),['\mbqoffset/\mbqunit=',num2str(gaps(gg))])
end
for gg = [2 3 1]
  % these are the "fake" torques calculated from the direct forces
  plot(rot,1000*ft(gg,:),'--')
end
colourplot(2);
draworigin;
moreticks;

xlabel('$z$ rotation $\mbqrotz$, deg.')
ylabel('$z$ torque $\mbqptorque$, \si{mN.m}')

matlabfrag('fig/mbq-rot-ex')

