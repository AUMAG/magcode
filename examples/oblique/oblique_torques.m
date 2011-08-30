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

if exist('addArrowhead')==2
  plot_schematics_bool = true;
else
  warning('You must install the `<a href="https://github.com/zprime/fletcher">fletcher</a>` Matlab package to draw schematics with arrow vectors. Only results in graph form will be shown.')
  plot_schematics_bool = false;  
end

%% Schematics for paper
%
% Note here the option (used only here) "plotextras" that includes
% far more geometric information on the schematic.
%
% This produces Figure 8.

if plot_schematics_bool
  
  figname = 'angle-schem';
  willfig(figname,'huge'); clf
  
  magoptions = {...
    'displ',[0;1;0]*0.01,...
    'magangle',30,...
    'gapratio',0.5,...
    'magratio',0.4...
    };
  
  [t f] = oblique_forces3(...
    'plot',true,'plotextras',true,...
    'plotvecscale',0.005,...
    'plottorquerratio',0.3,...
    'rotation',15,...
    magoptions{:}...
    );
  
  matlabfrag(['fig/mbq-',figname])
  
end

%% Individual schematics showing the small angle approximation
%
% This produces the graphics for Figure 9 in the paper.

if plot_schematics_bool
  
  figname = 'angle-schem-1';
  willfig(figname,'small'); clf
  
  magoptions = {...
    'displ',[0;1;0]*0.01,...
    'magangle',30,...
    'gapratio',0.25,...
    'magratio',0.4,...
    'leverratio',1.5 ...
    };
  
  [t f] = oblique_forces3(...
    'plot',true,'plotunrotated',true,...
    'plotforces',false,...
    'rotation',0, ...
    magoptions{:}...
    );
  matlabfrag(['fig/mbq-',figname])
  
  figname = 'angle-schem-2';
  willfig(figname,'small'); clf
  [t f] = oblique_forces3(...
    'plot',true,'plotmagnetrotation',true,'plotunrotated',true,...
    'plotforces',false,...
    'rotation',15, ...
    magoptions{:}...
    );
  matlabfrag(['fig/mbq-',figname])
  
  figname = 'angle-schem-3';
  willfig(figname,'small'); clf
  [t f] = oblique_forces3(...
    'plot',true,'plotmagnetrotation',false,'plotmagnetapprox',true,...
    'plotforces',false,...
    'rotation',15, ...
    magoptions{:}...
    );
  matlabfrag(['fig/mbq-',figname])
  
end

%% Animation of forces under vertical displacement
%
% Just for visualisation purposes. And a bit of fun.
  
displ_range = linspace(0.005,0.02,50);

willfig('plotanim2');
figuresize(6,20,'centimeters')

[f t] = oblique_forces3(...
  'displ',[0;1;0]*displ_range,...
  'plot',plot_schematics_bool,...
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
%
% This code calculates the results being processing in the next two cells;
% the first plots the schematics and the second draws the graph.

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

if plot_schematics_bool
  
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
  
end

%% Graphs of torque vs rotation
%
% This produces Figure 10.

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

xlabel('Rotation $\mbqrotz$, deg.')
ylabel('Torque $\mbqptorque$, \si{mN.m}')

matlabfrag('fig/mbq-rot-ex')

