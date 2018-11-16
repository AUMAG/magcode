%% Multipole array forces example with horizontal displacement
%

clear all
close all

xmax = 0.08;
xrange = linspace(-xmax,xmax,201);
zgap = repmat([0; 0; 0.015],[1 length(xrange)]);
displ = zgap + [1; 0; 0]*xrange;

%% Linear multipole array

fixed_array = ...
  struct(...
  'type','linear',    ...
  'align','y',        ...
  'face','up',        ...
  'msize',  [0.05 0.01 0.01],   ...
  'Nmag_per_wave', 4, ...
  'Nmag', 5,          ...
  'magn', 1,          ...
  'magdir_first', -90 ...
  );

float_array = fixed_array;
float_array.face = 'down';
float_array.magdir_first = 90;

datafile = 'data/linearhx_halbach.mat';
if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
else
  lforces = multipoleforces(fixed_array, float_array, displ);
  save(datafile,'lforces');
end



%% Planar multipole array
%
% Two "planar Halbach" arrays in opposition.
% Forces are plotted as a function of horizontal displacement
% with a fixed vertical gap.

fixed_array = ...
  struct(...
  'type','planar', ...
  'align', 'xy',   ...
  'face','up',        ...
  'msize',  [0.01 0.01 0.01],   ...
  'Nmag_per_wave', 4, ...
  'Nwaves', 1,   ...
  'magn', 1      ...
  );

float_array = fixed_array;
float_array.face = 'down';

datafile = 'data/planarhx_halbach.mat';
if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
else
  forces = multipoleforces(fixed_array, float_array, displ);
  save(datafile,'forces');
end


%% Quasi-Halbach multipole array
%
% Two "quasi Halbach" arrays in opposition.
% Forces are plotted as a function of horizontal displacement
% with a fixed vertical gap.

fixed_array = ...
  struct(...
  'type','quasi-halbach', ...
  'align', 'xy',   ...
  'face','up',        ...
  'msize',  [0.01 0.01 0.01],   ...
  'Nwaves', 1,   ...
  'magn', 1      ...
  );

float_array = fixed_array;
float_array.face = 'down';

datafile = 'data/planarhx_quasi.mat';
if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
else
  qforces = multipoleforces(fixed_array, float_array, displ);
  save(datafile,'qforces');
end


%% Plot

try
  willfig('planarhx-linear','tiny'); clf; hold on
catch
  figure
end

plot(1000*xrange,0.001*lforces([1 3],:))

set(gca,'box','on','ticklength',[0.02 0.05])
xlabel('Horiz.\ $\ax$ displacement, mm')
ylabel('Force, kN')
text(-40,0.250,'$F_z$','interpreter','latex');
text(-30,-0.120,'$F_x$','interpreter','latex');

try
  xlim(1000*[-xmax xmax])
  ylim(0.001*[-400 400])
  moreticks
  colourplot
  draworigin([0 0],'--')
  matlabfrag('fig/planarhx-linear','dpi',3200);
end


%% Plot

try
  willfig('planarhx-halbach','3across'); clf; hold on;
catch
  figure;
end

plot(xrange,forces([1 3],:));
set(gca,'box','on','ticklength',[0.02 0.05])

xlabel('Horizontal $y$-displacement, m')
%ylabel('Force, N')
set(gca,'yticklabel',[])
text( -0.02, 200,'$F_z$','interpreter','LaTeX');
text( -0.03,-100,'$F_y$','interpreter','LaTeX');

try
  xlim([-xmax xmax])
  ylim([-400 400])
  draworigin([0 0],'--')
  colourplot
  matlabfrag('fig/planarhx-halbach','dpi',3200);
end

%% Plot

figname = 'planarhx-quasi';
try
  willfig(figname,'3across'); clf; hold on;
catch
  figure;
end

plot(xrange,qforces([1 3],:));
set(gca,'box','on','ticklength',[0.02 0.05])

xlabel('Horizontal $y$-displacement, m')
%ylabel('Force, N')
pp = get(gca,'position');
set(gca,'yticklabel',[])
set(gca,'position',pp)
text( -0.018, 200,'$F_z$','interpreter','LaTeX');
text( -0.015,-220,'$F_y$','interpreter','LaTeX');

try
  xlim([-xmax xmax])
  ylim([-400 400])
  draworigin([0 0],'--')
  colourplot
  matlabfrag(['fig/',figname]);
end

