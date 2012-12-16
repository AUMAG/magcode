%% Multipole array forces example
%
% Reproducing similar linear Halbach array results of Allag, Yonnet & Latreche (2009).

%% Calculate

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

ymax = 0.08;
yrange = linspace(-ymax,ymax,201);
zgap = repmat([0; 0; 0.015],[1 length(yrange)]);
displ = zgap + [0; 1; 0]*yrange;

lforces = multipoleforces(fixed_array, float_array, displ,'y','z');

%% Plot

try
  willfig('allag-repro','small'); clf; hold on
catch
  figure
end

plot(yrange,lforces([2 3],:))

set(gca,'box','on','ticklength',[0.02 0.05])
xlabel('Horizontal $y$-displacement, m')
ylabel('Force, N')
text(-0.018,200,'$F_z$','interpreter','latex');
text( 0.015,250,'$F_y$','interpreter','latex');

try
  xlim([-ymax ymax])
  ylim([-400 400])
  colourplot
  draworigin([0 0],'--')
  matlabfrag('fig/allag-repro','dpi',3200);
end



%% Planar multipole array forces example
%
% Two "planar Halbach" arrays in opposition.
% Forces are plotted as a function of horizontal displacement
% with a fixed vertical gap.

%% Setup

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

displ_steps = 51;
ymax = 0.08;
yrange = linspace(-ymax,ymax,displ_steps);
zgap = repmat([0; 0; 0.015],[1 displ_steps]);
displ = zgap + [0; 1; 0]*yrange;

%% Calculate or load data

datafile = 'data/planar_multipole_example_data.mat';
if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
else
  forces = multipoleforces(fixed_array, float_array, displ);
  save(datafile,'forces','yrange');
end

%% Plot

try
  willfig('planar-halbach','small'); clf; hold on;
catch
  figure;
end

plot(yrange,forces([2 3],:));
set(gca,'box','on','ticklength',[0.02 0.05])

xlabel('Horizontal $y$-displacement, m')
ylabel('Force, N')
text( -0.02, 200,'$F_z$','interpreter','LaTeX');
text( -0.03,-100,'$F_y$','interpreter','LaTeX');

try
  xlim([-ymax ymax])
  ylim([-400 400])
  draworigin([0 0],'--')
  colourplot
  matlabfrag('fig/planar-halbach','dpi',3200);
end

