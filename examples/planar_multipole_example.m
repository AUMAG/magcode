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
yrange = linspace(-0.08,0.08,displ_steps);
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
  willfig('planar-halbach'); clf; hold on;
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
  axistight
  draworigin([0 0],'--')
  colourplot
  matlabfrag('fig/planar-halbach','dpi',3200);
end

