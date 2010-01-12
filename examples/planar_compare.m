%% Multipole array forces
%
% Calculating the force between multipole magnet arrays
% in a range of configurations.

%%

clc
datafile = 'data/planar_compare_data.mat';

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  simple_graph = true;
else
  simple_graph = false;
end

%% Array setup

array_height = 0.01;
array_length = 0.05;

displ_steps = 20;
zrange = array_height*(1+linspace(0.001,1,displ_steps));
displ = [0; 0; 1]*zrange;

%% Linear array

fixed_array = ...
  struct(...
  'type','linear',  ...
  'face','up',      ...
  'length', array_length,   ...
  'width',  array_length,   ...
  'height', array_height,   ...
  'Nmag_per_wave', 4, ...
  'Nwaves', 1,        ...
  'magn', 1           ...
  );

float_array = fixed_array;
float_array.face = 'down';

forces.linear = multipoleforces(fixed_array, float_array, displ);

%% Planar Halbach array

fixed_array = ...
  struct(...
  'type','planar',  ...
  'face','up',      ...
  'length', array_length,   ...
  'width',  array_length,   ...
  'height', array_height,   ...
  'Nmag_per_wave', 4, ...
  'Nwaves', 1,        ...
  'magn', 1           ...
  );

float_array = fixed_array;
float_array.face = 'down';

forces.halbach = multipoleforces(fixed_array, float_array, displ);

%%

save(datafile,'forces','zrange');

%%

if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
end

%% Plot

willfig('planar-compare'); clf; hold on;

plot(zrange,forces.linear(3,:),'Tag','Linear Halbach');
plot(zrange,forces.halbach(3,:),'Tag','Planar Halbach');
set(gca,'box','on','ticklength',[0.02 0.05])
axis_tight

xlabel('Vertical $z$ displacement, m')
ylabel('Vertical force, N')
text( -0.02, 200,'$F_z$');
text( -0.03,-100,'$F_y$');

if ~simple_graph
  [h1 h2] = draworigin;
  set([h1 h2],'linestyle','--');
  colourplot
  % labelplot
  matlabfrag('fig/planar-compare','dpi',3200);
end

