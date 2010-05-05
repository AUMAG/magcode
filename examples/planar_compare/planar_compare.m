%% Multipole array forces
%
% Calculating the force between multipole magnet arrays
% in a range of configurations.

%% Array setup

array_height = 0.01;
array_length = 0.05;

displ_steps = 30;
zrange = array_height*(1+linspace(0.001,1,displ_steps));
gap = repmat( [0; 0; 0], [1 displ_steps] );
displ = gap + [0; 0; 1]*zrange;

%%% Linear Halbach

linear1 = ...
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

linear2 = linear1;
linear2.face = 'down';

%%% Single magnet

single1 = ...
  struct(...
  'dim',[array_length; array_length; array_height],...
  'magdir', [0; 0; 1], ...
  'magn', 1  ...
  );

single2 = single1;
single2.magdir = [0;0;-1];

%%% Planar Halbach

halbach1 = ...
  struct(...
  'type','planar',  ...
  'face','up',      ...
  'msize', array_height*[1 1 1], ...
  'Nmag_per_wave', 4, ...
  'Nwaves', 1,        ...
  'magn', 1           ...
  );

halbach2 = halbach1;
halbach2.face = 'down';

%%% Patchwork

patchwork1 = ...
  struct(...
  'type','patchwork',  ...
  'face', 'up', ...
  'Nmag', 5, ...
  'msize', array_height*[1 1 1], ...
  'magn', 1           ...
  );

patchwork2 = patchwork1;
patchwork2.face = 'down';


%%% Quasi-Halbach

quasi1 = ...
  struct(...
  'type', 'quasi-halbach', ...
  'face', 'up', ...
  'Nwaves', 1, ...
  'msize', array_height*[1 1 1], ...
  'magn', 1 ...
  );

quasi2 = quasi1;
quasi2.face = 'down';


%% Load or calculate data

datafile = 'data/planar_compare_data.mat';
if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
else
  forces.single    = magnetforces(single1, single2, displ);
  forces.linear    = multipoleforces(linear1,    linear2,    displ);
  forces.halbach   = multipoleforces(halbach1,   halbach2,   displ);
  forces.patchwork = multipoleforces(patchwork1, patchwork2, displ);
  forces.quasi     = multipoleforces(quasi1,     quasi2,     displ);
  save(datafile,'forces','zrange');
end

%% Plot

try
  willfig('planar-compare'); clf; hold on
catch
  figure; hold on
end

plot(zrange,forces.linear(3,:),   '.-', 'Tag','Linear Halbach');
plot(zrange,forces.halbach(3,:),  '--', 'Tag','Planar Halbach');
plot(zrange,forces.quasi(3,:),          'Tag','Quasi Halbach' );
plot(zrange,forces.single(3,:),   '.' , 'Tag','Single magnet' );
plot(zrange,forces.patchwork(3,:),'.--','Tag','Patchwork'     );

set(gca,'box','on','ticklength',[0.02 0.05])
set(gca,'xlim',[0.0095 0.02]);

xlabel('Vertical displacement, m')
ylabel('Vertical force, N')

try
  draworigin([0.01 0],'v','--');
  colourplot(1,5:-1:1);
  labelplot('northeast','vertical')
  legendshrink
  matlabfrag('fig/planar-compare','dpi',3200);
end

