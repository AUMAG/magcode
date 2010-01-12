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


displ_steps = 30;
zrange = array_height*(1+linspace(0.001,1,displ_steps));
gap = repmat( [0; 0; 0], [1 displ_steps] );
displ = gap + [0; 0; 1]*zrange;

%% Linear

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

forces.linear  = multipoleforces(linear1,  linear2, displ, 'debug');

%% Single magnet

single1 = ...
  struct(...
  'dim',[array_length; array_length; array_height],...
  'magdir', [0; 0; 1], ...
  'magn', 1  ...
  );

single2 = single1;
single2.magdir = [0;0;-1];

forces.single  = magnetforces(single1,  single2, displ, 'debug');

%% Halbach

halbach1 = ...
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

halbach2 = halbach1;
halbach2.face = 'down';

forces.halbach = multipoleforces(halbach1, halbach2, displ);
disp('Halbach planar array calculations complete.')

%% Patchwork

patchwork1 = ...
  struct(...
  'type','generic',  ...
  'mcount',[5 5 1], ...
  'msize', array_height*[1 1 1],   ...
  'magdir_fn', @(ii,jj,kk) [0; 0; (-1)^(ii+jj+1)] , ...
  'magn', 1           ...
  );

patchwork2 = patchwork1;
patchwork2.magdir_fn = @(ii,jj,kk) [0; 0; (-1)^(ii+jj)];

forces.patchwork = multipoleforces(patchwork1, patchwork2, displ,'debug');

%% Quasi-Halbach

quasi1 = ...
  struct(...
  'type', 'generic', ...
  'mcount', [5 5 1], ...
  'msize', array_height*[1 1 1], ...
  'magn', 1, ...
  'magdir_fn', @(ii,jj,kk) [  sind(90*ii)*cosd(90*jj) ;
                              cosd(90*ii)*sind(90*jj) ;
                              sind(90*ii)*sind(90*jj) ] ...
  );

quasi2 = quasi1;
quasi2.magdir_fn = @(ii,jj,kk)  [  sind(90*ii)*cosd(90*jj) ;
                                   cosd(90*ii)*sind(90*jj) ;
                                  -sind(90*ii)*sind(90*jj) ] ;

forces.quasi = multipoleforces(quasi1,quasi2,displ,'debug');

%%

save(datafile,'forces','zrange');

%%

if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
end

%% Plot

willfig('planar-compare','large'); clf; hold on;

plot(zrange,forces.single(3,:),   ':' ,'Tag','Single magnet' );
plot(zrange,forces.patchwork(3,:),'-.','Tag','Patchwork'     );
plot(zrange,forces.linear(3,:),   '.-','Tag','Linear Halbach');
plot(zrange,forces.halbach(3,:),  '--','Tag','Planar Halbach');
plot(zrange,forces.quasi(3,:),         'Tag','Quasi Halbach' );

set(gca,'box','on','ticklength',[0.02 0.05])
set(gca,'xlim',[0.0095 0.02]);

xlabel('Vertical displacement, m')
ylabel('Vertical force, N')
text( -0.02, 200,'$F_z$');
text( -0.03,-100,'$F_y$');

if ~simple_graph
  h1 = draworigin([0.01 0],'v');
  set(h1,'linestyle','--');
  colourplot
  labelplot('northeast','vertical')
  legendshrink
  matlabfrag('fig/planar-compare','dpi',3200);
end

