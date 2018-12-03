%% Multipole array forces
%
% Calculating the force between multipole magnet arrays
% in a range of configurations.

%%

clc
datafile = 'data/planar_compare_x_data.mat';

%% Array setup

array_height = 0.01;
array_length = 0.05;

displ_steps = 31;
xrange = array_length*linspace(-1,1,displ_steps);
gap = repmat( [0; 0; 1.5*array_height], [1 displ_steps] );
displ = gap + [1; 0; 0]*xrange;

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

%% Single magnet

single1 = ...
  struct(...
  'dim',[array_length; array_length; array_height],...
  'magdir', [0; 0; 1], ...
  'magn', 1  ...
  );

single2 = single1;
single2.magdir = [0;0;-1];

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


%%

if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
else
  forces.linear  = multipoleforces(linear1,  linear2, displ, 'debug');
  forces.single  = magnetforces(single1,  single2, displ, 'debug');
  forces.halbach = multipoleforces(halbach1, halbach2, displ, 'debug');
  forces.patchwork = multipoleforces(patchwork1, patchwork2, displ,'debug');
  forces.quasi = multipoleforces(quasi1,quasi2,displ,'debug');
  save(datafile,'forces','xrange');
end

%% Plot

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  simple_graph = true;
else
  simple_graph = false;
end


for ii = [1 3]
  
  willfig(['planar-compare-x',num2str(ii)]); clf; hold on;
  
  plot(xrange,forces.linear(ii,:),   '.-','Tag','Linear Halbach');
  plot(xrange,forces.halbach(ii,:),  '--','Tag','Planar Halbach');
  plot(xrange,forces.quasi(ii,:),         'Tag','Quasi Halbach' );
  plot(xrange,forces.single(ii,:),   '.' ,'Tag','Single magnet' );
  plot(xrange,forces.patchwork(ii,:),'.--','Tag','Patchwork'     );
  
  set(gca,'box','on','ticklength',[0.02 0.05])
  
  xlabel('Displacement, m')
  ylabel('Force, N')
  
  if ~simple_graph
    draworigin
    colourplot(1,5:-1:1);
    labelplot('northeast','vertical')
    legendshrink
  end

end