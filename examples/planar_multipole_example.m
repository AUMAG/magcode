%% Planar multipole array forces example
%
% In which magnetisations are a function of magnet position in x and y.

%% Setup
%
% In case you don't have the various bits'n'pieces that I use to create
% my Matlab graphics (probably likely).

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  simple_graph = 1;
else
  simple_graph = 0;
end

%% Calculate

fixed_array = ...
  struct(...
  'type','planar-xy', ...
  'face','up',        ...
  'msize',  [0.01 0.01 0.01],   ...
  'Nmag_per_wave', [4 4], ...
  'Nmag', [5 5],          ...
  'magn', 1               ...
  );

float_array = fixed_array;
float_array.face = 'down';
float_array.magdir_first = [-90 -90];

displ_steps = 201;
yrange = linspace(-0.08,0.08,displ_steps);
zgap = repmat([0; 0; 0.015],[1 displ_steps]);
displ = zgap + [0; 1; 0]*yrange;

forces = multipoleforces(fixed_array, float_array, displ);

%% Plot

willfig('planar-patchwork'); clf; hold on;

plot(yrange,forces(2,:),'Tag','y');
plot(yrange,forces(3,:),'Tag','z');
set(gca,'box','on')
set(gca,'ticklength',[0.02 0.05])
axis tight

% We want the vertical axis 5% less tight:
ylim = get(gca,'ylim');
ylim_range = ylim(2)-ylim(1);
yp = 0.05;
ylim = [ylim(1)-yp*ylim_range ylim(2)+yp*ylim_range];
set(gca,'ylim',ylim);


xlabel('Horizontal $y$ displacement, m')
ylabel('Force, N')
text( -0.02, 200,'$F_z$');
text( -0.03,-100,'$F_y$');

if ~simple_graph
  [h1 h2] = draworigin;
  set([h1 h2],'linestyle','--');
  colourplot
  % labelplot
  matlabfrag('fig/planar-patchwork','dpi',3200);
end

% save data/multipole_planar_example_data forces drange
