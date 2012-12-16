
%% Array setup
%
% Linear Halbach array with half-cube magnets.

clear all

array_height = 0.01;
mag_length = array_height;

displ_steps = 50;
zrange = linspace(0.0001,array_height/2,displ_steps);
gap = repmat( [0; 0; array_height], [1 displ_steps] );
displ = gap + [0; 0; 1]*zrange;

linear1 = ...
  struct(...
  'type','linear',    ...
  'face','up',        ...
  'msize', [mag_length array_height array_height],   ...
  'Nmag_per_wave', 4, ...
  'Nmag',9,           ...
  'magn', 1,          ...
  'mgap',[0; 0; 0]    ...
  );

linear2 = linear1;
linear2.face = 'down';

single_mag = struct(...
  'dim',[mag_length array_height array_height],...
  'magn',1,...
  'magdir',[0 0 1]...
);
single_mag2 = single_mag;
single_mag2.magdir = [0 0 -1];

single_array = struct(...
  'dim',[9*mag_length array_height array_height],...
  'magn',1,...
  'magdir',[0 0 1]...
);
single_array2 = single_array;
single_array2.magdir = [0 0 -1];

gaps = mag_length*[0 0.1 0.25 0.5 1];

%% Data

datafile = 'data/multipole_gaps.mat';
if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
else
  tic
  for g = 1:numel(gaps)
    linear1.mgap(1) = gaps(g);
    linear2.mgap(1) = gaps(g);
    forces(g).linear = multipoleforces(linear1, linear2, displ);
  end
  single_force = magnetforces(single_mag,single_mag2,displ);
  array_force = magnetforces(single_array,single_array2,displ);
  toc
  save(datafile,'forces','zrange','gaps','single_force','array_force');
end

%% Plot

try
  willfig('halbach-gaps'); clf; hold on
catch
  figure; hold on
end

for g = 1:numel(gaps)
  plot(zrange/mag_length,forces(g).linear(3,:),   '-', 'Tag',['\num{',num2str(gaps(g)/mag_length),'}']);
end
%plot(zrange/mag_length,array_force(3,:),'k-','UserData','colourplot:ignore');
plot(zrange/mag_length,linear1.Nmag*single_force(3,:),'k--','UserData','colourplot:ignore');

set(gca,'box','on','ticklength',[0.02 0.05])
xlim([0 array_height/2]/mag_length);
ylim([0 350])

xlabel('Normalised vertical displacement $\mupvdispl/\mupheight$')
ylabel('Vertical force, N')

try
  hl=labelplot('east','vertical','$\mupfacegap/\mupmaglength$');
  lpos = get(hl,'Position');
  set(hl,'Position',lpos+[0.1 0.1 0 0]);
  
  colourplot
  legendshrink
  matlabfrag('fig/halbach-gaps');
end