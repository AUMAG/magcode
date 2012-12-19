%% Experimenting with a double-Halbach design
%

clear all
close all

%% Setup
%
% 5x5 array as usual

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

displ_steps = 201;
ymax = 0.1;
yrange = linspace(-ymax,ymax,displ_steps);
zgap = repmat([0; 0; 0.015],[1 displ_steps]);
displ = zgap + [0; 1; 0]*yrange;

%% Calculate or load data

datafile = 'data/double_halbach.mat';
if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
else
  tic
  dforces = multipoleforces(fixed_array, float_array, displ);
  dforcestoc = toc;
  fprintf('Time for calculation: %f\n',dforcestoc)
  save(datafile,'dforces','yrange','dforcestoc');
end

%% Plot

try
  willfig('planarh-halbach'); clf; hold on;
catch
  figure;
end

offset = 0.009;
offsetrange = find(round(1000*yrange)==round(1000*offset))-find(yrange==0);

indsum = find(yrange>=yrange(1)+offset & yrange<=yrange(end)-offset);
ind1 = indsum+offsetrange;
ind2 = indsum-offsetrange;

ddforces = dforces(:,ind1)+dforces(:,ind2);

%plot(yrange(indsum),dforces([2 3],ind1));
%plot(yrange(indsum),dforces([2 3],ind2));
plot(yrange(indsum),ddforces([2 3],:));
set(gca,'box','on','ticklength',[0.02 0.05])

xlabel('Horizontal $y$-displacement, m')
ylabel('Force, N')
text( -0.022, 300,'$F_z$','interpreter','LaTeX');
text( -0.031,-120,'$F_y$','interpreter','LaTeX');

try
  xlim(0.8*[-ymax ymax])
  ylim([-150 350])
  draworigin([0 0],'--')
  colourplot
  matlabfrag('fig/double-halbach');
end

