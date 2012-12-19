%% Multipole array forces example with horizontal displacement
%

clear all
close all

ymax = 0.08;
yrange = linspace(-ymax,ymax,201);
zgap = repmat([0; 0; 0.015],[1 length(yrange)]);
displ = zgap + [0; 1; 0]*yrange;

%% Linear halbach array

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

datafile = 'data/linearh_halbach.mat';
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

datafile = 'data/planarh_halbach.mat';
if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
else
  pforces = multipoleforces(fixed_array, float_array, displ);
  save(datafile,'pforces');
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

datafile = 'data/planarh_quasi.mat';
if exist(datafile,'file')
  % Delete (or rename) the data file to re-run the calculations
  load(datafile)
else
  qforces = multipoleforces(fixed_array, float_array, displ);
  save(datafile,'qforces');
end


%% Plot all

% linear halbach

try
  willfig('planarh-linear','tiny'); clf; hold on
catch
  figure
end

plot(1000*yrange,0.001*lforces([2 3],:))

set(gca,'box','on','ticklength',[0.02 0.05])
xlabel('$y$-displacement, $\muphdispl$, mm')
ylabel('Force, kN')
text(-20,0.250,'$F_z$','interpreter','latex');
text( 15,0.250,'$F_y$','interpreter','latex');

try
  xlim(1000*[-ymax ymax])
  ylim(0.001*[-400 400])
  colourplot
  moreticks
  draworigin([0 0],'--')
  matlabfrag('fig/planarh-linear');
end

% planar halbach

try
  willfig('planarh-halbach','tiny'); clf; hold on;
catch
  figure;
end

plot(1000*yrange,1000\pforces([2 3],:));
set(gca,'box','on','ticklength',[0.02 0.05])

xlabel('$y$-displacement, $\muphdispl$, mm')
ylabel('Force, kN')
text( -20, 0.200,'$F_z$','interpreter','LaTeX');
text( -30,-0.100,'$F_y$','interpreter','LaTeX');

try
  xlim(1000*[-ymax ymax])
  ylim(1000\[-400 400])
  draworigin([0 0],'--')
  colourplot
  moreticks
  matlabfrag('fig/planarh-halbach');
end

% quasi

figname = 'planarh-quasi';
try
  willfig(figname,'tiny'); clf; hold on;
catch
  figure;
end

plot(1000*yrange,1000\qforces([2 3],:));
set(gca,'box','on','ticklength',[0.02 0.05])

xlabel('$y$-displacement, $\muphdispl$, mm')
ylabel('Force, kN')
text( -18, 0.200,'$F_z$','interpreter','LaTeX');
text( -15,-0.220,'$F_y$','interpreter','LaTeX');

try
  xlim(1000*[-ymax ymax])
  ylim(1000\[-400 400])
  draworigin([0 0],'--')
  colourplot
  moreticks
  matlabfrag(['fig/',figname]);
end

