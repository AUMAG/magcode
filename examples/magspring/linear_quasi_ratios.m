%% Linear quasi-Halbach arrays with varying sizes of magnet
%
%
% In which the vertically oriented magnets have different sizes to the
% horizontally oriented magnets.
%
% Many calculations are made to optimise the exact ratio of widths
% to be used; hence slow.


%% Calculate

ratios_array = 0:0.025:1.5;
waves_array  = [1 2 4];

datafile = 'data/linear_quasi_ratios.mat';
if exist(datafile,'file')
  current_array = ratios_array;
  load(datafile);
  if isequal(current_array,ratios_array)
    recalculate = false;
  else
    ratios_array = current_array;
    recalculate = true;
  end
else
    recalculate = true;
end

%%

if recalculate
  
  displ_steps = 51;
  zrange = linspace(0.0101,0.02,displ_steps);
  displ = [0; 0; 1]*zrange;
  
  array_height = 0.01;
  array_length = 0.1;
  
  forces = repmat(NaN,...
    [3 displ_steps length(ratios_array) length(waves_array)]);
  
  for ww = 1:length(waves_array)
    for nn = 1:length(ratios_array)
      ratio = ratios_array(nn); 
      fixed_array = ...
        struct(...
        'type','linear-quasi',    ...
        'align','x',        ...
        'face','up',        ...
        'length', array_length,   ...
        'width',  array_height,   ...
        'height', array_height,   ...
        'Nwaves', waves_array(ww) ,  ...
        'ratio', ratio ...
        );
      float_array = fixed_array;
      float_array.face = 'down';
      forces(:,:,nn,ww) = multipoleforces(fixed_array, float_array, displ);
    end
  end
  save(datafile,'zrange','ratios_array','waves_array','array_height','array_length','forces');
else
  load(datafile);
end

%% Setup plot parameters
%
% When I create these plots I use some non-standard functions
% to make them pretty.

if ~exist('willfig','file')
  close all
  willfig = @(varargin) figure;
  simple_graph = true;
else
  simple_graph = false;
end

%% Plot integrals of force vs. displacement

figname = 'ratios-forcesum';
willfig(figname,'small'); clf; hold on

style = [0.5 1 1.5];

fs = repmat(NaN,[length(ratios_array) length(waves_array)]);
norm_ind = find(ratios_array==1);

for ww = 1:length(waves_array)
  fs(:,ww) = squeeze(sum(forces(3,:,:,ww),2))/squeeze(sum(forces(3,:,norm_ind,ww)));
  plot(ratios_array,fs(:,ww),...
    'Tag',num2str(waves_array(ww)),...
    'LineWidth',style(ww));
end

xlabel('Magnet size ratio $\mupqratio$','interpreter','none')
ylabel('Normalised force integral')

if ~simple_graph
  
  set(gca,'box','on','ticklength',[0.02 0.05])
  colourplot(1)
  H = labelplot('south','vertical','$\mupNwaves$');
  legendshrink
  
  axis tight
  axistight
  draworigin([1 1],':')
  
  matlabfrag(['fig/',figname],'dpi',3200);
  
end

