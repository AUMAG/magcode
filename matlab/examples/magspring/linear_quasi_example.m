%% Linear quasi-Halbach arrays with varying sizes of magnet
%
%
% In which the vertically oriented magnets have different sizes to the
% horizontally oriented magnets.


%% Calculate

ratios_array = [0 0.5 1];

datafile = 'data/linear_quasi_example.mat';
if exist(datafile,'file')
  
  current_array = ratios_array;
  load(datafile);
  if isequal(current_array,ratios_array)
    recalculate = false;
  else
    ratios_array = current_array;
    recalculate = true;
  end
end

if recalculate
  
  displ_steps = 50;
  zrange = linspace(0.0101,0.02,displ_steps);
  displ = [0; 0; 1]*zrange;
  
  array_height = 0.01;
  array_length = 0.1;
  
  forces = repmat(NaN,...
    [3 displ_steps length(ratios_array)]);
  
  
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
      'Nwaves', 2 ,  ...
      'ratio', ratio ...
      );
    
    float_array = fixed_array;
    float_array.face = 'down';
    
    forces(:,:,nn) = multipoleforces(fixed_array, float_array, displ);
    
  end
  
  save(datafile,'zrange','ratios_array','array_height','array_length','forces');

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

%% Plot results

style = {'.-','-','--'};

zdata = zrange/array_height;

figname = 'ratios-compare';
willfig(figname,'tiny'); clf; hold on;
  
for ii = 1:length(ratios_array)
      
  plot(zdata,squeeze(forces(3,:,ii)),...
      style{ii},...
      'linewidth',1,...
      'Tag',num2str( ratios_array(ii) ));

end

xlabel('Normalised vertical displacement')
ylabel('Vertical force, N')

if ~simple_graph
  
  set(gca,'box','on','ticklength',[0.02 0.05])
  
  H = labelplot('north','vertical','$\mupqratio$');
  pos = get(H,'position');
  set(H,'position',[0.65 0.5 pos(3:4)]);
  legendshrink(0.5)
  colourplot

  axistight(gca,0.1,'-x')
  set(gca,'ylim',[0 400]);
  
  draworigin([1 0],'v',':')
  
  matlabfrag(['fig/',figname],'dpi',3200);

end

%% Plot a few normalised

style = {'.-','-','--'};

zdata = zrange/array_height;

figname = 'ratios-compare-norm';
willfig(figname,'tiny'); clf; hold on;

for ii = 1:length(ratios_array)
      
  plot(zdata,squeeze(forces(3,:,ii))./squeeze(forces(3,:,3)),...
      style{ii},...
      'linewidth',1,...
      'Tag',num2str( ratios_array(nn) ));

end

xlabel('Normalised vertical displacement')
ylabel('Vertical force, N')

if ~simple_graph
  
  set(gca,'box','on','ticklength',[0.02 0.05])
  colourplot
  
  H = labelplot('north','vertical','$\mupqratio$');
  pos = get(H,'position');
  set(H,'position',[0.6 0.5 pos(3:4)]);
  legendshrink(0.5)
  
  matlabfrag(['fig/',figname],'dpi',3200)

end