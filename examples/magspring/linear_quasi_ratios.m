%% Linear quasi-Halbach arrays with varying sizes of magnet
%
%
% In which the vertically oriented magnets have different sizes to the
% horizontally oriented magnets.


%% Calculate

ratios_array = 0:0.025:1.5;

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

%% Plot all

N = length(ratios_array);

zdata = zrange/array_height;

figname = 'ratios_compare_all';
willfig(figname,'small'); clf; hold on;
  
for nn = 1:N
        
  plot(zdata,squeeze(forces(3,:,nn)),...
      'Tag',num2str( ratios_array(nn) ));

end

set(gca,'box','on','ticklength',[0.02 0.05])
colourplot

xlabel('Normalised vertical displacement')
ylabel('Vertical force, N')

H = labelplot('north','vertical','$\mupqratio$');
pos = get(H,'position');
set(H,'position',[0.6 0.5 pos(3:4)]);
legendshrink

matlabfrag(['fig/',figname],'dpi',3200);

%% Plot integral of force with displacement

figname = 'ratios_forcesum';
willfig(figname); clf; hold on

fs = squeeze(sum(forces(3,:,:),2));
 
plot(ratios_array,fs/fs(ratios_array==1));

set(gca,'box','on','ticklength',[0.02 0.05])
colourplot

axis tight
axistight
draworigin([1 1],'--')

xlabel('Magnet size ratio $\mupqratio$')
ylabel('Normalised force sum')

matlabfrag(['fig/',figname],'dpi',3200);

%% Plot a few

some = [0 0.5 1];
style = {'.-','-','--'};

zdata = zrange/array_height;

figname = 'ratios-compare';
willfig(figname); clf; hold on;
  
for ii = 1:length(some)
      
  nn = find(ratios_array==some(ii));
  plot(zdata,squeeze(forces(3,:,nn)),...
      style{ii},...
      'linewidth',1,...
      'Tag',num2str( ratios_array(nn) ));

end

set(gca,'box','on','ticklength',[0.02 0.05])
colourplot

xlabel('Normalised vertical displacement')
ylabel('Vertical force, N')

H = labelplot('north','vertical','$\mupqratio$');
pos = get(H,'position');
set(H,'position',[0.6 0.5 pos(3:4)]);
legendshrink(0.5)

matlabfrag(['fig/',figname],'dpi',3200);

%% Plot a few normalised

some = [0 0.5 1];
style = {'.-','-','--'};

zdata = zrange/array_height;

figname = 'ratios-compare-norm';
willfig(figname); clf; hold on;

mm = find(ratios_array==some(3));

for ii = 1:length(some)
      
  nn = find(ratios_array==some(ii));
  plot(zdata,squeeze(forces(3,:,nn))./squeeze(forces(3,:,mm)),...
      style{ii},...
      'linewidth',1,...
      'Tag',num2str( ratios_array(nn) ));

end

set(gca,'box','on','ticklength',[0.02 0.05])
colourplot

xlabel('Normalised vertical displacement')
ylabel('Vertical force, N')

H = labelplot('north','vertical','$\mupqratio$');
pos = get(H,'position');
set(H,'position',[0.6 0.5 pos(3:4)]);
legendshrink(0.5)

matlabfrag(['fig/',figname],'dpi',3200);