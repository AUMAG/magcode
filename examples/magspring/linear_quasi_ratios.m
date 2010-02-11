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
waves_array  = [2 4 8];

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
  
  displ_steps = 5;
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

end

%% Plot integrals of force vs. displacement

figname = 'ratios_forcesum';
willfig(figname); clf; hold on

fs = repmat(NaN,[length(ratios_array) length(waves_array)]);

for ww = 1:length(waves_array)
  fs(:,ww) = squeeze(sum(forces(3,:,:,ww),2));
  plot(ratios_array,fs(:,ww));
end

set(gca,'box','on','ticklength',[0.02 0.05])
colourplot

axis tight
axistight
draworigin([1 1],'--')

xlabel('Magnet size ratio $\mupqratio$','interpreter','none')
ylabel('Normalised force sum')

matlabfrag(['fig/',figname],'dpi',3200);

