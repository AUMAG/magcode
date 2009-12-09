%% Multipole array forces
%
% Calculating the force between multipole magnet arrays
% in a range of configurations.

%%

clc
timestamp
tic

%% Calculate

if exist('multipole_compare_data.mat')

  load multipole_compare_data.mat
  
else
  
  displ_steps = 50;
  zrange = linspace(0.0101,0.02,displ_steps);
  
  Nmag_per_wave_array = [1 2 4 8];
  Nwaves_array = [1 2 4];
  
  array_height = 0.01;
  array_length = 0.1;
  
  forces = repmat(NaN,...
    [length(Nwaves_array) length(Nmag_per_wave_array) displ_steps 3]);
  
  
  for nn = 1:length(Nwaves_array)
    
    Nwaves = Nwaves_array(nn);
    
    for ww = 1:length(Nmag_per_wave_array)
      
      Nmag_per_wave = Nmag_per_wave_array(ww);
      
      fixed_array = ...
        struct(...
        'type','linear-x',  ...
        'face','up',        ...
        'length', array_length,   ...
        'width',  array_height,   ...
        'height', array_height,   ...
        'Nmag_per_wave', Nmag_per_wave, ...
        'Nwaves', Nwaves,   ...
        'magn', 1,          ...
        'magdir_first', 90  ...
        );
      
      float_array = fixed_array;
      float_array.face = 'down';
      float_array.magdir_first = -90;
      
      for ii = 1:displ_steps
        displ = [0 0 zrange(ii)];
        forces(nn,ww,ii,:) = multipoleforces(fixed_array, float_array, displ);
      end
      
    end
    
  end
  
  toc
  
  save multipole_compare_data displ_steps zrange ...
    Nmag_per_wave_array Nwaves_array array_height array_length forces

end

%% Plot

for nn = 1:length(Nwaves_array)
  
  Nwaves = Nwaves_array(nn);
  
  figname = ['halbach-waves-Nwaves-',num2str(Nwaves)];
  willfig(figname); clf; hold on;
  
  for ww = 1:length(Nmag_per_wave_array)
    
    Nmag_per_wave = Nmag_per_wave_array(ww);
        
    plot(zrange/array_height,squeeze(forces(nn,ww,:,3)),'Tag',num2str(Nmag_per_wave));
    
  end
  
  set(gca,'xlim',[0.009 zrange(end)]/array_height);
  h1 = draworigin([1 0],'v'); set(h1,'linestyle','--');
  colourplot
  labelplot('northeast','vertical','$M$')
  legendshrink
  xlabel('Normalised vertical displacement')
  ylabel('Vertical force, N')
  
end


for ww = 1:length(Nmag_per_wave_array)
  
  Nmag_per_wave = Nmag_per_wave_array(ww);

  figname = ['halbach-discrete-Nmag-',num2str(Nmag_per_wave)];
  willfig(figname); clf; hold on;
  
  for nn = 1:length(Nwaves_array)
    
    Nwaves = Nwaves_array(nn);    
    plot(zrange/array_height,squeeze(forces(nn,ww,:,3)),'Tag',num2str(Nwaves));
    
  end
  
  set(gca,'xlim',[0.009 zrange(end)]/array_height);
  h1 = draworigin([1 0],'v'); set(h1,'linestyle','--');
  colourplot
  labelplot('northeast','vertical','$N$')
  legendshrink
  xlabel('Normalised vertical displacement')
  ylabel('Vertical force, N')
  
end
