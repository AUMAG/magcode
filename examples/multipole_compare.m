%% Multipole array forces
%
% Calculating the force between multipole magnet arrays
% in a range of configurations.

%%

clc
timestamp
tic

%% Common setup

displ_steps = 5;
zrange = linspace(0.0101,0.03,displ_steps);

Nmag_per_wave_array = [1 2 4 8];
Nwaves_array = [1 2 4];
height = 0.01;

forces = repmat(NaN,...
  [length(Nwaves_array) length(Nmag_per_wave_array) displ_steps 3]);

%% Calculate 

for nn = 1:length(Nwaves_array)
  
  Nwaves = Nwaves_array(nn);
    
  for ww = 1:length(Nmag_per_wave_array)
    
    Nmag_per_wave = Nmag_per_wave_array(ww);
    
    fixed_array = ...
      struct(...
      'type','linear-x',  ...
      'face','up',        ...
      'length', 0.10,     ...
      'width',  0.01,     ...
      'height', height,   ...
      'Nmag_per_wave', Nmag_per_wave, ...
      'Nwaves', Nwaves,   ...
      'magn', 1,          ...
      'magdir_first', 90  ...
      );
    
    float_array = fixed_array;
    float_array.face = 'down';
    float_array.magdir_first = -90;
        
    for ii = 1:N
      displ = [0 0 zrange(ii)];
      forces(nn,ww,ii,:) = multipoleforces(fixed_array, float_array, displ);
    end
        
  end
    
end

toc

%% Plot

for nn = 1:length(Nwaves_array)
  
  Nwaves = Nwaves_array(nn);
  
  figname = ['halbach-waves-Nwaves-',num2str(Nwaves)];
  willfig(figname); clf; hold on;
  
  for ww = 1:length(Nmag_per_wave_array)
    
    Nmag_per_wave = Nmag_per_wave_array(ww);
        
    plot(zrange/height,squeeze(forces(nn,ww,:,3)),'Tag',num2str(Nmag_per_wave));
    
  end
  
  set(gca,'xlim',[0.009 zrange(end)]);
  h1 = draworigin([0.01 0],'v');
  set(h1,'linestyle','--');
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
    plot(zrange,squeeze(forces(nn,ww,:,3)),'Tag',num2str(Nwaves));
    
  end
  
  set(gca,'xlim',[0.009 zrange(end)]);
  h1 = draworigin([0.01 0],'v');
  set(h1,'linestyle','--');
  colourplot
  labelplot('northeast','vertical','$N$')
  legendshrink
  xlabel('Normalised vertical displacement')
  ylabel('Vertical force, N')
  
end
