%% Multipole array forces
%
% Calculating the force between linear Halbach magnet arrays
% in a range of configurations.
%
% These calculations can take quite a while to perform.



%% Calculate

datafile = 'data/multipole_const_Nmag_data.mat';
if exist(datafile,'file')
  load(datafile);
else
  
  displ_steps = 80;
  zrange = [ linspace(0.0101,0.0115,displ_steps/2)...
             linspace(0.0116,0.015,displ_steps/2) ];
  displ = [0; 0; 1]*zrange;
  
  Nmag_per_wave_array = [2 4 8];
  
  array_height = 0.01;
  array_length = 0.1;
  
  forces = repmat(NaN,...
    [3 displ_steps length(Nmag_per_wave_array)]);
  
    
    for ww = 1:length(Nmag_per_wave_array)
      
      Nmag_per_wave = Nmag_per_wave_array(ww);
      
      fixed_array = ...
        struct(...
        'type','linear',    ...
        'align','x',        ...
        'face','up',        ...
        'length', array_length,   ...
        'width',  array_height,   ...
        'height', array_height,   ...
        'Nmag_per_wave', Nmag_per_wave, ...
        'Nmag', 50,   ...
        'magn', 1           ...
        );
      
      float_array = fixed_array;
      float_array.face = 'down';
      
      forces(:,:,ww) = multipoleforces(fixed_array, float_array, displ);
      
    end
  
  % One magnet:
  fixed_array = ...
    struct(...
    'type','linear',    ...
    'align','x',        ...
    'face','up',        ...
    'length', array_length,   ...
    'width',  array_height,   ...
    'height', array_height,   ...
    'Nmag_per_wave', 1, ...
    'Nwaves', 1,        ...
    'magn', 1           ...
    );
  
  float_array = fixed_array;
  float_array.face = 'down';
  
  one_mag_force = multipoleforces(fixed_array, float_array, displ);
  
  save(datafile,'zrange','one_mag_force','Nmag_per_wave_array','array_height','array_length','forces');
  
end

%% Plot

zdata = zrange/array_height;
th = [0.5 1 1.5];

figname = 'halbach-max-Nmag';
willfig(figname,'small'); clf; hold on;

for ww = length(Nmag_per_wave_array):-1:1
  
  Nmag_per_wave = Nmag_per_wave_array(ww);
  
  plot(zdata,squeeze(forces(3,:,ww)),...
    'LineWidth',th(ww),...
    'Tag',num2str(Nmag_per_wave));
  
end

plot(zdata,one_mag_force(3,:),'--k','UserData','colourplot:ignore');

set(gca,'xlim',[0.95*zdata(1) 1.5],...
  'xtick', 1:0.25:2 );

set(gca,'box','on','ticklength',[0.02 0.05])

draworigin([1 0],'v',':')
colourplot(1,[3 2 1])

xlabel('Normalised vertical displacement')

H = labelplot('north','vertical','$\mupmagperwave$');
pos = get(H,'position');
set(H,'position',[0.6 0.5 pos(3:4)]);
legendshrink

ylabel('Vertical force, N')

matlabfrag(['fig/',figname],'dpi',3200);


