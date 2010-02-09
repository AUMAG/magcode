%% Multipole array forces
%
% Calculating the force between linear Halbach magnet arrays
% in a range of configurations.
%
% These calculations can take quite a while to perform.



%% Calculate

Nmag_per_wave_array = [2 4 8];
Nwaves_array = [1 2 4];
  
datafile = 'data/multipole_compare_data.mat';
if exist(datafile,'file')
  load(datafile);
else
  
  displ_steps = 50;
  zrange = linspace(0.0101,0.02,displ_steps);
  displ = [0; 0; 1]*zrange;
  
  array_height = 0.01;
  array_length = 0.1;
  
  forces = repmat(NaN,...
    [3 displ_steps length(Nwaves_array) length(Nmag_per_wave_array)]);
  
  
  for nn = 1:length(Nwaves_array)
    
    Nwaves = Nwaves_array(nn);
    
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
        'Nwaves', Nwaves,   ...
        'magn', 1           ...
        );
      
      float_array = fixed_array;
      float_array.face = 'down';
      
      forces(:,:,nn,ww) = multipoleforces(fixed_array, float_array, displ);
      
    end
    
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
  
  save(datafile,'zrange','one_mag_force','Nmag_per_wave_array','Nwaves_array','array_height','array_length','forces');

end

%% Plot

zdata = zrange/array_height;
th = [0.5 1 1.5];

for nn = 1:length(Nwaves_array)
  
  Nwaves = Nwaves_array(nn);
  
  figname = ['halbach-waves-Nwaves-',num2str(Nwaves)];
  willfig(figname,'small'); clf; hold on;
  
  for ww = length(Nmag_per_wave_array):-1:1
    
    Nmag_per_wave = Nmag_per_wave_array(ww);
        
    plot(zdata,squeeze(forces(3,:,nn,ww)),...
      'LineWidth',th(ww),...
      'Tag',num2str(Nmag_per_wave));
    
  end
  
  plot(zdata,one_mag_force(3,:),'--k','UserData','colourplot:ignore');

  set(gca,'xlim',[0.9*zdata(1) zdata(end)],...
    'xtick', [1:0.25:2] );
  set(gca,'ylim',[0 550])
  
  set(gca,'box','on','ticklength',[0.02 0.05])
  
  draworigin([1 0],'v',':')
  colourplot(1,[3 2 1])
  
  xlabel('Normalised vertical displacement')
  title(['$\mupNwaves=',num2str(Nwaves),'$'])
  
  if nn == 3
    H = labelplot('north','vertical','$\mupmagperwave$');
    pos = get(H,'position');
    set(H,'position',[0.6 0.5 pos(3:4)]);
    legendshrink
  end
  
  if nn == 1
    ylabel('Vertical force, N')
    pos2 = get(gca,'position');
  else
    set(gca,'yticklabel',[]);
    set(gca,'position',pos2);
  end
  
  matlabfrag(['fig/',figname],'dpi',3200);
  
end

%%

for ww = 1:length(Nmag_per_wave_array)
  
  Nmag_per_wave = Nmag_per_wave_array(ww);

  figname = ['halbach-discrete-Nmag-',num2str(Nmag_per_wave)];
  willfig(figname,'small'); clf; hold on;
  
  for nn = length(Nwaves_array):-1:1
    
    Nwaves = Nwaves_array(nn);    
    plot(zdata,squeeze(forces(3,:,nn,ww)),...
      'LineWidth',th(nn),...    
      'Tag',num2str(Nwaves));
    
  end
  
  plot(zdata,one_mag_force(3,:),'--k','UserData','colourplot:ignore');

  set(gca,'xlim',[0.9*zdata(1) zdata(end)],'xtick',[1:0.25:2]);
  set(gca,'ylim',[0 550])
  
  set(gca,'box','on','ticklength',[0.02 0.05])
  
  draworigin([1 0],'v',':')
  colourplot(1,[3 2 1])
  
  xlabel('Normalised vertical displacement')
  title(['$\mupmagperwave=',num2str(Nmag_per_wave),'$'])
  
  if ww == 3
    H = labelplot('north','vertical','$\mupNwaves$');
    pos = get(H,'position');
    set(H,'position',[0.6 0.5 pos(3:4)]);
    legendshrink
  end
  
  if ww == 1
    ylabel('Vertical force, N')
    pos2 = get(gca,'position');
  else
    set(gca,'yticklabel',[]);
    set(gca,'position',pos2);
  end
  
  matlabfrag(['fig/',figname],'dpi',3200);
  
end

%%

figname = 'halbach-discrete-8v4-Nmag';
willfig(figname,'small'); clf; hold on;

for nn = 1:length(Nwaves_array)
  
  plot(zdata,squeeze(forces(3,:,nn,3)./forces(3,:,nn,2)))
  
end

axis tight
set(gca,'xlim',[0.9*zdata(1) zdata(end)]);
draworigin([1 1],':')
set(gca,'box','on','ticklength',[0.02 0.05])

s = '\providecommand\Nwaves{N_w}';
text(1.2, 1.19,[s '$\Nwaves=4$'],'interpreter','LaTeX')
text(1.5, 1.12,[s '$\Nwaves=2$'],'interpreter','LaTeX')
text(1.7, 1.05,[s '$\Nwaves=1$'],'interpreter','LaTeX')

colourplot

xlabel('Normalised vertical displacement')
ylabel('Forces ratio')


figname = 'halbach-discrete-4v2-Nmag';
willfig(figname,'small'); clf; hold on;

for nn = 1:length(Nwaves_array)
  
  plot(zdata,squeeze(forces(3,:,nn,2)./forces(3,:,nn,1)))
  
end

axis tight
set(gca,'xlim',[0.9*zdata(1) zdata(end)]);
draworigin([1 1],':')
set(gca,'box','on','ticklength',[0.02 0.05])

s = '\providecommand\Nwaves{N_w}';
text(1.2, 1.53,[s '$\Nwaves=4$'],'interpreter','LaTeX')
text(1.5, 1.33,[s '$\Nwaves=2$'],'interpreter','LaTeX')
text(1.7, 1.07,[s '$\Nwaves=1$'],'interpreter','LaTeX')

colourplot

xlabel('Normalised vertical displacement')
ylabel('Forces ratio')

