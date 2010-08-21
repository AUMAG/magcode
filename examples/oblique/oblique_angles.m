%% Oblique magnets for low stiffness
%
% Explanation forthcoming
%

%%

close all
clc
timestamp(mfilename)

m = 0.01;
angles = 0:5:90;
N_ngl = length(angles);

datafile = [mfilename,'.mat'];
if exist(datafile,'file')
  load(datafile)
else
  [ngl_yrange ngl_forces] = oblique_forces(...
    'magn',1,...
    'unitlength',m,...
    'magratio',0.5,...
    'dispratio',1,...
    'points',50,...
    'magangle',angles,...
    'gapratio',0 ...
    );
  save(datafile,'ngl_yrange','ngl_forces')
end

%% forces

willfig('forces by angle'); clf; hold on

for tt = N_ngl:-1:1
  ff = squeeze(ngl_forces(2,tt,:));
  yy = 1000*squeeze(ngl_yrange(tt,:));
  [Fmax,iFmax]=max(ff);
  
  if iFmax ~= 1
    plot(yy(1:iFmax),ff(1:iFmax),'--','color',0.7*[1 1 1],'UserData','colourplot:ignore')
    plot(yy(iFmax:end),ff(iFmax:end))
    plot(yy(iFmax),Fmax,'.','color',[0 0 0],'UserData','colourplot:ignore')
  else
    plot(yy,ff,'tag',num2str(angles(tt)))
  end
end

set(gca,'xtick',0:2.5:10)
xlim([0 10])
ylim([0 50])
xlabel('Displacement, mm')
ylabel('Force, N')
colourplot(1,N_ngl:-1:1)

H = annotation('arrow',[0.8 0.7],[0.6 0.2]);
set(H,'HeadStyle','cback3')
text(8.5,35,{'Increasing','angle'},'HorizontalAlignment','Center')

matlabfrag('fig/mbq-fvx-angle')


%% stiffnesses

willfig('stiffnesses by angle'); clf; hold on

for tt = 1:N_ngl
  
  ff = squeeze(ngl_forces(2,tt,:));
  yy = 1000*squeeze(ngl_yrange(tt,:));
  kk = -gradient(ff,yy(2)-yy(1));
  [Fmax,iFmax]=max(ff);
  
  if iFmax ~= 1
    plot(yy(1:iFmax),kk(1:iFmax),'--','color',0.7*[1 1 1],'UserData','colourplot:ignore')
    plot(yy(iFmax:end),kk(iFmax:end))
  else
    plot(yy,kk)
  end
  
  draworigin([0 0],'--')
  
end

xlim([0 10])
ylim([-5 20])
set(gca,'xtick',0:2.5:10)
xlabel('Displacement, mm')
ylabel('Stiffness, N')
colourplot

H = annotation('arrow',[0.24 0.35],[0.45 0.65]);
set(H,'HeadStyle','cback3')
text(2.8,12.8,{'Increasing','angle'})

matlabfrag('fig/mbq-kvx-angle')

%% stiffness v force

willfig('stiffness v force, by angle'); clf; hold on

for tt = 1:N_ngl
  
  ff = squeeze(ngl_forces(2,tt,:));
  yy = squeeze(ngl_yrange(tt,:));
  kk = -gradient(ff,yy(2)-yy(1));
  
  plot(ff,kk/1000)
  
end

xlim([0 45])
ylim([0 8])
colourplot

xlabel('Force, N')
ylabel('Stiffness, kN/m')

H = annotation('arrow',[0.8 0.5],[0.65 0.8]);
set(H,'HeadStyle','cback3')
text(15,7,{'Increasing','angle'},'HorizontalAlignment','Center')


matlabfrag('fig/mbq-kvf-angle')
