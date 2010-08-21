%% Oblique magnets for low stiffness
%
% Explanation forthcoming
%

close all
clc

timestamp(mfilename)

m = 0.01;
gaps = 0:0.1:1;
N_gaps = length(gaps);

datafile = [mfilename,'.mat'];
if exist(datafile,'file')
  load(datafile)
else
  [gap_yrange gap_forces] = oblique_forces(...
    'magn',1,...
    'unitlength',m,...
    'magratio',0.5,...
    'dispratio',1,...
    'points',200,...
    'magangle',45,...
    'gapratio',gaps ...
    );
  save(datafile,'gap_yrange','gap_forces')
end

%% forces

willfig('forces by gap'); clf; hold on

for tt = 1:N_gaps
  ff = squeeze(gap_forces(2,tt,:));
  yy = squeeze(gap_yrange(tt,:)-gap_yrange(tt,1));
  [Fmax,iFmax]=max(ff);
  
  if iFmax ~= 1
    plot(1000*yy(iFmax),Fmax,'.')
    plot(1000*yy(1:iFmax),ff(1:iFmax),'--','color',0.7*[1 1 1],'UserData','colourplot:ignore')
    plot(1000*yy(iFmax:end),ff(iFmax:end))
  else
    plot(1000*yy,ff)
  end
end

colourplot
ylim([0 25])
set(gca,'xtick',0:2.5:10)
xlabel('Displacement, mm')
ylabel('Force, N')

H = annotation('textarrow',[0.8 0.75],[0.52 0.2]);
set(H,'HeadStyle','cback3','String',{'Increasing','gap'})

matlabfrag('fig/mbq-fvx-gaps')



%% stiffnesses


willfig('stiffnesses by gap'); clf; hold on

for tt = 1:N_gaps
  
  ff = squeeze(gap_forces(2,tt,:));
  yy = 1000*squeeze(gap_yrange(tt,:)-gap_yrange(tt,1));
  kk = -gradient(ff,yy(2)-yy(1));
  [Fmax,iFmax]=max(ff);
  
  if iFmax ~= 1
    plot(yy(1:iFmax),kk(1:iFmax),'--','color',0.7*[1 1 1],'UserData','colourplot:ignore')
    plot(yy(iFmax:end),kk(iFmax:end))
  else
    plot(yy,kk)
  end
end

draworigin([0 0],'--')
colourplot
xlabel('Displacement, mm')
ylabel('Stiffness, kN/m')
ylim([-0.5 5])
set(gca,'xtick',0:2.5:10)

H = annotation('textarrow',[0.8 0.75],[0.52 0.25]);
set(H,'HeadStyle','cback3','String',{'Increasing','gap'})

matlabfrag('fig/mbq-kvx-gaps')


%% stiffness v force

willfig('stiffness v force, by gap'); clf; hold on

for tt = 1:N_gaps
  
  ff = squeeze(gap_forces(2,tt,:));
  yy = squeeze(gap_yrange(tt,:)-gap_yrange(tt,1));
  kk = -gradient(ff,yy(2)-yy(1));
  
  plot(ff,kk/1000)
  
end

ylim([0 6])
colourplot

xlabel('Force, N')
ylabel('Stiffness, kN/m')
H = annotation('textarrow',[0.77 0.5],[0.77 0.41]);
set(H,'String',{'Increasing','gap'},'HorizontalAlignment','Left')

matlabfrag('fig/mbq-kvf-gaps')
