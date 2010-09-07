%% Oblique magnets for low stiffness
%
% Explanation forthcoming
%

close all
clc
arrowsetup = {'HeadLength',6.5,'HeadWidth',4.5,'HeadStyle','cback3'};

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  colourplot = @(varargin) disp('');
  draworigin = @(varargin) disp('');
  matlabfrag = @(varargin) disp('');
end

%%

m = 0.01;
gaps = 0:0.05:0.5;
N_gaps = length(gaps);

if isempty(mfilename)
  error('This code chunk must be executed as an m-file')
end

datafile = [mfilename,'.mat'];
if exist(datafile,'file')
  load(datafile)
else
  [gap_yrange gap_forces] = oblique_forces(...
    'magn',1,...
    'unitlength',m,...
    'magratio',0.4,...
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
    plot(1000*yy(1:iFmax),ff(1:iFmax),'--','color',0.7*[1 1 1],'UserData','colourplot:ignore')
    plot(1000*yy(iFmax:end),ff(iFmax:end))
    plot(1000*yy(iFmax),Fmax,'.','color',[0 0 0],'UserData','colourplot:ignore')
  else
    plot(1000*yy,ff)
  end
end

colourplot;
ylim([0 25])
set(gca,'xtick',0:2.5:10)
xlabel('Displacement, mm')
ylabel('Force, N')

H = annotation('textarrow',[0.8 0.75],[0.52 0.2]);
set(H,'String',{'Increasing','gap'},arrowsetup{:})

matlabfrag('fig/mbq-fvx-gaps')



%% stiffnesses


willfig('stiffnesses by gap'); clf; hold on

for tt = 1:N_gaps
  
  ff = squeeze(gap_forces(2,tt,:));
  yy = squeeze(gap_yrange(tt,:)-gap_yrange(tt,1));
  kk = -gradient(ff,yy(2)-yy(1));
  [Fmax,iFmax]=max(ff);
  
  if iFmax ~= 1
    plot(1000*yy(1:iFmax),1000\kk(1:iFmax),'--','color',0.7*[1 1 1],'UserData','colourplot:ignore')
    plot(1000*yy(iFmax:end),1000\kk(iFmax:end))
  else
    plot(1000*yy,1000\kk)
  end
end

draworigin([0 0],'h','--')
colourplot;
xlabel('Displacement, mm')
ylabel('Stiffness, kN/m')
ylim([-0.5 5])
set(gca,'xtick',0:2.5:10)

H = annotation('textarrow',[0.8 0.75],[0.52 0.25]);
set(H,'String',{'Increasing','gap'},arrowsetup{:})

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
colourplot;

xlabel('Force, N')
ylabel('Stiffness, kN/m')

H = annotation('textarrow',[0.77 0.5],[0.77 0.41]);
set(H,'String',{'Increasing','gap'},'HorizontalAlignment','Left',arrowsetup{:})

matlabfrag('fig/mbq-kvf-gaps')


%% resonance v force

willfig('frequency v force, by gap'); clf; hold on

best_gap = gaps==0.05;

for tt = 1:N_gaps
  
  ff = squeeze(gap_forces(2,tt,:));
  yy = squeeze(gap_yrange(tt,:)-gap_yrange(tt,1));
  kk = -gradient(ff,yy(2)-yy(1));
  kk(kk<0)=NaN;
  
  ww = sqrt(kk./(ff/9.81)); % resonance frequency
  
  plot(ff,ww/(2*pi))
  
  if tt == find(best_gap)
    text(ff(1)+1,ww(1)/(2*pi),['\num{',num2str(gaps(tt)),'}'],'interpreter','none')
    plot(ff(1),ww(1)/(2*pi),'.','linewidth',2,'color',[0 0 0],'UserData','colourplot:ignore')
  end
  
end

colourplot;

xlim([0 50])

xlabel('Load force, N')
ylabel('Natural frequency, Hz')

H = annotation('textarrow',[0.38 0.2],[0.72 0.4]);
set(H,'String',{'Increasing','gap'},'HorizontalAlignment','Left',arrowsetup{:})

matlabfrag('fig/mbq-wvf-gaps')