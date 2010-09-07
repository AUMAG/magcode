%% Oblique magnets for low stiffness
%
% Explanation forthcoming
%

%%

close all
clc

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  colourplot = @(varargin) disp('');
  draworigin = @(varargin) disp('');
  matlabfrag = @(varargin) disp('');
end

%%

m = 0.01;
gaps = 0.05:0.05:0.5;
N_gaps = length(gaps);
magangles = 30:2:60;
N_ngl = length(magangles);

dd = 0.0001;

calc_f = @(offset) oblique_forces(...
    'magn',1,...
    'unitlength',m,...
    'magratio',0.4,...
    'dispratio',1,...
    'points',1000,...
    'magangle',magangles,...
    'gapratio',gaps, ...
    'dispoffset',offset...
    );

if isempty(mfilename)
  error('This code chunk must be executed as an m-file')
end

datafile = [mfilename,'.mat'];
if exist(datafile,'file')
  load(datafile)
else
  [ga_yrange1 ga_forces]   = calc_f([  0; 0;   0]);
  [ga_yrange2 ga_forcesX2] = calc_f([ dd; 0;   0]);
  [ga_yrange3 ga_forcesZ2] = calc_f([  0; 0;  dd]);
  ga_yrange = ga_yrange1;
  save(datafile,...
    'ga_yrange','ga_forces','ga_forcesX2','ga_forcesZ2'...
  )
end

%%

yy = ga_yrange(1,1,:);
yys = yy;
dy = yy(2)-yy(1);

ga_stiffness_Y = -squeeze((ga_forces(2,:,:,[3,3:end,end])-ga_forces(2,:,:,[1,1:end-2,end-2]))/(2*dy));
ga_stiffness_X = -squeeze(ga_forcesX2(1,:,:,:)/dd);
ga_stiffness_Z = -squeeze(ga_forcesZ2(3,:,:,:)/dd);

yys([1 end]) = [];
ga_stiffness_X(:,:,[1 end]) = [];
ga_stiffness_Y(:,:,[1 end]) = [];
ga_stiffness_Z(:,:,[1 end]) = [];

ga_k_ratio = nan(size(ga_stiffness_Y));

for nn = 1:N_ngl
for mm = 1:N_gaps
  pospos = ga_stiffness_Y(nn,mm,:)>0 & ga_stiffness_X(nn,mm,:)>0;
  ga_k_ratio(nn,mm,pospos) = ga_stiffness_Y(nn,mm,pospos)./ga_stiffness_X(nn,mm,pospos);
end
end

%%

for nn = 1:0%N_ngl
  
  willfig(['load/res at angle ',num2str(magangles(nn))]); clf; hold on
  
  nicekratio = inf;
  nicek = ga_k_ratio < nicekratio & ga_k_ratio > 1/nicekratio;
  
  for mm = 1:N_gaps
    
    ff = squeeze(ga_forces(2,nn,mm,2:end-1));
    
    kk = squeeze(ga_stiffness_Y(nn,mm,:));
    kk(kk<0) = NaN;
    ww = sqrt(kk./(ff/9.81)); % resonance frequency
    ww_nice = ww;
    ww_nice(~nicek(nn,mm,:)) = NaN;
    
    plot(ff,ww_nice/(2*pi))
    
    kk_not_nice = kk;
    kk_not_nice(nicek(nn,mm,:)) = NaN;
    ww_not_nice = sqrt(kk_not_nice./(ff/9.81)); % resonance frequency
    
    plot(ff,ww_not_nice/(2*pi),'-','color',0.8*[1 1 1],'userdata','colourplot:ignore')
    
  end
  
  N_displ = length(ff);
  d_incr = 0.001;
  d_first = find(yy>0.001,1);
  displ_step = round(d_incr/((ga_yrange(1,1,end)-ga_yrange(1,1,1))/N_displ));
  
  for dd = d_first:displ_step:N_displ
    
    ff = squeeze(ga_forces(2,nn,:,dd));
    kk = squeeze(ga_stiffness_Y(nn,:,dd))';
    
    ww = sqrt(kk./(ff/9.81)); % resonance frequency
    
    plot(ff,ww/(2*pi),':','color',0*[1 1 1],'userdata','colourplot:ignore')
    
  end
  
  colourplot;
  xlim([5 50])
  ylim([0 10])
  xlabel('Load force, N')
  ylabel('Natural frequency, Hz')

end

%%

for mm = 1:0%N_gaps
  
  willfig(['load/res at gap ',num2str(gaps(mm))]); clf; hold on
  
  nicekratio = inf;
  nicek = ga_k_ratio < nicekratio & ga_k_ratio > 1/nicekratio;
  
  for nn = 1:N_ngl
    
    ff = squeeze(ga_forces(2,nn,mm,2:end-1));
    
    kk = squeeze(ga_stiffness_Y(nn,mm,:));
    kk(kk<0) = NaN;
    ww = sqrt(kk./(ff/9.81)); % resonance frequency
    ww_nice = ww;
    ww_nice(~nicek(nn,mm,:)) = NaN;
    
    plot(ff,ww_nice/(2*pi))
    
    kk_not_nice = kk;
    kk_not_nice(nicek(nn,mm,:)) = NaN;
    ww_not_nice = sqrt(kk_not_nice./(ff/9.81)); % resonance frequency
    
    plot(ff,ww_not_nice/(2*pi),'-','color',0.8*[1 1 1],'userdata','colourplot:ignore')
    
  end
  
  N_displ = length(ff);
  d_incr = 0.001;
  d_first = find(yy>0.001,1);
  displ_step = round(d_incr/((ga_yrange(1,1,end)-ga_yrange(1,1,1))/N_displ));
  
  for dd = d_first:displ_step:N_displ
    
    ff = squeeze(ga_forces(2,:,mm,dd));
    kk = squeeze(ga_stiffness_Y(:,mm,dd))';
    
    ww = sqrt(kk./(ff/9.81)); % resonance frequency
    
    plot(ff,ww/(2*pi),':','color',0*[1 1 1],'userdata','colourplot:ignore')
    
  end

  colourplot;
  xlim([5 50])
  ylim([0 10])
  xlabel('Load force, N')
  ylabel('Natural frequency, Hz')  

end


%%

nn = find(magangles==40);

willfig(['paper graph @ angle ',num2str(magangles(nn))]); clf; hold on

nicekratio = inf;
nicek = ga_k_ratio < nicekratio & ga_k_ratio > 1/nicekratio;

taggap = [1 2 3];

ff_unstabl = nan([1 N_gaps]);
ww_unstabl = nan([1 N_gaps]);

for mm = 1:N_gaps
  
  ff = squeeze(ga_forces(2,nn,mm,2:end-1));
  
  kk = squeeze(ga_stiffness_Y(nn,mm,:));
  kk(kk<0) = NaN;
  ww = sqrt(kk./(ff/9.81)); % resonance frequency
  ww_nice = ww;
  ww_nice(~nicek(nn,mm,:)) = NaN;
  
  plot(ff,ww_nice/(2*pi))
  
  pp = find(~isnan(ww),1,'first');
  ff_unstabl(mm) = ff(pp);
  ww_unstabl(mm) = ww(pp);
  if any(taggap==mm)
    plot(ff(pp),ww(pp)/2/pi,'k.','userdata','colourplot:ignore')
    text(ff(pp),ww(pp)/2/pi,[num2str(gaps(mm)),'~'],'horizontalalignment','right')
  end
  
  kk = squeeze(ga_stiffness_Y(nn,mm,:));
  kk(kk<0)=NaN;
  kk(nicek(nn,mm,:)) = NaN;
  ww = sqrt(kk./(ff/9.81)); % resonance frequency
  
  plot(ff,ww/(2*pi),'-','color',0.85*[1 1 1],'userdata','colourplot:ignore')
  
end

% plot(ff_unstabl,ww_unstabl/2/pi,'k--')

N_displ = length(ff);
d_incr = 0.001;
d_first = find(yy>0.001,1);
displ_step = round(d_incr/((ga_yrange(1,1,end)-ga_yrange(1,1,1))/N_displ));

lbl = [1 2 3];

incr = 0;

for dd = d_first:displ_step:N_displ
  incr = incr+1;
  
  ff = squeeze(ga_forces(2,nn,:,dd));
  kk = squeeze(ga_stiffness_Y(nn,:,dd))';
  kk(kk<0)=NaN;
  
  ww = sqrt(kk./(ff/9.81)); % resonance frequency
  
  plot(ff,ww/(2*pi),':','color',0*[1 1 1],'LineWidth',1,'userdata','colourplot:ignore')
  
  if any(lbl==incr)
    if incr == 3, halign = 'center'; else halign = 'left'; end
    plot(ff(1),ww(1)/2/pi,'k.','userdata','colourplot:ignore')
    text(ff(1),ww(1)/2/pi,[' ', num2str(round(1000*ga_yrange(1,1,dd))),'\,mm'],...
      'Interpreter','none',...
      'HorizontalAlignment',halign,...
      'VerticalAlignment','baseline')
  end
  
end

colourplot;
xlim([5 37])
ylim([0 6.9])
xlabel('Load force, N')
ylabel('Natural frequency, Hz')

% annotation('arrow',[0.5 0.3],[0.9 0.9])
% annotation('arrow',[0.45 0.25],[0.2 0.2])

matlabfrag('fig/mbq-wvf-angle-stabl')
