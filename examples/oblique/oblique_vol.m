%% Oblique vs volume
%
% Explanatino forthcoming

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

msizes = [0.01:0.005:0.05];
N_msz = length(msizes);

dd = 0.0001;

calc_f = @(offset) oblique_forces(...
    'magn',1,...
    'unitlength',msizes,...
    'magratio',0.4,...
    'dispratio',1,...
    'points',50,...
    'magangle',40,...
    'gapratio',0.15, ...
    'dispoffset',offset...
    );

if isempty(mfilename)
  error('This code chunk must be executed as an m-file')
end

datafile = [mfilename,'.mat'];
if exist(datafile,'file')
  load(datafile)
else
  [vol_yrange1 vol_forces]   = calc_f([  0; 0;   0]);
  [vol_yrange2 vol_forcesX2] = calc_f([ dd; 0;   0]);
  [vol_yrange3 vol_forcesZ2] = calc_f([  0; 0;  dd]);
  vol_yrange = vol_yrange1;
  save(datafile,...
    'vol_yrange','vol_forces','vol_forcesX2','vol_forcesZ2'...
  )
end

%% Calculate stiffnesses

yy = vol_yrange(1,:);
yys = yy;
dy = yy(2)-yy(1);

vol_stiffness_Y = -squeeze((vol_forces(2,:,[3,3:end,end])-vol_forces(2,:,[1,1:end-2,end-2])))./(2*repmat(vol_yrange(:,2)-vol_yrange(:,1),[1 50]));
vol_stiffness_X = -squeeze(vol_forcesX2(1,:,:))./repmat(vol_yrange(:,2)-vol_yrange(:,1),[1 50]);
vol_stiffness_Z = -squeeze(vol_forcesZ2(3,:,:))./repmat(vol_yrange(:,2)-vol_yrange(:,1),[1 50]);

yys([1 end]) = [];
vol_stiffness_X(:,[1 end]) = [];
vol_stiffness_Y(:,[1 end]) = [];
vol_stiffness_Z(:,[1 end]) = [];

vol_k_ratio = nan(size(vol_stiffness_Y));

for tt = 1:N_msz
  pospos = vol_stiffness_Y(tt,:)>0 & vol_stiffness_X(tt,:)>0;
  vol_k_ratio(tt,pospos) = vol_stiffness_Y(tt,pospos)./vol_stiffness_X(tt,pospos);
end

nicekratio = inf;
nicek = vol_k_ratio < nicekratio & vol_k_ratio > 1/nicekratio;

%%


willfig('forces by volume'); clf; hold on

for tt = 1:N_msz
  ff = squeeze(vol_forces(2,tt,:));
  yy = squeeze(vol_yrange(tt,:)-vol_yrange(tt,1));
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
xlabel('Displacement, mm')
ylabel('Force, N')


%%

willfig('frequency v force, by volume'); clf; hold on


for tt = 1:2:N_msz
  
  ff = squeeze(vol_forces(2,tt,2:end-1));
  
  kk = squeeze(vol_stiffness_Y(tt,:))';
  kk(kk<0) = NaN;
  ww = sqrt(kk./(ff/9.81)); % resonance frequency
  ww_nice = ww;
  ww_nice(~nicek(tt,:)) = NaN;
  
  last_stable = find(~isnan(ww_nice),1,'last');
  
  plot(1000\ff,ww_nice/(2*pi))
  plot(1000\ff(last_stable),ww_nice(last_stable)/(2*pi),'k.','userdata','colourplot:ignore')
  
  text(1000\ff(last_stable),ww_nice(last_stable)/(2*pi),...
    [' \SI{', num2str(round(1000*msizes(tt))),'}{mm^3}'],...
      'Interpreter','none',...
      'HorizontalAlignment','left',...
      'VerticalAlignment','baseline')
  
  kk_not_nice = kk;
  kk_not_nice(nicek(tt,:)) = NaN;
  ww_not_nice = sqrt(kk_not_nice./(ff/9.81)); % resonance frequency
  
  plot(1000\ff,ww_not_nice/(2*pi),'-','color',0.85*[1 1 1],'userdata','colourplot:ignore')
  
  
end

colourplot;
xlim([0 0.58])
ylim([0 6.5])
set(gca,'xtick',[0:0.1:0.6])
xlabel('Load force, kN')
ylabel('Natural frequency, Hz')

matlabfrag('fig/mbq-wvf-vol')