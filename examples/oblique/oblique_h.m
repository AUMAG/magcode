%% Oblique magnets for low stiffness: vary angle; 3 DOF
%
% First investigation into 3 DOF behaviour while varying magnet angle.
% (|oblique_gaps_h.m| is the same but varying magnet gap.)
% Stiffnesses are analysed by calculated the numerical gradient of the
% forces in each direction.
%
% Many graphs are produced; few are used for publication. A graph is
% produced demonstrating the positive stiffness can be achieved in two
% orthogonal directions simultaneously.

%% setup

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

%% constants and calculations

m = 0.01;
angles = 0:5:90;
N_ngl = length(angles);
dd = 0.0001;

calc_f = @(offset) oblique_forces(...
    'magn',1,...
    'unitlength',m,...
    'magratio',0.4,...
    'dispratio',1,...
    'points',200,...
    'magangle',angles,...
    'gapratio',0.2, ...
    'dispoffset',offset...
    );

  
if isempty(mfilename)
  error('This code chunk must be executed as an m-file')
end

datafile = ['data/',mfilename,'.mat'];
if exist(datafile,'file')
  load(datafile)
else
  [hk_yrange1 hk_forces]   = calc_f([  0; 0;   0]);
  [hk_yrange2 hk_forcesX2] = calc_f([ dd; 0;   0]);
  [hk_yrange3 hk_forcesZ2] = calc_f([  0; 0;  dd]);
  hk_yrange = hk_yrange1;
  save(datafile,...
    'hk_yrange','hk_forces','hk_forcesX2','hk_forcesZ2'...
  )
end

%% calculate stiffnesses

yy = hk_yrange(1,:);
yys = yy;
dy = yy(2)-yy(1);

hk_stiffness_Y = -squeeze((hk_forces(2,:,[3,3:end,end])-hk_forces(2,:,[1,1:end-2,end-2]))/(2*dy));
hk_stiffness_X = -squeeze(hk_forcesX2(1,:,:)/dd);
hk_stiffness_Z = -squeeze(hk_forcesZ2(3,:,:)/dd);

yys([1 end]) = [];
hk_stiffness_X(:,[1 end]) = [];
hk_stiffness_Y(:,[1 end]) = [];
hk_stiffness_Z(:,[1 end]) = [];


%% X-offset forces
%
% Forces in xyz for offset in X

willfig('x force x'); clf; hold on

for mm = 1:N_ngl
  plot(yy,squeeze(hk_forcesX2(1,mm,:)))
end
colourplot;

willfig('x force y'); clf; hold on

for mm = 1:N_ngl
  plot(yy,squeeze(hk_forcesX2(2,mm,:)))
end
colourplot;

if ~all(hk_forces(3,mm,:)<1e-10);
  disp('WARNING: z forces not zero')
  willfig('x force z'); clf; hold on
  for mm = 1:N_ngl
    plot(yy,squeeze(hk_forcesX2(3,mm,:)))
  end
  colourplot;
end



%% XYZ stiffnesses
%
% These are the stiffnesses in each direction
% due to vertical displacement

willfig('x stiffness2'); clf; hold on

for mm = 1:N_ngl
  plot(yys,hk_stiffness_X(mm,:),'tag','X')
end
draworigin;
colourplot;

willfig('y stiffness2'); clf; hold on

for mm = 1:N_ngl
  plot(yys,hk_stiffness_Y(mm,:),'tag','Y')
end
draworigin;
colourplot;

willfig('z stiffness2'); clf; hold on

for mm = 1:N_ngl
  plot(yys,hk_stiffness_Z(mm,:),'tag','Z')
end
draworigin;
colourplot;

%%

willfig('compare stiffness');
figuresize(14,6,'centimeters')

% for colours only:
for mm = 1:N_ngl
  plot(1000*yys,1000\hk_stiffness_X(mm,:))
end
colours = colourplot;

subplot(1,2,1); clf; hold on
subplot(1,2,2); clf; hold on

for mm = 1:N_ngl
  
  pospos = hk_stiffness_Y(mm,:)>0 & hk_stiffness_X(mm,:)>0;
  
  subplot(1,2,1); hold on
  plot(1000*yys(pospos),1000\hk_stiffness_X(mm,pospos),'color',colours(mm,:))
  hkXneg = hk_stiffness_X(mm,:);
  hkXneg(pospos) = NaN;
  plot(1000*yys,1000\hkXneg,'-','color',0.85*[1 1 1],'userdata','colourplot:ignore')
  
  subplot(1,2,2); hold on
  plot(1000*yys(pospos),1000\hk_stiffness_Y(mm,pospos),'color',colours(mm,:))
  hkYneg = hk_stiffness_Y(mm,:);
  hkYneg(pospos) = NaN;
  plot(1000*yys,1000\hkYneg,'-','color',0.85*[1 1 1],'userdata','colourplot:ignore')
  
end

subplot(1,2,1);
xlim([0 9.9])
ylim([-1 7.9])
xlabel('Displ.\ $\mbqvdisp$, mm')
set(gca,'xtick',0:2:10)
ylabel('Stiffness, kN/m')
draworigin;
annotation('arrow',[0.2 0.3],[0.7 0.7],arrowsetup{:});
text(6,7,'Vertical')

subplot(1,2,2);
xlim([0 9.9])
ylim([-1 7.9])
xlabel('Displ.\ $\mbqvdisp$, mm')
set(gca,'xtick',0:2:10)
draworigin;
annotation('arrow',[0.65 0.8],[0.22 0.22],arrowsetup{:})
text(5,7,'Horizontal')

matlabfrag('fig/mbq-kvxy-gaps')

%%

willfig('stiffness ratio'); clf; hold on

hk_k_ratio = nan(size(hk_stiffness_Y));

for mm = 1:N_ngl
  pospos = hk_stiffness_Y(mm,:)>0 & hk_stiffness_X(mm,:)>0;
  hk_k_ratio(mm,pospos) = hk_stiffness_Y(mm,pospos)./hk_stiffness_X(mm,pospos);
  plot(yys,hk_k_ratio(mm,:))
end
colourplot;
ylim([0 10])

%%

willfig('load graph'); clf; hold on

nicekratio = inf;
nicek = hk_k_ratio < nicekratio & hk_k_ratio > 1/nicekratio;

for mm = 1:N_ngl
  
  ff = squeeze(hk_forces(2,mm,2:end-1));
  kk = squeeze(hk_stiffness_Y(mm,:)');
  kk(kk<0) = NaN;
  
  ww = sqrt(kk./(ff/9.81)); % resonance frequency
  ww_nice = ww;
  ww_nice(~nicek(mm,:)) = NaN;
  
  plot(ff,ww_nice/(2*pi))
  
  kk_not_nice = kk;
  kk_not_nice(nicek(mm,:)) = NaN;
  ww_not_nice = sqrt(kk_not_nice./(ff/9.81)); % resonance frequency
  
  plot(ff,ww_not_nice/(2*pi),'--','color',0.8*[1 1 1],'userdata','colourplot:ignore')

end

colourplot;
xlim([5 35])
ylim([0 8])
xlabel('Load force, N')
ylabel('Natural frequency, Hz')





