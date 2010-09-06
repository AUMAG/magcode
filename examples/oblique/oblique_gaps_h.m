%% Oblique magnets for low stiffness
%
% Explanation forthcoming
%

%%

close all
clc
timestamp(mfilename)

if isempty(mfilename)
  error('This code chunk must be executed as an m-file')
end

m = 0.01;
gaps = 0.05:0.05:0.5;
N_gaps = length(gaps);
magangle = 45;

dd = 0.0001;

calc_f = @(offset) oblique_forces(...
    'magn',1,...
    'unitlength',m,...
    'magratio',0.4,...
    'dispratio',1,...
    'points',200,...
    'magangle',magangle,...
    'gapratio',gaps, ...
    'dispoffset',offset...
    );

datafile = [mfilename,'.mat'];
if exist(datafile,'file')
  load(datafile)
else
  [hkg_yrange1 hkg_forces]   = calc_f([  0; 0;   0]);
  [hkg_yrange2 hkg_forcesX2] = calc_f([ dd; 0;   0]);
  [hkg_yrange3 hkg_forcesZ2] = calc_f([  0; 0;  dd]);
  hkg_yrange = hkg_yrange1;
  save(datafile,'magangle',...
    'hkg_yrange','hkg_forces','hkg_forcesX2','hkg_forcesZ2'...
  )
end

%%

yy = hkg_yrange(1,:);
yys = yy;
dy = yy(2)-yy(1);

hkg_stiffness_Y = -squeeze((hkg_forces(2,:,[3,3:end,end])-hkg_forces(2,:,[1,1:end-2,end-2]))/(2*dy));
hkg_stiffness_X = -squeeze(hkg_forcesX2(1,:,:)/dd);
hkg_stiffness_Z = -squeeze(hkg_forcesZ2(3,:,:)/dd);

yys([1 end]) = [];
hkg_stiffness_X(:,[1 end]) = [];
hkg_stiffness_Y(:,[1 end]) = [];
hkg_stiffness_Z(:,[1 end]) = [];



%% XYZ plain forces
%
% These are the forces in each direction without any offset

if ~all(hkg_forces(1,:,:)<1e-10);
  disp('WARNING: x forces not zero')  
  willfig('force x'); clf; hold on
  for mm = 1:N_gaps
    plot(yy,squeeze(hkg_forces(1,mm,:)))
  end
  colourplot
end

if false % shown below
  willfig('force y'); clf; hold on
  for mm = 1:N_gaps
    plot(yy,squeeze(hkg_forces(2,mm,:)))
  end
  colourplot
end

if ~all(hkg_forces(3,:,:)<1e-10);
  disp('WARNING: z forces not zero')
  willfig('force z'); clf; hold on
  for mm = 1:N_gaps
    plot(yy,squeeze(hkg_forces(3,mm,:)))
  end
  colourplot
end


%% X-offset forces
%
% Forces in xyz for offset in X

willfig('x force x'); clf; hold on

for mm = 1:N_gaps
  plot(yy,squeeze(hkg_forcesX2(1,mm,:)))
end
colourplot

willfig('x force y'); clf; hold on

for mm = 1:N_gaps
  plot(yy,squeeze(hkg_forcesX2(2,mm,:)))
end
colourplot

if ~all(hkg_forces(3,mm,:)<1e-10);
  disp('WARNING: z forces not zero')
  willfig('x force z'); clf; hold on
  for mm = 1:N_gaps
    plot(yy,squeeze(hkg_forcesX2(3,mm,:)))
  end
  colourplot
end



%% XYZ stiffnesses
%
% These are the stiffnesses in each direction
% due to vertical displacement

willfig('x stiffness2'); clf; hold on

for mm = 1:N_gaps
  plot(yys,hkg_stiffness_X(mm,:),'tag','X')
end
draworigin
colourplot

willfig('y stiffness2'); clf; hold on

for mm = 1:N_gaps
  plot(yys,hkg_stiffness_Y(mm,:),'tag','Y')
end
draworigin
colourplot

willfig('z stiffness2'); clf; hold on

for mm = 1:N_gaps
  plot(yys,hkg_stiffness_Z(mm,:),'tag','Z')
end
draworigin
colourplot

%%

willfig('compare stiffness');
subplot(1,2,1); clf; hold on
subplot(1,2,2); clf; hold on

for mm = 1:N_gaps
  
  pospos = hkg_stiffness_Y(mm,:)>0 & hkg_stiffness_X(mm,:)>0;
  
  subplot(1,2,1); hold on
  plot(yys( pospos),1000\hkg_stiffness_X(mm, pospos))
  hkXneg = hkg_stiffness_X(mm,:);
  hkXneg(pospos) = NaN;
  plot(yys,1000\hkXneg,'--','color',0.7*[1 1 1],'userdata','colourplot:ignore')
  
  subplot(1,2,2); hold on
  plot(yys( pospos),1000\hkg_stiffness_Y(mm, pospos))
  hkYneg = hkg_stiffness_Y(mm,:);
  hkYneg(pospos) = NaN;
  plot(yys,1000\hkYneg,'--','color',0.7*[1 1 1],'userdata','colourplot:ignore')
  
end

subplot(1,2,1);
axis tight
colourplot
ylim([-1 10])

subplot(1,2,2);
axis tight
colourplot
title('Y')
ylim([-1 10])

%%

willfig('stiffness ratio'); clf; hold on

hkg_k_ratio = nan(size(hkg_stiffness_Y));

for mm = 1:N_gaps
  pospos = hkg_stiffness_Y(mm,:)>0 & hkg_stiffness_X(mm,:)>0;
  hkg_k_ratio(mm,pospos) = hkg_stiffness_Y(mm,pospos)./hkg_stiffness_X(mm,pospos);
  plot(yys,hkg_k_ratio(mm,:))
end
colourplot
ylim([0 10])

%%

willfig(['load graph gaps at angle ',num2str(magangle)]); clf; hold on

nicekratio = inf;
nicek = hkg_k_ratio < nicekratio & hkg_k_ratio > 1/nicekratio;

for mm = 1:N_gaps
  
  ff = squeeze(hkg_forces(2,mm,2:end-1));
  
  kk = squeeze(hkg_stiffness_Y(mm,:)');
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

colourplot
xlim([5 40])
ylim([0 10])
xlabel('Load force, N')
ylabel('Natural frequency, Hz')





