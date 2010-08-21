%% Oblique magnets for low stiffness
%
% Explanation forthcoming
%

clear all
clc

timestamp(mfilename)

%% Constants

magn = 1;

% Magnet setup

mag = @(a,b,magdir) ...
  struct('dim',[b b a],'magn',magn,'magdir',[0 0 magdir]);


%% Calculations

m = 0.01;
volume = m^3;

% sizes
magratio = 0.1:0.1:1;
Nm = length(magratio);

% height ranges
ymax = 2*m;
Ny = 25;


% initialise variables
forces = repmat(NaN,[3 Nm Ny]);

ymin = 0.0001;
yrange = linspace(ymin,ymin+ymax,Ny);
  
for mm = 1:Nm
  
  a = (volume*magratio(mm)^2)^(1/3);
  b = (volume/magratio(mm))^(1/3);
  
  for yy = 1:Ny
    
    forces(:,mm,yy) = magnetforces(mag(a,b,1),mag(a,b,-1), [0;0; a+yrange(yy)] );
    
  end
end

%%

willfig('mag ratios'); clf; hold on

for mm = 1:Nm
  
  plot(yrange,squeeze(forces(3,mm,:)))
    
end

colourplot
axis tight


%%

willfig('mag ratios 2'); clf; hold on

ind = magratio==1;

for mm = 1:Nm
  
  ratio_array = repmat(magratio(mm),[1 length(yrange)]);
  plot(yrange,squeeze(forces(3,mm,:))./squeeze(forces(3,ind,:)))
    
end

colourplot
axis tight


%% Calculations

% sizes
Nm = 100;
magratio = linspace(0.1,1,Nm);

% initialise variables
forces = repmat(NaN,[3 Nm Ny]);
yrange = 0.001:0.001:m;
Ny = length(yrange);
  
for mm = 1:Nm
  a = (volume*magratio(mm)^2)^(1/3);
  b = (volume/magratio(mm))^(1/3);
  for yy = 1:Ny
    forces(:,mm,yy) = magnetforces(mag(a,b,1),mag(a,b,-1), [0;0; a+yrange(yy)] );
  end
end

willfig('mag ratios 3'); clf; hold on

for dd = 1:Ny
  
  plot(magratio,squeeze(forces(3,:,dd)))
  [yy ii] = max(squeeze(forces(3,:,dd)));
  plot(magratio(ii),yy,'.');
    
end

xlim([0.1 1])
colourplot
xlabel 'Magnet size ratio $\mbqmagratio$'
ylabel 'Force, N'

set(gca,'xtick',0:0.25:1)

H = annotation('textarrow',[0.77 0.77],[0.67 0.3]);
set(H,'HeadStyle','cback3', 'String',{'Increasing','displacement'},'HorizontalAlignment','Center')

matlabfrag('fig/mbq-magratios')

