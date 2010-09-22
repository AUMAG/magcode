%% Optimising force between two magnets with fixed volume
%

clear all
close all
clc

% In case you don't have the various bits'n'pieces that I use to create
% my Matlab graphics (probably likely).

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  simple_graph = true;
else
  simple_graph = false;
end

%% Calculations


m = 0.01;
volume = m^3;

mag = @(a,b,magdir) ...
  struct('dim',[b b a],'magn',1,'magdir',[0 0 magdir]);

% magnet ratio range
magratio = 0.1:0.1:1;
Nm = length(magratio);

% displacement range
ymax = 1*m;
Ny = 25;
ymin = 0.0001;
yrange = linspace(ymin,ymax,Ny);

% initialise variable
forces = repmat(NaN,[3 Nm Ny]);
  
for mm = 1:Nm
  
  a = (volume*magratio(mm)^2)^(1/3);
  b = (volume/magratio(mm))^(1/3);
  
  for yy = 1:Ny
    
    face_gap = a+yrange(yy);
    forces(:,mm,yy) = magnetforces(mag(a,b,1),mag(a,b,-1), [0;0; face_gap] );
    
  end
end

%% Force versus displacement figures

willfig('mag ratios1'); clf; hold on

for mm = 1:4
  
  plot(1000*yrange,squeeze(forces(3,mm,:))./squeeze(forces(3,end,:)),'tag',num2str(magratio(mm)))
    
end

axis tight
text(8,0.8,'$0.1$')
text(6,1.05,'$0.2$')
text(4,1.13,'$0.3$')
text(1.5,1.22,'$\gamma=0.4$')

ylim([0.7,1.3])

xlabel('Displacement $x$, mm')
ylabel('Normalised force $\bar F$')

if ~simple_graph
  colourplot
  matlabfrag('fig/magratio1')
end


willfig('mag ratios2'); clf; hold on
for mm = Nm-2:-1:4
  
  plot(1000*yrange,squeeze(forces(3,mm,:))./squeeze(forces(3,end,:)),'tag',num2str(magratio(mm)))
    
end

axis tight
text(7.9,1.3,'$\gamma=0.4$','HorizontalAlignment','Right')
text(8,1.255,'$0.5$')
text(8.25,1.21,'$0.6$')
text(8.5,1.155,'$0.7$')
text(8.75,1.1,'$0.8$')

xlabel('Displacement $x$, mm')
ylabel('Normalised force $\bar F$')

if ~simple_graph
  colourplot
  matlabfrag('fig/magratio2')
end


%% Recalculate for force versus magnet ratio figure


yrange = 0.001:0.001:m;
Ny = length(yrange);

% sizes
Nm = 100;
magratio = linspace(0.1,1,Nm);

% initialise variables
forces = repmat(NaN,[3 Nm Ny]);
  
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
xlabel 'Magnet size ratio $\mbqmagratio$'
ylabel 'Force, N'

set(gca,'xtick',0:0.2:1)

H = annotation('textarrow',[0.77 0.77],[0.67 0.3]);
set(H,...
  'HeadLength',6.5,'HeadWidth',4.5,'HeadStyle','cback3',...
  'String',{'Increasing','displacement'},'HorizontalAlignment','Center')


if ~simple_graph
  colourplot
  matlabfrag('fig/magratios')
end

