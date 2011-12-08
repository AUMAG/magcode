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

%% Calculations for cuboid magnets


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
forces = nan([3 Nm Ny]);
  
for mm = 1:Nm
  a = (volume*magratio(mm)^2)^(1/3);
  b = (volume/magratio(mm))^(1/3);
  
  for yy = 1:Ny
    forces(:,mm,yy) = magnetforces(mag(a,b,1),mag(a,b,-1), [0;0; a+yrange(yy)] );
  end
end


%% Force versus displacement figures

willfig('mag ratios1','small'); clf; hold on

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


willfig('mag ratios2','small'); clf; hold on
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
forces = nan([3 Nm Ny]);
  
for mm = 1:Nm
  a = (volume*magratio(mm)^2)^(1/3);
  b = (volume/magratio(mm))^(1/3);
  for yy = 1:Ny
    forces(:,mm,yy) = magnetforces(mag(a,b,1),mag(a,b,-1), [0;0; a+yrange(yy)] );
  end
end


willfig('mag ratios 3','small'); clf; hold on

for dd = 1:Ny %#ok<FORPF> :: plots should be in a certain order
  
  plot(magratio,squeeze(forces(3,:,dd)))
  [yy, ii] = max(squeeze(forces(3,:,dd)));
  plot(magratio(ii),yy,'.','MarkerSize',8);
    
end

xlim([0.1 1])
ylim([0 30])
xlabel 'Magnet size ratio $\mbqmagratio$'    'interpreter' 'none'
ylabel 'Force, N'                            'interpreter' 'none'

set(gca,'xtick',0:0.2:1)

H = annotation('textarrow',[0.77 0.77],[0.75 0.3]);
set(H,...
  'HeadLength',6.5,'HeadWidth',4.5,'HeadStyle','cback3',...
  'String',{'Increasing','displacement'},'HorizontalAlignment','Center')


if ~simple_graph
  colourplot
  matlabfrag('fig/magratios')
end



%% And the same thing for cylindrical magnets
%
% Here the volume is fixed at
%
% $$ V = \pi r^2 a $$
%
% with magnet ratio
%
% $$ \alpha = r/a $$
%
% So, the radius and side length, resp., are
%
% $$ r = \left( \frac{\alpha V}{\pi} \right)^{1/3} $$
%
% $$ a = \left( \frac{V}{\pi \alpha^2} \right)^{1/3} $$

m = 0.01;
volume = m^3;

cylmag = @(r,a,magdir) ...
  struct('dim',[r a],'magn',1,'magdir',[0 0 magdir]);

% magnet ratio range
cylmagratio = [0.6:0.1:1.0, 1.2:0.2:2.0];
Ncm = length(cylmagratio);
Ncrit = find(cylmagratio==1.0);

% displacement range
ymax = 1.5*m;
ymin = 0.0001;
Ny = 50;
yrange = linspace(ymin,ymax,Ny);

% initialise variable
cylforces = nan([3 Ncm Ny]);

for mm = 1:Ncm
  a = (volume/(pi*cylmagratio(mm)^2))^(1/3);
  r = (volume*cylmagratio(mm)/pi)^(1/3);
  
  for yy = 1:Ny
    cylforces(:,mm,yy) = magnetforces(cylmag(r,a,1),cylmag(r,a,-1), [0;0; a+yrange(yy)] );
  end
end

willfig('cylmag ratios1');
clf; hold on

yy = 8;
yyind = find(round(1000*yrange)==yy,1,'first');

normforce = squeeze(cylforces(3,cylmagratio==1.0,:));

for mm = 1:Ncrit-1
  plot(1000*yrange,squeeze(cylforces(3,mm,:))./normforce,'tag',num2str(cylmagratio(mm)))
  text(yy,squeeze(cylforces(3,mm,yyind))./normforce(yyind),num2str(cylmagratio(mm)) )
end

xlabel 'Displacement $x$, mm'
ylabel 'Normalised force $\bar F$' 'interpreter' 'none'

if ~simple_graph
  axistight
  colourplot
  matlabfrag('fig/magratio1-cyl')
end

willfig('cylmag ratios2');
clf; hold on

yy = 1;
yyind = find(round(1000*yrange)==yy,1,'first');

normforce = squeeze(cylforces(3,cylmagratio==1.0,:));

for mm = Ncrit+1:Ncm
  plot(1000*yrange,squeeze(cylforces(3,mm,:))./normforce,'tag',num2str(cylmagratio(mm)))
  text(yy,squeeze(cylforces(3,mm,yyind))./normforce(yyind),num2str(cylmagratio(mm)) )
end

xlabel 'Displacement $x$, mm'
ylabel 'Normalised force $\bar F$' 'interpreter' 'none'

if ~simple_graph
  axistight
  colourplot
  matlabfrag('fig/magratio2-cyl')
end


%% Recalculate cylindrical case for force versus magnet ratio figure

yrange = 0.001:0.001:m;
Ny = length(yrange);

% sizes
Nm = 100;
cylmagratio = linspace(0.1,2.5,Nm);

% initialise variables
cylforces = nan([3 Nm Ny]);
 
for mm = 1:Nm
  a = (volume/(pi*cylmagratio(mm)^2))^(1/3);
  r = (volume*cylmagratio(mm)/pi)^(1/3);
  for yy = 1:Ny
    cylforces(:,mm,yy) = magnetforces(cylmag(r,a,1),cylmag(r,a,-1), [0;0; a+yrange(yy)] );
  end
end

willfig('cylmag ratios 3','small'); clf; hold on

for dd = 1:Ny %#ok<FORPF> :: plots should be in a certain order

  plot(cylmagratio,squeeze(cylforces(3,:,dd)))
  [yy, ii] = max(squeeze(cylforces(3,:,dd)));
  plot(cylmagratio(ii),yy,'.','MarkerSize',8);

end

xlim([0 2.5])
xlabel 'Magnet size ratio $\ratioMag$'       'interpreter' 'none'
ylabel 'Force, N'                            'interpreter' 'none'

H = annotation('textarrow',[0.77 0.77],[0.75 0.3]);
set(H,...
  'HeadLength',6.5,'HeadWidth',4.5,'HeadStyle','cback3',...
  'String',{'Increasing','displacement'},'HorizontalAlignment','Center')


if ~simple_graph
  colourplot
  matlabfrag('fig/magratios-cyl')
end



