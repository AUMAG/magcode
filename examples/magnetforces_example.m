%% Setup
%
% In case you don't have the various bits'n'pieces that I use to create
% my Matlab graphics (probably likely).

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  simple_graph = 1;
else
  simple_graph = 0;
end


%% Parallel magnets example
%
% Reproducing the results of Akoun & Yonnet (1984)

magnet_fixed.dim = [0.02 0.012 0.006];
magnet_float.dim = [0.012 0.02 0.006];

magnet_fixed.magn = 0.38;
magnet_float.magn = 0.38;

magnet_fixed.magdir = [0  90]; % z
magnet_float.magdir = [0  90]; % z

offset = [-0.004 -0.004 0.008];
N = 1001;
displ = linspace(0, 0.03, N);
f1_xyz = repmat(NaN, [3 N]);

for ii = 1:N
  f1_xyz(:,ii) = magnetforces(magnet_fixed,magnet_float,offset+[displ(ii) 0 0]);
end

willfig('akoun'); clf; hold on
plot(displ,f1_xyz(1,:),'Tag','x')
plot(displ,f1_xyz(2,:),'Tag','y')
plot(displ,f1_xyz(3,:),'Tag','z')
xlabel('Horizontal $x$ displacement, m')
ylabel('Forces, N')
set(gca,'box','on');
set(gca,'ticklength',[0.02 0.05])
text(0.004,-0.5,'$F_x$')
text(0.004, 0.8,'$F_y$')
text(0.004,-1.7,'$F_z$')

if ~simple_graph
  h1 = draworigin([0 0],'h');
  set(h1,'linestyle','--');
  colourplot
  % labelplot
  matlabfrag('fig/akoun-repro','dpi',3200);
end

%% Orthogonal
%
% Replicating the results of Janssen et al. (2009)

magnet_fixed.dim = [0.01  0.026 0.014];
magnet_float.dim = [0.014 0.026 0.01 ];

magnet_fixed.magn = 1.23;
magnet_float.magn = 1.23;

magnet_fixed.magdir = [0  90]; % z
magnet_float.magdir = [0  0];  % x

offset = [0 -0.008 0.015];
N = 1001;
displ = linspace(-0.05, 0.05, N);
f_xyz = repmat(NaN, [3 N]);

for ii = 1:N
  f_xyz(:,ii) = magnetforces(magnet_fixed,magnet_float,offset+[displ(ii) 0 0]);
end

%% Plot

willfig('janssen'); clf; hold on
plot(displ,f_xyz(1,:),'Tag','x')
plot(displ,f_xyz(2,:),'Tag','y')
plot(displ,f_xyz(3,:),'Tag','z')
set(gca,'box','on');
set(gca,'ticklength',[0.02 0.05])
axis tight

% We want the vertical axis 5% less tight:
ylim = get(gca,'ylim');
ylim_range = ylim(2)-ylim(1);
yp = 0.05;
ylim = [ylim(1)-yp*ylim_range ylim(2)+yp*ylim_range];
set(gca,'ylim',ylim);

% add some more xticks:
set(gca,'xtick',[-0.05 -0.025 0 0.025 0.05]);

xlabel('Horizontal $x$ displacement, m')
ylabel('Forces, N')
text(-0.03,-10,'$F_x$')
text(0.01,   10,'$F_y$')
text(-0.025, 10,'$F_z$')

if ~simple_graph
  [h1 h2] = draworigin;
  set([h1 h2],'linestyle','--');
  colourplot
  % labelplot
  matlabfrag('fig/janssen-repro','dpi',3200);
end
