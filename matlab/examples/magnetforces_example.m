%% Setup
%
% In case you don't have the various bits'n'pieces that I use to create
% my Matlab graphics (probably likely).

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  simple_graph = true;
else
  simple_graph = false;
end

%% Parallel magnets example
%
% Reproducing the results of Akoun & Yonnet (1984)

magnet_fixed.dim = [0.02 0.012 0.006];
magnet_float.dim = [0.012 0.02 0.006];

magnet_fixed.magn = 0.38;
magnet_float.magn = 0.38;

magnet_fixed.magdir = [0 0 1]; % z
magnet_float.magdir = [0 0 1]; % z

N = 501;
offset = repmat([-0.004; -0.004; 0.008],[1 N]);
displ = linspace(0, 0.03, N);
displ_range = offset+[1; 0; 0]*displ;

f1_xyz = magnetforces(magnet_fixed,magnet_float,displ_range);

willfig('akoun'); clf; hold on
plot(displ,f1_xyz(1,:),'Tag','x')
plot(displ,f1_xyz(2,:),'Tag','y')
plot(displ,f1_xyz(3,:),'Tag','z')
xlabel('Horizontal $x$ displacement, m')
ylabel('Forces, N')
set(gca,'box','on','ticklength',[0.02 0.05])
text(0.004,-0.5,'$F_x$','interpreter','LaTeX')
text(0.004, 0.8,'$F_y$','interpreter','LaTeX')
text(0.004,-1.7,'$F_z$','interpreter','LaTeX')

if ~simple_graph
  draworigin([0 0],'h','--')
  colourplot
  matlabfrag('fig/akoun-repro','dpi',3200);
end

%% Orthogonal
%
% Replicating the results of Janssen et al. (2009)

magnet_fixed.dim = [0.01  0.026 0.014];
magnet_float.dim = [0.014 0.026 0.01 ];

magnet_fixed.magn = 1.23;
magnet_float.magn = 1.23;

magnet_fixed.magdir = [0 0 1]; % z
magnet_float.magdir = [1 0 0];  % x

N = 501;
offset = repmat([0; -0.008; 0.015],[1 N]);
displ = linspace(-0.05, 0.05, N);
displ_range = offset+[1; 0; 0]*displ;

f_xyz = magnetforces(magnet_fixed,magnet_float,displ_range);

%% Plot

willfig('janssen'); clf; hold on
plot(displ,f_xyz(1,:),'Tag','x')
plot(displ,f_xyz(2,:),'Tag','y')
plot(displ,f_xyz(3,:),'Tag','z')
set(gca,'box','on','ticklength',[0.02 0.05])

% add some more xticks:
set(gca,'xtick',[-0.05 -0.025 0 0.025 0.05]);

xlabel('Horizontal $x$ displacement, m')
ylabel('Forces, N')
text(-0.03,-10,'$F_x$','interpreter','LaTeX')
text( 0.01, 10,'$F_y$','interpreter','LaTeX')
text(-0.025,10,'$F_z$','interpreter','LaTeX')

if ~simple_graph
  axistight
  draworigin([0 0],'--')
  colourplot
  matlabfrag('fig/janssen-repro','dpi',3200);
end
