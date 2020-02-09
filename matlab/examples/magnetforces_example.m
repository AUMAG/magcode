%% Setup

%% Parallel magnets example
%
% Reproducing the results of Akoun & Yonnet (1984)

magnet_fixed = magnetdefine('type','cuboid','dim',[0.02 0.012 0.006],'magn',0.38,'magdir','z');
magnet_float = magnetdefine('type','cuboid','dim',[0.012 0.02 0.006],'position',[-0.004; -0.004; 0.008],'magn',0.38,'magdir','z');

N = 501;
displ = linspace(0, 0.03, N);
displ_range = [1; 0; 0]*displ;

f1_xyz = magnetforces(magnet_fixed,magnet_float,displ_range);

figure(1); clf; hold on
plot(displ,f1_xyz(1,:),'Tag','x')
plot(displ,f1_xyz(2,:),'Tag','y')
plot(displ,f1_xyz(3,:),'Tag','z')
xlabel('Horizontal $x$ displacement, m')
ylabel('Forces, N')
set(gca,'box','on','ticklength',[0.02 0.05])
text(0.004,-0.5,'$F_x$','interpreter','LaTeX')
text(0.004, 0.8,'$F_y$','interpreter','LaTeX')
text(0.004,-1.7,'$F_z$','interpreter','LaTeX')

yline(0,'--');


%% Orthogonal
%
% Replicating the results of Janssen et al. (2009)

magnet_fixed = magnetdefine('type','cuboid','dim',[0.01  0.026 0.014],'magn',1.23,'magdir','z');
magnet_float = magnetdefine('type','cuboid','dim',[0.014 0.026 0.01 ],'position',[0; -0.008; 0.015],'magn',1.23,'magdir','x');

N = 501;
displ = linspace(-0.05, 0.05, N);
displ_range = [1; 0; 0]*displ;

f_xyz = magnetforces(magnet_fixed,magnet_float,displ_range);

% Plot

figure(2); clf; hold on
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

axis tight
yline(0,'--');
