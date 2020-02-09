%% An example of drawing the magnets in 3D
%
% Replicating the Akoun and Yonnet geometry


magnet_fixed = magnetdefine('type','cuboid',...
  'dim',[0.02 0.012 0.006],...
  'magn',0.38,'magdir','z',...
  'color',[0.8 0 0]);
magnet_float = magnetdefine('type','cuboid',...
  'dim',[0.012 0.02 0.006],...
  'position',[-0.004; -0.004; 0.008],...
  'magn',0.38,'magdir','z',...
  'color',[0 0.8 0]);

N = 501;
displ = linspace(0, 0.03, N);
displ_range = [1; 0; 0]*displ;

figure(1); clf
f1_xyz = magnetforces(magnet_fixed,magnet_float,displ_range,'draw',true,'drawpath',true,'markpath',{'d','o','s'},'markpathN',2);
axis equal
view(3)
grid on

%%

figure(2); clf; hold on; box on
yline(0,'--');
h1 = plot(displ,f1_xyz(1,:),'linewidth',2);
h2 = plot(displ,f1_xyz(2,:),'linewidth',2);
h3 = plot(displ,f1_xyz(3,:),'linewidth',2);
xlabel('Horizontal $x$ displacement, m','interpreter','latex')
ylabel('Force $F_i$, N','interpreter','latex')
set(gca,'ticklength',[0.02 0.05])
ii = 150;
text(displ(ii),f1_xyz(1,ii),'$F_x$','interpreter','LaTeX','backgroundcolor','white','color',h1.Color,'edgecolor',h1.Color,'horizontalalignment','center')
text(displ(ii),f1_xyz(2,ii),'$F_y$','interpreter','LaTeX','backgroundcolor','white','color',h2.Color,'edgecolor',h2.Color,'horizontalalignment','center')
text(displ(ii),f1_xyz(3,ii),'$F_z$','interpreter','LaTeX','backgroundcolor','white','color',h3.Color,'edgecolor',h3.Color,'horizontalalignment','center')


