%% Demo of magnetic field from axially-magnetised cylinder

%% Field calculations

M = 1;

W = 10;
H = 10;

R = 1;
L = 2;

N = 1000;
rho_range = W*linspace(-1,1,N);
z_range   = H*linspace(-1,1,N);

[rr,zz] = meshgrid(rho_range,z_range);

[B_rho, B_z] = cylinder_field_axial(M,R,L,rr,zz);
Bmag = sqrt(B_rho.^2+B_z.^2);

figure(1); clf; hold on

h = surf(rho_range,z_range,Bmag);
h.EdgeColor = 'none';
view(2)
axis equal
pbaspect([1 1 1])
xlabel('RHO')
ylabel('Z')

plot([-R -R R R -R],[L -L -L L L],'k-','linewidth',3)


%% Magnetic field lines (example only!)

% non-connected
N = 50;
t = linspace(pi,2*pi,N);
xx = W*cos(t);
yy = H*sin(t);

%plot(xx,yy,'.','markersize',20)

XY = stream2(rr,zz,B_rho,B_z,xx,yy,[1,10000]);
h = streamline(XY);
set(h(:),'color','white','linewidth',1)

axis([-W W -H H]/1.5)


% connected
startx2 = linspace(-0.9*W,-1.1*R,5);
starty2 = zeros(size(startx2));
XY = stream2(rr,zz,B_rho,B_z,startx2,starty2,[1,5000]);
h = streamline(XY);
set(h(:),'color','white','linewidth',1)