%% Demo of magnetic field from axially-magnetised cylinder

%% Field calculations (side)

M = 1;

W = 10;
H = 10;

R = 1;
L = 2;

N = 1000;
rho_range = W*linspace(-1,1,N);
z_range   = H*linspace(-1,1,N);
phi = 0;

[rr,zz] = meshgrid(rho_range,z_range);

[B_rho, B_phi, B_z] = cylinder_field_transverse(M,R,L,rr,phi,zz);
Bmag = sqrt(B_rho.^2+B_phi.^2+B_z.^2);

figure(1); clf; hold on

h = surf(rho_range,z_range,Bmag);
h.EdgeColor = 'none';
view(2); axis equal; pbaspect([1 1 1])
xlabel('RHO')
ylabel('Z')

plot([-R -R R R -R],[L -L -L L L],'k-','linewidth',3)

%% Evenly spaced streamlines to visualise flux lines
%
% See: <https://github.com/keithfma/evenly_spaced_streamlines>

figure(2); clf; hold on

h = surf(rho_range,z_range,Bmag);
h.EdgeColor = 'none';
view(2); axis equal; pbaspect([1 1 1])
xlabel('RHO')
ylabel('Z')

plot([-R -R R R -R],[L -L -L L L],'k-','linewidth',3)

even_stream_arrow(rr,zz,B_rho,B_z,0.5,10,'Color','w','lineWidth',1,'ArrowLength',7);



%% Field calculations (top)

M = 1;

W = 10;
H = 10;

R = 1;
L = 2;

N = 1000;
xx_range = W*linspace(-1,1,N);
yy_range = H*linspace(-1,1,N);

[xx,yy] = meshgrid(xx_range,yy_range);

rr = sqrt(xx.^2+yy.^2);
pp = atan2(yy,xx);
zz = 0;

[B_rho, B_phi, B_z] = cylinder_field_transverse(M,R,L,rr,pp,zz);

B_x = B_rho.*cos(pp)-B_phi.*sin(pp);
B_y = B_rho.*sin(pp)+B_phi.*cos(pp);

Bmag = sqrt(B_x.^2+B_y.^2+B_z.^2);

figure(3); clf; hold on

h = surf(xx_range,yy_range,Bmag);
h.EdgeColor = 'none';
view(2); axis equal; pbaspect([1 1 1])
xlabel('X')
ylabel('Y')

plot(R*cos(linspace(0,2*pi)),R*sin(linspace(0,2*pi)),'k-','Linewidth',3)


%% Evenly spaced streamlines to visual flux lines
%
% See: <https://github.com/keithfma/evenly_spaced_streamlines>

figure(4); clf; hold on

h = surf(xx_range,yy_range,Bmag);
h.EdgeColor = 'none';
view(2); axis equal; pbaspect([1 1 1])

plot(R*cos(linspace(0,2*pi)),R*sin(linspace(0,2*pi)),'k-','Linewidth',3)

even_stream_arrow(xx,yy,B_x,B_y,0.5,20,'Color','w','lineWidth',1,'ArrowLength',7);
