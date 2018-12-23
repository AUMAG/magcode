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


%% Evenly spaced streamlines to visual flux lines
%
% See: <https://github.com/keithfma/evenly_spaced_streamlines>

figure(4); clf; hold on

h = surf(rho_range,z_range,Bmag);
h.EdgeColor = 'none';
view(2)
axis equal
pbaspect([1 1 1])
xlabel('RHO')
ylabel('Z')

plot([-R -R R R -R],[L -L -L L L],'k-','linewidth',3)

even_stream_arrow(rr,zz,B_rho,B_z,0.5,10,'Color','w','lineWidth',1,'ArrowLength',7);
