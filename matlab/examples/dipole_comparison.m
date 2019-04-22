%% Dipole comparison example


%%

clc
clear all
close all

JM = 1;
magJa = [0;0;1];
magJb = [0;0;1];

a = 0.02;
coil = magnetdefine('type','cuboid','magn',JM,'magdir',magJa,'dim',[a 2*a a]);
mag  = magnetdefine('type','cuboid','magn',JM,'magdir',magJb,'dim',[2*a a a]);

N = 50;
displ = [0; 0; 1]*linspace(2*a,5*a,N)+[2*a;a;0];

coil_moment = 1/(4*pi*1e-7)*magJa*coil.volume;
mag_moment  = 1/(4*pi*1e-7)*magJb*mag.volume;

[coilmagforce,coilmagtorque] = magnetforces(coil,mag,displ,'force','torque');
[dipoleforce,dipoletorque]   = magnetforces(coil,mag,displ,'force','torque','method','dipole');

figure(1); clf; hold on
subplot(2,1,1); cla; hold on
h = {0;0;0};
for ii = 1:3
  h{ii} = plot(displ(3,:),coilmagforce(ii,:));
end
for ii = 1:3
  plot(displ(3,:),dipoleforce(ii,:),'--','color',h{ii}.Color)
end
legend('X','Y','Z')
xlabel('Displacement, m')
ylabel('Force, N')
title('Cuboid exact (solid) vs dipole (dashed)')

subplot(2,1,2); cla; hold on
h = {0;0;0};
for ii = 1:3
  h{ii} = plot(displ(3,:),coilmagtorque(ii,:));
end
for ii = 1:3
  plot(displ(3,:),dipoletorque(ii,:),'--','color',h{ii}.Color)
end
legend('X','Y','Z')
xlabel('Displacement, m')
ylabel('Torque, Nm')
title('Cuboid exact (solid) vs dipole (dashed)')


