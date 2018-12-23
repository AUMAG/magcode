%% Demo of cuboid torque equations

%% ZZ
%
% Geometries and axis ordering to match Janssen (2010).
% DOI: 10.1109/TMAG.2010.2043224
%
% Note the magnitude and shape of the graphs does not appear to
% published results exactly.

size1  = [ 10;  26;  14 ]/2/1000;
size2  = [ 10;  26;  14 ]/2/1000;
lever  = [  0;   0; -47 ]/1000;
offset = [  0;  -8;  15 ]/1000;
J  = 1.23;
J1 = [0;0;+J];
J2 = [0;0;-J];

N = 100;
displ_range  = linspace(0,50/1000,N);
displ = [1;0;0]*displ_range;

torque_zz = cuboid_torque_z_z(size1,size2,offset+displ,lever,J1,J2);

figure(1);
plot(displ_range*1000,torque_zz)
title('ZZ')
legend('TX','TY','TZ')
xlabel('Displacement, mm')
ylabel('Torque, N.m')


%% ZY
%
% Geometries and axis ordering to match Janssen (2011).
% DOI: 10.1109/TMAG.2011.2154315
%
% Note the magnitude and shape of the graphs does not appear to
% published results exactly.

size1  = [ 13;  5;   7 ]/1000;
size2  = [ 13;  7;   5 ]/1000;
lever  = [  0;  0; -47 ]/1000;
offset = [ -8;  0;  15 ]/1000;
J  = 1.23;
JZ = [0; 0;J];
JY = [0;-J;0];

N = 100;
displ_range  = linspace(0,50/1000,N);
displ = [0;1;0]*displ_range;

torque_zy = cuboid_torque_z_y(size1,size2,offset+displ,lever,JZ,JY);

figure(1);
plot(displ_range*1000,torque_zy)
title('ZY')
legend('TY','TX','TZ')
xlabel('Displacement, mm')
ylabel('Torque, N.m')

  