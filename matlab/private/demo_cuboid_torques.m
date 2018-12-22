%% Demo of cuboid torque equations

%% ZZ

size1  = [ 10;  26;  14 ]/2/1000;
size2  = [ 14;  26;  10 ]/2/1000;
lever  = [  0;   0; -47 ]/1000;
offset = [  0;  -8;  15 ]/1000;
J = 1.23;
JZ = [0;0;J];

N = 100;
displ_range  = linspace(0,50/1000,N);
displ = [1;0;0]*displ_range;

torque_zz = cuboid_torque_z_z(size1,size2,offset+displ,lever,JZ,-JZ);

figure(1);
plot(displ_range,torque_zz)
title('ZZ')
legend('TX','TY','TZ')


%% ZY

size1  = [ 10;  26;  14 ]/2/1000;
size2  = [ 14;  26;  10 ]/2/1000;
lever  = [  0;   0; -47 ]/1000;
offset = [  0;  -8;  15 ]/1000;
J = 1.23;
JZ = [0;0;J];
JY = [0;J;0];

N = 100;
displ_range  = linspace(0,50/1000,N);
displ = [1;0;0]*displ_range;

torque_zy = cuboid_torque_z_y(size1,size2,offset+displ,lever,JZ,JY);

figure(1);
plot(displ_range,torque_zy)
title('ZY')
legend('TX','TY','TZ')
