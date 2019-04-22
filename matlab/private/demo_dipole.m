%% Demo of dipole equations


%% Verify Yung 1998 example (Table 1)

ma = [100;0;0];
mb = [100;0;0];
rab = [1;0;0];
f1 = dipole_forcetorque(ma,mb,rab)

ma = [0;-100;0];
mb = [100;0;0];
rab = [1;0;0];
dipole_forcetorque(ma,mb,rab)
f2 = dipole_forcetorque(ma,mb,rab)

ma = [100/sqrt(2);0;100/sqrt(2)];
mb = [0;0;100];
rab = [1;0;0];
dipole_forcetorque(ma,mb,rab)
f3 = dipole_forcetorque(ma,mb,rab)
f3*sqrt(2)

ma = [0;0;100];
mb = [0;0;100];
rab = [1;0;0];
dipole_forcetorque(ma,mb,rab)
f4 = dipole_forcetorque(ma,mb,rab)

%% Offset dipoles in repulsion

ma  = [0;0;+1];
mb  = [0;0;+1];
z = linspace(0,0.5);
rab = [0;0;1]*z + [0.1;0;0];
[f5,t5] = dipole_forcetorque(ma,mb,rab);

figure(1); hold on; clf

yyaxis left
plot(z,f5(3,:))
xlabel('Displacement')
ylabel('Force')

yyaxis right
plot(z,t5(2,:))
ylabel('Torque')


%% Rotating dipoles

ma = [0;0;1];

R = 0.5;
t = linspace(0,pi);

mb  = [-cos(t);zeros(size(t));sin(t)];
rab = [R*cos(t);zeros(size(t));R*sin(t)];

[f6,t6] = dipole_forcetorque(ma,mb,rab);
figure(1); hold on; clf

yyaxis left
plot(t,f6(3,:))
ylabel('Force')

yyaxis right
plot(t,t6(2,:))
ylabel('Torque')

xlabel('Angle, deg.')
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlim(t([1 end]))

