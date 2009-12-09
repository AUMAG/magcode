
%% Parallel

magnet_fixed.dim = [0.04 0.04 0.04];
magnet_float.dim =  magnet_fixed.dim;

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.magdir = [0  90]; % z
magnet_float.magdir = [0  90]; % z

offset = [0 0 0.12];
N = 51;
displ = linspace(-0.12, 0.12, N);

%%

f_xyz = repmat(NaN, [3 N]);
for ii = 1:N

  f_xyz(:,ii) = magnetforces(magnet_fixed,magnet_float,offset+[0 displ(ii) 0]);

end

%%
willfig('test'); clf; hold on
plot(displ,f_xyz(1,:),'UserData','x')
plot(displ,f_xyz(2,:),'UserData','y')
plot(displ,f_xyz(3,:),'UserData','z')
colourplot
labelplot



%% Orthogonal
%
% Replicating Janssen's plot

! ~/bin/mtangle magnetforces

magnet_fixed.dim = [0.01  0.026 0.014];
magnet_float.dim = [0.014 0.026 0.01 ];

magnet_fixed.magn = 1.23;
magnet_float.magn = 1.23;

magnet_fixed.magdir = [0  90]; % z
magnet_float.magdir = [0  0];  % x

offset = [0 -0.008 0.015];
N = 201;
displ = linspace(-0.05, 0.05, N);
f_xyz = repmat(NaN, [3 N]);

for ii = 1:N

  f_xyz(:,ii) = magnetforces(magnet_fixed,magnet_float,offset+[displ(ii) 0 0]);

end

willfig('janssen test2'); clf; hold on
plot(displ,f_xyz(1,:),'Tag','x')
plot(displ,f_xyz(2,:),'Tag','y')
plot(displ,f_xyz(3,:),'Tag','z')
colourplot
labelplot

%% Stiffnesses of the above

! ~/bin/mtangle magnetforces

magnet_fixed.dim = [0.01  0.026 0.014];
magnet_float.dim = [0.014 0.026 0.01 ];

magnet_fixed.magn = 1.23;
magnet_float.magn = 1.23;

magnet_fixed.magdir = [0  90]; % z
magnet_float.magdir = [0   0];  % x

offset = [0.02 -0.008 0.015];
N = 101;
displ = linspace(-0.05, 0.05, N);
dx = displ(2)-displ(1);
f_xyz = repmat(NaN, [3 N]);
k_xyz = repmat(NaN, [3 N]);

for ii = 1:N

  [f_xyz(:,ii) k_xyz(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[displ(ii) 0 0],'force','stiffness');

end

fkx = -gradient(f_xyz(1,:),dx);

willfig('test2'); clf; hold on
subplot(3,2,1);
plot(displ,f_xyz(1,:))
title('Force')
ylabel('x')
subplot(3,2,2);
plot(displ,k_xyz(1,:))
plot(displ,fkx,'.-');
title('Stiffness')


for ii = 1:N

  [f_xyz(:,ii) k_xyz(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[0 displ(ii) 0],'force','stiffness');

end

fky = -gradient(f_xyz(2,:),dx);

subplot(3,2,3);
plot(displ,f_xyz(2,:))
ylabel('y')
subplot(3,2,4);
plot(displ,k_xyz(2,:))
plot(displ,fky,'.-');


for ii = 1:N

  [f_xyz(:,ii) k_xyz(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[0 0 displ(ii)],'force','stiffness');

end

fkz = -gradient(f_xyz(3,:),dx);

subplot(3,2,5);
plot(displ,f_xyz(1,:))
ylabel('z')
subplot(3,2,6);
plot(displ,k_xyz(3,:))
plot(displ,fkz,'.-');