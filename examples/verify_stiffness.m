%% Stiffnesses verification

%% Parallel stiffnesses

magnet_fixed.dim = [0.01  0.026 0.014];
magnet_float.dim = [0.014 0.026 0.01 ];

magnet_fixed.magn = 1.23;
magnet_float.magn = 1.23;

magnet_fixed.magdir = [0  90]; % z
magnet_float.magdir = [0  90]; % z

offset = [0.02 -0.008 0.015];
N = 51;
displ = linspace(-0.05, 0.05, N);
dd = displ(2)-displ(1);

f_x = repmat(NaN, [3 N]);
k_x = repmat(NaN, [3 N]);
f_y = repmat(NaN, [3 N]);
k_y = repmat(NaN, [3 N]);
f_z = repmat(NaN, [3 N]);
k_z = repmat(NaN, [3 N]);

for ii = 1:N

  [f_x(:,ii) k_x(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[displ(ii) 0 0],'force','stiffness');
  [f_y(:,ii) k_y(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[0 displ(ii) 0],'force','stiffness');
  [f_z(:,ii) k_z(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[0 0 displ(ii)],'force','stiffness');

end

fkx = -gradient(f_x(1,:),dd);
fky = -gradient(f_y(2,:),dd);
fkz = -gradient(f_z(3,:),dd);

willfig('parallel-stiffness'); clf; hold on
plot(displ,k_x(1,:),displ,k_y(2,:),displ,k_z(3,:),'Tag','z')
plot(displ,fkx,'.',displ,fky,'.',displ,fkz,'.');
legend('x','y','z')

%% Orthogonal stiffnessnes

magnet_fixed.dim = [0.01  0.026 0.014];
magnet_float.dim = [0.014 0.026 0.01 ];

magnet_fixed.magn = 1.23;
magnet_float.magn = 1.23;

magnet_fixed.magdir = [0 0 1];
magnet_float.magdir = [1 0 0];

N = 51;
offset = [0.02 -0.008 0.015];
displ = linspace(-0.05, 0.05, N);
dd = displ(2)-displ(1);

f_x = repmat(NaN, [3 N]);
k_x = repmat(NaN, [3 N]);
f_y = repmat(NaN, [3 N]);
k_y = repmat(NaN, [3 N]);
f_z = repmat(NaN, [3 N]);
k_z = repmat(NaN, [3 N]);

for ii = 1:N

  [f_x(:,ii) k_x(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[displ(ii) 0 0],'force','stiffness');
  [f_y(:,ii) k_y(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[0 displ(ii) 0],'force','stiffness');
  [f_z(:,ii) k_z(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[0 0 displ(ii)],'force','stiffness');

end

fkx = -gradient(f_x(1,:),dd);
fky = -gradient(f_y(2,:),dd);
fkz = -gradient(f_z(3,:),dd);

willfig('orth-stiffness-x'); clf; hold on
plot(displ,k_x(1,:),displ,k_y(2,:),displ,k_z(3,:),'Tag','z')
plot(displ,fkx,'.',displ,fky,'.',displ,fkz,'.');
legend('x','y','z')


%% Orthogonal stiffnessnes again
%

magnet_fixed.dim = [0.01  0.026 0.014];
magnet_float.dim = [0.014 0.026 0.01 ];

magnet_fixed.magn = 1.23;
magnet_float.magn = 1.23;

magnet_fixed.magdir = [0 0 1];
magnet_float.magdir = [0 1 0];

N = 51;
offset = [0.02 -0.008 0.015];
displ = linspace(-0.05, 0.05, N);
dd = displ(2)-displ(1);

f_x = repmat(NaN, [3 N]);
k_x = repmat(NaN, [3 N]);
f_y = repmat(NaN, [3 N]);
k_y = repmat(NaN, [3 N]);
f_z = repmat(NaN, [3 N]);
k_z = repmat(NaN, [3 N]);

for ii = 1:N

  [f_x(:,ii) k_x(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[displ(ii) 0 0],'force','stiffness');
  [f_y(:,ii) k_y(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[0 displ(ii) 0],'force','stiffness');
  [f_z(:,ii) k_z(:,ii)] = magnetforces(magnet_fixed,magnet_float,offset+[0 0 displ(ii)],'force','stiffness');

end

fkx = -gradient(f_x(1,:),dd);
fky = -gradient(f_y(2,:),dd);
fkz = -gradient(f_z(3,:),dd);

willfig('orth-stiffness-y'); clf; hold on
plot(displ,k_x(1,:),displ,k_y(2,:),displ,k_z(3,:),'Tag','z')
plot(displ,fkx,'.',displ,fky,'.',displ,fkz,'.');
legend('x','y','z')


