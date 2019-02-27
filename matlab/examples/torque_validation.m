%% Torque validation
%
% Using Akoun & Yonnet's geometry with both parallel and orthogonal
% magnets, with reference torques calculated using the "O'Connell" method.

%% Parallel magnetisation

magnet_fixed = magnetdefine(    ...
  'type'  , 'cuboid',           ...
  'dim'   , [0.02 0.012 0.006], ...
  'magn'  , 0.38,               ...
  'magdir', [0 0 +1]            ...
);
magnet_float = magnetdefine(    ...
  'type'  , 'cuboid',           ...
  'dim'   , [0.012 0.02 0.006], ...
  'magn'  , 0.38,               ...
  'magdir', [0 0 +1]            ...
);

N = 101;
offset = repmat([-0.004; -0.004; 0.008],[1 N]);
displ = linspace(0, 0.03, N);
displ_range = offset+[1; 0; 0]*displ;

torq_para = magnetforces(magnet_fixed,magnet_float,displ_range,'torque');

figure(1); clf; hold on

h(1) = plot(displ*1000,torq_para(1,:));
h(2) = plot(displ*1000,torq_para(2,:));
h(3) = plot(displ*1000,torq_para(3,:));

load paratorques
plot(x*1000,T(:,1),'--','color',h(1).Color,'linewidth',2)
plot(x*1000,T(:,2),'--','color',h(2).Color,'linewidth',2)
plot(x*1000,T(:,3),'--','color',h(3).Color,'linewidth',2)

legend('x','y','z')
xlabel('Displacement, mm')
ylabel('Torques, Nm')

%% Orthogonal magnetisation

magnet_fixed = magnetdefine(    ...
  'type'  , 'cuboid',           ...
  'dim'   , [0.02 0.012 0.006], ...
  'magn'  , 0.38,               ...
  'magdir', [0 0 +1]            ...
);
magnet_float = magnetdefine(    ...
  'type'  , 'cuboid',           ...
  'dim'   , [0.012 0.02 0.006], ...
  'magn'  , 0.38,               ...
  'magdir', [0 +1 0]            ...
);

N = 51;
offset = repmat([-0.004; -0.004; 0.008],[1 N]);
displ = linspace(0, 0.03, N);
displ_range = offset+[1; 0; 0]*displ;

torq_orth = magnetforces(magnet_fixed,magnet_float,displ_range,'torque');

figure(1); clf; hold on

h(1) = plot(displ*1000,torq_orth(1,:),'.-');
h(2) = plot(displ*1000,torq_orth(2,:),'.-');
h(3) = plot(displ*1000,torq_orth(3,:),'.-');

load orthtorques
plot(x*1000,T(:,1),'--','color',h(1).Color)
plot(x*1000,T(:,2),'--','color',h(2).Color)
plot(x*1000,T(:,3),'--','color',h(3).Color)

legend('x','y','z')
xlabel('Displacement, mm')
ylabel('Torques, Nm')
