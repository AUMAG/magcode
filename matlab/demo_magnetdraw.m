%% Demo for magnetdraw

function demo_magnetdraw()

%% Plot setup

figure(1); clf; hold on
view(3);
axis equal
camlight

title('Cuboid and cylinder magnet')
xlabel('x');
ylabel('y');
zlabel('z');


%% Cuboid

mag_cuboid = magnetdefine(...
  'type','cuboid',...
  'dim',[0.1 0.2 0.3],...
  'magdir','z',...
  'magn',1);

magnetdraw(mag_cuboid,[0; 0; 0]);


%% Cylinder

mag_cyl = magnetdefine(...
  'type','cylinder',...
  'dim',[0.1 0.2],...
  'magdir','z',...
  'magn',1);

magnetdraw(mag_cyl,[0.5; 0; 0]);

