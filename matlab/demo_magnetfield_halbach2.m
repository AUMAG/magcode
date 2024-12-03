%% Demo for magnetdraw

function demo_magnetfield_cuboid()

clc

%% Plot setup

figure(1); clf; hold on
view(3);
axis equal
camlight

title('Cuboid magnets')
xlabel('x');
ylabel('y');
zlabel('z');

mag_x = 0.1;
array_h = 0.1;
zplane_scale = 1.5;
array_d = 0.2;
xyplane_ex = 2;

mag_cuboid1 = magnetdefine(...
  'type','cuboid',...
  'dim',[mag_x array_d array_h],...
  'position',[0;0;0],...
  'rotation',[0, 0, 0],...
  'magn',1,'magdir',[0 0 -1]);

mag_cuboid2 = magnetdefine(...
  'type','cuboid',...
  'dim',[mag_x array_d array_h],...
  'position',[mag_x;0;0],...
  'rotation',[0, 0, 0],...
  'magn',1,'magdir',[1 0 0]);

mag_cuboid3 = magnetdefine(...
  'type','cuboid',...
  'dim',[mag_x array_d array_h],...
  'position',[-mag_x;0;0],...
  'rotation',[0, 0, 0],...
  'magn',1,'magdir',[-1 0 0]);

mag_cuboid4 = magnetdefine(...
  'type','cuboid',...
  'dim',[mag_x array_d array_h],...
  'position',[-2*mag_x;0;0],...
  'rotation',[0, 0, 0],...
  'magn',1,'magdir',[0 0 1]);

mag_cuboid5 = magnetdefine(...
  'type','cuboid',...
  'dim',[mag_x array_d array_h],...
  'position',[2*mag_x;0;0],...
  'rotation',[0, 0, 0],...
  'magn',1,'magdir',[0 0 1]);

magnetdraw(mag_cuboid1);
magnetdraw(mag_cuboid2);
magnetdraw(mag_cuboid3);
magnetdraw(mag_cuboid4);
magnetdraw(mag_cuboid5);


%% bottom 

N = 50;
[ptx,pty] = meshgrid( ...
   linspace(-xyplane_ex*array_d,xyplane_ex*array_d,N), ...
   linspace(-xyplane_ex*array_d,xyplane_ex*array_d,N));
ptz = repmat(-zplane_scale*array_h/2,size(ptx));

magB1 = magnetfield(mag_cuboid1,[ptx(:),pty(:),ptz(:)].');
magB2 = magnetfield(mag_cuboid2,[ptx(:),pty(:),ptz(:)].');
magB3 = magnetfield(mag_cuboid3,[ptx(:),pty(:),ptz(:)].');
magB4 = magnetfield(mag_cuboid4,[ptx(:),pty(:),ptz(:)].');
magB5 = magnetfield(mag_cuboid5,[ptx(:),pty(:),ptz(:)].');
magB = magB1 + magB2 + magB3 + magB4 + magB5;
Bmag = sqrt(magB(1,:).^2+magB(2,:).^2+magB(3,:).^2);

h = surf(ptx,pty,ptz,reshape(Bmag,size(ptx)));
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;

cmax = 0.5;
caxis([0 cmax])

%% top 

N = 50;
[ptx,pty] = meshgrid( ...
   linspace(-xyplane_ex*array_d,xyplane_ex*array_d,N), ...
   linspace(-xyplane_ex*array_d,xyplane_ex*array_d,N));
ptz = repmat(zplane_scale*array_h/2,size(ptx));

magB1 = magnetfield(mag_cuboid1,[ptx(:),pty(:),ptz(:)].');
magB2 = magnetfield(mag_cuboid2,[ptx(:),pty(:),ptz(:)].');
magB3 = magnetfield(mag_cuboid3,[ptx(:),pty(:),ptz(:)].');
magB4 = magnetfield(mag_cuboid4,[ptx(:),pty(:),ptz(:)].');
magB5 = magnetfield(mag_cuboid5,[ptx(:),pty(:),ptz(:)].');
magB = magB1 + magB2 + magB3 + magB4 + magB5;
Bmag = sqrt(magB(1,:).^2+magB(2,:).^2+magB(3,:).^2);

h = surf(ptx,pty,ptz,reshape(Bmag,size(ptx)));
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;

caxis([0 cmax])


%% middle 

N = 50;
pmin = -0.4;
pmax = 0.4;
[ptz,ptx] = meshgrid(linspace(pmin,pmax,N),linspace(pmin,pmax,N));
pty = repmat(0,size(ptz));

magB1 = magnetfield(mag_cuboid1,[ptx(:),pty(:),ptz(:)].');
magB2 = magnetfield(mag_cuboid2,[ptx(:),pty(:),ptz(:)].');
magB3 = magnetfield(mag_cuboid3,[ptx(:),pty(:),ptz(:)].');
magB4 = magnetfield(mag_cuboid4,[ptx(:),pty(:),ptz(:)].');
magB5 = magnetfield(mag_cuboid5,[ptx(:),pty(:),ptz(:)].');
magB = magB1 + magB2 + magB3 + magB4 + magB5;
Bmag = sqrt(magB(1,:).^2+magB(2,:).^2+magB(3,:).^2);

h = surf(ptx,pty,ptz,reshape(Bmag,size(ptx)));
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;

caxis([0 cmax])





