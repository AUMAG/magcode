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

mag_cuboid = magnetdefine(...
  'type','cuboid',...
  'dim',[0.1 0.2 0.3],...
  'position',[0;-0.2;0],...
  'rotation',[0, pi/8, 0],...
  'magn',1,'magdir',[1 0 1]);

magnetdraw(mag_cuboid);


%% bottom 

N = 50;
[ptx,pty] = meshgrid(linspace(-0.2,0.2,N),linspace(-0.2,0.2,N));
ptz = repmat(-0.2,size(ptx));

magB = magnetfield(mag_cuboid,[ptx(:),pty(:),ptz(:)].');
Bmag = sqrt(magB(1,:).^2+magB(2,:).^2+magB(3,:).^2);

h = surf(ptx,pty,ptz,reshape(Bmag,size(ptx)));
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;

%% curve

phi = linspace(0,pi);
r = 0.2;
y = linspace(-0.2,0.1);

[pp,yy] = meshgrid(phi,y);

xx = r*cos(pp);
zz = r*sin(pp);

magB = magnetfield(mag_cuboid,[xx(:),yy(:),zz(:)].');
Bmag = sqrt(magB(1,:).^2+magB(2,:).^2+magB(3,:).^2);

h = surf(xx,yy,zz,reshape(Bmag,size(pp)));
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;


