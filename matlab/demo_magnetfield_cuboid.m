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


%% Cuboid

mag_cuboid = magnetdefine(...
  'type','cuboid',...
  'dim',[0.1 0.2 0.3],...
  'magn',1,'magdir',[5 0 5]);

magnetdraw(mag_cuboid,[0; 0; 0]);



%% bottom 

N = 50;
[ptx,pty] = meshgrid(linspace(-0.2,0.2,N),linspace(-0.2,0.2,N));
ptz = repmat(-0.2,size(ptx));

[Bx,By,Bz] = magnetfield(mag_cuboid,ptx,pty,ptz);
Bmag = sqrt(Bx.^2+By.^2+Bz.^2);

h = surf(ptx,pty,ptz,Bmag);
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;

%% curve

phi = linspace(0,pi);
r = 0.2;
y = linspace(-0.1,0.1);

[pp,yy] = meshgrid(phi,y);

xx = r*cos(pp);
zz = r*sin(pp);

[Bx,By,Bz] = magnetfield(mag_cuboid,xx,yy,zz);
Bmag = sqrt(Bx.^2+By.^2+Bz.^2);

h = surf(xx,yy,zz,Bmag);
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;


