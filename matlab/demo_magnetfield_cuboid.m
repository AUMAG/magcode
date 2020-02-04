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
  'magn',1,'magdir',[0 0 0]);

mag_cuboid.magdir = [0 0 5];
magnetdraw(mag_cuboid,[0; 0; 0]);


%% top 

N = 50;
[ptx,pty] = meshgrid(linspace(-0.2,0.2,N),linspace(-0.2,0.2,N));
ptz = repmat(0.2,size(ptx));

[Bx,By,Bz,ptx,pty,ptz] = magnetfield(mag_cuboid,ptx,pty,ptz);
Bmag = sqrt(Bx.^2+By.^2+Bz.^2);

h = surf(ptx,pty,ptz,Bmag);
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;

%% bottom 

N = 50;
[ptx,pty] = meshgrid(linspace(-0.2,0.2,N),linspace(-0.2,0.2,N));
ptz = repmat(-0.2,size(ptx));

[Bx,By,Bz,ptx,pty,ptz] = magnetfield(mag_cuboid,ptx,pty,ptz);
Bmag = sqrt(Bx.^2+By.^2+Bz.^2);

h = surf(ptx,pty,ptz,Bmag);
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;

