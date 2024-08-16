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

Nmag = 5;
mag_cuboid = cell(Nmag,1);

for ii = 1:5

mag_cuboid{ii} = magnetdefine(...
  'type','cuboid',...
  'dim',[0.1 0.1 0.1],...
  'position',[ii*0.12;-0.2;0],...
  'rotation',[0, pi/2+pi/2*ii, 0],...
  'magn',1,'magdir',[0 0 1]);

magnetdraw(mag_cuboid{ii});

end

%% bottom 

N = 50;
[ptx,pty] = meshgrid(linspace(-0.2,1,N),linspace(-0.4,0.2,N));
ptz = repmat(-0.1,size(ptx));

magB = nan(Nmag,3,N*N);
Bmag = nan(Nmag,N*N);

for ii = 1:5

  magB(ii,:,:) = magnetfield(mag_cuboid{ii},[ptx(:),pty(:),ptz(:)].');
  Bmag(ii,:) = sqrt(magB(ii,1,:).^2+magB(ii,2,:).^2+magB(ii,3,:).^2);

end

Bmag = sum(Bmag,1);

h = surf(ptx,pty,ptz,reshape(Bmag,size(ptx)));
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;

%% curve

phi = linspace(0,pi);
r = 0.2;
y = linspace(-0.4,0.1);

[pp,yy] = meshgrid(phi,y);

xx = r*cos(pp);
zz = r*sin(pp);

magB = magnetfield(mag_cuboid,[xx(:),yy(:),zz(:)].');
Bmag = sqrt(magB(1,:).^2+magB(2,:).^2+magB(3,:).^2);

h = surf(xx,yy,zz,reshape(Bmag,size(pp)));
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;


