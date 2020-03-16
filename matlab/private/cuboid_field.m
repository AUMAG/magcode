function [g_magB] = cuboid_field(mag,g_points)

J = mag.magn*mag.magdir;

points = transpose(mag.rotation)*(g_points - mag.position);

% Set up variables
vertx = mag.vertices(1,:);
verty = mag.vertices(2,:);
vertz = mag.vertices(3,:);

% Mesh for vectorisation
[ptx,vertx] = meshgrid(points(1,:),vertx);
[pty,verty] = meshgrid(points(2,:),verty);
[ptz,vertz] = meshgrid(points(3,:),vertz);

% Calculate components
if abs(J(1)) > eps
  [Bxx,Bxy,Bxz] = cuboid_field_x(ptx,pty,ptz,vertx,verty,vertz,J(1));
else
  Bxx = zeros(size(points(1,:))); Bxy = Bxx; Bxz = Bxx;
end

if abs(J(2)) > eps
  [Byx,Byy,Byz] = cuboid_field_y(ptx,pty,ptz,vertx,verty,vertz,J(2));
else
  Byx = zeros(size(points(1,:))); Byy = Byx; Byz = Byx;
end

if abs(J(3)) > eps
  [Bzx,Bzy,Bzz] = cuboid_field_z(ptx,pty,ptz,vertx,verty,vertz,J(3));
else
  Bzx = zeros(size(points(1,:))); Bzy = Bzx; Bzz = Bzx;
end

% Finish
l_magB = [Bxx+Byx+Bzx;
          Bxy+Byy+Bzy;
          Bxz+Byz+Bzz];

g_magB = mag.rotation*l_magB;
      
end


function [Bx,By,Bz] = cuboid_field_x(x,y,z,xprime,yprime,zprime,J)

DD = J/(4*pi)*(-1).^(0:7)';
D = repmat(DD,1,size(x,2));

% Solve field
zeta = sqrt((x-xprime).^2+(y-yprime).^2+(z-zprime).^2);
Bxc = D.*atan((y-yprime).*(z-zprime)./((x-xprime).*zeta));
Byc = D.*log(-z+zprime+zeta);
Bzc = D.*log(-y+yprime+zeta);

% Solve singularities:
if any(any((abs(y-yprime)<eps | abs(z-zprime)<eps) & abs(x-xprime)<eps))
    index = (abs(y-yprime)<eps | abs(z-zprime)<eps) & abs(x-xprime)<eps;
    Bxc(index) = 0;
end
if any(any(abs(zeta-z+zprime)<eps))
    index = abs(zeta-z+zprime)<eps;
    Byc(index) = D(index).*log(1./zeta(index));
end
if any(any(abs(zeta-y+yprime)<eps))
    index = abs(zeta-y+yprime)<eps;
    Bzc(index) = D(index).*log(1./zeta(index));
end

Bx = sum(Bxc,1);
By = sum(Byc,1);
Bz = sum(Bzc,1);

end

function [Bx,By,Bz] = cuboid_field_y(x,y,z,xprime,yprime,zprime,J)

DD = J/(4*pi)*(-1).^(0:7)';
D = repmat(DD,1,size(x,2));

% Solve field
zeta = sqrt((x-xprime).^2+(y-yprime).^2+(z-zprime).^2);
Bxc = D.*log(-z+zprime+zeta);
Byc = D.*atan((x-xprime).*(z-zprime)./((y-yprime).*zeta));
Bzc = D.*log(-x+xprime+zeta);

% Solve singularities:
if any(any(abs(zeta-z+zprime)<eps))
    index = abs(zeta-z+zprime)<eps;
    Bxc(index) = D(index).*log(1./zeta(index));
end
if any(any(abs(zeta-x+xprime)<eps))
    index = abs(zeta-x+xprime)<eps;
    Bzc(index) = D(index).*log(1./zeta(index));
end
if any(any((abs(x-xprime)<eps | abs(z-zprime)<eps) & abs(y-yprime)<eps))
    index = (abs(x-xprime)<eps | abs(z-zprime)<eps) & abs(y-yprime)<eps;
    Byc(index) = 0;
end

Bx = sum(Bxc,1);
By = sum(Byc,1);
Bz = sum(Bzc,1);

end

function [Bx,By,Bz] = cuboid_field_z(x,y,z,xprime,yprime,zprime,J)

DD = J/(4*pi)*(-1).^(0:7)';
D = repmat(DD,1,size(x,2));

% Solve field
zeta = sqrt((x-xprime).^2+(y-yprime).^2+(z-zprime).^2);
Bxc = D.*log(-y+yprime+zeta);
Byc = D.*log(-x+xprime+zeta);
Bzc = D.*atan((x-xprime).*(y-yprime)./((z-zprime).*zeta));

% Solve singularities:
if any(any(abs(zeta-y+yprime)<eps))
    index = abs(zeta-y+yprime)<eps;
    Bxc(index) = D(index).*log(1./zeta(index));
end
if any(any(abs(zeta-x+xprime)<eps))
    index = abs(zeta-x+xprime)<eps;
    Byc(index) = D(index).*log(1./zeta(index));
end
if any(any((abs(x-xprime)<eps | abs(y-yprime)<eps) & abs(z-zprime)<eps))
    index = (abs(x-xprime)<eps | abs(y-yprime)<eps) & abs(z-zprime)<eps;
    Bzc(index) = 0;
end

Bx = sum(Bxc,1);
By = sum(Byc,1);
Bz = sum(Bzc,1);

end
