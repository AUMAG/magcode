function magB = cuboid_field(mag,xyz)

J = mag.magn*mag.magdir;

if size(xyz,1) == 3
elseif size(xyz,2) == 3
  warning('xyz should be column vectors of position stacked along rows')
  xyz = transpose(xyz);
else
  error('xyz funny size')
end

% Set up variables
xprime = mag.vertices(:,1);
yprime = mag.vertices(:,2);
zprime = mag.vertices(:,3);
x = xyz(1,:);
y = xyz(2,:);
z = xyz(3,:);

% Mesh everything for vectorisation
[x,xprime] = meshgrid(x,xprime);
[y,yprime] = meshgrid(y,yprime);
[z,zprime] = meshgrid(z,zprime);

Bxcontribution = cuboid_field_x(x,y,z,xprime,yprime,zprime,J);
Bycontribution = cuboid_field_y(x,y,z,xprime,yprime,zprime,J);
Bzcontribution = cuboid_field_z(x,y,z,xprime,yprime,zprime,J);

magB = Bxcontribution+Bycontribution+Bzcontribution;

end

function B = cuboid_field_x(x,y,z,xprime,yprime,zprime,J)

Jx = J(1);

% Solve field
D = repmat(Jx/(4*pi)*(-1).^(0:7)',1,length(x));
zeta = sqrt((x-xprime).^2+(y-yprime).^2+(z-zprime).^2);
Bx = D.*atan((y-yprime).*(z-zprime)./((x-xprime).*zeta));
By = D.*log(-z+zprime+zeta);
Bz = D.*log(-y+yprime+zeta);

% Solve singularities:
if any(any((abs(y-yprime)<eps | abs(z-zprime)<eps) & abs(x-xprime)<eps))
    index = (abs(y-yprime)<eps | abs(z-zprime)<eps) & abs(x-xprime)<eps;
    Bx(index) = 0;
end
if any(any(abs(zeta-z+zprime)<eps))
    index = abs(zeta-z+zprime)<eps;
    By(index) = D(index).*log(1./zeta(index));
end
if any(any(abs(zeta-y+yprime)<eps))
    index = abs(zeta-y+yprime)<eps;
    Bz(index) = D(index).*log(1./zeta(index));
end

B = [sum(Bx);sum(By);sum(Bz)];

end

function B = cuboid_field_y(x,y,z,xprime,yprime,zprime,J)

Jy = J(2);

% Solve field
D = repmat(Jy/(4*pi)*(-1).^(0:7)',1,length(x));
zeta = sqrt((x-xprime).^2+(y-yprime).^2+(z-zprime).^2);
Bx = D.*log(-z+zprime+zeta);
By = D.*atan((x-xprime).*(z-zprime)./((y-yprime).*zeta));
Bz = D.*log(-x+xprime+zeta);

% Solve singularities:
if any(any(abs(zeta-z+zprime)<eps))
    index = abs(zeta-z+zprime)<eps;
    Bx(index) = D(index).*log(1./zeta(index));
end
if any(any(abs(zeta-x+xprime)<eps))
    index = abs(zeta-x+xprime)<eps;
    Bz(index) = D(index).*log(1./zeta(index));
end
if any(any((abs(x-xprime)<eps | abs(z-zprime)<eps) & abs(y-yprime)<eps))
    index = (abs(x-xprime)<eps | abs(z-zprime)<eps) & abs(y-yprime)<eps;
    By(index) = 0;
end

B = [sum(Bx);sum(By);sum(Bz)];

end

function B = cuboid_field_z(x,y,z,xprime,yprime,zprime,J)

Jz = J(3);

% Solve field
D = repmat(Jz/(4*pi)*(-1).^(0:7)',1,length(x));
zeta = sqrt((x-xprime).^2+(y-yprime).^2+(z-zprime).^2);
Bx = D.*log(-y+yprime+zeta);
By = D.*log(-x+xprime+zeta);
Bz = D.*atan((x-xprime).*(y-yprime)./((z-zprime).*zeta));

% Solve singularities:
if any(any(abs(zeta-y+yprime)<eps))
    index = abs(zeta-y+yprime)<eps;
    Bx(index) = D(index).*log(1./zeta(index));
end
if any(any(abs(zeta-x+xprime)<eps))
    index = abs(zeta-x+xprime)<eps;
    By(index) = D(index).*log(1./zeta(index));
end
if any(any((abs(x-xprime)<eps | abs(y-yprime)<eps) & abs(z-zprime)<eps))
    index = (abs(x-xprime)<eps | abs(y-yprime)<eps) & abs(z-zprime)<eps;
    Bz(index) = 0;
end

B = [sum(Bx);sum(By);sum(Bz)];

end