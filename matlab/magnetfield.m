function magB = magnetfield(mag,xyz,varargin)
%MAGNETFIELD Calculate magnetic field from a magnet source

switch mag.type
  
  case 'cuboid'
  
    magB = calc_cuboid_field(mag,xyz);
    
  case 'cylinder'
    
    if isequal(mag.magdir,[0;0;1]) && isequal(mag.dir,[0;0;1])
        magB = calc_cylinder_axial_field(mag,xyz);
    end
    
end

end


function magB = calc_cuboid_field(mag,xyz)

% Br = mag.magn;

J = mag.magn*mag.magdir;

% u0 = 4*pi*10^-7 ;   


if size(xyz,1) == 3
elseif size(xyz,2) == 3
  warning('xyz should be column vectors of position stacked along rows')
  xyz = transpose(xyz);
else
  error('xyz funny size')
end

% Set up variables
n = size(xyz,2);
Bx = zeros(1,n);
By = zeros(1,n);
Bz = zeros(1,n);

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


% N = size(X);
% 
% x_m = [-mag.dim(1)/2 +mag.dim(1)/2];
% y_m = [-mag.dim(2)/2 +mag.dim(2)/2];
% z_m = [-mag.dim(3)/2 +mag.dim(3)/2];
% 
% [ii,jj,kk]=meshgrid(1:2,1:2,1:2);
% %reshape index array for 8 summations and size of field array
% ii=repmat(reshape(ii,[1,1,8]),N);
% jj=repmat(reshape(jj,[1,1,8]),N);
% kk=repmat(reshape(kk,[1,1,8]),N);
% %reshape field array to match index array
% x=ones(N(1),N(2),8).*X;
% y=ones(N(1),N(2),8).*Y;
% z=ones(N(1),N(2),8).*Z;
% X_m = reshape(x_m(ii),N(1),N(2),8);
% Y_m = reshape(y_m(jj),N(1),N(2),8);
% Z_m = reshape(z_m(kk),N(1),N(2),8);
% %define equation constant D  - ravaud 2009
% D=((-1).^(ii+jj+kk));
% %define equation constant zeta  - ravaud 2009
% zeta = sqrt((x-X_m).^2+(y-Y_m).^2+(z-Z_m).^2);
% 
% % Set up zero arrays
% Bx = zeros(size(X));
% By = zeros(size(Y));
% Bz = zeros(size(Z));
% 
% % X magnetisation
% Jx = Br*mag.magdir(1);
% Bx = Bx + (Jx/(4*pi))*sum(D.*(atan(((y-Y_m).*(z-Z_m))./((x-X_m).*zeta))),3);
% By = By + (Jx/(4*pi))*sum(D.*-real(log((z-Z_m)+zeta)),3);
% Bz = Bz + (Jx/(4*pi))*sum(D.*-real(log((y-Y_m)+zeta)),3);
% 
% % Y magnetisation
% Jy = Br*mag.magdir(2);
% Bx = Bx + (Jy/(4*pi))*sum(D.*-real(log((z-Z_m)+zeta)),3);
% By = By + (Jy/(4*pi))*sum(D.*(atan(((x-X_m).*(z-Z_m))./((y-Y_m).*zeta))),3);
% Bz = Bz + (Jy/(4*pi))*sum(D.*-real(log((x-X_m)+zeta)),3);
% 
% % Z magnetisation
% Jz = Br*mag.magdir(3);
% Bx = Bx + (Jz/(4*pi))*sum(D.*-real(log((y-Y_m)+zeta)),3);
% By = By + (Jz/(4*pi))*sum(D.*-real(log((x-X_m)+zeta)),3);
% Bz = Bz + (Jz/(4*pi))*sum(D.*(atan(((x-X_m).*(y-Y_m))./((z-Z_m).*zeta))),3);
% 
% magB = sum([Bx;By;Bz],3);

end

function B = cuboid_field_x(x,y,z,xprime,yprime,zprime,J)

Jx = J(1);

D = Jx/(4*pi)*(-1).^(0:7)';
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
    By(index) = D(any(index')).*log(1./zeta(index));
end
if any(any(abs(zeta-y+yprime)<eps))
    index = abs(zeta-y+yprime)<eps;
    Bz(index) = D(any(index')).*log(1./zeta(index));
end

B = [sum(Bx);sum(By);sum(Bz)];

end

function B = cuboid_field_y(x,y,z,xprime,yprime,zprime,J)

Jy = J(2);

D = Jy/(4*pi)*(-1).^(0:7)';
zeta = sqrt((x-xprime).^2+(y-yprime).^2+(z-zprime).^2);
Bx = D.*log(-z+zprime+zeta);
By = D.*atan((x-xprime).*(z-zprime)./((y-yprime).*zeta));
Bz = D.*log(-x+xprime+zeta);

% Solve singularities:
if any(any(abs(zeta-z+zprime)<eps))
    index = abs(zeta-z+zprime)<eps;
    Bx(index) = D(any(index')).*log(1./zeta(index));
end
if any(any(abs(zeta-x+xprime)<eps))
    index = abs(zeta-x+xprime)<eps;
    Bz(index) = D(any(index')).*log(1./zeta(index));
end
if any(any((abs(x-xprime)<eps | abs(z-zprime)<eps) & abs(y-yprime)<eps))
    index = (abs(x-xprime)<eps | abs(z-zprime)<eps) & abs(y-yprime)<eps;
    By(index) = 0;
end

B = [sum(Bx);sum(By);sum(Bz)];

end

function B = cuboid_field_z(x,y,z,xprime,yprime,zprime,J)

Jz = J(3);

D = Jz/(4*pi)*(-1).^(0:7)';
zeta = sqrt((x-xprime).^2+(y-yprime).^2+(z-zprime).^2);
Bx = D.*log(-y+yprime+zeta);
By = D.*log(-x+xprime+zeta);
Bz = D.*atan((x-xprime).*(y-yprime)./((z-zprime).*zeta));

% Solve singularities:
if any(any(abs(zeta-y+yprime)<eps))
    index = abs(zeta-y+yprime)<eps;
    Bx(index) = D(any(index')).*log(1./zeta(index));
end
if any(any(abs(zeta-x+xprime)<eps))
    index = abs(zeta-x+xprime)<eps;
    By(index) = D(any(index')).*log(1./zeta(index));
end
if any(any((abs(x-xprime)<eps | abs(y-yprime)<eps) & abs(z-zprime)<eps))
    index = (abs(x-xprime)<eps | abs(y-yprime)<eps) & abs(z-zprime)<eps;
    Bz(index) = 0;
end

B = [sum(Bx);sum(By);sum(Bz)];

end


function magB = calc_cylinder_axial_field(mag,xyz)
%field equations from Caciagli (2018)

%set up variables
M = mag.magn;
R = mag.dim(1);
L = mag.dim(2)/2;

mu0 = 4*pi*10^(-7);
xyz = xyz';

rho = sqrt(xyz(:,1).^2+xyz(:,2).^2);
Z = xyz(:,3);

%compute equations
zeta = [Z+L,Z-L];
alpha = 1./sqrt(zeta.^2+(rho+R).^2);
beta = zeta.*alpha;
gamma = (rho-R)./(rho+R);
ksq = (zeta.^2+(rho-R).^2)./(zeta.^2+(rho+R).^2);

[K,E,P] = ellipkepi(1-[gamma,gamma].^2,1-ksq);

P1 = K - 2./(1-ksq).*(K-E);
P2 = -gamma./(1-gamma.^2).*(P-K)-1./(1-gamma.^2).*(gamma.^2.*P-K);

%evaluate the singularities at rho = 0
P1(rho==0,:) = 0;
P2(rho==0,:) = pi/2;

Brho = M*R/pi*(alpha(:,1).*P1(:,1)-alpha(:,2).*P1(:,2));
Bz = M*R./(pi*(rho+R)).*(beta(:,1).*P2(:,1)-beta(:,2).*P2(:,2));

%evaluate the z-field at rho = R (Ravaud 2010)
index = abs(rho-R)<eps;
if any(index)
    Bz(index) = imag(sum((-1).^[0,1].*alpha(index,:).*zeta(index,:).*(ellipticF(-asin(1./beta(index,:).^2),beta(index,:).^2)+ellipticK(beta(index,:).^2)),2))*M/2/pi;
end
    
%convert to Cartesian coordinates
d = sqrt(xyz(:,1).^2+xyz(:,2).^2)+eps;
Bx = Brho.*xyz(:,1)./d;
By = Brho.*xyz(:,2)./d;

magB = [Bx';By';Bz'];

end


