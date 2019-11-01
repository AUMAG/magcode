function magB = magnetfield(mag,xyz,varargin)
%MAGNETFIELD Calculate magnetic field from a magnet source

switch mag.type
  
  case 'cuboid'
  
    magB = calc_cuboid_field(mag,xyz);
    
  case 'cylinder'
    
    
    
end

end


function magB = calc_cuboid_field(mag,xyz)

Br = mag.magn;
u0 = 4*pi*10^-7 ;   


if size(xyz,1) == 3
elseif size(xyz,2) == 3
  warning('xyz should be column vectors of position stacked along rows')
  xyz = transpose(xyz);
else
  error('xyz funny size')
end

X = xyz(1,:);
Y = xyz(2,:);
Z = xyz(3,:);
N = size(X);

x_m = [-mag.dim(1)/2 +mag.dim(1)/2];
y_m = [-mag.dim(2)/2 +mag.dim(2)/2];
z_m = [-mag.dim(3)/2 +mag.dim(3)/2];

[ii,jj,kk]=meshgrid(1:2,1:2,1:2);
%reshape index array for 8 summations and size of field array
ii=repmat(reshape(ii,[1,1,8]),N);
jj=repmat(reshape(jj,[1,1,8]),N);
kk=repmat(reshape(kk,[1,1,8]),N);
%reshape field array to match index array
x=ones(N(1),N(2),8).*X;
y=ones(N(1),N(2),8).*Y;
z=ones(N(1),N(2),8).*Z;
X_m = reshape(x_m(ii),N(1),N(2),8);
Y_m = reshape(y_m(jj),N(1),N(2),8);
Z_m = reshape(z_m(kk),N(1),N(2),8);
%define equation constant D  - ravaud 2009
D=((-1).^(ii+jj+kk));
%define equation constant zeta  - ravaud 2009
zeta = sqrt((x-X_m+eps).^2+(y-Y_m+eps).^2+(z-Z_m+eps).^2);
zeta(isnan(zeta))=0;

% Set up zero arrays
Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));

% X magnetisation
Jx = Br*mag.magdir(1);
Bx = Bx + (Jx/(4*pi))*sum(D.*(atan(((y-Y_m+eps).*(z-Z_m+eps))./((x-X_m+eps).*zeta))),3);
By = By + (Jx/(4*pi))*sum(D.*-real(log((z-Z_m+eps)+zeta)),3);
Bz = Bz + (Jx/(4*pi))*sum(D.*-real(log((y-Y_m+eps)+zeta)),3);

% Y magnetisation
Jy = Br*mag.magdir(2);
Bx = Bx + (Jy/(4*pi))*sum(D.*-real(log((z-Z_m+eps)+zeta)),3);
By = By + (Jy/(4*pi))*sum(D.*(atan(((x-X_m+eps).*(z-Z_m+eps))./((y-Y_m+eps).*zeta))),3);
Bz = Bz + (Jy/(4*pi))*sum(D.*-real(log((x-X_m+eps)+zeta)),3);

% Z magnetisation
Jz = Br*mag.magdir(3);
Bx = Bx + (Jz/(4*pi))*sum(D.*-real(log((y-Y_m+eps)+zeta)),3);
By = By + (Jz/(4*pi))*sum(D.*-real(log((x-X_m+eps)+zeta)),3);
Bz = Bz + (Jz/(4*pi))*sum(D.*(atan(((x-X_m+eps).*(y-Y_m+eps))./((z-Z_m+eps).*zeta))),3);

magB = sum([Bx;By;Bz],3);


end