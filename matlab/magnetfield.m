function magB = magnetfield(mag,xyz,varargin)
%MAGNETFIELD Calculate magnetic field from a magnet source

switch mag.type
  
  case 'cuboid'
  
    magB = calc_cuboid_field(mag,xyz);
    
  case 'cylinder'
    
    
end

end


function magB = calc_cuboid_field(mag,xyz)

Br = 1;
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
%define equation constant D  - ravaud 2009
D=((-1).^(ii+jj+kk));
%define equation constant zeta  - ravaud 2009
zeta = sqrt((x-x_m(ii)+eps).^2+(y-y_m(jj)+eps).^2+(z-z_m(kk)+eps).^2);
zeta(isnan(zeta))=0;

% +Z magnetisation

Bx = (Br*u0/(4*pi))*sum(D.*-real(log((y-y_m(jj)+eps)+zeta)),3);
By = (Br*u0/(4*pi))*sum(D.*-real(log((x-x_m(ii)+eps)+zeta)),3);
Bz = (Br*u0/(4*pi))*sum(D.*(atan(((x-x_m(ii)+eps).*(y-y_m(jj)+eps))./((z-z_m(kk)+eps).*zeta))),3);

magB = sum([Bx;By;Bz],2);


end