%Matthew Forbes - matthew.forbes@adelaide.edu.au

%Eqn from Ravaud 2009  - MAGNETIC FIELD PRODUCED BY A PARALLELEPI-
%                        PEDIC MAGNET OF VARIOUS AND UNIFORM POLAR-
%                        ZATION

function H = cuboid_field(sigma,phi,theta,x_m,y_m,z_m,X,Y,Z)

%theta      - CCW angle from the +x-axis to the polarisation vector 
%           (XY-plane // about Z)
%phi        - CW angle from the +z-axis to the polarisation vector
%x1,y1,z1   - reference corner of magnet
%x2,y2,z2   - opposite corner of magnet, defining direction of coordinate
%           system and dimension of magnet
%           x_m = [x1,x2] ...
%sigma      - magnetic field strength, tesla
%X,Y,Z      - field points

%unit conversion
M0 = sigma*8*10^5; % A/m

%magnetic permeability of free space
u0 = 4*pi*10^-7 ;  %H/m

%summation of function over i,j,k 1:2
[i,j,k]=meshgrid(1:2,1:2,1:2);
d=size(X);

%reshape index array for 8 summations and size of field array
i=repmat(reshape(i,[1,1,8]),d);
j=repmat(reshape(j,[1,1,8]),d);
k=repmat(reshape(k,[1,1,8]),d);

%reshape field array to match index array
x=ones(d(1),d(2),8).*X;
y=ones(d(1),d(2),8).*Y;
z=ones(d(1),d(2),8).*Z;

%define equation constant D 
D=((-1).^(i+j+k));

%define equation constant zeta 
zeta = sqrt((x-x_m(i)+eps).^2+(y-y_m(j)+eps).^2+(z-z_m(k)+eps).^2);
zeta(isnan(zeta))=0;

%calculate magnetic field strength
Hx = ((M0*sin(phi)*cos(theta))/(4*pi))*sum(D.*(atan(((y-y_m(j)+eps)...
                        .*(z-z_m(k)+eps))./((x-x_m(i)+eps).*zeta))),3)...
    +((M0*sin(phi)*sin(theta))/(4*pi))*sum(D.*-real(log((z-z_m(k)+eps)...
                        +zeta)),3)...
    +((M0*cos(phi))/(4*pi))*sum(D.*-real(log((y-y_m(j)+eps)+zeta)),3);

Hy = ((M0*sin(phi)*cos(theta))/(4*pi))*sum(D.*-real(log((z-z_m(k)+eps)...
                                                            +zeta)),3)...
    +((M0*sin(phi)*sin(theta))/(4*pi))*sum(D.*(atan(((x-x_m(i)+eps)...
                        .*(z-z_m(k)+eps))./((y-y_m(j)+eps).*zeta))),3)...
    +((M0*cos(phi))/(4*pi))*sum(D.*-real(log((x-x_m(i)+eps)+zeta)),3);

Hz = ((M0*sin(phi)*cos(theta))/(4*pi))*sum(D.*-real(log((y-y_m(j)+eps)...
                                                            +zeta)),3)...
    +((M0*sin(phi)*sin(theta))/(4*pi))*sum(D.*-real(log((x-x_m(i)+eps)...
                                                            +zeta)),3)...
    +((M0*cos(phi))/(4*pi))*sum(D.*(atan(((x-x_m(i)+eps)...
                           .*(y-y_m(j)+eps))./((z-z_m(k)+eps).*zeta))),3);
                       
H = cat(3,Hx,Hy,Hz);

%magnetic flux density
B = u0.*H;

end

