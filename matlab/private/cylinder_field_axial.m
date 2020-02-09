function magB = cylinder_field_axial(mag,xyz_g)

% Set up variables
M = mag.magn;
R = mag.dim(1);
L = mag.dim(2)/2;

xyz = transpose(mag.rotation)*xyz_g-mag.position;

xyz = xyz';
rho = sqrt(xyz(:,1).^2+xyz(:,2).^2);
Z = xyz(:,3);

% Solve field equations (Caciagli 2018)
zeta = [Z+L,Z-L];
alpha = 1./sqrt(zeta.^2+(rho+R).^2);
beta = zeta.*alpha;
gamma = (rho-R)./(rho+R);
ksq = (zeta.^2+(rho-R).^2)./(zeta.^2+(rho+R).^2);

[K,E,P] = ellipkepi(1-[gamma,gamma].^2,1-ksq);

P1 = K - 2./(1-ksq).*(K-E);
P2 = -gamma./(1-gamma.^2).*(P-K)-1./(1-gamma.^2).*(gamma.^2.*P-K);

% Evaluate the numeric singularities at rho = 0
P1(rho==0,:) = 0;
P2(rho==0,:) = pi/2;

Brho = M*R/pi*(alpha(:,1).*P1(:,1)-alpha(:,2).*P1(:,2));
Bz = M*R./(pi*(rho+R)).*(beta(:,1).*P2(:,1)-beta(:,2).*P2(:,2));

% Evaluate the z-field at rho = R (Ravaud 2010)
index = abs(rho-R)<eps;
if any(index)
    Bz(index) = imag(sum((-1).^[0,1].*alpha(index,:).*zeta(index,:).*(ellipticF(-asin(1./beta(index,:).^2),beta(index,:).^2)+ellipticK(beta(index,:).^2)),2))*M/2/pi;
end
    
% Convert to Cartesian coordinates
d = sqrt(xyz(:,1).^2+xyz(:,2).^2)+eps;
Bx = Brho.*xyz(:,1)./d;
By = Brho.*xyz(:,2)./d;

magB = [Bx';By';Bz'];

end

