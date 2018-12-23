function [B_rho, B_phi, B_z] = cylinder_field_transverse(M,R,L,rho,phi,z)
% CYLINDER_FIELD_TRANSVERSE  Magnetic field of transversally-magnetised cylinder
%
% These equations are from Caciagli (2018).

xi_1 = z+L;
xi_2 = z-L;
alpha_1 = (xi_1.^2+(rho+R).^2).^(-1/2);
alpha_2 = (xi_2.^2+(rho+R).^2).^(-1/2);
beta_1 = xi_1.*alpha_1;
beta_2 = xi_2.*alpha_2;
gamma = (rho-R)./(rho+R);
gg = 1-gamma.^2;
m_1 = (xi_1.^2+(rho-R).^2)./(xi_1.^2+(rho+R).^2);
m_2 = (xi_2.^2+(rho-R).^2)./(xi_2.^2+(rho+R).^2);

[K_1,E_1,PI_1] = ellipkepi(gg,(1-m_1));
[K_2,E_2,PI_2] = ellipkepi(gg,(1-m_2));

P11 = K_1-2./(1-m_1).*(K_1-E_1);
P12 = K_2-2./(1-m_2).*(K_2-E_2);

P31 = 1./(1-m_1).*(K_1-E_1) - gamma.^2./gg.*(PI_1-K_1);
P32 = 1./(1-m_2).*(K_2-E_2) - gamma.^2./gg.*(PI_2-K_2);

P41 = gamma./gg.*(PI_1.*(1+gamma.^2) - 2*K_1) - P11;
P42 = gamma./gg.*(PI_2.*(1+gamma.^2) - 2*K_2) - P12;

B_rho = 2e-7*M*R*cos(phi)./rho.*(beta_1.*P41-beta_2.*P42);
B_phi = 4e-7*M*R*sin(phi)./rho.*(beta_1.*P31-beta_2.*P32);
B_z   = 4e-7*M*R*cos(phi).*(alpha_1.*P11-alpha_2.*P12);

end

