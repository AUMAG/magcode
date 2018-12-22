function [B_rho, B_z] = cylinder_field_axial(M,R,L,rho,z)
% CYLINDER_FIELD_AXIAL  Magnetic field of axially-magnetised cylinder
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

P21 = -gamma./gg.*(PI_1-K_1)-1./gg.*(gamma.^2.*PI_1-K_1);
P22 = -gamma./gg.*(PI_2-K_2)-1./gg.*(gamma.^2.*PI_2-K_2);

B_rho = 4e-7*M*R*(alpha_1.*P11-alpha_2.*P12);
B_z   = 4e-7*M*R./(rho+R).*(beta_1.*P21-beta_2.*P22);

end

