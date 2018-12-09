%% cylinder_force_coaxial(J1,J2,r1,r2,h1,h2,displ)
%
% Calculates axial force between two cylindrical coaxial magnets.
%
% INPUTS
%        J = magnetisation strength [T]
%        r = magnet radii [m]
%        h = magnet height [m]
%    displ = displacement(s) between magnet centres [m]
%
% OUTPUTS
%    force = axial force on the second magnet [N]
%
% Code is vectorised.

function force_axial = cylinder_force_coaxial(J1,J2,r1,r2,h1,h2,displ_input)

displ    = displ_input(:);
K        = zeros(numel(displ),4);
E        = zeros(numel(displ),4);
PI_term  = zeros(numel(displ),4);

z = [h1 -h1 h2 -h2]/2;
ii = [1, 2, 1, 2];
jj = [3, 3, 4, 4];

a1   = z(ii)-z(jj)-displ;
a1sq = a1.^2;
a2   = 1+(r2-r1).^2./a1sq;
a3sq = (r1+r2).^2 + a1sq ;
a3   = sqrt(a3sq);
a4   = 4*r1.*r2./a3sq;

% singularity at a1=0 (i.e., f_z = 0 for coincident faces)
ind_zero = abs(a1)<eps;
if any(ind_zero(:))
  a2(ind_zero) = eps; % need to ensure E/a2 = 0
  K(ind_zero) = 0;
  E(ind_zero) = 0;
  PI_term(ind_zero) = 0;
end

% singularity at a2=1 (i.e., equal radii)
ind_singu = ( a2 == 1 | isnan(a2) );
if any(ind_singu(:))
  [K(ind_singu), E(ind_singu)] = ellipke(a4(ind_singu));
  PI_term(ind_singu) = 0;
end

% all remaining cases
ind = ~ind_zero & ~ind_singu;
if any(ind(:))
  [K(ind), E(ind), PI] = ellipkepi( a4(ind)./(1-a2(ind)) , a4(ind) );
  PI_term(ind) = (1-a1sq(ind)./a3sq(ind)).*PI;
end

C_d = (-1).^(ii+jj).*a1.*a2.*a3.*(K-E./a2-PI_term);
f_z = sum(C_d,2);

force_axial = J1*J2/(8*pi*1e-7)*f_z;
force_axial = reshape(force_axial,size(displ));

end

