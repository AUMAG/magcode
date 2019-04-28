%% cuboid_force_z_x
%
% Calculate the forces between two parallel cuboid magnets, magnetised in the
% z- and x-directions respectively.

% \START

% \begin{mfunction}{cuboid_force_z_x}
%
%  This is a translation of cuboid_force_z_y into a rotated coordinate system.
%
% \begin{center}
% \begin{tabular}{lll}
% Inputs:
% & |size1|=|(a,b,c)| & the half dimensions of the fixed magnet \\
% & |size2|=|(A,B,C)| & the half dimensions of the floating magnet \\
% & |offset|=|(dx,dy,dz)| & distance between magnet centres \\
% & |(J1,J2)| & magnetisation vectors of the magnets \\
% Outputs:
% & |forces_xyz|=|(Fx,Fy,Fz)| & Forces of the second magnet \\
% \end{tabular}
% \end{center}

function forces_xyz = cuboid_force_z_x(size1,size2,offset,J1,J2)

J1 = J1(3);
J2 = J2(1);

if ( abs(J1)<eps || abs(J2)<eps )
  forces_xyz  =  [0; 0; 0];
  return;
end

component_x = 0;
component_y = 0;
component_z = 0;

for ii = [1 -1]
  for jj = [1 -1]
    for kk = [1 -1]
      for ll = [1 -1]
        for pp = [1 -1]
          for qq = [1 -1]

            ind_sum = ii*jj*kk*ll*pp*qq;

            u = -offset(2,:) + size2(2)*jj - size1(2)*ii;
            v =  offset(1,:) + size2(1)*ll - size1(1)*kk;
            w =  offset(3,:) + size2(3)*qq - size1(3)*pp;
            r = sqrt(u.^2+v.^2+w.^2);

            tmp = zeros(size(u));
            atan_term_u = tmp; atan_term_v = tmp; atan_term_w = tmp;
            log_ru      = tmp; log_rw      = tmp; log_rv      = tmp;
            
            ind = abs(u) > eps;
            if any(ind), atan_term_u(ind) = atan(v(ind).*w(ind)./(r(ind).*u(ind))); end
            ind = abs(v) > eps;
            if any(ind), atan_term_v(ind) = atan(u(ind).*w(ind)./(r(ind).*v(ind))); end
            ind = abs(w) > eps;
            if any(ind), atan_term_w(ind) = atan(u(ind).*v(ind)./(r(ind).*w(ind))); end

            ind = abs(r-u) > eps;
            if any(ind), log_ru(ind) = log(r(ind)-u(ind)); end
            ind = abs(r+w) > eps;
            if any(ind), log_rw(ind) = log(r(ind)+w(ind)); end
            ind = abs(r+v) > eps;
            if any(ind), log_rv(ind) = log(r(ind)+v(ind)); end

            cx = ...
              + v.*w.*log_ru ...
              - v.*u.*log_rw ...
              - u.*w.*log_rv ...
              + 0.5*u.^2.*atan_term_u ...
              + 0.5*v.^2.*atan_term_v ...
              + 0.5*w.^2.*atan_term_w;

            cy = ...
              - 0.5*(u.^2-v.^2).*log_rw ...
              + u.*w.*log_ru ...
              + u.*v.*atan_term_v ...
              + 0.5*w.*r;

            cz = ...
              - 0.5*(u.^2-w.^2).*log_rv ...
              + u.*v.*log_ru ...
              + u.*w.*atan_term_w ...
              + 0.5*v.* r;

            component_x = component_x + ind_sum.*cx;
            component_y = component_y + ind_sum.*cy;
            component_z = component_z + ind_sum.*cz;

          end
        end
      end
    end
  end
end

forces_xyz = J1*J2/(4*pi*(4*pi*1e-7))*[ component_y; -component_x; component_z ];

end

% \end{mfunction}

