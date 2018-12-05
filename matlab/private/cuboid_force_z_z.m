%\begin{mfunction}{forces_calc_z_z}
%
%The expressions here follow directly from \textcite{akoun1984}.
%
%\begin{center}
%\begin{tabular}{lll}
% Inputs:
% & |size1|=|(a,b,c)| & the half dimensions of the fixed magnet \\
% & |size2|=|(A,B,C)| & the half dimensions of the floating magnet \\
% & |offset|=|(dx,dy,dz)| & distance between magnet centres \\
% &     |(J,J2)| & magnetisations of the magnet in the z-direction \\
% Outputs:
% & |forces_xyz|=|(Fx,Fy,Fz)| & Forces of the second magnet \\
%\end{tabular}
%\end{center}

function forces_xyz = cuboid_force_z_z(size1,size2,offset,J1,J2)

magconst = 1/(4*pi*(4*pi*1e-7));

J1 = J1(3);
J2 = J2(3);

if ( abs(J1)<eps || abs(J2)<eps )
  forces_xyz = [0; 0; 0];
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
            
            u = offset(1) + size2(1)*jj - size1(1)*ii;
            v = offset(2) + size2(2)*ll - size1(2)*kk;
            w = offset(3) + size2(3)*qq - size1(3)*pp;
            r = sqrt(u.^2+v.^2+w.^2);
            
            if w == 0
              atan_term = 0;
            else
              atan_term = atan(u.*v./(r.*w));
            end
            if abs(r-u) < eps
              log_ru = 0;
            else
              log_ru = log(r-u);
            end
            if abs(r-v) < eps
              log_rv = 0;
            else
              log_rv = log(r-v);
            end
            
            cx = ...
              + 0.5*(v.^2-w.^2).*log_ru ...
              + u.*v.*log_rv ...
              + v.*w.*atan_term...
              + 0.5*r.*u;
            
            cy = ...
              + 0.5*(u.^2-w.^2).*log_rv ...
              + u.*v.*log_ru ...
              + u.*w.*atan_term ...
              + 0.5*r.*v;
            
            cz = ...
              - u.*w.*log_ru ...
              - v.*w.*log_rv ...
              + u.*v.*atan_term ...
              - r.*w;
            
            ind_sum = ii*jj*kk*ll*pp*qq;
            component_x = component_x + ind_sum.*cx;
            component_y = component_y + ind_sum.*cy;
            component_z = component_z + ind_sum.*cz;
            
          end
        end
      end
    end
  end
end

forces_xyz = J1*J2*magconst.*[component_x; component_y; component_z];

end

%\end{mfunction}