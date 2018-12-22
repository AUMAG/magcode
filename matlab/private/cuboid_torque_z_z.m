%% cuboid_torque_z_z
%
% Calculate the torques between two parallel cuboid magnets, both magnetised in the
% z-direction.
%
% In the Janssen code, the torque is calculated on the first magnet and the
% lever arm ("torque reference point") is also w.r.t. the first magnet.

% \START

% \begin{mfunction}{cuboid_torque_z_z}
%
%  The expressions here follow directly from \textcite{janssen2010-ietm}.
%  The code below was largely written by Allan Liu; thanks!
%  We have checked it against Janssen's own Matlab code and the two give
%  identical output.
%
%  Note that despite this verification this code produces results which are
%  inconsistent with the graph in the \textcite{janssen2010-ietm} paper.
%  This appears to have been an oversight in the publication.
%
% \begin{center}
% \begin{tabular}{lll}
%  Inputs:
%  & |size1|=|(a1,b1,c1)| & the half dimensions of the fixed magnet \\
%  & |size2|=|(a2,b2,c2)| & the half dimensions of the floating magnet \\
%  & |displ|=|(a,b,c)| & distance between magnet centres \\
%  & |lever|=|(d,e,f)| & distance between floating magnet and its centre of rotation \\
%  &     |(J,J2)| & magnetisations of the magnet in the z-direction \\
%  Outputs:
%  & |forces_xyz|=|(Fx,Fy,Fz)| & Forces of the second magnet \\
% \end{tabular}
% \end{center}

function torque_zz = cuboid_torque_z_z(size1,size2,offset,lever,J1,J2)

br1 = J1(3);
br2 = J2(3);

if br1==0 || br2==0
  torque_zz = 0*offset;
  return
end

Tx = zeros([1 size(offset,2)]);
Ty = Tx;
Tz = Tx;

for ii=[0,1]
  for jj=[0,1]
    for kk=[0,1]
      for ll=[0,1]
        for mm=[0,1]
          for nn=[0,1]
            
            Cu = (-1)^ii.*size1(1) - lever(1,:);
            Cv = (-1)^kk.*size1(2) - lever(2,:);
            Cw = (-1)^mm.*size1(3) - lever(3,:);
            
            u = offset(1,:) - (-1)^ii.*size1(1) + (-1)^jj.*size2(1);
            v = offset(2,:) - (-1)^kk.*size1(2) + (-1)^ll.*size2(2);
            w = offset(3,:) - (-1)^mm.*size1(3) + (-1)^nn.*size2(3);
            
            s = sqrt(u.^2+v.^2+w.^2);
            
            Ex=(1/8).*(...
              Cw.*(8.*v.*u-2*s.^2-4.*v.*s)-...
              w.*(-8.*v.*u+s.^2+8.*Cv.*s+6.*v.*s)+...
              2.*(2.*Cw+w).*(u.^2+w.^2).*log(v+s)+...
              4.*(...
                2.*Cv.*u.*w.*acoth(u./s) + ...
                w.*(v.^2+2.*Cv.*v-w.*(2.*Cw+w)).*acoth(v./s) - ...
                u.*(...
                  2*w.*(Cw+w).*atan(v./w) + ...
                  2*v.*(Cw+w).*log(s-u) + ...
                  (w.^2+2.*Cw.*w-v.*(2.*Cv+v)).*atan( u.*v./(w.*s) ) ...
                )...
              )...
              );
            
            Ey=(1/8)*...
              ((2.*Cw+w).*u.^2-8.*u.*v.*(Cw+w)+8.*u.*v.*(Cw+w).*log(s-v)...
              +4.*Cw.*u.*s+6.*w.*s.*u+(2.*Cw+w).*(v.^2+w.^2)+...
              4.*w.*(w.^2+2.*Cw.*w-u.*(2.*Cu+u)).*acoth(u./s)+...
              4.*v.*(-2.*Cu.*w.*acoth(v./s)+2.*w.*(Cw+w).*atan(u./w)...
              +(w.^2+2.*Cw.*w-u.*(2.*Cu+u)).*atan(u.*v./(w.*s)))...
              -2.*(2.*Cw+w).*(v.^2+w.^2).*log(u+s)+8.*Cu.*w.*s);
            
            Ez=(1/36).*(-u.^3-18.*v.*u.^2-6.*u.*(w.^2+6.*Cu...
              .*v-3.*v.*(2.*Cv+v)+3.*Cv.*s)+v.*(v.^2+6.*(w.^2+...
              3.*Cu.*s))+6.*w.*(w.^2-3.*v.*(2.*Cv+v)).*atan(u./w)...
              -6.*w.*(w.^2-3.*u.*(2.*Cu+u)).*atan(v./w)-9.*...
              (2.*(v.^2+2.*Cv.*v-u.*(2.*Cu+u)).*w.*atan(u.*v./(w.*s))...
              -2.*u.*(2.*Cu+u).*v.*log(s-u)-(2.*Cv+v).*(v.^2-w.^2)...
              .*log(u+s)+2.*u.*v.*(2.*Cv+v).*log(s-v)+(2.*Cu+...
              u).*(u.^2-w.^2).*log(v+s)));
            
            Tx = Tx + (-1)^(ii+jj+kk+ll+mm+nn)*Ex;
            Ty = Ty + (-1)^(ii+jj+kk+ll+mm+nn)*Ey;
            Tz = Tz + (-1)^(ii+jj+kk+ll+mm+nn)*Ez;
            
          end
        end
      end
    end
  end
end

torque_zz = real([Tx; Ty; Tz].*br1*br2/(16*pi^2*1e-7));

end
% \end{mfunction}
