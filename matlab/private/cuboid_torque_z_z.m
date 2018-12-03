
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

a1 = size1(1);
b1 = size1(2);
c1 = size1(3);

a2 = size2(1);
b2 = size2(2);
c2 = size2(3);

a = offset(1,:);
b = offset(2,:);
c = offset(3,:);

d = a+lever(1,:);
e = b+lever(2,:);
f = c+lever(3,:);

Tx = zeros([1 size(offset,2)]);
Ty = Tx;
Tz = Tx;

for ii=[0,1]
  for jj=[0,1]
    for kk=[0,1]
      for ll=[0,1]
        for mm=[0,1]
          for nn=[0,1]
            
            Cu = (-1)^ii.*a1 - d;
            Cv = (-1)^kk.*b1 - e;
            Cw = (-1)^mm.*c1 - f;
            
            u = a-(-1)^ii.*a1 + (-1)^jj.*a2;
            v = b-(-1)^kk.*b1 + (-1)^ll.*b2;
            w = c-(-1)^mm.*c1 + (-1)^nn.*c2;
            
            s=sqrt(u.^2+v.^2+w.^2);
            
            Ex=(1/8).*(...
              -2.*Cw.*(-4.*v.*u+s.^2+2.*v.*s)-...
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
            
            Tx=Tx+(-1)^(ii+jj+kk+ll+mm+nn)*Ex;
            Ty=Ty+(-1)^(ii+jj+kk+ll+mm+nn)*Ey;
            Tz=Tz+(-1)^(ii+jj+kk+ll+mm+nn)*Ez;
            
          end
        end
      end
    end
  end
end

torque_zz = real([Tx; Ty; Tz].*br1*br2/(16*pi^2*1e-7));

end
% \end{mfunction}
