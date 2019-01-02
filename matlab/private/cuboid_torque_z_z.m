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

Tx = zeros([1, size(offset,2)]);
Ty = Tx;
Tz = Tx;

for ii=[0,1]
  for jj=[0,1]
    for kk=[0,1]
      for ll=[0,1]
        for mm=[0,1]
          for nn=[0,1]
            
            Cu = (-1)^ii.*size1(1) - lever(1);
            Cv = (-1)^kk.*size1(2) - lever(2);
            Cw = (-1)^mm.*size1(3) - lever(3);
            
            u = offset(1,:) - (-1)^ii.*size1(1) + (-1)^jj.*size2(1);
            v = offset(2,:) - (-1)^kk.*size1(2) + (-1)^ll.*size2(2);
            w = offset(3,:) - (-1)^mm.*size1(3) + (-1)^nn.*size2(3);
            
            u2 = u.^2;
            v2 = v.^2;
            w2 = w.^2;
            s2 = u2+v2+w2;
            s = sqrt(s2);
            
            % find indexes where cuboid faces align
            a = (u2<eps) & (v2<eps);
            b = (u2<eps) & (w2<eps);
            c = (v2<eps) & (w2<eps);
            
            % and all those that do not
            d = ~a & ~b & ~c;
            
            Ex = nan(1,size(offset,2));
            Ey = nan(1,size(offset,2));
            Ez = nan(1,size(offset,2));
            
            if any(a)
%              error a
              Ex(a) = 1/8*w(a).*(-w2(a)-2*Cw.*w(a)-8*Cv.*abs(w(a))+w(a).*(2*Cw+w(a)).*log(w2(a)));
              Ey(a) = 1/8*w(a).*(+w2(a)+2*Cw.*w(a)+8*Cu.*abs(w(a))-w(a).*(2*Cw+w(a)).*log(w2(a)));
              Ez(a) = 1/4*(Cu-Cv)*w2(a).*log(w2(a));
            end
            
            if any(b)
              Ex(b) = -1/4*Cw*v(b).*(v(b)+2*abs(v(b)));
              Ey(b) = -1/4*Cw*v(b).^2.*(log(v2(b))-1);
              Ez(b) = 1/72*v(b).*(2*v(b).^2+36*Cu*abs(v(b))+9*v(b).*(2*Cv+v(b)).*log(v2(b)));
            end
                        
            if any(c)
              Ex(c) = 1/4*Cw*u2(c).*(log(u2(c))-1);
              Ey(c) = 1/4*Cw.*(u2(c)+2*abs(u(c)).*u(c));
              Ez(c) = -1/72*u(c).*(2*u2(c)+36*Cv*abs(u(c))+9*u(c).*(2*Cu+u(c)).*log(u2(c)));
            end
            
            if any(d)
              
              Ex(d) = 1/8.*(...
                - s2(d).*(2.*Cw+w(d)) + ...
                - 2.*s(d).*( 2.*Cw.*v(d) + 4.*w(d).*Cv + 3.*w(d).*v(d) ) + ...
                + 8.*Cv.*u(d).*w(d).*acoth(s(d)./u(d)) + ...
                - 8.*u(d).*v(d).*(Cw+w(d)).*log(s(d)-u(d)) + ...
                + 2.*(2.*Cw+w(d)).*(u2(d)+w2(d)).*log(v(d)+s(d)) + ...
                + 4.*w(d).*(v2(d)-w2(d)+2.*Cv.*v(d)-2.*w(d).*Cw).*acoth(s(d)./v(d)) + ...
                - 4.*u(d).*(w2(d)-v2(d)+2.*Cw.*w(d)-2.*v(d).*Cv).*atan(u(d).*v(d)./(w(d).*s(d))) ...
                           );
              
              Ey(d) = (1/8)*...
                ((2.*Cw+w(d)).*u(d).^2-8.*u(d).*v(d).*(Cw+w(d))+8.*u(d).*v(d).*(Cw+w(d)).*log(s(d)-v(d))...
                +4.*Cw.*u(d).*s(d)+6.*w(d).*s(d).*u(d)+(2.*Cw+w(d)).*(v(d).^2+w(d).^2)+...
                4.*w(d).*(w(d).^2+2.*Cw.*w(d)-u(d).*(2.*Cu+u(d))).*acoth(s(d)./u(d))+...
                4.*v(d).*(-2.*Cu.*w(d).*acoth(s(d)./v(d))+2.*w(d).*(Cw+w(d)).*atan(u(d)./w(d))...
                +(w(d).^2+2.*Cw.*w(d)-u(d).*(2.*Cu+u(d))).*atan(u(d).*v(d)./(w(d).*s(d))))...
                -2.*(2.*Cw+w(d)).*(v(d).^2+w(d).^2).*log(u(d)+s(d))+8.*Cu.*w(d).*s(d));
              
              Ez(d) = (1/36).*(-u(d).^3-18.*v(d).*u(d).^2-6.*u(d).*(w(d).^2+6.*Cu...
                .*v(d)-3.*v(d).*(2.*Cv+v(d))+3.*Cv.*s(d))+v(d).*(v(d).^2+6.*(w(d).^2+...
                3.*Cu.*s(d)))+6.*w(d).*(w(d).^2-3.*v(d).*(2.*Cv+v(d))).*atan(u(d)./w(d))...
                -6.*w(d).*(w(d).^2-3.*u(d).*(2.*Cu+u(d))).*atan(v(d)./w(d))-9.*...
                (2.*(v(d).^2+2.*Cv.*v(d)-u(d).*(2.*Cu+u(d))).*w(d).*atan(u(d).*v(d)./(w(d).*s(d)))...
                -2.*u(d).*(2.*Cu+u(d)).*v(d).*log(s(d)-u(d))-(2.*Cv+v(d)).*(v(d).^2-w(d).^2)...
                .*log(u(d)+s(d))+2.*u(d).*v(d).*(2.*Cv+v(d)).*log(s(d)-v(d))+(2.*Cu+...
                u(d)).*(u(d).^2-w(d).^2).*log(v(d)+s(d))));
              
            end
            
            Tx = Tx + (-1)^(ii+jj+kk+ll+mm+nn)*Ex;
            Ty = Ty + (-1)^(ii+jj+kk+ll+mm+nn)*Ey;
            Tz = Tz + (-1)^(ii+jj+kk+ll+mm+nn)*Ez;
            
          end
        end
      end
    end
  end
end

torque_zz = [Tx; Ty; Tz].*br1*br2/(16*pi^2*1e-7);

end
% \end{mfunction}
