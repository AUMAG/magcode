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

Txyz = zeros([3, size(offset,2)]);

for ii=[0,1]
  for jj=[0,1]
    for kk=[0,1]
      for ll=[0,1]
        for mm=[0,1]
          for nn=[0,1]
            
            Cu = (-1)^ii.*size1(1) - offset(1,:) - lever(1,:);
            Cv = (-1)^kk.*size1(2) - offset(2,:) - lever(2,:);
            Cw = (-1)^mm.*size1(3) - offset(3,:) - lever(3,:);
            
            u = offset(1,:) - (-1)^ii.*size1(1) + (-1)^jj.*size2(1);
            v = offset(2,:) - (-1)^kk.*size1(2) + (-1)^ll.*size2(2);
            w = offset(3,:) - (-1)^mm.*size1(3) + (-1)^nn.*size2(3);
            
            Cuu = 2*Cu + u;
            Cvv = 2*Cv + v;
            Cww = 2*Cw + w;
            
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
              Ex(a) = 1/8*w(a).*(-w2(a)-2*Cw(a).*w(a)-8*Cv(a).*abs(w(a))+w(a).*Cww(a).*log(w2(a)));
              Ey(a) = 1/8*w(a).*(+w2(a)+2*Cw(a).*w(a)+8*Cu(a).*abs(w(a))-w(a).*Cww(a).*log(w2(a)));
              Ez(a) = 1/4*(Cu(a)-Cv(a))*w2(a).*log(w2(a));
            end
            
            if any(b)
              Ex(b) = -1/4*Cw(b).*v(b).*(v(b)+2*abs(v(b)));
              Ey(b) = -1/4*Cw(b).*v2(b).*(log(v2(b))-1);
              Ez(b) =  1/72*v(b).*(2*v2(b)+36*Cu(b).*abs(v(b))+9*v(b).*Cvv(b).*log(v2(b)));
            end
                        
            if any(c)
              Ex(c) =  1/4*Cw(c).*u2(c).*(log(u2(c))-1);
              Ey(c) =  1/4*Cw(c).*(u2(c)+2*abs(u(c)).*u(c));
              Ez(c) = -1/72*u(c).*(2*u2(c)+36*Cv(c).*abs(u(c))+9*u(c).*Cuu(c).*log(u2(c)));
            end
            
            if any(d)
              
              Ex(d) = 1/8.*(...
                - Cww(d).*s2(d) + ...
                - 2.*s(d).*(v(d).*Cww(d)+2.*Cvv(d).*w(d)) + ...
                ...
                + 2.*Cww(d).*(s2(d)-v2(d)).*log(s(d)+v(d)) + ...
                - 8.*u(d).*v(d).*(Cw(d)+w(d)).*log(s(d)-u(d)) + ...
                ...
                + 8.*Cv(d).*u(d).*w(d).*acoth(s(d)./u(d)) + ...
                + 4.*w(d).*(v(d).*Cvv(d)-w(d).*Cww(d)).*acoth(s(d)./v(d)) + ...
                + 4.*u(d).*(v(d).*Cvv(d)-w(d).*Cww(d)).*atan(u(d).*v(d)./(w(d).*s(d))) + ...
                           0);
              
              Ey(d) = 1/8*(...
                + Cww(d).*s2(d) + ...
                + 2.*s(d).*(u(d).*Cww(d)+2.*Cuu(d).*w(d)) ...
                ...
                - 2.*Cww(d).*(s2(d)-u2(d)).*log(s(d)+u(d)) + ...
                + 8.*u(d).*v(d).*(Cw(d)+w(d)).*(log(s(d)-v(d))-1) + ...
                ...
                - 8.*Cu(d).*v(d).*w(d).*acoth(s(d)./v(d)) + ...
                - 4.*w(d).*(u(d).*Cuu(d)-w(d).*Cww(d)).*acoth(s(d)./u(d)) + ...
                - 4.*v(d).*(u(d).*Cuu(d)-w(d).*Cww(d)).*atan(u(d).*v(d)./(w(d).*s(d))) + ...
                ...
...                + 8.*v(d).*w(d).*(Cw(d)+w(d)).*atan(u(d)./w(d)) + ...
                         0);
              
              Ez(d) = 1/36.*(...
                - u(d).^3 + ...
                + v(d).^3 + ...
                +  6.*w2(d).*(v(d)-u(d)) + ...
                + 18.*s(d).*(Cu(d).*v(d)-u(d).*Cv(d)) + ...
                ...
                + 18.*u(d).*v(d).*( Cuu(d).*(log(s(d)-u(d))-1) - Cvv(d).*(log(s(d)-v(d))-1) ) + ...
                +  9.*Cvv(d).*(v2(d)-w2(d)).*log(s(d)+u(d)) + ...
                -  9.*Cuu(d).*(u2(d)-w2(d)).*log(s(d)+v(d)) ...
                ...
                +  6.*w(d).*(w2(d)-3.*v(d).*Cvv(d)).*atan(u(d)./w(d)) + ...
                -  6.*w(d).*(w2(d)-3.*u(d).*Cuu(d)).*atan(v(d)./w(d)) + ...
                - 18.*w(d).*(Cvv(d).*v(d)-u(d).*Cuu(d)).*atan(u(d).*v(d)./(w(d).*s(d))) + ...
                            0);
              
            end
            
            Txyz = Txyz + (-1)^(ii+jj+kk+ll+mm+nn)*[Ex; Ey; Ez];
            
          end
        end
      end
    end
  end
end

torque_zz = Txyz.*br1*br2/(16*pi^2*1e-7);

end
% \end{mfunction}
