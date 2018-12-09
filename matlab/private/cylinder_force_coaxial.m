% \begin{mfunction}{cylinder_force_coaxial}
  function force_axial = cylinder_force_coaxial(size1,size2,h_gap,J1,J2)
    
    % inputs
    
    r1 = size1(1);
    r2 = size2(1);
    
    % implicit
    
    z(4,:) = h_gap + size2(2)/2;
    z(3,:) = h_gap - size2(2)/2;
    z(2,:) =  size1(2)/2;
    z(1,:) = -size1(2)/2;
    
    C_d = zeros(size(h_gap));
    
    for ii = [1 2]
      for jj = [3 4]
        
        a1 = z(ii,:) - z(jj,:);
        
        if abs(a1)<eps
          continue % singularity at a1=0 (i.e., f_z = 0 for coincident faces)
        end
        
        a2 = 1 + ( (r1-r2)./a1 ).^2;
        a3 = sqrt( (r1+r2).^2 + a1.^2 );
        a4 = 4*r1.*r2./( (r1+r2).^2 + a1.^2 );
                
        if ( a2 == 1 || isnan(a2) ) % singularity at a2=1 (i.e., equal radii)
          [K, E] = ellipke(a4);
          PI_term = 0;
        else
          [K, E, PI] = ellipkepi( a4./(1-a2) , a4 );
          PI_term = (1-a1.^2./a3.^2).*PI;
        end
               
        f_z = a1.*a2.*a3.*( K - E./a2 - PI_term );
        
        C_d = C_d + (-1)^(ii+jj).*f_z;
        
      end
    end
    
    force_axial = J1*J2/(8*pi*1e-7)*C_d;
    
  end
% \end{mfunction}
