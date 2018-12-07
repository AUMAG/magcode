% \begin{mfunction}{cylinder_force_coaxial}
  function calc_out = cylinder_force_coaxial(size1,size2,h_gap,J1,J2)
    
    % inputs
    
    r1 = size1(1);
    r2 = size2(1);
    
    % implicit
    
    z = nan(4,length(h_gap));
    z(1,:) = -size1(2)/2;
    z(2,:) =  size1(2)/2;
    z(3,:) = h_gap - size2(2)/2;
    z(4,:) = h_gap + size2(2)/2;
    
    C_d = zeros(size(h_gap));
    
    for ii = [1 2]
      for jj = [3 4]
        
        a1 = z(ii,:) - z(jj,:);
        a2 = 1 + ( (r1-r2)./a1 ).^2;
        a3 = sqrt( (r1+r2).^2 + a1.^2 );
        a4 = 4*r1.*r2./( (r1+r2).^2 + a1.^2 );
                
        a2_ind = ( a2 == 1 | isnan(a2) );
        if a2_ind % singularity at a2=1 (i.e., equal radii)
          [K, E] = ellipkepi( a4./(1-a2) , a4 );
          PI_term = 0;
        else
          [K, E, PI] = ellipkepi( a4./(1-a2) , a4 );
          PI_term = (1-a1.^2./a3.^2).*PI;
        end
        
        if abs(a1)<eps
          f_z = 0; % singularity at a1=0 (i.e., coincident faces)
        else
          f_z = a1.*a2.*a3.*( K - E./a2 - PI_term );
        end        
        
        C_d = C_d + (-1)^(ii+jj).*f_z;
        
      end
    end
    
    calc_out = J1*J2/(8*pi*1e-7)*C_d;
    
  end
% \end{mfunction}

% \begin{mfunction}{ellipkepi}
% Complete elliptic integrals calculated with the arithmetric-geometric mean
% algorithms contained here: \url{http://dlmf.nist.gov/19.8}.
% Valid for $a<=1$ and $m<=1$.

  function [k,e,PI] = ellipkepi(a,m)
    
    a0 = 1;
    g0 = sqrt(1-m);
    s0 = m;
    nn = 0;
    
    p0 = sqrt(1-a);
    Q0 = 1;
    QQ = Q0;
    
    Q1 = 1;
    w1 = 1;
    
    while max(w1(:)) > eps % || max(Q1(:)) > eps %% <- this is probably correct but I need to test more thoroughly
      
      % for Elliptic I
      a1 = (a0+g0)/2;
      g1 = sqrt(a0.*g0);
      
      % for Elliptic II
      nn = nn + 1;
      c1 = (a0-g0)/2;
      w1 = 2^nn*c1.^2;
      s0 = s0 + w1;
      
      % for Elliptic III
      rr = p0.^2+a0.*g0;
      p1 = rr./(2.*p0);
      Q1 = 0.5*Q0.*(p0.^2-a0.*g0)./rr;
      QQ = QQ+Q1;
      
      a0 = a1;
      g0 = g1;
      Q0 = Q1;
      p0 = p1;
      
    end
    
    k = pi./(2*a1);
    e = k.*(1-s0/2);
    PI = pi./(4.*a1).*(2+a./(1-a).*QQ);
    
    im = find(m == 1);
    if ~isempty(im)
      k(im) = inf;
      e(im) = ones(length(im),1);
      PI(im) = inf;
    end
    
  end
% \end{mfunction}