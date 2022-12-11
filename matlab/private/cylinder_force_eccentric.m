% Based on CONWAY 2013: https://doi.org/10.1109/TMAG.2013.2251652
%
% \begin{mfunction}{cylinder_force_eccentric}
  function calc_out = cylinder_force_eccentric(size1,size2,h_gap,e_displ,J1,J2)
    
    r1 = size1(1);
    r2 = size2(1);
    
    z1 = -size1(2)/2;
    z2 =  size1(2)/2;
    z3 = h_gap - size2(2)/2;
    z4 = h_gap + size2(2)/2;
        
    h = [z4-z2; z3-z2; z4-z1; z3-z1];
    
    fn = @(t) [xdir(t,r1,r2,h,e_displ), zdir(t,r1,r2,h,e_displ)];
    fn_int = integral(fn,0,pi,'ArrayValued',true,'AbsTol',1e-6);
    
    calc_out = -1e7*J1*J2*r1*r2*fn_int/4/pi/pi;
    
    function gx = xdir(t,r,R,h,p)
      
      X = sqrt(r^2+R^2-2*r*R*cos(t));
      hh = h.^2;
      ff = (p+X)^2+hh;
      gg = (p-X)^2+hh;
      f = sqrt(ff);
      g = sqrt(gg);
      m = 1-gg./ff;  % equivalent to $m = 4pX/f^2$
      
      [KK, EE] = ellipke(m);
      [F2, E2] = arrayfun(@elliptic12,asin(h./g),1-m);
      
      Ta = f.*EE;
      Tb = (p^2-X^2).*KK./f;
      Tc = sign(p-X)*h.*( F2.*(EE-KK) + KK.*E2 - 1 );
      Td = -pi/2*h;
      
      T = cos(t)/p*(Ta+Tb+Tc+Td);
      gx = -T(1)+T(2)+T(3)-T(4);
      
    end
    
    function gz = zdir(t,r,R,h,p)
      
      XX = p^2+R^2-2*p*R*cos(t);
      rr = r.^2;
      X = sqrt(XX);
      hh = h.^2;
      ff = (r+X)^2+hh;
      gg = (r-X)^2+hh;
      f = sqrt(ff);
      g = sqrt(gg);
      m = 1-gg./ff;
      
      [KK, EE] = ellipke(m);
      [F2, E2] = arrayfun(@elliptic12,asin(h./g),1-m);
      
      Ta = +h.*f.*(EE-KK);
      Tb = -h.*KK.*(r-X)^2./f;
      Tc = abs(rr-XX).*( F2.*(EE-KK) + KK.*E2 - 1 );
      Td = 4/pi.*min(rr,XX);   % note $r^2+X^2 - \lvert r^2-X^2\rvert = 2\min(r^2,X^2)$
      
      T = (R-p.*cos(t))./(2.*r.*XX).*(Ta+Tb+Tc+Td);
      gz = -T(1)+T(2)+T(3)-T(4);
      
    end
    
  end
% \end{mfunction}


% \begin{mfunction}{elliptic12}
  function [F,E] = elliptic12(u,m)
    % ELLIPTIC12 evaluates the value of the Incomplete Elliptic Integrals
    % of the First, Second Kind.
    % GNU GENERAL PUBLIC LICENSE Version 2, June 1991
    % Copyright (C) 2007 by Moiseev Igor.
    
    % EDITED BY WSPR to optimise for numel(u)=numel(m)=1
    % TODO: re-investigate vectorising once the wrapper code is properly in place
    
    tol = eps; % making this 1e-6 say makes it slower??
        
    F = zeros(size(u)); E = F; Z = E;
    
    m(m<eps) = 0;
    
    I = uint32( find(m ~= 1 & m ~= 0) );
    if ~isempty(I)
      signU = sign(u(I));
      
      % pre-allocate space and augment if needed
      chunk = 7;
      a = zeros(chunk,1);
      c = a;
      b = a;
      a(1,:) = 1;
      c(1,:) = sqrt(m);
      b(1,:) = sqrt(1-m);
      n = uint32( zeros(1,1) );
      i = 1;
      while any(abs(c(i,:)) > tol) % Arithmetic-Geometric Mean of A, B and C
        i = i + 1;
        if i > size(a,1)
          a = [a; zeros(2,1)];
          b = [b; zeros(2,1)];
          c = [c; zeros(2,1)];
        end
        a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
        b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
        c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
        in = uint32( find((abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol)) );
        if ~isempty(in)
          [mi,ni] = size(in);
          n(in) = ones(mi,ni)*(i-1);
        end
      end
      
      mmax = length(I);
      mn = double(max(n));
      phin = zeros(1,mmax);     C  = zeros(1,mmax);
      Cp = C;  e  = uint32(C);  phin(:) = signU.*u(I);
      i = 0;   c2 = c.^2;
      while i < mn % Descending Landen Transformation
        i = i + 1;
        in = uint32(find(n > i));
        if ~isempty(in)
          phin(in) = atan(b(i)./a(i).*tan(phin(in))) + ...
            pi.*ceil(phin(in)/pi - 0.5) + phin(in);
          e(in) = 2.^(i-1) ;
          C(in) = C(in)  + double(e(in(1)))*c2(i);
          Cp(in)= Cp(in) + c(i+1).*sin(phin(in));
        end
      end
      
      Ff = phin ./ (a(mn).*double(e)*2);
      F(I) = Ff.*signU;                        % Incomplete Ell. Int. of the First Kind
      E(I) = (Cp + (1 - 1/2*C) .* Ff).*signU;  % Incomplete Ell. Int. of the Second Kind
    end
    
    % Special cases: m == {0, 1}
    m0 = find(m == 0);
    if ~isempty(m0), F(m0) = u(m0); E(m0) = u(m0); end
    
    m1 = find(m == 1);
    um1 = abs(u(m1));
    if ~isempty(m1)
      N = floor( (um1+pi/2)/pi );
      M = find(um1 < pi/2);
      
      F(m1(M)) = log(tan(pi/4 + u(m1(M))/2));
      F(m1(um1 >= pi/2)) = Inf.*sign(u(m1(um1 >= pi/2)));
      
      E(m1) = ((-1).^N .* sin(um1) + 2*N).*sign(u(m1));
    end
  end
% \end{mfunction}
