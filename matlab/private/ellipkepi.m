% \begin{mfunction}{ellipkepi}
% Complete elliptic integrals calculated with the arithmetric-geometric mean
% algorithms contained here: \url{http://dlmf.nist.gov/19.8}.
% Valid for $a<=1$ and $m<=1$.

function [K,E,PI] = ellipkepi(a,m)

a1 = 1;
g1 = sqrt(1-m);
s1 = m;
nn = 0;

p1 = sqrt(1-a);
q1 = 1;
qq = 1;
w1 = 1;

while max(w1(:)) > eps || max(q1(:)) > eps
  
  % Update from previous loop
  a0 = a1;
  g0 = g1;
  p0 = p1;
  q0 = q1;
  
  % for Elliptic I
  a1 = (a0+g0)/2;
  g1 = sqrt(a0.*g0);
  
  % for Elliptic II
  nn = nn + 1;
  c1 = (a0-g0)/2;
  w1 = 2^nn*c1.^2;
  s1 = s1 + w1;
  
  % for Elliptic III
  rr = p0.^2+a0.*g0;
  p1 = rr./(2.*p0);
  q1 = 0.5*q0.*(p0.^2-a0.*g0)./rr;
  qq = qq + q1;
  
end

K  = 1./a1*pi/2;
E  = K.*(1-s1/2);
PI = K.*(1+a./(2-2*a).*qq);

im = find(m == 1);
if ~isempty(im)
  K(im) = inf;
  E(im) = ones(length(im),1);
  PI(im) = inf;
end

end
% \end{mfunction}