%% ellipkepi
%
% Calculate complete integrals of the first three kinds.
% If only the first two are needed, use the built-in Matlab function |ellipke| instead.


% \START
% \begin{mfunction}{ellipkepi}
% Complete elliptic integrals calculated with the arithmetric-geometric mean
% algorithms contained here: \url{http://dlmf.nist.gov/19.8}.
% Valid for $0\le a\le 1$ and $0\le m\le 1$.

function [K,E,PI] = ellipkepi(a,m)

a1 = 1;
g1 = sqrt(1-m);
p1 = sqrt(1-a);
q1 = 1;
w1 = 1;

nn = 0;
qq = 1;
ww = m;

while max(abs(w1(:))) > eps || max(abs(q1(:))) > eps

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
  d1 = (a0-g0)/2;
  w1 = 2^nn*d1.^2;
  ww = ww + w1;

  % for Elliptic III
  rr = p0.^2+a0.*g0;
  p1 = rr./p0/2;
  q1 = q0.*(p0.^2-a0.*g0)./rr/2;
  qq = qq + q1;

end

K  = 1./a1*pi/2;
E  = K.*(1-ww/2);
PI = K.*(1+a./(2-2*a).*qq);

im = find(m == 1);
if ~isempty(im)
  K(im) = inf;
  E(im) = ones(length(im),1);
  PI(im) = inf;
end

end

% \end{mfunction}
