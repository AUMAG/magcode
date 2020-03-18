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
ww = m;
qq = 1;

while max(abs([w1(:);q1(:)])) > eps

  % Update from previous loop
  a0 = a1;
  g0 = g1;
  p0 = p1;
  q0 = q1;

  p0p0 = p0.^2;
  ag0  = a0.*g0;
  rr   = ag0./p0p0;

  % for Elliptic I
  a1 = (a0+g0)/2;
  g1 = sqrt(ag0);

  % for Elliptic II
  nn = nn + 1;
  w1 = 2^(nn-2)*(a0-g0).^2;
  ww = ww + w1;

  % for Elliptic III
  p1 = p0.*(1+rr)/2;
  q1 = q0.*(1-rr)./(1+rr)/2;
  qq = qq + q1;

end

K  = 1./a1*pi/2;
E  = K.*(1-ww/2);
PI = K.*(1+a./(1-a).*qq/2);

im = find(m == 1);
if ~isempty(im)
  K(im) = inf;
  E(im) = ones(length(im),1);
  PI(im) = inf;
end

end

% \end{mfunction}
