%% Example of ELLIPKEPI

clear all
close all
clc

format compact
format long g

%% Once-off

disp('Testing singularity')

J1 = rand(1);
J2 = rand(1);
r1 = rand(1);
h1 = rand(1);
r2 = rand(1);
h2 = rand(1);
displ = h1/2+h2/2;

fz1 = cylinder_force_coaxial(J1,J2,r1,r2,h1,h2,displ)
fz1 = cylinder_force_coaxial(J1,J2,r1,r2,h1,h2,[0 displ])


%% Once-off

disp('Testing inputs')

J1 = rand(1);
J2 = rand(1);
r1 = rand(1);
h1 = rand(1);
r2 = rand(1);
h2 = rand(1);
displ = h1+h2;

fz1 = cylinder_force_coaxial(J1,J2,r1,r2,h1,h2,displ)
fz2 = cylinder_force_coaxial(J1,J2,r1,r2,h1,h2,[displ;displ*2])




%% Vectoring
%
% The first iteration doesn't show much improvement between a for-loop and
% a vectorised calculation. Subsequent iterations with vectorised code
% are much faster.

disp('Testing vectorisation')

N = 100;

for jj = 1:2
  
  J1 = rand(1);
  J2 = rand(1);
  r1 = rand(1);
  h1 = rand(1);
  r2 = rand(1);
  h2 = rand(1);
  displ = linspace(0,h1/2+h2/2,N);

  fz1 = nan(N,1);
  
  tic
  for ii = 1:N
    fz1(ii) = cylinder_force_coaxial(J1,J2,r1,r2,h1,h2,displ(ii));
  end
  t1 = toc;
  
  tic
  fz2 = cylinder_force_coaxial(J1,J2,r1,r2,h1,h2,displ);
  t2 = toc;
  
  fprintf('For loop (N=%4i) : %0.9f\n',N,t1)
  fprintf('Vector, no loop   : %0.9f\n',t2)
  
  assert(all(abs(fz1-fz2)<eps),'fz')
  fprintf("%2.2f%% improvement, iteration %i\n\n",100*(t1-t2)/t1,jj)
  
end