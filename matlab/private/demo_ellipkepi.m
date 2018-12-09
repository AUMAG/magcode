%% Example of ELLIPKEPI

clear all
close all
clc

format compact
format long g

%% Once-off

[K1,E1] = ellipke(0.5);% in-built into Matlab
[K2,E2,P2] = ellipkepi(0.5,0.5);

K = [K1; K2]
E = [E1; E2]
P2


%% Vectoring
%
% The first iteration doesn't show much improvement between a for-loop and
% a vectorised calculation. Subsequent iterations with vectorised code
% are much faster.

N = 1000;

for jj = 1:2
  
  a = rand(N,1);
  m = rand(N,1);
  
  K3 = nan(N,1);
  E3 = nan(N,1);
  P3 = nan(N,1);
  
  tic
  for ii = 1:N
    [K3(ii),E3(ii),P3(ii)] = ellipkepi(a(ii),m(ii));
  end
  t1 = toc;
  
  tic
  [K4,E4,P4] = ellipkepi(a,m);
  t2 = toc;
  
  fprintf('For loop (N=%4i) : %0.9f\n',N,t1)
  fprintf('Vector, no loop   : %0.9f\n',t2)
  
  assert(all(abs(K3-K4)<eps),'K')
  assert(all(abs(E3-E4)<eps),'E')
  assert(all(abs(P3-P4)<eps),'PI, %i of %i inconsistent',sum(abs(P3-P4)>=eps),N)
  
  fprintf("%2.2f%% improvement, iteration %i\n\n",100*(t1-t2)/t1,jj)
  
end