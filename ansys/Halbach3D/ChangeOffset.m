clear all;close all;clc;
for k=2:2
    n=0+.001*(k-1);       %offset
    halbach3d(n)
    
   k=k+1; 
end
% 
% n=0.02;
% 
% halbach3d(n)