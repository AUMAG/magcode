
close all;clear all;clc;
for k=1:8
    N=4+(k-1)*4;
   
    halbach2d(N);
    load('polardata.mat','E','t','A4');
    B(1)='b';
    B(2)='r';
    B(3)='k';
    B(4)='g';
    B(5)='b';
    B(6)='r';
    B(7)='k';
    B(8)='g';
    M(k)=N;
    l=linspace(1,1000,1000);
    q(k)=trapz(l,A4(E(l)));
%     plot(N,q,'o');
    %polar(t,A4(E(l))',B(k)); 
%     s{k}=sprintf('%s%g','N=',N);
%     legend(s)
%     hold on;
%     title('Magnetic BSUM between rings');
    k=k+1;    
end
