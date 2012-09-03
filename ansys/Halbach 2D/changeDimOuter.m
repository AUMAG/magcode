%% This script plots a polar figure of halbach2D magnets with varying outer ring magnet dimensions
clc;clear all;
for k=1:5

%     e = .03;       % Outer magnet radial length
%     f = .03;       % Outer magnet tangential length
%     g = .02;       % Inner magnet radial length
%     h = .02;       % Inner magnet tangential length 

   
    g=.02;
    h=.02;
    %f=.03-.005*(k-1);
    f=.04-.005*(k-1);
    e=.02;
    
    halbach2d(e,f,g,h);
    load('polardata.mat','E','t','A4');
    B(1)='b';
    B(2)='r';
    B(3)='k';
    B(4)='g';
    B(5)='m';
    p(k)=A4(E(1));
    v(k)=e;
    l=linspace(1,1000,1000);
    polar(t,A4(E(l))',B(k));  
    s{k}=sprintf('%s%g','outermagnetdimension=',f);
    legend(s,'Location','Northeast')
    hold on;
    title('Magnetic BSUM between rings with varying outer ring dimension');
    k=k+1;    
end
