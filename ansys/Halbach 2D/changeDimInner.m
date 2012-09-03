%% This script plots a polar figure of halbach2D magnets with varying inner ring magnet dimensions

for k=1:20
%     e = .02;       % Outer magnet radial length
%     f = .02;       % Outer magnet tangential length
%     g = .02;       % Inner magnet radial length
%     h = .02;       % Inner magnet tangential length 

   
    e=.02;
    f=.02;
    g=.02-.0005*(k-1);
    h=g;
    
    halbach2d(e,f,g,h);
    load('polardata.mat','E','t','A4');
    B(1)='b';
    B(2)='r';
    B(3)='k';
    B(4)='g';
    B(5)='m';
    
    l=linspace(1,1000,1000);
    p(k)=A4(E(1));
    v(k)=g;
    %polar(t,A4(E(l))',B(k));
    %s{k}=sprintf('%s%g','innermagnetdimension=',g);
    %legend(s)
    %hold on;
    %title('Magnetic BSUM between rings with varying inner ring dimension');
    k=k+1;    
end
