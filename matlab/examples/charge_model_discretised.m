%% Compare with Furlani discrete charge model

%% Setup

clear all

mu0 = 4*pi*1e-7;

magn = 1;
a = 10/1000;
b = 10/1000;
c = 10/1000;

M1 = [0;0;+1];
M2 = [0;0;-1];
magM1 = magn*M1;
magM2 = magn*M2;

normn{1} = [0;0;-1];
normn{2} = [0;0;+1];

Nd = 20;
drange = linspace(1.2*c,2*c,Nd);

%% Akoun & Yonnet calculation

mag1 = magnetdefine('type','cuboid','dim',[a b c],'magdir',M1,'magn',magn);
mag2 = magnetdefine('type','cuboid','dim',[a b c],'magdir',M2,'magn',magn);

disp('Akoun & Yonnet calculation')
tic
F1 = magnetforces(mag1,mag2,[0;0;1]*drange);
toc

FZ1 = F1(3,:); % only compare Z forces

%% Discrete charge model

disp('Furlani calculation')
tic
for dd = 1:Nd
  
  displ = [0; 0; drange(dd)];
      
  N = 20;
  xx = linspace(-a/2,a/2,N+1)+a/N/2; xx(end) = [];
  yy = linspace(-b/2,b/2,N+1)+b/N/2; yy(end) = [];
  zz = [-c/2 c/2];
  
  [x1g,y1g,x2g,y2g] = ndgrid(xx,yy,xx+displ(1),yy+displ(2));
  
  dx = xx(2)-xx(1);
  dy = yy(2)-yy(1);
  dA = dx*dy;
  
  FF = 0;
  
  for k1 = 1:2
    for k2 = 1:2
      
      q1 = dot(magM1,normn{k1})/mu0;
      q2 = dot(magM2,normn{k2})/mu0;
      
      x1 = x1g(:);
      y1 = y1g(:);
      x2 = x2g(:);
      y2 = y2g(:);
      z1 = zz(k1);
      z2 = zz(k2)+displ(3);
      
      distvec = sqrt( (x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2 );
      fx = (x2-x1)./(distvec.^3);
      fy = (y2-y1)./(distvec.^3);
      fz = (z2-z1)./(distvec.^3);
      
      FF = FF + mu0/4/pi*dA*dA*q1*q2*[fx,fy,fz];
      
    end
  end
  
  F2 = transpose(sum(FF));
  FZ2(dd) = F2(3);
    
end
toc

%% Plot

figure(1); clf;

ystretch = 0.2;

subplot(3,1,1)
plot(drange,FZ1,'o-')
xlim([drange(1) drange(end)]+[-1 1]*(drange(2)-drange(1)))
ylim(ylim()+ystretch*[-1 1].*diff(ylim()))
set(gca,'xticklabel',{})
ylabel('Force, N')
title('Akoun & Yonnet model')

subplot(3,1,2)
plot(drange,FZ2,'o-')
xlim([drange(1) drange(end)]+[-1 1]*(drange(2)-drange(1)))
ylim(ylim()+ystretch*[-1 1].*diff(ylim()))
set(gca,'xticklabel',{})
xlabel('Displacement, z direction')
ylabel('Force, N')
title('Furlani discrete charge model')

subplot(3,1,3)
plot(drange,100*(FZ1-FZ2)./FZ1,'ro-')
xlim([drange(1) drange(end)]+[-1 1]*(drange(2)-drange(1)))
ylim(ylim()+ystretch*[-1 1].*diff(ylim()))
xlabel('Displacement, z direction')
ylabel('Error, %')
title('Error')
