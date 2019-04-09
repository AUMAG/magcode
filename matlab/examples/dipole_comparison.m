%% Dipole comparison example

%% Check dpt/cross product vectorisation

a = [1;2;3];
b = [4;5;6];
c = [7;8;9];
t1 = [cross(a,b),cross(a,c)];
t2 = cross([a,a],[b,c]);
assert(all(t1(:)==t2(:)))

t1 = [dot(a,b),dot(a,c)];
t2 = dot([a,a],[b,c]);
assert(all(t1(:)==t2(:)))


%% Verify Yung 1998 example (Table 1)

ma = [0;0;100];
mb = [0;0;100];
rab = [0;0;1];
calcdipoleforcetorque(ma,mb,rab)

ma = [0;-100;0];
mb = [0;0;100];
rab = [0;0;1];
calcdipoleforcetorque(ma,mb,rab)

%%

clc

a = 0.02;
coil = magnetdefine('type','cuboid','magn',+1,'magdir',[0;0;1],'dim',[a 2*a a]);
mag  = magnetdefine('type','cuboid','magn',-1,'magdir',[0;0;1],'dim',[2*a a a]);

N = 50;
displ = [0; 0; 1]*linspace(2*a,5*a,N)+[2*a;a;0];

[coilmagforce,coilmagtorque] = magnetforces(coil,mag,displ,'force','torque');
[dipoleforce,dipoletorque]  = calcdipoleforcetorque(coil.dipolemoment,mag.dipolemoment,displ);

figure(1); clf; hold on
subplot(2,1,1); cla; hold on
h = {0;0;0};
for ii = 1:3
  h{ii} = plot(displ(3,:),coilmagforce(ii,:));
end
for ii = 1:3
  plot(displ(3,:),dipoleforce(ii,:),'--','color',h{ii}.Color)
end
legend('X','Y','Z')
xlabel('Displacement, m')
ylabel('Force, N')
set(gca,'box','on')

subplot(2,1,2); cla; hold on
h = {0;0;0};
for ii = 1:3
  h{ii} = plot(displ(3,:),coilmagtorque(ii,:));
end
for ii = 1:3
  plot(displ(3,:),dipoletorque(ii,:),'--','color',h{ii}.Color)
end
legend('X','Y','Z')
xlabel('Displacement, m')
ylabel('Torque, Nm')
set(gca,'box','on')



%% Dipole functions

function [f,t] = calcdipoleforcetorque(m_a,m_b,r_ab)

% m = 1/(4*pi*1e-7)*Br*V

% ensure dipole moment vectors are replicated to the same length as the
% displacement vector:
N = size(r_ab,2);
if size(m_a,2) == 1
  m_a = repmat(m_a,[1,N]);
  m_b = repmat(m_b,[1,N]);
else
  error('Magnetic dipole moment must be defined as a 3x1 vector')
end

R_ab = sqrt(r_ab(1,:).^2+r_ab(2,:).^2+r_ab(3,:).^2);
ur   = r_ab./R_ab;

drma = dot(ur,m_a);
drmb = dot(ur,m_b);

f = 3e-7./(R_ab.^4).*(...
   + ur.*dot(m_a,m_b) ...
   + m_a.*drmb ...
   + m_b.*drma ...
   - 5*ur.*drma.*drmb ...
  );

t = 1e-7./R_ab.^5.*( ...
  + cross(3*m_b,dot(m_a,r_ab).*r_ab) ...  
  - r_ab.^2.*cross(m_a,m_b) ...
  );

end

