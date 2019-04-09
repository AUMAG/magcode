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
clear all
close all

JM = 1;
magJa = [0;0;1];
magJb = [0;0;1];

a = 0.02;
coil = magnetdefine('type','cuboid','magn',JM,'magdir',magJa,'dim',[a 2*a a]);
mag  = magnetdefine('type','cuboid','magn',JM,'magdir',magJb,'dim',[2*a a a]);

N = 50;
displ = [0; 0; 1]*linspace(2*a,5*a,N)+[2*a;a;0];

coil_moment = 1/(4*pi*1e-7)*magJa*coil.volume;
mag_moment  = 1/(4*pi*1e-7)*magJb*mag.volume;

[coilmagforce,coilmagtorque] = magnetforces(coil,mag,displ,'force','torque');
[dipoleforce,dipoletorque]  = calcdipoleforcetorque(coil_moment,mag_moment,displ);

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
title('Cuboid exact (solid) vs dipole (dashed)')

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
title('Cuboid exact (solid) vs dipole (dashed)')



%% Dipole function

function [f,t] = calcdipoleforcetorque(m_a,m_b,r_ab)
% [F,T] = CALCDIPOLEFORCETORQUE(MA,MB,R)
%
% Calculates the force and torque on magnetic dipole MB due to magnetic
% dipole MA with distance vector R.
%
% Magnetic dipole moments can be calculated for a hard magnet using the equation
%           m = 1/(4*pi*1e-7)*Br*V
% where Br is the remanence magnetisation vector in Tesla and V is the magnet volume.
%
% For a coil, the magnetic dipole moment is
%           m = N*I*A*n
% where I is the current, N the number of turns, A the cross-sectional area
% of the current loop(s), and n is the normal vector to the loop(s).
%
% These equations have been adapted from the following pair of papers:
% * Yung 1998:      http://doi.org/10.1155/1998/79537
% * Landecker 1999: http://doi.org/10.1155/1999/97902
%

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
rnorm   = r_ab./R_ab;

dot_rma = dot(rnorm,m_a);
dot_rmb = dot(rnorm,m_b);

f = 3e-7./R_ab.^4.*(...
   + rnorm.*dot(m_a,m_b) ...
   + m_a.*dot_rmb ...
   + m_b.*dot_rma ...
   - 5*rnorm.*dot_rma.*dot_rmb ...
  );

t = 1e-7./R_ab.^3.*( ...
  + cross(3*m_b,dot_rma.*rnorm) ...  
  - cross(m_a,m_b) ...
  );

end

