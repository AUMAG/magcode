%% Dipole comparison example

%% verify

ma = [0;0;100];
mb = [0;0;100];
rab = [0;0;1];
calcdipoleforce(rab,ma,mb)

ma = [0;-100;0];
mb = [0;0;100];
rab = [0;0;1];
calcdipoleforce(rab,ma,mb)

%%

a = 0.02;
coil = magnetdefine('type','cuboid','magn',1,'magdir',[0;0;1],'dim',[a a a]);
mag  = magnetdefine('type','cuboid','magn',1,'magdir',[0;0;1],'dim',[a a a]);

N = 50;
displ = [0; 0; 1]*linspace(5*a,20*a,N);

coilmagforce = magnetforces(coil,mag,displ);

dipoleforce = nan(3,N);
for ii = 1:N
  dipoleforce(:,ii) = calcdipoleforce(displ(:,ii),coil.dipolemoment,mag.dipolemoment);
end

figure(1); clf; hold on
plot(displ(3,:),coilmagforce(3,:))
plot(displ(3,:),dipoleforce(3,:))
legend('Exact','Dipole')
xlabel('Axial displacement, m')
ylabel('Axial force, N')
set(gca,'box','on','ticklength',[0.02 0.05],'xscale','log','yscale','log')

function f = calcdipoleforce(r_ab,m_a,m_b)

f = 3e-7/(norm(r_ab)^4)*(...
      + cross(cross(r_ab,m_a),m_b) + ...
      + cross(cross(r_ab,m_b),m_a) + ...
      - 2*r_ab.*(dot(m_a,m_b)) + ...
      + 5*r_ab.*(dot(cross(r_ab,m_a),cross(r_ab,m_b))) + ...
    0);

end