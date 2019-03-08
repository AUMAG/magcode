%% Dipole comparison example


%% Parallel magnets example

coil = magnetdefine('type','cylinder','turns',100,'current',1,'dim',[0.02 0.04]);
mag  = magnetdefine('type','cylinder','dim',[0.005 0.005],'magn',1);

N = 50;
displ = [0; 0; 1]*linspace(0.03,0.1,N);

coilmagforce = magnetforces(coil,mag,displ);

dipoleforce = nan(3,N);
for ii = 1:N
  dipoleforce(:,ii) = calcdipoleforce(displ(:,ii),coil.dipolemoment,mag.dipolemoment);
end

figure(1); clf; hold on
plot(displ(3,:),coilmagforce(3,:))
plot(displ(3,:),dipoleforce(3,:))
xlabel('Axial displacement, m')
ylabel('Axial force, N')
set(gca,'box','on','ticklength',[0.02 0.05])

function f = calcdipoleforce(r_ab,m_a,m_b)

f = 3e-7/(norm(r_ab)^4)*(...
      + cross(cross(r_ab,m_a),m_b) + ...
      + cross(cross(r_ab,m_b),m_a) + ...
      - 2*r_ab.*(dot(m_a,m_b)) + ...
      + 5*r_ab.*(dot(cross(r_ab,m_a),cross(r_ab,m_b))) + ...
    0);

end