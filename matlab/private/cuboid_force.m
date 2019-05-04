%% cuboid_force
%
% Calculate the forces between two parallel cuboid magnets.

% \START

% \begin{mfunction}{cuboid_force}
function force_out = cuboid_force(magnet_fixed,magnet_float,displ)

hsize1 = magnet_fixed.dim(:)/2;
hsize2 = magnet_float.dim(:)/2;
J1 = magnet_fixed.magM;
J2 = magnet_float.magM;

force_out = zeros(size(displ));
force_out = force_out + cuboid_force_x_x(hsize1,hsize2,displ,J1,J2);
force_out = force_out + cuboid_force_x_y(hsize1,hsize2,displ,J1,J2);
force_out = force_out + cuboid_force_x_z(hsize1,hsize2,displ,J1,J2);
force_out = force_out + cuboid_force_y_x(hsize1,hsize2,displ,J1,J2);
force_out = force_out + cuboid_force_y_y(hsize1,hsize2,displ,J1,J2);
force_out = force_out + cuboid_force_y_z(hsize1,hsize2,displ,J1,J2);
force_out = force_out + cuboid_force_z_x(hsize1,hsize2,displ,J1,J2);
force_out = force_out + cuboid_force_z_y(hsize1,hsize2,displ,J1,J2);
force_out = force_out + cuboid_force_z_z(hsize1,hsize2,displ,J1,J2);

end
% \end{mfunction}