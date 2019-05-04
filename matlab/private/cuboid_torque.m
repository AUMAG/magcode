%% cuboid_torque
%
% Calculate the torques between two parallel cuboid magnets.

% \START


% \begin{mfunction}{cuboid_torque}
function torques_out = cuboid_torque(size1,size2,displ,lever,J1,J2)

swap_x_z = @(vec) vec([3 2 1],:);
swap_y_z = @(vec) vec([1 3 2],:);

rotate_z_to_x = @(vec) [  vec(3,:);  vec(2,:); -vec(1,:) ] ; % Ry( 90)
rotate_x_to_z = @(vec) [ -vec(3,:);  vec(2,:);  vec(1,:) ] ; % Ry(-90)

rotate_y_to_z = @(vec) [  vec(1,:); -vec(3,:);  vec(2,:) ] ; % Rx( 90)
rotate_z_to_y = @(vec) [  vec(1,:);  vec(3,:); -vec(2,:) ] ; % Rx(-90)

size1_x = swap_x_z(size1);
size2_x = swap_x_z(size2);
J1_x    = rotate_x_to_z(J1);
J2_x    = rotate_x_to_z(J2);

size1_y = swap_y_z(size1);
size2_y = swap_y_z(size2);
J1_y    = rotate_y_to_z(J1);
J2_y    = rotate_y_to_z(J2);

torque_components = nan([size(displ) 9]);

d_x  = rotate_x_to_z(displ);
d_y  = rotate_y_to_z(displ);
d_z  = displ;

l_x = rotate_x_to_z(lever);
l_y = rotate_y_to_z(lever);
l_z = lever;

torque_components(:,:,9) = cuboid_torque_z_z( size1,size2,d_z,l_z,J1,J2 );

torque_components(:,:,8) = cuboid_torque_z_y( size1,size2,d_z,l_z,J1,J2 );

torque_components(:,:,7) = torques_calc_z_x( size1,size2,d_z,l_z,J1,J2 );

torque_components(:,:,1) = ...
  rotate_z_to_x( cuboid_torque_z_z(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );

torque_components(:,:,2) = ...
  rotate_z_to_x( cuboid_torque_z_y(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );

torque_components(:,:,3) = ...
  rotate_z_to_x( torques_calc_z_x(size1_x,size2_x,d_x,l_x,J1_x,J2_x) );

torque_components(:,:,4) = ...
  rotate_z_to_y( torques_calc_z_x(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );

torque_components(:,:,5) = ...
  rotate_z_to_y( cuboid_torque_z_z(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );

torque_components(:,:,6) = ...
  rotate_z_to_y( cuboid_torque_z_y(size1_y,size2_y,d_y,l_y,J1_y,J2_y) );

torques_out = sum(torque_components,3);
end

% \end{mfunction}


% \begin{mfunction}{torques_calc_z_x}
function calc_out = torques_calc_z_x(size1,size2,offset,lever,J1,J2)

if J1(3)~=0 && J2(1)~=0
  error('Torques cannot be calculated for orthogonal magnets yet.')
end

calc_out = 0*offset;

end
% \end{mfunction}

