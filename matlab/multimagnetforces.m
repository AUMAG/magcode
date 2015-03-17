function F = multimagnetforces(magnet_float, displ, magnets_fixed, pos)
% F = multimagnetforces(magnet_float, displ, magnets_fixed, pos). 
% Calculates the forces acting on a floating magnet with multiple fixed 
% magnets n. Note: only for forces still.
%
% magnet_float  = floating/levitated magnet data
% displ         = displacement of floating magnet wrt its center of gravity
% magnets_fixed = fixed magnets data (nx1 structures)
% pos           = position of the fixed magnet wrt the floating magnet (nx3)
%
% Output F      = Matrix with the forces acting on the floating magnet in
%                 each direction for m different displacements (3xm)
%
% Note: floating magnet has position xyz = [0 0 0];
%

figure
view(3)
grid on

% Draw floating magnet in green color
drawmagnet(magnet_float,[0 0 0],[0 1 0])

% Check input parameters
Nmagnets = length(magnets_fixed);
if size(pos,2)==3
elseif size(pos,1)==3
    pos = pos';
else
    error('Position matrix does not correspond with fixed magnets data');
end

% Calculate forces
for i = 1:Nmagnets
    magnet_fixed = magnets_fixed(i);
    % Draw fixed magnet in red color
    drawmagnet(magnet_fixed,pos(i,:),[1 0 0])
    % Invert position to make it with respect to the fixed magnet
    pos(i,:) = -1*pos(i,:);
    % Check displacement input and add positions
    if size(displ,2)==3
        displ2 = [displ(:,1)+pos(i,1) displ(:,2)+pos(i,2) displ(:,3)+pos(i,3)];
    else
        displ = displ';
        displ2 = [displ(:,1)+pos(i,1) displ(:,2)+pos(i,2) displ(:,3)+pos(i,3)];
    end
    % Calculate the forces
    if i==1
        F = magnetforces(magnet_fixed, magnet_float, displ2);
    else
        F(:,:) = F(:,:) + magnetforces(magnet_fixed, magnet_float, displ2);
    end
end

axis equal
end