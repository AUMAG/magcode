function drawmagnet(magnet_struct,pos,color)
% drawmagnets(magnet_struct,pos,color). Draws the magnet in 3D by using the
% same structure used in other magnet scripts.
%
% magnet_struct = magnet data
% pos = position of the magnet (1x3)
% color = color of the magnet (1x3)
%

%% Cylindrical magnets
if length(magnet_struct.dim)==2
    r = magnet_struct.dim(1);
    h = magnet_struct.dim(2);
    n = 100;
    [X,Y,Z] = cylinder(r,n);
    Z(2,:) = h*Z(2,:);
    X = X+pos(1);
    Y = Y+pos(2);
    Z = Z-0.5*h+pos(3);
    % Draw cylinder sides
    for i = 1:(size(X,2)-1);
        patch([X(1,i) X(1,i+1) X(2,i+1) X(2,i)],...
            [Y(1,i) Y(1,i+1) Y(2,i+1) Y(2,i)],...
            [Z(1,i) Z(1,i+1) Z(2,i+1) Z(2,i)],color);
    end
    % Bottom cover
    patch(X(1,:),Y(1,:),Z(1,:),color);
    % Top cover
    patch(X(2,:),Y(2,:),Z(2,:),color);
end

%% Cube magnets
if length(magnet_struct.dim)==3
    hdim = magnet_struct.dim/2;
    
    % Define all points of cube
    % Top plate: no. 1, 2, 3, 4
    % Bottom plate: no. 5, 6, 7, 8
    %
    crd = zeros(8,3);
    crd(1,:) = [-hdim(1) -hdim(2) hdim(3)]+pos;
    crd(2,:) = [hdim(1) -hdim(2) hdim(3)]+pos;
    crd(3,:) = [hdim(1) hdim(2) hdim(3)]+pos;
    crd(4,:) = [-hdim(1) hdim(2) hdim(3)]+pos;
    crd(5,:) = [-hdim(1) -hdim(2) -hdim(3)]+pos;
    crd(6,:) = [hdim(1) -hdim(2) -hdim(3)]+pos;
    crd(7,:) = [hdim(1) hdim(2) -hdim(3)]+pos;
    crd(8,:) = [-hdim(1) hdim(2) -hdim(3)]+pos;
    
    % Draw cube
    %
    % Draw top patch
    patch(crd(1:4,1),crd(1:4,2),crd(1:4,3),color);
    % Draw bottom patch
    patch(crd(5:8,1),crd(5:8,2),crd(5:8,3),color);
    % Draw side patches
    patch(crd([5,6,2,1],1),crd([5,6,2,1],2),crd([5,6,2,1],3),color);
    patch(crd([8,7,3,4],1),crd([8,7,3,4],2),crd([8,7,3,4],3),color);
    patch(crd([6,7,3,2],1),crd([6,7,3,2],2),crd([6,7,3,2],3),color);
    patch(crd([5,8,4,1],1),crd([5,8,4,1],2),crd([5,8,4,1],3),color);
end

end