function magnetdraw(magnet,pos,varargin)
% drawmagnets(magnet_struct,pos,color). Draws the magnet in 3D by using the
% same structure used in other magnet scripts.
%
% magnet_struct = magnet data
% pos = position of the magnet (1x3)
% color = color of the magnet (1x3)
%

color = [0.7 0 0.2];

switch magnet.type
  case 'cuboid',  draw_cube(magnet,pos);
  case 'cylinder', draw_cyl(magnet,pos);
  otherwise, error(['Cannot draw magnet of type "',magnet.type,'".'])
end

%% Sub-functions

function draw_cube(magnet,pos)

    hdim = magnet.dim/2;
    
    % Define vertices of cube
    % Top plate:    no. 1, 2, 3, 4
    % Bottom plate: no. 5, 6, 7, 8
    
    vrtc = zeros(8,3);
    vrtc(1,:) = [-hdim(1); -hdim(2);  hdim(3)] + pos;
    vrtc(2,:) = [ hdim(1); -hdim(2);  hdim(3)] + pos;
    vrtc(3,:) = [ hdim(1);  hdim(2);  hdim(3)] + pos;
    vrtc(4,:) = [-hdim(1);  hdim(2);  hdim(3)] + pos;
    vrtc(5,:) = [-hdim(1); -hdim(2); -hdim(3)] + pos;
    vrtc(6,:) = [ hdim(1); -hdim(2); -hdim(3)] + pos;
    vrtc(7,:) = [ hdim(1);  hdim(2); -hdim(3)] + pos;
    vrtc(8,:) = [-hdim(1);  hdim(2); -hdim(3)] + pos;
    
    % Draw top & bottom patches
    patch(vrtc(1:4,1),vrtc(1:4,2),vrtc(1:4,3),color);
    patch(vrtc(5:8,1),vrtc(5:8,2),vrtc(5:8,3),color);
    
    % Draw side patches
    patch(vrtc([5,6,2,1],1),vrtc([5,6,2,1],2),vrtc([5,6,2,1],3),color);
    patch(vrtc([8,7,3,4],1),vrtc([8,7,3,4],2),vrtc([8,7,3,4],3),color);
    patch(vrtc([6,7,3,2],1),vrtc([6,7,3,2],2),vrtc([6,7,3,2],3),color);
    patch(vrtc([5,8,4,1],1),vrtc([5,8,4,1],2),vrtc([5,8,4,1],3),color);

    % TEST: figure(1); clf; magnetdraw(struct('type','cuboid','dim',[0.1 0.2 0.3],'magdir','z'),[0; 0; 0]); view(3); xlabel('x'); ylabel('y'); zlabel('z')
end

function draw_cyl(magnet,pos)

    r = magnet.dim(1);
    h = magnet.dim(2);
    n = 50;
    
    [X,Y,Z] = cylinder(r,n);
    Z(2,:) = h*Z(2,:);
    X = X + pos(1);
    Y = Y + pos(2);
    Z = Z-0.5*h + pos(3);
    
    % Draw cylinder sides
    for ii = 1:n
        patch([X(1,ii) X(1,ii+1) X(2,ii+1) X(2,ii)],...
              [Y(1,ii) Y(1,ii+1) Y(2,ii+1) Y(2,ii)],...
              [Z(1,ii) Z(1,ii+1) Z(2,ii+1) Z(2,ii)],...
              color,'edgecolor','none');
    end
    
    % Bottom & Top cover
    patch(X(1,:),Y(1,:),Z(1,:),color);
    patch(X(2,:),Y(2,:),Z(2,:),color);
    
    % TEST: figure(1); clf; magnetdraw(struct('type','cylinder','dim',[0.1 0.2],'magdir','z'),[0; 0; 0]); view(3); xlabel('x'); ylabel('y'); zlabel('z')
end

end