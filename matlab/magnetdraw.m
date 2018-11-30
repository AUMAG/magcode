function magnetdraw(magnet,pos,varargin)
% drawmagnets(magnet_struct,pos,color). Draws the magnet in 3D by using the
% same structure used in other magnet scripts.
%
% magnet_struct = magnet data
% pos = position of the magnet (1x3)
% color = color of the magnet (1x3)
%

if nargin == 0
  disp('No arguments; running MAGNETDRAW demo.')
  demo_magnetdraw;
  return
end

color = [1 0 0.2];

switch magnet.type
  case 'cuboid',  draw_cube(magnet,pos);
  case 'cylinder', draw_cyl(magnet,pos);
  otherwise, error(['Cannot draw magnet of type "',magnet.type,'".'])
end

%% Sub-functions

function draw_cube(magnet,pos)

  pos = transpose(pos(:));
hdim = magnet.dim/2;

vrtc = zeros(8,3);
vrtc(1,:) = [-hdim(1); -hdim(2);  hdim(3)]; % top plate
vrtc(2,:) = [ hdim(1); -hdim(2);  hdim(3)];
vrtc(3,:) = [ hdim(1);  hdim(2);  hdim(3)];
vrtc(4,:) = [-hdim(1);  hdim(2);  hdim(3)];
vrtc(5,:) = [-hdim(1); -hdim(2); -hdim(3)]; % bottom plate
vrtc(6,:) = [ hdim(1); -hdim(2); -hdim(3)];
vrtc(7,:) = [ hdim(1);  hdim(2); -hdim(3)];
vrtc(8,:) = [-hdim(1);  hdim(2); -hdim(3)];

faces(1,:) = [1,2,3,4]; % top
faces(2,:) = [5,6,7,8]; % bottom
faces(3,:) = [5,6,2,1]; % sides
faces(4,:) = [8,7,3,4]; %
faces(5,:) = [6,7,3,2]; %
faces(6,:) = [5,8,4,1]; %

[vrtc_p, faces_p] = split_patches(vrtc,faces,+magnet.magdir);
[vrtc_n, faces_n] = split_patches(vrtc,faces,-magnet.magdir);

patch('Faces',faces_p,'Vertices',vrtc_p+pos,'FaceColor',color)
patch('Faces',faces_n,'Vertices',vrtc_n+pos,'FaceColor',color/2)

end

function draw_cyl(magnet,pos)

pos = transpose(pos(:));

r = magnet.dim(1);
h = magnet.dim(2);
n = 50;

[X,Y,Z] = cylinder(r,n);
Z(2,:) = h*Z(2,:);
X = X;
Y = Y;
Z = Z-0.5*h;

vrtc = [X(:), Y(:), Z(:)];
faces = nan(n,4);
for ii = 1:n
  faces(ii,:) = 2*(ii-1)+[1 3 4 2];
end

% Sides

[vrtc_p, faces_p] = split_patches(vrtc,faces,+magnet.magdir);
[vrtc_n, faces_n] = split_patches(vrtc,faces,-magnet.magdir);

patch('Faces',faces_p,'Vertices',vrtc_p+pos,'FaceColor',color,'EdgeColor','none')
patch('Faces',faces_n,'Vertices',vrtc_n+pos,'FaceColor',color/2,'EdgeColor','none')

% Bottom & Top cover

faces = [1:2:2*(n+1);2:2:2*(n+1)];
[vrtc_p, faces_p] = split_patches(vrtc,faces,+magnet.magdir);
[vrtc_n, faces_n] = split_patches(vrtc,faces,-magnet.magdir);

patch('Faces',faces_p,'Vertices',vrtc_p+pos,'FaceColor',color,'EdgeColor','none')
patch('Faces',faces_n,'Vertices',vrtc_n+pos,'FaceColor',color/2,'EdgeColor','none')

end

end

function [vrtc_new, faces_new] = split_patches(vrtc,faces,norm)

Nfaces = size(faces,1);
Nvrtc  = size(faces,2);

faces_pn = calc_face_vertex_side(faces,vrtc);

% remove faces on the wrong side of the plane
faces(all(faces_pn<0,2),:) = [];
Nfaces = size(faces,1);

faces_pn = calc_face_vertex_side(faces,vrtc);

vrtc_new  = vrtc;
faces_new = faces;

for ff = 1:Nfaces
  line_ind = [faces(ff,:),faces(ff,1)];
  for vv = 1:Nvrtc
    p2 = vrtc(line_ind(vv+1),:);
    p1 = vrtc(line_ind(vv),:);
    line_vec = p2-p1;
    
    s = -dot(norm,p1)/(dot(norm,line_vec));
    c = p1+s*line_vec;
    
    if s>0 && s<=1
      if faces_pn(ff,vv) > 0
        vrtc_new(line_ind(vv+1),:) = c;
      else
        vrtc_new(line_ind(vv  ),:) = c;
      end
    end
    
  end
end

  function faces_pn = calc_face_vertex_side(faces,vrtc)
    
    faces_pn = nan(size(faces));
    
    % calculate whether face vertices are pos or neg w.r.t. the plane
    for fff = 1:Nfaces
      for vvv = 1:Nvrtc
        pp = vrtc(faces(fff,vvv),:);
        faces_pn(fff,vvv) = sign(dot(norm,pp));
      end
    end

  end

end