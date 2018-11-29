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

[vrtc_p, vrtc_n] = split_patches(vrtc,faces,magnet.magdir);

patch('Faces',faces,'Vertices',vrtc_p+pos,'FaceColor',color)
patch('Faces',faces,'Vertices',vrtc_n+pos,'FaceColor',color/2)

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

[vrtc_p, vrtc_n, faces_p, faces_n] = split_patches(vrtc,faces,magnet.magdir);

patch('Faces',faces_p,'Vertices',vrtc_p+pos,'FaceColor',color,'EdgeColor','none')
patch('Faces',faces_n,'Vertices',vrtc_n+pos,'FaceColor',color/2,'EdgeColor','none')

% Bottom & Top cover

faces = [1:2:2*(n+1);2:2:2*(n+1)];
[vrtc_p, vrtc_n, faces_p, faces_n] = split_patches(vrtc,faces,magnet.magdir);

patch('Faces',faces_p,'Vertices',vrtc_p+pos,'FaceColor',color,'EdgeColor','none')
patch('Faces',faces_n,'Vertices',vrtc_n+pos,'FaceColor',color/2,'EdgeColor','none')

end

end

function [vrtc_p, vrtc_n, faces_p, faces_n] = split_patches(vrtc,faces,norm)

Nfaces = size(faces,1);
Nedges = size(faces,2);

vrtc_p = vrtc;
vrtc_n = vrtc;
faces_p = faces;
faces_n = faces;

for ff = 1:Nfaces
  line_ind = [faces(ff,:),faces(ff,1)];
  for ll = 1:Nedges
    p2 = vrtc(line_ind(ll+1),:);
    p1 = vrtc(line_ind(ll),:);
    line_vec = p2-p1;
    
    s = -dot(norm,p1)/(dot(norm,line_vec));
    c = p1+s*line_vec;
    pn = sign(dot(norm,p1));
    
    if s>0 && s<=1
      if pn > 0
        vrtc_n(line_ind(ll  ),:) = c;
        vrtc_p(line_ind(ll+1),:) = c;
      else
        vrtc_n(line_ind(ll+1),:) = c;
        vrtc_p(line_ind(ll  ),:) = c;
      end
    else
      if pn > 0
%        faces_n(ff,:) = NaN;
      elseif pn < 0
%        faces_p(ff,:) = NaN;
      end
    end
    
  end
end

end