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
%[vrtc_n, faces_n] = split_patches(vrtc,faces,-magnet.magdir);

vrtc_p = vrtc_p+pos;
%vrtc_n = vrtc_n+pos;

NN = size(faces_p,1);
for ii = 1:NN
patch('Faces',faces_p(ii,:),'Vertices',vrtc_p,'FaceColor',[ii/NN 0 (NN-ii)/NN],'FaceAlpha',1)
end
%patch('Faces',faces_p,'Vertices',vrtc_p,'FaceColor',color,'FaceAlpha',0.5)
%patch('Faces',faces_n,'Vertices',vrtc_n,'FaceColor',color/2,'FaceAlpha',0.5)

for ii = 1:size(vrtc_p,1)
  plot3(vrtc_p(ii,1),vrtc_p(ii,2),vrtc_p(ii,3),'.r','markersize',20)
  text(vrtc_p(ii,1),vrtc_p(ii,2),vrtc_p(ii,3),num2str(ii),'color','red','fontsize',20)
end
for ii = []%1:size(vrtc_n,1)
  plot3(vrtc_n(ii,1),vrtc_n(ii,2),vrtc_n(ii,3),'.b','markersize',20)
  text(vrtc_n(ii,1),vrtc_n(ii,2),vrtc_n(ii,3),num2str(ii),'color','blue','fontsize',20)
end

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

function [vrtc_new, faces_new] = split_patches(vrtc,faces,normm)

vrtc_new  = vrtc;
faces_new = faces(:,[1:end,1]);

% remove faces on the wrong side of the plane
faces_pn = calc_face_vertex_side(faces_new,vrtc_new);
faces_new(all(faces_pn<1,2),:) = [];
faces_pn = calc_face_vertex_side(faces_new,vrtc_new);

Nfaces = size(faces_new,1);
Nedges = size(faces_new,2)-1;

for ff = 1:Nfaces
  
  this_face = faces_new(ff,:);
  new_face = [];
  vrtc_to_delete = [];
  faces_update = nan(1,Nfaces-1);
  
  for vv = 1:Nedges
    
    v1 = vv;
    v2 = vv+1;
    faces_update(vv) = faces_new(ff,vv);
    
    fprintf('Face %i, Edge %i-%i\n',ff,this_face(v1),this_face(v2))

    p1 = vrtc_new(this_face(v1),:);
    p2 = vrtc_new(this_face(v2),:);
    line_vec = p2-p1;
    
    pside1 = faces_pn(ff,v1);
    pside2 = faces_pn(ff,v2);
    
    if pside1 == 0 && pside2 == 0
      disp('Delete this edge.')
      vrtc_to_delete = [vrtc_to_delete,faces_new(ff,v1),faces_new(ff,v2)];
    elseif pside1 > 0 && pside2 > 0
      disp('Keep this edge.')
      new_face = [new_face,v1,v2];
    else
      disp('Cut this edge.')
      
      s = -dot(normm,p1)/(dot(normm,line_vec));
      fprintf('s: %2.2f\n',s)
      c = p1+s*line_vec;
      plot3(c(1),c(2),c(3),'k.','markersize',20)
      
      assert( s>0 && s<=1 , 'Cut not calculated correctly?')
      vvnew = vv+faces_pn(ff,vv);
      Nnew = size(vrtc_new,1)+1;
      vrtc_new(Nnew,:) = c;
      faces_update(vvnew) = Nnew;
      fprintf('Edge cut; new vertex %i created.\n',Nnew)
      
  end
  
  
end
faces_new(ff,:) = faces_update;
  for ii = 1:numel(vrtc_to_delete)
    ind = faces_new(ff,:) == vrtc_to_delete(ii);
    faces_new(ff,ind) = NaN;
  end
end

  function faces_pn = calc_face_vertex_side(faces,vrtc)
    
    Nfaces = size(faces,1);
    Nedges = size(faces,2);
    
    faces_pn = nan(size(faces));
    
    % calculate whether face vertices are pos or neg w.r.t. the plane
    for fff = 1:Nfaces
      for vvv = 1:Nedges
        pp = vrtc(faces(fff,vvv),:);
        faces_pn(fff,vvv) = (1+sign(dot(normm,pp)))/2;
      end
    end

  end

end