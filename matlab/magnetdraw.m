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
    
    vrtc_p = vrtc_p+pos;
    vrtc_n = vrtc_n+pos;
    
    patch('Faces',faces_p,'Vertices',vrtc_p,'FaceColor',color,'FaceAlpha',1)
    patch('Faces',faces_n,'Vertices',vrtc_n,'FaceColor',color/2,'FaceAlpha',1)
    
    for ii = []%1:size(vrtc_p,1)
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
    Z = h*Z-0.5*h;
    
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
  this_face(isnan(this_face)) = [];
  new_face = [];
  
  Nfaceedges = numel(this_face)-1;
  
  for vv = 1:Nfaceedges
    
    ind1 = vv;
    ind2 = vv+1;
    v1 = this_face(ind1);
    v2 = this_face(ind2);
    
%    fprintf('Face %i, Edge %i-%i\n',ff,v1,v2)
    
    p1 = vrtc_new(v1,:);
    p2 = vrtc_new(v2,:);
    line_vec = p2-p1;
    
    pside1 = faces_pn(ff,ind1);
    pside2 = faces_pn(ff,ind2);
    
    if pside1 == 0 && pside2 == 0
%      disp('Drop this edge.')
    elseif pside1 > 0 && pside2 > 0
%      disp('Keep this edge.')
      new_face = [new_face,v1,v2];
    else
%      disp('Cut this edge.')
      
      s = -dot(normm,p1)/(dot(normm,line_vec));
      c = p1+s*line_vec;
      % plot3(c(1),c(2),c(3),'k.','markersize',20)
      
      assert( s>0 && s<=1 , 'Cut not calculated correctly?')
      
      findc = find(all(abs(vrtc_new-c)<eps,2));
      if numel(findc) == 0
        Nnew = size(vrtc_new,1)+1;
        vrtc_new(Nnew,:) = c;
%        fprintf('Edge cut; new vertex %i created.\n',Nnew)
      else
        if numel(findc) > 1
          findc = findc(1);
        end
        Nnew = findc;
%        fprintf('Edge cut; vertex %i used.\n',Nnew)
      end
      if pside1 > 0
        new_face = [new_face,v1,Nnew];
      else
        new_face = [new_face,Nnew,v2];
      end
      
    end
    
    new_face(isnan(new_face)) = [];
    Nedges_new = numel(new_face)-1;
    
    if Nedges_new > Nedges
      faces_new = [faces_new nan(Nfaces,Nedges_new-Nedges)];
      Nedges = Nedges_new;
    elseif Nedges_new < Nedges
      tmp = nan(1,Nedges+1);
      tmp(1:(Nedges_new+1)) = new_face;
      new_face = tmp;
    end
    faces_new(ff,:) = new_face;
    
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