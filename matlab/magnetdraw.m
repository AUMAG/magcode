function magnetdraw(magnet,arg_pos,varargin)
% MAGNETDRAW(mag,pos) Draws a magnet in 3D.
%
% # Mandatory arguments
%
% mag = magnet created by magnetdefine()
% pos = position of the magnet (1x3)
%
% # Optional arguments
%
% 'color', color1 : color of the magnet (1x3)
% 'color2',color2 : color of the "alternate face" of the magnet (1x3)
% 'alpha', alpha  : transparency of the magnet
%
% By default color2 is 0.5*color1.
%

if nargin == 0
  disp('No arguments; running MAGNETDRAW demo.')
  demo_magnetdraw;
  return
elseif nargin == 1
  arg_pos = [0;0;0];
end

if ~isfield(magnet,'fndefined')
  magnet = magnetdefine(magnet);
end

if isfield(magnet,'alpha')
  default_alpha = magnet.alpha;
else
  default_alpha  = 0.8;
end
if isfield(magnet,'color')
  default_color = magnet.color;
else
  default_color  = [1 0 0.2];
end
if isfield(magnet,'color2')
  default_color2 = magnet.color2;
else
  default_color2 = [nan nan nan];
end

p = inputParser;
p.addRequired('pos',@(x) size(x,1)==3 );
p.addParameter('color',default_color);
p.addParameter('color2',default_color2);
p.addParameter('alpha',default_alpha);
p.parse(arg_pos,varargin{:});

pos    = p.Results.pos;
Ndispl = size(pos,2);
color1 = p.Results.color;
color2 = p.Results.color2;
alpha  = p.Results.alpha;
if all(isnan(color2))
  color2 = color1/2;
end
patch_opts1 = {'FaceColor',color1,'FaceAlpha',alpha,'EdgeColor',color1/4};
patch_opts2 = {'FaceColor',color2,'FaceAlpha',alpha,'EdgeColor',color1/4};
patch_opts3 = {'FaceColor',color1,'FaceAlpha',alpha,'EdgeColor','none'};
patch_opts4 = {'FaceColor',color2,'FaceAlpha',alpha,'EdgeColor','none'};

for ii = 1:Ndispl
  switch magnet.type
    case 'cuboid',  draw_cube(magnet,pos(:,ii));
    case 'cylinder', draw_cyl(magnet,pos(:,ii));
    otherwise, error(['Cannot draw magnet of type "',magnet.type,'".'])
  end
end

%% Sub-functions

  function draw_cube(magnet,pos)
    
    pos = transpose(pos(:));
    hdim = transpose(magnet.dim(:))/2;
    
    vrtc = [-1 -1 +1; % top plate
            +1 -1 +1;
            +1 +1 +1;
            -1 +1 +1;
            -1 -1 -1; % bottom plate
            +1 -1 -1;
            +1 +1 -1;
            -1 +1 -1].*hdim;
    
    faces = [1,2,3,4;  % top
             5,6,7,8;  % bottom
             5,6,2,1;  % sides
             8,7,3,4;  %
             6,7,3,2;  %
             5,8,4,1]; %
    
    [vrtc_p, faces_p] = split_patches(vrtc,faces,+magnet.magdir);
    [vrtc_n, faces_n] = split_patches(vrtc,faces,-magnet.magdir);
    
    vrtc_p = transpose(magnet.rotation*transpose(vrtc_p) + magnet.position);
    vrtc_n = transpose(magnet.rotation*transpose(vrtc_n) + magnet.position);
    
    vrtc_p = vrtc_p+pos;
    vrtc_n = vrtc_n+pos;
    
    patch('Faces',faces_p,'Vertices',vrtc_p,patch_opts1{:})
    patch('Faces',faces_n,'Vertices',vrtc_n,patch_opts2{:})
    
  end

  function draw_cyl(magnet,pos)
    
    pos = transpose(pos(:));
    
    r = magnet.dim(1);
    h = magnet.dim(2);
    n = 50;
    
    [X,Y,Z] = cylinder(r,n);
    vrtc = [X(:), Y(:), h*(Z(:)-0.5)];

    faces = nan(n,4);
    for ii = 1:n
      faces(ii,:) = 2*(ii-1)+[1 3 4 2];
    end
    
    % Sides
    [vrtc_p, faces_p] = split_patches(vrtc,faces,+magnet.magdir);
    [vrtc_n, faces_n] = split_patches(vrtc,faces,-magnet.magdir);
    
    vrtc_p = transpose(magnet.rotation*transpose(vrtc_p) + magnet.position);
    vrtc_n = transpose(magnet.rotation*transpose(vrtc_n) + magnet.position);

    patch('Faces',faces_p,'Vertices',vrtc_p+pos,patch_opts3{:})
    patch('Faces',faces_n,'Vertices',vrtc_n+pos,patch_opts4{:})
    
    % Bottom & Top cover
    faces = [1:2:2*(n+1);2:2:2*(n+1)];
    [vrtc_p, faces_p] = split_patches(vrtc,faces,+magnet.magdir);
    [vrtc_n, faces_n] = split_patches(vrtc,faces,-magnet.magdir);

    vrtc_p = transpose(magnet.rotation*transpose(vrtc_p) + magnet.position);
    vrtc_n = transpose(magnet.rotation*transpose(vrtc_n) + magnet.position);

    patch('Faces',faces_p,'Vertices',vrtc_p+pos,patch_opts1{:})
    patch('Faces',faces_n,'Vertices',vrtc_n+pos,patch_opts2{:})
    
  end

end


function [vrtc_new, faces_new] = split_patches(vrtc,faces,normm)

vrtc_new  = vrtc;
faces_new = faces(:,[1:end,1]);

% remove faces on the wrong side of the plane
faces_pn = calc_face_vertex_side(faces_new,vrtc_new,normm);
faces_new(all(faces_pn<1,2),:) = [];
faces_pn = calc_face_vertex_side(faces_new,vrtc_new,normm);

Nfaces = size(faces_new,1);
Nedges = size(faces_new,2)-1;

for ff = 1:Nfaces
  
  this_face = faces_new(ff,:);
  this_face(isnan(this_face)) = [];
  new_face = [];
  
  for vv = 1:numel(this_face)-1
    
    ind1 = vv;
    ind2 = vv+1;
    v1 = this_face(ind1);
    v2 = this_face(ind2);
    p1 = vrtc_new(v1,:);
    p2 = vrtc_new(v2,:);
    pside1 = faces_pn(ff,ind1);
    pside2 = faces_pn(ff,ind2);
%    fprintf('Face %i, Edge %i-%i\n',ff,v1,v2)
    
    if pside1 < 0 && pside2 < 0    %disp('Drop this edge.')
    elseif pside1 > 0 && pside2 > 0  %disp('Keep this edge.')
      new_face = [new_face,v1,v2];
    else                             %disp('Cut this edge.')
      
      s = -dot(normm,p1)/(dot(normm,p2-p1));
      c = p1+s*(p2-p1); % plot3(c(1),c(2),c(3),'k.','markersize',20)
      assert( s>0 && s<=1 , 'Cut not calculated correctly?')
      
      findc = find(all(abs(vrtc_new-c)<eps,2));
      if numel(findc) == 0
        Nnew = size(vrtc_new,1)+1;
        vrtc_new(Nnew,:) = c;
      else
        Nnew = findc(1);
      end
      if pside1 > 0
        new_face = [new_face,v1,Nnew];
      else
        new_face = [new_face,Nnew,v2];
      end
      
    end
    
    % Combine new face back into face vertices array
    
    new_face(isnan(new_face)) = [];
    Nedges_new = numel(new_face)-1;
    
    if Nedges_new > Nedges
      % enlarge face vertices array with nan to fit face
      faces_new = [faces_new, nan(Nfaces,Nedges_new-Nedges)];
      Nedges = Nedges_new;
    elseif Nedges_new < Nedges
      % pad face with nan to fit into face vertices array
      tmp = nan(1,Nedges+1);
      tmp(1:(Nedges_new+1)) = new_face;
      new_face = tmp;
    end
    faces_new(ff,:) = new_face;
    
  end
  
end

% for ii = 1:size(vrtc_new,1)
%   plot3(vrtc_new(ii,1),vrtc_new(ii,2),vrtc_new(ii,3),'.r','markersize',20)
%   text(vrtc_new(ii,1),vrtc_new(ii,2),vrtc_new(ii,3),num2str(ii),'color','red','fontsize',20)
% end

end



function faces_pn = calc_face_vertex_side(faces,vrtc,normm)

Nfaces = size(faces,1);
Nedges = size(faces,2);

faces_pn = nan(size(faces));

% calculate whether face vertices are pos or neg w.r.t. the plane
for fff = 1:Nfaces
  for vvv = 1:Nedges
    pp = vrtc(faces(fff,vvv),:);
    faces_pn(fff,vvv) = sign(dot(normm,pp));
  end
end

end