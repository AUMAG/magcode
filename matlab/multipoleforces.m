
function [varargout] = multipoleforces(fixed_array, float_array, displ, varargin)


%% MULTIPOLEFORCES  Calculate forces between two multipole arrays of magnets
%
% Finish this off later. Please read the PDF documentation instead for now.
%



debug_disp = @(str) disp([]);
calc_force_bool = false;
calc_stiffness_bool = false;
calc_torque_bool = false;

% Undefined calculation flags for the three directions:
calc_xyz = [-1; -1; -1];

for ii = 1:length(varargin)
  switch varargin{ii}
    case 'debug',      debug_disp = @(str) disp(str);
    case 'force',      calc_force_bool     = true;
    case 'stiffness',  calc_stiffness_bool = true;
    case 'torque',     calc_torque_bool    = true;
    case 'x',  calc_xyz(1) = 1;
    case 'y',  calc_xyz(2) = 1;
    case 'z',  calc_xyz(3) = 1;
    otherwise
      error(['Unknown calculation option ''',varargin{ii},''''])
  end
end

% If none of |'x'|, |'y'|, |'z'| are specified, calculate all.
if all( calc_xyz == -1 )
  calc_xyz = [1; 1; 1];
end

calc_xyz( calc_xyz == -1 ) = 0;

if ~calc_force_bool && ~calc_stiffness_bool && ~calc_torque_bool
  varargin{end+1} = 'force';
  calc_force_bool = true;
end


if size(displ,1) == 3
  % all good
elseif size(displ,2) == 3
  displ = transpose(displ);
else
  error(['Displacements matrix should be of size (3, D) ',...
         'where D is the number of displacements.'])
end

Ndispl = size(displ,2);

if calc_force_bool
  forces_out = repmat(NaN,[3 Ndispl]);
end

if calc_stiffness_bool
  stiffnesses_out = repmat(NaN,[3 Ndispl]);
end

if calc_torque_bool
  torques_out = repmat(NaN,[3 Ndispl]);
end


part = @(x,y) x(y);

fixed_array = complete_array_from_input(fixed_array);
float_array = complete_array_from_input(float_array);

if calc_force_bool
  array_forces = repmat(NaN,[3 Ndispl fixed_array.total float_array.total]);
end

if calc_stiffness_bool
  array_stiffnesses = repmat(NaN,[3 Ndispl fixed_array.total float_array.total]);
end

displ_from_array_corners = displ ...
  + repmat(fixed_array.size/2,[1 Ndispl]) ...
  - repmat(float_array.size/2,[1 Ndispl]);


for ii = 1:fixed_array.total

  fixed_magnet = struct(...
        'dim',    fixed_array.dim(ii,:), ...
        'magn',   fixed_array.magn(ii), ...
        'magdir', fixed_array.magdir(ii,:) ...
  );

  for jj = 1:float_array.total

    float_magnet = struct(...
      'dim',    float_array.dim(jj,:), ...
      'magn',   float_array.magn(jj), ...
      'magdir', float_array.magdir(jj,:) ...
    );

    mag_displ = displ_from_array_corners ...
                  - repmat(fixed_array.magloc(ii,:)',[1 Ndispl]) ...
                  + repmat(float_array.magloc(jj,:)',[1 Ndispl]) ;

    if calc_force_bool && ~calc_stiffness_bool
      array_forces(:,:,ii,jj) = ...
          magnetforces(fixed_magnet, float_magnet, mag_displ,varargin{:});
    elseif calc_stiffness_bool && ~calc_force_bool
      array_stiffnesses(:,:,ii,jj) = ...
          magnetforces(fixed_magnet, float_magnet, mag_displ,varargin{:});
    else
      [array_forces(:,:,ii,jj) array_stiffnesses(:,:,ii,jj)] = ...
          magnetforces(fixed_magnet, float_magnet, mag_displ,varargin{:});
    end

  end
end

if calc_force_bool
  forces_out = sum(sum(array_forces,4),3);
end

if calc_stiffness_bool
  stiffnesses_out = sum(sum(array_stiffnesses,4),3);
end


varargout = {};

for ii = 1:length(varargin)
  switch varargin{ii}
    case 'force'
      varargout{end+1} = forces_out;

    case 'stiffness'
      varargout{end+1} = stiffnesses_out;

    case 'torque'
      varargout{end+1} = torques_out;
  end
end





function array = complete_array_from_input(array)

if ~isfield(array,'type')
  array.type = 'generic';
end


if ~isfield(array,'face')
  array.face = 'undefined';
end

linear_index = 0;
planar_index = [0 0];

switch array.type
  case 'generic'
  case 'linear',            linear_index = 1;
  case 'linear-quasi',        linear_index = 1;
  case 'planar',            planar_index = [1 2];
  case 'quasi-halbach',     planar_index = [1 2];
  case 'patchwork',         planar_index = [1 2];
  otherwise
    error(['Unknown array type ''',array.type,'''.'])
end

if ~isequal(array.type,'generic')
  if linear_index == 1
    if ~isfield(array,'align')
      array.align = 'x';
    end
    switch array.align
      case 'x', linear_index = 1;
      case 'y', linear_index = 2;
      case 'z', linear_index = 3;
    otherwise
      error('Alignment for linear array must be ''x'', ''y'', or ''z''.')
    end
  else
    if ~isfield(array,'align')
      array.align = 'xy';
    end
    switch array.align
      case 'xy', planar_index = [1 2];
      case 'yz', planar_index = [2 3];
      case 'xz', planar_index = [1 3];
    otherwise
      error('Alignment for planar array must be ''xy'', ''yz'', or ''xz''.')
    end
  end
end

switch array.face
  case {'+x','-x'},   facing_index = 1;
  case {'+y','-y'},   facing_index = 2;
  case {'up','down'}, facing_index = 3;
  case {'+z','-z'},   facing_index = 3;
  case 'undefined',   facing_index = 0;
end

if linear_index ~= 0
  if linear_index == facing_index
    error('Arrays cannot face into their alignment direction.')
  end
elseif ~isequal( planar_index, [0 0] )
  if any( planar_index == facing_index )
    error('Planar-type arrays can only face into their orthogonal direction')
  end
end


switch array.type
  case 'linear'

array = extrapolate_variables(array);

array.mcount = ones(1,3);
array.mcount(linear_index) = array.Nmag;

  case 'linear-quasi'


if isfield(array,'ratio') && isfield(array,'mlength')
  error('Cannot specify both ''ratio'' and ''mlength''.')
elseif ~isfield(array,'ratio') && ~isfield(array,'mlength')
  error('Must specify either ''ratio'' or ''mlength''.')
end


array.Nmag_per_wave = 4;
array.magdir_rotate = 90;

if isfield(array,'Nwaves')
  array.Nmag = array.Nmag_per_wave*array.Nwaves+1;
else
  error('''Nwaves'' must be specified.')
end

if isfield(array,'mlength')
  if numel(array.mlength) ~=2
    error('''mlength'' must have length two for linear-quasi arrays.')
  end
  array.ratio = array.mlength(2)/array.mlength(1);
else
  if isfield(array,'length')
    array.mlength(1) = 2*array.length/(array.Nmag*(1+array.ratio)+1-array.ratio);
    array.mlength(2) = array.mlength(1)*array.ratio;
  else
    error('''length'' must be specified.')
  end
end

array.mcount = ones(1,3);
array.mcount(linear_index) = array.Nmag;

array.msize = repmat(NaN,[array.mcount 3]);

[sindex_x sindex_y sindex_z] = ...
  meshgrid(1:array.mcount(1), 1:array.mcount(2), 1:array.mcount(3));

%% Because the array is linear, the |sindex| terms will be linear also.

all_indices = [1 1 1];
all_indices(linear_index) = 0;
all_indices(facing_index) = 0;
width_index = find(all_indices);

for ii = 1:array.Nmag
  array.msize(sindex_x(ii),sindex_y(ii),sindex_z(ii),linear_index) = ...
    array.mlength(mod(ii-1,2)+1);
  array.msize(sindex_x(ii),sindex_y(ii),sindex_z(ii),facing_index) = ...
    array.height;
  array.msize(sindex_x(ii),sindex_y(ii),sindex_z(ii),width_index) = ...
    array.width;
end


  case 'planar'

if isfield(array,'length')
  if length(array.length) == 1
    if isfield(array,'width')
      array.length = [ array.length array.width ];
    else
      array.length = [ array.length array.length ];
    end
  end
end

if isfield(array,'mlength')
  if length(array.mlength) == 1
    if isfield(array.mwidth)
      array.mlength = [ array.mlength array.mwidth ];
    else
      array.mlength = [ array.mlength array.mlength ];
    end
  end
end

var_names = {'length','mlength','wavelength','Nwaves',...
             'Nmag','Nmag_per_wave','magdir_rotate'};

tmp_array1 = struct();
tmp_array2 = struct();
var_index = zeros(size(var_names));

for iii = 1:length(var_names)
  if isfield(array,var_names(iii))
    tmp_array1.(var_names{iii}) = array.(var_names{iii})(1);
    tmp_array2.(var_names{iii}) = array.(var_names{iii})(end);
  else
    var_index(iii) = 1;
  end
end

tmp_array1 = extrapolate_variables(tmp_array1);
tmp_array2 = extrapolate_variables(tmp_array2);

for iii = find(var_index)
  array.(var_names{iii}) = [tmp_array1.(var_names{iii}) tmp_array2.(var_names{iii})];
end

array.width  = array.length(2);
array.length = array.length(1);

array.mwidth  = array.mlength(2);
array.mlength = array.mlength(1);

array.mcount = ones(1,3);
array.mcount(planar_index) = array.Nmag;

  case 'quasi-halbach'

if isfield(array,'mcount')
  if numel(array.mcount) ~=3
    error('''mcount'' must always have three elements.')
  end
elseif isfield(array,'Nwaves')
  if numel(array.Nwaves) > 2
    error('''Nwaves'' must have one or two elements only.')
  end
  array.mcount(facing_index) = 1;
  array.mcount(planar_index) = 4*array.Nwaves+1;
elseif isfield(array,'Nmag')
  if numel(array.Nmag) > 2
    error('''Nmag'' must have one or two elements only.')
  end
  array.mcount(facing_index) = 1;
  array.mcount(planar_index) = array.Nmag;
else
  error('Must specify the number of magnets (''mcount'' or ''Nmag'') or wavelengths (''Nwaves'')')
end

  case 'patchwork'

if isfield(array,'mcount')
  if numel(array.mcount) ~=3
    error('''mcount'' must always have three elements.')
  end
elseif isfield(array,'Nmag')
  if numel(array.Nmag) > 2
    error('''Nmag'' must have one or two elements only.')
  end
  array.mcount(facing_index) = 1;
  array.mcount(planar_index) = array.Nmag;
else
  error('Must specify the number of magnets (''mcount'' or ''Nmag'')')
end

end


array.total = prod(array.mcount);

if ~isfield(array,'msize')
  array.msize = [NaN NaN NaN];
  if linear_index ~=0
    array.msize(linear_index) = array.mlength;
    array.msize(facing_index) = array.height;
    array.msize(isnan(array.msize)) = array.width;
  elseif ~isequal( planar_index, [0 0] )
    array.msize(planar_index) = [array.mlength array.mwidth];
    array.msize(facing_index) = array.height;
  else
    error('The array property ''msize'' is not defined and I have no way to infer it.')
  end
elseif numel(array.msize) == 1
  array.msize = repmat(array.msize,[3 1]);
end

if numel(array.msize) == 3
  array.msize_array = ...
      repmat(reshape(array.msize,[1 1 1 3]), array.mcount);
else
  if isequal([array.mcount 3],size(array.msize))
    array.msize_array = array.msize;
  else
    error('Magnet size ''msize'' must have three elements (or one element for a cube magnet).')
  end
end
array.dim = reshape(array.msize_array, [array.total 3]);

if ~isfield(array,'mgap')
  array.mgap = [0; 0; 0];
elseif length(array.mgap) == 1
  array.mgap = repmat(array.mgap,[3 1]);
end



if ~isfield(array,'magn')
  array.magn = 1;
end

if length(array.magn) == 1
  array.magn = repmat(array.magn,[array.total 1]);
else
  error('Magnetisation magnitude ''magn'' must be a single value.')
end



if ~isfield(array,'magdir_fn')

  if ~isfield(array,'face')
    array.face = '+z';
  end

  switch array.face
    case {'up','+z','+y','+x'},   magdir_rotate_sign =  1;
    case {'down','-z','-y','-x'}, magdir_rotate_sign = -1;
  end

  if ~isfield(array,'magdir_first')
    array.magdir_first = magdir_rotate_sign*90;
  end

  magdir_fn_comp{1} = @(ii,jj,kk) 0;
  magdir_fn_comp{2} = @(ii,jj,kk) 0;
  magdir_fn_comp{3} = @(ii,jj,kk) 0;

  switch array.type
  case 'linear'
    magdir_theta = @(nn) ...
      array.magdir_first+magdir_rotate_sign*array.magdir_rotate*(nn-1);

    magdir_fn_comp{linear_index} = @(ii,jj,kk) ...
      cosd(magdir_theta(part([ii,jj,kk],linear_index)));

    magdir_fn_comp{facing_index} = @(ii,jj,kk) ...
      sind(magdir_theta(part([ii,jj,kk],linear_index)));

  case 'linear-quasi'

    magdir_theta = @(nn) ...
      array.magdir_first+magdir_rotate_sign*90*(nn-1);

    magdir_fn_comp{linear_index} = @(ii,jj,kk) ...
      cosd(magdir_theta(part([ii,jj,kk],linear_index)));

    magdir_fn_comp{facing_index} = @(ii,jj,kk) ...
      sind(magdir_theta(part([ii,jj,kk],linear_index)));

  case 'planar'

    magdir_theta = @(nn) ...
      array.magdir_first(1)+magdir_rotate_sign*array.magdir_rotate(1)*(nn-1);

    magdir_phi = @(nn) ...
      array.magdir_first(end)+magdir_rotate_sign*array.magdir_rotate(end)*(nn-1);

    magdir_fn_comp{planar_index(1)} = @(ii,jj,kk) ...
      cosd(magdir_theta(part([ii,jj,kk],planar_index(2))));

    magdir_fn_comp{planar_index(2)} = @(ii,jj,kk) ...
      cosd(magdir_phi(part([ii,jj,kk],planar_index(1))));

    magdir_fn_comp{facing_index} = @(ii,jj,kk) ...
      sind(magdir_theta(part([ii,jj,kk],planar_index(1)))) ...
      + sind(magdir_phi(part([ii,jj,kk],planar_index(2))));

  case 'patchwork'

    magdir_fn_comp{planar_index(1)} = @(ii,jj,kk) 0;

    magdir_fn_comp{planar_index(2)} = @(ii,jj,kk) 0;

    magdir_fn_comp{facing_index} = @(ii,jj,kk) ...
      magdir_rotate_sign*(-1)^( ...
            part([ii,jj,kk],planar_index(1)) ...
            + part([ii,jj,kk],planar_index(2)) ...
            + 1 ...
          );

  case 'quasi-halbach'

    magdir_fn_comp{planar_index(1)} = @(ii,jj,kk) ...
      sind(90*part([ii,jj,kk],planar_index(1))) ...
      * cosd(90*part([ii,jj,kk],planar_index(2)));

    magdir_fn_comp{planar_index(2)} = @(ii,jj,kk) ...
      cosd(90*part([ii,jj,kk],planar_index(1))) ...
      * sind(90*part([ii,jj,kk],planar_index(2)));

    magdir_fn_comp{facing_index} = @(ii,jj,kk) ...
      magdir_rotate_sign ...
      * sind(90*part([ii,jj,kk],planar_index(1))) ...
      * sind(90*part([ii,jj,kk],planar_index(2)));

  otherwise
    error('Array property ''magdir_fn'' not defined and I have no way to infer it.')
  end

  array.magdir_fn = @(ii,jj,kk)   ...
    [ magdir_fn_comp{1}(ii,jj,kk) ...
      magdir_fn_comp{2}(ii,jj,kk) ...
      magdir_fn_comp{3}(ii,jj,kk) ];

end





array.magloc = repmat(NaN,[array.total 3]);
array.magdir = array.magloc;
arrat.magloc_array = repmat(NaN,[array.mcount(1) array.mcount(2) array.mcount(3) 3]);

nn = 0;
for iii = 1:array.mcount(1)
  for jjj = 1:array.mcount(2)
    for kkk = 1:array.mcount(3)
      nn = nn + 1;
      array.magdir(nn,:) = array.magdir_fn(iii,jjj,kkk);
    end
  end
end

magsep_x = zeros(size(array.mcount(1)));
magsep_y = zeros(size(array.mcount(2)));
magsep_z = zeros(size(array.mcount(3)));

magsep_x(1) = array.msize_array(1,1,1,1)/2;
magsep_y(1) = array.msize_array(1,1,1,2)/2;
magsep_z(1) = array.msize_array(1,1,1,3)/2;

for iii = 2:array.mcount(1)
  magsep_x(iii) = array.msize_array(iii-1,1,1,1)/2 ...
                + array.msize_array(iii  ,1,1,1)/2 ;
end
for jjj = 2:array.mcount(2)
  magsep_y(jjj) = array.msize_array(1,jjj-1,1,2)/2 ...
                + array.msize_array(1,jjj  ,1,2)/2 ;
end
for kkk = 2:array.mcount(3)
  magsep_z(kkk) = array.msize_array(1,1,kkk-1,3)/2 ...
                + array.msize_array(1,1,kkk  ,3)/2 ;
end

magloc_x = cumsum(magsep_x);
magloc_y = cumsum(magsep_y);
magloc_z = cumsum(magsep_z);

for iii = 1:array.mcount(1)
  for jjj = 1:array.mcount(2)
    for kkk = 1:array.mcount(3)
      array.magloc_array(iii,jjj,kkk,:) = ...
        [magloc_x(iii); magloc_y(jjj); magloc_z(kkk)] ...
        + [iii-1; jjj-1; kkk-1].*array.mgap;
    end
  end
end
array.magloc = reshape(array.magloc_array,[array.total 3]);

array.size = squeeze( array.magloc_array(end,end,end,:) ...
           - array.magloc_array(1,1,1,:) ...
           + array.msize_array(1,1,1,:)/2 ...
           + array.msize_array(end,end,end,:)/2 );

debug_disp('Magnetisation directions')
debug_disp(array.magdir)

debug_disp('Magnet locations:')
debug_disp(array.magloc)


end



function array_out = extrapolate_variables(array)

var_names = {'wavelength','length','Nwaves','mlength',...
             'Nmag','Nmag_per_wave','magdir_rotate'};

if isfield(array,'Nwaves')
  mcount_extra  =  1;
else
  mcount_extra  =  0;
end

if isfield(array,'mlength')
  mlength_adjust = false;
else
  mlength_adjust = true;
end

variables = repmat(NaN,[7 1]);

for iii = 1:length(var_names);
  if isfield(array,var_names(iii))
    variables(iii) = array.(var_names{iii});
  end
end

var_matrix = ...
    [1,  0,  0, -1,  0, -1,  0;
     0,  1,  0, -1, -1,  0,  0;
     0,  0,  1,  0, -1,  1,  0;
     0,  0,  0,  0,  0,  1,  1];

var_results = [0 0 0 log(360)]';
variables = log(variables);

idx = ~isnan(variables);
var_known = var_matrix(:,idx)*variables(idx);
var_calc = var_matrix(:,~idx)\(var_results-var_known);
variables(~idx) = var_calc;
variables = exp(variables);

for iii = 1:length(var_names);
  array.(var_names{iii}) = variables(iii);
end

array.Nmag = round(array.Nmag) + mcount_extra;
array.Nmag_per_wave = round(array.Nmag_per_wave);

if mlength_adjust
  array.mlength  =  array.mlength * (array.Nmag-mcount_extra)/array.Nmag;
end

array_out = array;

end



end

