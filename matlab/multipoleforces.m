  
function [varargout]  =  multipoleforces(fixed_array, float_array, displ, varargin) 
 
  
%% MULTIPOLEFORCES  Calculate forces between two multipole arrays of magnets 
% 
% Finish this off later. 
% 
 
 
 
  
if size(displ,1) == 3 
  % all good 
elseif size(displ,2) == 3 
  displ  =  transpose(displ); 
else 
  error('Displacements matrix should be of size (3, D) where D is the number of displacements.') 
end 
 
Ndispl  =  size(displ,2); 
 
Nvargin  =  length(varargin); 
debug_disp  =  @(str) disp([]); 
calc_force_bool  =  false; 
calc_stiffness_bool  =  false; 
 
for ii  =  1:Nvargin 
  switch varargin{ii} 
    case 'debug' 
      debug_disp  =  @(str) disp(str); 
    case 'force' 
      calc_force_bool  =  true; 
    case 'stiffness' 
      calc_stiffness_bool  =  true; 
    otherwise 
      error(['Unknown calculation option ''',varargin{ii},'''']) 
  end 
end 
 
if ~calc_force_bool && ~calc_stiffness_bool 
  calc_force_bool  =  true; 
end 
 
if calc_force_bool 
  forces_out  =  repmat(NaN,[3 Ndispl]); 
end 
if calc_stiffness_bool 
  stiffnesses_out  =  repmat(NaN,[3 Ndispl]); 
end 
 
 
  
fixed_array  =  complete_array_from_input(fixed_array); 
float_array  =  complete_array_from_input(float_array); 
 
if calc_force_bool 
  array_forces  =  repmat(NaN,[fixed_array.total float_array.total 3 Ndispl]); 
end 
 
if calc_stiffness_bool 
  array_stiffnesses  =  repmat(NaN,[fixed_array.total float_array.total 3 Ndispl]); 
end 
 
displ_from_array_corners  =  displ  ... 
  + repmat(fixed_array.size/2,[1 Ndispl])  ... 
  - repmat(float_array.size/2,[1 Ndispl]); 
 
 
  
for ii  =  1:fixed_array.total 
 
  fixed_magnet  =  struct( ... 
 'dim',    fixed_array.dim(ii,:),  ... 
 'magn',   fixed_array.magn(ii),  ... 
 'magdir', fixed_array.magdir(ii,:)  ... 
  ); 
 
  for jj  =  1:float_array.total 
 
    float_magnet  =  struct( ... 
      'dim',    float_array.dim(jj,:),  ... 
      'magn',   float_array.magn(jj),  ... 
      'magdir', float_array.magdir(jj,:)  ... 
    ); 
 
    for dd  =  1:Ndispl 
 
      mag_displ  =  displ_from_array_corners(:,dd)  ... 
                  - fixed_array.magloc(ii,:)' + float_array.magloc(jj,:)' ; 
 
      if calc_force_bool 
        array_forces(ii,jj,:,dd)  =   ... 
            magnetforces(fixed_magnet, float_magnet, mag_displ,'force'); 
      end 
 
      if calc_stiffness_bool 
        array_stiffnesses(ii,jj,:,dd)  =   ... 
            magnetforces(fixed_magnet, float_magnet, mag_displ,'stiffness'); 
      end 
 
    end 
 
  end 
end 
 
if calc_force_bool 
  forces_out  =  squeeze(sum(sum(array_forces,1),2)); 
end 
 
if calc_stiffness_bool 
  stiffnesses_out  =  squeeze(sum(sum(array_stiffnesses,1),2)); 
end 
 
 
  
varargout{1}  =  forces_out; 
for ii  =  1:Nvargin 
  switch varargin{ii} 
    case 'force' 
      varargout{ii}  =  forces_out; 
    case 'stiffness' 
      varargout{ii}  =  stiffnesses_out; 
  end 
end 
 
 
 
 
  
  
function array_out  =  complete_array_from_input(array) 
 
if ~isfield(array,'type') 
  array.type  =  'generic'; 
end 
 
linear_index  =  0; 
planar_index  =  [0 0]; 
 
switch array.type 
  case 'generic' 
  case 'linear',    linear_index  =  1; 
  case 'linear-x',  linear_index  =  1; 
  case 'linear-y',  linear_index  =  2; 
  case 'linear-z',  linear_index  =  3; 
  case 'planar',    planar_index  =  [1 2]; 
  case 'planar-xy', planar_index  =  [1 2]; 
  case 'planar-yz', planar_index  =  [2 3]; 
  case 'planar-xz', planar_index  =  [1 3]; 
  otherwise 
    error(['Unknown array type ''',array.type,'''.']) 
end 
 
switch array.face 
  case {'+x','-x'},   facing_index  =  1; 
  case {'+y','-y'},   facing_index  =  2; 
  case {'up','down'}, facing_index  =  3; 
  case {'+z','-z'},   facing_index  =  3; 
end 
 
if linear_index ~= 0 
  if linear_index == facing_index 
    error('Arrays cannot face into their alignment direction.') 
  end 
  
array  =  extrapolate_variables(array); 
 
array.mcount  =  ones(1,3); 
array.mcount(linear_index)  =  array.Nmag; 
 
 
elseif ~isequal( planar_index, [0 0] ) 
  if any( planar_index == facing_index ) 
    error('Planar arrays can only face into their orthogonal direction') 
  end 
  
var_names  =  {'length','mlength','wavelength','Nwaves', ... 
             'Nmag','Nmag_per_wave','magdir_rotate'}; 
 
% In the ['length'] direction 
tmp_array1  =  struct(); 
tmp_array2  =  struct(); 
var_index  =  zeros(size(var_names)); 
 
for iii  =  1:length(var_names) 
  if isfield(array,var_names(iii)) 
    tmp_array1.(var_names{iii})  =  array.(var_names{iii})(1); 
    tmp_array2.(var_names{iii})  =  array.(var_names{iii})(2); 
  else 
    var_index(iii)  =  1; 
  end 
end 
 
tmp_array1  =  extrapolate_variables(tmp_array1); 
tmp_array2  =  extrapolate_variables(tmp_array2); 
 
for iii  =  find(var_index) 
  array.(var_names{iii})  =  [tmp_array1.(var_names{iii}) tmp_array2.(var_names{iii})]; 
end 
 
array.depth   =  array.length(2); 
array.length  =  array.length(1); 
 
array.mdepth   =  array.mlength(2); 
array.mlength  =  array.mlength(1); 
 
array.mcount  =  ones(1,3); 
array.mcount(planar_index)  =  array.Nmag; 
 
 
 
end 
 
 
  
array.total  =  prod(array.mcount); 
 
if ~isfield(array,'msize') 
  array.msize  =  [NaN NaN NaN]; 
  if linear_index ~=0 
    array.msize(linear_index)  =  array.mlength; 
    array.msize(facing_index)  =  array.height; 
    array.msize(isnan(array.msize))  =  array.depth; 
  elseif ~isequal( planar_index, [0 0] ) 
    array.msize(planar_index)  =  [array.mlength array.mdepth]; 
    array.msize(facing_index)  =  array.height; 
  else 
    error('The array property ''msize'' is not defined and I have no way to infer it.') 
  end 
elseif numel(array.msize) == 1 
  array.msize  =  repmat(array.msize,[3 1]); 
end 
 
if numel(array.msize) == 3 
  array.msize_array  =   ... 
      repmat(reshape(array.msize,[1 1 1 3]), array.mcount); 
  array.dim  =  reshape(array.msize_array, [array.total 3]); 
else 
  error('Magnet size ''msize'' must have three elements (or one element for a cube magnet).') 
end 
 
if ~isfield(array,'mgap') 
  array.mgap  =  [0; 0; 0]; 
elseif length(array.mgap) == 1 
  array.mgap  =  repmat(array.mgap,[3 1]); 
end 
 
 
 
  
if length(array.magn) == 1 
  array.magn  =  repmat(array.magn,[array.total 1]); 
else 
  error('Magnetisation magnitude ''magn'' must be a single value.') 
end 
 
 
 
 
  
part  =  @(x,y) x(y); 
 
if ~isfield(array,'magdir_fn') 
 
  if ~isfield(array,'face') 
    array.face  =  '+z'; 
  end 
 
  switch array.face 
    case {'up','+z','+y','+x'},   magdir_rotate_sign  =   1; 
    case {'down','-z','-y','-x'}, magdir_rotate_sign  =  -1; 
  end 
 
  magdir_fn_comp{1}  =  @(ii,jj,kk) 0; 
  magdir_fn_comp{2}  =  @(ii,jj,kk) 0; 
  magdir_fn_comp{3}  =  @(ii,jj,kk) 0; 
 
  if linear_index ~= 0 
    magdir_theta  =  @(nn)  ... 
      array.magdir_first+magdir_rotate_sign * array.magdir_rotate * (nn-1); 
 
    magdir_fn_comp{linear_index}  =  @(ii,jj,kk)  ... 
      cosd(magdir_theta(part([ii,jj,kk],linear_index))); 
 
    magdir_fn_comp{facing_index}  =  @(ii,jj,kk)  ... 
      sind(magdir_theta(part([ii,jj,kk],linear_index))); 
 
  elseif ~isequal( planar_index, [0 0] ) 
 
    magdir_theta  =  @(nn)  ... 
      array.magdir_first(1)+magdir_rotate_sign * array.magdir_rotate(1) * (nn-1); 
 
    magdir_phi  =  @(nn)  ... 
      array.magdir_first(2)+magdir_rotate_sign * array.magdir_rotate(2) * (nn-1); 
 
    magdir_fn_comp{planar_index(1)}  =  @(ii,jj,kk)  ... 
      cosd(magdir_theta(part([ii,jj,kk],planar_index(2)))); 
 
    magdir_fn_comp{planar_index(2)}  =  @(ii,jj,kk)  ... 
      cosd(magdir_phi(part([ii,jj,kk],planar_index(1)))); 
 
    magdir_fn_comp{facing_index}  =  @(ii,jj,kk)  ... 
      sind(magdir_theta(part([ii,jj,kk],planar_index(1))))  ... 
      + sind(magdir_phi(part([ii,jj,kk],planar_index(2)))); 
 
  else 
    error('Array property ''magdir_fn'' not defined and I have no way to infer it.') 
  end 
 
  array.magdir_fn  =  @(ii,jj,kk)    ... 
    [ magdir_fn_comp{1}(ii,jj,kk)  ... 
      magdir_fn_comp{2}(ii,jj,kk)  ... 
      magdir_fn_comp{3}(ii,jj,kk) ]; 
 
end 
 
 
 
 
  
array.magloc  =  repmat(NaN,[array.total 3]); 
array.magdir  =  array.magloc; 
arrat.magloc_array  =  repmat(NaN,[array.mcount(1) array.mcount(2) array.mcount(3) 3]); 
 
nn  =  0; 
for iii  =  1:array.mcount(1) 
  for jjj  =  1:array.mcount(2) 
    for kkk  =  1:array.mcount(3) 
      nn  =  nn + 1; 
      array.magdir(nn,:)  =  array.magdir_fn(iii,jjj,kkk); 
    end 
  end 
end 
 
magsep_x  =  zeros(size(array.mcount(1))); 
magsep_y  =  zeros(size(array.mcount(2))); 
magsep_z  =  zeros(size(array.mcount(3))); 
 
magsep_x(1)  =  array.msize_array(1,1,1,1)/2; 
magsep_y(1)  =  array.msize_array(1,1,1,2)/2; 
magsep_z(1)  =  array.msize_array(1,1,1,3)/2; 
 
for iii  =  2:array.mcount(1) 
  magsep_x(iii)  =  array.msize_array(iii-1,1,1,1)/2  ... 
                + array.msize_array(iii  ,1,1,1)/2 ; 
end 
for jjj  =  2:array.mcount(2) 
  magsep_y(jjj)  =  array.msize_array(1,jjj-1,1,2)/2  ... 
                + array.msize_array(1,jjj  ,1,2)/2 ; 
end 
for kkk  =  2:array.mcount(3) 
  magsep_z(kkk)  =  array.msize_array(1,1,kkk-1,3)/2  ... 
                + array.msize_array(1,1,kkk  ,3)/2 ; 
end 
 
magloc_x  =  cumsum(magsep_x); 
magloc_y  =  cumsum(magsep_y); 
magloc_z  =  cumsum(magsep_z); 
 
for iii  =  1:array.mcount(1) 
  for jjj  =  1:array.mcount(2) 
    for kkk  =  1:array.mcount(3) 
      array.magloc_array(iii,jjj,kkk,:)  =   ... 
        [magloc_x(iii); magloc_y(jjj); magloc_z(kkk)]  ... 
        + [iii-1; jjj-1; kkk-1].*array.mgap; 
    end 
  end 
end 
array.magloc  =  reshape(array.magloc_array,[array.total 3]); 
 
array.size  =  squeeze( array.magloc_array(end,end,end,:)  ... 
           - array.magloc_array(1,1,1,:)  ... 
           + array.msize_array(1,1,1,:)/2  ... 
           + array.msize_array(end,end,end,:)/2 ); 
 
debug_disp('Magnetisation directions') 
debug_disp(array.magdir) 
 
debug_disp('Magnet locations:') 
debug_disp(array.magloc) 
 
 
 
array_out  =  array; 
end 
 
 
 
  
function array_out  =  extrapolate_variables(array) 
 
var_names  =  {'wavelength','length','Nwaves','mlength', ... 
             'Nmag','Nmag_per_wave','magdir_rotate'}; 
 
mcount_extra  =  0; 
if isfield(array,'Nwaves') 
  mcount_extra  =  1; 
end 
 
variables  =  repmat(NaN,[7 1]); 
 
for iii  =  1:length(var_names); 
  if isfield(array,var_names(iii)) 
    variables(iii)  =  array.(var_names{iii}); 
  end 
end 
 
var_matrix  =   ... 
    [1,  0,  0, -1,  0, -1,  0; 
     0,  1,  0, -1, -1,  0,  0; 
     0,  0,  1,  0, -1,  1,  0; 
     0,  0,  0,  0,  0,  1,  1]; 
 
var_results  =  [0 0 0 log(360)]'; 
variables  =  log(variables); 
 
idx  =  ~isnan(variables); 
var_known  =  var_matrix(:,idx) * variables(idx); 
var_calc  =  var_matrix(:,~idx)\(var_results-var_known); 
variables(~idx)  =  var_calc; 
variables  =  exp(variables); 
 
for iii  =  1:length(var_names); 
  array.(var_names{iii})  =  variables(iii); 
end 
 
array.Nmag  =  round(array.Nmag) + mcount_extra; 
array.Nmag_per_wave  =  round(array.Nmag_per_wave); 
array.mlength  =  array.mlength * (array.Nmag-mcount_extra)/array.Nmag; 
 
array_out  =  array; 
 
end 
 
 
 
 
 
 
end 
 

