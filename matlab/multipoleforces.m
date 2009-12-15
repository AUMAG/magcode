  
function [varargout]  =  multipoleforces(fixed_array, float_array, displ, varargin) 
 
  
%% MULTIPOLEFORCES  Calculate forces between two multipole arrays of magnets 
% 
% Finish this off later. 
% 
 
 
 
  
Nvargin  =  length(varargin); 
 
if ( Nvargin ~=0 && Nvargin ~= nargout ) 
  error('Must have as many outputs as calculations requested.') 
end 
 
calc_force_bool  =  false; 
calc_stiffness_bool  =  false; 
 
if Nvargin == 0 
  calc_force_bool  =  true; 
else 
  for ii  =  1:Nvargin 
    switch varargin{ii} 
      case 'force' 
        calc_force_bool  =  true; 
      case 'stiffness' 
        calc_stiffness_bool  =  true; 
      otherwise 
        error(['Unknown calculation option ''',varargin{ii},'''']) 
    end 
  end 
end 
 
 
 
  
fixed_array  =  complete_array_from_input(fixed_array); 
float_array  =  complete_array_from_input(float_array); 
 
if calc_force_bool 
  array_forces  =  repmat(NaN,[fixed_array.total float_array.total 3]); 
end 
 
if calc_stiffness_bool 
  array_stiffnesses  =  repmat(NaN,[fixed_array.total float_array.total 3]); 
end 
 
displ  =  reshape(displ,[3 1]); 
 
displ_from_array_corners  =  displ + fixed_array.size/2 - float_array.size/2; 
 
for mm  =  1:fixed_array.total 
 
  fixed_magnet  =  struct( ... 
 'dim',    fixed_array.dim(mm,:),  ... 
 'magn',   fixed_array.magn(mm),  ... 
 'magdir', fixed_array.magdir(mm,:)  ... 
  ); 
 
  for nn  =  1:float_array.total 
 
    mag_displ  =  displ_from_array_corners  ... 
                - fixed_array.magloc(mm,:)' + float_array.magloc(nn,:)' ; 
 
    float_magnet  =  struct( ... 
      'dim',    float_array.dim(nn,:),  ... 
      'magn',   float_array.magn(nn),  ... 
      'magdir', float_array.magdir(nn,:)  ... 
    ); 
 
    if calc_force_bool 
      array_forces(mm,nn,:)  =   ... 
          magnetforces(fixed_magnet, float_magnet, mag_displ,'force'); 
    end 
 
    if calc_stiffness_bool 
      array_stiffnesses(mm,nn,:)  =   ... 
          magnetforces(fixed_magnet, float_magnet, mag_displ,'stiffness'); 
    end 
 
  end 
end 
 
debug_disp('Forces:') 
debug_disp(reshape(array_forces,[],3)) 
 
if calc_force_bool 
  forces_out  =  squeeze(sum(sum(array_forces,1),2)); 
end 
 
if calc_stiffness_bool 
  stiffnesses_out  =  squeeze(sum(sum(array_stiffnesses,1),2)); 
end 
 
 
 
  
if Nvargin == 0 
  varargout{1}  =  forces_out; 
else 
  for ii  =  1:Nvargin 
    switch varargin{ii} 
      case 'force' 
        varargout{ii}  =  forces_out; 
      case 'stiffness' 
        varargout{ii}  =  stiffnesses_out; 
    end 
  end 
end 
 
 
 
 
end 
 
  
  
function array_out  =  complete_array_from_input(array) 
 
if ~isfield(array,'type') 
  array.type  =  'generic'; 
end 
 
switch array.type 
  case 'generic' 
  case 'linear-x',  linear_index  =  1; linear_vect  =  [1 0 0]; 
  case 'linear-y',  linear_index  =  2; linear_vect  =  [0 1 0]; 
  case 'linear-z',  linear_index  =  3; linear_vect  =  [0 0 1]; 
  otherwise 
    error(['Unknown array type ''',array.type,'''.']) 
end 
 
switch array.face 
  case {'+x','-x'},   facing_index  =  1; 
  case {'+y','-y'},   facing_index  =  2; 
  case {'up','down'}, facing_index  =  3; 
  case {'+z','-z'},   facing_index  =  3; 
end 
 
if strncmp(array.type,'linear',6) 
    
var_names  =  {'wavelength','length','Nwaves','mlength', ... 
             'Nmag','Nmag_per_wave','magdir_rotate'}; 
 
mcount_extra  =  0; 
if isfield(array,'Nwaves') 
  mcount_extra  =  1; 
end 
 
variables  =  repmat(NaN,[7 1]); 
 
for ii  =  1:length(var_names); 
  if isfield(array,var_names(ii)) 
    variables(ii)  =  array.(var_names{ii}); 
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
 
for ii  =  1:length(var_names); 
  array.(var_names{ii})  =  variables(ii); 
end 
 
array.Nmag  =  round(array.Nmag) + mcount_extra; 
array.Nmag_per_wave  =  round(array.Nmag_per_wave); 
 
array.mlength  =  array.mlength * (array.Nmag-mcount_extra)/array.Nmag; 
 
array.mcount  =  ones(1,3); 
array.mcount(linear_index)  =  array.Nmag; 
 
 
 
end 
 
  
array.total  =  prod(array.mcount); 
 
if ~isfield(array,'msize') 
  array.msize  =  [NaN NaN NaN]; 
  array.msize(linear_index)  =  array.mlength; 
  array.msize(facing_index)  =  array.height; 
  array.msize(isnan(array.msize))  =  array.depth; 
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
 
  magdir_theta  =  @(nn)  ... 
    array.magdir_first+magdir_rotate_sign * array.magdir_rotate * (nn-1); 
 
  magdir_fn_comp{linear_index}  =  @(ii,jj,kk)  ... 
    cosd(magdir_theta(part([ii,jj,kk],linear_index))); 
 
  magdir_fn_comp{facing_index}  =  @(ii,jj,kk)  ... 
    sind(magdir_theta(part([ii,jj,kk],linear_index))); 
 
  array.magdir_fn  =  @(ii,jj,kk)    ... 
    [ magdir_fn_comp{1}(ii,jj,kk)  ... 
      magdir_fn_comp{2}(ii,jj,kk)  ... 
      magdir_fn_comp{3}(ii,jj,kk) ]; 
 
end 
 
 
 
 
  
array.magloc  =  repmat(NaN,[array.total 3]); 
array.magdir  =  array.magloc; 
magloc_array  =  repmat(NaN,[array.mcount(1) array.mcount(2) array.mcount(3) 3]); 
 
ii  =  0; 
for xx  =  1:array.mcount(1) 
  for yy  =  1:array.mcount(2) 
    for zz  =  1:array.mcount(3) 
      ii  =  ii + 1; 
      array.magdir(ii,:)  =  array.magdir_fn(xx,yy,zz); 
    end 
  end 
end 
 
magsep_x  =  zeros(size(array.mcount(1))); 
magsep_y  =  zeros(size(array.mcount(2))); 
magsep_z  =  zeros(size(array.mcount(3))); 
 
magsep_x(1)  =  array.msize_array(1,1,1,1)/2; 
magsep_y(1)  =  array.msize_array(1,1,1,2)/2; 
magsep_z(1)  =  array.msize_array(1,1,1,3)/2; 
 
for ii  =  2:array.mcount(1) 
  magsep_x(ii)  =  array.msize_array(ii-1,1,1,1)/2  ... 
               + array.msize_array(ii  ,1,1,1)/2 ; 
end 
for jj  =  2:array.mcount(2) 
  magsep_y(jj)  =  array.msize_array(1,jj-1,1,2)/2  ... 
               + array.msize_array(1,jj  ,1,2)/2 ; 
end 
for kk  =  2:array.mcount(3) 
  magsep_z(kk)  =  array.msize_array(1,1,kk-1,3)/2  ... 
               + array.msize_array(1,1,kk  ,3)/2 ; 
end 
 
magloc_x  =  cumsum(magsep_x); 
magloc_y  =  cumsum(magsep_y); 
magloc_z  =  cumsum(magsep_z); 
 
for ii  =  1:array.mcount(1) 
  for jj  =  1:array.mcount(2) 
    for kk  =  1:array.mcount(3) 
      array.magloc_array(ii,jj,kk,:)  =   ... 
        [magloc_x(ii); magloc_y(jj); magloc_z(kk)]  ... 
        + [ii-1; jj-1; kk-1].*array.mgap; 
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
 
 
 
 
function debug_disp(str) 
  %disp(str) 
end 
 
 
 
 

