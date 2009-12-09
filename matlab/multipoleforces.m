  
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
 
 
for mm  =  1:fixed_array.total 
 
  fixed_magnet  =  struct( ... 
 'dim',    fixed_array.dim(mm,:),  ... 
 'magn',   fixed_array.magn(mm),  ... 
 'magdir', fixed_array.magdir(mm,:)  ... 
  ); 
 
  for nn  =  1:float_array.total 
 
    float_magnet  =  struct( ... 
      'dim',    float_array.dim(nn,:),  ... 
      'magn',   float_array.magn(nn),  ... 
      'magdir', float_array.magdir(nn,:)  ... 
    ); 
 
    mag_displ  =  displ - fixed_array.loc(mm,:) + float_array.loc(nn,:); 
 
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
 
  case 'linear-x' 
    linear_index  =  1; 
  otherwise 
    error(['Unknown array type ''',array.type,'''.']) 
end 
 
array.mcount_extra  =  0; 
if isfield(array,'Nwaves') 
  array.mcount_extra  =  1; 
end 
 
if strncmp(array.type,'linear',6) 
 
  var_names  =  {'wavelength','length','Nwaves','mlength', ... 
               'Nmag','Nmag_per_wave','magdir_rotate'}; 
 
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
 
  array.Nmag  =  round(array.Nmag) + array.mcount_extra; 
  array.Nmag_per_wave  =  round(array.Nmag_per_wave); 
 
  array.mlength  =  array.mlength * (array.Nmag-array.mcount_extra)/array.Nmag; 
 
  array.mcount  =  ones(1,3); 
  array.mcount(linear_index)  =  array.Nmag; 
 
end 
 
array.total  =  prod(array.mcount); 
 
  
if ~isfield(array,'msize') 
  array.msize  =  [array.mlength array.width array.height]; 
end 
 
if length(array.msize) == 3 
  array.dim_array  =   ... 
      repmat(reshape(array.msize,[1 1 1 3]), array.mcount); 
  array.dim  =  reshape(array.dim_array, [array.total 3]); 
else 
  error('Not yet implemented.') 
end 
 
 
 
  
if length(array.magn) == 1 
  array.magn  =  repmat(array.magn,[array.total 1]); 
else 
  error('Not yet implemented.') 
end 
 
 
 
  
if ~isfield(array,'mgap') 
  array.mgap  =  [0; 0; 0]; 
end 
 
if length(array.mgap) == 1 
  array.mgap  =  repmat(array.mgap, [3 1]); 
end 
 
 
 
  
array.magdir  =  repmat(NaN,[array.total 2]); 
 
switch array.face 
  case 'up' 
    array.magdir_rotate_sign  =   1; 
  case 'down' 
    array.magdir_rotate_sign  =  -1; 
  otherwise 
    if ~isfield(array,'magdir_rotate_sign') 
      disp('huh?') 
      array.magdir_rotate_sign  =  1; 
    else 
      error('huh?') 
    end 
end 
 
if ~isfield(array,'magdir_fn') 
  switch array.type 
    case 'linear-x' 
      array.magdir_fn  =  @(ii,jj,kk)   ... 
        [0 array.magdir_first+array.magdir_rotate_sign * array.magdir_rotate * (ii-1)]; 
  end 
end 
 
ii  =  0; 
for xx  =  1:array.mcount(1) 
  for yy  =  1:array.mcount(2) 
    for zz  =  1:array.mcount(3) 
      ii  =  ii + 1; 
      array.magdir(ii,:)  =  array.magdir_fn(xx,yy,zz); 
    end 
  end 
end 
 
debug_disp('Magnetisation directions') 
debug_disp(mod(array.magdir,360)) 
 
 
 
  
array.loc  =  repmat(NaN,[array.total 3]); 
 
ii  =  0; 
for xx  =  1:array.mcount(1) 
  for yy  =  1:array.mcount(2) 
    for zz  =  1:array.mcount(3) 
      ii  =  ii + 1; 
      array.loc(ii,:)  =   ... 
        [xx-1; yy-1; zz-1].* ... 
        (squeeze(array.dim_array(xx,yy,zz,:)) + array.mgap); 
    end 
  end 
end 
 
debug_disp('Magnet locations:') 
debug_disp(array.loc) 
 
 
 
 
array_out  =  array; 
 
end 
 
 
 
function debug_disp(str) 
  disp(str) 
end 
 
 
 
 

