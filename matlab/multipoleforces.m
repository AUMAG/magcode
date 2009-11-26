  
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
 
    displ  =  displ - fixed_array.loc(mm,:) + float_array.loc(nn,:); 
 
    if calc_force_bool 
      array_forces(mm,nn,:)  =   ... 
          magnetforces(fixed_magnet, float_magnet, displ,'force'); 
    end 
 
    if calc_stiffness_bool 
      array_stiffnesses(mm,nn,:)  =   ... 
          magnetforces(fixed_magnet, float_magnet, displ,'stiffness'); 
    end 
 
  end 
end 
 
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
 
 
 
 
  
  
function array_out  =  complete_array_from_input(array) 
 
if ~isfield(array,'type') 
  array.type  =  'generic'; 
end 
 
array.total  =  prod(array.mcount); 
array.loc  =  repmat(NaN,[array.total 3]); 
array.magdir  =  repmat(NaN,[array.total 2]); 
 
  
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
 
 
 
  
if length(array.mgap) == 3 
  array.gaps  =  array.mgap; 
elseif length(array.mgap) == 1 
  array.gaps  =  repmat(array.mgap, [3 1]); 
else 
  error('Not yet implemented.') 
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
 
 
 
  
ii  =  0; 
for xx  =  1:array.mcount(1) 
  for yy  =  1:array.mcount(2) 
    for zz  =  1:array.mcount(3) 
      ii  =  ii + 1; 
      array.loc(ii,:)  =   ... 
        [xx-1; yy-1; zz-1].*(squeeze(array.dim_array(xx,yy,zz,:))+array.gaps); 
    end 
  end 
end 
 
debug_disp('Magnet locations:') 
debug_disp(array.loc) 
 
 
 
 
array_out  =  array; 
 
end 
 
 
 
function debug_disp(str) 
  %disp(str) 
end 
 
 
 
end 
 

