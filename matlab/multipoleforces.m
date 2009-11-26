  
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
 
 
 
  
if ~isfield(fixed_array,'type') 
  fixed_array.type  =  'generic'; 
end 
 
switch fixed_array.type 
 
end 
 
fixed_array.total  =  prod(fixed_array.mcount); 
fixed_array.loc  =  repmat(NaN,[fixed_array.total 3]); 
fixed_array.magdir  =  repmat(NaN,[fixed_array.total 2]); 
 
float_array.total  =  prod(float_array.mcount); 
float_array.loc  =  repmat(NaN,[float_array.total 3]); 
float_array.magdir  =  repmat(NaN,[float_array.total 2]); 
 
 
  
if length(fixed_array.msize) == 3 
  fixed_array.dim_array  =   ... 
      repmat(reshape(fixed_array.msize,[1 1 1 3]), fixed_array.mcount); 
  fixed_array.dim  =  reshape(fixed_array.dim_array, [fixed_array.total 3]); 
else 
  error('Not yet implemented.') 
end 
 
if length(float_array.msize) == 3 
  float_array.dim_array  =   ... 
      repmat(reshape(float_array.msize,[1 1 1 3]), float_array.mcount); 
  float_array.dim  =  reshape(float_array.dim_array, [float_array.total 3]); 
else 
  error('Not yet implemented.') 
end 
 
 
  
if length(fixed_array.magn) == 1 
  fixed_array.magn  =  repmat(fixed_array.magn,[fixed_array.total 1]); 
else 
  error('Not yet implemented.') 
end 
 
if length(float_array.magn) == 1 
  float_array.magn  =  repmat(float_array.magn,[float_array.total 1]); 
else 
  error('Not yet implemented.') 
end 
 
  
if length(fixed_array.mgap) == 3 
  fixed_gaps  =  fixed_array.mgap; 
elseif length(fixed_array.mgap) == 1 
  fixed_gaps  =  repmat(fixed_array.mgap, [3 1]); 
else 
  error('Not yet implemented.') 
end 
 
if length(float_array.mgap) == 3 
  float_gaps  =  float_array.mgap; 
elseif length(float_array.mgap) == 1 
  float_gaps  =  repmat(float_array.mgap, [3 1]); 
else 
  error('Not yet implemented.') 
end 
 
 
  
ii  =  0; 
for xx  =  1:fixed_array.mcount(1) 
  for yy  =  1:fixed_array.mcount(2) 
    for zz  =  1:fixed_array.mcount(3) 
      ii  =  ii + 1; 
      fixed_array.magdir(ii,:)  =  fixed_array.magdir_fn(xx,yy,zz); 
    end 
  end 
end 
 
debug_disp('Fixed magnetisation directions') 
debug_disp(mod(fixed_array.magdir,360)) 
 
ii  =  0; 
for xx  =  1:float_array.mcount(1) 
  for yy  =  1:float_array.mcount(2) 
    for zz  =  1:float_array.mcount(3) 
      ii  =  ii + 1; 
      float_array.magdir(ii,:)  =  float_array.magdir_fn(xx,yy,zz); 
    end 
  end 
end 
 
debug_disp('Float magnetisation directions') 
debug_disp(mod(float_array.magdir,360)) 
 
 
  
ii  =  0; 
for xx  =  1:fixed_array.mcount(1) 
  for yy  =  1:fixed_array.mcount(2) 
    for zz  =  1:fixed_array.mcount(3) 
      ii  =  ii + 1; 
      fixed_array.loc(ii,:)  =   ... 
        [xx-1; yy-1; zz-1].*(squeeze(fixed_array.dim_array(xx,yy,zz,:))+fixed_gaps); 
    end 
  end 
end 
 
debug_disp('Fixed magnet locations:') 
debug_disp(fixed_array.loc) 
 
ii  =  0; 
for xx  =  1:float_array.mcount(1) 
  for yy  =  1:float_array.mcount(2) 
    for zz  =  1:float_array.mcount(3) 
      ii  =  ii + 1; 
      float_array.loc(ii,:)  =   ... 
        [xx-1; yy-1; zz-1].*(squeeze(float_array.dim_array(xx,yy,zz,:))+float_gaps); 
    end 
  end 
end 
 
debug_disp('Float magnet locations:') 
debug_disp(float_array.loc) 
 
 
 
  
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
 
 
 
 
function debug_disp(str) 
  %disp(str) 
end 
 
end 
 

