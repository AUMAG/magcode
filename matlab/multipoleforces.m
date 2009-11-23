  
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
 
 
 
  
M  =  prod(fixed_array.mcount); 
N  =  prod(float_array.mcount); 
 
fixed_magnet_loc  =  repmat(NaN,[M 3]); 
float_magnet_loc  =  repmat(NaN,[N 3]); 
 
fixed_magnet_magdir  =  repmat(NaN,[M 2]); 
float_magnet_magdir  =  repmat(NaN,[N 2]); 
 
  
if length(fixed_array.msize) == 3 
  fixed_magnet_dim_array  =   ... 
      repmat(reshape(fixed_array.msize,[1 1 1 3]), fixed_array.mcount); 
  fixed_magnet_dim  =  reshape(fixed_magnet_dim_array, [M 3]); 
else 
  error('Not yet implemented.') 
end 
 
if length(float_array.msize) == 3 
  float_magnet_dim_array  =   ... 
      repmat(reshape(float_array.msize,[1 1 1 3]), float_array.mcount); 
  float_magnet_dim  =  reshape(float_magnet_dim_array, [N 3]); 
else 
  error('Not yet implemented.') 
end 
 
  
if length(fixed_array.magn) == 1 
  fixed_magnet_magn  =  repmat(fixed_array.magn,[M 1]); 
else 
  error('Not yet implemented.') 
end 
 
if length(float_array.magn) == 1 
  float_magnet_magn  =  repmat(float_array.magn,[N 1]); 
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
      fixed_magnet_loc(ii,:)  =   ... 
        [xx-1; yy-1; zz-1].*(squeeze(fixed_magnet_dim_array(xx,yy,zz,:))+fixed_gaps); 
    end 
  end 
end 
 
ii  =  0; 
for xx  =  1:float_array.mcount(1) 
  for yy  =  1:float_array.mcount(2) 
    for zz  =  1:float_array.mcount(3) 
      ii  =  ii + 1; 
      float_magnet_loc(ii,:)  =   ... 
        [xx-1; yy-1; zz-1].*(squeeze(float_magnet_dim_array(xx,yy,zz,:))+float_gaps); 
    end 
  end 
end 
 
  
ii  =  0; 
for xx  =  1:fixed_array.mcount(1) 
  for yy  =  1:fixed_array.mcount(2) 
    for zz  =  1:fixed_array.mcount(3) 
      ii  =  ii + 1; 
      fixed_magnet_magdir(ii,:)  =  fixed_array.magdir_fn(xx,yy,zz); 
    end 
  end 
end 
 
 
ii  =  0; 
for xx  =  1:float_array.mcount(1) 
  for yy  =  1:float_array.mcount(2) 
    for zz  =  1:float_array.mcount(3) 
      ii  =  ii + 1; 
      float_magnet_magdir(ii,:)  =  float_array.magdir_fn(xx,yy,zz); 
    end 
  end 
end 
 
 
 
 
  
if calc_force_bool 
  array_forces  =  repmat(NaN,[M N 3]); 
end 
 
if calc_stiffness_bool 
  array_stiffnesses  =  repmat(NaN,[M N 3]); 
end 
 
for mm  =  1:M 
 
  fixed_magnet  =  struct( ... 
 'dim',    fixed_magnet_dim(mm,:),  ... 
 'magn',   fixed_magnet_magn(mm),  ... 
 'magdir', fixed_magnet_magdir(mm,:)  ... 
  ); 
 
  for nn  =  1:N 
 
    float_magnet  =  struct( ... 
      'dim',    float_magnet_dim(nn,:),  ... 
      'magn',   float_magnet_magn(nn),  ... 
      'magdir', float_magnet_magdir(nn,:)  ... 
    ); 
 
    displ  =  displ - fixed_magnet_loc(mm,:) + float_magnet_loc(nn,:); 
 
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
 
 
 
 

