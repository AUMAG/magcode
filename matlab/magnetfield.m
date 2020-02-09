function [magB] = magnetfield(mag,points,varargin)
%MAGNETFIELD Calculate magnetic field from a magnet source

switch mag.type
  
  case 'cuboid'
    
    magB = cuboid_field(mag,points);
    
  case 'cylinder'
    
    if isequal(mag.magdir,[0;0;1]) && isequal(mag.dir,[0;0;1])
      magB = cylinder_field_axial(mag,points);
    end
    
end

end





