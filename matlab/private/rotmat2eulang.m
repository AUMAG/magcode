%\begin{mfunction}{rotmat2eulang}
function eul = rotmat2eulang(R,order)
% ROTMAT2EULANG Calculates fixed-euler angles from rotation matrix
%
% * `R`: Rotation matrix (3x3)
% * `order`: Order of Euler angles to return; any of:
%              {'XYZ','XZY','YXZ','YZX','ZYX','ZXY'}
% * `EUL`: Euler angles in specified order

switch order

  case 'XYZ'
    eul(1) = atan2(-R(2,3),R(3,3));
    eul(2) = asin(R(1,3));
    eul(3) = atan2(-R(1,2),R(1,1));
    
  case 'XZY'
    eul(1) = atan2(R(3,2),R(2,2));
    eul(2) = -asin(R(1,2));
    eul(3) = atan2(R(1,3),R(1,1));
    
  case 'YXZ'
    eul(1) = atan2(R(1,3),R(3,3));
    eul(2) = -asin(R(2,3));
    eul(3) = atan2(R(2,1),R(2,2));
    
  case 'YZX'
    eul(1) = atan2(-R(3,1),R(1,1));
    eul(2) = asin(R(2,1));
    eul(3) = atan2(-R(2,3),R(2,2));
    
  case 'ZXY'
    eul(1) = atan2(-R(1,2),R(2,2));
    eul(2) = asin(R(3,2));
    eul(3) = atan2(-R(3,1),R(3,3));

  case 'ZYX'
    eul(1) = atan2(R(2,1),R(1,1));
    eul(2) = -asin(R(3,1));
    eul(3) = atan2(R(3,2),R(3,3));
    
  otherwise
    error('Don''t know specified rotation order %s', order)

end

end
%\end{mfunction}