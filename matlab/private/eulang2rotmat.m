function rotM = eulang2rotmat(eul,order)
%
% Calculate rotation matrix `rotM` from Euler angles `eul` (3x1 or 1x3).
% Angles are specified in radians.
%
% `order` specifies the axes order to apply the Euler angles, any of:
%     {'XYZ','XZY','YXZ','YZX','ZYX','ZXY'}
%
% E.g., for planar rotation matrices {Rx, Ry, Rz}, eulang2rotmat([a b c],'ZYX')
% creates rotation matrix Rz(a)*Ry(b)*Rx(c).

oo = order - 'X' + 1; % 'X' -> 1, 'Y' -> 2, 'Z' -> 3
reord = [find(order=='X'),find(order=='Y'),find(order=='Z')];

c = cos(eul(reord));
s = sin(eul(reord));

RxRyRz(:,:,1) = [1 0 0; 0 c(1) -s(1); 0 s(1) c(1)];
RxRyRz(:,:,2) = [c(2) 0 s(2); 0 1 0; -s(2) 0 c(2)];
RxRyRz(:,:,3) = [c(3) -s(3) 0; s(3) c(3) 0; 0 0 1];

rotM = RxRyRz(:,:,oo(1))*RxRyRz(:,:,oo(2))*RxRyRz(:,:,oo(3));

end

% Licence included in README.