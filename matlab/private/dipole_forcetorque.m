%% dipole_forcetorque
%
% Calculate the forces and torques between two magnetic dipoles.

% \START

% \begin{mfunction}{dipole_forcetorque}

function [f,t] = dipole_forcetorque(m_a,m_b,r_ab)
% [F,T] = CALCDIPOLEFORCETORQUE(MA,MB,R)
%
% Calculates the force and torque on magnetic dipole MB due to magnetic
% dipole MA with distance vector RAB.
%
% MA and MB must have equal size and be 3x1 or 3xM vector arrays.
% RAB 
%
% Magnetic dipole moments can be calculated for a hard magnet using the equation
%           m = 1/(4*pi*1e-7)*Br*V
% where Br is the remanence magnetisation vector in Tesla and V is the magnet volume.
%
% For a coil, the magnetic dipole moment is
%           m = N*I*A*n
% where I is the current, N the number of turns, A the cross-sectional area
% of the current loop(s), and n is the normal vector to the loop(s).
%
% These equations have been adapted from the following pair of papers:
% * Yung 1998:      http://doi.org/10.1155/1998/79537
% * Landecker 1999: http://doi.org/10.1155/1999/97902


assert( size(r_ab,1)==3 , "Displacement vector RAB must be 3xM size")
assert( size(m_a,1)==3 && size(m_b,1)==3 , "Dipole moment vectors MA and MB must be 3xM size")
assert( size(m_a,2)==size(m_b,2) , "Dipole moment vectors MA and MB must be equal sizes")
assert( size(m_a,2)==size(r_ab,2) , "Dipole moment vectors MA and MB and displacement vector RAB must be equal sizes")

% replicate dipole moment vectors to the same length as the displacement vector:
N = size(r_ab,2);
if size(m_a,2) == 1
  m_a = repmat(m_a,[1,N]);
  m_b = repmat(m_b,[1,N]);
end

R_ab = sqrt(r_ab(1,:).^2+r_ab(2,:).^2+r_ab(3,:).^2);
rnorm   = r_ab./R_ab;

dot_rma = dot(rnorm,m_a);
dot_rmb = dot(rnorm,m_b);

f = 3e-7./R_ab.^4.*(...
   + rnorm.*dot(m_a,m_b) ...
   + m_a.*dot_rmb ...
   + m_b.*dot_rma ...
   - 5*rnorm.*dot_rma.*dot_rmb ...
  );

t = 1e-7./R_ab.^3.*( ...
  + cross(3*m_b,dot_rma.*rnorm) ...  
  - cross(m_a,m_b) ...
  );

end

% \end{mfunction}

