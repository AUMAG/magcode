%% An incorrect example of calculating the torques between magnets
%
% This is the accompanying code to a letter written Nov. 2009,
% commenting on a recent publication by Allag & Yonnet in the Oct.
% issue of the IEEE Transactions on Magnetics: "3-D Analytical
% Calculation of the Torque and Force Exerted Between Two Cuboidal Magnets"
%     <http://dx.doi.org/10.1109/TMAG.2009.2025047>
%
% Copyright 2009 Will Robertson
%
% Released under the Apache License v2.0:
%    <http://www.apache.org/licenses/LICENSE-2.0>
% This means, in essense, that you may freely modify and distribute this
% code provided that you acknowledge your changes to the work and retain
% my copyright. See the License text for the specific language governing
% permissions and limitations under the License.

%% Initialise

clear all
clc

%% Variables
% Cube magnets with z-displacement only and both with positive-z
% magnetisation.

a = 0.015/2;
b = 0.015/2;
c = 0.015/2;
A = 0.015/2;
B = 0.015/2;
C = 0.015/2;

J1 = 1;
J2 = 1;

offset = [0.04 0.04 0.04];

%% Now calculate the forces and torques.

% Pre-allocate variables:
torque_xyz = zeros(3,1);
torque_xyz2 = zeros(3,1);
f_comp_x = zeros(2,2,2);
f_comp_y = zeros(2,2,2);
f_comp_z = zeros(2,2,2);
l_comp_x = zeros(2,2,2);
l_comp_y = zeros(2,2,2);
l_comp_z = zeros(2,2,2);
indices = zeros(2,2,2,3);

% Avoid singularities:
sum_without_NaN = @(x) sum(x(~isnan(x)));
            
% Loop through the dimensions of the second magnet:
for jj = [0, 1]
  for ll = [0, 1]
    for qq = [0, 1]
      
      % These are the lever arms:
      lx = A*(-1)^jj;
      ly = B*(-1)^ll;
      lz = C*(-1)^qq;
      
      % Save them for later:
      l_comp_x(jj+1,ll+1,qq+1) = lx;
      l_comp_y(jj+1,ll+1,qq+1) = ly;
      l_comp_z(jj+1,ll+1,qq+1) = lz;
      
      indices(jj+1,ll+1,qq+1,:) = [jj, ll, qq];
      
      % Loop through the dimensions of the first magnet:
      for ii = [0, 1]
        for kk = [0, 1]
          for pp = [0, 1]

            % The distance between the corners:
            u = offset(1) + lx - a*(-1)^ii;
            v = offset(2) + ly - b*(-1)^kk;
            w = offset(3) + lz - c*(-1)^pp;
            r = sqrt(u^2+v^2+w^2);
            
            % These are the corner forces:
            f_x = sum_without_NaN([
              + 0.5*(v^2-w^2)*log(r-u);
              + u*v*log(r-v);
              + v*w*atan(u*v/(r*w));
              + 0.5*r*u
              ]);
            
            f_y = sum_without_NaN([
              + 0.5*(u^2-w^2)*log(r-v);
              + u*v*log(r-u);
              + u*w*atan(u*v/(r*w));
              + 0.5*r*v
              ]);
            
            f_z = sum_without_NaN([
              - u*w*log(r-u);
              - v*w*log(r-v);
              + u*v*atan(u*v/(r*w));
              - r*w
              ]);
            
            % These are the corner torques:
            t_x = f_y*(lz-w/2) - f_z*(ly-v/2);
            t_y = f_z*(lx-u/2) - f_x*(lz-w/2);
            t_z = f_x*(ly-v/2) - f_y*(lx-u/2);
            

            index_sum = (-1)^(ii+jj+kk+ll+pp+qq);

            % This is the total torque:
            torque_xyz = torque_xyz + cross( [lx; ly; lz] , [f_x; f_y; f_z] )*index_sum;
            torque_xyz2 = torque_xyz2 + [t_x;t_y;t_z]*index_sum
            
            % These are the total corner forces from the first magnet onto
            % a single corner of the second:
            f_comp_x(jj+1,ll+1,qq+1) = f_comp_x(jj+1,ll+1,qq+1) + f_x*index_sum;
            f_comp_y(jj+1,ll+1,qq+1) = f_comp_y(jj+1,ll+1,qq+1) + f_y*index_sum;
            f_comp_z(jj+1,ll+1,qq+1) = f_comp_z(jj+1,ll+1,qq+1) + f_z*index_sum;

          end
        end
      end
    end
  end
end

% The total torque with correct units:
magconst = 1/(4*pi*(4*pi*1e-7));
torque_xyz = J1*J2*magconst*torque_xyz;
torque_xyz2 = J1*J2*magconst*torque_xyz2;

% Calculate and reshape the results into single column vectors:
f_corner_x = J1*J2*magconst*reshape(f_comp_x,[8 1]);
f_corner_y = J1*J2*magconst*reshape(f_comp_y,[8 1]);
f_corner_z = J1*J2*magconst*reshape(f_comp_z,[8 1]);
l_corner_x = reshape(l_comp_x,[8 1]);
l_corner_y = reshape(l_comp_y,[8 1]);
l_corner_z = reshape(l_comp_z,[8 1]);
corners = reshape(indices,[8 3]);

% Calculate total force:
fx = sum(f_corner_x);
fy = sum(f_corner_y);
fz = sum(f_corner_z);

% Components of and total x-torque:
tx_corner_y =  l_corner_y.*f_corner_z;
tx_corner_z = -l_corner_z.*f_corner_y;
tx_corner   = tx_corner_y+tx_corner_z;
tx = sum(tx_corner);

% Components of and total y-torque:
ty_corner_x =  l_corner_z.*f_corner_x;
ty_corner_z = -l_corner_x.*f_corner_z;
ty_corner   = ty_corner_x+ty_corner_z;
ty = sum(ty_corner);

% Components of and total z-torque:
tz_corner_x =  l_corner_x.*f_corner_y;
tz_corner_y = -l_corner_y.*f_corner_x;
tz_corner   = tz_corner_x+tz_corner_y;
tz = sum(tz_corner);

%% Display the results

format short g

fprintf('\nTotal forces:\n')
disp([fx fy fz])
disp('These forces are correct for magnets with z-displacement only.')

fprintf('\nTotal torques calculated with three different methods:\n')
disp([tx ty tz])
disp(torque_xyz')
disp(torque_xyz2')
disp('These torques must be wrong; they should all be zero for z-displacement only.')

fprintf('\nBreakdown of y-torque into components of the x-force and z-force for each corner:\n\n')
disp('   y-torque-x   y-torque-z')
disp([ty_corner_x ty_corner_z])
disp('SUM  --------     --------')
disp([sum(ty_corner_x) sum(ty_corner_z)])
fprintf('\nTherefore, only the z-force needs examining.\n')

fprintf('\nBreakdown of the y-torque into components per corner of z-force:\n\n')
disp('                                 corner      z-force   y-torque-z')
disp([corners round(f_corner_z,3) round(ty_corner_z,3)])

fprintf(['\nThese forces are not symmetrical and therefore invalidate the argument\n',...
             'that the corner forces can be used for calculating torques.\n'])

           
disp(torque_xyz2)
magnetforces(struct('dim',[2*a 2*b 2*c],'magdir',[0 0 1],'magn',1),struct('dim',[2*A 2*B 2*C],'magdir',[0 0 1],'magn',1),offset,'torque')

