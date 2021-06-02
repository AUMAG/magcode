%% cuboid_stiffness
%
% Calculate the stiffnesses between two parallel cuboid magnets with
% magnetisations along one of the primary axes.

% \START

% \begin{mfunction}{cuboid_stiffness}
%
% Cuboid forces are written out explicitly for speed. Try something
% differently here to shorten the code.

function stiffness_out = cuboid_stiffness(size1,size2,offset,J1,J2)

stiffness_out = zeros(size(offset));

stiffness_out = stiffness_out + stiffnesses_calc_z_z(size1,size2,offset,J1(1),J2(1),[-3;  2; 1], [ 3;  2;  1]); % xx
stiffness_out = stiffness_out + stiffnesses_calc_z_y(size1,size2,offset,J1(1),J2(2),[-3;  2; 1], [ 3;  2;  1]); % xy
stiffness_out = stiffness_out + stiffnesses_calc_z_y(size1,size2,offset,J1(1),J2(3),[-2; -3; 1], [-3; -1; -2]); % xz
stiffness_out = stiffness_out + stiffnesses_calc_z_y(size1,size2,offset,J1(2),J2(1),[ 3;  1; 2], [ 2;  3;  1]); % yx
stiffness_out = stiffness_out + stiffnesses_calc_z_z(size1,size2,offset,J1(2),J2(2),[ 1; -3; 2], [ 1;  3;  2]); % yy
stiffness_out = stiffness_out + stiffnesses_calc_z_y(size1,size2,offset,J1(2),J2(3),[ 1; -3; 2], [-1; -3; -2]); % yz
stiffness_out = stiffness_out + stiffnesses_calc_z_y(size1,size2,offset,J1(3),J2(1),[-2;  1; 3], [ 2;  1;  3]); % zx
stiffness_out = stiffness_out + stiffnesses_calc_z_y(size1,size2,offset,J1(3),J2(2),[ 1;  2; 3], [ 1;  2;  3]); % zy
stiffness_out = stiffness_out + stiffnesses_calc_z_z(size1,size2,offset,J1(3),J2(3),[ 1;  2; 3], [ 1;  2;  3]); % zz

end

% \end{mfunction}

% \begin{mfunction}{stiffnesses_calc_z_z}

function calc_out = stiffnesses_calc_z_z(size1,size2,offset,J1,J2,di,oi)

if (abs(J1)<eps || abs(J2)<eps)
  calc_out = zeros(size(offset));
  return
end

si = abs(di);

comps_xyz = zeros(size(offset));

for ii = [0,1]
  for jj = [0,1]
    for kk = [0,1]
      for ll = [0,1]
        for pp = [0,1]
          for qq = [0,1]
            
            u = sign(di(1))*offset(si(1),:) + size2(si(1))*(-1).^jj - size1(si(1))*(-1).^ii;
            v = sign(di(2))*offset(si(2),:) + size2(si(2))*(-1).^ll - size1(si(2))*(-1).^kk;
            w = sign(di(3))*offset(si(3),:) + size2(si(3))*(-1).^qq - size1(si(3))*(-1).^pp;
            r = sqrt(u.^2+v.^2+w.^2);
            
            component_x = - r - (u.^2 .*v)./(u.^2+w.^2) - v.*log(r-v) ;
            component_y = - r - (v.^2 .*u)./(v.^2+w.^2) - u.*log(r-u) ;
            component_z = - component_x - component_y;
            
            comps_xyz = comps_xyz + (-1)^(ii+jj+kk+ll+pp+qq).*[component_x; component_y; component_z];
            
          end
        end
      end
    end
  end
end

calc_out = J1*J2/(4*pi*(4*pi*1e-7))*sign(oi).*comps_xyz(si,:);

end

% \end{mfunction}


% \begin{mfunction}{stiffnesses_calc_z_y}

function calc_out = stiffnesses_calc_z_y(size1,size2,offset,J1,J2,di,oi)

if (abs(J1)<eps || abs(J2)<eps)
  calc_out = zeros(size(offset));
  return
end

si = abs(di);
comps_xyz = zeros(size(offset));

for ii = [0,1]
  for jj = [0,1]
    for kk = [0,1]
      for ll = [0,1]
        for pp = [0,1]
          for qq = [0,1]
            
            u = sign(di(1))*offset(si(1),:) + size2(si(1))*(-1).^jj - size1(si(1))*(-1).^ii;
            v = sign(di(2))*offset(si(2),:) + size2(si(2))*(-1).^ll - size1(si(2))*(-1).^kk;
            w = sign(di(3))*offset(si(3),:) + size2(si(3))*(-1).^qq - size1(si(3))*(-1).^pp;
            r = sqrt(u.^2+v.^2+w.^2);
            
            component_x =  ((u.^2 .*v)./(u.^2 + v.^2)) + (u.^2 .*w)./(u.^2 + w.^2) ...
              - u.*atan1(v.*w,r.*u) + multiply_x_log_y( w , r + v ) + ...
              + multiply_x_log_y( v , r + w );

            component_y = - v/2 + (u.^2 .*v)./(u.^2 + v.^2) - (u.*v.*w)./(v.^2 + w.^2) ...
              -  u.*atan1(u.*w,r.*v) - multiply_x_log_y( v , r + w );

            component_z = - component_x - component_y;
            
            comps_xyz = comps_xyz + (-1)^(ii+jj+kk+ll+pp+qq).*[component_x; component_y; component_z];
            
          end
        end
      end
    end
  end
end

calc_out = J1*J2/(4*pi*(4*pi*1e-7))*sign(oi).*comps_xyz(abs(oi),:);

end

% \end{mfunction}

% \subsubsection{Helpers}
%
%  The equations contain two singularities. Specifically, the equations
%  contain terms of the form $x \log(y)$, which becomes |NaN| when both $x$
%  and $y$ are zero since $\log(0)$ is negative infinity.
%
% \begin{mfunction}{multiply_x_log_y}
%  This function computes $x \log(y)$, special-casing the singularity to output
%  zero, instead. (This is indeed the value of the limit.)
function out = multiply_x_log_y(x,y)
out = x.*log(y);
out(~isfinite(out))=0;
end
% \end{mfunction}

% \begin{mfunction}{atan1}
% We're using |atan| instead of |atan2| (otherwise the wrong results
%	are calculated --- I guess I don't totally understand that), which becomes
%	a problem when trying to compute |atan(0/0)| since |0/0| is |NaN|.
function out = atan1(x,y)
out = zeros(size(x));
ind = x~=0 & y~=0;
out(ind) = atan(x(ind)./y(ind));
end
% \end{mfunction}
