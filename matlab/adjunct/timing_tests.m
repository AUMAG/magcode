function timing_tests()

magconst = 1/(4*pi*(4*pi*1e-7));

[index_i, index_j, index_k, index_l, index_p, index_q] = ndgrid([0 1]);

index_sum = (-1).^(index_i+index_j+index_k+index_l+index_p+index_q);

calc_xyz = [true true true];

tol = 1e-9;

checkloop = 6;

s1 = rand(1,3); s2 = rand(1,3); dd = rand(1,3); J1 = rand(1); J2 = rand(1);
F1 = forces_calc_sumsum(s1,s2,dd,J1,J2);
F2 = forces_calc_loop(  s1,s2,dd,J1,J2);
F3 = forces_calc_loop2( s1,s2,dd,J1,J2);
[F1 F2 F3]
assert(all(abs(F2-F1)<tol));
assert(all(abs(F3-F1)<tol));

s1 = [1 1 1]; s2 = [1 1 1]; dd = 2*[1 1 1]; J1 = 1; J2 = 1;
F1 = forces_calc_sumsum(s1,s2,dd,J1,J2);
F2 = forces_calc_loop(  s1,s2,dd,J1,J2);
F3 = forces_calc_loop2( s1,s2,dd,J1,J2);
[F1 F2 F3]
assert(all(abs(F2-F1)<tol));
assert(all(abs(F3-F1)<tol));

%%

N = 1000;
tic
for noop = 1:N
  F = forces_calc_sumsum(rand(1,3),rand(1,3),rand(1,3),rand(1),rand(1));
end
toc
tic
for noop = 1:N
  F = forces_calc_loop(rand(1,3),rand(1,3),rand(1,3),rand(1),rand(1));
end
toc
tic
for noop = 1:N
  F = forces_calc_loop2(rand(1,3),rand(1,3),rand(1,3),rand(1),rand(1));
end
toc

%%

  function calc_out = forces_calc_loop(size1,size2,offset,J1,J2)
    
    component_x = 0;
    component_y = 0;
    component_z = 0;
    
    count = 0;
    for ii = [1 -1]
      for jj = [1 -1]
        for kk = [1 -1]
          for ll = [1 -1]
            for pp = [1 -1]
              for qq = [1 -1]
                count = count + 1;
                
                u = offset(1) + size2(1)*jj - size1(1)*ii;
                v = offset(2) + size2(2)*ll - size1(2)*kk;
                w = offset(3) + size2(3)*qq - size1(3)*pp;
                r = sqrt(u.^2+v.^2+w.^2);
                
                cx = ...
                  + multiply_x_log_y( 0.5*(v.^2-w.^2), r-u ) ...
                  + multiply_x_log_y( u.*v, r-v ) ...
                  + v.*w.*atan1(u.*v,r.*w) ...
                  + 0.5*r.*u;
                
                cy = ...
                  + multiply_x_log_y( 0.5*(u.^2-w.^2), r-v ) ...
                  + multiply_x_log_y( u.*v, r-u ) ...
                  + u.*w.*atan1(u.*v,r.*w) ...
                  + 0.5*r.*v;
                
                cz = ...
                  - multiply_x_log_y( u.*w, r-u ) ...
                  - multiply_x_log_y( v.*w, r-v ) ...
                  + u.*v.*atan1(u.*v,r.*w) ...
                  - r.*w;
                
                ind_sum = ii*jj*kk*ll*pp*qq;
                component_x = component_x + ind_sum.*cx;
                component_y = component_y + ind_sum.*cy;
                component_z = component_z + ind_sum.*cz;
                
              end
            end
          end
        end
      end
    end
    
    calc_out = J1*J2*magconst.*[component_x; component_y; component_z];
    
    function out = multiply_x_log_y(x,y)
      out = x.*log(y);
      out(~isfinite(out))=0;
    end
    
    %%
    
    function out = atan1(x,y)
      out = zeros(size(x));
      ind = x~=0 & y~=0;
      out(ind) = atan(x(ind)./y(ind));
    end
  end


  function calc_out = forces_calc_loop2(size1,size2,offset,J1,J2)
    
    component_x = 0;
    component_y = 0;
    component_z = 0;
    
    count = 0;
    for ii = [1 -1]
      for jj = [1 -1]
        for kk = [1 -1]
          for ll = [1 -1]
            for pp = [1 -1]
              for qq = [1 -1]
                
                count = count + 1;
                
                u = offset(1) + size2(1)*jj - size1(1)*ii;
                v = offset(2) + size2(2)*ll - size1(2)*kk;
                w = offset(3) + size2(3)*qq - size1(3)*pp;
                r = sqrt(u.^2+v.^2+w.^2);
                
                if w == 0
                  atan_term = 0;
                else
                  atan_term = atan(u.*v./(r.*w));
                end
                if abs(r-u) < eps
                  log_ru = 0;
                else
                  log_ru = log(r-u);
                end
                if abs(r-v) < eps
                  log_rv = 0;
                else
                  log_rv = log(r-v);
                end
                
                cx = ...
                  + 0.5*(v.^2-w.^2).*log_ru ...
                  + u.*v.*log_rv ...
                  + v.*w.*atan_term...
                  + 0.5*r.*u;
                
                cy = ...
                  + 0.5*(u.^2-w.^2).*log_rv ...
                  + u.*v.*log_ru ...
                  + u.*w.*atan_term ...
                  + 0.5*r.*v;
                
                cz = ...
                  - u.*w.*log_ru ...
                  - v.*w.*log_rv ...
                  + u.*v.*atan_term ...
                  - r.*w;

                ind_sum = ii*jj*kk*ll*pp*qq;
                component_x = component_x + ind_sum.*cx;
                component_y = component_y + ind_sum.*cy;
                component_z = component_z + ind_sum.*cz;
                
              end
            end
          end
        end
      end
    end
    
    calc_out = J1*J2*magconst.*[component_x; component_y; component_z];
    
  end

%%

  function calc_out = forces_calc_sumsum(size1,size2,offset,J1,J2)
    
    u = offset(1) + size2(1)*(-1).^index_j - size1(1)*(-1).^index_i;
    v = offset(2) + size2(2)*(-1).^index_l - size1(2)*(-1).^index_k;
    w = offset(3) + size2(3)*(-1).^index_q - size1(3)*(-1).^index_p;
    r = sqrt(u.^2+v.^2+w.^2);
    
    if calc_xyz(1)
      component_x = ...
        + multiply_x_log_y( 0.5*(v.^2-w.^2), r-u ) ...
        + multiply_x_log_y( u.*v, r-v ) ...
        + v.*w.*atan1(u.*v,r.*w) ...
        + 0.5*r.*u;
    end
    
    if calc_xyz(2)
      component_y = ...
        + multiply_x_log_y( 0.5*(u.^2-w.^2), r-v ) ...
        + multiply_x_log_y( u.*v, r-u ) ...
        + u.*w.*atan1(u.*v,r.*w) ...
        + 0.5*r.*v;
    end
    
    if calc_xyz(3)
      component_z = ...
        - multiply_x_log_y( u.*w, r-u ) ...
        - multiply_x_log_y( v.*w, r-v ) ...
        + u.*v.*atan1(u.*v,r.*w) ...
        - r.*w;
    end
    
    component_x = index_sum.*component_x;
    component_y = index_sum.*component_y;
    component_z = index_sum.*component_z;
    
    calc_out = J1*J2*magconst .* ...
      [ sum(component_x(:)) ; sum(component_y(:)) ; sum(component_z(:)) ] ;
    
    function out = multiply_x_log_y(x,y)
      out = x.*log(y);
      out(~isfinite(out))=0;
    end
    
    function out = atan1(x,y)
      out = zeros(size(x));
      ind = x~=0 & y~=0;
      out(ind) = atan(x(ind)./y(ind));
    end
    
  end

end