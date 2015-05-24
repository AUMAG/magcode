function timing_tests()

debug_disp = @(s) disp([])

magconst = 1/(4*pi*(4*pi*1e-7));

[index_i, index_j, index_k, index_l, index_p, index_q] = ndgrid([0 1]);

index_sum = (-1).^(index_i+index_j+index_k+index_l+index_p+index_q);

calc_xyz = [true true true];

N = 1000;
tic
for ii = 1:N
  F = forces_calc_sumsum(rand(1,3),rand(1,3),rand(1,3),rand(1,3),rand(1,3));
end
toc


  function calc_out = forces_calc_sumsum(size1,size2,offset,J1,J2)
    
    J1 = J1(3);
    J2 = J2(3);
    
    if (J1==0 || J2==0)
      debug_disp('Zero magnetisation.')
      calc_out = [0; 0; 0];
      return;
    end
    
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
    
    
    if calc_xyz(1)
      component_x = index_sum.*component_x;
    else
      component_x = 0;
    end
    
    if calc_xyz(2)
      component_y = index_sum.*component_y;
    else
      component_y = 0;
    end
    
    if calc_xyz(3)
      component_z = index_sum.*component_z;
    else
      component_z = 0;
    end
    
    calc_out = J1*J2*magconst .* ...
      [ sum(component_x(:)) ;
        sum(component_y(:)) ;
        sum(component_z(:)) ] ;
    
    debug_disp(calc_out')
    
  end

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