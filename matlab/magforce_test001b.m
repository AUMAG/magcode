 
 
disp('=================') 
fprintf('TEST 001b: ') 
 
magnet_fixed.dim  =  [0.04 0.04 0.04]; 
magnet_float.dim  =   magnet_fixed.dim; 
 
magnet_fixed.magn  =  1.3; 
magnet_float.magn  =  1.3; 
 
 
 
fzyz  =  []; 
 
for ii  =  [1, -1] 
  for jj  =  [1, -1] 
    for kk  =  [1, -1] 
 
      magnet_fixed.magdir  =  ii * [0  90];  % $\pm z$ 
      magnet_float.magdir  =  jj * [90 0];  % $\pm y$ 
      displ  =  kk * [0 0 0.1];  % $\pm z$ 
      fzyz(:, end+1)  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
    end 
  end 
end 
 
dirforces  =  chop( fzyz(2,:), 8 ); 
otherforces  =  fzyz([1 3],:); 
 
 
 
 
assert (  ... 
      all( abs( otherforces(:) ) < 1e-11 ) ,  ... 
      'Orthogonal forces should be zero'  ... 
    ) 
assert (  ... 
      all( abs(dirforces) == abs(dirforces(1)) ) ,  ... 
      'Force magnitudes should be equal'  ... 
    ) 
assert (  ... 
      all( dirforces(1:4) == -dirforces(5:8) ) ,  ... 
      'Forces should be opposite with reversed fixed magnet magnetisation'  ... 
    ) 
assert (  ... 
      all( dirforces([1 3 5 7]) == -dirforces([2 4 6 8]) ) ,  ... 
      'Forces should be opposite with reversed float magnet magnetisation'  ... 
    ) 
 
 
 
 
 
 
fzxz  =  []; 
 
for ii  =  [1, -1] 
  for jj  =  [1, -1] 
    for kk  =  [1, -1] 
 
      magnet_fixed.magdir  =  ii * [0  90];  % $\pm z$ 
      magnet_float.magdir  =  [90+jj * 90 0];  % $\pm x$ 
      displ  =  kk * [0.1 0 0];  % $\pm x$ 
      fzxz(:, end+1)  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
    end 
  end 
end 
 
dirforces  =  chop( fzxz(3,:), 8 ); 
otherforces  =  fzxz([1 2],:); 
 
 
 
 
assert (  ... 
      all( abs( otherforces(:) ) < 1e-11 ) ,  ... 
      'Orthogonal forces should be zero'  ... 
    ) 
assert (  ... 
      all( abs(dirforces) == abs(dirforces(1)) ) ,  ... 
      'Force magnitudes should be equal'  ... 
    ) 
assert (  ... 
      all( dirforces(1:4) == -dirforces(5:8) ) ,  ... 
      'Forces should be opposite with reversed fixed magnet magnetisation'  ... 
    ) 
assert (  ... 
      all( dirforces([1 3 5 7]) == -dirforces([2 4 6 8]) ) ,  ... 
      'Forces should be opposite with reversed float magnet magnetisation'  ... 
    ) 
 
 
 
 
 
 
fzxx  =  []; 
 
for ii  =  [1, -1] 
  for jj  =  [1, -1] 
    for kk  =  [1, -1] 
 
      magnet_fixed.magdir  =  ii * [0  90];  % $\pm z$ 
      magnet_float.magdir  =  [90+jj * 90 0];  % $\pm x$ 
      displ  =  kk * [0 0 0.1];  % $\pm z$ 
      fzxx(:, end+1)  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
    end 
  end 
end 
 
dirforces  =  chop( fzxx(1,:), 8 ); 
otherforces  =  fzxx([2 3],:); 
 
 
 
 
 
 
assert (  ... 
      all( abs( otherforces(:) ) < 1e-11 ) ,  ... 
      'Orthogonal forces should be zero'  ... 
    ) 
assert (  ... 
      all( abs(dirforces) == abs(dirforces(1)) ) ,  ... 
      'Force magnitudes should be equal'  ... 
    ) 
assert (  ... 
      all( dirforces(1:4) == -dirforces(5:8) ) ,  ... 
      'Forces should be opposite with reversed fixed magnet magnetisation'  ... 
    ) 
assert (  ... 
      all( dirforces([1 3 5 7]) == -dirforces([2 4 6 8]) ) ,  ... 
      'Forces should be opposite with reversed float magnet magnetisation'  ... 
    ) 
 
 
 
 
 
 
fzyy  =  []; 
 
for ii  =  [1, -1] 
  for jj  =  [1, -1] 
    for kk  =  [1, -1] 
 
      magnet_fixed.magdir  =  ii * [0  90];  % $\pm z$ 
      magnet_float.magdir  =  jj * [90 0];  % $\pm y$ 
      displ  =  kk * [0 0.1 0];  % $\pm y$ 
      fzyy(:, end+1)  =  magnetforces(magnet_fixed,magnet_float,displ); 
 
    end 
  end 
end 
 
dirforces  =  chop( fzyy(3,:), 8 ); 
otherforces  =  fzyy([1 2],:); 
 
 
 
 
assert (  ... 
      all( abs( otherforces(:) ) < 1e-11 ) ,  ... 
      'Orthogonal forces should be zero'  ... 
    ) 
assert (  ... 
      all( abs(dirforces) == abs(dirforces(1)) ) ,  ... 
      'Force magnitudes should be equal'  ... 
    ) 
assert (  ... 
      all( dirforces(1:4) == -dirforces(5:8) ) ,  ... 
      'Forces should be opposite with reversed fixed magnet magnetisation'  ... 
    ) 
assert (  ... 
      all( dirforces([1 3 5 7]) == -dirforces([2 4 6 8]) ) ,  ... 
      'Forces should be opposite with reversed float magnet magnetisation'  ... 
    ) 
 
 
 
 
fprintf('passed\n') 
disp('=================') 
 

