
clear all

magnet_fixed.dim = [0.04 0.04 0.04];
magnet_float.dim = magnet_fixed.dim;

magnet_fixed.type = 'cuboid';
magnet_float.type = 'cuboid';

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;
offset = 0.1;

f = [];

for ii = [1, -1]
  magnet_fixed.magdir = [0 0 ii];
  for jj = [1, -1]
    magnet_float.magdir = [0 0 jj];
    for kk = [1, -1]
      displ = kk*[0 0 offset];
      f(:,end+1) = magnetforces(magnet_fixed,magnet_float,displ);
    end
  end
end

dirforces = round( f(3,:), 8 );
otherforces = f([1 2],:);


assert ( ...
             all( abs( otherforces(:) ) < 1e-11 ) , ...
             'Orthogonal forces should be zero' ...
           )
assert ( ...
             all( abs(dirforces) == abs(dirforces(1)) ) , ...
             'Force magnitudes should be equal' ...
           )
assert ( ...
             all( dirforces(1:4) == -dirforces(5:8) ) , ...
             'Forces should be opposite with reversed fixed magnet magnetisation' ...
           )
assert ( ...
             all( dirforces([1 3 5 7]) == -dirforces([2 4 6 8]) ) , ...
             'Forces should be opposite with reversed float magnet magnetisation' ...
           )




f = [];

for ii = [1, -1]
  magnet_fixed.magdir = [ii 0 0]; % $\pm x$
  for jj = [1, -1]
    magnet_float.magdir = [jj 0 0];
    for kk = [1, -1]
      displ = kk*[offset 0 0];
      f(:,end+1) = magnetforces(magnet_fixed,magnet_float,displ);
    end
  end
end

dirforces = round( f(1,:), 8 );
otherforces = f([2 3],:);


assert ( ...
             all( abs( otherforces(:) ) < 1e-11 ) , ...
             'Orthogonal forces should be zero' ...
           )
assert ( ...
             all( abs(dirforces) == abs(dirforces(1)) ) , ...
             'Force magnitudes should be equal' ...
           )
assert ( ...
             all( dirforces(1:4) == -dirforces(5:8) ) , ...
             'Forces should be opposite with reversed fixed magnet magnetisation' ...
           )
assert ( ...
             all( dirforces([1 3 5 7]) == -dirforces([2 4 6 8]) ) , ...
             'Forces should be opposite with reversed float magnet magnetisation' ...
           )




f = [];

for ii = [1, -1]
  magnet_fixed.magdir = [0 ii 0]; % $\pm y$
  for jj = [1, -1]
    magnet_float.magdir = [0 jj 0];
    for kk = [1, -1]
      displ = kk*[0 offset 0];
      f(:,end+1) = magnetforces(magnet_fixed,magnet_float,displ);
    end
  end
end

dirforces = round( f(2,:), 8 );
otherforces = f([1 3],:);



assert ( ...
             all( abs( otherforces(:) ) < 1e-11 ) , ...
             'Orthogonal forces should be zero' ...
           )
assert ( ...
             all( abs(dirforces) == abs(dirforces(1)) ) , ...
             'Force magnitudes should be equal' ...
           )
assert ( ...
             all( dirforces(1:4) == -dirforces(5:8) ) , ...
             'Forces should be opposite with reversed fixed magnet magnetisation' ...
           )
assert ( ...
             all( dirforces([1 3 5 7]) == -dirforces([2 4 6 8]) ) , ...
             'Forces should be opposite with reversed float magnet magnetisation' ...
           )




