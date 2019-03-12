
clear all

magnet_fixed.dim = [0.04 0.04 0.04];
magnet_float.dim =  magnet_fixed.dim;

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.type = 'cuboid';
magnet_float.type = 'cuboid';

f = [];

for ii = [-1 1]
  for jj = [-1 1]
    for xx = 0.12*[-1, 1]
      for yy = 0.12*[-1, 1]
        for zz = 0.12*[-1, 1]

              magnet_fixed.magdir = [0 0 ii];
              magnet_float.magdir = [0 0 jj];
              displ = [xx yy zz];
              f(:,end+1) = magnetforces(magnet_fixed,magnet_float,displ);

        end
      end
    end
  end
end

f = round(f,8);

uniquedir = f(3,:);
otherdir  = f([1 2],:);

test1 = abs(diff(abs(f(1,:))))<1e-10 ;
test2 = abs(diff(abs(f(2,:))))<1e-10 ;
test3 = abs(diff(abs(f(3,:))))<1e-10 ;
assert ( all(test1) && all(test2) && all(test3) , ...
         'All forces in a single direction should be equal' )

test = abs(diff(abs(otherdir))) < 1e-11;
assert ( all(test) , 'Orthogonal forces should be equal' )

test1 = f(:,1:8) == f(:,25:32);
test2 = f(:,9:16) == f(:,17:24);
assert ( all( test1(:) ) && all( test2(:)) , ...
             'Reverse magnetisation shouldn''t make a difference' )


f = [];

for ii = [-1 1]
  for jj = [-1 1]
    for xx = 0.12*[-1, 1]
      for yy = 0.12*[-1, 1]
        for zz = 0.12*[-1, 1]

              magnet_fixed.magdir = [0 0 ii];
              magnet_float.magdir = [0 jj 0];
              displ = [xx yy zz];
              f(:,end+1) = magnetforces(magnet_fixed,magnet_float,displ);

        end
      end
    end
  end
end

f = round( f , 8 );

uniquedir = f(1,:);
otherdir  = f([2 3],:);


test1 = abs(diff(abs(f(1,:))))<1e-10 ;
test2 = abs(diff(abs(f(2,:))))<1e-10 ;
test3 = abs(diff(abs(f(3,:))))<1e-10 ;
assert ( all(test1) && all(test2) && all(test3) , ...
         'All forces in a single direction should be equal' )

test = abs(diff(abs(otherdir))) < 1e-11;
assert ( all(test) , 'Orthogonal forces should be equal' )

test1 = f(:,1:8) == f(:,25:32);
test2 = f(:,9:16) == f(:,17:24);
assert ( all( test1(:) ) && all( test2(:)) , ...
             'Reverse magnetisation shouldn''t make a difference' )




