
clear all

% Fixed parameters

fixed_array = ...
  struct(...
        'length', [0.10 0.10], ...
        'width',  0.10, ...
        'height', 0.01, ...
        'Nmag_per_wave', [4 4], ...
        'Nwaves', [1 1], ...
        'magn', 1, ...
        'magdir_first', [90 90] ...
  );

float_array = fixed_array;
float_array.magdir_first = [-90 -90];

f = nan([3 0]);

% The varying calculations

fixed_array.type = 'planar';
float_array.type = fixed_array.type;
fixed_array.align = 'xy';
float_array.align = fixed_array.align;
fixed_array.face = 'up';
float_array.face = 'down';
displ = [0 0 0.02];
f(:,end+1) = multipoleforces(fixed_array, float_array, displ);

fixed_array.type = 'planar';
float_array.type = fixed_array.type;
fixed_array.align = 'yz';
float_array.align = fixed_array.align;
fixed_array.face = '+x';
float_array.face = '-x';
displ = [0.02 0 0];
f(:,end+1) = multipoleforces(fixed_array, float_array, displ);

fixed_array.type = 'planar';
float_array.type = fixed_array.type;
fixed_array.align = 'xz';
float_array.align = fixed_array.align;
fixed_array.face = '+y';
float_array.face = '-y';
displ = [0 0.02 0];
f(:,end+1) = multipoleforces(fixed_array, float_array, displ);

ind = [3 4 8];

assert( all(round(f(ind) * 100)/100==589.05) ,  ...
  'Arrays aligned in different directions should produce consistent results.');

assert( all(f(~ind)<1e-10) ,  ...
  'These forces should all be (essentially) zero.');


