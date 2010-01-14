%% Multipole array forces example
%
% Reproducing the linear Halbach array results of
% Allag, Yonnet & Latreche (2009).

%% Calculate

fixed_array = ...
  struct(...
  'type','linear-y',  ...
  'face','up',        ...
  'msize',  [0.01 0.01 0.01],   ...
  'Nmag_per_wave', 4, ...
  'Nmag', 5,          ...
  'magn', 1,          ...
  'magdir_first', -90 ...
  );

float_array = fixed_array;
float_array.face = 'down';
float_array.magdir_first = 90;

yrange = linspace(-0.08,0.08,201);
zgap = repmat([0; 0; 0.015],[1 length(yrange)]);
displ = zgap + [0; 1; 0]*yrange;

forces = multipoleforces(fixed_array, float_array, displ);

%% Plot

try
  willfig('allag-repro'); clf; hold on
catch
  figure
end

plot(yrange,forces([2 3],:))

set(gca,'box','on','ticklength',[0.02 0.05])
xlabel('Horizontal $y$-displacement, m')
ylabel('Force, N')
text(-0.018,40,'$F_z$','interpreter','latex');
text( 0.015,30,'$F_y$','interpreter','latex');

try
  colourplot
  axistight
  draworigin([0 0],'--')
  matlabfrag('fig/allag-repro','dpi',3200);
end
