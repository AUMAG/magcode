%% Setup
%
% In case you don't have the various bits'n'pieces that I use to create
% my Matlab graphics (probably likely).

if ~exist('willfig','file')
  close all
  willfig = @(str) figure;
  simple_graph = true;
else
  simple_graph = false;
end

%% Varying magnet aspect ratios

magnet_fixed.magn = 1;
magnet_float.magn = 1;

magnet_fixed.magdir = [0 0 1];
magnet_float.magdir = [0 0 -1];

displ_steps = 51;

%%

volume = 0.01^3;
ratios = [0.25 0.5 1 2];

height = @(a)    (volume*a^2)^(1/3);
face_size = @(a) (volume/a)^(1/3);
magsize = @(a) ...
  [ face_size(a) face_size(a) height(a) ];


willfig('mag_ratios'); clf; hold on

for this_ratio = ratios

  magnet_fixed.dim = magsize(this_ratio);
  magnet_float.dim = magsize(this_ratio);
  
  displ = 2*height(this_ratio)+ linspace(0, 0.015, displ_steps);
  displ_range = [0; 0; 1]*displ;
  
  these_forces = magnetforces(magnet_fixed,magnet_float,displ_range);
  
  face_gap = displ - displ(1);
  plot(face_gap,these_forces(3,:),'Tag',num2str(this_ratio));

end

xlabel('Face gap, m')
ylabel('Vertical force, N')
set(gca,'box','on','ticklength',[0.02 0.05])
axis tight

%%

if ~simple_graph
  colourplot
  labelplot('east','vertical','Ratio $\hwratio$')
  legendshrink
  matlabfrag('fig/mag_ratio','dpi',3200);
end
