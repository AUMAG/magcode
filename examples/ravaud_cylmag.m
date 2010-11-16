%% Matlab example for forces between cylindrical magnets/coils
%
% The script below uses the "magnetforces" function in the "matlab/"
% sibling directory. This function may be used to calculate the force
% between cylindrical and cuboid magnets with a common and consistent
% interface.


%% Define geometry and properties
%
% The magnet moves inside the coil and we want to calculate the force from
% the coil on the magnet.

% Magnet:
r2 = 0.015; % metres
h2 = 0.015; % metres
J2 = 1; % Tesla

% Coil:
r1 = 0.02; % metres
h1 = 0.02; % metres

Nturns = 100;
current = 1; % Amperes
mu0 = 4*pi*10^(-7);
J1 = mu0*Nturns*current/h1;

% Displacement:
NN = 45;
displ_range = 0.045;
displ = linspace(-displ_range,displ_range,NN);

% Calculate forces:
fcyl = magnetforces(...
  struct('turns',Nturns,'current',current,'dim',[r1 h1],'dir',[0 0  1]),...
  struct('magn',J2,'dim',[r2 h2],'dir',[0 0 1]),...
  displ'*[0 0 1]...
);

%% Plot output

try
  willfig('cylmag','tiny'); clf; hold on
catch
  figure
end

plot(1000*displ,fcyl(3,:))

xlabel('Displacement, mm')
ylabel('Force, N')
title('Matlab')

try
  colourplot
  axistight
  draworigin
  matlabfrag('fig/cylmag-matlab')
end

%% Save output to file

f = fopen('data/magcyl-matlab.txt','w');
fprintf(f,'%f %f\n',[1000*displ' fcyl(3,:)']');