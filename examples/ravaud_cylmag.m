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
  struct('turns',Nturns,'current',current,'dim',[r1 h1],'dir',[0 0 1]),...
  struct('magn',J2,'dim',[r2 h2],'dir',[0 0 1]),...
  displ'*[0 0 1]...
);

% Plot output

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
  % matlabfrag('fig/cylmag-matlab')
end

%% Save output to file

f = fopen('data/magcyl-matlab.txt','w');
fprintf(f,'%f %f\n',[1000*displ' fcyl(3,:)']');

%% Two magnets



% Magnet 1:
r1 = 0.02; % radius, metres
h1 = 0.02; % height, metres
grade1 = 'N35';

% Magnet 2:
r2 = 0.015; % radius, metres
h2 = 0.015; % height, metres
grade2 = 'N42';

% Displacement:
gap_max = 0.05; % metres
face_gap = linspace(0,gap_max);
displ = face_gap + h1/2 + h2/2;

% Calculate forces:
fcyl = magnetforces(...
  struct('grade',grade1,'dim',[r1 h1],'dir',[0 0 -1]),...
  struct('grade',grade2,'dim',[r2 h2],'dir',[0 0 1]),...
  displ'*[0 0 1]...
);

figure(1); clf; hold on
plot(1000*face_gap,fcyl(3,:))

xlabel('Face gap, mm')
ylabel('Force, N')



%% Benchmark tests
%
% First do a number of tests in a loop; this will be much slower
% (than Mathematica) because of the overheads in my magnetforces().

clc

NN = 1000;

Nturns = round(500*rand);
current = 10*rand;
J1 = 2*rand;
J2 = 2*rand;
r1 = 0.1*rand;
r2 = 0.1*rand;
h1 = 0.1*rand;
h2 = 0.1*rand;
displ = 0.1*rand([1 NN]);

f_rand = nan([NN 3]);

coil = struct('turns',Nturns,'current',current,'dim',[r1 h1],'dir',[0 0 1]);
magnet = struct('magn',J2,'dim',[r2 h2],'dir',[0 0 1]);

tic
for nn = 1:NN
  f_rand(nn,:) = magnetforces( coil, magnet, [0 0 displ(nn)] );
end
elapsedtime = toc;
disp(['Average time: ',num2str(1000*elapsedtime/NN),' ms (loop)']);

% But now compare this to the vectorised code, which avoids the overheads
% of processing magnetforces() multiple times. With vector operations,
% the performance is much improved.

tic
f_rand2 = magnetforces( coil, magnet, displ'*[0 0 1] );
elapsedtime = toc;
disp(['Average time: ',num2str(1000*elapsedtime/NN),' ms (no loop)']);

close_enough = @(x) round(x*1e6);
assert(all(close_enough(f_rand(:,3)) == close_enough(f_rand2(3,:))'),'loop and vector code must be equal')

%% Vectorised benchmark
%
% The vectorised code is so much faster that we should draw a graph

clc

Ntotal = 100;
Ndisp = round(logspace(0,5,Ntotal));

Ttime = nan([Ntotal 1]);

for nn = 1:Ntotal
  
  displ = 0.1*rand([1 Ndisp(nn)]);
  f_rand3 = nan([Ndisp(nn) 3]);
  
  tic
  f_rand3 = magnetforces( coil, magnet, displ'*[0 0 1] );
  elapsedtime = toc;
  
  Ttime(nn) = 1000*elapsedtime/Ndisp(nn);
  
end


try
  willfig('cylmag-bench','large'); clf; hold on
catch
  figure
end
plot(Ndisp,Ttime)
set(gca,'xscale','log','yscale','log')
xlabel('Number of calculations')
ylabel('Execution time per number of calculations, ms')


f = fopen('data/magcyl-matlab-benchmark.txt','w');
fprintf(f,'%f %f\n',[Ndisp' Ttime]');