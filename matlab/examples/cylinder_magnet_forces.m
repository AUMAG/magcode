%% Examples of cylindrical magnets/coils

%% Init

close all
clear all

%% Geometry

% Coil:

r1 = 0.02; % metres
h1 = 0.02; % metres

Nturns = 100;
current = 1; % Amperes
mu0 = 4*pi*10^(-7);
J1 = mu0*Nturns*current/h1;

coil_cyl_z = magnetdefine('type','cylinder','turns',Nturns,'current',current,'dim',[r1 h1],'dir',[0 0 1]);

% Magnet:

r2 = 0.015; % metres
h2 = 0.015; % metres
J2 = 1; % Tesla

mag_cyl_z = magnetdefine('type','cylinder','magn',J2,'dim',[r2 h2],'dir',[0 0 1]);

% "Equivalent" cuboid magnet:

x2 = r2*sqrt(pi);
y2 = r2*sqrt(pi);
h2 = 0.015; % metres

V1 = pi*r2^2*h2;
V2 = x2*y2*h2;
assert(abs(V1-V2)<eps,'Volume of cyl and cuboid not same.')

J2 = 1; % Tesla
mag_cuboid_z = magnetdefine('type','cuboid','magn',J2,'dim',[x2 y2 h2],'magdir',[0 0 1]);


%% 1. Coaxial example
%
% The magnet moves inside the coil and we want to calculate the force from
% the fixed coil on the moving magnet.
% Coaxial cylindrical calculations are relatively fast.

% Displacement:
NN = 50;
displ_max = 0.045;
displ_range = linspace(-displ_max,displ_max,NN);
displ = displ_range'*[0 0 1];

% Calculate forces:
tic
fcyl = magnetforces(coil_cyl_z,mag_cyl_z,displ);
toc

% Plot output
figure(1); clf; hold on
plot(1000*displ_range,fcyl(3,:))
box on
xlabel('Z displacement, mm')
ylabel('Force, N')


%% 2. Fixed eccentricity
%
% The magnet moves with a small offset inside the coil and we want to
% calculate the force from the fixed coil on the moving magnet.
% Non-coaxial cylindrical calculations are quite slow. 

% Displacement:
NN = 50;
displ_max = 0.045;
ecc_offset = (r1-r2)/2;
displ_range = linspace(-displ_max,displ_max,NN);
displ = displ_range'*[0 0 1];
ecc = repmat(ecc_offset,[NN,1])*[1 0 0];

% Calculate forces:
tic
fcyl2 = magnetforces(coil_cyl_z,mag_cyl_z,displ+ecc);
toc

% Plot output
figure(1);
plot(1000*displ_range,fcyl2(3,:))
plot(1000*displ_range,fcyl2(1,:))
box on
xlabel('Z displacement, mm')
ylabel('Force, N')

legend('Coaxial, Z force','Eccentric, Z force','Eccentric, X force')


%% 3. Two z magnets, x displacement
%
% Now look at two magnets magnetised in the Z direction with an offset
% between then in the Z direction. Then displace the floating magnet in the
% X direction.
% The cylindrical magnet calcs are around 100 times slower than the cuboid.

% Displacement:
NN = 50;
displ_max = 0.045;
displ_range = linspace(-displ_max,displ_max,NN);
displ = displ_range'*[1 0 0];
offset = repmat(h2*1.5,[NN,1])*[0 0 1];

% Calculate forces:
tic
fcyl3 = magnetforces(mag_cyl_z,mag_cyl_z,displ+offset);
toc

% Plot output
figure(3); clf; hold on

plot(1000*displ_range,fcyl3(3,:))
plot(1000*displ_range,fcyl3(1,:))
box on
xlabel('X displacement, mm')
ylabel('Force, N')

tic
fcyl4 = magnetforces(mag_cuboid_z,mag_cuboid_z,displ+offset);
toc
plot(1000*displ_range,fcyl4(3,:),'--')
plot(1000*displ_range,fcyl4(1,:),'--')

legend('Cyl Z', 'Cyl X','Cuboid Z', 'Cuboid X')
