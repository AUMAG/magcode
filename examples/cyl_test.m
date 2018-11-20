%% Matlab example for forces between cylindrical magnets/coils
%
% The script below uses the "magnetforces" function in the "matlab/"
% sibling directory. This function may be used to calculate the force
% between cylindrical and cuboid magnets with a common and consistent
% interface.


%% Geometry one

r1 = 0.1/2; % metres
h1 = 0.03; % metres
J1 = 1.3; % Tesla
mag1 = struct('magn',J1,'dim',[r1 h1],'magdir',[0 0  1]);

r2 = 0.075/2; % metres
h2 = 0.03; % metres
J2 = 1.3; % Tesla
mag2 = struct('magn',J2,'dim',[r2 h2],'magdir',[0 0 -1]);

% Displacement:
NN = 51;
displ_range = 0.1;
face_gap = linspace(0,displ_range,NN);
displ = h1/2+h2/2+face_gap;

% Calculate a single force:
L = 0.05;
fcyl = magnetforces(mag1,mag2,[0 0 h1/2+h2/2+L]...
);
fprintf('With face gap %g mm, force is %2.2f N\n',1000*L,fcyl(3));

% Calculate over a range:
fcyl = magnetforces(mag1,mag2,displ'*[0 0 1]);

figure(1)
plot(1000*(face_gap),fcyl(3,:))
xlabel('Face gap, mm')
ylabel('Force, N')



%% Geometry two

r1 = 0.15/2; % metres
h1 = 0.06; % metres
J1 = 1.3; % Tesla
mag1 = struct('magn',J1,'dim',[r1 h1],'magdir',[0 0  1]);

r2 = 0.07/2; % metres
h2 = 0.06; % metres
J2 = 1.3; % Tesla
mag2 = struct('magn',J2,'dim',[r2 h2],'magdir',[0 0 -1]);

% Displacement:
NN = 51;
displ_range = 0.1;
face_gap = linspace(0,displ_range,NN);
displ = h1/2+h2/2+face_gap;

% Calculate a single force:
L = 0.05;
fcyl = magnetforces(mag1,mag2,[0 0 h1/2+h2/2+L]);
fprintf('With face gap %g mm, force is %2.2f N\n',1000*L,fcyl(3));

% Calculate over a range:
fcyl = magnetforces(mag1,mag2,displ'*[0 0 1]);

figure(2)
plot(1000*(face_gap),fcyl(3,:))

xlabel('Face gap, mm')
ylabel('Force, N')