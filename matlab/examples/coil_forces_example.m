
% Coils:
r1 = 0.02; % metres
h1 = 0.02; % metres

r2 = 0.015; % metres
h2 = 0.015; % metres

Nturns = 100;
current = 1; % Amperes

coil1 = magnetdefine('type','coil','turns',Nturns,'current',current,'dim',[r1 h1],'dir',[0 0 1]);
coil2 = magnetdefine('type','coil','turns',Nturns,'current',current,'dim',[r2 h2],'dir',[0 0 1]);

% Displacement:
NN = 45;
displ_range = 0.045;
displ = linspace(-displ_range,displ_range,NN);

% Calculate forces:
fcyl = magnetforces(coil1,coil2,[0;0;1]*displ);

% Plot output
figure(1);
plot(1000*displ,fcyl(3,:))
xlabel('Displacement, mm')
ylabel('Force, N')
