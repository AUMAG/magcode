%% An example of force between two ring magnets

%% Setup

clear all

%% Geometry
%
% Of course these variables can be inserted directly into the
% |magnetdefine| function in the next section.

% magnet 1
r11 = 0.035; % inner radius
r12 = 0.045; % outer radius
h1  = 0.04;  % height

% magnet 2
r21 = 0.02;
r22 = 0.03;
h2  = 0.05;

% displacement vector
zmax = 0.15;
zdispl = linspace(-zmax,zmax);

%% Define magnets and calculate forces

ringmag1 = magnetdefine('type','cylinder','radius',[r11 r12],'height',h1,'grade','N42','dir','+z');
ringmag2 = magnetdefine('type','cylinder','radius',[r21 r22],'height',h2,'grade','N42','dir','-z');
 
Fring = magnetforces(ringmag1,ringmag2,[0; 0; 1]*zdispl);

%% Plot

figure(1); clf; hold on
box on

plot(1000*zdispl,Fring(3,:))
plot(xlim(),[0 0],'k--')
plot([0 0],ylim(),'k--')

xlabel('Axial displacement, mm')
ylabel('Axial force, N')
