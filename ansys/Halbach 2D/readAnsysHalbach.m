function []=readAnsysHalbach(c,d,e,f,g,h,theta,M)
%% This script reads and plots Ansys results from NodeData.txt
% clear all;
%close all;clc;
fid=fopen('NodeData.txt');


A1 = [];    % Node number
A2 = [];    % Node x position
A3 = [];    % Node y position
A4 = [];    % Node magnetic B sum
A5 = [];    % Node magnetic B x direction
A6 = [];    % Node magnetic B y direction
A7 = [];    % Magnetic vector potential

A_data=textscan(fid,'%f %f %f %f %f %f %f'); % Read NodeData.txt columns
A1 = [A1;A_data{1}];
A2 = [A2;A_data{2}];
A3 = [A3;A_data{3}];
A4 = [A4;A_data{4}];
A5 = [A5;A_data{5}];
A6 = [A6;A_data{6}];
A7 = [A7;A_data{7}];

% Remove lines with zero 'B'
condition=A4(:,1)==0;
A1(condition,:)=[];
A2(condition,:)=[];
A3(condition,:)=[];
A4(condition,:)=[];
A5(condition,:)=[];
A6(condition,:)=[];
A7(condition,:)=[];

%% Determine magnetic field B between inner and outer ring
N = 1000;            % Number of iterations

nodes = [A2 A3];
t = linspace(0,2*pi,N);
R = ((c-e/2)+(d+g/2))/2;
p(:,1) = R*cos(t);
p(:,2) = R*sin(t);
D = zeros(N,1);
E = zeros(N,1);

for i=1:N
    D(i,:) = norm(nodes-repmat(p(i,:),[length(nodes) 1]));
    displ= nodes-repmat(p(i,:),[length(nodes) 1]);
    D = sqrt(displ(:,1).^2+displ(:,2).^2);
    E(i,:) = find(D==min(D));
end
i = i+1;
save('polardata.mat','E','t','A4');


%{
% Arrows(c,e,f,d,g,h,theta,M);
    %% Plot of BSUM exactly between inner and outer ring

figure(1);
plot(t,A4(E));
figure(2);
l=linspace(1,N,N);
polar(t,A4(E(l))');

%% Scatter plot


figure(3);
hold on;

Arrows(c,e,f,d,g,h,theta,M);
scatter(A2,A3,[],A4,'filled');

axis equal;


title('Magnetic flux B on points');
xlabel('x,m');
ylabel('y,m');
hold off;

%% Contour plot of magnetic BSUM

figure(4);
x1=linspace(min(A2),max(A2),1000); % Don't make number of steps too big to avoid memory problems
y1=linspace(min(A3),max(A3),1000);
[X1 Y1]=meshgrid(x1,y1);
Z1=griddata(A2,A3,A4,X1,Y1);
Z2=griddata(A2,A3,A7,X1,Y1);
contour(X1,Y1,Z1,30);
axis equal;
title('Contour plot of magnetic flux B');
xlabel('x,m');
ylabel('y,m');

%% Contour plot of magnetic flux potential (2D flux lines)
figure(5)
Z2=griddata(A2,A3,A7,X1,Y1);
contour(X1,Y1,Z2,30);
axis equal;
title('Contour plot of 2D flux lines (AZ)');
xlabel('x,m');
ylabel('y,m');

% Velocity plot
figure(6);
quiver(A2,A3,A5,A6,3);
axis equal;
title('Magnetic flux B velocity vectors');
xlabel('x,m');
ylabel('y,m');


% % Velocity plot with colored arrows
% figure(7);
% title('Magnetic flux B velocity vectors');
% xlabel('x,m');
% ylabel('y,m');
% quiverc(A2,A3,A5,A6,3);
% colorbar;
% axis equal
%}

end