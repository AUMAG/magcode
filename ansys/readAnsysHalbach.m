function []=readAnsysHalbach()
%% This script reads and plots Ansys results from NodeData.txt
clear all;close all;clc;
fid=fopen('NodeData.txt');

A1 = [];    % Node number
A2 = [];    % Node x position
A3 = [];    % Node y position
A4 = [];    % Node magnetic B sum
A5 = [];    % Node magnetic B x direction
A6 = [];    % Node magnetic B y direction

A_data=textscan(fid,'%f %f %f %f %f %f'); % Read NodeData.txt columns
A1 = [A1;A_data{1}];
A2 = [A2;A_data{2}];
A3 = [A3;A_data{3}];
A4 = [A4;A_data{4}];
A5 = [A5;A_data{5}];
A6 = [A6;A_data{6}];

% Remove lines with zero 'B'
condition=A4(:,1)==0;
A1(condition,:)=[];
A2(condition,:)=[];
A3(condition,:)=[];
A4(condition,:)=[];
A5(condition,:)=[];
A6(condition,:)=[];

% Scatter plot
figure(1);
caxis([min(A2) max(A2)]);
scatter(A2,A3,[],A4,'filled');
title('Magnetic flux B on points');
xlabel('x,m');
ylabel('y,m');


% Contour plot
figure(2);
x1=linspace(min(A2),max(A2),1000); % Don't make number of steps too big to avoid memory problems
y1=linspace(min(A3),max(A3),1000);
[X1 Y1]=meshgrid(x1,y1);
Z1=griddata(A2,A3,A4,X1,Y1);
contour(X1,Y1,Z1,30);
title('Contour plot of magnetic flux B');
xlabel('x,m');
ylabel('y,m');

% Velocity plot
figure(3);
figure(3);
quiver(A2,A3,A5,A6,3);
title('Magnetic flux B velocity vectors');
xlabel('x,m');
ylabel('y,m');


% Velocity plot with colored arrows
figure(4);
title('Magnetic flux B velocity vectors');
xlabel('x,m');
ylabel('y,m');figure(4);
quiverc(A2,A3,A5,A6,3);
colorbar;

end