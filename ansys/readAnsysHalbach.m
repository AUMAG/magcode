function []=readAnsysHalbach()
clear all;close all;clc;
fid=fopen('NodeData.txt');

A1 = [];
A2 = [];
A3 = [];
A4 = [];

A_data=textscan(fid,'%f %f %f %f'); % Read NodeData.txt columns
A1 = [A1;A_data{1}];
A2 = [A2;A_data{2}];
A3 = [A3;A_data{3}];
A4 = [A4;A_data{4}];

% Remove lines with zero 'B'
condition=A4(:,1)==0;
A1(condition,:)=[];
A2(condition,:)=[];
A3(condition,:)=[];
A4(condition,:)=[];

% Scatter plot
caxis([min(A2) max(A2)]);
scatter(A2,A3,[],A4,'filled');

% Contour plot
figure(2);
x1=linspace(min(A2),max(A2),1000); % Don't make number of steps too big to avoid memory problems
y1=linspace(min(A3),max(A3),1000);
[X1 Y1]=meshgrid(x1,y1);
Z1=griddata(A2,A3,A4,X1,Y1);
contourf(X1,Y1,Z1,10);
end