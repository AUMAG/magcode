function []=readAnsysHalbach(c,d,e,f,g,h,theta,M)
%% This script reads and plots Ansys results from NodeData.txt
% clear all;
close all;%clc;
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
N = 200;            % Number of iterations

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

%% Plot magnet configuration
% Outside magnets
R1 = c-e/2;                 % Radius from center to midpoint inside of magnet
R2 = c+e/2;                 % Radius from center to midpoint outside of magnet
L1 = sqrt(R1^2+(f/2)^2);    % Radius from center to corner inside of magnet
L2 = sqrt(R2^2+(f/2)^2);    % Radius from center to corner outside of magnet
phi = atan((f/2)/R1);       % Angle of inside corner
alfa = atan((f/2)/R2);      % Angle of outside corner
% Inside magnets
R1in = d-g/2;                 % Radius from center to midpoint inside of magnet
R2in = d+g/2;                 % Radius from center to midpoint outside of magnet
L1in = sqrt(R1in^2+(h/2)^2);    % Radius from center to corner inside of magnet
L2in = sqrt(R2in^2+(h/2)^2);    % Radius from center to corner outside of magnet
delta = .01;
beta = .005;
phiin = atan((h/2)/R1in);       % Angle of inside corner
alfain = atan((h/2)/R2in);      % Angle of outside corner

gamma = atan(delta/c);
R2out = sqrt(c^2+delta^2);
ksi = atan(beta/c);
R3out = sqrt(c^2+beta^2);
    for i=1:M

    U = [L1*cos(theta*(i-1)-phi) L1*sin(theta*(i-1)-phi) ; L2*cos(theta*(i-1)-alfa) L2*sin(theta*(i-1)-alfa) ; L2*cos(theta*(i-1)+alfa) L2*sin(theta*(i-1)+alfa) ; L1*cos(theta*(i-1)+phi) L1*sin(theta*(i-1)+phi)];
    V = [L1in*cos(theta*(i-1)-phiin) L1in*sin(theta*(i-1)-phiin) ; L2in*cos(theta*(i-1)-alfain) L2in*sin(theta*(i-1)-alfain) ; L2in*cos(theta*(i-1)+alfain) L2in*sin(theta*(i-1)+alfain) ; L1in*cos(theta*(i-1)+phiin) L1in*sin(theta*(i-1)+phiin)];
  
% Determine positions circle outer ring

W(i,:)=c.*[cos(theta*(i-1)) sin(theta*(i-1))];  % Outer ring magnet positions
T(i,:)=d.*[cos(theta*(i-1)) sin(theta*(i-1))];  % Inner ring magnet positions
figure(10);
     hold on;
     patch(U(:,1),U(:,2),'w');
     patch(V(:,1),V(:,2),'w');
  
   hold off;
    axis equal;
    i=i+1;

    end
    %% Plot outer ring arrows
   for i=1:M
   arrowLength=e/2;
   strAnnotationType='arrow';
   afStartingPoint=[W(i,1),W(i,2)];
   afEndingPoint=[arrowLength*cos((i-1)*theta*(M/4+1))+W(i,1),arrowLength*sin((i-1)*theta*(M/4+1))+W(i,2)];
    
   set(gcf,'Units','normalized');
   afXAxisLimits=get(gca,'XLim');
   afYAxisLimits=get(gca,'YLim');
   afAxesDimensionsAndPositions=get(gca,'Position');
   fXAxisPosition=afAxesDimensionsAndPositions(1);
   fYAxisPosition=afAxesDimensionsAndPositions(2);
   fXAxisLength=afAxesDimensionsAndPositions(3);
   fYAxisLength=afAxesDimensionsAndPositions(4);
   fXonYaxesRatio=fXAxisLength/fYAxisLength;
   afFigurePosition=get(gcf,'Position');
   fXonYDimensionRatio=afFigurePosition(3);
   afStartingPoint_FU(1)=(afStartingPoint(1)-afXAxisLimits(1))/(afXAxisLimits(2)-afXAxisLimits(1))*fXAxisLength+fXAxisPosition;
   afStartingPoint_FU(2)=(afStartingPoint(2)-afYAxisLimits(1))/(afYAxisLimits(2)-afYAxisLimits(1))*fYAxisLength+fYAxisPosition;
   afEndingPoint_FU(1)=(afEndingPoint(1)-afXAxisLimits(1))/(afXAxisLimits(2)-afXAxisLimits(1))*fXAxisLength+fXAxisPosition;
   afEndingPoint_FU(2)=(afEndingPoint(2)-afYAxisLimits(1))/(afYAxisLimits(2)-afYAxisLimits(1))*fYAxisLength+fYAxisPosition;

   handleToAnnotation1=annotation(...
       strAnnotationType,...
       [afStartingPoint_FU(1) afEndingPoint_FU(1)],...
       [afStartingPoint_FU(2) afEndingPoint_FU(2)]);
   i=i+1;
   end
        % Plot inner ring arrows
   for i=1:M
   arrowLength=g/2;
   strAnnotationType='arrow';
   afStartingPoint=[T(i,1),T(i,2)];
   afEndingPoint=[arrowLength*cos(-(i-1)*theta*(M/4-1))+T(i,1),arrowLength*sin(-(i-1)*theta*(M/4-1))+T(i,2)];
    
   set(gcf,'Units','normalized');
   afXAxisLimits=get(gca,'XLim');
   afYAxisLimits=get(gca,'YLim');
   afAxesDimensionsAndPositions=get(gca,'Position');
   fXAxisPosition=afAxesDimensionsAndPositions(1);
   fYAxisPosition=afAxesDimensionsAndPositions(2);
   fXAxisLength=afAxesDimensionsAndPositions(3);
   fYAxisLength=afAxesDimensionsAndPositions(4);
   fXonYaxesRatio=fXAxisLength/fYAxisLength;
   afFigurePosition=get(gcf,'Position');
   fXonYDimensionRatio=afFigurePosition(3);
   afStartingPoint_FU(1)=(afStartingPoint(1)-afXAxisLimits(1))/(afXAxisLimits(2)-afXAxisLimits(1))*fXAxisLength+fXAxisPosition;
   afStartingPoint_FU(2)=(afStartingPoint(2)-afYAxisLimits(1))/(afYAxisLimits(2)-afYAxisLimits(1))*fYAxisLength+fYAxisPosition;
   afEndingPoint_FU(1)=(afEndingPoint(1)-afXAxisLimits(1))/(afXAxisLimits(2)-afXAxisLimits(1))*fXAxisLength+fXAxisPosition;
   afEndingPoint_FU(2)=(afEndingPoint(2)-afYAxisLimits(1))/(afYAxisLimits(2)-afYAxisLimits(1))*fYAxisLength+fYAxisPosition;

   handleToAnnotation1=annotation(...
       strAnnotationType,...
       [afStartingPoint_FU(1) afEndingPoint_FU(1)],...
       [afStartingPoint_FU(2) afEndingPoint_FU(2)]);
   i=i+1;
   end
    %% Plot of BSUM exactly between inner and outer ring
figure(1);
plot(t,A4(E));
figure(2);
l=linspace(1,200,200);
polar(t,A4(E(l))');

%% Scatter plot
figure(3);
caxis([min(A2) max(A2)]);
scatter(A2,A3,[],A4,'filled');
axis equal;
title('Magnetic flux B on points');
xlabel('x,m');
ylabel('y,m');


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
figure(5);
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


end