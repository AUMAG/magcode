function []=Arrows(c,e,f,d,g,h,theta,M)
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

phiin = atan((h/2)/R1in);       % Angle of inside corner
alfain = atan((h/2)/R2in);      % Angle of outside corner

for i=1:M
    % U and V contain all the magnet coordinates [x1 y1;x2 y2;..,;xM yM]
    U = [L1*cos(theta*(i-1)-phi) L1*sin(theta*(i-1)-phi) ; L2*cos(theta*(i-1)-alfa) L2*sin(theta*(i-1)-alfa) ; L2*cos(theta*(i-1)+alfa) L2*sin(theta*(i-1)+alfa) ; L1*cos(theta*(i-1)+phi) L1*sin(theta*(i-1)+phi)];
    V = [L1in*cos(theta*(i-1)-phiin) L1in*sin(theta*(i-1)-phiin) ; L2in*cos(theta*(i-1)-alfain) L2in*sin(theta*(i-1)-alfain) ; L2in*cos(theta*(i-1)+alfain) L2in*sin(theta*(i-1)+alfain) ; L1in*cos(theta*(i-1)+phiin) L1in*sin(theta*(i-1)+phiin)];
  
    % Determine positions circle outer ring
    W(i,:)=c.*[cos(theta*(i-1)) sin(theta*(i-1))];  % Outer ring magnet positions
    T(i,:)=d.*[cos(theta*(i-1)) sin(theta*(i-1))];  % Inner ring magnet positions
    figure(99);
    hold on;
    % Plot outer and inner magnets
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
   arrowLength=g/2;                 % Arrow length
   strAnnotationType='arrow';       % Annotation type
   afStartingPoint=[T(i,1),T(i,2)]; % Arrow starting points
   afEndingPoint=[arrowLength*cos(-(i-1)*theta*(M/4-1))+T(i,1),arrowLength*sin(-(i-1)*theta*(M/4-1))+T(i,2)];   % Arrow ending points
   
   % Code to transform to correct coordinate system
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

   % Plot arrows
   handleToAnnotation1=annotation(...
   strAnnotationType,...
   [afStartingPoint_FU(1) afEndingPoint_FU(1)],...
   [afStartingPoint_FU(2) afEndingPoint_FU(2)]);
   i=i+1;
   end
  
   saveas(gcf, 'configuration.fig', 'fig');






end