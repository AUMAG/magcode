%% An example of drawing the magnets in 3D
%
% Akoun and Yonnet geometry with increased displacement.
% Drawing the magnets and forces overlaid onto the 3D geometry.

magnet_fixed = magnetdefine('type','cuboid',...
  'dim',[0.02 0.012 0.006],'magn',0.38,...
  'magdir','z','color',[0.1 0.1 0.8],'alpha',0.6);

magnet_float = magnetdefine('type','cuboid',...
  'dim',[0.012 0.02 0.006],'magn',0.38,...
  'magdir','x','color',[0 0.7 0.2],'alpha',0.4);

N = 501;
offset = repmat([-0.008; -0.008; 0.016],[1 N]);
displ = linspace(0, 0.03, N);
displ_range = offset+[1; 0; 0]*displ;

figure(1); clf; box on
figuresize(12,9,'cm')

f1_xyz = magnetforces(magnet_fixed,magnet_float,displ_range,...
  'draw',true,...
  'drawpath',true,...
  'drawforce',true,...
  'drawforcescale',0.025,...
  'markpath',{'d','o','s'},...
  'markpathN',5,...
  'markpathopt',{'color','black','markerfacecolor','black','markersize',10});

axis equal
view(-13,21)
grid on

xticklabels([])
yticklabels([])
zticklabels([])

set(gca,'position',[0.05 0.05 0.9 0.95],'projection','perspective')

print -dpng -r300 magforce.png


%% Magnetic fields

%%%% TODO BEFORE MOENA: proper user interface for magnetic fields!!!


%% Field calculations (side)

M = 1;

W = 5;
H = 5;


R = 1;
L = 4;
cylmag = magnetdefine('type','cylinder',...
  'radius',R,'height',L,'magn',M,...
  'dir','z','magdir',[1 0 0],...
  'alpha',1);

N = 1000;
rho_range = W*linspace(-1,1,N);
z_range   = H*linspace(-1,1,N);
phi = 0;

[rr,zz] = meshgrid(rho_range,z_range);

[B_rho, B_phi, B_z] = cylinder_field_transverse(M,R,L/2,rr,phi,zz);
Bmag = sqrt(B_rho.^2+B_phi.^2+B_z.^2);

figure(1); clf; hold on
figuresize(10,10,'cm')

h = surf(rr,zeros(size(Bmag)),zz,Bmag);
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;

% See: <https://github.com/keithfma/evenly_spaced_streamlines>
even_stream_arrow(rr,zz,B_rho,B_z,0.5,10,'Color','w','lineWidth',2,'Plane','xz');


N = 1000;
xx_range = W*linspace(-1,1,N);
yy_range = H*linspace(-1,1,N);

[xx,yy] = meshgrid(xx_range,yy_range);

rr = sqrt(xx.^2+yy.^2);
pp = atan2(yy,xx);
zz = 0;

[B_rho, B_phi, B_z] = cylinder_field_transverse(M,R,L/2,rr,pp,zz);

B_x = B_rho.*cos(pp)-B_phi.*sin(pp);
B_y = B_rho.*sin(pp)+B_phi.*cos(pp);

Bmag = sqrt(B_x.^2+B_y.^2+B_z.^2);

h = surf(xx,yy,zeros(size(xx)),Bmag);
h.EdgeColor = 'none';
h.FaceAlpha = 0.8;

plot(R*cos(linspace(0,2*pi)),R*sin(linspace(0,2*pi)),'k-','Linewidth',3)

even_stream_arrow(xx,yy,B_x,B_y,0.5,20,'Color','w','lineWidth',2,'Plane','xy');

%plot3([-W W W -W -W],[H H -H -H H],[0 0 0 0 0],'k-','linewidth',3)
%plot3([-W W W -W -W],[0 0 0 0 0],[H H -H -H H],'k-','linewidth',3)
%plot3([-W W],[-0.1 -0.1],[0.1 0.1],'k-','linewidth',3)

magnetdraw(cylmag,[0;0;0])

view(3);
axis equal;
axis off;
pbaspect([1 1 1])
set(gca,'position',[0.01 0.01 0.99 0.99])

print -dpng -r300 magfield.png

%%

function [B_rho, B_phi, B_z] = cylinder_field_transverse(M,R,L,rho,phi,z)
% CYLINDER_FIELD_TRANSVERSE  Magnetic field of transversally-magnetised cylinder
%
% These equations are from Caciagli (2018).

xi_1 = z+L;
xi_2 = z-L;
alpha_1 = (xi_1.^2+(rho+R).^2).^(-1/2);
alpha_2 = (xi_2.^2+(rho+R).^2).^(-1/2);
beta_1 = xi_1.*alpha_1;
beta_2 = xi_2.*alpha_2;
gamma = (rho-R)./(rho+R);
gg = 1-gamma.^2;
m_1 = (xi_1.^2+(rho-R).^2)./(xi_1.^2+(rho+R).^2);
m_2 = (xi_2.^2+(rho-R).^2)./(xi_2.^2+(rho+R).^2);

[K_1,E_1,PI_1] = ellipkepi(gg,(1-m_1));
[K_2,E_2,PI_2] = ellipkepi(gg,(1-m_2));

P11 = K_1-2./(1-m_1).*(K_1-E_1);
P12 = K_2-2./(1-m_2).*(K_2-E_2);

P31 = 1./(1-m_1).*(K_1-E_1) - gamma.^2./gg.*(PI_1-K_1);
P32 = 1./(1-m_2).*(K_2-E_2) - gamma.^2./gg.*(PI_2-K_2);

P41 = gamma./gg.*(PI_1.*(1+gamma.^2) - 2*K_1) - P11;
P42 = gamma./gg.*(PI_2.*(1+gamma.^2) - 2*K_2) - P12;

B_rho = 2e-7*M*R*cos(phi)./rho.*(beta_1.*P41-beta_2.*P42);
B_phi = 4e-7*M*R*sin(phi)./rho.*(beta_1.*P31-beta_2.*P32);
B_z   = 4e-7*M*R*cos(phi).*(alpha_1.*P11-alpha_2.*P12);

end


%% ellipkepi
%
% Calculate complete integrals of the first three kinds.
% If only the first two are needed, use the built-in Matlab function |ellipke| instead.


% \START
% \begin{mfunction}{ellipkepi}
% Complete elliptic integrals calculated with the arithmetric-geometric mean
% algorithms contained here: \url{http://dlmf.nist.gov/19.8}.
% Valid for $0\le a\le 1$ and $0\le m\le 1$.

function [K,E,PI] = ellipkepi(a,m)

a1 = 1;
g1 = sqrt(1-m);
p1 = sqrt(1-a);
q1 = 1;
w1 = 1;

nn = 0;
qq = 1;
ww = m;

while max(abs(w1(:))) > eps || max(abs(q1(:))) > eps

  % Update from previous loop
  a0 = a1;
  g0 = g1;
  p0 = p1;
  q0 = q1;

  % for Elliptic I
  a1 = (a0+g0)/2;
  g1 = sqrt(a0.*g0);

  % for Elliptic II
  nn = nn + 1;
  d1 = (a0-g0)/2;
  w1 = 2^nn*d1.^2;
  ww = ww + w1;

  % for Elliptic III
  rr = p0.^2+a0.*g0;
  p1 = rr./p0/2;
  q1 = q0.*(p0.^2-a0.*g0)./rr/2;
  qq = qq + q1;

end

K  = 1./a1*pi/2;
E  = K.*(1-ww/2);
PI = K.*(1+a./(2-2*a).*qq);

im = find(m == 1);
if ~isempty(im)
  K(im) = inf;
  E(im) = ones(length(im),1);
  PI(im) = inf;
end

end

% \end{mfunction}





function [xy, hl, ha] = even_stream_arrow(varargin)
% Compute and plot evenly-spaced streamlines for a vector field with
% arrow glyphs to indicate the flow direction. Uses the 'arrow' package by
% Dr. Erik A. Johnson from the Mathworks File Exchange.
%
% even_stream_arrow(xx, yy, uu, vv)
% even_stream_arrow(xx, yy, uu, vv, min_density)
% even_stream_arrow(xx, yy, uu, vv, min_density, max_density)
% even_stream_arrow(xy)
% even_stream_arrow(..., Name, Value)
% [xy, hh] = even_stream_arrow(...)
%
% Arguments:
%   xx, yy: Matrices or vectors, x-coord and y-coord. If matrices, the
%       size must match uu and vv. If vectors, xx must match the number of
%       columns in uu and vv, and yy must match the number of rows.
%   uu, vv: Matrices, vector field x-component and y-component
%   min_density: Scalar, specifies the minimum density (spacing) of
%       streamlines. From the streamslice() documentation: "modifies the
%       automatic spacing of the streamlines. DENSITY must be greater than
%       0. The default value is 1; higher values will produce more
%       streamlines on each plane. For example, 2 will produce
%       approximately twice as many streamlines while 0.5 will produce
%       approximately half as many."
%   max_density: Scalar, " " maximum " ", default is 2
%   xy: Matrix, [x, y] coordinates for stream line points, each row is a
%       point, individual lines are separated by NaNs.
%   hl = Graphics object for streamlines
%   ha = Vector of graphics objects for arrows
%
% Parameters (Name, Value):
%   'LineStyle': line style as in plot(), default = '-'
%   'LineWidth': line width as in plot(), default = 0.5
%   'Color': line color as in plot(), default = 'b'
%   'ArrowLength': arrow head length in pixels, default = 20
%   'ArrowTipAngle': arrow head tip angle in degrees, default = 20
%   'ArrowBaseAngle': arrow head base angle in degrees, default = 10
%   'ArrowDensity': arrow head density in arbitrary units, default
%       is 1, higher values will produce more closey spaced arrows
%  
% References: 
% [1] Jobard, B., & Lefer, W. (1997). Creating Evenly-Spaced Streamlines of
%   Arbitrary Density. In W. Lefer & M. Grave (Eds.), Visualization in
%   Scientific Computing ?97: Proceedings of the Eurographics Workshop in
%   Boulogne-sur-Mer France, April 28--30, 1997 (pp. 43?55). inbook,
%   Vienna: Springer Vienna. http://doi.org/10.1007/978-3-7091-6876-9_5
%
% Example: Plot and replot streamlines
%   [xx, yy] = meshgrid(0:0.2:2, 0:0.2:2);
%   uu = cos(xx).*yy;
%   vv = sin(xx).*yy;
%   subplot(1,2,1)
%   xy = even_stream_arrow(xx, yy, uu, vv, 1, 2, 'Color', 'b');
%   subplot(1,2,2)
%   even_stream_arrow(xy, 'Color', 'r', 'ArrowDensity', 3);
% %

% handle inputs
parser = inputParser;
parser.CaseSensitive = false;
parser.PartialMatching = false;
parser.KeepUnmatched = false;

parser.addRequired('xx_or_xy');
parser.addOptional('yy', []);
parser.addOptional('uu', []);
parser.addOptional('vv', []);
parser.addOptional('min_density', []);
parser.addOptional('max_density', []);
parser.addParameter('LineStyle', '-');
parser.addParameter('LineWidth', 0.5);
parser.addParameter('Color', 'b');
parser.addParameter('ArrowLength', 7);
parser.addParameter('ArrowTipAngle', 20);
parser.addParameter('ArrowBaseAngle', 45);
parser.addParameter('ArrowDensity', 1);
parser.addParameter('Plane', 'xy');

parser.parse(varargin{:});
xx_or_xy = parser.Results.xx_or_xy;
yy = parser.Results.yy;
uu = parser.Results.uu;
vv = parser.Results.vv;
min_density = parser.Results.min_density;
max_density = parser.Results.max_density;
line_style = parser.Results.LineStyle;
line_width = parser.Results.LineWidth;
color = parser.Results.Color;
arrow_length = parser.Results.ArrowLength;
arrow_tip_angle = parser.Results.ArrowTipAngle;
arrow_base_angle = parser.Results.ArrowBaseAngle;
arrow_density = parser.Results.ArrowDensity;

% get stream lines
optional = ~[isempty(yy), isempty(uu), isempty(vv), isempty(min_density), isempty(max_density)];
if all(optional)
    xx = xx_or_xy;
    xy = even_stream_data(xx, yy, uu, vv, min_density, max_density);
elseif ~any(optional)
    xy = xx_or_xy;
else
    error('Invalid input arguments, see help %s', mfilename)
end

% plot lines
switch parser.Results.Plane
  case 'xy'
    hl = plot3(xy(:,1), xy(:,2), zeros(size(xy(:,2))), ...
      'LineStyle', line_style, 'LineWidth', line_width, 'Color', color);

  case 'yz'
    hl = plot3(zeros(size(xy(:,2))), xy(:,1), xy(:,2), ...
      'LineStyle', line_style, 'LineWidth', line_width, 'Color', color);
    
  case 'xz'
    hl = plot3(xy(:,1), zeros(size(xy(:,2))), xy(:,2), ...
      'LineStyle', line_style, 'LineWidth', line_width, 'Color', color);
    
  otherwise
    error('Only "xy","yz","xz" allowed.')
end
    
end