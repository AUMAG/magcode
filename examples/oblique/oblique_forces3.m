function [ forces, torques ] = oblique_forces3( varargin )
%OBLIQUE_FORCES
%
% This function calculates the forces generated on a magnetic spring
% consisting of two pairs of inclined permanent magnets.
%
% INPUTS:         DEFAULT  DESCRIPTION
%
%  'magn'          1T      magnet magnetisation
%  'unitlength'  * 0.01    cube root of the magnet volume, metres
%  'displ    '   * [0;0;0] displacement
%  'dispoffset'    [0;0;0] offset added to each displacement
%  'magangle'    * 45      inclination angle(s) of the magnets, degrees
%  'magratio'    * 0.4     ratio(s) between magnet height/weight
%  'gapratio'    * 0       horizontal offset(s) at zero displacement
%                            in terms of 'unitlength' 
%  'leverratio'  * 2       horizontal offset(s) between magnet centres and
%                            centre of rotation of the floating section
%  'rotation'    * 0       rotation of the floating section in degrees
%                            (around midpoint between magnet centres)
%  'plot'          false   whether to plot the geometry
%  'plotvecscale'  0.01    scaling factor for force vectors on the plot
%  'plotextras'    false   whether to show labels etc on the plot
%
% Starred inputs may take vector inputs.
%
% OUTPUTS:
%
%  DISPL   max size: [  <magangles> <gaps> <mag ratios> <length ratio> <rotation> <length displ>]
%  FORCES  max size: [3 <magangles> <gaps> <mag ratios> <length ratio> <rotation> <length displ>]
%
% Outputs are SQUEEZEd, so non-vector input will reduce the number
% of dimensions in the output arrays.
%
% Coordinate system:
%    X - "horizontal" in the plane of the magnet rotations
%    Y - "vertical"
%    Z - "out of plane"
%

%% Inputs and initial setup

p = inputParser;
p.addParamValue('magn',1);
p.addParamValue('unitlength',0.01);
p.addParamValue('magratio',0.4);
p.addParamValue('displ',[0; 0; 0]);
p.addParamValue('magangle',45);
p.addParamValue('gapratio',0);
p.addParamValue('dispoffset',[0; 0; 0]);
p.addParamValue('leverratio',2);
p.addParamValue('rotation',0);
p.addParamValue('plot',false);
p.addParamValue('plotsize',0);
p.addParamValue('plotextras',false);
p.addParamValue('plotvecscale',0.01);
p.addParamValue('plottorquescale',0);
p.addParamValue('plottorquearc',pi/4);
p.addParamValue('plottorquerratio',0.5);
p.parse(varargin{:});

unitlength = p.Results.unitlength;
magratio   = p.Results.magratio;
magn       = p.Results.magn;
displ      = p.Results.displ;
magangle   = p.Results.magangle;
gapratio   = p.Results.gapratio;
dispoffset = p.Results.dispoffset;
leverratio = p.Results.leverratio;
rotation   = p.Results.rotation;
plot_bool  = p.Results.plot;
plotsize   = p.Results.plotsize;
plot_extras_bool  = p.Results.plotextras;
plotvecscale      = p.Results.plotvecscale;
plot_torque_scale      = p.Results.plottorquescale;
plot_torque_arc      = p.Results.plottorquearc;
plot_torque_radius_ratio      = p.Results.plottorquerratio;

%% Theory

% Rotation matrices
%
% R21 rotates a vector in frame 1 into frame 2.
% (Consider it a CCW rotation.) E.g., if frame 2 is +45° from frame 1,
% a vector at 60° in frame 1 is at 15° in frame 2.
%
% R12 is obviously the opposite, but it can also be considered to be simply
% a regular rotation matrix; a vector at 15° applied to R12(45°) will then
% be at 60°.

R21 = @(T) [cosd(T) sind(T) 0;-sind(T) cosd(T) 0; 0 0 1];
R12 = @(T) R21(-T);
R31 = @(T) R21(180-T);
R13 = @(T) R21(T-180);

% Magnet setup:
% In the local coordinate system of the magnets, their magnetisation is
% in the x direction.

mag = @(a,b,magdir) ...
  struct('dim',[a b b],'magn',magn,'magdir',[magdir 0 0]);

magl = @(a,b,magdir,ll) ...
  struct('dim',[a b b],'magn',magn,'magdir',[magdir 0 0],'lever',ll);

% Displacements of the magnets due to rotation in the global coordinate system

r11 = @(l,r)  R12(r)*[-l; 0; 0] - [-l; 0; 0];
r21 = @(l,r)  R12(r)*[ l; 0; 0] - [ l; 0; 0];

% Displacements between magnets in the magnets' coordinate system

D1 = @(Dxyz,T,a,d,l,r)  [ a; 0; 0 ] + R21(T)*(Dxyz+[ d; 0; 0]+r11(l,r)+dispoffset);
D2 = @(Dxyz,T,a,d,l,r)  [ a; 0; 0 ] + R31(T)*(Dxyz+[-d; 0; 0]+r21(l,r)+dispoffset);

% Displacements between floating magnets and centre of rotation (leverarms)
% in the coordinate system of the magnets
% (Here R12 is being used simply as a CCW rotation matrix.)

L1  = @(T,l,r) R21(T)*( -R12(r)*[-l; 0; 0] );
L2  = @(T,l,r) R31(T)*( -R12(r)*[ l; 0; 0] );

% Forces

% f1 = @(Dxyz,T,a,b,d,l,r)  magnetforces(mag(a,b,1),mag(a,b,-1),D1(Dxyz,T,a,d,l,r));
% f2 = @(Dxyz,T,a,b,d,l,r)  magnetforces(mag(a,b,1),mag(a,b,-1),D2(Dxyz,T,a,d,l,r));

rotate_z_to_x = @(vec) [  vec(3,:);  vec(2,:); -vec(1,:) ] ; % Ry( 90)
rotate_x_to_z = @(vec) [ -vec(3,:);  vec(2,:);  vec(1,:) ] ; % Ry(-90)

forces_calc_x_x = @(a,b,c,A,B,C,d,J1,J2) ...
  rotate_z_to_x( forces_calc_z_z(c,b,a,C,B,A, rotate_x_to_z(d), J1,J2) );

f1 = @(Dxyz,T,a,b,d,l,r) ...
  forces_calc_x_x(a/2,b/2,b/2,a/2,b/2,b/2,D1(Dxyz,T,a,d,l,r),magn,-magn);
f2 = @(Dxyz,T,a,b,d,l,r) ...
  forces_calc_x_x(a/2,b/2,b/2,a/2,b/2,b/2,D2(Dxyz,T,a,d,l,r),magn,-magn);

% t1 = @(Dxyz,T,a,b,d,l,r)  magnetforces(mag(a,b,1),magl(a,b,-1,L1(T,l,r)),D1(Dxyz,T,a,d,l,r),'torque');
% t2 = @(Dxyz,T,a,b,d,l,r)  magnetforces(mag(a,b,1),magl(a,b,-1,L2(T,l,r)),D2(Dxyz,T,a,d,l,r),'torque');

t1 = @(Dxyz,T,a,b,d,l,r) ...
...%  torques_calc_x_x(a/2,b/2,b/2,a/2,b/2,b/2,...
  ztorque_calc_x_x(a/2,b/2,b/2,a/2,b/2,b/2,...
    D1(Dxyz,T,a,d,l,r),L1(T,l,r),magn,-magn);

t2 = @(Dxyz,T,a,b,d,l,r) ...
...%  torques_calc_x_x(a/2,b/2,b/2,a/2,b/2,b/2,...
  ztorque_calc_x_x(a/2,b/2,b/2,a/2,b/2,b/2,...
    D2(Dxyz,T,a,d,l,r),L2(T,l,r),magn,-magn);

%% Setup

Ns = length(unitlength);
Ng = length(gapratio);
Nm = length(magratio);
NT = length(magangle);
NL = length(leverratio);
NR = length(rotation);

points = size(displ,2);

% initialise variables
forces  = repmat(NaN,[3 Ns NT Ng Nm NL NR points]);
torques = repmat(NaN,[3 Ns NT Ng Nm NL NR points]);

%% Forces calc

if plot_bool
  plot_h = gcf;
  set(plot_h,'defaultlinelinewidth',0.7,'defaulttextinterpreter','none');
  pvec = @(t) repmat(t,[1 5]);
end

for uu = 1:Ns
  volume = unitlength(uu)^3;
  for mm = 1:Nm
    a = (volume*magratio(mm)^2)^(1/3);
    b = (volume/magratio(mm))^(1/3);
    for gg = 1:Ng
      d = gapratio(gg)*unitlength(uu);
      for tt = 1:NT
        t = magangle(tt);
        for ll = 1:NL
          l = leverratio(ll)*unitlength(uu);
          for rr = 1:NR
            r = rotation(rr);
            for yy = 1:points
              
              f11 = R12(t)*f1(displ(:,yy),t,a,b,d,l,r);
              f21 = R13(t)*f2(displ(:,yy),t,a,b,d,l,r);
              forces(:,uu,tt,gg,mm,ll,rr,yy)  = f11 + f21;
              
              if nargout == 2
                t11 = t1(displ(:,yy),t,a,b,d,l,r);
                t21 = t2(displ(:,yy),t,a,b,d,l,r);
                torques(:,uu,tt,gg,mm,ll,rr,yy) = t11 + t21;
              end

              if plot_bool
                
                cla
                hold on

                rotdisp = repmat(r11(l,r),[1 5]);
                rotdisp(3,:) = [];
                
                rot = @(t) [cosd(t) -sind(t);sind(t) cosd(t)];
                mag1 = rot(t)*[[-a;-b] [a;-b] [a;b] [-a;b] [-a;-b]]/2;
                mag2 = pvec(2*[d+l;0]+[2*a*cosd(t);0]) + ...
                       rot(-t)*([[-a;-b] [a;-b] [a;b] [-a;b] [-a;-b]]/2);
                     
                mag3u = mag1 + pvec(displ([1 2],yy)+[ d;0]+[ a*cosd(t);a*sind(t)]);
                mag3 = mag3u + rotdisp;
                mag3c = mean(mag3(:,1:4),2);
                mag3r = rot(r)*(mag3 - pvec(mag3c)) + pvec(mag3c);
                
                mag4u = mag2 + pvec(displ([1 2],yy)+[-d;0]+[-a*cosd(t);a*sind(t)]);
                mag4 = mag4u - rotdisp;
                mag4c = mean(mag4(:,1:4),2);
                mag4r = rot(r)*(mag4 - pvec(mag4c)) + pvec(mag4c);
                
                rotc = [mean([mag3c(1), mag4c(1)]),mean([mag3c(2), mag4c(2)])]';

                plot([mag3c(1) mag4c(1)],[mag3c(2) mag4c(2)],'k')

                plot(mag1(1,:),mag1(2,:),'k')
                plot(mag2(1,:),mag2(2,:),'k')
                
                if plot_extras_bool
                  plot([mag1(1,4) mag2(1,3)],[mag1(2,1) mag2(2,2)],'k--')

                  plot(mag3(1,:),mag3(2,:),'k')
                  plot(mag4(1,:),mag4(2,:),'k')
                  if r ~= 0
                    grey = 0.7*[1 1 1];
                    plot(mag3r(1,:),mag3r(2,:),'--','color',grey)
                    plot(mag4r(1,:),mag4r(2,:),'--','color',grey)
                  end
                else
                  plot(mag3r(1,:),mag3r(2,:),'k')
                  plot(mag4r(1,:),mag4r(2,:),'k')
                end
                
                if plot_extras_bool
                  
                  text(mag3c(1)+[0.25 0.75]*(mag4c(1)-mag3c(1)),mag3c(2)+[0.25 0.75]*(mag4c(2)-mag3c(2)),...
                    '$\mbqlever$',...
                    'VerticalAlignment','bottom')
                  
                  plot_angle_label( rotc, r, a, '$\mbqrotz$')
                  plot_angle_label( mag1(:,1), t, a, '$\theta$')
                end
                
                axis tight
                set(gca,'xtick',[],'ytick',[])
                if plotsize ~= 0
                  xl = xlim;
                  xlim(mean(xl)+[-plotsize(1) plotsize(1)]/2)
                  axis equal
                  yl = ylim;
                  ylim([-0.1*plotsize(1) yl(2)-yl(1)]);
                else
                  axis equal
                end
                axis off
                
                plot_vec(mag3c,f11,'facecolor',[0 0 1],'edgecolor',[0 0 1])
                plot_vec(mag4c,f21,'facecolor',[0 0 1],'edgecolor',[0 0 1])
                plot_vec(rotc,f11+f21,'facecolor',[1 0 0],'edgecolor',[1 0 0],'linewidth',1.4)
                
                plot(rotc(1),rotc(2),'k.','markersize',12)
                
                this_torque = t11(3,:)+t21(3,:);
                if abs(this_torque) > 1e-6
                  
                  torque_scale = 10;
                  if plot_torque_scale == 0
                    torque_rot = plot_torque_arc;
                    disp(['Torque scale of ',num2str(abs(pi/4/this_torque))])
                  else
                    torque_rot = abs(this_torque)*plot_torque_scale;
                  end
                  
                  torque_radius = l*plot_torque_radius_ratio;
                  
                  torque_points = 50;
                  if this_torque > 0
                    torque_angles = linspace(-torque_rot+r*pi/180,torque_rot+r*pi/180,torque_points);
                  else
                    torque_angles = fliplr(linspace(pi-torque_rot+r*pi/180,pi+torque_rot+r*pi/180,torque_points));                    
                  end
                  
                  torque_arc = [rotc(1)+torque_radius*cos(torque_angles);
                                rotc(2)+torque_radius*sin(torque_angles)];
                  
                  torque_arrow_colour = [0 0.8 0.2];
                  plot(torque_arc(1,:),torque_arc(2,:),...
                    'color',torque_arrow_colour,'linewidth',1.4)
                  
                  arrow(torque_arc(:,end-1),torque_arc(:,end),'Length',12,'facecolor',torque_arrow_colour,'edgecolor',torque_arrow_colour)
                  
                end

                drawnow
                
              end
              
            end
          end
        end
      end
    end
  end
end

forces = squeeze(forces);
if nargout == 2
  torques = squeeze(torques);
end

  function plot_vec(r,v,varargin)
    
    arrow(r(:),r(:)+plotvecscale*v(1:2),'length',10,varargin{:})
    
  end

  function plot_angle_label(cp,t,l,s)
    plot(cp(1)+[0,l],cp(2)+[0,0],'k--')
    ap = cp+rot(t/2)*[2*l/3; 0];
    lp = ap+[0; -l/2];
    plot([ap(1) lp(1)],[ap(2) lp(2)],'color',0.6*[1 1 1])
    % plot(ap(1),ap(2),'k.')
    text(lp(1),lp(2),['\,',s])
  end

end



function [force] = forces_calc_z_z(a,b,c,A,B,C,displ,J1,J2)
%% forces between cuboid magnets with parallel z-direction magnetisation
%
% The theory of Akoun & Yonnet (1984)
%
% Inputs:     (a,b,c) - the HALF dimensions of the
%                       fixed magnet
%             (A,B,C) - the HALF dimensions of the
%                       floating magnet
%  (alpha,beta,gamma) - distance between magnet centres
%              J[,J2] - magnetisations of the
%                       magnet(s) in the z-direction
%
% Outputs: (Fx,Fy,Fz) - Forces on the second magnet
%
%
% You probably want to call
%   warning off MATLAB:divideByZero
%   warning off MATLAB:log:logOfZero

alpha = displ(1);
beta = displ(2);
gamma = displ(3);

[index_h, index_j, index_k, index_l, index_p, index_q] = ndgrid([0 1]);
index_sum = (-1).^(index_h+index_j+index_k+index_l+index_p+index_q);

% (Using this method is actually LESS efficient than using six for
% loops for h..q over [0 1], but it looks a bit nicer, huh?)

u = alpha + A*(-1).^index_j - a*(-1).^index_h;
v = beta  + B*(-1).^index_l - b*(-1).^index_k;
w = gamma + C*(-1).^index_q - c*(-1).^index_p;
r = sqrt(u.^2+v.^2+w.^2);

magconst = J1*J2/(4*pi*(4*pi*1e-7));

f_x = ...
  + 0.5*(v.^2-w.^2).*log(r-u) ...
  + u.*v.*log(r-v) ...
  + v.*w.*atan(u.*v./r./w) ...
  + 0.5*r.*u;
f_y = ...
  + 0.5*(u.^2-w.^2).*log(r-v) ...
  + u.*v.*log(r-u) ...
  + u.*w.*atan(u.*v./r./w)...
  + 0.5*r.*v;
f_z = ...
  - u.*w.*log(r-u) ...
  - v.*w.*log(r-v) ...
  + u.*v.*atan(u.*v./r./w) ...
  - r.*w;

fx = index_sum.*f_x;
fy = index_sum.*f_y;
fz = index_sum.*f_z;

Fx = magconst*sum(fx(:));
Fy = magconst*sum(fy(:));
Fz = magconst*sum(fz(:));

force = [Fx; Fy; Fz];

end



function torques = ztorque_calc_x_x(a1,b1,c1,a2,b2,c2,offset,lever,J1,J2)
%% forces between cuboid magnets with parallel z-direction magnetisation
%
% coordinates are transformed from x-magnetisation to z and then back again

a = -offset(3,:);
b =  offset(2,:);
c =  offset(1,:);

d = a-lever(3,:);
e = b+lever(2,:);
f = c+lever(1,:);

Tx=zeros([1 size(offset,2)]);
Ty=Tx;
Tz=Tx;

for ii=[0,1]
  for jj=[0,1]
    for kk=[0,1]
      for ll=[0,1]
        for mm=[0,1]
          for nn=[0,1]
                        
            Cu=(-1)^ii.*c1-d;
            Cv=(-1)^kk.*b1-e;
            Cw=(-1)^mm.*a1-f;
            
            u=a-(-1)^ii.*c1+(-1)^jj.*c2;
            v=b-(-1)^kk.*b1+(-1)^ll.*b2;
            w=c-(-1)^mm.*a1+(-1)^nn.*a2;
            
            s=sqrt(u.^2+v.^2+w.^2);

            Ex=(1/8).*(...
              -2.*Cw.*(-4.*v.*u+s.^2+2.*v.*s)-...
              w.*(-8.*v.*u+s.^2+8.*Cv.*s+6.*v.*s)+...
              2.*(2.*Cw+w).*(u.^2+w.^2).*log(v+s)+...
              4.*(...
                2.*Cv.*u.*w.*acoth(u./s) + ...
                w.*(v.^2+2.*Cv.*v-w.*(2.*Cw+w)).*acoth(v./s) - ...
                u.*(...
                  2*w.*(Cw+w).*atan(v./w) + ...
                  2*v.*(Cw+w).*log(s-u) + ...
                  (w.^2+2.*Cw.*w-v.*(2.*Cv+v)).*atan( u.*v./(w.*s) ) ...
                )...
              )...
            );

            Tx=Tx+(-1)^(ii+jj+kk+ll+mm+nn)*Ex;
            
          end
        end
      end
    end
  end
end

torques = real([0; 0; -Tx].*J1*J2/(16*pi^2*1e-7));
                        
end