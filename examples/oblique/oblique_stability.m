%% Stability in quasi-3DOF for the oblique magnet spring
%
% For now a hodge-podge of examples.

%%

close all



%% schematic for paper

willfig('angle-schem','huge'); clf

[f t] = oblique_forces3(...
  'displ',[0;1;0]*0.01,...
  'plot',true,'plotextras',true,...
  'plotvecscale',0.005,...
  'plottorquerratio',0.3,...
  'magangle',30,...
  'gapratio',0.5,...
  'rotation',15 ...
);

matlabfrag('fig/mbq-angle-schem')



%% vertical displacement

displ_range = linspace(0.005,0.02,50);

willfig('plotanim2');
figuresize(6,20,'centimeters')

[f t] = oblique_forces3(...
  'displ',[0;1;0]*displ_range,...
  'plot',true,...
  'plotvecscale',0.005,...
  'magangle',30,...
  'gapratio',0,...
  'rotation',0....
);

willfig('y displ; y force2')
plot(1000*displ_range,squeeze(f(2,:)))
xlabel('Y DISPL, mm'); ylabel('Y FORCE, N')

%%

clc
rot = 0:0.5:10;
vert = 0.01;
lever = [2 3];
gaps = [0.25 0.5 1];

[f t] = oblique_forces3(...
  'displ',[0;1;0]*vert,...
  'magangle',30,...
  'leverratio',lever,...
  'gapratio',gaps,...
  'rotation',rot...
);

tx = t(1,:); ty = t(2,:); tz = 1000*squeeze(t(3,:,:,:));
assert( all(tx(:)<1e-10) && all(ty(:)<1e-10) ,...
  'out-of-plane torques should equal zero')

%%

ta = [pi/2 pi/4 pi/4];
for ii = 1:3
willfig(['show-rot-',num2str(ii)]); clf; hold on
[f t] = oblique_forces3(...
  'displ',[0;vert;0],...
  'leverratio',lever(1),...
  'plot',true,...
  'plotsize',0.1,...
  'plottorquearc',ta(ii),...
  'magangle',30,...
  'gapratio',gaps(ii),...
  'rotation',rot(round(end/2))...
);
matlabfrag(['fig/mbq-rot-ex-diag-',num2str(ii)])
end

%%

willfig('show-rot','huge'); clf

for jj = 1:2
for ii = 1:3
subplot(2,3,(jj-1)*3+ii); hold on
[f t] = oblique_forces3(...
  'displ',[0;vert;0],...
  'leverratio',lever(jj),...
  'plot',true,...
  'magangle',30,...
  'gapratio',gaps(ii),...
  'rotation',rot(round(end/2))...
);
end
end





%% Paper figure (rotation example)

ltext = @(x,y,s) text(x,y,num2str(s),'horizontalalignment','c','userdata',...
  ['matlabfrag:\fboxsep=1pt\colorbox{white}{$',num2str(s),'$}']);

willfig('rotation; z torque 2','small'); clf; hold on
ind = 17*ones(size(gaps));
for gg = 2:length(gaps)
  plot(rot,squeeze(tz(gg,1,:)))
  ltext(rot(ind(gg)),tz(gg,1,ind(gg)),num2str(gaps(gg)))
end
for gg = 1
  plot(rot,squeeze(tz(gg,1,:)))
  ltext(rot(ind(gg)),tz(gg,1,ind(gg)),['\mbqoffset/\mbqunit=',num2str(gaps(gg))])
end
colourplot
draworigin

xlabel('$z$ rotation $\mbqrotz$, deg.')
moreticks
ylabel('$z$ torque $\mbqptorque$, \si{mN.m}')

matlabfrag('fig/mbq-rot-ex')

%%

willfig('rotation; z torque - 2'); clf; hold on

ind = 17*ones(size(gaps));
for gg = 1:length(gaps)
  plot(rot,squeeze(tz(gg,2,:)),'--')
  ltext(rot(ind(gg)),tz(gg,2,ind(gg)),num2str(gaps(gg)))
end
colourplot

xlabel('$z$ rotation, deg.')
moreticks
ylabel('$z$ torque, \si{mN.m}')

matlabfrag('fig/mbq-rot-ex-2')




%% Dynamic analysis (only vertical)
%
% Just to show the dynamic simulation works and gives reasonable results.

tmax = 10;

[T1, X1] = oblique_dynamics(...
  'mass',3,...
  'dampingratio',[0.2 0.2 0.5],...
  'tmax',4,...
  'RelTol',1e-5,'AbsTol',1e-5,...
  'perturb',[0 0.001 0]...
  );

willfig('dyn-y'); clf; hold on
plot(T1,1000*X1(:,1))
xlabel('Time, s'); ylabel('Vertical displ., mm')

%% Dynamic analysis (vertical & horizontal)
%
% Specifically choose parameters such that the horizontal stiffness is
% positive.

[T2, X2, param] = oblique_dynamics(...
  'mass',3,...
  'dampingratio',[0.2 0.2 0.5],...
  'tmax',12,...
  'RelTol',1e-5,'AbsTol',1e-5,...
  'perturb',[0.001 0.001 0]...
  );

y0 = param.y0;

willfig('dyn-xy'); clf; hold on
plot(T2,1000*X2(:,1),T2,1000*(X2(:,3)-y0))
legend('x','y')

willfig('map-xy'); clf;
plot(1000*X2(:,1),1000*(X2(:,3)-y0))
xlabel('Horiz., mm')
ylabel('Vert., mm')
axis equal

%% Vertical and rotational

[T3, X3, param] = oblique_dynamics(...
  'mass',3,...
  'dampingratio',[0.2 0.2 0.5],...
  'tmax',12,...
  'RelTol',1e-6,'AbsTol',1e-6,...
  'perturb',[0 0.001 1]...
  );

%%

willfig('dyn-yr'); clf; hold on
plot(T3,1000*(X3(:,1)-param.y0),T3,X3(:,3))
legend('y','r')





%% yann check


facesize = .0254 ; %mm
thickness = .0127; %mm
m = (facesize^2*thickness)^(1/3);

figure
[f t] = oblique_forces3(...
  'magn',1.3,...
  'unitlength',m,...
  'magratio',0.5,...
  'displ',[0;0.01702;0],...
  'plot',true,...
  'magangle',85,...
  'gapratio',0.25,...
  'rotation',0....
);
axis tight
2*f(2)



