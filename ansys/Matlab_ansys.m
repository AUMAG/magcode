%%This script compares analytical Matlab with numerical ANSYS results

clear all;close all;clc;
!del output.txt
!del total.txt
% Magnet 1 is upper magnet
a=0.3;          %air gap width
b=0.6;          %air gap height
c=0.05;         %radius magnet 2
e=0.14;         %distance between magnet 2 and bottom of airgap
d=0.03;         %height of magnet 2
%x=0.03;        %distance between 2 magnets, not needed in this script
w=0.03;         %height of magnet 1
u=0.035;        %raduis magnet 1
B1=1.2;         %magnetisation magnet 1 (T)
B2=1.2;         %magnetisation magnet 2(T)
mu=1.05;        %permeability (N/A^2)
mu0=4*pi*10^-7; %permeability of space (N/A^2)
H1=B1/(mu*mu0); %coercivity magnet 1 (A/m)
H2=B2/(mu*mu0); %coercivity magnet 2 (A/m)
mesh=0.005;     %meshsize numerical solver

%% analytical solution
% Dimensions [radius height]
magnet_upper.dim=[u w];
magnet_lower.dim=[c d];

% Magnetisation
magnet_upper.magn=B1;
magnet_lower.magn=-B2;                          %IMPORTANT, change for repulsion(-) or attraction

% Displacement
N=10;              %number of calculation steps
offset=repmat([0;0;0.03],[1 N]);
displ=linspace(0,.4,N);
displ_range=offset+[0;0;1]*displ;

% Calculate forces
forces=magnetforces(magnet_upper,magnet_lower,displ_range);

% Plot output
plot((displ),forces(3,:))
title('Force between two cylindrical magnets');
xlabel('Displacement, m');
ylabel('Force, N');

% Numerical solution
for i=1:N-1
    x=.4/N*i;           % Variable distance between magnets
A_str=num2str(a);
B_str=num2str(b);
C_str=num2str(c);
D_str=num2str(d);
E_str=num2str(e);
ED_str=num2str(e+d);
X_str=num2str(x);
W_str=num2str(w);
EDX_str=num2str(e+d+x);
EDXW_str=num2str(e+d+x+w);
U_str=num2str(u);
Mu_str=num2str(mu);
B1_str=num2str(B1);
B2_str=num2str(B2);
H1_str=num2str(H1);
H2_str=num2str(-H2);                % IMPORTANT, change for repulsion(-) or attraction
MESH_str=num2str(mesh);
semi='''';


% Print ansys batch code to text file total.txt
fid1=fopen('total.txt','wt');
fprintf(fid1,'/PREP7\n');
fprintf(fid1,'/TITLE,numerical solution\n');
fprintf(fid1,'l=ARG1\n');
fprintf(fid1,'KEYW,PR_SET,1\n');
fprintf(fid1,'KEYW,PR_ELMAG,1\n');
fprintf(fid1,'KEYW,MAGNOD,1\n');
fprintf(fid1,'ET,1,PLANE53\n');
fprintf(fid1,'KEYOPT,1,3,1\n');
fprintf(fid1,'MP,MURX,1,1\n');
fprintf(fid1,'MP,MURX,2,');
fprintf(fid1,'%s',B2_str);
fprintf(fid1,'\n');
fprintf(fid1,'MP,MGYY,2,');
fprintf(fid1,'%s',H1_str);
fprintf(fid1,'\n');
fprintf(fid1,'MP,MURX,3,');
fprintf(fid1,'%s',B1_str);
fprintf(fid1,'\n');
fprintf(fid1,'MP,MGYY,3,');
fprintf(fid1,'%s',H2_str);
fprintf(fid1,'\n');
fprintf(fid1,'RECTNG,0,');
fprintf(fid1,'%s',A_str);
fprintf(fid1,',0,');
fprintf(fid1,'%s',B_str);
fprintf(fid1,'\n');
fprintf(fid1,'RECTNG,0,');  %+l?
fprintf(fid1,'%s',C_str);
fprintf(fid1,',');
fprintf(fid1,'%s',E_str);
fprintf(fid1,',');
fprintf(fid1,'%s',ED_str);
fprintf(fid1,'\n');
fprintf(fid1,'RECTNG,0,');
fprintf(fid1,'%s',U_str);
fprintf(fid1,',');
fprintf(fid1,'%s',EDX_str);
fprintf(fid1,',');
fprintf(fid1,'%s',EDXW_str);
fprintf(fid1,'\n');

fprintf(fid1,'ASEL,S,AREA,,2\n');
fprintf(fid1,'AATT,2,,1\n');
fprintf(fid1,'ASEL,S,AREA,,3\n');
fprintf(fid1,'AATT,3,,1\n');

fprintf(fid1,'ASEL,ALL\n');
fprintf(fid1,'ASBA,1,ALL,,,KEEP\n');
fprintf(fid1,'ASEL,U,AREA,,2,3\n');
fprintf(fid1,'AATT,1,,1\n');
%boundary conditions
fprintf(fid1,'ALLSEL\n');
fprintf(fid1,'LSEL,S,EXT\n');
fprintf(fid1,'DL,ALL,,AZ,0\n');
%mesh
fprintf(fid1,'ALLSES\n');
fprintf(fid1,'AESIZE,ALL,');
fprintf(fid1,'%s',MESH_str);
fprintf(fid1,'\n');
fprintf(fid1,'AMESH,ALL\n');
%apply force flags
fprintf(fid1,'ESEL,S,MAT,,3\n');
fprintf(fid1,'CM,mag1,ELEM\n');
fprintf(fid1,'FMAGBC,');
fprintf(fid1,'%s',semi);
fprintf(fid1,'mag1');
fprintf(fid1,'%s',semi);
fprintf(fid1,'\n');
fprintf(fid1,'PRIM\n');

fprintf(fid1,'FINISH\n');
fprintf(fid1,'ALLSEL\n');
fprintf(fid1,'/SOL\n');
fprintf(fid1,'MAGSOLV,2\n');
fprintf(fid1,'FINISH\n');

fprintf(fid1,'/POST1\n');
fprintf(fid1,'PLF2D\n');
fprintf(fid1,'fmagsum,');
fprintf(fid1,'%s',semi);
fprintf(fid1,'mag1');
fprintf(fid1,'%s',semi);
fprintf(fid1,'\n');

fprintf(fid1,'ESEL,S,MAT,,3\n');
fprintf(fid1,'NSEL,S,EXT\n');
fprintf(fid1,'ESLN\n');
fprintf(fid1,'ESEL,U,MAT,,3\n');
fprintf(fid1,'SSUM\n');

fprintf(fid1,'etab,fmgx1,fmag,x\n');
fprintf(fid1,'etab,fmgy1,fmag,y\n');
fprintf(fid1,'etab.fmgz1,fmag,z\n');

fprintf(fid1,'SSUM\n');

fprintf(fid1,'*dim,ff,ARRAY,4\n');
fprintf(fid1,'*get,ff(1),ssum,fmx_x\n');
fprintf(fid1,'*get,ff(2),ssum,fvw_x\n');
fprintf(fid1,'*get,ff(3),ssum,fmx_y\n');
fprintf(fid1,'*get,ff(4),ssum,fvw_y\n');

fprintf(fid1,'/PAGE,5000,5000,5000,5000\n');
fprintf(fid1,'/HEADER,OFF,OFF,OFF,OFF,OFF,OFF\n');

fprintf(fid1,'/OUTPUT,');
fprintf(fid1,'%s',semi);
fprintf(fid1,'output'); %name can be changed
fprintf(fid1,'%s',semi);
fprintf(fid1,',');
fprintf(fid1,'%s',semi);
fprintf(fid1,'txt');
fprintf(fid1,'%s',semi);
%fprintf(fid1,',,\n');
fprintf(fid1,',,APPEND\n'); %changes with previous line to append to output

fprintf(fid1,'*VWRITE,ff(1),ff(2),ff(3),ff(4)\n');
fprintf(fid1,'(E10.3)\n');
fprintf(fid1,'/OUTPUT\n');
fprintf(fid1,'/EOF');
fclose(fid1);


% Run the ANSYS batch file
!"C:\Program Files\ANSYS Inc\v140\ansys\bin\intel\ANSYS140" -b -i C:\Code\matlab\total.txt -o C:\Code\\matlab\result.out
    
    i=i+1;
end

% Load ANSYS results    
load output.txt;
Q=output;
clear output;


% Plot ANSYS results
for k=1:N-1
    hold on;
    plot(displ(k),Q((k*4)-1),'g+:');
    plot(displ(k),Q((k*4)),'b+:');
    title('numerical results');
    xlabel('displacement, m');
    ylabel('Force, N');
    k=k+1;
end