%% This script compares analytical results with numerical results
clear all; close all; clc;

hold on;
    
mu=1.05;                        %permeability (N/A^2)
mesh=0.005;                     %meshsize numerical solver

% Dimensions [radius height (thickness)]
magnet_fixed.dim=[0.06 0.03];  % magnet 1
magnet_float.dim=[0.04 0.03];  % magnet 2

% Add or remove comment sign in order to use ring magnets
magnet_fixed2.dim=[.03 .03]; % Used for upper ring magnet
%magnet_float2.dim=[.03 .03]; % Used for lower ring magnet

% Check for existence of magnet dimensions
A=length(magnet_fixed.dim);
B=exist('magnet_fixed2','var');
C=exist('magnet_float2','var');

if A==2 && B==0 && C==0
    magnettype='axisymmetric';
    disp('Matlab: axisymmetric magnets');
end
if A==2 && B==1 && C==0
    magnettype='axiring';
    disp('Matlab: cylinder and ring magnet');
end
if A==2 && B==1 && C==1
    magnettype='rings';
    disp('Matlab: ring magnets');
end
if (length(magnet_fixed.dim))==3
    magnettype='planar';
    disp('Matlab: planar magnets');
end


% Magnetisation
magnet_fixed.magn=1.2;
magnet_float.magn=-1.2;                 % IMPORTANT, change for repulsion(-) or attraction
magnet_fixed2.magn=-magnet_fixed.magn;  % Upper, inner magnet used for superposition
magnet_float2.magn=-magnet_float.magn;  % Lower, inner magnet used for superposition

% Displacement
displ_max = 0.05;
N=20;                           % Number of calculation steps
            switch magnettype
                case 'axisymmetric'
                    offset=repmat([0;0;0.03],[1 N*10]);
                    displ=linspace(0,.3,N*10);
                    displ_range=offset+[0;0;1]*displ;
                    % Calculate forces
                    forces=magnetforces(magnet_fixed,magnet_float,displ_range);
                    % Plot output
                    plot((displ),forces(3,:))
                    title('Force between two cylindrical magnets');
                    xlabel('Displacement, m');
                    ylabel('Force, N');
                    %ansys(N,displ,mu,mesh,magnet_fixed,magnet_float);
                case 'planar'
                    magnet_fixed.magdir=[0 1 0];
                    magnet_float.magdir=[0 1 0];
                    offset=repmat([0;0.03;0],[1 N*10]);
                    displ=linspace(eps,displ_max,N*10);
                    displ_range=offset+[0;1;0]*displ;
                    % Calculate forces
                    forces=magnetforces(magnet_fixed,magnet_float,displ_range);
                    % Plot output
                    plot((displ),forces(2,:))
                    title('Force between two planar magnets');
                    xlabel('Displacement, m');
                    ylabel('Force, N');
                    %ansys(N,displ,mu,mesh,magnet_fixed,magnet_float);
                case 'axiring'
                    offset=repmat([0;0;0.03],[1 N*10]);
                    displ=linspace(eps,displ_max,N*10);
                    displ_range=offset+[0;0;1]*displ; 
                    % Calculate forces
                    forces=magnetforces(magnet_fixed,magnet_float,displ_range)+magnetforces(magnet_fixed2,magnet_float,displ_range);
                    % Plot output
                    plot((displ),forces(3,:))
                    title('Force between cylinder and ring magnet');
                    xlabel('Displacement, m');
                    ylabel('Force, N');
                    ansys(N,displ,mu,mesh,magnet_fixed,magnet_float,magnet_fixed2);
                case 'rings'
                    offset=repmat([0;0;0.03],[1 N*10]);
                    displ=linspace(0,.3,N*10);
                    displ_range=offset+[0;0;1]*displ; 
                    % Calculate forces
                    forces=magnetforces(magnet_fixed,magnet_float,displ_range)+magnetforces(magnet_fixed2,magnet_float,displ_range);   
                    plot((displ),forces(3,:))
                    title('Force between ring magnets');
                    xlabel('Displacement, m');
                    ylabel('Force, N');
                    %ansys(N,displ,mu,mesh,magnet_fixed,magnet_float,magnet_fixed2,magnet_float2);
                otherwise
            end

% Load ANSYS results    
load output.txt;
Q=output;
clear output;

Qmx = Q(3:4:end);
Qvw = Q(4:4:end);
ansys_displ = linspace(0,displ_max,numel(Qvw));

% Plot ANSYS results
switch magnettype
    case 'axisymmetric'
        plot(ansys_displ,Qmx,'g+:');
        plot(ansys_displ,Qvw,'b+:');
    case 'planar'
        plot(ansys_displ,Qmx*magnet_fixed.dim(1,3),'g+:');
        plot(ansys_displ,Qvw*magnet_fixed.dim(1,3),'b+:');
    case 'axiring'
        plot(ansys_displ,Qmx,'g+:');
        plot(ansys_displ,Qvw,'b+:');
    case 'rings'
        plot(ansys_displ,Qmx,'g+:');
        plot(ansys_displ,Qvw,'b+:');
    otherwise
end