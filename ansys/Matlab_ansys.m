%% This script compares analytical results with numerical results
clear all; close all; clc;

figure(77)
hold on;
    
mu=1.05;                        %permeability (N/A^2)
mesh=0.005;                     %meshsize numerical solver

% Dimensions [radius height (thickness)]
magnet_fixed.dim=[0.035 0.03 .1];  % magnet 1
magnet_float.dim=[0.05  0.03 .1];  % magnet 2

if (length(magnet_fixed.dim))==2
    magnettype='axisymmetric';
     disp('Matlab: axisymmetric magnets');
end
if (length(magnet_fixed.dim))==3
    magnettype='planar';
    disp('Matlab: planar magnets');
end


% Magnetisation
magnet_fixed.magn=1.2;
magnet_float.magn=-1.2;         %IMPORTANT, change for repulsion(-) or attraction

% Displacement
displ_max = 0.05;
N=10;                           %number of calculation steps
            switch magnettype
                case 'axisymmetric'
                    offset=repmat([0;0;0.03],[1 N]);
                    displ=linspace(0,.3,N);
                    displ_range=offset+[0;0;1]*displ; 
                case 'planar'
                    magnet_fixed.magdir=[0 1 0];
                    magnet_float.magdir=[0 1 0];
                    offset=repmat([0;0.03;0],[1 N]);
                    displ=linspace(eps,displ_max,N);
                    displ_range=offset+[0;1;0]*displ;
                otherwise
            end

% Calculate forces
forces=magnetforces(magnet_fixed,magnet_float,displ_range);

        switch magnettype
            case 'axisymmetric'
                % Plot output
                plot((displ),forces(3,:))
                title('Force between two cylindrical magnets');
                xlabel('Displacement, m');
                ylabel('Force, N');
            case 'planar'
                % Plot output
                plot((displ),forces(2,:))
                title('Force between two planar magnets');
                xlabel('Displacement, m');
                ylabel('Force, N');
            otherwise
        end

% Calculate using Ansys
ansys(N,displ,mu,mesh,magnet_fixed,magnet_float);

% Load ANSYS results    
load output.txt;
Q=output;
clear output;

Qvw = Q(3:4:end);
Qmx = Q(4:4:end);
ansys_displ = linspace(0,displ_max,numel(Qvw));

% Plot ANSYS results
switch magnettype
    case 'axisymmetric'
        plot(ansys_displ,Qvw,'g+:');
        plot(ansys_displ,Qmx,'b+:');
        
    case 'planar'
        plot(ansys_displ,Qvw*magnet_fixed.dim(1,3),'g+:');
        plot(ansys_displ,Qmx*magnet_fixed.dim(1,3),'b+:');
        
    otherwise
end