%% This script compares analytical results with numerical results
clear all; close all; clc;

hold on;
    
mu = 1.05;                                  % Permeability (N/A^2)
mesh = 0.005;                               % Meshsize numerical solver

% Dimensions [radius height (thickness)]
magnet_fixed.dim = [0.06 0.03];             % Magnet 1
magnet_float.dim = [0.06 0.03];             % Magnet 2

% Add or remove comment sign in order to use ring magnets or not
magnet_fixed2.dim = [.03 .03];              % Used for upper ring magnet
magnet_float2.dim = [.03 .03];              % Used for lower ring magnet

% Check for existence of magnet dimensions
A = length(magnet_fixed.dim);
B = exist('magnet_fixed2','var');
C = exist('magnet_float2','var');

% Determine magnet configuration
if A==2 && B==0 && C==0
    magnettype='axisymmetric';              % Two axisymmetric (cylindrical) magnets
    disp('Matlab: axisymmetric magnets');
end
if A==2 && B==1 && C==0
    magnettype='axiring';                   % An axisymmetric (cylindrical) magnet and a ring magnet
    disp('Matlab: cylinder and ring magnet');
end
if A==2 && B==1 && C==1
    magnettype='rings';                     % Two ring magnets
    disp('Matlab: ring magnets');
end
if (length(magnet_fixed.dim))==3
    magnettype='planar';                    % Two planar (cubic) magnets
    disp('Matlab: planar magnets');
end

% Magnetisation, change for repulsion or attraction
magnet_fixed.magn = 1.2;                    % Upper magnet
magnet_float.magn = -1.2;                   % Lower magnet
magnet_fixed2.magn = -magnet_fixed.magn;    % Upper, inner magnet used for superposition for ring magnet
magnet_float2.magn = -magnet_float.magn;    % Lower, inner magnet used for superposition for ring magnet

% Displacement
displ_max = 0.05;
N = 10;                                     % Number of calculation steps
            switch magnettype
                case 'axisymmetric'
                    offset = repmat([0;0;0.03],[1 N*10]);
                    displ = linspace(0,displ_max,N*10);
                    displ_range = offset+[0;0;1]*displ;
                    % Calculate forces
                    forces = magnetforces(magnet_fixed,magnet_float,displ_range);
                    % Plot output
                    plot((displ),forces(3,:),'k')
                    title('Force between two cylindrical magnets');
                    xlabel('Displacement, m');
                    ylabel('Force, N');
                    %ansys(N,displ,mu,mesh,magnet_fixed,magnet_float);
                case 'planar'
                    magnet_fixed.magdir = [0 1 0];
                    magnet_float.magdir = [0 1 0];
                    offset = repmat([0;0.03;0],[1 N*10]);
                    displ = linspace(eps,displ_max,N*10);
                    displ_range = offset+[0;1;0]*displ;
                    % Calculate forces
                    forces = magnetforces(magnet_fixed,magnet_float,displ_range);
                    % Plot output
                    plot((displ),forces(2,:),'k')
                    title('Force between two planar magnets');
                    xlabel('Displacement, m');
                    ylabel('Force, N');
                    %ansys(N,displ,mu,mesh,magnet_fixed,magnet_float);
                case 'axiring'
                    offset = repmat([0;0;0.03],[1 N*10]);
                    displ = linspace(eps,displ_max,N*10);
                    displ_range = offset+[0;0;1]*displ; 
                    % Calculate forces
                    forces = magnetforces(magnet_fixed,magnet_float,displ_range)+magnetforces(magnet_fixed2,magnet_float,displ_range);
                    % Plot output
                    plot((displ),forces(3,:),'k')
                    title('Force between cylinder and ring magnet');
                    xlabel('Displacement, m');
                    ylabel('Force, N');
                    %ansys(N,displ,mu,mesh,magnet_fixed,magnet_float,magnet_fixed2);
                case 'rings'
                    offset = repmat([0;0;0.03],[1 N*10]);
                    displ = linspace(0,displ_max,N*10);
                    displ_range = offset+[0;0;1]*displ; 
                    % Calculate forces
                    forces = magnetforces(magnet_fixed,magnet_float,displ_range)+magnetforces(magnet_fixed2,magnet_float,displ_range)+magnetforces(magnet_fixed,magnet_float2,displ_range)+magnetforces(magnet_fixed2,magnet_float2,displ_range);
                    plot((displ),forces(3,:),'k')
                    title('Force between ring magnets');
                    xlabel('Displacement, m');
                    ylabel('Force, N');
                    %ansys(N,displ,mu,mesh,magnet_fixed,magnet_float,magnet_fixed2,magnet_float2);
                otherwise
            end

% Load ANSYS results    
load output.txt;
Q = output;
clear output;

% Select correct Ansys values
Qmx = Q(3:4:end);                               % Maxwell force
Qvw = Q(4:4:end);                               % Virtual work force
ansys_displ = linspace(0,displ_max,numel(Qvw));

% Plot ANSYS results
switch magnettype
    case 'axisymmetric'
        plot(ansys_displ,Qmx,'g+:');
        plot(ansys_displ,Qvw,'r+:');
        legend('Matlab','Maxwell','Virtual work');
    case 'planar'
        plot(ansys_displ,Qmx*magnet_fixed.dim(1,3),'g+:');
        plot(ansys_displ,Qvw*magnet_fixed.dim(1,3),'r+:');
        legend('Matlab','Maxwell','Virtual work');
    case 'axiring'
        plot(ansys_displ,Qmx,'g+:');
        plot(ansys_displ,Qvw,'r+:');
        legend('Matlab','Maxwell','Virtual work');
    case 'rings'
        plot(ansys_displ,Qmx,'g+:');
        plot(ansys_displ,Qvw,'r+:');
        legend('Matlab','Maxwell','Virtual work');
    otherwise
end