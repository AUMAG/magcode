function []=ansys(N,displ,mu,mesh,magnet_fixed,magnet_float,magnet_fixed2,magnet_float2)
% Numerical solution of two magnets
!del output.txt
!del total.txt

%% Determine magnettype
A=length(magnet_fixed.dim);
B=exist('magnet_fixed2','var');
C=exist('magnet_float2','var');

if A==2 && B==0 && C==0
    magnettype='axisymmetric';
    disp('Ansys: axisymmetric magnets');
end
if A==2 && B==1 && C==0
    magnettype='axiring';

    disp('Ansys: cylinder and ring magnet');
end
if A==2 && B==1 && C==1
    magnettype='rings';
    disp('Ansys: ring magnets');
end
if A==3
    magnettype='planar';
    disp('Ansys: planar magnets');
end

        
%% Write Ansys code
    for i=0:N
    % Upper (fixed) magnet is magnet 1
    % Lower (float) magnet is magnet 2
 
    g=max(displ)-min(displ);
    x=min(displ)+g/N*i;             % Variable distance between magnets
    disp(x);
    mu0=4*pi*10^-7;                 %permeability of space (N/A^2)
    h1=magnet_fixed.magn/(mu*mu0);  %coercivity magnet 1 (A/m)
    h2=magnet_float.magn/(mu*mu0);  %coercivity magnet 2 (A/m)

    b=(magnet_fixed.dim(1,2)+magnet_float.dim(1,2))*10;     % Air gap height, equals total magnet height*10
    a=magnet_float.dim(1,1)*10;                                                  % Air gap width
    c=magnet_float.dim(1,1);        % Lower magnet (outer) radius
    d=magnet_float.dim(1,2);        % Lower magnet height
    u=magnet_fixed.dim(1,1);        % Upper magnet (outer) radius
    e=b/2-(x/2)-d;                  % Distance from bottom of airgap to bottom of lower magnet
    ed=e+d;                         % Distance from bottom of airgap to top of lower magner
    w=magnet_fixed.dim(1,2);        % Upper magnet height
    edx=e+d+x;                      % Distance from bottom of airgap to bottom of upper magnet
    edxw=e+d+x+w;                   % Distance from bottom of airgap to top of upper magnet

    b1=magnet_fixed.magn;
    b2=magnet_float.magn*-1;
    semi='''';
    
    fid1=fopen('total.txt','wt');
    fprintf(fid1,'/PREP7\n');
    fprintf(fid1,'/TITLE,numerical solution\n');
    %fprintf(fid1,'l=ARG1\n');
    fprintf(fid1,'KEYW,PR_SET,1\n');
    fprintf(fid1,'KEYW,PR_ELMAG,1\n');
    fprintf(fid1,'KEYW,MAGNOD,1\n');
    fprintf(fid1,'ET,1,PLANE53\n');
       
        switch magnettype
            case 'axisymmetric'
                fprintf(fid1,'KEYOPT,1,3,1\n');
            case 'planar'
                f=a/2-c/2;          % Distance from left side to left side lower magnet
                h=a/2-u/2;          % Distance from left side to left side upper magnet
            case 'axiring'
                fprintf(fid1,'KEYOPT,1,3,1\n');
                h=magnet_fixed2.dim(1,1); % Distance from left side to inner radius upper magnet
            case 'rings'
                fprintf(fid1,'KEYOPT,1,3,1\n');
                h=magnet_fixed2.dim(1,1); % Distance from left side to inner radius upper magnet
                j=magnet_float2.dim(1,1); % Distance from left side to inner radius lower magnet
        end

    fprintf(fid1,'MP,MURX,1,1\n');
    fprintf(fid1,'MP,MURX,2,');
    fprintf(fid1,'%f',b2);
    fprintf(fid1,'\n');
    fprintf(fid1,'MP,MGYY,2,');
    fprintf(fid1,'%f',h1);
    fprintf(fid1,'\n');
    fprintf(fid1,'MP,MURX,3,');
    fprintf(fid1,'%f',b1);
    fprintf(fid1,'\n');
    fprintf(fid1,'MP,MGYY,3,');            
    fprintf(fid1,'%f',h2);
    fprintf(fid1,'\n');
    fprintf(fid1,'RECTNG,0,');
    fprintf(fid1,'%f',a);
    fprintf(fid1,',0,');
    fprintf(fid1,'%f',b);
    fprintf(fid1,'\n');
    
        switch magnettype
            case 'axisymmetric'
            fprintf(fid1,'RECTNG,0,');
            fprintf(fid1,'%f',c);
            fprintf(fid1,',');
            fprintf(fid1,'%f',e);
            fprintf(fid1,',');
            fprintf(fid1,'%f',ed);
            fprintf(fid1,'\n');
            fprintf(fid1,'RECTNG,0,');
            fprintf(fid1,'%f',u);
            fprintf(fid1,',');
            fprintf(fid1,'%f',edx);
            fprintf(fid1,',');
            fprintf(fid1,'%f',edxw);
            fprintf(fid1,'\n');   
            case 'planar'
            fprintf(fid1,'RECTNG,');
            fprintf(fid1,'%f',f);
            fprintf(fid1,',');
            fprintf(fid1,'%f',f+c);
            fprintf(fid1,',');
            fprintf(fid1,'%f',e);
            fprintf(fid1,',');
            fprintf(fid1,'%f',ed);
            fprintf(fid1,'\n');
            fprintf(fid1,'RECTNG,');
            fprintf(fid1,'%f',h);
            fprintf(fid1,',');
            fprintf(fid1,'%f',h+u);
            fprintf(fid1,',');
            fprintf(fid1,'%f',edx);
            fprintf(fid1,',');
            fprintf(fid1,'%f',edxw);
            fprintf(fid1,'\n');  
            case 'axiring'
            fprintf(fid1,'RECTNG,0,');
            fprintf(fid1,'%f',c);
            fprintf(fid1,',');
            fprintf(fid1,'%f',e);
            fprintf(fid1,',');
            fprintf(fid1,'%f',ed);
            fprintf(fid1,'\n');
            fprintf(fid1,'RECTNG,');
            fprintf(fid1,'%f',h);
            fprintf(fid1,',');
            fprintf(fid1,'%f',u);
            fprintf(fid1,',');
            fprintf(fid1,'%f',edx);
            fprintf(fid1,',');
            fprintf(fid1,'%f',edxw);
            fprintf(fid1,'\n');
            case 'rings'
            fprintf(fid1,'RECTNG,');
            fprintf(fid1,'%f',j);
            fprintf(fid1,',');
            fprintf(fid1,'%f',c);
            fprintf(fid1,',');
            fprintf(fid1,'%f',e);
            fprintf(fid1,',');
            fprintf(fid1,'%f',ed);
            fprintf(fid1,'\n');
            fprintf(fid1,'RECTNG,');
            fprintf(fid1,'%f',h);
            fprintf(fid1,',');
            fprintf(fid1,'%f',u);
            fprintf(fid1,',');
            fprintf(fid1,'%f',edx);
            fprintf(fid1,',');
            fprintf(fid1,'%f',edxw);
            fprintf(fid1,'\n');
        end
    
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
    fprintf(fid1,'%f',mesh);
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
    !"C:\Program Files\ANSYS Inc\v140\ansys\bin\intel\ANSYS140" -b -i C:\magcode\ansys\total.txt -o C:\magcode\ansys\result.out

        i=i+1;
    end
    %delete 'unnecessary' files
    !del file.BCS
    !del file.emat
    !del file.err
    !del file.esav
    !del file.full
    !del file.rmg
    !del file.stat
end
