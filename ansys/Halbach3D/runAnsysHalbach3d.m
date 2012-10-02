function [] = runAnsysHalbach3d(a,b,k,l,m,n,N,h1,c,theta,angle,e,f,d,h2,g,h,mesh)
quote = '''';
fid1 = fopen('halbach3d.txt','wt');
fprintf(fid1,'/PREP7\n');
fprintf(fid1,'/TITLE, HALBACH CUBES\n');
fprintf(fid1,'ET,1,SOLID98\n');
fprintf(fid1,'MP,MURX,1,1\n');

% Draw air
fprintf(fid1,'BLOCK,');
fprintf(fid1,'%f',-a/2);
fprintf(fid1,',');
fprintf(fid1,'%f',a/2);
fprintf(fid1,',');
fprintf(fid1,'%f',-b/2);
fprintf(fid1,',');
fprintf(fid1,'%f',b/2);
fprintf(fid1,',');
fprintf(fid1,'%f',-c/2);
fprintf(fid1,',');
fprintf(fid1,'%f',c/2);
fprintf(fid1,'\n');

% Draw outer 'cubes'
for i = 1:N
    if i == 1
        % For first magnet MGYY=0, so is not used
        fprintf(fid1,'MP,MURX,');
        fprintf(fid1,'%f',i+1);
        fprintf(fid1,',1.05\n');
        fprintf(fid1,'MP,MGXX,');
        fprintf(fid1,'%f',i+1);
        fprintf(fid1,',');
        fprintf(fid1,'%f',h1);
        fprintf(fid1,'\n');

        fprintf(fid1,'LOCAL,');
        fprintf(fid1,'%f',10+i);
        fprintf(fid1,',0,');
        fprintf(fid1,'%f',d);
        fprintf(fid1,',0,0\n');
    else 
        % Rest of the magnets
        fprintf(fid1,'MP,MURX,');
        fprintf(fid1,'%f',i+1);
        fprintf(fid1,',1.05\n');
        fprintf(fid1,'MP,MGXX,');
        fprintf(fid1,'%f',i+1);
        fprintf(fid1,',');
        fprintf(fid1,'%f',h1*cos((i-1)*(N/4+1)*theta));
        fprintf(fid1,'\n');
        fprintf(fid1,'MP,MGYY,');
        fprintf(fid1,'%f',i+1);
        fprintf(fid1,',');
        fprintf(fid1,'%f',h1*sin((i-1)*(N/4+1)*theta));
        fprintf(fid1,'\n');
        fprintf(fid1,'LOCAL,');
        fprintf(fid1,'%f',10+i);
        fprintf(fid1,',0,');
        fprintf(fid1,'%f',d*cos((i-1)*theta));
        fprintf(fid1,',');
        fprintf(fid1,'%f',d*sin((i-1)*theta));
        fprintf(fid1,',0,');
        fprintf(fid1,'%f',((i-1)*angle));
        fprintf(fid1,'\n');
    end

    fprintf(fid1,'WPCSYS,1,');
    fprintf(fid1,'%f',10+i);
    fprintf(fid1,'\n');
    fprintf(fid1,'BLOCK,');
    fprintf(fid1,'%f',-f/2);
    fprintf(fid1,',');
    fprintf(fid1,'%f',f/2);
    fprintf(fid1,',');
    fprintf(fid1,'%f',-g/2);
    fprintf(fid1,',');
    fprintf(fid1,'%f',g/2);
    fprintf(fid1,',');
    fprintf(fid1,'%f',-h/2);
    fprintf(fid1,',');
    fprintf(fid1,'%f',h/2);
    fprintf(fid1,'\n');
    i = i+1;  
%     end
end

% Draw inner 'cubes'
for i = 1:N
    if i == 1
        % For first magnet MGYY=0, so is not used
        fprintf(fid1,'MP,MURX,');
        fprintf(fid1,'%f',i+1+N);
        fprintf(fid1,',1.05\n');
        fprintf(fid1,'MP,MGXX,');
        fprintf(fid1,'%f',i+1+N);
        fprintf(fid1,',');
        fprintf(fid1,'%f',h2);
        fprintf(fid1,'\n');
        
        fprintf(fid1,'LOCAL,');
        fprintf(fid1,'%f',10+i);
        fprintf(fid1,',0,');
        fprintf(fid1,'%f',e);
        fprintf(fid1,',0,0\n');
    else
        % Rest of the magnets
        fprintf(fid1,'MP,MURX,');
        fprintf(fid1,'%f',i+1+N);
        fprintf(fid1,',1.05\n');
        fprintf(fid1,'MP,MGXX,');
        fprintf(fid1,'%f',i+1+N);
        fprintf(fid1,',');
        fprintf(fid1,'%f',h2*cos(-(i-1)*(N/4-1)*(theta)));
        fprintf(fid1,'\n');
        fprintf(fid1,'MP,MGYY,');
        fprintf(fid1,'%f',i+1+N);
        fprintf(fid1,',');
        fprintf(fid1,'%f',h2*sin(-(i-1)*(N/4-1)*(theta)));
        fprintf(fid1,'\n');
        
        fprintf(fid1,'LOCAL,');
        fprintf(fid1,'%f',10+i);
        fprintf(fid1,',0,');
        fprintf(fid1,'%f',e*cos((i-1)*theta));
        fprintf(fid1,',');
        fprintf(fid1,'%f',e*sin((i-1)*theta));
        fprintf(fid1,',0,');
        fprintf(fid1,'%f',((i-1)*angle));
        fprintf(fid1,'\n');
    end
    fprintf(fid1,'WPCSYS,1,');
    fprintf(fid1,'%f',10+i);
    fprintf(fid1,'\n');
    fprintf(fid1,'BLOCK,');
    fprintf(fid1,'%f',-k/2);
    fprintf(fid1,',');
    fprintf(fid1,'%f',k/2);
    fprintf(fid1,',');
    fprintf(fid1,'%f',-l/2);
    fprintf(fid1,',');
    fprintf(fid1,'%f',l/2);
    fprintf(fid1,',');
    fprintf(fid1,'%f',-m/2+n);
    fprintf(fid1,',');
    fprintf(fid1,'%f',m/2+n);
    fprintf(fid1,'\n');
    i = i+1;
end

fprintf(fid1,'WPCSYS,-1,0\n');          % Back to original working plane

% Assign material properties outer ring
for i = 1:N

    fprintf(fid1,'VSEL,S,VOLU,,');
    fprintf(fid1,'%f',i+1);
    fprintf(fid1,'\n');
    fprintf(fid1,'VATT,');
    fprintf(fid1,'%f',i+1);
    fprintf(fid1,',,');
    fprintf(fid1,'%f',1);
    fprintf(fid1,'\n');

    i=i+1;
end

% Assign material properties inner ring
for i = 1:N
    fprintf(fid1,'VSEL,S,VOLU,,');
    fprintf(fid1,'%f',i+1+N);
    fprintf(fid1,'\n');
    fprintf(fid1,'VATT,');
    fprintf(fid1,'%f',i+1+N);
    fprintf(fid1,',,');
    fprintf(fid1,'%f',1);
    fprintf(fid1,'\n');

    i = i+1;
end

fprintf(fid1,'VSEL,ALL\n');
fprintf(fid1,'VSBV,1,ALL,,,KEEP\n');
fprintf(fid1,'VSEL,U,VOLU,,2,');
fprintf(fid1,'%f',2*N+1);
fprintf(fid1,'\n');
fprintf(fid1,'VATT,1,,1\n');

% Boundary conditions
fprintf(fid1,'ALLSEL\n');
fprintf(fid1,'ASEL,S,EXT\n');
fprintf(fid1,'DL,ALL,,AZ,0\n');

fprintf(fid1,'ALLSEL\n');
fprintf(fid1,'SMRTSIZE,1\n');
% fprintf(fid1,'ESIZE,,');
% fprintf(fid1,'%f',mesh);
% fprintf(fid1,'\n');
fprintf(fid1,'VMESH,ALL\n');

% Apply force flags
fprintf(fid1,'ESEL,S,MAT,,');
fprintf(fid1,'%f',N+1);
fprintf(fid1,',');
fprintf(fid1,'%f',2*N);
fprintf(fid1,'\n');
fprintf(fid1,'CM,mag1,ELEM\n');
fprintf(fid1,'FMAGBC,');
fprintf(fid1,'%s',quote);
fprintf(fid1,'mag1');
fprintf(fid1,'%s',quote);
fprintf(fid1,'\n');

fprintf(fid1,'ALLSEL\n');
fprintf(fid1,'/SOL\n');
fprintf(fid1,'MAGSOLV\n');
fprintf(fid1,'FINISH\n');


%% Code for calculating forces
fprintf(fid1,'/POST1\n');
fprintf(fid1,'fmagsum,');
fprintf(fid1,'%s',quote);
fprintf(fid1,'mag1');
fprintf(fid1,'%s',quote);
fprintf(fid1,'\n');

% fprintf(fid1,'etab,fmgx1,fmag,x\n');
% fprintf(fid1,'etab,fmgy1,fmag,y\n');
% fprintf(fid1,'etab,fmgz1,fmag,z\n');

fprintf(fid1,'SSUM\n');

fprintf(fid1,'*dim,ff,ARRAY,6\n');
fprintf(fid1,'*get,ff(1),ssum,fmx_x\n');
fprintf(fid1,'*get,ff(2),ssum,fvw_x\n');
fprintf(fid1,'*get,ff(3),ssum,fmx_y\n');
fprintf(fid1,'*get,ff(4),ssum,fvw_y\n');
fprintf(fid1,'*get,ff(5),ssum,fmx_z\n');
fprintf(fid1,'*get,ff(6),ssum,fvw_z\n');

fprintf(fid1,'/PAGE,5000,5000,5000,5000\n');
fprintf(fid1,'/HEADER,OFF,OFF,OFF,OFF,OFF,OFF\n');
fprintf(fid1,'/OUTPUT,');
fprintf(fid1,'%s',quote);
fprintf(fid1,'Halbach3Ddata');
fprintf(fid1,'%s',quote);
fprintf(fid1,',txt,,APPEND\n');

% Write data

fprintf(fid1,'*VWRITE,ff(1),ff(2),ff(3),ff(4),ff(5),ff(6)\n');
fprintf(fid1,'(E14.6,2X,E14.8,2X,E14.8,2X,E14.8,2X,E14.8,2X,E14.8)\n');
fprintf(fid1,'/OUTPUT\n');


fprintf(fid1,'/EOF\n');
fclose(fid1);

%%
% Run the ANSYS batch file, locations must be correct
% !"C:\Program Files\ANSYS Inc\v140\ansys\bin\intel\ANSYS140" -b -i C:\magcode\ansys\Halbach3D\halbach3d.txt -o C:\magcode\ansys\Halbach3D\halbach3dresult.out

% Delete 'unnecessary' files
!del file.BCS
!del file.emat
!del file.err
!del file.esav
!del file.full
!del file.rst
!del file.stat



end