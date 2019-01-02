%% Test torque zy singularities

%% Test setup

clear all
testname = [mfilename,'.m'];
separator = repmat('=',[1,6+numel(testname)]);
disp(separator);
disp(['== ',testname,' =='])
disp(separator);

%% Magnet setup

magnet_fixed = magnetdefine('type','cuboid','magn',1.3,'magdir',[0 0 1],'dim',[0.03 0.04 0.05]);
magnet_float = magnetdefine('type','cuboid','magn',1.3,'magdir',[0 1 0],'dim',[0.04 0.05 0.06]);

smidge = 1e-6;
prec   = 1e3;

%% Test 1

fprintf('TEST cuboid torque zy singularity UV: ')

displ = [0.005; 0.005; 0.12];
T1 = magnetforces(magnet_fixed,magnet_float,displ,'torque');
T2 = magnetforces(magnet_fixed,magnet_float,displ+smidge,'torque');

check = round([T1,T2]*prec);
assert( all(~isnan(check(:))) , 'UV no nans' )
assert( all(check(:,1)==check(:,2)) , 'UV singularity consistent' )

fprintf('passed\n')

%% Test 2

fprintf('TEST cuboid torque zy singularity UW: ')

displ = [0.005; 0.1; 0.005];
T1 = magnetforces(magnet_fixed,magnet_float,displ,'torque');
T2 = magnetforces(magnet_fixed,magnet_float,displ+smidge,'torque');

check = round([T1,T2]*prec);
assert( all(~isnan(check(:))) , 'UW no nans' )
assert( all(check(:,1)==check(:,2)) , 'UW singularity consistent' )

fprintf('passed\n')

%% Test 3

fprintf('TEST cuboid torque zy singularity WV: ')

displ = [0.08; 0.005; 0.005];
T1 = magnetforces(magnet_fixed,magnet_float,displ,'torque');
T2 = magnetforces(magnet_fixed,magnet_float,displ+smidge,'torque');

check = round([T1,T2]*prec);
assert( all(~isnan(check(:))) , 'UW no nans' )
assert( all(check(:,1)==check(:,2)) , 'UW singularity consistent' )

fprintf('passed\n')

%% End

disp(separator);




