
clear all
disp('=================')
fprintf('TEST cylinder eccentric force consistency #1: ')

magnet_fixed = magnetdefine('type','cylinder',...
  'dim',[0.02 0.02],'magn',1.3,'dir',[0;0;1],'magdir',[0;0;1]);
magnet_float = magnetdefine('type','cylinder',...
  'dim',[0.02 0.02],'magn',1.3,'dir',[0;0;1],'magdir',[0;0;1]);

zdispl = 0.05;
F1 = magnetforces(magnet_fixed,magnet_float,[    0; 0; zdispl]);
F2 = magnetforces(magnet_fixed,magnet_float,[ 1e-9; 0; zdispl]);

check = round(1e6*[F1(3) F2(3)]);
assert( check(1) == check(2) , 'Inconsistent force with negligible offset (coaxial %i ~= eps offset %i)',check(1),check(2));

fprintf('passed\n')


fprintf('TEST cylinder eccentric force consistency #2: ')

magnet_fixed = magnetdefine('type','cylinder',...
  'dim',[0.02 0.02],'magn',1.3,'dir',[0;0;1],'magdir',[0;0;1]);
magnet_float = magnetdefine('type','cylinder',...
  'dim',[0.015 0.015],'magn',1.3,'dir',[0;0;1],'magdir',[0;0;1]);

zdispl = 0.05;
F1 = magnetforces(magnet_fixed,magnet_float,[    0; 0; zdispl]);
F2 = magnetforces(magnet_fixed,magnet_float,[ 1e-9; 0; zdispl]);

check = round(1e6*[F1(3) F2(3)]);
assert( check(1) == check(2) , 'Inconsistent force with negligible offset (coaxial %i ~= eps offset %i)',check(1),check(2));

fprintf('passed\n')
disp('=================')