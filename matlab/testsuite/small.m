%% Test torque zz singularities

displ = {...
  [0.005; 0.005; 0.12 ], ...
  [0.005; 0.1  ; 0.005], ...
  [0.08 ; 0.005; 0.005], ...
  };
displ = displ{1};
magdir1 = {'+x','+y','+z'};
magdir2 = {'+x','+y','+z'};

magnet_fixed = magnetdefine('type','cuboid','magn',1.3,'magdir',magdir1{2},'dim',[0.03 0.04 0.05]);
magnet_float = magnetdefine('type','cuboid','magn',1.3,'magdir',magdir2{1},'dim',[0.04 0.05 0.06]);

smidge = 1e-6;
F1 = magnetforces(magnet_fixed,magnet_float,displ,        'force')
F2 = magnetforces(magnet_fixed,magnet_float,displ+smidge, 'force')
F3 = magnetforces(magnet_fixed,magnet_float,[2*displ,displ],'force');
F3 = F3(:,2)
