%% Test torque zz singularities

classdef test_cuboid_dipole_01 < matlab.unittest.TestCase
  
  properties(TestParameter)
    a = {0.05,0.05};
    displ = {...
      [0; 0; 0.1 ], ...
      [0; 0; 1 ], ...
    };
    dim = {
      [0.02 0.02 0.02], ...
      [0.2  0.2  0.2 ], ...
    };
    roundn = {3,1};
  end
  
  methods(Test, ParameterCombination='sequential')
    
    function test_force(testCase,displ,dim,roundn)
      
      magnet_fixed = magnetdefine('type','cuboid','magn',+1,'magdir',[0 0 1],'dim',dim);
      magnet_float = magnetdefine('type','cuboid','magn',-1,'magdir',[0 0 1],'dim',dim);
      
      f1 = magnetforces(magnet_fixed,magnet_float,displ,'method','auto'  );
      f2 = magnetforces(magnet_fixed,magnet_float,displ,'method','dipole');
      
      check = round([f1(3),f2(3)],roundn);
      assert(check(1)==check(2), "Cuboid force and dipole force should be pretty close here: "+f1(3)+" ~= "+f2(3) );
      
    end
    
    function test_torque(testCase,a,displ,dim,roundn)
      
      magnet_fixed = magnetdefine('type','cuboid','magn',+1,'magdir',[0 0 1],'dim',dim);
      magnet_float = magnetdefine('type','cuboid','magn',-1,'magdir',[0 0 1],'dim',dim);
      
      f1 = magnetforces(magnet_fixed,magnet_float,displ+[a;0;0],'torque','method','auto'  );
      f2 = magnetforces(magnet_fixed,magnet_float,displ+[a;0;0],'torque','method','dipole');
      
      check = round([f1(2),f2(2)],roundn+1);
      assert(check(1)==check(2), "Cuboid torque and dipole torque should be pretty close here: "+f1(2)+" ~= "+f2(2) );
      
    end
    
    function test_sym(testCase,a,displ,dim,roundn)
      
      magnet_fixed = magnetdefine('type','cuboid','magn',+1,'magdir',[0 0 1],'dim',dim);
      magnet_float = magnetdefine('type','cuboid','magn',-1,'magdir',[0 0 1],'dim',dim);
      
      [f1,t1] = magnetforces(magnet_fixed,magnet_float,displ-[a;0;0],'force','torque','method','dipole');
      [f2,t2] = magnetforces(magnet_fixed,magnet_float,displ+[a;0;0],'force','torque','method','dipole');
      
      check = [f1(1),f2(1)];
      assert(abs(check(1)+check(2))<eps, "Odd symmetry for forces" );
      
      check = [f1(3),f2(3)];
      assert(abs(check(1)-check(2))<eps, "Even symmetry for forces" );
      
      check = [t1(2),t2(2)];
      assert(abs(check(1)+check(2))<eps, "Odd symmetry for torques" );
      
    end
    
  end
end

