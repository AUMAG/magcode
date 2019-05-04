
classdef test_cuboid_stiffness < matlab.unittest.TestCase
  
  properties(TestParameter)
    magdir1 = {...
      [1;0;0],...
      [0;1;0],...
      [0;0;1],...
      [0.55;0.32;0.76],...
    }
    magdir2 = {...
      [1;0;0],...
      [0;1;0],...
      [0;0;1],...
      [0.35;0.61;0.71],...
    }
    displ = {...
      [0.1;  0.09; 0.11],...
      [0.05; 0.01; 0.01],...
    };
  end
  
  methods(Test)
    
    function test_stiffness(testCase,magdir1,magdir2,displ)
                
      magnet_fixed = magnetdefine('type','cuboid','magn',1,'magdir',magdir1,'dim',[0.03  0.04  0.05]);
      magnet_float = magnetdefine('type','cuboid','magn',1,'magdir',magdir2,'dim',[0.055 0.045 0.035]);
      
      incr = 1e-6;
      xi = [incr;0;0];
      yi = [0;incr;0];
      zi = [0;0;incr];
      
      fi = magnetforces(magnet_fixed,magnet_float,[displ-xi,displ+xi,displ-yi,displ+yi,displ-zi,displ+zi]);
      k1x = -(fi(1,2)-fi(1,1))/incr/2;
      k1y = -(fi(2,4)-fi(2,3))/incr/2;
      k1z = -(fi(3,6)-fi(3,5))/incr/2;
      k_man = [k1x;k1y;k1z];
      k_all = magnetforces(magnet_fixed,magnet_float,displ,'stiffness');
      
      testCase.verifyEqual( round(k_all,3), round(k_man,3) , ...
        "stiffnesses not consistent with manually calculated values from forces" );

      
    end
    
  end
  
end




