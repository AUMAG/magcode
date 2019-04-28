%% Test torque zz singularities

classdef test_cuboid_force_singu < matlab.unittest.TestCase
  
  properties(TestParameter)
    displ = {...
      [0.005; 0.005; 0.12 ], ...
      [0.005; 0.1  ; 0.005], ...
      [0.08 ; 0.005; 0.005], ...
    };
    magdir1 = {'+x','+y','+z'};
    magdir2 = {'+x','+y','+z'};
  end
  
  methods(Test)
    
    function test_singularities(testCase,displ,magdir1,magdir2)
      
      magnet_fixed = magnetdefine('type','cuboid','magn',1.3,'magdir',magdir1,'dim',[0.03 0.04 0.05]);
      magnet_float = magnetdefine('type','cuboid','magn',1.3,'magdir',magdir2,'dim',[0.04 0.05 0.06]);
      
      smidge = 1e-6;
      F1 = magnetforces(magnet_fixed,magnet_float,displ,        'force');
      F2 = magnetforces(magnet_fixed,magnet_float,displ+smidge, 'force');
      F3 = magnetforces(magnet_fixed,magnet_float,[2*displ,displ],'force');
      F3 = F3(:,2);
      
      tol = 1e-3;
      testCase.verifyThat([F1,F2,F3], ~HasNaN, ...
        "singularity: should not have NaN" );
      testCase.verifyLessThan(abs(F1-F2),tol*abs(F1), ...
        "singularity case is not consistent with nearby value" );
      testCase.verifyEqual(F1,F3, ...
        "singularity case not vectorised" );
      
    end
  end
end

% "import" function
function c = HasNaN(varargin)
c = matlab.unittest.constraints.HasNaN(varargin{:});
end