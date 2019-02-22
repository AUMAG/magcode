%% Test torque zz singularities

classdef testcuboidtorque03 < matlab.unittest.TestCase
  
  properties(TestParameter)
    singularity = {'UV','UW','VW'};
    displ = {...
      [0.005; 0.005; 0.12 ], ...
      [0.005; 0.1  ; 0.005], ...
      [0.08 ; 0.005; 0.005], ...
    };
  end
  
  methods(Test, ParameterCombination='sequential')
    
    function test_torque_zz_singularity(testCase,displ,singularity)
      
      magnet_fixed = magnetdefine('type','cuboid','magn',1.3,'magdir',[0 0 1],'dim',[0.03 0.04 0.05]);
      magnet_float = magnetdefine('type','cuboid','magn',1.3,'magdir',[0 0 1],'dim',[0.04 0.05 0.06]);
      
      smidge = 1e-6;
      T1 = magnetforces(magnet_fixed,magnet_float,displ,       'torque');
      T2 = magnetforces(magnet_fixed,magnet_float,displ+smidge,'torque');
      
      check = round([T1,T2],3);
      testCase.verifyThat(check, ~HasNaN, ...
        "singularity "+singularity+": should not have NaN" );
      testCase.verifyEqual(check(:,1), check(:,2), ...
        "singularity "+singularity+" is not consistent with nearby value" );
      
    end
  end
end

% "import" function
function c = HasNaN(varargin)
c = matlab.unittest.constraints.HasNaN(varargin{:});
end