%% Run test suite

choice = 2;

%% Plain execution
%
% This code uses the Matlab Unit Test framework.
% Very little wrapper code is needed to run the test suite, but for
% Continuous Integration especially it's useful to have the following
% to ensure that the path is correctly set up.

if choice == 1
  addpath(pwd)
  results = runtests('testsuite/');
  assert(all(~[results.Failed]));
end


%% Coverage report

if choice == 2
  import matlab.unittest.TestSuite
  import matlab.unittest.TestRunner
  import matlab.unittest.plugins.CodeCoveragePlugin
  import matlab.unittest.plugins.codecoverage.CoberturaFormat
  
  suite  = TestSuite.fromFolder('testsuite/');
  runner = TestRunner.withTextOutput;
  plugin = CodeCoveragePlugin.forFolder('testsuite/',...
    'Producing',CoberturaFormat('CoverageResults.xml'));
  
  runner.addPlugin(plugin);
  results = runner.run(suite);
  assert(all(~[results.Failed]));
end