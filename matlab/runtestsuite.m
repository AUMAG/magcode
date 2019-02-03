%% Wrapper script to run the MAGCODE test suite

addpath(pwd)

choice = 2;

%% Plain execution
%
% This code uses the Matlab Unit Test framework.
% Very little wrapper code is needed to run the test suite, but for
% Continuous Integration especially it's useful to have the following
% to ensure that the path is correctly set up.

if choice == 1
  results = runtests('testsuite/');
end


%% Coverage report
%
% This uses the far more clunky interface to run tests but produces
% something very interesting...

if choice == 2
  import matlab.unittest.TestSuite
  import matlab.unittest.TestRunner
  import matlab.unittest.plugins.CodeCoveragePlugin
  import matlab.unittest.plugins.codecoverage.CoberturaFormat
  
  suite  = TestSuite.fromFolder('testsuite/');
  runner = TestRunner.withTextOutput;
  plugin = CodeCoveragePlugin.forFolder(pwd,'IncludingSubfolders',true,...
    'Producing',CoberturaFormat('coverage.xml'));
  
  runner.addPlugin(plugin);
  results = runner.run(suite);
end

%% Finish up

assert(all(~[results.Failed]));
