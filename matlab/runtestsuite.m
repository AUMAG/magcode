%% Wrapper script to run the MAGCODE test suite
%
% This code uses the Matlab Unit Test framework.
% Very little wrapper code is needed to run the test suite, but for
% Continuous Integration especially it's useful to have the following
% to ensure that the path is correctly set up.

addpath(pwd)

%% Which code path?
%
% By default (for CI) we want to run the profiler and code coverage
% plug-ins. For development purposes or testing you may wish to switch
% to |choice=1| to keep things simple and speed things up.

choice = 2;

%% Plain execution
%
% Nice and simple.

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
  plugin = CodeCoveragePlugin.forFolder(pwd+["","/private"],...
    'Producing',CoberturaFormat('coverage.xml'));
  
  runner.addPlugin(plugin);
  results = runner.run(suite);
end

%% Coverage report debugging

if choice == 3
  import matlab.unittest.TestSuite
  import matlab.unittest.TestRunner
  import matlab.unittest.plugins.CodeCoveragePlugin
  import matlab.unittest.plugins.codecoverage.CoberturaFormat
  
  suite  = TestSuite.fromFile('testsuite/testcuboidtorque03.m');
  runner = TestRunner.withTextOutput;
  plugin = CodeCoveragePlugin.forFolder(pwd+["","/private"],...
    'Producing',CoberturaFormat('coverage-debug.xml'));
  
  runner.addPlugin(plugin);
  results = runner.run(suite);
end

%% Finish up
%
% Without this line we don't get the correct exit code if results failed.

assert(all(~[results.Failed]));
