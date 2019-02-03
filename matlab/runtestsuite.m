%% Run test suite
%
% This code uses the Matlab Unit Test framework.
% Very little wrapper code is needed to run the test suite, but for
% Continuous Integration especially it's useful to have the following
% to ensure that the path is correctly set up.

addpath(pwd)

results = runtests('IncludeSubfolders', true);

assert(all(~[results.Failed]));
