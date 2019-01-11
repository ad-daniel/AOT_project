function [pval, stat] = indtest_template(X, Y, Z, pars)
% function [pval, stat] = indtest_template(X, Y, Z, pars)
%
% Template for writing indtest_* functions that perform a (conditional) independence test.
% Its sole purpose is to define a common interface for the various indtest_* functions;
% it does not in itself perform any useful function.
%
%
% INPUT:
%   X         Nxd1 matrix of samples (N data points, d1 dimensions)
%   Y         Nxd2 matrix of samples (N data points, d2 dimension)
%   Z         Nxd3 matrix of samples (N data points, d3 dimensions)
%   pars      structure containing parameters for the independence test
%
% OUTPUT:
%   pval      p value of the test
%   stat      test statistic


  help('indtest_template');

return
