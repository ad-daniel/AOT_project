function [mKurt, expKurt, delta] = mardiaKurtosis(X)
% INPUT
%   X - multivariate data matrix (T x n) where T = datapoints, n = dimensionality 
% OUTPUT
%   mKurt   : mardia kurtosis
%   expKurt : expected kurtosis of pure gaussian
%   delta   : distance between actual and expected kurtosis

[T,n] = size(X);

if(n > T)
   error('data provided in the wrong form')
end

%% covariance matrix normalized by n
COVm = cov(X,1); % by default, it normalizes by (n-1)

%% compute squared Mahalonobis distance
% remove mean
dfm = zeros(T,n); 
for i = 1:n
   dfm(:,i) = X(:,i) - mean(X(:,i));
end

D = dfm * inv(COVm) * dfm'; % squared Mahalonobis distance

% mSkew = (sum(sum(D.^3))) / (T^2);  % multivariate skewness
mKurt   = trace(D.^2) / T;           % multivariate kurtosis
expKurt = n*(n + 2);                 % kurtosis of a pure gaussian
delta   = abs(mKurt-expKurt);   % gaussianity distance

