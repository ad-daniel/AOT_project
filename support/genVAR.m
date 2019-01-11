function [fw_data, bw_data] = genVAR(phi, X0, eps_t, p, T, n)
% Generates an (V)AR[p] process with a given noise vector
% INPUT:
%   phi     : AR coefficients (n x n x p)
%   X0      : intercept (n x 1)
%   eps_t   : noise vector (n x T)
%   p       : AR order
%   T       : time series datapoints
%   n       : dimensionality (1 for univariate)
% OUTPUT: 
%   fw_data : time series in the chronologically direct direction
%   bw_data : time series in the chronologically reverse direction

%% Sanitize
if(size(X0,1) ~= n)
   error('Intercept - wrong dimensionality')
end
if (size(eps_t,1) ~= n || size(eps_t,2) ~= T)
   error('Noise vector - wrong size')
end
if (size(phi,1) ~= n || size(phi,2) ~= n || size(phi,3) ~= p)
   error('AR coefficients - wrong size')
end

%% Preallocate and init
X = zeros(n,T);
X(:,1:p) = eps_t(:,1:p) + repmat(X0,1,p);

for i = p+1:T
   % add intercept
   X(:,i) = X(:,i) + X0;
   
   % add lags
   for m = 1:p
      X(:,i) = X(:,i) + phi(:,:,m)*X(:,i-m);
   end
   % add noise vector
   X(:,i) = X(:,i) + eps_t(:,i);
end

%% Check data sanity
if(nnz(isnan(X)) || nnz(isinf(X)))
   error('Time series generated is divergent')
end

%% Return time series
fw_data = X';
bw_data = flipud(fw_data);