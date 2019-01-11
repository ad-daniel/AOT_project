function eps_t = genNoise(T, n, NOISE_TYPE, r, MU, SIGMA)

settings;

eps_t = zeros(n,T);

switch(NOISE_TYPE)
   case NOISE_UNIF
      for nn = 1:n
         eps_t(nn,:) = rand(1,T)-0.5*ones(1,T);
         eps_t(nn,:) = eps_t(nn,:) ./ std(eps_t(nn,:));      % unit std
      end
   case NOISE_SWING
      if(nargin < 4)
         error('Provide r value for NOISE_SWING')
      end
       
      for nn = 1:n
         Z = randn(1,T);
         eps_t(nn,:) = sign(Z) .* abs(Z) .^ r;   
         eps_t(nn,:) = eps_t(nn,:) ./ std(eps_t(nn,:));      % unit std
      end
   case NOISE_GAUSS
      for nn = 1:n
         eps_t(nn,:) = randn(1,T);
         eps_t(nn,:) = eps_t(nn,:) ./ std(eps_t(nn,:));      % unit std
      end
      
   case NOISE_GAUSS_SPEC
         % https://stats.stackexchange.com/questions/120179/generating-data-with-a-given-sample-covariance-matrix      
         % create innovations that have same covariance matrix as residuals
         eps_t = randn(T, n); % normal noise
         eps_t = bsxfun(@minus, eps_t, mean(eps_t));
         eps_t = eps_t * inv(chol(cov(eps_t)));
         eps_t = eps_t * chol(SIGMA); % now cov(eps_t) = SIGMA
         eps_t = eps_t';
         %eps_t = eps_t + (MU - mean(eps_t, 2)); % now eps_t has mean MU
         
   otherwise
      error('Noise type unknown')
end
