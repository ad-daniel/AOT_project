function [sign,nullset] = genNullset(vartype, ref, p, T, n, rep)

settings;

switch(vartype)
   case 'uvar'
      nullset = zeros(rep,1);

      for a = 1:rep
         if(~mod(a,25))
            disp(['> nullset a = ', num2str(a), '/', num2str(rep)])
         end
         % generate noise vector
         eps_t = genNoise(T, n, NOISE_GAUSS); % NB: it's unit std
         % adapt noise variance to that of the reference series
         eps_t = std(ref.res) * eps_t;
         % generate time series
         [fw_data, bw_data] = genVAR(ref.phi, ref.X0, eps_t, p, T, n);
         % fit AR process 
         [~,~,~,res_fw] = CBIG_RL2017_ar_mls(fw_data, p);
         [~,~,~,res_bw] = CBIG_RL2017_ar_mls(bw_data, p);
         % gaussianity measure
         [~, ~, deltaKurt_fw(a)] = mardiaKurtosis(res_fw');
         [~, ~, deltaKurt_bw(a)] = mardiaKurtosis(res_bw'); 
         
         nullset(a) = deltaKurt_bw(a) - deltaKurt_fw(a);
      end
      
      % compare nullset to reference
      LI = abs(nullset) > abs(ref.delta_bwfw);
      sign = (rep - sum(LI)) / rep;
      
   case 'mvar'
      nullset = zeros(rep,1);

      for a = 1:rep
         if(~mod(a,25))
            disp(['> nullset a = ', num2str(a), '/', num2str(rep)])
         end         
         MU = mean(ref.res,2);
         SIGMA = cov(ref.res');
         % generate noise vector
         eps_t = genNoise(T, n, NOISE_GAUSS_SPEC, 0,  MU, SIGMA);
         % generate time series
         [fw_data, bw_data] = genVAR(ref.phi, ref.X0, eps_t, p, T, n);
         % fit AR process 
         [~,~,~,res_fw] = CBIG_RL2017_ar_mls(fw_data, p);
         [~,~,~,res_bw] = CBIG_RL2017_ar_mls(bw_data, p);
         % gaussianity measure
         [~, ~, deltaKurt_fw(a)] = mardiaKurtosis(res_fw');
         [~, ~, deltaKurt_bw(a)] = mardiaKurtosis(res_bw'); 
         
         nullset(a) = deltaKurt_bw(a) - deltaKurt_fw(a);
      end
      
      % compare nullset to reference
      LI = abs(nullset) > abs(ref.delta_bwfw);
      sign = (rep - sum(LI)) / rep;
     
   otherwise
      error('Unknown nullset type')
end