clear all; close all; clc;

%% Include
addpath(genpath('../support/'));
settings;

%% Simulated AR(p) process
rep = 100;  % repetitions for accuracy measurement
T = 8000;   % time series length
p = 1;      % lag order
nmax = 10;  % dimensionality
r = 0.5;    % noise type

%% Settings

NOISE_TYPE = NOISE_SWING;
%NOISE_TYPE = NOISE_UNIF;
%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

switch(NOISE_TYPE)
   case NOISE_SWING
      filename = 'mvar_Ngaussians_gauss_swing';
   case NOISE_UNIF
      filename = 'mvar_Ngaussians_gauss_unif';
end

%% preallocate
results = zeros(length(0:nmax), 2);

%% MAIN
if(F_TASK == F_RECOMPUTE)
   for ndims = 0 : nmax
      deltaKurt_fw = zeros(1,rep);
      deltaKurt_bw = zeros(1,rep);

      for a = 1 : rep 
         disp(['doing #n = ', num2str(ndims), ' : rep ', num2str(a), '/', num2str(rep)]);

         % generate noise vector
         eps_t = zeros(nmax,T);
         for nn = 1:nmax                        
            if(nn <= ndims)            
               eps_t(nn,:) = genNoise(T,1,NOISE_GAUSS);
            else
               eps_t(nn,:) = genNoise(T,1,NOISE_TYPE, r);
            end
         end         

         % generate AR coefficients
         phi = genCoeff(nmax,p);
         % generate time series
         [fw_data, bw_data] = genVAR(phi, zeros(nmax,1), eps_t, p, T, nmax);
         % fit AR process 
         [~,~,~,res_fw] = CBIG_RL2017_ar_mls(fw_data, p);
         [~,~,~,res_bw] = CBIG_RL2017_ar_mls(bw_data, p);
         % gaussianity measure
         [~, ~, deltaKurt_fw(a)] = mardiaKurtosis(res_fw');
         [~, ~, deltaKurt_bw(a)] = mardiaKurtosis(res_bw');
      end

      LI = deltaKurt_fw > deltaKurt_bw;
      
      results(ndims+1, 1) = sum(LI);
      results(ndims+1, 2) = rep - sum(LI);
   end
   
   save(['results/', filename, '.mat'], 'results');
else
   % load last results 
   load(['results/', filename, '.mat'])
end

%% Plot
set(0,'defaultAxesFontSize',16)

fig = figure;
bar(0:nmax, results,'stacked', 'FaceColor','flat'); axis tight;
xticks(0:nmax)
xlabel('# of gaussian dimensions'); ylabel('accuracy');
legend('forward', 'backward', 'Location', 'SouthEast');

%% Print
print(fig, '-depsc2', ['images/', filename])

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end