clear all; close all; clc;


%% Include
addpath(genpath('../support/'));
settings;

%% Simulated AR(p) process
rep = 100;  % repetitions for accuracy measurement
T = 5000;   % time series length
n = 1;      % dimensionality
r = 0.25;   % noise type

%% Settings
pVals = 1:1:15;

NOISE_TYPE = NOISE_SWING;
%NOISE_TYPE = NOISE_UNIF;
%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

switch(NOISE_TYPE)
   case NOISE_SWING
      filename = 'uvar_ARorder_gauss_swing';
   case NOISE_UNIF
      filename = 'uvar_ARorder_gauss_unif';
end

%% preallocate
results = zeros(length(pVals), 2);

%% MAIN
if(F_TASK == F_RECOMPUTE)
   for pp = 1 : length(pVals)
      p = pVals(pp);

      deltaKurt_fw = zeros(1,rep);
      deltaKurt_bw = zeros(1,rep);

      for a = 1 : rep 
         disp(['doing p = ', num2str(p), ' : rep ', num2str(a), '/', num2str(rep)]);

         % generate noise vector
         eps_t = genNoise(T, n, NOISE_TYPE, r);
         % generate AR coefficients
         phi = genCoeff(n,p);
         % generate time series
         [fw_data, bw_data] = genVAR(phi, zeros(n,1), eps_t, p, T, n);
         % fit AR process 
         [~,~,~,res_fw] = CBIG_RL2017_ar_mls(fw_data, p);
         [~,~,~,res_bw] = CBIG_RL2017_ar_mls(bw_data, p);
         % gaussianity measure
         [~, ~, deltaKurt_fw(a)] = mardiaKurtosis(res_fw');
         [~, ~, deltaKurt_bw(a)] = mardiaKurtosis(res_bw');
      end

      LI = deltaKurt_fw > deltaKurt_bw;
      
      results(pp, 1) = sum(LI);
      results(pp, 2) = rep - sum(LI);
   end
   
   save(['results/', filename, '.mat'], 'results');
else
   % load last results 
   load(['results/', filename, '.mat'])
end

%% Plot
set(0,'defaultAxesFontSize',16)
fig = figure;
bar(pVals(1:12), results(1:12,:),'stacked', 'FaceColor','flat'); axis tight;
xticks(pVals(1:12))
xlabel('p'); ylabel('accuracy');
legend('forward', 'backward', 'Location', 'SouthEast');

%% Print
print(fig, '-depsc2', ['images/', filename])

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end