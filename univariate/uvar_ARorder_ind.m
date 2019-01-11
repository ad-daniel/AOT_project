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
pVals = 1 : 1 : 15; % noise type
sig1 = 0.1;
sig2 = 0.05;

NOISE_TYPE = NOISE_SWING;
%NOISE_TYPE = NOISE_UNIF;
%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

switch(NOISE_TYPE)
   case NOISE_SWING
      filename = 'uvar_ARorder_ind_swing';
   case NOISE_UNIF
      filename = 'uvar_ARorder_ind_unif';
end

%% preallocate
results = zeros(length(pVals), 3);

%% MAIN
if(F_TASK == F_RECOMPUTE)
   for pp = 1 : length(pVals)
      p = pVals(pp);
      
      decision = -1*ones(1,rep);
      
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
         % independance measure
         [pval_fw, ~] = indtest_hsic(res_fw', fw_data(1:end-p), [], []);
         [pval_bw, ~] = indtest_hsic(res_bw', bw_data(1:end-p), [], []);
         
         decision(a) = indtest_decision(pval_fw, pval_bw, sig1, sig2);
      end
      
      results(pp, 1) = sum( decision == DIR_FW );
      results(pp, 2) = sum( decision == DIR_BW );
      results(pp, 3) = sum( decision == DIR_UNKNOWN );
   end
   
   save(['results/', filename, '.mat'], 'results');
else
   % load last results 
   load(['results/', filename, '.mat'])
end

%% Plot
set(0,'defaultAxesFontSize',16)

fig = figure;
bar(pVals(1:12), results(1:12,:), 'stacked', 'FaceColor','flat'); axis tight;
xticks(pVals(1:12))
xlabel('p'); ylabel('accuracy');
legend('forward', 'backward', 'undecided', 'Location', 'NorthEast');

%% Print
print(fig, '-depsc2', ['images/', filename])

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end