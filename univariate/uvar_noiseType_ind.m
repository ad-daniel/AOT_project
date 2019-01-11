clear all; close all; clc;

%% Include
addpath(genpath('../support/'));
settings;

%% Simulated AR(p) process
rep = 100;  % repetitions for accuracy measurement
T = 5000;   % time series length
p = 1;      % lag order
n = 1;      % dimensionality

%% Settings
rVals = 0.1 : 0.1 : 2; % noise type
sig1 = 0.1;
sig2 = 0.05;

NOISE_TYPE = NOISE_SWING;
%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

filename = 'uvar_noiseType_ind_swing';

%% preallocate
results = zeros(length(rVals), 3);

%% MAIN

if(F_TASK == F_RECOMPUTE)
   for rr = 1 : length(rVals)
      r = rVals(rr);
      
      decision = -1*ones(1,rep);
      
      for a = 1 : rep 
         disp(['doing r = ', num2str(r), ' : rep ', num2str(a), '/', num2str(rep)]);

         % generate noise vector
         eps_t = genNoise(T, n, NOISE_TYPE, r);
         % generate AR coefficients
         %phi = 0.9*rand(n,n,p)-0.45*ones(n,n,p); % avoid poles on unit circle
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
      
      results(rr, 1) = sum( decision == DIR_FW );
      results(rr, 2) = sum( decision == DIR_BW );
      results(rr, 3) = sum( decision == DIR_UNKNOWN );
   end
   
   save(['results/', filename, '.mat'], 'results');
else
   % load last results 
   load(['results/', filename, '.mat'])
end

%% Plot
set(0,'defaultAxesFontSize',14)

fig = figure;
bar(rVals, results, 'stacked', 'FaceColor','flat'); axis tight;
xticks([0.1,0.5,1.0,1.5,2.0] )
xticklabels({'0.1','0.5','1.0', '1.5', '2.0'})
xlabel('r'); ylabel('accuracy');
legend('forward', 'backward', 'undecided', 'Location', 'NorthEast');

%% Print
print(fig, '-depsc2', ['images/', filename])

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end