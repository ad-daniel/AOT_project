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

NOISE_TYPE = NOISE_SWING;

%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

filename = 'uvar_noiseType_gauss_swing';

%% preallocate
results = zeros(length(rVals), 2);

%% MAIN

if(F_TASK == F_RECOMPUTE)
   for rr = 1 : length(rVals)
      r = rVals(rr);

      deltaKurt_fw = zeros(1,rep);
      deltaKurt_bw = zeros(1,rep);

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
         % gaussianity measure
         [~, ~, deltaKurt_fw(a)] = mardiaKurtosis(res_fw');
         [~, ~, deltaKurt_bw(a)] = mardiaKurtosis(res_bw');
      end

      LI = deltaKurt_fw > deltaKurt_bw;
      
      results(rr, 1) = sum(LI);
      results(rr, 2) = rep - sum(LI);
   end
   
   save(['results/', filename, '.mat'], 'results');
else
   % load last results 
   load(['results/', filename, '.mat'])
end

%% Plot
set(0,'defaultAxesFontSize',14)

fig = figure;
bar(rVals, results,'stacked', 'FaceColor','flat'); axis tight;
xticks([0.1,0.5,1.0,1.5,2.0] )
xticklabels({'0.1','0.5','1.0', '1.5', '2.0'})
xlabel('r'); ylabel('accuracy');
legend('forward', 'backward', 'Location', 'SouthEast');

%% Print
print(fig, '-depsc2', ['images/', filename])

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end