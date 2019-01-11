clear all; close all; clc;

%% Include
addpath(genpath('../support/'));
settings;

%% Simulated AR(p) process
rep = 100;  % repetitions for accuracy measurement
T = 8000;   % time series length
p = 1;      % lag order
r = 0.5;    % noise type

%% Settings
nVals = [5,10,15,20,30,50,100];

NOISE_TYPE = NOISE_SWING;
%NOISE_TYPE = NOISE_UNIF;
%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

switch(NOISE_TYPE)
   case NOISE_SWING
      filename = 'mvar_dimensionality_gauss_swing';
   case NOISE_UNIF
      filename = 'mvar_dimensionality_gauss_unif';
end

%% preallocate
results = zeros(length(nVals), 2);

%% MAIN
if(F_TASK == F_RECOMPUTE)
   for nn = 1 : length(nVals)
      n = nVals(nn);

      deltaKurt_fw = zeros(1,rep);
      deltaKurt_bw = zeros(1,rep);

      for a = 1 : rep 
         disp(['doing n = ', num2str(n), ' : rep ', num2str(a), '/', num2str(rep)]);

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
      
      results(nn, 1) = sum(LI);
      results(nn, 2) = rep - sum(LI);
   end
   
   save(['results/', filename, '.mat'], 'results');
else
   % load last results 
   load(['results/', filename, '.mat'])
end

%% Plot
set(0,'defaultAxesFontSize',16)

fig = figure;
bar(results,'stacked', 'FaceColor','flat'); axis tight;
xticklabels(nVals)
xlabel('n dimensional'); ylabel('accuracy');
legend('forward', 'backward', 'Location', 'SouthEast');

%% Print
print(fig, '-depsc2', ['images/', filename])

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end