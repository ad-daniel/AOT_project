clear all; close all; clc;

%% Include
addpath(genpath('../support/'));
settings;

%% Simulated AR(p) process
rep = 10000;  % repetitions for accuracy measurement
T = 5000;   % time series length
p = 1;      % lag order
n = 1;      % dimensionality

NOISE_TYPE = NOISE_UNIF;
filename = 'uvar_nullset_unif';

%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

%% MAIN
if(F_TASK == F_RECOMPUTE)   
   %% Generate reference time series
   % generate uniform noise vector
   eps_t = genNoise(T, n, NOISE_TYPE);
   % generate AR coefficients
   phi = genCoeff(n,p);
   % generate time series
   X0 = zeros(n,1); % intercept
   [fw_data, bw_data] = genVAR(phi, X0, eps_t, p, T, n);
   % fit AR process 
   [~,B,~,res_fw] = CBIG_RL2017_ar_mls(fw_data, p);
   [~,~,~,res_bw] = CBIG_RL2017_ar_mls(bw_data, p);
   % gaussianity measure
   [~, ~, deltaKurt_fw] = mardiaKurtosis(res_fw');
   [~, ~, deltaKurt_bw] = mardiaKurtosis(res_bw');
   
   delta_bwfw = deltaKurt_bw - deltaKurt_fw;
   
   % store reference info in structure
   ref.X0 = B(:,1);
   phi_id = zeros(n,n,p);
   for pp = 1:p
      phi_id(:,:,pp) = B(:,2+n*(pp-1):n*(pp-1)+n+1);
   end
   ref.phi = phi_id;
   ref.delta_bwfw = delta_bwfw;
   ref.res = res_fw;
   
   % generate nullset
   [sign,nullset] = genNullset('uvar', ref, p, T, n, rep);
   
   % save data for backup
   results.ref = ref;
   results.nullset = nullset;
   results.sign = sign;
   
   save(['results/', filename, '.mat'], 'results');
else
   % load last results 
   load(['results/', filename, '.mat'])
end

%% Plot
fig = figure;
stem(results.nullset, 'b.'); grid minor; hold on;
line([0, rep], [results.ref.delta_bwfw,results.ref.delta_bwfw], 'Color', 'r');
%line([0, rep], [-results.ref.delta_bwfw,-results.ref.delta_bwfw], 'Color', 'r');hold off; axis tight;
legend('nullset', 'original TS','Location','SouthEast');
xlabel('attempt'); ylabel('k4_{bw-fw}')
%title(['Result more significant than ', num2str(results.sign*100), '% of the nullset'])

%% Print
print(fig, '-depsc2', ['images/', filename])

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end