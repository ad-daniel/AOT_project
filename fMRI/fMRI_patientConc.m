clc; clear all; close all;

addpath(genpath('../support/'));
addpath('dataset');
settings;

%% Load rest state
load X_RS_21subjs.mat;
p = 1;

rep = 100;     % nullset dimension

%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

filename = 'fMRI_patientConc';


if(F_TASK == F_RECOMPUTE)  
   for npat = 1:21
      data = zX_RS(:,:,1:npat); data = data(:,:);
      [n, T, ~] = size(data);
      disp(['patient = 1:', num2str(npat)])
      disp(['doing original TS']);

      % prepare data
      fw_data = data';
      bw_data = flipud(fw_data);
      % fit AR process 
      [~,B,~,res_fw] = CBIG_RL2017_ar_mls(fw_data, p);
      [~,~,~,res_bw] = CBIG_RL2017_ar_mls(bw_data, p);
      % gaussianity measure
      [~, ~, deltaKurt_fw] = mardiaKurtosis(res_fw');
      [~, ~, deltaKurt_bw] = mardiaKurtosis(res_bw');

      delta_bwfw = deltaKurt_bw - deltaKurt_fw;
      results{npat}.ref = delta_bwfw;

      % store reference info in structure
      ref.X0 = B(:,1);
      phi_id = zeros(n,n,p);
      for pp = 1:p
         phi_id(:,:,pp) = B(:,2+n*(pp-1):n*(pp-1)+n+1);
      end
      ref.phi = phi_id;
      ref.delta_bwfw = delta_bwfw;
      ref.res = res_fw;         

      %% Create nullset
      disp(['doing nullset']);
      [sign,nullset] = genNullset('mvar', ref, p, T, n, rep);

      results{npat}.nullset = nullset;
      results{npat}.sign = sign;
      results{npat}.rep = rep;

      save(['results/', filename, '.mat'], 'results');
   end
else
   % load last results 
   load(['results/', filename, '.mat'])
end

%% Plot
fig = figure; hold on;
for npat = 1:21 
   if(results{npat}.sign < 0.95)
      L(1) = stem(npat, results{npat}.ref, 'r.');
   else
      L(2) = stem(npat, results{npat}.ref, 'b.'); 
   end
   lbs{npat} = npat;
end
grid minor; hold off; axis tight;
xticks(1:npat)
xticklabels(lbs)
xlabel('patient considered (concatenated)'); ylabel('k4_{bw-fw}'); 
legend(L, {'non-significant', 'significant'}, 'Location', 'SouthEast');

%% Print
print(fig, '-depsc2', ['images/', filename])

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end