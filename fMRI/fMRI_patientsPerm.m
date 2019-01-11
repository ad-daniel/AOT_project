clc; clear all; close all;

addpath(genpath('../support/'));
addpath('dataset');
settings;

%% Load rest state
load X_RS_21subjs.mat;
p = 1;

rep = 100;     % nullset dimension
npat_max = 21; % only pick up to npat_max random concatenated patients
do_times = 3;  % how many different combinations to use 

%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

filename = 'fMRI_patientsPerm';


if(F_TASK == F_RECOMPUTE)  
   for npat = 1:npat_max % taking up to npat_max random patients
      for dt = 1:do_times
         perm = randperm(21,npat); % get npat patients from the 21
         % concatenate data from patients in the permutation
         data = [];
         for j = 1:length(perm)
            data = cat(2, data, zX_RS(:,:,perm(j):perm(j)));
         end

         [n, T, ~] = size(data);
         disp(['patients considered = ', num2str(npat), ' do times: ', num2str(dt)])
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
         results{npat}.attempt{dt}.ref = delta_bwfw;

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
         
         results{npat}.attempt{dt}.perm = perm;
         results{npat}.attempt{dt}.nullset = nullset;
         results{npat}.attempt{dt}.sign = sign;
         results{npat}.attempt{dt}.rep = rep;

         save(['results/', filename, '.mat'], 'results');
      end
   end
else
   % load last results 
   load(['results/', filename, '.mat'])
end

%% Plot
set(0,'defaultAxesFontSize',12)

plot_attempt = 3;
modif = ['v', num2str(plot_attempt)];

fig = figure; hold on;
for npat = 1:npat_max 
   if(results{npat}.attempt{plot_attempt}.sign < 0.95)
      L(1) = stem(npat, results{npat}.attempt{plot_attempt}.ref, 'r.');
   else
      L(2) = stem(npat, results{npat}.attempt{plot_attempt}.ref, 'b.'); 
   end
   lbs{npat} = mat2str(results{npat}.attempt{plot_attempt}.perm);
end
grid minor; hold off; axis tight;
xticks(1:npat_max)
ylabel('k4_{bw-fw}');xlabel('# patients considered in the permutation')
legend(L, {'non-significant', 'significant'}, 'Location', 'SouthEast');


%% Print
print(fig, '-depsc2', ['images/', filename, modif])

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end