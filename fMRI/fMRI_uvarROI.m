clc; clear all; close all;

addpath(genpath('../support/'));
addpath('dataset');
settings;

p = 40;
n = 360;
rep = 100;

%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

filename = 'fMRI_uvarROI';

if(F_TASK == F_RECOMPUTE)
   for i = 1:2
      ROI_deltaKurt = zeros(n,1);
      sign = zeros(n,1);
      
      if(i == 1)
         %% Load rest state
         load X_RS_21subjs.mat;
         data = zX_RS(:,:);
         [n, T] = size(data);
         % standardize
         for k = 1:n
            data(k,:) = (data(k,:) - mean(data(k,:)))/std(data(k,:));
         end
      else
         %% Load task state
         load X_TASK_21subjs.mat;
         data = zX_TASK(:,:);
         [n, T] = size(data);
         % standardize
         for k = 1:n
            data(k,:) = (data(k,:) - mean(data(k,:)))/std(data(k,:));
         end         
      end
      
      % prepare data
      fw_data = data'; % prepare as Txn
      bw_data = flipud(fw_data);

      % Check each ROI separatedly using multivar residuals
      for nn = 1:n 
         disp(['Doing data i = ', num2str(i), ' ROI n = ', num2str(nn)])
         % apply AoT univar to each ROI
         TS_fw = fw_data(:,nn);
         TS_bw = bw_data(:,nn); % already flipped
         
         [~,B,~,res_fw] = CBIG_RL2017_ar_mls(TS_fw, p);
         [~,~,~,res_bw] = CBIG_RL2017_ar_mls(TS_bw, p);
         
         [~, ~, deltaKurt_fw] = mardiaKurtosis(res_fw');
         [~, ~, deltaKurt_bw] = mardiaKurtosis(res_bw');
         
         ROI_deltaKurt(nn) = deltaKurt_bw - deltaKurt_fw;    
         
         % store reference info in structure
         ref.X0 = B(:,1);
         phi_id = zeros(1,1,p);
         for pp = 1:p
            phi_id(:,:,pp) = B(:,2+(pp-1):(pp-1)+2);
         end
         ref.phi = phi_id;
         ref.delta_bwfw = ROI_deltaKurt(nn);
         ref.res = res_fw;

         % generate nullset
         [sign(nn),~] = genNullset('uvar', ref, p, T, 1, rep);         
      end
      
      if(i == 1)
         results.RS.sign = sign;
         results.RS.deltaKurt = ROI_deltaKurt;
      else
         results.TS.sign = sign;
         results.TS.deltaKurt = ROI_deltaKurt;
      end
      
      save(['results/', filename, '.mat'], 'results');
   end
else
   % load last results 
   load(['results/', filename, '.mat'])
end

%% Plot
set(0,'defaultAxesFontSize',15)

fig1 = figure;
stem(1:n,results.RS.deltaKurt, 'b.'); grid minor; axis tight;
xlabel('ROI'); ylabel('k4_{bw-fw}');


fig2 = figure;
stem(1:n,results.TS.deltaKurt, 'b.'); grid minor; axis tight;
xlabel('ROI'); ylabel('k4_{bw-fw}');


%% Print
print(fig1, '-depsc2', 'images/fMRI_uvarROI_RS')
print(fig2, '-depsc2', 'images/fMRI_uvarROI_TS')

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end