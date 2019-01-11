clc; clear all; close all;

addpath(genpath('../support/'));
addpath('dataset');
settings;

p = 1;
n = 360;

%F_TASK = F_RECOMPUTE;
F_TASK = F_RELOAD;

filename = 'fMRI_mvarROI';

if(F_TASK == F_RECOMPUTE)
   for i = 1:2
      ROI_deltaKurt = zeros(n,1);
      
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
      % fit AR process 
      [~,~,~,res_fw] = CBIG_RL2017_ar_mls(fw_data, p);
      [~,~,~,res_bw] = CBIG_RL2017_ar_mls(bw_data, p);

      % Check each ROI separatedly by using multivar residuals
      for nn = 1:n 
         disp(['Doing data i = ', num2str(i), ' ROI n = ', num2str(nn)])
         % compute delta for each separate ROI
         [~, ~, deltaKurt_fw] = mardiaKurtosis(res_fw(nn,:)');
         [~, ~, deltaKurt_bw] = mardiaKurtosis(res_bw(nn,:)');
         
         ROI_deltaKurt(nn) = deltaKurt_bw - deltaKurt_fw;        
      end
      
      if(i == 1)
         results.RS.ROI_deltaKurt = ROI_deltaKurt;
      else
         results.TS.ROI_deltaKurt = ROI_deltaKurt;
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
stem(1:n,results.RS.ROI_deltaKurt, 'b.'); grid minor; axis tight;
xlabel('ROI'); ylabel('k4_{bw-fw}');


fig2 = figure;
stem(1:n,results.TS.ROI_deltaKurt, 'b.'); grid minor; axis tight;
xlabel('ROI'); ylabel('k4_{bw-fw}');


%% Print
print(fig1, '-depsc2', 'images/fMRI_uvarROI_RS')
print(fig2, '-depsc2', 'images/fMRI_uvarROI_TS')

[returnCode, hostName]=system('hostname');
if(strcmp(deblank(hostName),'miplabsrv3'))
   exit
end