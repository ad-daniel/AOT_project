clc; clear all; close all;

addpath(genpath('../support/'));
settings;

res_mvar = load('results/fMRI_mvarROI'); res_mvar = res_mvar.results;
res_uvar = load('results/fMRI_uvarROI'); res_uvar = res_uvar.results;

%% plot results from each method
set(0,'defaultAxesFontSize',15)

fig = figure;
stem(res_uvar.RS.deltaKurt, 'b.');
axis tight; grid minor;
xlabel('ROI'), ylabel('k4_{bw-fw}')
print(fig, '-depsc2', 'images/fMRI_ROI_prop1')

fig = figure;
stem(res_mvar.RS.ROI_deltaKurt, 'b.');
axis tight; grid minor;
xlabel('ROI'), ylabel('k4_{bw-fw}')
print(fig, '-depsc2', 'images/fMRI_ROI_prop2')


%% plot crosscorrelation between methods
set(0,'defaultAxesFontSize',15)

fig = figure;
crosscorr(res_uvar.RS.deltaKurt,res_mvar.RS.ROI_deltaKurt, 20);
title(''); ylabel('cross correlation')
print(fig, '-depsc2', 'images/fMRI_ROI_corr')

%% Plot significant results (uvar method only)
% Adapt pval to repeated tests
ROI = 1:360;
pStar = 0.05;

for i = 1:2
   if(i==1)
      imagename = 'images/RS_ROIsign';
      ROI_fw_listname = 'results/RS_fw_list.mat'; 
      ROI_bw_listname = 'results/RS_bw_list.mat';
      
      deltaKurt = res_uvar.RS.deltaKurt;  
      p = 1 - res_uvar.RS.sign;           % express significance to pvalue
      [pSorted,idx] = sort(p, 'ascend');
   else
      imagename = 'images/TS_ROIsign';
      ROI_fw_listname = 'results/TS_fw_list.mat';
      ROI_bw_listname = 'results/TS_bw_list.mat';
 
      deltaKurt = res_uvar.TS.deltaKurt;
      p = 1 - res_uvar.TS.sign;           % express significance as pvalue
      [pSorted,idx] = sort(p, 'ascend');      
   end
   
   ROI_fw_list = [];
   ROI_bw_list = [];
   
   fig = figure; hold on; grid minor;
   for nn = 1:360      
      if( pSorted(nn) <= (pStar * nn) / 360 )
         L(1) = stem(idx(nn), deltaKurt(idx(nn)), 'b.');
         if(deltaKurt(idx(nn)) < 0)
            ROI_fw_list = [ROI_fw_list, idx(nn)];
         else
            ROI_bw_list = [ROI_bw_list, idx(nn)];
         end

      else
         L(2) = stem(idx(nn), deltaKurt(idx(nn)), 'r.');
      end
   end
   hold off; axis tight;
   xlabel('ROI considered'); ylabel('k4_{bw-fw}'); 
   legend(L, {'significant', 'non-significant'}, 'Location', 'NorthWest');
   
   print(fig, '-depsc2', imagename)
   save(ROI_fw_listname, 'ROI_fw_list');
   save(ROI_bw_listname, 'ROI_bw_list');
end
