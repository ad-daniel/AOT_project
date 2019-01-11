clear all; close all; clc;


%% Include
addpath(genpath('../support/'));
settings;

% dataset choice
DATASET_ID = 1;
n = 1;
rep = 10000;

switch(DATASET_ID)
    case 1
        % Monthly numbers of new fish hatching in a region of the Pacific Ocean
        dataset = textread('datasets/Recruit.tsm', '%f');
        sTitle = 'Monthy number of fish hatchings in the Pacific';
        refname = 'fish';
        pbest = 2;
    case 2
        % NIKKEI avg
        % Time Period: 01 jan 1992 - 08 oct 2018, Frequency:Daily
        % Nikkei 255, OPEN price
        load('datasets/nikkei_255.mat');
        sTitle = 'Nikkei 255 open price daily';  
        dataset = nikkei_255(end-500:end);
        refname = 'nikkei_long';
        pbest = 3;
end

% normalize
dataset = (dataset - mean(dataset))/std(dataset);
data = iddata(dataset, [], 1);

% detrend
if(DATASET_ID == 2)
    data = detrend(data, 1, 10);
else
    data = detrend(data);
end

%% AR model
for i = 1 : 6 % vary model order
   % prepare data
   fw_data = data.y;
   bw_data = flip(fw_data);
   T = length(fw_data);
   
   % fit AR process 
   [~,B,~,res_fw] = CBIG_RL2017_ar_mls(fw_data, i);
   [~,~,~,res_bw] = CBIG_RL2017_ar_mls(bw_data, i);
   % gaussianity measure
   [~, ~, deltaKurt_fw] = mardiaKurtosis(res_fw');
   [~, ~, deltaKurt_bw] = mardiaKurtosis(res_bw');    
   
   delta_bwfw(i) = deltaKurt_bw - deltaKurt_fw;
   
   % store reference info in structure
   ref.X0 = B(:,1);
   phi_id = zeros(n,n,i);
   for pp = 1:i
      phi_id(:,:,pp) = B(:,2+n*(pp-1):n*(pp-1)+n+1);
   end
   ref.phi = phi_id;
   ref.delta_bwfw = delta_bwfw(i);
   ref.res = res_fw;
   
   % generate nullset
   [sign(i),nullset{i}] = genNullset('uvar', ref, i, T, n, rep);
   
   % compute model quality
   opt = arOptions;
   opt.Window = 'now';
   mdl_fw{i} = ar(fw_data, i, opt);
   mdl_bw{i} = ar(bw_data, i, opt);
    
   AIC_metric(i) = aic(mdl_fw{i}, 'aic');
   BIC_metric(i) = aic(mdl_bw{i}, 'bic');
end

%% Plot data
fig = figure;
set(0,'defaultAxesFontSize',11)

subplot(2,2,1)
plot(data.y, 'b'); grid minor; axis tight;

subplot(2,2,2)
stem(nullset{pbest}, 'b.'); hold on; grid
line([0, rep], [delta_bwfw(pbest),delta_bwfw(pbest)], 'Color', 'r');
axis tight; grid minor; hold off;
xlabel('nullset'); ylabel('k4_{bw-fw}')
legend('nullset', 'original TS','Location', 'NorthEast', 'Orientation', 'vertical');

subplot(2,2,3)
plot(AIC_metric, 'r'); hold on;
plot(BIC_metric, 'b'); axis tight; grid minor; hold off;
xlabel('order'); ylabel('criterion'); 
legend('AIC', 'BIC', 'Location', 'NorthEast')

subplot(2,2,4); hold on;
for i = 1:6
   if(sign(i) > 0.95)
      L(1) = stem(i,delta_bwfw(i), 'b.'); 
   else
      L(2) = stem(i,delta_bwfw(i), 'r.'); 
   end
end
axis tight; grid minor; hold off;
legend(L, {'sign.', 'non-sign.'}, 'Location', 'NorthEast', 'Orientation', 'vertical')
xlabel('order'); ylabel('k4_{bw-fw}')

%suptitle(['for p = ', num2str(pbest), ' more significant than ', num2str(sign(pbest)*100), '%'])

%% Print
filename = ['images/uvar_realdata_', num2str(refname)];
print(fig, '-depsc2', filename)