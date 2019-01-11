clear all; close all; clc;


%% Include
addpath(genpath('../support/'));
settings;

% dataset choice
DATASET_ID = 1;

switch(DATASET_ID)
    case 1
        % Monthly numbers of new fish hatching in a region of the Pacific Ocean
        dataset = textread('datasets/Recruit.tsm', '%f');
        sTitle = 'Monthy number of fish hatchings in the Pacific';
        filename = 'images/uvar_fish_causality';
    case 2
        % NIKKEI avg
        % Time Period: 01 jan 1992 - 08 oct 2018, Frequency:Daily
        % Nikkei 255, OPEN price
        load('datasets/nikkei_255.mat');
        sTitle = 'Nikkei 255 open price daily';  
        dataset = nikkei_255(end-500:end);
        filename = 'images/uvar_nikkei_causality';
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

%% Test causality
%dataset = randn(300,1);
dset_shuffled = data.y;

for i = 1:2000
    dset_shuffled = dset_shuffled(randperm(length(dset_shuffled)));
end

% 95% confidence level
L = length(data.y);
vcrit = sqrt(2)*erfinv(0.95);
low_conf = -vcrit/sqrt(L); up_conf = vcrit/sqrt(L);

fig = figure;
set(0,'defaultAxesFontSize',15)

stem(autocorr(dataset), 'r.'); hold on; grid minor;
stem(autocorr(dset_shuffled), 'b.'); hold off; axis tight;
hline = refline([0 low_conf]); hline.Color = 'g';
hline = refline([0 up_conf]); hline.Color = 'g';
xlabel('Lag'); ylabel('autocorrelation')
legend('unshuffled', 'shuffled', '95% conf', 'Location', 'NorthEast');

%% Print
print(fig, '-depsc2', filename)