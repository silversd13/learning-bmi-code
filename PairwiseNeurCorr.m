%% Calculate Correlations Btw 2 Neurons
clear, clc

%% load data
datadir = '/Volumes/data/RatBMI/data';
animal = 'S34';
session = '302-308';
datafile = fullfile(datadir,sprintf('data_block_%s_Cat%s',animal,session));

fid = load(datafile);

channels = [10,29];
gain = 1;

%% manually find time period of spontaneous activity 
% spontaneous = baseline, not sleep, and not task

% lfp
lfp = fid.data;

% time
dt = 1/fid.Fs_lfp;
time_lfp = (dt:dt:dt*size(lfp,2)) - dt;

% idx
idx0 = 14e5;
idxf = 19e5;

t0 = time_lfp(idx0);
tf = time_lfp(idxf);

% verify with LFP
strips(lfp(1:5:32,idx0:idxf)')

%% collect # of spikes for each neuron in 100ms bins and plot

width = 10;
bins = t0:width:tf;
spikes = [];
for i=1:length(channels),
    ch = channels(i);
    spikeTimes = [];
    for j=2:3,
        spikeTimes = [spikeTimes ...
            fid.TimeStamps{ch,j}(fid.TimeStamps{ch,j}>t0 & fid.TimeStamps{ch,j}<tf)];
    end
    [spikes(i,:),~] = histcounts(spikeTimes,bins);
end
R = corr(spikes');
r = R(1,2);

scatter(spikes(1,:)',spikes(2,:)'), hold on
title(sprintf('r = %.3f',r))



%% plot time to target across trials

% get first block
idx0 = find(fid.wave > 1, 1) - 1000;
idxf = 3.2e6;

[tstart,treward] = get_behavior(fid.wave(idx0:idxf));
timeToTarget = treward - tstart;

figure; hold on
for i=1:length(timeToTarget),
    plot(1:length(timeToTarget),timeToTarget,'v')
end

%% fit exponential learning curve to behavior
ls_opt = optimset('Display','off');
learningfun = @(beta,x) (beta(1) - beta(3))*exp(-1*x*beta(2)) + beta(3);
Beta = lsqcurvefit(learningfun,[5,1,2],1:length(timeToTarget),timeToTarget,[],[],ls_opt);

x = 1:length(timeToTarget);
yhat = learningfun(Beta,x);
plot(x,yhat,'--k')
title(sprintf('y = (%.2f)*exp{-%.2f*x} + %.2f',Beta(1)-Beta(3),Beta(2),Beta(3)))







