function [rta_early,rta_late] = reward_triggered_average(timeStamps,channels,binSize,win,events)
%% reward_triggered_average(timeStamps,channels,binSize,window)
%   Calculates the average firing rate on specified channels aligned to the
%   time of reward. Also creates a plot.
% 
%   timeStamps should be a cell array of spike times {channels x units}
%   channels is a vector of channels (2 x 1 or 1 x 2)
%   binSize is a scalar value of bin width for measuring rates (sec)
%   win is a (1x2) vector of time (secs) to calc rates e.g., [.1,.2]
%   calculates rates 100ms before to 200ms after the reward times
%   events vector of event times (secs) to trigger on
% 
% rta = reward_triggered_average(timeStamps,channels,binSize,epoch_idx)
%   returns firing rates for each channe

%% inputs
narginchk(1,5)
if ~iscell(timeStamps),
    error('timeStamps should be a cell array')
elseif ~isvector(channels),
    error('channels should be a vector of 2 channels')
elseif ~isscalar(binSize),
    error('binSize should be a scalar')
elseif ~isvector(win) || length(win)~=2,
    error('window should be a (1x2) vector')
elseif ~isvector(events),
    error('events should be a vector of event times (secs)')
end

%% use gaussian window to smooth
w = gausswin(15);
w = w / sum(w);

%% collect number of spikes on each channel
bins = (-1*win(1):binSize:win(2));
nBins = length(bins)-1;
E = length(events);
fifth = 20;round(E/5);

rta_early = zeros(length(channels),nBins);
rta_late = zeros(length(channels),nBins);
for i=1:length(channels),
    spike_mat = zeros(E,nBins);
    ch = channels(i);
    spikeTimes = [];
    for j=2:3,
        spikeTimes = cat(1, spikeTimes, timeStamps{ch,j}');
    end
    for e=1:E,
        [spike_mat(e,:),~] = histcounts(spikeTimes,bins+events(e));
    end
    rta_early(i,:) = mean(spike_mat(1:fifth,:)) / binSize;
    rta_early(i,:) = conv(rta_early(i,:),w,'same');
    rta_late(i,:) = mean(spike_mat(end-fifth+1:end,:)) / binSize;
    rta_late(i,:) = conv(rta_late(i,:),w,'same');
end

%% plot
cc = get(groot,'defaultAxesColorOrder');

hold on
for i=1:length(channels),
    plot(bins(2:end),rta_early(i,:),'--','Color',cc(i,:))
    plot(bins(2:end),rta_late(i,:),'-','Color',cc(i,:))
    vline(0,'-k')
end
title('reward triggered average')
xlabel('time wrt reward (secs)')
ylabel('firing rate')
legend('early','late','Location','NorthWest')


