function [trial_rates,r] = spike_correlation_by_trial(timeStamps,channels,win,events)
%% spike_correlation_by_trial(timeStamps,channels,binSize,window)
%   Calculates the spiking on 2 channels in a window (defined by the user)
%   with respect to the timing in the input events vector. This creates two
%   timeseries vectors (one for each channel). This function plots those
%   timeseries and returns a sliding average of their correlation  from
%   early in training to late in learning 
% 
%   timeStamps should be a cell array of spike times {channels x units}
%   channels is a vector of channels (2 x 1 or 1 x 2)
%   win is a (1x2) vector of time (secs) to calc rates e.g., [.1,.2]
%   calculates rates 100ms before to 200ms after the reward times
%   events vector of event times (secs) to trigger on
% 
% [trial_rates,r] = spike_correlation_by_trial(timeStamps,channels,binSize,epoch_idx)
%   returns timeseries for each channel (trial x channel)
%   returns early and late correlation btw channels

%% inputs
narginchk(1,4)
if ~iscell(timeStamps),
    error('timeStamps should be a cell array')
elseif ~isvector(channels) || length(channels)~=2,
    error('channels should be a vector of 2 channels')
elseif ~isvector(win) || length(win)~=2,
    error('window should be a (1x2) vector')
elseif ~isvector(events),
    error('events should be a vector of event times (secs)')
end

%% define gaussian window for smoothing
w = gausswin(15);
w = w / sum(w);

%% collect number of spikes on each channel
E = length(events);
T = 20;
W = win(1) + win(2);

trial_rates = zeros(E,length(channels));
for i=1:length(channels),
    spike_mat = zeros(E,1);
    ch = channels(i);
    spikeTimes = [];
    for j=2:3,
        spikeTimes = cat(1, spikeTimes, timeStamps{ch,j}');
    end
    for e=1:E,
        spike_mat(e,:) = sum(spikeTimes>=events(e)-win(1) & spikeTimes<events(e)+win(2));
    end
    trial_rates(:,i) = spike_mat' / W;
end

%% correlation
tstart = [1,E-T];
r = zeros(length(tstart),1);
for i=1:length(tstart),
    t = tstart(i);
    idx = t:t+T;
    tmp = corr(trial_rates(idx,:));
    r(i) = tmp(1,2);
end

%% plot
hold on
for i=1:length(tstart),
    t = tstart(i);
    idx = t:t+T;
    jitterx = 0 + .25*randn(size(trial_rates(idx,1)));
    jittery = 0 + .25*randn(size(trial_rates(idx,1)));
    scatter(trial_rates(idx,1)+jitterx,trial_rates(idx,2)+jittery)
end
title('trial x trial correlation')
xlabel(sprintf('Spikes (ch%i)',channels(1)))
ylabel(sprintf('Spikes (ch%i)',channels(2)))
legend(sprintf('early: r=%.02f',r(1)),sprintf('late: r=%.02f',r(2)))


