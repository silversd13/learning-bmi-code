function varargout = spike_spike_corr(timeStamps,channels,binSize,epoch_idx)
%% spike_spike_corr(timeStamps,channels,binSize,epoch_idx)
%   Calculates the neuron-to-neuron spiking correlations -> makes scatter
%   plot of rate1 vs. rate2 with correlation value
% 
%   timeStamps should be a cell array of spike times {channels x units}
%   channels is a vector of channels (2 x 1 or 1 x 2)
%   binSize is a scalar value corresponding to the time of an epoch in secs
%   epoch_idx is a logical vector that limits correlation measure to
%       specified bins
% 
% r = spike_spike_corr(timeStamps,channels,binSize,epoch_idx)
%   instead of making scatter plot, returns correlation value

%% inputs
narginchk(1,5)
if ~iscell(timeStamps),
    error('timeStamps should be a cell array')
elseif length(channels)~=2 || ~isvector(channels),
    error('channels should be a vector of 2 channels')
elseif ~isscalar(binSize),
    error('binSize should be a scalar')
elseif ~islogical(epoch_idx) || ~isvector(epoch_idx),
    error('epoch_idx should be a logical vector of indices')
end

%% outputs
nargoutchk(0,1)

%% collect number of spikes on each channel
bins = (0:binSize:length(epoch_idx)*binSize)';
nSpikes = [];
for i=1:length(channels),
    ch = channels(i);
    spikeTimes = [];
    for j=2:3,
        spikeTimes = cat(2, spikeTimes, timeStamps{ch,j});
    end
    [nSpikes(:,i),~] = histcounts(spikeTimes,bins);
end

%% select bins with epoch_idx
nSpikes = nSpikes(epoch_idx,:);

%% calc correlation
R = corr(nSpikes);
r = R(1,2);

%% plot/output
figure
scatter(nSpikes(:,1),nSpikes(:,2)), hold on
title(sprintf('r = %.3f',r))
xlabel(sprintf('Spikes (ch%02i)',channels(1)))
ylabel(sprintf('Spikes (ch%02i)',channels(2)))
if nargout==1,
    varargout{1} = r;
end
