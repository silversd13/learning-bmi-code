function varargout = neuron_stats(timeStamps,channels,binSize,epoch_idx)
%% neuron_stats(timeStamps,channels,binSize,epoch_idx)
%   Calculates the spiking statistics of specified channels. Also creates a
%   plot of spiking stats for the first two channels.
% 
%   timeStamps should be a cell array of spike times {channels x units}
%   channels is a vector of channels (2 x 1 or 1 x 2)
%   binSize is a scalar value corresponding to the time of an epoch in secs
%   epoch_idx is a logical vector that limits correlation measure to
%       specified bins, 
%       alternatively, you can input a scalar which is the time of the
%       lfp signal
% 
% FR = neuron_stats(timeStamps,channels,binSize,epoch_idx)
%   returns firing rates for each channe
% 
% [FR,FF] = neuron_stats(timeStamps,channels,binSize,epoch_idx)
%   returns firing rates and fano factors for each channels
% 
% [FR,FF,R] = neuron_stats(timeStamps,channels,binSize,epoch_idx)
%   also returns correlation matrix

%% inputs
narginchk(1,5)
if ~iscell(timeStamps),
    error('timeStamps should be a cell array')
elseif ~isvector(channels),
    error('channels should be a vector of 2 channels')
elseif ~isscalar(binSize),
    error('binSize should be a scalar')
elseif isscalar(epoch_idx),
elseif ~islogical(epoch_idx) || ~isvector(epoch_idx),
    error('epoch_idx should be a logical vector of indices')
end

%% collect number of spikes on each channel
if ~isscalar(epoch_idx),
    bins = (0:binSize:length(epoch_idx)*binSize)';
else % scalar (epoch_idx is time of lfp signal)
    bins = (0:binSize:epoch_idx)';
end

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
if ~isscalar(epoch_idx),
    nSpikes = nSpikes(epoch_idx,:);
end

%% calc statistics
FR = mean(nSpikes);
FF = var(nSpikes) ./ mean(nSpikes);
R = corr(nSpikes);

%% plot
if length(channels) > 1,
    scatter(nSpikes(:,1),nSpikes(:,2)), hold on
    text(.8*max(nSpikes(:,1)),.8*max(nSpikes(:,2)),...
        sprintf('r = %.3f',R(1,2)),'FontSize',18)
    xlabel(sprintf('Spikes (ch%02i)',channels(1)))
    ylabel(sprintf('Spikes (ch%02i)',channels(2)))
end

%% output
varargout{1} = FR;
if nargout>=2,
    varargout{2} = FF;
    if nargout==3,
        varargout{3} = R;
    end
end

