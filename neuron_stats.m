function varargout = neuron_stats(timeStamps,channels)
%% neuron_stats(timeStamps,channels)
%   Calculates the first order spiking statistics on specified channels. 
% 
%   timeStamps should be a cell array of spike times {channels x units}
%   channels is a vector of channels
% 
% FR = neuron_stats(timeStamps,channels)
%   returns firing rates for each channel
% 
% [FR,CV] = neuron_stats(timeStamps,channels)
%   also returns coefficient of variation for each channel
% 
% [FR,FF,ISI] = neuron_stats(timeStamps,channels)
%   also returns interspike intervals for each channel

%% inputs
narginchk(1,2)
assert(iscell(timeStamps),'timeStamps should be a cell array');
assert(isvector(channels),'channels should be a vector of channels');

%% collect spike times on each channel (each ch may have multiple units)
spikeTimes = cell(length(channels),1);
for i=1:length(channels),
    ch = channels(i);
    for j=2:3,
        idx = ~isnan(timeStamps{ch,j});
        spikeTimes{i} = cat(1, spikeTimes{i}, timeStamps{ch,j}(idx)');
    end
end

%% spike train statistics
for i=1:length(channels),
    ISI{i} = diff(spikeTimes{i}); % distribution
    FR{i} = 1 / mean(ISI{i});
    CV{i} = mean(ISI{i}) / std(ISI{i});
end

%% output
varargout{1} = FR;
if nargout>1,
    varargout{2} = CV;
    if nargout>2,
        varargout{3} = ISI;
    end
end

