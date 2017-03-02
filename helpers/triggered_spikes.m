function [raster,time] = triggered_spikes(timeStamps,events,win,binSize)
%% [spks,rates]=triggered_spikes(timeStamps,events,win)
%   collects waves triggered on times in events vector
%       timeStamps should be a cell array of spike times {ch x units}
%       events is a cell array of event times for each channel
%       win is a vector of length 2 (time window around events)
%       binSize is a scalar value of bin width for measuring rates (sec)
%
%   outputs a cell array, where each cell is a cell! and has the following
%   fields
%       unit - cell
%       ch - channel of unit
%       spks - a matrix of binned spikes (events x time) for each unit

%% deal with inputs
narginchk(1,4)
assert(ismatrix(timeStamps),'timeStamps should be a a cell array {ch x units}')
assert(iscell(events),'events should be a cell array')
assert(isvector(win) && length(win)==2,'win should be a vector of length 2')
assert(isscalar(binSize),'binSize should be a scalar')

%% collect number of spikes for each unit in bins defined relative to events
bins = (-1*win(1):binSize:win(2));
nBins = length(bins)-1;
time = bins(1:nBins) + binSize;

unit = 0;
raster = cell(0);
for ch=1:length(events), % each channel
    E = length(events{ch});
    for j=2:3, % each unit
        spikeTimes = timeStamps{ch,j};
        spikeTimes = spikeTimes(~isnan(spikeTimes));
        if length(spikeTimes)<10, % skip if empty
            continue,
        else
            spks = zeros(E,nBins);
            for e=1:E, % bin around each event
                [spks(e,:),~] = histcounts(spikeTimes,bins+events{ch}(e));
            end
            % output
            unit = unit + 1;
            raster{unit}.unit = unit;
            raster{unit}.ch = ch;
            raster{unit}.spks = spks;
        end
    end
end



