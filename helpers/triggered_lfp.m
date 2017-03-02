function [waves,time] = triggered_lfp(lfp,Fs,events,win,varargin)
%% triggered_lfp(lfp,Fs,events,win)
%   collects waves triggered on times in events vector
%       lfp should be a vector or matrix (samples x channels)
%       Fs is the sample rate of lfp (default is 24414.0625 / 24)
%       events is a cell array of event times for each channel
%       win is a vector of length 2 (time window around events)
% triggered_lfp(lfp,Fs,events,win,fpass)
%       fpass is vector of freqs to filter at, default is [.1,200]
% 
% [waves,time] = triggered_lfp(lfp,...)
%   returns cell array of waves (events x time) for each channel
%   returns time vector

%% deal with inputs
narginchk(1,5)
assert(ismatrix(lfp),'lfp should be a vector or matrix (samples x channels)')
assert(isscalar(Fs),'Fs should be a scalar')
assert(iscell(events),'events should be a cell array')
assert(isvector(win) && length(win)==2,'win should be a scalar')
if nargin==4,
    fpass = varargin{1};
    assert(isvector(fpass) && length(fpass)==2,...
        'fpass should be a vector of filter freqs')
else
    fpass = [.1,200];
end
    
%% filter 
[b,a] = butter(2,fpass/(Fs/2));
lfp = filtfilt(b,a,lfp);

%% go through each channel and collect waves triggered on events
waves = cell(size(lfp,2),1);
for ch=1:size(lfp,2),
    waves{ch} = createdatamatc(lfp(:,ch),events{ch},Fs,win);
end

%% time vector
time = linspace(-win(1),win(2),size(waves{1},1))';

