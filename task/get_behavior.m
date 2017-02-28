function varargout = get_behavior(wave,varargin)
%% get_behavior(wave)
%   wave should be a vector of pulses that correspond to task times
% 
% get_behavior(wave,Fs)
%   sampling rate for wave, default is Fs = 24414.0625 / 24
% 
% time_start = get_behavior(wave,...)
%   returns the start times for each trial
% 
% [time_start,time_reward] = get_behavior(wave,...)
%   also returns reward times

%% inputs
narginchk(1,2)
if nargin==1,
    Fs = 24414.0625 / 24;
else
    Fs = varargin{1};
end
if ~isvector(wave),
    error('Input should be a vector of real values')    
end
if ~isscalar(Fs),
    error('Fs should be a scalar')    
end

%% grab leading edge of each pulse
time = (0:length(wave)-1) / Fs;
i = 1:length(wave)-1;
pulse_times = time(wave(i)<1 & wave(i+1)>=1); 
pulse_times = [pulse_times,inf];

%% separate into single and double pulses
single_pulse = [];
for i=1:length(pulse_times)-1,
    if (pulse_times(i+1)-pulse_times(i)) > .5,
        single_pulse(end+1) = pulse_times(i);
    end
end

%% if there are any unsuccessful trials mark reward as NaN
trial_start = [];
trial_reward = [];
for i=1:length(single_pulse),
    pulse = single_pulse(i);
    trial_start(i,1) = pulse;
    % find next pulse time
    next_pulse = pulse_times(find(pulse_times>pulse,1));
    % is is a single pulse or a double pulse
    if any(next_pulse==single_pulse),
        trial_reward(end+1,1) = nan;
    elseif isinf(next_pulse),
        trial_reward(end+1,1) = nan;
    else
        trial_reward(end+1,1) = next_pulse;
    end
end

%% output
if nargout>=1,
    varargout{1} = trial_start;
    if nargout==2,
        varargout{2} = trial_reward;
    end
end
