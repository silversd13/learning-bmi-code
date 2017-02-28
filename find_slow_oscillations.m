function varargout = find_slow_oscillations(lfp,Fs,window,sleep_idx)
%% find_slow_oscillations(lfp)
%   lfp should be a vector or matrix (samples x channels)
%   produces a plot of the lfp data where epochs are classified as sleep or
%   awake
%   classification is done by looking at the avg pwr in delta and gamma
%   bands (via hilbert method) then clustered using K-Means alg
% 
% find_slow_oscillations(lfp,Fs)
%   Fs is the sample rate of lfp (default is 24414.0625 / 24)
% 
% find_slow_oscillations(lfp,Fs,window)
%   win is the size of the epochs in secs (default is 4 secs)
% 
% idx = find_slow_oscillations(lfp,...)
%   also returns indices 
%       1 = sleep
%       0 = awake

%% deal with inputs
narginchk(1,4)
assert(ismatrix(lfp),'lfp should be a vector or matrix (samples x channels)')
assert(isscalar(Fs),'Fs should be a scalar')
assert(isscalar(window),'window should be a scalar')
assert(islogical(sleep_idx) && isvector(sleep_idx),...
    'sleep_idx should be a logical vector')

%% sizing info
samples = round(Fs * window);
N = size(lfp,1) - mod(size(lfp,1),samples);
rows = samples;
cols = N/rows;

%% create time vector (important since concatenating)
dt = 1/Fs;
time = ((1:N)/Fs)' - dt;

% time during sleep
time = reshape(time(1:N),rows,cols)';
time = time(sleep_idx,:)';
time = time(:)';
    
%% filter in delta bands
% delta band
fpass = [.5,4];

% filter for delta
[b,a] = butter(2,fpass/(Fs/2));
lfp_delta = filtfilt(b,a,lfp);

%% go through each channel and find slow oscillations (pks in delta)
for ch=1:size(lfp_delta,2),
    % limit analysis to times when rat is asleep
    delta = reshape(lfp_delta(1:N,ch),rows,cols)';
    sleep = delta(sleep_idx,:)';
    sleep = sleep(:);
    

    % set threshold for oscillations
    thresh = prctile(sleep,90);

    % find peaks on each channel
    [~,idx{ch}] = findpeaks(sleep,time,'MinPeakHeight',thresh);
end
%% output
if nargout==1,
    varargout{1} = idx;
end

