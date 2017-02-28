function varargout = get_sleep_epochs(lfp,varargin)
%% get_sleep_epochs(lfp)
%   lfp should be a vector or matrix (samples x channels)
%   produces a plot of the lfp data where epochs are classified as sleep or
%   awake
%   classification is done by looking at the avg pwr in delta and gamma
%   bands (via hilbert method) then clustered using K-Means alg
% 
% get_sleep_epochs(lfp,Fs)
%   Fs is the sample rate of lfp (default is 24414.0625 / 24)
% 
% get_sleep_epochs(lfp,Fs,win)
%   win is the size of the epochs in secs (default is 4 secs)
% 
% sleep_idx = get_sleep_epochs(lfp,...)
%   also returns indices 
%       1 = sleep
%       0 = awake

%% deal with inputs
narginchk(1,3)
if nargin==1,
    Fs = 24414.0625 / 24;
    window = 4;
elseif nargin==2,
    Fs = varargin{1};
    window = 4;
elseif nargin==3,
    Fs = varargin{1};
    window = varargin{2};
end

if ~ismatrix(lfp),
    error('lfp should be a vector or matrix (samples x channels)')
end
if ~isscalar(Fs),
    error('Fs should be a scaler')
end
if ~isscalar(window),
    error('win should be a scaler')
end

%% sizing info
samples = round(Fs * window);
N = size(lfp,1) - mod(size(lfp,1),samples);
rows = samples;
cols = N/rows;

%% zscore lfp so that channels are comparable
lfp = zscore(lfp);

%% get pwr in delta and gamma bands
% frequency bands
fpass = {
    [.5,4]  % delta
    [40,70] % gamma
    };

% filter for delta, hilbert, and avg across chans
[b,a] = butter(2,fpass{1}/(Fs/2));
lfp_delta = filtfilt(b,a,lfp);
Xr_delta = mean(abs(hilbert(lfp_delta)),2);

% filter for gamma, hilbert, and avg across chans
[b,a] = butter(2,fpass{2}/(Fs/2));
lfp_gamma = filtfilt(b,a,lfp);
Xr_gamma = mean(abs(hilbert(lfp_gamma)),2);

%% average pwr in each epoch
delta = mean(reshape(Xr_delta(1:N),rows,cols))';
gamma = mean(reshape(Xr_gamma(1:N),rows,cols))';

%% do K-means to get indices of sleep
obj = fitgmdist([delta,gamma],2);
ctr = obj.mu;
idx = cluster(obj,[delta,gamma]);

if ctr(1,1)>ctr(2,1), % higher delta power
    sleep_idx = idx==1;
else
    sleep_idx = idx==2;
end

%% plot
gscatter(delta,gamma,sleep_idx), hold on
legend('awake','sleep')

%% output
if nargout==1,
    varargout{1} = sleep_idx;
end


