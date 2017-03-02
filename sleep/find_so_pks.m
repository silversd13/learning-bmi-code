function varargout = find_so_pks(lfp,Fs,win,sleep_idx)
%% find_so_pks(lfp,Fs,window,sleep_idx)
%   finds peaks of slow oscillations in delta band 
%       lfp should be a vector or matrix (samples x channels)
%       Fs is the sample rate of lfp (default is 24414.0625 / 24)
%       win is the size of the epochs in secs (default is 4 secs)
%       sleep_idx is a logical vector (1 - sleep, 0 - awake)
% 
% [pks,prv,nxt] = find_so_pks(lfp,...)
%   returns time of pks
%   returns time of previous troughs (for each peak)
%   returns time of next troughs (for each peak)

%% deal with inputs
narginchk(1,4)
assert(ismatrix(lfp),'lfp should be a vector or matrix (samples x channels)')
assert(isscalar(Fs),'Fs should be a scalar')
assert(isscalar(win),'win should be a scalar')
assert(islogical(sleep_idx) && isvector(sleep_idx),...
    'sleep_idx should be a logical vector')

%% sizing info
samples = round(Fs * win);
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
pks_idx = cell(size(lfp_delta,2),1);
prv_idx = cell(size(lfp_delta,2),1);
nxt_idx = cell(size(lfp_delta,2),1);
for ch=1:size(lfp_delta,2),
    % limit analysis to times when rat is asleep
    delta = reshape(lfp_delta(1:N,ch),rows,cols)';
    sleep = delta(sleep_idx,:)';
    sleep = sleep(:);
    
    % set threshold for oscillations
    thresh = prctile(sleep,90);

    % find peaks on each channel
    [~,pks_idx{ch}] = findpeaks(sleep,'MinPeakHeight',thresh);
    
    % find troughs in either direction
    nxt_idx{ch} = nan(size(pks_idx{ch}));
    prv_idx{ch} = nan(size(pks_idx{ch}));
    for i=1:length(pks_idx{ch}),
        idx = 2:pks_idx{ch}(i);
        prv_idx{ch}(i) = max([1,idx(1)+find(sleep(idx-1)>sleep(idx),1,'last')+1]);
        idx = pks_idx{ch}(i):length(sleep)-1;
        nxt_idx{ch}(i) = min([idx(1)+find(sleep(idx+1)>sleep(idx),1)-1,length(sleep)]);
    end
end
%% output
if nargout>0,
    varargout{1} = time(pks_idx);
    if nargout>1,
        varargout{2} = time(prv_idx);
        if nargout>2,
            varargout{3} = time(nxt_idx);
        end
    end
end

