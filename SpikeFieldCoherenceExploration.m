%% explore whether there is spike-field coherence at all, if so which freqs
% which phases of freqs, etc...

% load data
clear, clc

datadir = '/Volumes/data/RatBMI/data';
animal = 'S34';
session = '303';
datafile = fullfile(datadir,sprintf('data_block_%s_%s',animal,session));

fid = load(datafile,'data','Fs_lfp','TimeStamps');
binSize = 4;

%% sleep vs. awake
sleep_idx = get_sleep_epochs(fid.data',fid.Fs_lfp,binSize);

figure, hold on
for ch=1:32,
    leg{ch} = num2str(ch);
    
    %% collect lfp for each channel
    lfp = fid.data(ch,:)';
    samples = round(fid.Fs_lfp * binSize);
    N = length(lfp) - mod(length(lfp),samples);
    rows = samples;
    cols = N/rows;
    lfp = reshape(lfp(1:N),rows,cols);
    sleeplfp = lfp(:, sleep_idx);
    awakelfp = lfp(:,~sleep_idx);

    %% collect spike times from all units

    % all spikes
    spikeTimes = [];
    for j=2:3,
        spikeTimes = cat(1, spikeTimes, fid.TimeStamps{ch,j}');
    end

    % separate spikes by sleep vs. awake
    sleepSpikes = struct('spikeTimes',[]);
    awakeSpikes = struct('spikeTimes',[]);
    s_cnt = 1;
    a_cnt = 1;
    for i=1:length(sleep_idx),
        t0 = (i-1)*binSize;
        tf = i*binSize;
        idx = spikeTimes>=t0 & spikeTimes<tf;
        if sleep_idx(i),
            sleepSpikes(s_cnt).spikeTimes = spikeTimes(idx) - t0;
            s_cnt = s_cnt + 1;
        else
            awakeSpikes(a_cnt).spikeTimes = spikeTimes(idx) - t0;
            a_cnt = a_cnt + 1;
        end
    end

    %% calc spike-field coherence for awake vs. asleep
    T = binSize;
    W = 1;
    p = 1;
    params.tapers = [T,W,p];
    params.pad = -1;
    params.Fs = fid.Fs_lfp;
    params.fpass = [1,60];
    [C,phi,S12,S1,S2,f,zerosp]=coherencycpt(sleeplfp,sleepSpikes,params,[],[]);
%     [C,phi,S12,S1,S2,f,zerosp]=coherencycpt(awakelfp,awakeSpikes,params,[],[]);
    plot(f,mean(C(:,~zerosp),2))
end
legend(leg)
