%% sleep analysis
clear, clc

% Xu, Min, et al. "Basal forebrain circuit for sleep-wake control." 
% Nature neuroscience 18.11 (2015): 1641-1647.
% 
% Polysomnographic recordings and analysis.
% The signals were recorded and amplified using AM-Systems 
% amplifiers, filtered (0.1-1,000 Hz or 10-1,000 Hz for EEG and EMG 
% recordings, respectively) and digitized at 600 Hz using LabView. 
% Spectral analysis was carried out using fast Fourier transform (FFT)
% and NREM, REM and wake states were semi-automatically classified 
% using a sleep analysis software (SleepSign for Animal, Kissei Comtec
% America) for each 10-s epoch (wake: desynchronized EEG and high EMG 
% activity; NREM sleep: synchronized EEG with high power at 0.5-4 Hz 
% and low EMG activity; REM sleep: desynchronized EEG with high power 
% at theta frequencies (6-9 Hz) and low EMG activity).

%% load data
datadir = '/Volumes/data/RatBMI/data';
animal = 'S34';
% session = '302-308';
% datafile = fullfile(datadir,sprintf('data_block_%s_Cat%s',animal,session));
session = '302';
datafile = fullfile(datadir,sprintf('data_block_%s_%s',animal,session));

fid = load(datafile,'data','Fs_lfp');
lfp = fid.data(1,:);
Fs = fid.Fs_lfp;

%% manually get indices of SWS and not SWS / SWS = Slow Wave Sleep

plot(lfp)
title('Choose two timepoints that encapsulate a period of sleep')
[x,~] = ginput(2);
sleep_idx = round(x);
title('Choose two timepoints that encapsulate a period of awakeness')
[x,~] = ginput(2);
awake_idx = round(x);
close all

% make sure that they are both the same length
Nsleep = sleep_idx(2)-sleep_idx(1)+1;
Nawake = awake_idx(2)-awake_idx(1)+1;
if Nsleep < Nawake,
    awake_idx(2) = awake_idx(1) + Nsleep - 1;
else
    sleep_idx(2) = sleep_idx(1) + Nsleep - 1;
end

% assign LFP to sleep or awake
sleep = lfp(:,sleep_idx(1):sleep_idx(2));
awake = lfp(:,awake_idx(1):awake_idx(2));

%% spectral content of LFP

for ch=1:size(fid.data,1),
    % load  channels
    lfp = fid.data(ch,:);
    Fs = fid.Fs_lfp;
    
    % assign LFP to sleep or awake
    sleep = lfp(:,sleep_idx(1):sleep_idx(2));
    awake = lfp(:,awake_idx(1):awake_idx(2));

    % spectral content
    [Psleep,F] = pwelch(sleep,[],[],[],Fs);
    [Pawake,~] = pwelch(awake,[],[],[],Fs);

    % PLOT

    fpass = [.1,120];
    fidx = F>fpass(1) & F<fpass(2);

    clf;
    ax(1) = subplot(2,2,1);
    plot(sleep)
    title('sleep')

    ax(2) = subplot(2,2,2);
    plot(awake)
    title('awake')

    linkaxes(ax,'xy')

    subplot(2,2,[3,4]), hold on
    plot(F(fidx),smooth(10*log10(Psleep(fidx)),20))
    plot(F(fidx),smooth(10*log10(Pawake(fidx)),20))
    legend('sleep','awake')
    waitforbuttonpress;
end

%% power in entire lfp signal

for ch=1:size(fid.data,1),
    % load  channels
    lfp = fid.data(ch,:);
    Fs = fid.Fs_lfp;
    
    % assign LFP to sleep or awake
    sleep = lfp(:,sleep_idx(1):sleep_idx(2));
    awake = lfp(:,awake_idx(1):awake_idx(2));

    % filter in delta band
    [b,a] = butter(2,[.5,4]/(Fs/2));
    sleep_d = filtfilt(b,a,sleep);
    awake_d = filtfilt(b,a,awake);
    
    % filter in gamma band
    [b,a] = butter(2,[40,70]/(Fs/2));
    sleep_g = filtfilt(b,a,sleep);
    awake_g = filtfilt(b,a,awake);
    
    % power 
    [Psleep] = abs(hilbert(sleep_d)) ./ abs(hilbert(sleep_g));
    [Pawake] = abs(hilbert(awake_d)) ./ abs(hilbert(awake_g));

    % PLOT

    fpass = [.1,120];
    fidx = F>fpass(1) & F<fpass(2);

    clf;
    ax(1) = subplot(2,2,1);
    plot(sleep)
    title('sleep')

    ax(2) = subplot(2,2,2);
    plot(awake)
    title('awake')

    linkaxes(ax,'xy')

    subplot(2,2,[3,4]), hold on
    plot(smooth(Psleep))
    plot(smooth(Pawake))
    legend('sleep','awake')
    waitforbuttonpress;
end

%% in 10 sec windows, get pwr in different frequency bands

% load data for all channels
fid = load(datafile,'data','Fs_lfp');
lfp = fid.data;
Fs = fid.Fs_lfp;
sleep = lfp(:,sleep_idx(1):sleep_idx(2))';
awake = lfp(:,awake_idx(1):awake_idx(2))';
clear fid lfp

% frequency bands
fpass = {
    [.5,5]  % delta
    [5,12]   % theta
    [12,20] % beta
    [40,60] % gamma
    [60,250]% high gamma
    [.5 20] % low freq
    [40,250]% high freq
    };
freqs = {
    'delta'
    'theta'
    'beta'
    'gamma'
    'high gamma'
    'low freqs'
    'high freqs'
    };

window = 4;
samples = round(Fs * window);
Ps = [];
Pa = [];
for i=1:samples:size(sleep,1)-samples,
    idx = i:i+samples-1;
    [Pss,F] = pwelch(sleep(idx,:),[],[],[],Fs);
    [Paa,~] = pwelch(awake(idx,:),[],[],[],Fs);
    
    Ps(end+1,:) = zeros(1,length(fpass));
    Pa(end+1,:) = zeros(1,length(fpass));
    for fp=1:length(fpass),
        fidx = F>=fpass{fp}(1) & F<fpass{fp}(2);
        Ps(end,fp) = mean(mean(10*log10(Pss(fidx,:))));
        Pa(end,fp) = mean(mean(10*log10(Paa(fidx,:))));
    end
end

%% visualize each freq band
figure
for i=1:length(fpass),
    subplot(length(fpass),2,i), hold on
    histogram(Ps(:,i),'BinWidth',1,'Normalization','pdf','DisplayStyle','stairs')
    histogram(Pa(:,i),'BinWidth',1,'Normalization','pdf','DisplayStyle','stairs')
    title(freqs{i})
end

subplot(length(fpass),2,i+1:length(fpass)*2), hold on
scatter(Ps(:,end-1),Ps(:,end))
scatter(Pa(:,end-1),Pa(:,end))
legend('sleep','awake')

%% load all data and go through it in the same way
% load data for all channels
fid = load(datafile,'data','Fs_lfp');
lfp = zscore(fid.data');
Fs = fid.Fs_lfp;
clear fid

% frequency bands
fpass = {
    [.5,5]  % delta
    [40,60] % gamma
    };
freqs = {
    'delta'
    'gamma'
    };

window = 4;
samples = round(Fs * window);
P = [];
for i=1:samples:size(lfp,1)-samples,
    idx = i:i+samples-1;
    [Pxx,F] = pwelch(lfp(idx,:),[],[],[],Fs);
    
    P(end+1,:) = zeros(1,length(fpass));
    for fp=1:length(fpass),
        fidx = F>=fpass{fp}(1) & F<fpass{fp}(2);
        P(end,fp) = mean(mean(10*log10(Pxx(fidx,:))));
    end
end

% do K-means to get indices of sleep
KMM = kmeans(P,2);

% plot
figure
gscatter(P(:,1),P(:,2),KMM)
title('pwelch method')



%% do the same thing, except calc everything with hilbert
% load data for all channels
fid = load(datafile,'data','Fs_lfp');
lfp = zscore(fid.data');
Fs = fid.Fs_lfp;
clear fid

% frequency bands
fpass = {
    [.5,4]  % delta
    [40,70] % gamma
    };
freqs = {
    'delta'
    'gamma'
    };

[b,a] = butter(2,fpass{1}/(Fs/2));
lfp_delta = filtfilt(b,a,lfp);
Xr_delta = mean(abs(hilbert(lfp_delta)),2);

[b,a] = butter(2,fpass{2}/(Fs/2));
lfp_gamma = filtfilt(b,a,lfp);
Xr_gamma = mean(abs(hilbert(lfp_gamma)),2);
% Xr = mean( Xr_delta ./ Xr_gamma , 2);


% smoothing
window = 4; % secs
samples = round(Fs * window);
N = length(Xr_delta) - mod(length(Xr_delta),samples);
rows = samples;
cols = N/rows;
Xr_delta_avg = mean(reshape(Xr_delta(1:N),rows,cols))';
Xr_gamma_avg = mean(reshape(Xr_gamma(1:N),rows,cols))';

% do K-means to get indices of sleep
KMM = kmeans([Xr_delta_avg,Xr_gamma_avg],2);

% plot
figure
gscatter(Xr_delta_avg,Xr_gamma_avg,KMM)
title('hilbert method')



