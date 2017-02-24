%% Go through analysis for each session
clear, clc, close all
epochSize = 4; % seconds
datadir = '/Volumes/data/RatBMI/data';

experiments = get_experiments();
for experiment=11%:length(experiments),
    
    % make sure directories for saving exists
    savedir = fullfile('Users','daniel','Projects','LearningBMI','processed');
    if ~exist(savedir,'dir'),
        mkdir(savedir)
    end
    
    % session info
    rat = experiments(experiment).rat;
    channels = experiments(experiment).channels;
    name = sprintf('experiment_%02i',experiment);

    % output progress
    fprintf('\n%s: (expt %i of %i)',rat,experiment,length(experiments))
    clf, set(gcf,'Units','Normalized','Position',[0,0,1,1],...
        'Numbertitle','off','Name',name)
    
    %% sleep vs. awake
    % load sleep session to look at sleep vs awake neural activity
    session = experiments(experiment).sleep_session;
    if ~isempty(session),
        datafile = fullfile(datadir,sprintf('data_block_%s_%s',rat,session));
        fid = load(datafile,'data','Fs_lfp');

        subplot(3,3,1)
        sleep_idx = get_sleep_epochs(fid.data',fid.Fs_lfp,epochSize);

        % calc firing rates, fano factors, and correlation btw neurons
        fid = load(datafile,'TimeStamps');
        subplot(3,3,4)
        [FR_sleep,FF_sleep,R_sleep] = neuron_stats(fid.TimeStamps,channels,epochSize, sleep_idx);
        title('Sleep Correlation')
        subplot(3,3,7)
        [FR_awake,FF_awake,R_awake] = neuron_stats(fid.TimeStamps,channels,epochSize,~sleep_idx);
        title('Awake Correlation')
        if length(channels) > 1,
            r_sleep = R_sleep(1,2);
            r_awake = R_awake(1,2);
        else
            r_sleep = R_sleep;
            r_awake = R_awake;
        end
    else
        sleep_idx = [];
        r_sleep = [];
        r_awake = [];
        FR_sleep = [];
        FR_awake = [];
        FF_sleep = [];
        FF_awake = [];
    end
    
    %% baseline activity
    % load baseline session to look at spontaneous activity
    session = experiments(experiment).baseline_session;
    datafile = fullfile(datadir,sprintf('data_block_%s_%s',rat,session));
    fid = load(datafile,'data','Fs_lfp');
    N = round(size(fid.data,2) / fid.Fs_lfp);
    
    % calc firing rates, fano factors, and correlation btw neurons
    fid = load(datafile,'TimeStamps');
    subplot(3,3,2)
    [FR_base,FF_base,R_base] = neuron_stats(fid.TimeStamps,channels,epochSize, N);
    title('Baseline Correlation')
    if length(channels) > 1,
        r_base = R_base(1,2);
    else
        r_base = R_base;
    end

    % calc spike-field coherence during sleep & spontaneous activity
    % might add this later, but sfc seems pretty random across freqs, and
    % consistent across channels so it's not clear if it will be useful for
    % prediction

    %% task activity and learning
    % get behavior and calc learning rate
    session = experiments(experiment).task_session;
    datafile = fullfile(datadir,sprintf('data_block_%s_%s',rat,session));
    fid = load(datafile,'wave');
    [tstart,treward] = get_behavior(fid.wave);
    idx = ~isnan(treward);
    subplot(3,3,8:9)
    [beta,learning_fun] = calc_learning_rate(treward(idx)-tstart(idx));
    
    % use behavior and timestamps to get firing rates timed to reward
    fid = load(datafile,'TimeStamps');
    subplot(3,3,5)
    [rta_early{1},rta_late{1}] = reward_triggered_average(fid.TimeStamps,channels(1),.01,[4,1],treward(idx));
    if length(channels)==2,
        subplot(3,3,6)
        [rta_early{2},rta_late{2}] = reward_triggered_average(fid.TimeStamps,channels(2),.01,[4,1],treward(idx));
        subplot(3,3,3)
        [trial_rates,r_task] = spike_correlation_by_trial(fid.TimeStamps,channels,[1,0],treward(idx));
    else
        trial_rates = [];
        r_task = [];
    end
    
    % save relevant variables in mat file
    saveas(gcf,fullfile(savedir,name),'fig')
    print(gcf,fullfile(savedir,name),'-dpdf','-bestfit')
    save(fullfile(savedir,name),...
        'tstart','treward','beta','learning_fun','channels','epochSize',...
        'sleep_idx','r_sleep','r_awake','r_base',...
        'FR_sleep','FF_sleep','FR_awake','FF_awake','FR_base','FF_base',...
        'r_task','rta_early','rta_late','trial_rates')
end % experiments



