%% Sleep Analysis

epochSize = 4; % seconds
datadir = '/Volumes/data/RatBMI/data';

experiments = get_experiments();
for experiment=4%:length(experiments),
    
    % make sure directories for saving exists
    savedir = fullfile('~','Projects','LearningBMI','processed');
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
        pks = find_slow_oscillations(fid.data',fid.Fs_lfp,epochSize,sleep_idx);
                
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
    
%     % save relevant variables in mat file
%     saveas(gcf,fullfile(savedir,name),'fig')
%     print(gcf,fullfile(savedir,name),'-dpdf','-bestfit')
%     save(fullfile(savedir,name),...
%         'tstart','treward','beta','learning_fun','channels','epochSize',...
%         'sleep_idx','r_sleep','r_awake','r_base',...
%         'FR_sleep','FF_sleep','FR_awake','FF_awake','FR_base','FF_base',...
%         'r_task','rta_early','rta_late','trial_rates')
    
end % experiments

