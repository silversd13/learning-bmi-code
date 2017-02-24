%% BMI Analysis

epochSize = 4; % seconds
datadir = '/Volumes/data/RatBMI/data';

experiments = get_experiments();
for experiment=1%:length(experiments),
    
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
