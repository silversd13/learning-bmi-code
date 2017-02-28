%% Go through analysis for each session
clear, clc, close all
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
        
    %% baseline activity
    % load baseline session to look at spontaneous activity
    session = experiments(experiment).baseline_session;
    datafile = fullfile(datadir,sprintf('data_block_%s_%s',rat,session));
    fid = load(datafile,'wave','Fs_wave');
    N = floor(length(fid.wave)/fid.Fs_wave);
    win = [0,4];
    t = (0:win(2):N-win(2))';
    
    % calc firing rates, fano factors, and correlation btw neurons
    fid = load(datafile,'TimeStamps');
    subplot(3,3,2)
    [FR_base,FF_base,ISI_base] = neuron_stats(fid.TimeStamps,channels,.1,t,win);
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
    
%     % save relevant variables in mat file
%     saveas(gcf,fullfile(savedir,name),'fig')
%     print(gcf,fullfile(savedir,name),'-dpdf','-bestfit')
%     save(fullfile(savedir,name),...
%         'tstart','treward','beta','learning_fun','channels','epochSize',...
%         'sleep_idx','r_sleep','r_awake','r_base',...
%         'FR_sleep','FF_sleep','FR_awake','FF_awake','FR_base','FF_base',...
%         'r_task','rta_early','rta_late','trial_rates')

end % experiments



