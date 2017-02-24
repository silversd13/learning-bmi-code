%% Summary Script
clear, clc, close all

%% grab summary data, if does not exist for a session, fill in with nan
datadir = fullfile(pwd,'processed');
experiments = get_experiments();

X = zeros(length(experiments),15);
Y = zeros(length(experiments),3);
for e=1:length(experiments),
    % load
    datafile = sprintf('experiment_%02i',e);
    load(fullfile(datadir,datafile));
    
    % replace empties with nans
    if isempty(r_base),
        r_base = nan;
        FR_base = nan;
        FF_base = nan;
    end
    if isempty(r_sleep),
        r_sleep = nan;
        FR_sleep = nan(1,2);
        FF_sleep = nan(1,2);
    end
    if isempty(r_awake),
        r_awake = nan;
        FR_awake = nan(1,2);
        FF_awake = nan(1,2);
    end
    if length(FR_base)==1,
        FR_base = [FR_base,nan];
        FF_base = [FF_base,nan];
    end
    if length(FR_sleep)==1,
        FR_sleep = [FR_sleep,nan];
        FF_sleep = [FF_sleep,nan];
    end
    if length(FR_awake)==1,
        FR_awake = [FR_awake,nan];
        FF_awake = [FF_awake,nan];
    end
    
    % store
    X(e,:) = [...
        r_base,FR_base,FF_base,...
        r_sleep,FR_sleep,FF_sleep,...
        r_awake,FR_awake,FF_awake];
    Y(e,:) = beta;
end

%% start with linear regression
B = regress(Y(:,1),addones(X));









