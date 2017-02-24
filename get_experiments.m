function expts = get_experiments()
%% expts = get_experiments()
%   function to grab metadata for a few BMI learning experiments
%   returns a structure with the following fields
%       'rat' (string)
%       'baseline_session' (number)
%       'task_session' (vector of one or two number)
%       'channels' (vector of one or two numbers)

expts(1).rat = 'S23';
expts(1).baseline_session = '64';
expts(1).task_session = 'Cat65-70Task';
expts(1).sleep_session = '';
expts(1).channels = [15,21];

expts(2).rat = 'S23';
expts(2).baseline_session = '79';
expts(2).task_session = 'Cat80-84Task';
expts(2).sleep_session = '';
expts(2).channels = [26,29];

expts(3).rat = 'S23';
expts(3).baseline_session = '92';
expts(3).task_session = 'Cat94-98Task';
expts(3).sleep_session = 'Cat92-102Sleep';
expts(3).channels = [26,29];

expts(4).rat = 'S32';
expts(4).baseline_session = '341';
expts(4).task_session = 'Cat343-344Task';
expts(4).sleep_session = 'Cat341-345Sleep1';
expts(4).channels = [19,31];

expts(5).rat = 'S32';
expts(5).baseline_session = '364';
expts(5).task_session = 'Cat365-366Task';
expts(5).sleep_session = 'Cat364-370Sleep';
expts(5).channels = [22,31];

expts(6).rat = 'S34';
expts(6).baseline_session = '297';
expts(6).task_session = 'Cat298-299Task';
expts(6).sleep_session = '';
expts(6).channels = [27,6];

expts(7).rat = 'S34';
expts(7).baseline_session = '316';
expts(7).task_session = 'Cat317Task';
expts(7).sleep_session = 'Cat315-318Sleep';
expts(7).channels = [6,27];

expts(8).rat = 'T10';
expts(8).baseline_session = '462';
expts(8).task_session = 'Cat463-464Task';
expts(8).sleep_session = 'Cat462-465Sleep';
expts(8).channels = [5,7];

expts(9).rat = 'T10';
expts(9).baseline_session = '475';
expts(9).task_session = 'Cat476Task';
expts(9).sleep_session = '';
expts(9).channels = [5,7];

expts(10).rat = 'T12';
expts(10).baseline_session = '497';
expts(10).task_session = 'Cat498-499Task';
expts(10).sleep_session = 'Cat496-501Sleep';
expts(10).channels = [16,19];

expts(11).rat = 'T12';
expts(11).baseline_session = '528';
expts(11).task_session = 'Cat529-530Task';
expts(11).sleep_session = 'Cat528-532Sleep';
expts(11).channels = 17;

expts(12).rat = 'T12';
expts(12).baseline_session = '543';
expts(12).task_session = 'Cat544-545Task';
expts(12).sleep_session = 'Cat542-546Sleep';
expts(12).channels = [25,12];

% expts(13).rat = 'T17';
% expts(13).baseline_session = '606';
% expts(13).task_session = 'Cat607-608Task';
% expts(13).sleep_session = 'Cat606-611Sleep_Early';
% expts(13).channels = 31;
% 
% expts(14).rat = 'T17';
% expts(14).baseline_session = '611';
% expts(14).task_session = 'Cat612Task';
% expts(14).sleep_session = 'Cat606-611Sleep_Early';
% expts(14).channels = 6;



