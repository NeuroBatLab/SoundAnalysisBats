%% This script run through files that have been manually cured to apply what_calls.m
%% These are specific to the dataset and computer
BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
BaseCodeDir = 'C:\Users\tobias\Documents\GitHub\';
% WorkingDir = 'C:\Users\tobias\Documents\DeafWhoWorkDir\';
WhoLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWho.txt');
WhatLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWhat.txt');

%% Paths to code
% you should have pulled from github the last versions of
% LMC,LoggerDataProcessing and SoundAnalysisBats
addpath(genpath(fullfile(BaseCodeDir,'LMC')))
addpath(genpath(fullfile(BaseCodeDir, 'LoggerDataProcessing')))
addpath(genpath(fullfile(BaseCodeDir,'SoundAnalysisBats')))

%% Parameters of WhatCalls



%% Get the name of the next experiment that needs to be manually curated
ee=0;
% these are the list of files to run (replace with WhoLog)
if ~exist(WhatLog, 'file')
    error('Cannot find the list of file to run in: %s \n', WhatLog);
else
    FidWhat = fopen(WhatLog, 'r');
    Header = textscan(FidWhat,'%s\t%s\t%s\t%s\t%s\n',1);
    ListOfFilesToDo = textscan(FidWhat,'%s\t%s\t%s\t%.1f\t%d');
    fclose(FidWhat);
end

% This is the list of files that has been processed (to replace with the new
% WhatLog)
if ~exist(WhatLog, 'file') %% Here and below (down to line 44) you want to replace by WhatLog, a file that you are going to create to keep track of which file went through the Biosound analysis
    FidWhat = fopen(WhatLog, 'a'); % Rename h                        ere and below (down to line 44) FidWho for FidWhat
    fprintf(FidWhat, 'Subject\tDate\tTime\n');
    DoneListWhat = [];
else
    FidWhat = fopen(WhatLog, 'r');
    Header = textscan(FidWhat,'%s\t%s\t%s\n',1);
    DoneListWhat = textscan(FidWhat,'%s\t%s\t%s');
    fclose(FidWhat);
    FidWhat = fopen(WhatLog, 'a');
end

%% Run what_calls on the selected data
% 1- Do a loop through all the files that need to run
% 1.1- Get the number of experiments from the content of WhoLog that was
% read into ListOfFilesToDo
NToDo = length(ListOfFilesToDo{1});
% 1.2- Loop through these files
for ff=1:NToDo
    % Get the BatID the Date and the time of the ff experiment (use
    % ListOfFilesToDo)
    BatID = ListOfFilesToDo{1}{ff}
    Date = ListOfFilesToDo{2}{ff}
    ExpStartTime = ListOfFilesToDo{3}{ff}
    
    % 2- Now Check that this file/experiment listed in WhoLog is not already listed in
    % WhatLog (if they are listed in WhatLog it means they've already been
    % processed by WhatCalls and we might not want to run them and overwrite
    % the data)
    
    if ~isempty(WhatLog) %% WhatLog is a file path and not the content of the file that is read line 41 you want to replace WhatLog by DoneListWhat
        Done = sum(contains(DoneListWhat{1},BatID) .* contains(DoneListWhat{2},Date) .* contains(DoneListWhat{3},ExpStartTime));
    else
        Done = 0;
    end
    
    if Done % Biosound (what_calls) has already been applied to this experiment, proceed to the next ff
        continue
    else
        % Give the path to the data for that ff experiment
        ParamFile = dir(fullfile(BaseDataDir, ['20' Date], 'audio', sprintf('%s_%s_%s*RecOnly_param.txt', BatID, Date, ExpStartTime)));
        Logger_dir = fullfile(ParamFile.folder(1:(strfind(ParamFile.folder, 'audio')-1)), 'audiologgers');
        % Apply what_calls
        % 3- Once you've identified an experimental date that has not been run
        % through whatcalls, apply whatcalls
        % what_calls(Logger_dir, Date, ExpStartTime)
        fprintf(1, 'This is where I would call what_calls with Logger_dir being:\n %s\n Date: %s\n ExpStartTime:%s\n', Logger_diir, Date, ExpStartTime)
        pause()
        
    end
    
    
    
end

