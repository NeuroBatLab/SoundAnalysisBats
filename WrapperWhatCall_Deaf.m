%% This script run through files that have been manually cured to apply what_calls.m
Redo = 0; % 0: do not rerun things that have already been run; 1 rerun all 
%% These are specific to the dataset and computer
% BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
% BaseCodeDir = 'C:\Users\tobias\Documents\GitHub\';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
% WorkingDir = 'C:\Users\tobias\Documents\DeafWhoWorkDir\';
WhoLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWho.txt');
WhatLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWhat.txt');

%% Paths to code
% you should have pulled from github the last versions of
% LMC,LoggerDataProcessing and SoundAnalysisBats
addpath(genpath(fullfile(BaseCodeDir,'LMC')))
addpath(genpath(fullfile(BaseCodeDir, 'LoggerDataProcessing')))
addpath(genpath(fullfile(BaseCodeDir,'SoundAnalysisBats')))
addpath(genpath(fullfile(BaseCodeDir,'General')))

%% Parameters of WhatCalls



%% Get the name of the next experiment that needs to be manually curated
ee=0;
% these are the list of files to run
if ~exist(WhoLog, 'file')
    error('Cannot find the list of file to run in: %s \n', WhoLog);
else
    FidWho = fopen(WhoLog, 'r');
    Header = textscan(FidWho,'%s\t%s\t%s\t%s\t%s\n',1);
    ListOfFilesToDo = textscan(FidWho,'%s\t%s\t%s\t%.1f\t%d');
    fclose(FidWho);
end

% This is the list of files that has been processed
if ~exist(WhatLog, 'file')
    FidWhat = fopen(WhatLog, 'a');
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
    BatID = ListOfFilesToDo{1}{ff};
    Date = ListOfFilesToDo{2}{ff};
    ExpStartTime = ListOfFilesToDo{3}{ff};
    
    % 2- Now Check that this file/experiment listed in WhoLog is not already listed in
    % WhatLog (if they are listed in WhatLog it means they've already been
    % processed by WhatCalls and we might not want to run them and overwrite
    % the data)
    if ~(str2double(Date)>200122)  % these recordings where done not in the final set up and/or had issue with the extraction pipeline
        continue
    end
    if ~isempty(DoneListWhat)
        Done = sum(contains(DoneListWhat{1},BatID) .* contains(DoneListWhat{2},Date) .* contains(DoneListWhat{3},ExpStartTime));
    else
        Done = 0;
    end
    
    if Done && ~Redo % Biosound (what_calls) has already been applied to this experiment, proceed to the next ff
        continue
    else
        % Give the path to the data for that ff experiment
        ParamFile = dir(fullfile(BaseDataDir, ['20' Date], 'audio', sprintf('%s_%s_%s*RecOnly_param.txt', BatID, Date, ExpStartTime)));
        Logger_dir = fullfile(ParamFile.folder(1:(strfind(ParamFile.folder, 'audio')-1)), 'audiologgers');
        % 3- Once you've identified an experimental date that has not been run
        % through whatcalls, apply what_calls
        fprintf(1, '***********************************************\n* Running what_calls on:\n %s\n Date: %s\n ExpStartTime:%s *\n***********************************************\n', Logger_dir, Date, ExpStartTime)
        what_calls(Logger_dir, Date, ExpStartTime,0,0,1,1)
        %
        
      % Add code that records results in whatlog  
      if ~Redo
        fprintf(FidWhat, '%s\t%s\t%s\n', BatID, Date, ExpStartTime);
        %pause()
      end
      
      
    end
    
end

fclose(FidWhat);
