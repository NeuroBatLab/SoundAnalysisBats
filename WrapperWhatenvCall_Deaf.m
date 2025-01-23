%% This script run through files that have been manually cured and had run...
% through the old version of what_calls that was not gathering the amplitude...
% envelopes to apply what_calls_patch4env.m to gather the amplitude
% envelopes
Redo = 0; % 0: do not rerun things that have already been run; 1 rerun all 
%% These are specific to the dataset and computer
% BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
% BaseCodeDir = 'C:\Users\tobias\Documents\GitHub\';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
% WorkingDir = 'C:\Users\tobias\Documents\DeafWhoWorkDir\';
WhatLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWhat.txt');
% WhatenvLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWhatenv.txt');

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
if ~exist(WhatLog, 'file')
    error('Cannot find the list of file to run in: %s \n', WhatLog);
else
    FidWhat = fopen(WhatLog, 'r');
    Header = textscan(FidWhat,'%s\t%s\t%s\t%s\t%s\n',1);
    ListOfFilesToDo = textscan(FidWhat,'%s\t%s\t%s\t%.1f\t%d');
    fclose(FidWhat);
end

% % This is the list of files that has been processed
% if ~exist(WhatenvLog, 'file')
%     FidWhatenv = fopen(WhatenvLog, 'a');
%     fprintf(FidWhatenv, 'Subject\tDate\tTime\n');
%     DoneListWhatenv = [];
% else
%     FidWhatenv = fopen(WhatenvLog, 'r');
%     Header = textscan(FidWhatenv,'%s\t%s\t%s\n',1);
%     DoneListWhatenv = textscan(FidWhatenv,'%s\t%s\t%s');
%     fclose(FidWhatenv);
%     FidWhatenv = fopen(WhatenvLog, 'a');
% end

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
    if (str2double(Date)==200123) || (str2double(Date)==200124) || (str2double(Date)==200127) ||(str2double(Date)==200129)  % these recordings are done
        continue
    end
%     if ~isempty(DoneListWhatenv)
%         Done = sum(contains(DoneListWhatenv{1},BatID) .* contains(DoneListWhatenv{2},Date) .* contains(DoneListWhatenv{3},ExpStartTime));
%     else
%         Done = 0;
%     end
%     
%     if Done && ~Redo % Biosound (what_calls) has already been applied to this experiment, proceed to the next ff
%         continue
%     else
        % Give the path to the data for that ff experiment
        ParamFile = dir(fullfile(BaseDataDir, ['20' Date], 'audio', sprintf('%s_%s_%s*RecOnly_param.txt', BatID, Date, ExpStartTime)));
        Logger_dir = fullfile(ParamFile.folder(1:(strfind(ParamFile.folder, 'audio')-1)), 'audiologgers');
        % 3- Once you've identified an experimental date that has not been run
        % through whatcalls, apply what_calls
        fprintf(1, '***********************************************\n* Running what_calls_patch4env on:\n %s\n Date: %s\n ExpStartTime:%s *\n***********************************************\n', Logger_dir, Date, ExpStartTime)
        what_calls_patch4env(Logger_dir, Date, ExpStartTime)
        %
        
%       % Add code that records results in whatlog  
%       if ~Redo
%         fprintf(FidWhatenv, '%s\t%s\t%s\n', BatID, Date, ExpStartTime);
%         %pause()
%       end
      
      
%     end
    
end

% fclose(FidWhatenv);
