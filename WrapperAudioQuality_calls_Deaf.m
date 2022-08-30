%% This script run through files that have been manually cured to apply audioQuality_calls_deaf.m
Redo = 1; % 0: do not rerun dates that have already been run; 1 rerun all 
Redo2=1; % 0: do not rerun sets that have already been run and continue from where we left; 1 rerun all 
%% These are specific to the dataset and computer
% BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
% BaseCodeDir = 'C:\Users\tobias\Documents\GitHub\';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
% WorkingDir = 'C:\Users\tobias\Documents\DeafWhoWorkDir\';
WhatLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWhat.txt');
QualityLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalQuality.txt');


%% Paths to code
% you should have pulled from github the last versions of
% LMC,LoggerDataProcessing and SoundAnalysisBats
addpath(genpath(fullfile(BaseCodeDir,'LMC')))
addpath(genpath(fullfile(BaseCodeDir, 'LoggerDataProcessing')))
addpath(genpath(fullfile(BaseCodeDir,'SoundAnalysisBats')))
addpath(genpath(fullfile(BaseCodeDir,'General')))

%% Get the name of the next experiment that needs to be manually curated
ee=0;
% these are the list of files to run
if ~exist(WhatLog, 'file')
    error('Cannot find the list of file to run in: %s \n', WhatLog);
else
    FidWhat = fopen(WhatLog, 'r');
    Header = textscan(FidWhat,'%s\t%s\t%s\n',1);
    ListOfFilesToDo = textscan(FidWhat,'%s\t%s\t%s');
    fclose(FidWhat);
end

% This is the list of files that has been processed
if ~exist(QualityLog, 'file')
    FidQuality = fopen(QualityLog, 'a');
    fprintf(FidQuality, 'Subject\tDate\tTime\n');
    DoneListQuality = [];
else
    FidQuality = fopen(QualityLog, 'r');
    Header = textscan(FidQuality,'%s\t%s\t%s\n',1);
    DoneListQuality = textscan(FidQuality,'%s\t%s\t%s');
    fclose(FidQuality);
    FidQuality = fopen(QualityLog, 'a');
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
    
    % 2- Now Check that this file/experiment listed in WhatLog is not already listed in
    % QualityLog (if they are listed in QualityLog it means they've already been
    % processed by audioQuality_calls_deaf and we might not want to run them and overwrite
    % the data)
    if ~(str2double(Date)>200122)  % these recordings where done not in the final set up and/or had issue with the extraction pipeline
        continue
    end
    if ~isempty(DoneListQuality)
        Done = sum(contains(DoneListQuality{1},BatID) .* contains(DoneListQuality{2},Date) .* contains(DoneListQuality{3},ExpStartTime));
    else
        Done = 0;
    end
    
    if Done && ~Redo % the manual evaluation (audioQuality_calls_deaf) has already been applied to this experiment, proceed to the next ff
        continue
    else
        % Give the path to the data for that ff experiment
        ParamFile = dir(fullfile(BaseDataDir, ['20' Date], 'audio', sprintf('%s_%s_%s*RecOnly_param.txt', BatID, Date, ExpStartTime)));
        Logger_dir = fullfile(ParamFile.folder(1:(strfind(ParamFile.folder, 'audio')-1)), 'audiologgers');
        % 3- Once you've identified an experimental date that has not been run
        % apply audioQuality_calls_deaf
        fprintf(1, '***********************************************\n* Running what_calls on:\n %s\n Date: %s\n ExpStartTime:%s *\n***********************************************\n', Logger_dir, Date, ExpStartTime)
        audioQuality_calls_Deafs(Logger_dir, Date, ExpStartTime,Redo2)
        %
        
      % Records results in Qualitylog  
      if ~Redo
        fprintf(FidQuality, '%s\t%s\t%s\n', BatID, Date, ExpStartTime);
        %pause()
      end
      
      
    end
    
end

fclose(FidWhat);
