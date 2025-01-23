%% This script run through files that have been manually cured with ...
% audioQuality_calls_deaf.m to select onset and offset of true...
% vocalizations in IndVocStartPiezo and IndVocStopPiezo

%% These are specific to the dataset and computer
% BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';

QualityLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalQuality.txt');


%% Paths to code
% you should have pulled from github the last versions of
% LMC,LoggerDataProcessing and SoundAnalysisBats
addpath(genpath(fullfile(BaseCodeDir,'LMC')))
addpath(genpath(fullfile(BaseCodeDir, 'LoggerDataProcessing')))
addpath(genpath(fullfile(BaseCodeDir,'SoundAnalysisBats')))
addpath(genpath(fullfile(BaseCodeDir,'General')))

%% Get the name of the next experiment that needs to be corrected
ee=0;
% these are the list of files to run
if  ~exist(QualityLog, 'file')
    error('Cannot find the list of file to run in: %s \n', QualityLog);
else
    FidQuality = fopen(QualityLog, 'r');
    Header = textscan(FidQuality,'%s\t%s\t%s\n',1);
    ListOfFilesToDo = textscan(FidQuality,'%s\t%s\t%s');
    fclose(FidQuality);
end

%% Run audioQuality_correct_logger_data_voc on the selected data
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

    if ~(str2double(Date)>200122)  % these recordings where done not in the final set up and/or had issue with the extraction pipeline
        fprintf(1, 'Wrong dataset, Date = %s   Start Time = %s\n', Date,ExpStartTime)
        continue
    end


    % Give the path to the data for that ff experiment
    ParamFile = dir(fullfile(BaseDataDir, ['20' Date], 'audio', sprintf('%s_%s_%s*RecOnly_param.txt', BatID, Date, ExpStartTime)));
    Logger_dir = fullfile(ParamFile.folder(1:(strfind(ParamFile.folder, 'audio')-1)), 'audiologgers');
    % 2- Once you've identified an experimental date
    % apply audioQuality_correct_logger_data_voc
    fprintf(1, '***********************************************\n* Running audioQuality_correct_logger_data_voc on:\n %s\n Date: %s\n ExpStartTime:%s *\n***********************************************\n', Logger_dir, Date, ExpStartTime)
    audioQuality_correct_logger_data_voc(Logger_dir, Date, ExpStartTime)

end



