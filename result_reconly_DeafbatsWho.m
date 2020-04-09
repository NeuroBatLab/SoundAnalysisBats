BaseDataDir = 'X:\JulieE\DeafSalineGroup151\';
BaseCodeDir = 'C:\Users\Batman\Documents\Code\';
Path2RecordingTable = 'C:\Users\Batman\GoogleDriveNeuroBatGroup\JuvenileRecordings\DeafRecordingsNWAF155_Log.xlsx';
TTLFolder = 'C:\Users\Batman\GoogleDriveNeuroBatGroup\JuvenileRecordings';

% BaseDataDir = '/Volumes/Julie4T/JuvenileRecordings151/';
% BaseCodeDir = '/Users/elie/Documents/CODE';
% %Path2RecordingTable = '/Users/elie/Google Drive/JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx';
% Path2RecordingTable = '/Users/elie/Google Drive/JuvenileRecordings/JuvenileRecordingsNWAF155_Log.xlsx';
% TTLFolder = '/Users/elie/Documents/zero_playback_12h';

addpath(genpath(fullfile(BaseCodeDir,'LMC')))
addpath(genpath(fullfile(BaseCodeDir, 'LoggerDataProcessing')))
addpath(genpath(fullfile(BaseCodeDir,'SoundAnalysisBats')))

DatesDir = dir(fullfile(BaseDataDir,'20*'));
NDates = length(DatesDir);
ExpLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSal.txt');
WhoLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWho.txt');

if ~exist(ExpLog, 'file')
    error('Cannot find the list of file to run in: %s \n',ExpLog);
else
    FidExp = fopen(ExpLog, 'r');
    Header = textscan(FidExp,'%s\t%s\t%s\t%s\t%s\n');
    DoneListDetect = textscan(FidExp,'%s\t%s\t%s\t%.1f\t%d');
    fclose(FidExp);
end

if ~exist(WhoLog, 'file')
    FidWho = fopen(WhoLog, 'a');
    fprintf(FidWho, 'Subject\tDate\tTime\tNCalls\n');
    DoneListWho = [];
else
    FidWho = fopen(WhoLog, 'r');
    Header = textscan(FidWho,'%s\t%s\t%s\t%s\n');
    DoneListWho = textscan(FidWho,'%s\t%s\t%s\t%d');
    fclose(FidWho);
    FidWho = fopen(WhoLog, 'a');
end

NExpe = length(DoneListDetect{1});

for ee=1:NExpe
    BatsID = DoneListDetect{1}{ee};
    Date = DoneListDetect{2}{ee};
    Time = DoneListDetect{3}{ee};
    ParamFile = dir(fullfile(BaseDataDir,['20' Date],'audio',sprintf('%s_%s_%s*RecOnly_param.txt', BatsID, Date, Time)));
    fprintf(1, '\n\n\n Date: %s, experiment %d/%d\n%s\n', Date,dd,NExpe,ParamFile.name)
    % Check that the file was not already set aside or done
    
    if ~isempty(DoneListWho)
        Done = sum(contains(DoneListWho{1},BatsID) .* contains(DoneListWho{2},Date) .* contains(DoneListWho{3},Time));
    else
        Done=0;
    end
    if Done
        fprintf(1, '   -> Data already processed\n')
        continue
    end
    Filepath = fullfile(ParamFile.folder, ParamFile.name);
    NCalls = result_reconly_DbatsWho(Filepath,Path2RecordingTable,TTLFolder);
    fprintf(FidWho, '%s\t%s\t%s\t%d\n',ParamFile.name(1:4),ParamFile.name(6:11),ParamFile.name(13:16),NCalls);
end
close(FidWho)

%% INTERNAL FUNCTION
function [NCalls] = result_reconly_DbatsWho(Path2ParamFile, Logger_dir)

ForceWhoID = 0; % In case the identification of bats was already done but you want to re-do it again
close all

% Get the recording date
[AudioDataPath, DataFile ,~]=fileparts(Path2ParamFile);
Date = DataFile(6:11);

if nargin<2
    % Set the path to logger data
    Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'audiologgers');
        
end

%% Identify who is calling
fprintf('\n*** Identify who is calling ***\n')
WhoCall_dir = dir(fullfile(Logger_dir, sprintf('*%s_%s*whocalls*', Date, ExpStartTime)));
if isempty(WhoCall_dir) || ForceWhoID
    [IndVocStartRawMerged,~]=who_calls(AudioDataPath,Logger_dir,Date, ExpStartTime,200,1,1,0);
else
    fprintf('\n*** ALREADY DONE: Identify who is calling ***\n')
end
NCalls = sum(cellfun('length',IndVocStartRawMerged));

end
