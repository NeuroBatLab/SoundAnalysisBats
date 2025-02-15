% BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
% BaseCodeDir = 'C:\Users\BatLab\Documents\GitHub\';
% WorkingDir = 'C:\Users\BatLab\Documents\DeafWhoWorkDir\';
% Path2RecordingTable = 'C:\Users\BatLab\Google Drive\JuvenileRecordings\DeafRecordingsNWAF155_Log.xlsx';
% TTLFolder = 'C:\Users\BatLab\Google Drive\JuvenileRecordings';

% BaseDataDir = '/Volumes/Julie4T/JuvenileRecordings151/';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub';
WorkingDir = '/Users/elie/DeafWhoWorkDir/';
% %Path2RecordingTable = '/Users/elie/Google Drive/JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx';
% Path2RecordingTable = '/Users/elie/Google Drive/JuvenileRecordings/JuvenileRecordingsNWAF155_Log.xlsx';
% TTLFolder = '/Users/elie/Documents/zero_playback_12h';

addpath(genpath(fullfile(BaseCodeDir,'LMC')))
addpath(genpath(fullfile(BaseCodeDir, 'LoggerDataProcessing')))
addpath(genpath(fullfile(BaseCodeDir,'SoundAnalysisBats')))
addpath(genpath(fullfile(BaseCodeDir,'GeneralCode')))

DatesDir = dir(fullfile(BaseDataDir,'20*'));
NDates = length(DatesDir);
ExpLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSal.txt');
WhoLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWho.txt');
AlliLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalAllignement.txt');

if ~exist(ExpLog, 'file')
    error('Cannot find the list of file to run in: %s \n',ExpLog);
else
    FidExp = fopen(ExpLog, 'r');
    Header = textscan(FidExp,'%s\t%s\t%s\t%s\t%s\n',1);
    DoneListDetect = textscan(FidExp,'%s\t%s\t%s\t%.1f\t%d');
    fclose(FidExp);
end

if ~exist(WhoLog, 'file')
    FidWho = fopen(WhoLog, 'a');
    fprintf(FidWho, 'Subject\tDate\tTime\tNCalls\n');
    DoneListWho = [];
else
    FidWho = fopen(WhoLog, 'r');
    Header = textscan(FidWho,'%s\t%s\t%s\t%s\n',1);
    DoneListWho = textscan(FidWho,'%s\t%s\t%s\t%d');
    fclose(FidWho);
    FidWho = fopen(WhoLog, 'a');
end

if ~exist(AlliLog, 'file')
    FidAlli = fopen(AlliLog, 'a');
    fprintf(FidAlli, 'Subject\tDate\tTime\tAllignement\n');
    ListAlliOk = [];
else
    FidAlli = fopen(AlliLog, 'r');
    Header = textscan(FidAlli,'%s\t%s\t%s\t%s\n',1);
    ListAlliOk = textscan(FidAlli,'%s\t%s\t%s\t%d');
    fclose(FidAlli);
    FidAlli = fopen(AlliLog, 'a');
end

NExpe = length(DoneListDetect{1});

for ee=17:NExpe
    BatsID = DoneListDetect{1}{ee};
    Date = DoneListDetect{2}{ee};
    Time = DoneListDetect{3}{ee};
    ParamFile = dir(fullfile(BaseDataDir,['20' Date],'audio',sprintf('%s_%s_%s*RecOnly_param.txt', BatsID, Date, Time)));
    fprintf(1, '\n\n\n Date: %s, experiment %d/%d\n%s\n', Date,ee,NExpe,ParamFile.name)
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

    if ~isempty(ListAlliOk)
        AlliOkInd = find(contains(ListAlliOk{1},BatsID) .* contains(ListAlliOk{2},Date) .* contains(ListAlliOk{3},Time));
        AlliOk = ListAlliOk{4}(AlliOkInd);
    else
        AlliOk=[];
    end

    Filepath = fullfile(ParamFile.folder, ParamFile.name);
    NCalls = result_reconly_DbatsWho(Filepath, WorkingDir, AlliOk, FidAlli);
    fprintf(FidWho, '%s\t%s\t%s\t%d\n',ParamFile.name(1:4),ParamFile.name(6:11),ParamFile.name(13:16),NCalls);
end
fclose(FidWho);

%% INTERNAL FUNCTION
function [NCalls] = result_reconly_DbatsWho(Path2ParamFile,WorkingDir,AlliOk,FidAlli, Logger_dir)
ForceWhoID = 1; % In case the identification of bats was already done but you want to re-do it again

close all

% Get the recording date
[AudioDataPath, DataFile ,~]=fileparts(Path2ParamFile);
BatsID = DataFile(1:4);
Date = DataFile(6:11);
ExpStartTime = DataFile(13:16);

if nargin<5
    % Set the path to logger data
    Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'audiologgers');
        
end

if isempty(AlliOk)
    fprintf(1, '   -> Starting a new session\n')
    fprintf(1,'*** Check the clock drift correction of the logger ***\n')
    LoggersDir = dir(fullfile(Logger_dir, 'logger*'));
    Check = zeros(length(LoggersDir)+1,1);
    for ll=1:length(LoggersDir)
        FigCD = open(fullfile(LoggersDir(ll).folder, LoggersDir(ll).name,'extracted_data','CD_correction0.fig'));
        %                 fprintf(1, 'Go in %s\n',fullfile(BaseDir,sprintf('box%d',BoxOfInterest(bb)),'piezo',Date,'audiologgers','loggerxx','extracted_data'))
        %                 fprintf(1,'Open CD_correction0\n')
        Check(ll) = input('Is everything ok? (yes ->1, No -> 0): ');
        fprintf('\n')
        close(FigCD)
    end
    fprintf(1,'*** Check the allignement of the TTL pulses ***\n')
    AllignmentPath = fullfile(AudioDataPath,sprintf('%s_%s_CD_correction_audio_piezo.fig', Date, ExpStartTime));
    FigAP = open(AllignmentPath);
    %                 fprintf(1, 'Go in %s\n',fullfile(BaseDir,sprintf('box%d',BoxOfInterest(bb)),'bataudio'))
    %                 fprintf(1,'Search for %s_%s_CD_correction_audio_piezo\n', Date, Time)
    Check(length(LoggersDir)+1) = input('Is everything ok? (yes ->1, No -> 0): ');
    fprintf('\n')
    close(FigAP)
    if any(~Check)
        NCalls = nan;
        AlliOk=0;
        fprintf(FidAlli, '%s\t%s\t%s\t%d\n',BatsID,Date,ExpStartTime,AlliOk);
        fprintf(1,'\n****** Error in allignement reported ******\n')
        return
    else
        AlliOk=1;
        fprintf(FidAlli, '%s\t%s\t%s\t%d\n',BatsID,Date,ExpStartTime,AlliOk);
        fprintf(1,'\n****** Allignement reported as good! ******\n')
    end
elseif ~AlliOk
    fprintf(1, '   -> Session flagged as not alligned correctly\n')
    return
else
    fprintf(1, '   -> Starting from where we left on this session\n')
end
%% Identify who is calling
fprintf('\n*** Identify who is calling ***\n')
WhoCall_dir = dir(fullfile(Logger_dir, sprintf('*%s_%s*whocalls*', Date, ExpStartTime)));
if isempty(WhoCall_dir) || ForceWhoID
%     [IndVocStartRawMerged,~]=who_calls_playless(AudioDataPath,Logger_dir,Date, ExpStartTime,200,1,1,0,'UseResNetNoiseFilter',1, 'Sorter', 'JulieE','Working_dir',WorkingDir);
    [IndVocStartRawMerged,~]=who_calls_playless_ResNet(AudioDataPath,Logger_dir,Date, ExpStartTime,200,1,0, 'Sorter', 'JulieE','Working_dir',WorkingDir);
else
    fprintf('\n*** ALREADY DONE: Identify who is calling ***\n')
end
NCalls = sum(cellfun('length',IndVocStartRawMerged));

end
