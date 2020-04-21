BaseDataDir = 'Z:\JulieE\DeafSalineGroup151\';
BaseCodeDir = 'C:\Users\Batman\Documents\Code\';
Path2RecordingTable = 'C:\Users\Batman\Documents\GoogleDriveNeuroBatGroup\JuvenileRecordings\DeafRecordingsNWAF155_Log.xlsx';
TTLFolder = 'C:\Users\Batman\Documents\GoogleDriveNeuroBatGroup\JuvenileRecordings';

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

if ~exist(ExpLog, 'file')
    Fid = fopen(ExpLog, 'a');
    fprintf(Fid, 'Subject\tDate\tTime\tDuration(s)\tLoggerData\n');
    DoneList = [];
else
    Fid = fopen(ExpLog, 'r');
    Header = textscan(Fid,'%s\t%s\t%s\t%s\t%s\n',1);
    DoneList = textscan(Fid,'%s\t%s\t%s\t%.1f\t%d');
    fclose(Fid);
    Fid = fopen(ExpLog, 'a');
end

for dd=1:NDates
    ParamFile = dir(fullfile(DatesDir(dd).folder, DatesDir(dd).name,'audio','*RecOnly_param.txt'));
    fprintf(1, '\n\n\n Date: %s, experiment %d/%d\n%s\n', DatesDir(dd).name,dd,NDates,ParamFile.name)
    % Check that the file was not already treated
    BatsID = ParamFile.name(1:4);
    Date = ParamFile.name(6:11);
    Time = ParamFile.name(13:16);
    if ~isempty(DoneList)
        Done = sum(contains(DoneList{1},BatsID) .* contains(DoneList{2},Date) .* contains(DoneList{3},Time));
    else
        Done=0;
    end
    if Done
        fprintf(1, '   -> Data already processed\n')
        continue
    end
    Filepath = fullfile(ParamFile.folder, ParamFile.name);
    % check that the experiment has data!~
    fid = fopen(Filepath);
    data = textscan(fid,'%s','Delimiter', '\t');
    fclose(fid);
    
    % FIND THE LINE of your data
    IndexLine = find(contains(data{1}, 'Task stops at'));
    if ~isempty(IndexLine)
        IndexChar = strfind(data{1}{IndexLine},'after');
        IndexChar2 = strfind(data{1}{IndexLine},'seconds');
        
        % find the data into that line
        Temp = str2double(data{1}{IndexLine}((IndexChar + 6):(IndexChar2-2)));
        if Temp<600
            fprintf(1, '   -> Data too short\n')
            continue
        end
    end
    ProcessedOK = result_reconly_Dbats(Filepath,Path2RecordingTable,TTLFolder);
    Ind_ = strfind(ParamFile.name, '_param');
    fprintf(Fid, '%s\t%s\t%s\t%.1f\t%d\n',ParamFile.name(1:4),ParamFile.name(6:11),ParamFile.name(13:16),Temp,ProcessedOK);
end
fclose(Fid);

%% INTERNAL FUNCTION
function [Processed] = result_reconly_Dbats(Path2ParamFile, Path2RecordingTable, TTLFolder,Logger_dir)

TranscExtract = 1; % set to 1 to extract logger data and transceiver time
ForceExtract = 0; % set to 1 to redo the extraction of loggers otherwise the calculations will use the previous extraction data
ForceAllign = 0; % In case the TTL pulses allignment was already done but you want to do it again, set to 1
ForceVocExt1 = 0; % In case the localization on raw files of vocalizations that were manually extracted was already done but you want to do it again set to 1
ForceVocExt2 = 0; % In case the localization on Loggers of vocalizations that were manually extracted was already done but you want to do it again set to 1
ReAllignment = 0; % Incase we don't have a logger on all animals, it's better not to reallign the vocal data by cross correlation between the Microphone and the loggers
close all

% Get the recording date
[AudioDataPath, DataFile ,~]=fileparts(Path2ParamFile);
Date = DataFile(6:11);


if TranscExtract && nargin<2
    % Set the path to the recording log
    if contains(Path2ParamFile, 'LMC')
        Path2RecordingTable = '/Users/elie/Google Drive/BatmanData/RecordingLogs/recording_logs.xlsx';
    elseif contains(Path2ParamFile, 'Juvenile')
        Path2RecordingTable = '/Users/elie/Google Drive/JuvenileRecordings/JuvenileRecordingsNWAF155_Log.xlsx';
    elseif contains(Path2ParamFile, 'DeafSaline')
        Path2RecordingTable = '/Users/elie/Google Drive/JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx';
    end
end
if nargin<3
    % Set the path to the param file containing the TTL parameters
    TTLFolder = Path2RecordingTable;
end

if TranscExtract && nargin<4
    % Set the path to logger data
    if contains(Path2RecordingTable, 'BatmanData')
        Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'logger',['20' Date]);
    elseif contains(Path2RecordingTable, 'JuvenileRecording') || contains(Path2RecordingTable, 'DeafRecordings')
        Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'audiologgers');
    end
        
end

% Set the path to a working directory on the computer so logger data are
% transfered there and directly accessible for calculations
if TranscExtract
    WorkDir = ['~' filesep 'WorkingDirectory'];
end

%% Extracting sound events
% The samplestamp given by sound mex is not really reliable, so for each
% sound snippet, you want to find its exact location in the continuous
% recording files, then using TTL pulses, retrieve the time it correspond
% to in Deuteron, if requested.

% Checking what we have in terms of vocalization localization/extraction
ExpStartTime = DataFile(13:16);
VocExt_dir = dir(fullfile(AudioDataPath,sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)));

% Then run the logger extraction, allignment, and vocalization extraction
fprintf(1,'*** Extract Logger data if not already done ***\n');
% Find the ID of the recorded bats
if contains(Path2ParamFile, 'LMC')
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:P200','basic');
    RowData = find((cell2mat(RecTableData(2:end,1))== str2double(Date))) +1;
elseif contains(Path2ParamFile, 'Juvenile')
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:Q123','basic');
    RowData = find((cell2mat(RecTableData(4:end,1))== str2double(['20' Date]))) +3;
elseif contains(Path2ParamFile, 'DeafSaline')
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:AE41','basic');
    RowData = find((cell2mat(RecTableData(2:end,1))== str2double(['20' Date]))) +1;
end
DataInfo = RecTableData(RowData,:);
Header = RecTableData(1,:);
BatIDCol = find(contains(Header, 'Bat'));

% extract logger data if not already done
All_loggers_dir = dir(fullfile(Logger_dir, '*ogger*'));
DirFlags = [All_loggers_dir.isdir];
% Extract only those that are directories.
All_loggers_dir = All_loggers_dir(DirFlags);
TransceiverReset = struct(); % These are possible parameters for dealing with change of transceiver or sudden transceiver clock change. Set to empty before the first extraction
LoggerName = cell(length(All_loggers_dir),1);
BatID = cell(length(All_loggers_dir),1);
for ll=1:length(All_loggers_dir)
    Logger_i = fullfile(Logger_dir,All_loggers_dir(ll).name);
    Ind = strfind(All_loggers_dir(ll).name, 'r');
    Logger_num = str2double(All_loggers_dir(ll).name((Ind+1):end));
    NLogCol = find(contains(Header, 'NL'));% Columns of the neural loggers
    ALogCol = find(contains(Header, 'AL'));% Columns of the audio loggers
    LogCol = NLogCol(find(cell2mat(DataInfo(NLogCol))==Logger_num)); %#ok<FNDSB>
    if isempty(LogCol) % This is an audiologger and not a neural logger
        LogCol = ALogCol(find(cell2mat(DataInfo(ALogCol))==Logger_num)); %#ok<FNDSB>
        LoggerName{ll} = ['AL' num2str(Logger_num)];
    else
        LoggerName{ll} = ['NL' num2str(Logger_num)];
    end
    BatID{ll} = DataInfo{BatIDCol(find(BatIDCol<LogCol,1,'last'))};
    ParamFiles = dir(fullfile(Logger_i,'extracted_data','*extract_logger_data_parameters*mat'));
    if isempty(ParamFiles) || ForceExtract
        fprintf(1,'-> Extracting %s\n',All_loggers_dir(ll).name);
        
        % Bring data back on the computer
        Logger_local = fullfile(WorkDir, All_loggers_dir(ll).name);
        fprintf(1,'Transferring data from the server %s\n on the local computer %s\n', Logger_i, Logger_local);
        mkdir(Logger_local)
        [s,m,e]=copyfile(Logger_i, Logger_local, 'f');
        if ~s
            m %#ok<NOPRT>
            e %#ok<NOPRT>
            error('File transfer did not occur correctly for %s\n', Logger_i);
        end
        
        % run extraction
        if Logger_num==16 && str2double(Date)<190501
            % extract_logger_data(Logger_local, 'BatID', num2str(BatID), 'ActiveChannels', [0 1 2 3 4 5 6 7 8 9 10 12 13 14 15], 'AutoSpikeThreshFactor',5,'TransceiverReset',TransceiverReset)
            extract_logger_data(Logger_local, 'BatID', num2str(BatID{ll}), 'ActiveChannels', [0 1 2 3 4 5 6 7 8 9 10 12 13 14 15],'TransceiverReset',TransceiverReset,'AutoSpikeThreshFactor',4)
        else
            %extract_logger_data(Logger_local, 'BatID', num2str(BatID),'TransceiverReset',TransceiverReset)
            extract_logger_data(Logger_local, 'BatID', num2str(BatID{ll}),'TransceiverReset',TransceiverReset,'AutoSpikeThreshFactor',4)
        end
        
        % Keeps value of eventual clock reset
        Filename=fullfile(Logger_local, 'extracted_data', sprintf('%s_20%s_EVENTS.mat', num2str(BatID{ll}),Date));
        NewTR = load(Filename, 'TransceiverReset');
        if ~isempty(fieldnames(NewTR.TransceiverReset))% this will be used in the next loop!
            TransceiverReset = NewTR.TransceiverReset;
        end
        
        % Bring back data on the server
        fprintf(1,'Transferring data from the local computer %s\n back on the server %s\n', Logger_i, Logger_local);
        Remote_dir = fullfile(Logger_i, 'extracted_data');
        mkdir(Remote_dir)
        [s,m,e]=copyfile(fullfile(Logger_local, 'extracted_data'), Remote_dir, 'f');
        if ~s
            TicTransfer = tic;
            while toc(TicTransfer)<30*60
                [s,m,e]=copyfile(fullfile(Logger_local, 'extracted_data'), Remote_dir, 'f');
                if s
                    return
                end
            end
            if ~s
                s %#ok<NOPRT>
                m %#ok<NOPRT>
                e %#ok<NOPRT>
                error('File transfer did not occur correctly for %s\n Although we tried for 30min\n', Remote_dir);
            else
                fprintf('Extracted data transfered back on server in:\n%s\n',  Remote_dir);
            end
        else
            fprintf('Extracted data transfered back on server in:\n%s\n',  Remote_dir);
        end
        if s  %erase local data
            [sdel,mdel,edel]=rmdir(WorkDir, 's');
            if ~sdel
                TicErase = tic;
                while toc(TicErase)<30*60
                    [sdel,mdel,edel]=rmdir(WorkDir, 's');
                    if sdel
                        return
                    end
                end
            end
            if ~sdel
                sdel %#ok<NOPRT>
                mdel %#ok<NOPRT>
                edel %#ok<NOPRT>
                error('File erase did not occur correctly for %s\n Although we tried for 30min\n', WorkDir);
            end
        end
        
    else
        fprintf(1,'-> Already done for %s\n',All_loggers_dir(ll).name);
    end
end

if contains(Path2RecordingTable, 'BatmanData')
    % Get the serial numbers of the audiologgers that the two implanted
    % bats wear
    NLCol = find(contains(Header, 'NL'));
    ALThroatCol = find(contains(Header, 'AL-throat'));
    SerialNumberAL = nan(length(NLCol),1);
    SerialNumberNL = nan(length(NLCol),1);
    for dd=1:length(NLCol)
        SerialNumberAL(dd) = DataInfo{ALThroatCol(find(ALThroatCol<NLCol(dd),1,'last'))};
        SerialNumberNL(dd) = DataInfo{NLCol(dd)};
    end
end

% Alligning TTL pulses between soundmexpro and Deuteron
% for the RecOnly session

TTL_dir = dir(fullfile(AudioDataPath,sprintf( '%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
if isempty(TTL_dir) || ForceAllign
    fprintf(1,'\n*** Alligning TTL pulses for the RecOnly session ***\n');
    if contains(Path2RecordingTable, 'BatmanData')
        %         align_soundmexAudio_2_logger(AudioDataPath, Logger_dir, ExpStartTime,'TTL_pulse_generator','Avisoft','Method','risefall', 'Session_strings', {'all voc reward stop', 'rec only stop'}, 'Logger_list', [SerialNumberAL; SerialNumberNL]);
        align_soundmexAudio_2_logger(AudioDataPath, Logger_dir, ExpStartTime,'TTL_pulse_generator','Avisoft','Method','risefall', 'Session_strings', {'rec only start', 'rec only stop'}, 'Logger_list', [SerialNumberAL; SerialNumberNL]);
    elseif contains(Path2RecordingTable, 'JuvenileRecordings')
        align_soundmexAudio_2_logger(AudioDataPath, Logger_dir, ExpStartTime,'TTL_pulse_generator','Avisoft','Method','risefall', 'Session_strings', {'rec only start', 'rec only stop'}, 'TTLFolder',TTLFolder);
    end
    close all
else
    fprintf(1,'\n*** ALREADY DONE: Alligning TTL pulses for the free session ***\n');
end


%% Identify vocalizations using the piezo
if (isempty(VocExt_dir) || ForceVocExt1)
    fprintf(1,'\n*** Localizing and extracting vocalizations using the piezos ***\n');
    voc_localize_using_piezo(Logger_dir, AudioDataPath,Date, ExpStartTime)
elseif ~(isempty(VocExt_dir) || ForceVocExt1)
    fprintf(1,'\n*** ALREADY DONE: Localizing and extracting vocalizations using the piezos ***\n');
end

%% Identify the same vocalizations on the piezos and save sound extracts, onset and offset times
LogVoc_dir = dir(fullfile(Logger_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)));
if isempty(LogVoc_dir) || ForceVocExt1 || ForceVocExt2
    fprintf('\n*** Localizing vocalizations on piezo recordings ***\n')
    get_logger_data_voc(AudioDataPath, Logger_dir,Date, ExpStartTime, 'ReAllignment',ReAllignment);
else
    fprintf('\n*** ALREADY DONE: Localizing vocalizations on piezo recordings ***\n')
    
end

% Save the ID of the bat for each logger
Filename_ID = fullfile(Logger_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, 200));
if isfile(Filename_ID)
    save(Filename_ID, 'BatID','LoggerName','-append')
else
    save(Filename_ID, 'BatID','LoggerName')
end

Processed=1;

end
