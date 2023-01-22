
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
WhatLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWhat.txt');
%% Paths to code
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
    Header = textscan(FidWhat,'%s\t%s\t%s\t%s\t%s\n',1);
    ListOfFilesToDo = textscan(FidWhat,'%s\t%s\t%s\t%.1f\t%d');
    fclose(FidWhat);
end


%% Run the correction of BatName on the selected data
NToDo = length(ListOfFilesToDo{1});
%Loop through these files
for ff=1:NToDo
    % Get the BatID the Date and the time of the ff experiment (use
    % ListOfFilesToDo)
    BatID = ListOfFilesToDo{1}{ff};
    Date = ListOfFilesToDo{2}{ff};
    ExpStartTime = ListOfFilesToDo{3}{ff};
    
    if ~(str2double(Date)>200122) || str2double(Date)==200123  % these recordings where done not in the final set up and/or had issue with the extraction pipeline or have already run
        continue
    end
        % Give the path to the data for that ff experiment
        ParamFile = dir(fullfile(BaseDataDir, ['20' Date], 'audio', sprintf('%s_%s_%s*RecOnly_param.txt', BatID, Date, ExpStartTime)));
        Logger_dir = fullfile(ParamFile.folder(1:(strfind(ParamFile.folder, 'audio')-1)), 'audiologgers');
        % apply BatIDNameCorrect
        fprintf(1, '***********************************************\n* Running BatIDNameCorrect on:\n %s\n Date: %s\n ExpStartTime:%s *\n***********************************************\n', Logger_dir, Date, ExpStartTime)
        BatIDNameCorrect(Logger_dir, Date, ExpStartTime)
end

%% INTERNAL FUNCTION
function BatIDNameCorrect(Loggers_dir, Date, ExpStartTime)
Error=0; %counter for the number of vocalizations misslabeled
% Find the google drive folder
GGFolder = '/Users/elie/Google Drive/My Drive/';
if ~(exist(GGFolder, 'file')==7)
    GGFolder = '/Users/elie/Google Drive/Mon Drive/';
end
if ~(exist(GGFolder, 'file')==7)
    warning('cannot find GGFolder at %s\n', GGFolder)
    GGFolder = input('Please enter the GGFolder path: ', 's');
    keyboard
end

% Set the path to the recording log
if contains(Loggers_dir, 'LMC')
    Path2RecordingTable = fullfile(GGFolder,'/BatmanData/RecordingLogs/recording_logs.xlsx');
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:P200','basic');
    RowData = find((cell2mat(RecTableData(2:end,1))== str2double(Date))) +1;
elseif contains(Loggers_dir, 'Juvenile')
    Path2RecordingTable = fullfile(GGFolder,'/JuvenileRecordings/JuvenileRecordingsNWAF155_Log.xlsx');
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:Q123','basic');
    RowData = find((cell2mat(RecTableData(4:end,1))== str2double(['20' Date]))) +3;
elseif contains(Loggers_dir, 'DeafSaline')
    Path2RecordingTable = fullfile(GGFolder,'/JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx');
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:AE41','basic');
    RowData = find((cell2mat(RecTableData(2:end,1))== str2double(['20' Date]))) +1;
end
% Find the ID of the recorded bats
DataInfo = RecTableData(RowData,:);
Header = RecTableData(1,:);
BatIDCol = find(contains(Header, 'Bat'));
All_loggers_dir = dir(fullfile(Loggers_dir, '*ogger*'));
DirFlags = [All_loggers_dir.isdir];
% Extract only those that are directories.
All_loggers_dir = All_loggers_dir(DirFlags);
LoggerName = cell(length(All_loggers_dir),1);
BatID = cell(length(All_loggers_dir),1);
for ll=1:length(All_loggers_dir)
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
end


% Load data
Data1 = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractDat*.mat', Date, ExpStartTime)));
% select the correct files
Gdf = zeros(length(Data1),1);
for df=1:length(Data1)
    if length(strfind(Data1(df).name, '_'))==2
        Gdf(df)=1;
    end
end

if sum(Gdf)==0
    warning('No vocalization data extracted by who_calls.m or get_logger_data_voc.m')
else
    
    DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractDat*_*.mat', Date, ExpStartTime)));
    if length(DataFile)~=sum(Gdf)
        warning('The number of files generated by get_logger_data_voc (Data1) is not the same as the number generated by who_calls (DataFile)')
        keyboard
    end
    % Loop through the datafiles
    for df=1:sum(Gdf) %1
        % bringing the file back on the local computer (we're going to write
        % pretty often to it)
        Data1 = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData%d.mat', Date, ExpStartTime, df)));
        if isempty(Data1)
            Data1 = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)));
        end
        DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData%d_*.mat', Date, ExpStartTime,df)));
        if isempty(DataFile) % who calls was the earlier format
            DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date, ExpStartTime)));
        end
        fprintf(1,'Set %d/%d\nwith file %s and %s\n', df, sum(Gdf), Data1.name, DataFile.name)
        % Save the BatID and Logger names
        save(fullfile(Loggers_dir, DataFile.name), 'BatID', 'LoggerName', '-append');
        
        load(fullfile(Loggers_dir, DataFile.name),'BioSoundFilenames');
        if ~exist('BioSoundFilenames', 'var')
            continue
        end
        %% Loop through calls, correct wavfile and pdf names 
        % Turn off warning notifications for python 2 struct conversion
        NV = size(BioSoundFilenames,1);
        if isempty(BioSoundFilenames)
            continue
        end
        
        for NVocFile=1:NV
            OldBSName1 = BioSoundFilenames{NVocFile,1};
            OldBSName2=BioSoundFilenames{NVocFile,2};
            OldBSName = OldBSName1;
            if isempty(BioSoundFilenames{NVocFile,1})
                if isempty(BioSoundFilenames{NVocFile,2})
                    continue
                else
                    OldBSName = OldBSName2;
                end
            end
            ALInd = strfind(OldBSName, 'AL');
            ElmtInd = strfind(OldBSName, 'Elmt');
            ALNum = OldBSName((ALInd+2):(ElmtInd-2));

            % ID of the bat
            ALIndex = strcmp(LoggerName, ['AL' ALNum]);
            if sum(ALIndex)~=1
                keyboard
            end
            BatID_local =BatID{ALIndex};
            BatInd = strfind(OldBSName, 'Bat');
            BatID_old = OldBSName((BatInd+3):(ALInd-2));
            if ~strcmp(BatID_old, num2str(BatID_local))
                Error = Error+1;
                if ~isempty(OldBSName1)
                    BioSoundFilenames{NVocFile,1}((BatInd+3):(ALInd-2)) = num2str(BatID_local);
                    [Status,MSG]=movefile(OldBSName1, BioSoundFilenames{NVocFile,1});
                    if ~Status
                        warning('Issue with renaming file, error is:')
                        MSG %#ok<*NOPRT>
                        keyboard
                    end
                    OldPDFName1 = [OldBSName1(1:end-3) 'pdf'];
                    [Status,MSG]=movefile(OldPDFName1, [BioSoundFilenames{NVocFile,1}(1:end-3) 'pdf']);
                    if ~Status
                        warning('Issue with renaming file, error is:')
                        MSG
                        keyboard
                    end
                end
                
                BioSoundFilenames{NVocFile,2}((BatInd+3):(ALInd-2)) = num2str(BatID_local);
                [Status,MSG]=movefile(OldBSName2, BioSoundFilenames{NVocFile,2});
                if ~Status
                    warning('Issue with renaming file, error is:')
                    MSG
                    keyboard
                end
                OldPDFName2 = [OldBSName2(1:end-3) 'pdf'];
                [Status,MSG]=movefile(OldPDFName2, [BioSoundFilenames{NVocFile,2}(1:end-3) 'pdf']);
                if ~Status
                    warning('Issue with renaming file, error is:')
                    MSG
                    keyboard
                end
            end   
        end
        % save the values!
        save(fullfile(Loggers_dir, DataFile.name), 'BioSoundFilenames','-append');
        clear BioSoundFilenames
    end
end
fprintf(1,'\n**** %d vocalizations were mislabeled with the wrong bat ID ****\n', Error)
end


