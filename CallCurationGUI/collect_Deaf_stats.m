% This script gather all the Deaf experiments
%% These are specific to the dataset and computer
BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
Today = datetime;
diary(fullfile(BaseDataDir, sprintf('StatsDiary_%s_%d%d.txt', date,Today.Hour,Today.Minute)));

DatesDir = dir(fullfile(BaseDataDir,'20*'));
NDates = length(DatesDir);
MinDur = 10;%(in min)


%% This first part go through all files and count all experiments
AllExpCount = 0;
ExpMissMicDataCount = 0;
ExpTooShortCount = 0;
ExpMissLogDataCount = 0;


for dd=1:NDates
    ParamFile = dir(fullfile(DatesDir(dd).folder, DatesDir(dd).name,'audio','*RecOnly_param.txt'));
    Nsessions = length(ParamFile);
    for nn=1:Nsessions
        fprintf(1, '\n\n\n Date: %s, experiment %d/%d, session %d/%d\n%s\n', DatesDir(dd).name,dd,NDates,nn,Nsessions,ParamFile(nn).name)
        Filepath = fullfile(ParamFile(nn).folder, ParamFile(nn).name);
        
        % check that the experiment has data!
        MicFiles = dir(fullfile(ParamFile(nn).folder, [ParamFile(nn).name(1:25) 'mic1*']));
        fid = fopen(Filepath);
        data = textscan(fid,'%s','Delimiter', '\t');
        fclose(fid);
        IndexLine = find(contains(data{1}, 'Task stops at'));
        if ~isempty(IndexLine)
            IndexChar = strfind(data{1}{IndexLine},'after');
            IndexChar2 = strfind(data{1}{IndexLine},'seconds');
            
            % find the data into that line
            Temp = round(str2double(data{1}{IndexLine}((IndexChar + 6):(IndexChar2-2)))/60);% duration in min
        else %No stop line, estimate the duration of the experiment by looking at the number of microphone files
            Temp = (length(MicFiles)-1)*10;
        end
        MicData = isempty(MicFiles);
        
        % Check that the experiment has logger data
        Logger_dir = [ParamFile(nn).folder(1:(strfind(ParamFile(nn).folder, 'audio')+4)) 'loggers'];
        All_loggers_dir = dir(fullfile(Logger_dir, '*ogger*'));
        LogData = isempty(All_loggers_dir);
        
        AllExpCount = AllExpCount +1;
        ExpMissMicDataCount = ExpMissMicDataCount + MicData;
        if ~MicData % There are mic data
            ExpTooShortCount = ExpTooShortCount + (Temp<=MinDur) ; % count if too short
            if Temp>MinDur % There are mic data and long enough
                ExpMissLogDataCount = ExpMissLogDataCount + LogData; % count if no Logger data
            end
        end
        
    end
end

fprintf(1,'Total # of experiments: %d\n', AllExpCount)
fprintf(1, 'Experiments with missing microphone data: %d/%d, %d%%\n', ExpMissMicDataCount,AllExpCount, round(ExpMissMicDataCount*100/AllExpCount))
fprintf(1, 'Experiments with Mic data that are too short: %d/%d, %d%%\n', ExpTooShortCount, AllExpCount-ExpMissMicDataCount, round(ExpTooShortCount*100/(AllExpCount-ExpMissMicDataCount)))
fprintf(1, 'Experiments with Mic data, longer that %d min that have no logger data: %d/%d, %d%%\n',MinDur,ExpMissLogDataCount, AllExpCount-ExpMissMicDataCount-ExpTooShortCount, round(ExpMissLogDataCount*100/(AllExpCount-ExpMissMicDataCount-ExpTooShortCount)))
TotCleanExp = AllExpCount-ExpMissMicDataCount-ExpTooShortCount-ExpMissLogDataCount;
fprintf(1, 'total number of clean experiments to process: %d\n', TotCleanExp)


%% Stats of number fo experiments that went through vocalization detection
BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
BaseCodeDir = 'C:\Users\BatLab\Documents\GitHub\';
ExpLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSal.txt'); % List of experiments for which vocalizations where detected on piezo
WhoLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWho.txt');% List of experiments on which manual curation was performed
AlliLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalAllignement.txt');

% Files that have their putative vocalizations extracted by
% result_reconly_Deafbats.m
FidExp = fopen(ExpLog, 'r');
Header = textscan(FidExp,'%s\t%s\t%s\t%s\t%s\n',1);
DoneListDetect = textscan(FidExp,'%s\t%s\t%s\t%.1f\t%d');
fclose(FidExp);

% Extract list of files that have been manually cured using CallCura
FidWho = fopen(WhoLog, 'r');
HeaderWho = textscan(FidWho,'%s\t%s\t%s\t%s\n',1);
DoneListWho = textscan(FidWho,'%s\t%s\t%s\t%d');
fclose(FidWho);

% Files that belong to the list of DoneListDetect and DoneListWho have
% their sequences extracted
ExtractedExp.Date = [DoneListDetect{2}; DoneListWho{2}];
ExtractedExp.Time = [DoneListDetect{3}; DoneListWho{3}];
% eliminate data duplicates
Dup = zeros(1,length(ExtractedExp.Date));
for ee=1:length(ExtractedExp.Date)
    Duplicate = find(contains(ExtractedExp.Date((ee+1):end), ExtractedExp.Date{ee}) .* contains(ExtractedExp.Time((ee+1):end), ExtractedExp.Time{ee}));
    if ~isempty(Duplicate)
        Dup(ee+Duplicate) = 1;
        Duplicate = [];
    end
end
ExtractedExp.Date(logical(Dup)) = [];
ExtractedExp.Time(logical(Dup)) = [];
fprintf(1, 'Experiments with piezo detected vocalizations: %d/%d, %d%%\n', length(ExtractedExp.Date), TotCleanExp, round(length(ExtractedExp.Date)/TotCleanExp*100))

% Experiments that have been manually curated
CuratedExp.Date = DoneListWho{2};
CuratedExp.Time = DoneListWho{3};
% eliminate data duplicates
Dup = zeros(1,length(CuratedExp.Date));
for ee=1:length(CuratedExp.Date)
    Duplicate = find(contains(CuratedExp.Date((ee+1):end), CuratedExp.Date{ee}) .* contains(CuratedExp.Time((ee+1):end), CuratedExp.Time{ee}));
    if ~isempty(Duplicate)
        Dup(ee+Duplicate) = 1;
        Duplicate = [];
    end
end
CuratedExp.Date(logical(Dup)) = [];
CuratedExp.Time(logical(Dup)) = [];
fprintf(1, 'Experiments with manually curated vocalizations: %d/%d, %d%%\n', length(CuratedExp.Date), TotCleanExp, round(length(CuratedExp.Date)/TotCleanExp*100))

% Add the current set being manually curated



% Gather the number of vocalizations manually curated for each experiment
CuratedExp.NumVoc = zeros(length(CuratedExp.Date),1);
CuratedExp.NumSeq = zeros(length(CuratedExp.Date),1);
CuratedExp.NumFullSeq = zeros(length(CuratedExp.Date),1); % Sequences with vocalizations
CuratedExp.VocPerBat = cell(length(CuratedExp.Date),1);
CuratedExp.BatID = cell(length(CuratedExp.Date),1);
CuratedExp.BatStatus = cell(length(CuratedExp.Date),1);
CuratedExp.UniqueSorterNames = {};
CuratedExp.SorterSeqNum = [];
CuratedExp.UniqueBatNames = [];
CuratedExp.BatVocNum = [];

for ee=1:length(CuratedExp.Date)
    ManuFiles = dir(fullfile(BaseDataDir,['20' CuratedExp.Date{ee}], 'audiologgers', sprintf('%s_%s_VocExtractData*_*', CuratedExp.Date{ee}, CuratedExp.Time{ee})));
    
        
    if ~strcmp(ManuFiles(1).name(27), '_')
        % re-order files to Find the first file that should contain bat ID and loggerID
        IndFile = nan(length(ManuFiles),1);
        for ff=1:length(ManuFiles)
            IndData = strfind(ManuFiles(ff).name, 'Data')+4;
            Ind_ = strfind(ManuFiles(ff).name, '_')-1;
            IndFile(ff) = str2double(ManuFiles(ff).name(IndData:Ind_(end)));
        end
        [~,OrdInd ] = sort(IndFile);
        ManuFiles=ManuFiles(OrdInd);
    else
        ManuFiles = dir(fullfile(BaseDataDir,['20' CuratedExp.Date{ee}], 'audiologgers', sprintf('%s_%s_VocExtractData_*', CuratedExp.Date{ee}, CuratedExp.Time{ee})));
    end
    for ff=1:length(ManuFiles)
        load(fullfile(ManuFiles(ff).folder, ManuFiles(ff).name), 'IndVocStartRaw_merged', 'SorterName')
        if ff==1
            load(fullfile(ManuFiles(ff).folder, ManuFiles(ff).name), 'BatID', 'LoggerName')
            if ~exist('BatID', 'var')
                keyboard
            end
            UniqueBatNames_local = unique([cell2mat(BatID); CuratedExp.UniqueBatNames]);
            CallCount_local = zeros(length(UniqueBatNames_local),1);
            if ~isempty(CuratedExp.UniqueBatNames) % Get previous data
                for bb=1:length(UniqueBatNames_local)
                    if sum(CuratedExp.UniqueBatNames == UniqueBatNames_local(bb))
                        CallCount_local(bb) = CuratedExp.BatVocNum(CuratedExp.UniqueBatNames == UniqueBatNames_local(bb));
                    end
                end
            end
            CuratedExp.UniqueBatNames = UniqueBatNames_local;
        end
        if exist('IndVocStartRaw_merged', 'var')
            CuratedExp.NumSeq(ee) = length(IndVocStartRaw_merged)+CuratedExp.NumSeq(ee);
            NumCall = nan(length(IndVocStartRaw_merged),1);
            for cc=1:length(IndVocStartRaw_merged)
                if ~isempty(IndVocStartRaw_merged{cc})
                    NumCall(cc)=length([IndVocStartRaw_merged{cc}{:}]);
                    for aa=1:length(BatID)
                        IndexBat = find(CuratedExp.UniqueBatNames == BatID{aa});
                        CallCount_local(IndexBat) = CallCount_local(IndexBat) + length(IndVocStartRaw_merged{cc}{aa});
                    end
                else
                    NumCall(cc) = 0;
                end
            end
            CuratedExp.NumFullSeq(ee) = sum(NumCall>0)+CuratedExp.NumFullSeq(ee);
            CuratedExp.NumVoc(ee) = sum(NumCall)+ CuratedExp.NumVoc(ee);
            clear IndVocStartRaw_merged
        else
            Ind_ = strfind(ManuFiles(ff).name, '_');
            load(fullfile(ManuFiles(ff).folder, [ManuFiles(ff).name(1:(Ind_(end)-1)) '.mat']), 'VocFilename')
            CuratedExp.NumMissedSeq(ee) = length(VocFilename);
            clear VocFilename
        end
        if exist('SorterName', 'var')
            USorterName = unique([SorterName CuratedExp.UniqueSorterNames]);
            SorterNumSeq = nan(length(USorterName),1);
            for sn=1:length(USorterName)
                SorterNumSeq(sn) = CuratedExp.SorterSeqNum(contains(CuratedExp.UniqueSorterNames, USorterName{sn})) + sum(contains(SorterName,USorterName{sn}));
            end
            CuratedExp.SorterSeqNum = SorterNumSeq;
            CuratedExp.UniqueSorterNames = USorterName;
            clear SorterName
        end
    end
    CuratedExp.BatVocNum = CallCount_local;
    clear BatID LoggerName
end
fprintf(1, 'Total number of manually curated sequences %d\n',sum(CuratedExp.NumSeq))
fprintf(1, 'Total number of sequences with vocalizations %d/%d, %d%%\n', sum(CuratedExp.NumFullSeq),sum(CuratedExp.NumSeq),round(sum(CuratedExp.NumFullSeq)*100/sum(CuratedExp.NumSeq)))
fprintf(1, 'Total number of vocalizations %d\n', sum(CuratedExp.NumVoc))
diary OFF
%% Some figures
BatStatus.name = [11648 14461 14463 14464 65696 71043 71047 71351 71353 71354];
BatStatus.sex = {'F' 'M' 'F' 'M' 'F' 'M' 'F' 'M' 'F' 'F'};
BatStatus.deaf = [0 0 1 0 0 1 1 1 0 1];

% Plot num Voc per batID
figure()
subplot(1,2,1)
BAR = bar(CuratedExp.BatVocNum);
ylabel('# Vocalizations')
xlabel('Bat Name')
BAR.Parent.XTickLabel = CuratedExp.UniqueBatNames;
BAR.FaceColor = 'flat';
for bb=1:length(CuratedExp.UniqueBatNames)
    if strcmp(BatStatus.sex(BatStatus.name==CuratedExp.UniqueBatNames(bb)), 'F')
        BAR.CData(bb,1) = 1;
    end
end

subplot(1,2,2)
BAR = bar(CuratedExp.BatVocNum);
ylabel('# Vocalizations')
xlabel('Bat Name')
BAR.Parent.XTickLabel = CuratedExp.UniqueBatNames;
BAR.FaceColor = 'flat';
for bb=1:length(CuratedExp.UniqueBatNames)
    if BatStatus.deaf(BatStatus.name==CuratedExp.UniqueBatNames(bb))
        BAR.CData(bb,:) = [1 1 1];
    else
        BAR.CData(bb,:) = [0 0 0];
    end
end


% Plot number of sequences per curator
figure()
BAR = bar(CuratedExp.SorterSeqNum);
ylabel('# sequences')
xlabel('Curator')
BAR.Parent.XTickLabel = CuratedExp.UniqueSorterNames;








