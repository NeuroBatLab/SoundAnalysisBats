%% UMAP on MPS and some acoustic analysis of calls manually curated in detail!
% Playing around with all experimental dates manually checked so far

BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
QualityLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalQuality.txt');

% these are the list of files to analyse
if ~exist(QualityLog, 'file')
    error('Cannot find the list of file to anayse in: %s \n', QualityLog);
else
    FidQuality = fopen(QualityLog, 'r');
    Header = textscan(FidQuality,'%s\t%s\t%s\n',1);
    DoneListQuality = textscan(FidQuality,'%s\t%s\t%s');
    fclose(FidQuality);
end

% Find out the number of files to analyse
NToDo = length(DoneListQuality{1});
NF = nan(NToDo,1);
DataFiles = cell(1,NToDo);
for ff=1:NToDo
    % Get the BatID the Date and the time of the ff experiment 
    Date = DoneListQuality{2}{ff};
    ExpStartTime = DoneListQuality{3}{ff};
    Logger_dir = fullfile(BaseDataDir,['20' Date], 'audiologgers');
    DataFiles_local = dir(fullfile(Logger_dir, sprintf('%s_%s_VocExtractDat*_*.mat', Date, ExpStartTime)));
    NF(ff) = length(DataFiles_local);
    DataFiles{ff} = DataFiles_local';
end
DataFiles = [DataFiles{:}]';
%% 
% Gather data from all calls elmts

NF = length(DataFiles);

BatID = cell(NF*100,1);
CallType = cell(NF*100,1);
MicAudioGood = nan(NF*100,1);
Duration = nan(NF*100,1);
RMS_mean = nan(NF*100,1);
Sal_mean = nan(NF*100,1);
F0_mean = nan(NF*100,1);
Spect_mean = nan(NF*100,1);
Spect_std = nan(NF*100,1);
Spect_kurt = nan(NF*100,1);
Spect_skew = nan(NF*100,1);
Spect_ent = nan(NF*100,1);
Time_mean = nan(NF*100,1);
Time_std = nan(NF*100,1);
Time_kurt = nan(NF*100,1);
Time_skew = nan(NF*100,1);
Time_ent = nan(NF*100,1);
SpectralMean_mean = nan(NF*100,1);
Q1_mean = nan(NF*100,1);
Q2_mean = nan(NF*100,1);
Q3_mean = nan(NF*100,1);
MPS = cell(NF*100,1);
SpectMic_mean = nan(NF*100,1);
SpectMic_std = nan(NF*100,1);
SpectMic_kurt = nan(NF*100,1);
SpectMic_skew = nan(NF*100,1);
SpectMic_ent = nan(NF*100,1);
Q1Mic_mean = nan(NF*100,1);
Q2Mic_mean = nan(NF*100,1);
Q3Mic_mean = nan(NF*100,1);
VocCount = nan(NF,1);
% VocUID = nan(NF*100,1);
warning('off', 'MATLAB:Python:UnsupportedLoad')
for nf = 1:NF
    fprintf(1,'File %d/%d\n', nf, NF);
    DataFile = load(fullfile(DataFiles(nf).folder, DataFiles(nf).name));
    if isempty(DataFile.BioSoundCalls)
        fprintf(1,'No Data for this file %s\n', DataFiles(nf).name)
        continue
    end
    VVCount = 0;
    for vv=1:size(DataFile.BioSoundFilenames,1)
        if isempty(DataFile.BioSoundFilenames{vv,2}) % That was probably an noise that have not been kept
            continue
        end
        [~,Name,~] = fileparts(DataFile.BioSoundFilenames{vv,2});
        NC_local = sum(DataFile.ManualAnnotationOK{vv}==1);
        IndCall = find(DataFile.ManualAnnotationOK{vv}==1);
        if iscell(DataFile.ManualCallType{vv}) && sum(~cellfun('isempty',DataFile.ManualCallType{vv}))~= NC_local
            error('There is an issue with the identification of syllables')
        elseif ~iscell(DataFile.ManualCallType{vv}) && ~isempty(DataFile.ManualCallType{vv})~= NC_local
            error('There is an issue with the identification of syllables')
        end
        BatInd = strfind(Name,'Bat');
        for ee=1:NC_local
            cc = IndCall(ee);
            VVCount = VVCount +1;
            BatID{VVCount+sum(VocCount, 'omitnan')} = Name(BatInd+(3:7));
            CallType{VVCount+sum(VocCount, 'omitnan')} = DataFile.ManualCallType{vv}{cc};
            if ~isempty(DataFile.AudioGood{vv}) % Issue here! AudioGood was not saved for part of the data!
                MicAudioGood(VVCount+sum(VocCount, 'omitnan')) = DataFile.AudioGood{vv}(cc);
            else
                warning('Looks like AudioGood was not saved for %s vv=%d', DataFiles(nf).name, vv)
                MicAudioGood(VVCount+sum(VocCount, 'omitnan')) = 0;
            end
            Duration(VVCount+sum(VocCount, 'omitnan')) = DataFile.Duration(vv);
            RMS_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.RMS(vv);
            Sal_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.MeanSal(cc);
            F0_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.MeanF0(cc);
            Spect_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.meanspect(cc);
            Spect_std(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.stdspect(cc);
            Spect_kurt(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.kurtosisspect(cc);
            Spect_skew(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.skewspect(cc);
            Spect_ent(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.entropyspect(cc);
            Time_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.meantime(cc);
            Time_std(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.stdtime(cc);
            Time_kurt(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.kurtosistime(cc);
            Time_skew(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.skewtime(cc);
            Time_ent(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.entropytime(cc);
            Q1_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.MeanQ1t(cc);
            Q2_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.MeanQ2t(cc);
            Q3_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.MeanQ3t(cc);
            SpectralMean_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,2}.MeanSpectralMean(cc);
            if cc<=length(DataFile.BioSoundCalls{vv,1}.meanspect) % for some reason sometimes the last call is not prperly anayzed for microphone data under What_calls  -> Need to fix that at some point!!
                SpectMic_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,1}.meanspect(cc);
                SpectMic_std(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,1}.stdspect(cc);
                SpectMic_kurt(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,1}.kurtosisspect(cc);
                SpectMic_skew(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,1}.skewspect(cc);
                SpectMic_ent(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,1}.entropyspect(cc);
                Q1Mic_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,1}.MeanQ1t(cc);
                Q2Mic_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,1}.MeanQ2t(cc);
                Q3Mic_mean(VVCount+sum(VocCount, 'omitnan')) = DataFile.BioSoundCalls{vv,1}.MeanQ3t(cc);
            end
            if isfield(DataFile.BioSoundCalls{vv,2}, 'wf_elmts')
                MPSRows = round(DataFile.BioSoundCalls{vv,2}.wf_elmts{cc}*10^5)>=0;
                MPSCol = logical((round(DataFile.BioSoundCalls{vv,2}.wt_elmts{cc}/10)>=-30) .* (round(DataFile.BioSoundCalls{vv,2}.wt_elmts{cc}/10)<=30));
                MPS{VVCount+sum(VocCount, 'omitnan')} = reshape(DataFile.BioSoundCalls{vv,2}.mps_elmts{cc}(MPSRows,MPSCol),sum(MPSRows)*sum(MPSCol),1);
            else
                MPSRows = round(DataFile.BioSoundCalls{vv,2}.wf*10^5)>=0;
                MPSCol = logical((round(DataFile.BioSoundCalls{vv,2}.wt/10)>=-30) .* (round(DataFile.BioSoundCalls{vv,2}.wt/10)<=30));
                MPS{VVCount+sum(VocCount, 'omitnan')} = reshape(DataFile.BioSoundCalls{vv,2}.mps(MPSRows,MPSCol),sum(MPSRows)*sum(MPSCol),1);
            end
%             VocUID(VVCount+sum(VocCount, 'omitnan')) = nf*1000 + vv;
            if VVCount+sum(VocCount, 'omitnan')==1
                MPS_wf = DataFile.BioSoundCalls{vv,2}.wf(MPSRows);
                MPS_wt = DataFile.BioSoundCalls{vv,2}.wt(MPSCol);
            elseif (~isfield(DataFile.BioSoundCalls{vv,2}, 'wf_elmts')) && (any(round(MPS_wf*10^5) ~= round(DataFile.BioSoundCalls{vv,2}.wf(MPSRows)*10^5)) || any(round(MPS_wt/10) ~= round(DataFile.BioSoundCalls{vv,2}.wt(MPSCol)/10)))
                keyboard
            elseif (isfield(DataFile.BioSoundCalls{vv,2}, 'wf_elmts')) && (any(round(MPS_wf*10^5) ~= round(DataFile.BioSoundCalls{vv,2}.wf_elmts{cc}(MPSRows)*10^5)) || any(round(MPS_wt/10) ~= round(DataFile.BioSoundCalls{vv,2}.wt_elmts{cc}(MPSCol)/10)))
                keyboard
            end
        end

    end
    VocCount(nf) = VVCount;
    fprintf(1, 'NCalls = %d\n', VVCount)
end
warning('on', 'MATLAB:Python:UnsupportedLoad')
%%
NCalls = sum(VocCount, 'omitnan');
fprintf(1, 'Total number of detected and cut-out call elements: %d\n', NCalls)
BatID = BatID(1:NCalls);
CallType = CallType(1:NCalls);
MicAudioGood = MicAudioGood(1:NCalls);
Duration = Duration(1:NCalls);
RMS_mean = RMS_mean(1:NCalls);
Sal_mean = Sal_mean(1:NCalls);
F0_mean = F0_mean(1:NCalls);
Spect_mean = Spect_mean(1:NCalls);
Spect_std = Spect_std(1:NCalls);
Spect_kurt = Spect_kurt(1:NCalls);
Spect_skew = Spect_skew(1:NCalls);
Spect_ent = Spect_ent(1:NCalls);
Time_mean = Time_mean(1:NCalls);
Time_std = Time_std(1:NCalls);
Time_kurt = Time_kurt(1:NCalls);
Time_skew = Time_skew(1:NCalls);
Time_ent = Time_ent(1:NCalls);
SpectralMean_mean = SpectralMean_mean(1:NCalls);
Q1_mean = Q1_mean(1:NCalls);
Q2_mean = Q2_mean(1:NCalls);
Q3_mean = Q3_mean(1:NCalls);
SpectMic_mean = SpectMic_mean(1:NCalls);
SpectMic_std = SpectMic_std(1:NCalls);
SpectMic_kurt = SpectMic_kurt(1:NCalls);
SpectMic_skew = SpectMic_skew(1:NCalls);
SpectMic_ent = SpectMic_ent(1:NCalls);
Q1Mic_mean = Q1Mic_mean(1:NCalls);
Q2Mic_mean = Q2Mic_mean(1:NCalls);
Q3Mic_mean = Q3Mic_mean(1:NCalls);
MPS = [MPS{1:NCalls}]';
% VocUID = VocUID(1:sum(VocCount));
save(fullfile(BaseDataDir, 'Data4_DeafBats_CatCalls.mat'))
%% 
% Get color vectors ready

% Get the color vector ready for call type
ColorCode = [get(groot, 'DefaultAxesColorOrder'); 0 1 1; 0.5 0.5 0.5; 1 0 0 ; 0 1 0 ; 0 0 1; 0 0 0];
UCT = unique(CallType);
UCTFull = UCT;
UCTFull{strcmp(UCT, 'Tr')} = 'Trill';
UCTFull{strcmp(UCT, 'Ba')} = 'Bark';
UCTFull{strcmp(UCT, 'Pi')} = 'Pitchy Call';
UCTFull{strcmp(UCT, 'Bu')} = 'Low Buzz';
UCTFull{strcmp(UCT, 'Pa')} = 'Panting';
UCTFull{strcmp(UCT, 'LT')} = 'Low Tuck';
UCTFull{strcmp(UCT, 'Sq')} = 'Squeal';
UCTFull{strcmp(UCT, 'Ra')} = 'Rattle';
UCTFull{strcmp(UCT, 'Ch')} = 'Chuckles';
UCTFull{strcmp(UCT, 'BB')} = 'Bark Buzz';
UCTFull{strcmp(UCT, 'PB')} = 'Pitchy Call Buzz';
UCTFull{strcmp(UCT, 'SB')} = 'Squeal Buzz';
UCTFull{strcmp(UCT, 'Un')} = 'Unknown';
CCT = nan(length(CallType),3);
for ct=1:length(UCT)
    CCT(contains(CallType,UCT(ct)),:) = repmat(ColorCode(ct,:),sum(contains(CallType,UCT(ct))),1);
end
%%
% Get the color vector ready for Batname
UBat = unique(BatID);
CBat = nan(length(BatID),3);
for bat=1:length(UBat)
    CBat(contains(BatID,UBat(bat)),:) = repmat(ColorCode(bat,:),sum(contains(BatID,UBat(bat))),1);
end
%%
% Get the color vector ready for Sex and Deafness
Path2RecordingTable = '/Volumes/GoogleDrive/My Drive/JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx';
[~,~,RecTableData]=xlsread(Path2RecordingTable,2,'A1:k3','basic');
BatName = cell2mat(RecTableData(3,2:11));
CSexDeaf =  nan(length(BatID),3);
SexDeaf =  cell(length(BatID),1);
Sex =  cell(length(BatID),1);
Deaf = cell(length(BatID),1);
for bat = 1:length(BatName)
    ColID = 1+ bat;
    if strcmp(RecTableData(1,ColID), 'K') && strcmp(RecTableData(2,ColID), 'M')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(1,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'DM'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'M'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'D'};
    elseif strcmp(RecTableData(1,ColID), 'S') && strcmp(RecTableData(2,ColID), 'M')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(6,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'HM'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'M'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'H'};
    elseif strcmp(RecTableData(1,ColID), 'S') && strcmp(RecTableData(2,ColID), 'F')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(3,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'HF'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'F'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'H'};
    elseif strcmp(RecTableData(1,ColID), 'K') && strcmp(RecTableData(2,ColID), 'F')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(2,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'DF'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'F'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'D'};
    end

end
USexDeaf = {'Deaf Male' 'Hearing Male' 'Deaf Female' 'Hearing Female'};
UCSexDeaf = ColorCode([1 6 2 3], :);
%% 
% Construct the table of PAF for statistical purposes

Tbl =  table(BatID, CallType, Sex,Deaf,log10(Duration),RMS_mean,Sal_mean,F0_mean,Spect_mean,Spect_std,Spect_kurt,Spect_skew,Spect_ent,...
    Time_mean,Time_std,Time_kurt,Time_skew,Time_ent,SpectralMean_mean,Q1_mean,Q2_mean,Q3_mean,SpectMic_mean,SpectMic_std,...
    SpectMic_kurt,SpectMic_skew,SpectMic_ent,Q1Mic_mean,Q2Mic_mean,Q3Mic_mean, 'VariableNames', {'BatID', 'CT','Sex','HorD',...
    'Duration','RMS_mean','Sal_mean','F0_mean','Spect_mean','Spect_std','Spect_kurt','Spect_skew','Spect_ent',...
    'Time_mean','Time_std','Time_kurt','Time_skew','Time_ent','SpectralMean_mean','Q1_mean','Q2_mean','Q3_mean','SpectMic_mean','SpectMic_std',...
    'SpectMic_kurt','SpectMic_skew','SpectMic_ent','Q1Mic_mean','Q2Mic_mean','Q3Mic_mean'});
Tbl.HorD = categorical(Tbl.HorD, {'H' 'D'});
%% Explore the number of calls per bat

histogram(categorical(BatID))
xlabel('Bat')
ylabel('# calls')
%% Scatter of various acoustic features

tiledlayout(1,3)
nexttile
scatter(log10(Duration), RMS_mean,12,CCT,"filled")
xlabel('Duration (ms, log10 scale)')
ylabel('RMS')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
XLabels = cellfun(@str2num,(get(gca, 'XTickLabels')));
set(gca, 'XTickLabels', 10.^XLabels)

nexttile
scatter(log10(Duration), RMS_mean,12,CBat,"filled")
xlabel('Duration (ms, log10 scale)')
ylabel('RMS')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
XLabels = cellfun(@str2num,(get(gca, 'XTickLabels')));
set(gca, 'XTickLabels', 10.^XLabels)
nexttile
scatter(log10(Duration), RMS_mean,12,CSexDeaf,"filled")
xlabel('Duration (ms, log10 scale)')
ylabel('RMS')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
XLabels = cellfun(@str2num,(get(gca, 'XTickLabels')));
set(gca, 'XTickLabels', 10.^XLabels)
%% 
% *Stats for Duration*


lme = fitlme(Tbl, 'Duration ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Duration ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Duration ~ Sex + (1|BatID)')
lme4 = fitlme(Tbl, 'Duration ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Stats for RMS*

lme = fitlme(Tbl, 'RMS_mean ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'RMS_mean ~ Sex + HorD + (1|BatID)')
lme3 = fitlme(Tbl, 'RMS_mean ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'RMS_mean ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* Deaf animals significantly produce louder calls (p=0.022), Males 
% produce louder calls than females (p<10^-4)

tiledlayout(1,3)
nexttile
scatter(Sal_mean, F0_mean,12,CCT,"filled")
xlabel('Saliency')
ylabel('Fundamental (Hz)')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(Sal_mean, F0_mean,12,CBat,"filled")
xlabel('Saliency')
ylabel('Fundamental(Hz)')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(Sal_mean, F0_mean,12,CSexDeaf,"filled")
xlabel('Saliency')
ylabel('Fundamental(Hz)')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for F0*

lme = fitlme(Tbl, 'F0_mean ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'F0_mean ~ Sex + HorD + (1|BatID)')
lme3 = fitlme(Tbl, 'F0_mean ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'F0_mean ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Stats for Saliency*

lme = fitlme(Tbl, 'Sal_mean ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Sal_mean ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Sal_mean ~ Sex + (1|BatID)')
lme4 = fitlme(Tbl, 'Sal_mean ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* Deaf animals produce significantly higher pitch calls (p=0.044), 
% Males produce higher pitch and more pitchy calls than females (p<10^-2, p<10^-2)

tiledlayout(1,3)
nexttile
scatter(Spect_mean, Spect_std,12,CCT,"filled")
xlabel('Spectral mean')
ylabel('Spectral std')
YLim = get(gca, 'ylim');
XLim = [0 5000];
set(gca, 'xlim', XLim);
for ct=1:length(UCT)
    text(diff(XLim)*4/5+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(Spect_mean, Spect_std,12,CBat,"filled")
xlabel('Spectral mean')
ylabel('Spectral std')
YLim = get(gca, 'ylim');
XLim = [0 5000];
set(gca, 'xlim', XLim);
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(Spect_mean, Spect_std,12,CSexDeaf,"filled")
xlabel('Spectral mean')
ylabel('Spectral std')
YLim = get(gca, 'ylim');
XLim = [0 5000];
set(gca, 'xlim', XLim);
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Spectral mean*

lme = fitlme(Tbl, 'Spect_mean ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Spect_mean ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Spect_mean ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Spect_mean ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Stats for Spectral STD*

lme = fitlme(Tbl, 'Spect_std ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Spect_std ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Spect_std ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Spect_std ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* No effect of deafening or sex on the spectral mean and spectral 
% STD as recorded from the piezo

tiledlayout(1,3)
nexttile
scatter(Spect_ent, Spect_skew,12,CCT,"filled")
xlabel('Spectral entropy')
ylabel('Spectral skewness')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(Spect_ent, Spect_skew,12,CBat,"filled")
xlabel('Spectral entropy')
ylabel('Spectral skewness')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(Spect_ent, Spect_skew,12,CSexDeaf,"filled")
xlabel('Spectral entropy')
ylabel('Spectral skewness')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Spectral skewness*

lme = fitlme(Tbl, 'Spect_skew ~ Sex * HorD + (1|BatID)')
lme2 = fitlme(Tbl, 'Spect_skew ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Spect_skew ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Spect_skew ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Stats for Spectral entropy*

lme = fitlme(Tbl, 'Spect_ent ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Spect_ent ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Spect_ent ~ Sex + (1|BatID)')
lme4 = fitlme(Tbl, 'Spect_ent ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* Deafened and Male individuals produce calls with higher spectral 
% skweness, in particular Male deafened (HorD p<10^-3, Sex p<10-3, HorD*Sex p=0.05). 
% Male bats produce calls with lower spectral entropy (p<10^-2)

tiledlayout(1,3)
nexttile
scatter(Spect_ent, Spect_kurt,12,CCT,"filled")
xlabel('Spectral entropy')
ylabel('Spectral kurtosis')
YLim = [0 100];set(gca, 'ylim',YLim);
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(Spect_ent, Spect_kurt,12,CBat,"filled")
xlabel('Spectral entropy')
ylabel('Spectral kurtosis')
YLim = [0 100];set(gca, 'ylim', YLim);
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(Spect_ent, Spect_kurt,12,CSexDeaf,"filled")
xlabel('Spectral entropy')
ylabel('Spectral kurtosis')
YLim = [0 100];set(gca, 'ylim', YLim);
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Spectral kurtosis*

lme = fitlme(Tbl, 'Spect_kurt ~ Sex * HorD + (1|BatID)')
lme2 = fitlme(Tbl, 'Spect_kurt ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Spect_kurt ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Spect_kurt ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* Deafened and Male individuals produce calls with higher spectral 
% kurtosis, in particular Male deafened (HorD p<10^-2, Sex p<10-2, HorD*Sex p=0.03).

MicAudioGood = logical(MicAudioGood);
tiledlayout(1,3)
nexttile
scatter(SpectMic_mean(MicAudioGood), SpectMic_std(MicAudioGood),12,CCT(MicAudioGood,:),"filled")
xlabel('Mic Spectral mean')
ylabel('Mic Spectral std')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*4/5+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(SpectMic_mean(MicAudioGood), SpectMic_std(MicAudioGood),12,CBat(MicAudioGood,:),"filled")
xlabel('Mic Spectral mean')
ylabel('Mic Spectral std')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(SpectMic_mean(MicAudioGood), SpectMic_std(MicAudioGood),12,CSexDeaf(MicAudioGood,:),"filled")
xlabel('Mic Spectral mean')
ylabel('Mic Spectral std')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Mic Spectral mean*

TblMic = Tbl(MicAudioGood,:);
lme = fitlme(TblMic, 'SpectMic_mean ~ Sex * HorD + (1|BatID)')
lme2 = fitlme(TblMic, 'SpectMic_mean ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(TblMic, 'SpectMic_mean ~ Sex + (1|BatID)');
lme4 = fitlme(TblMic, 'SpectMic_mean ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Stats for Mic Spectral std*

lme = fitlme(TblMic, 'SpectMic_std ~ Sex * HorD + (1|BatID)')
lme2 = fitlme(TblMic, 'SpectMic_std ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(TblMic, 'SpectMic_std ~ Sex + (1|BatID)');
lme4 = fitlme(TblMic, 'SpectMic_std ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* Deafened and Male individuals and only them produce calls with 
% higher spectral std as recorded from the microphone (HorD*Sex p<10^-2).

tiledlayout(1,3)
nexttile
scatter(SpectMic_ent(MicAudioGood), SpectMic_skew(MicAudioGood),12,CCT(MicAudioGood,:),"filled")
xlabel('Mic Spectral entropy')
ylabel('Mic Spectral skewness')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(SpectMic_ent(MicAudioGood), SpectMic_skew(MicAudioGood),12,CBat(MicAudioGood,:),"filled")
xlabel('Mic Spectral entropy')
ylabel('Mic Spectral skewness')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(SpectMic_ent(MicAudioGood), SpectMic_skew(MicAudioGood),12,CSexDeaf(MicAudioGood,:),"filled")
xlabel('Mic Spectral entropy')
ylabel('Mic Spectral skewness')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Mic Spectral entropy*

lme = fitlme(TblMic, 'SpectMic_ent ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(TblMic, 'SpectMic_ent ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(TblMic, 'SpectMic_ent ~ Sex + (1|BatID)')
lme4 = fitlme(TblMic, 'SpectMic_ent ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Stats for Mic Spectral skewness*

lme = fitlme(TblMic, 'SpectMic_skew ~ Sex * HorD + (1|BatID)')
lme2 = fitlme(TblMic, 'SpectMic_skew ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(TblMic, 'SpectMic_skew ~ Sex + (1|BatID)');
lme4 = fitlme(TblMic, 'SpectMic_skew ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* Decrease spectral entropy for Males (p=0.03). Decrease spectral 
% skewness only for deafened male individuals (HorD*Sex p=0.018).

tiledlayout(1,3)
nexttile
scatter(SpectMic_ent(MicAudioGood), SpectMic_kurt(MicAudioGood),12,CCT(MicAudioGood,:),"filled")
xlabel('Mic Spectral entropy')
ylabel('Mic Spectral kurtosis')
YLim = [0 100];set(gca, 'ylim',YLim);
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(SpectMic_ent(MicAudioGood), SpectMic_kurt(MicAudioGood),12,CBat(MicAudioGood,:),"filled")
xlabel('Mic Spectral entropy')
ylabel('Mic Spectral kurtosis')
YLim = [0 100];set(gca, 'ylim', YLim);
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(SpectMic_ent(MicAudioGood), SpectMic_kurt(MicAudioGood),12,CSexDeaf(MicAudioGood,:),"filled")
xlabel('Mic Spectral entropy')
ylabel('Mic Spectral kurtosis')
YLim = [0 100];set(gca, 'ylim', YLim);
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Mic Spectral kurtosis*

lme = fitlme(TblMic, 'SpectMic_kurt ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(TblMic, 'SpectMic_kurt ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(TblMic, 'SpectMic_kurt ~ Sex + (1|BatID)')
lme4 = fitlme(TblMic, 'SpectMic_kurt ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* Increase spectral kurtosis for Males (p=0.021).

tiledlayout(1,3)
nexttile
scatter(Time_mean, Time_std,12,CCT,"filled")
xlabel('Temporal mean')
ylabel('Temporal std')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*4/5+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(Time_mean, Time_std,12,CBat,"filled")
xlabel('Temporal mean')
ylabel('Temporal std')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(Time_mean, Time_std,12,CSexDeaf,"filled")
xlabel('Temporal mean')
ylabel('Temporal std')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Temporal mean*

lme = fitlme(Tbl, 'Time_mean ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Time_mean ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Time_mean ~ Sex + (1|BatID)')
lme4 = fitlme(Tbl, 'Time_mean ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Stats for Temporal std*

lme = fitlme(Tbl, 'Time_std ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Time_std ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Time_std ~ Sex + (1|BatID)')
lme4 = fitlme(Tbl, 'Time_std ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* Decrease time mean and time std for Males (p=0.034; p=0.02) 
% -> shorter calls with sharper attack

tiledlayout(1,3)
nexttile
scatter(Time_ent, Time_kurt,12,CCT,"filled")
xlabel('Temporal entropy')
ylabel('Temporal kurtosis')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*5/6+XLim(1),YLim(2)-diff(YLim)*(ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(Time_ent, Time_kurt,12,CBat,"filled")
xlabel('Temporal entropy')
ylabel('Temporal kurtosis')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(Time_ent, Time_kurt,12,CSexDeaf,"filled")
xlabel('Temporal entropy')
ylabel('Temporal kurtosis')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Temporal entropy*

lme = fitlme(Tbl, 'Time_ent ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Time_ent ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Time_ent ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Time_ent ~ HorD + (1|BatID)')
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Stats for Temporal kurtosis*

lme = fitlme(Tbl, 'Time_kurt ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Time_kurt ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Time_kurt ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Time_kurt ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* Increase time entropy for deafened bats (p=0.016) -> less time 
% variation in power

tiledlayout(1,3)
nexttile
scatter(Time_ent, Time_skew,12,CCT,"filled")
xlabel('Temporal entropy')
ylabel('Temporal skewness')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*5/6+XLim(1),YLim(2)-diff(YLim)*(ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(Time_ent, Time_skew,12,CBat,"filled")
xlabel('Temporal entropy')
ylabel('Temporal skewness')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(Time_ent, Time_skew,12,CSexDeaf,"filled")
xlabel('Temporal entropy')
ylabel('Temporal skewness')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Temporal entropy*

lme = fitlme(Tbl, 'Time_ent ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Time_ent ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Time_ent ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Time_ent ~ HorD + (1|BatID)')
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Stats for Temporal skewness*

lme = fitlme(Tbl, 'Time_skew ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Time_skew ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Time_skew ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Time_skew ~ HorD + (1|BatID)')
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* Decrease time skewness and increase time entropy for Deaf bats 
% (p=0.033; p=0.016)

tiledlayout(1,3)
nexttile
scatter(SpectralMean_mean, Spect_mean,12,CCT,"filled")
xlabel('Spectral Mean')
ylabel('Spect mean')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*5/6+XLim(1),YLim(2)-diff(YLim)*(ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(SpectralMean_mean, Spect_mean,12,CBat,"filled")
xlabel('Spectral Mean')
ylabel('Spect mean')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(SpectralMean_mean, Spect_mean,12,CSexDeaf,"filled")
xlabel('Spectral Mean')
ylabel('Spect mean')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Spectral Mean*

lme = fitlme(Tbl, 'SpectralMean_mean ~ Sex * HorD + (1|BatID)')
lme2 = fitlme(Tbl, 'SpectralMean_mean ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'SpectralMean_mean ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'SpectralMean_mean ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* No effect of Mean spectral mean

tiledlayout(1,3)
nexttile
scatter(Q1_mean, Q2_mean,12,CCT,"filled")
xlabel('First Spectral Quartile')
ylabel('Spectral median')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*5/6+XLim(1),YLim(2)-diff(YLim)*(ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(Q1_mean, Q2_mean,12,CBat,"filled")
xlabel('First Spectral Quartile')
ylabel('Spectral median')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(Q1_mean, Q2_mean,12,CSexDeaf,"filled")
xlabel('First Spectral Quartile')
ylabel('Spectral median')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for First spectral quartile*

lme = fitlme(Tbl, 'Q1_mean ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Q1_mean ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Q1_mean ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Q1_mean ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Stats for Spectral median*

lme = fitlme(Tbl, 'Q2_mean ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Q2_mean ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Q2_mean ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Q2_mean ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* No Effect of Q1 and Q2

tiledlayout(1,3)
nexttile
scatter(Q3_mean, Q2_mean,12,CCT,"filled")
xlabel('Third Spectral Quartile')
ylabel('Spectral median')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UCT)
    text(diff(XLim)*5/6+XLim(1),YLim(2)-diff(YLim)*(ct)/30, UCTFull(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end

nexttile
scatter(Q3_mean, Q2_mean,12,CBat,"filled")
xlabel('Third Spectral Quartile')
ylabel('Spectral median')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(UBat)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, UBat(ct), 'Color',ColorCode(ct,:), 'FontWeight','bold' )
end
nexttile
scatter(Q3_mean, Q2_mean,12,CSexDeaf,"filled")
xlabel('Third Spectral Quartile')
ylabel('Spectral median')
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
for ct=1:length(USexDeaf)
    text(diff(XLim)*2/3+XLim(1),YLim(2)-diff(YLim)*(2+ct)/30, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
%% 
% *Stats for Third spectral quartile*

lme = fitlme(Tbl, 'Q3_mean ~ Sex * HorD + (1|BatID)');
lme2 = fitlme(Tbl, 'Q3_mean ~ Sex + HorD + (1|BatID)');
lme3 = fitlme(Tbl, 'Q3_mean ~ Sex + (1|BatID)');
lme4 = fitlme(Tbl, 'Q3_mean ~ HorD + (1|BatID)');
TestLME_Int=compare(lme2, lme)
TestLME_HD=compare(lme3, lme2)
TestLME_Sex=compare(lme4, lme2)
%% 
% *Conclusion:* No effect on third spectral quartile
%% Plot the average MPS for each call category
% First average MPS over all call

tiledlayout(1,1)
MPS_mean = reshape(mean(MPS)',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_mean, MPS_wf, MPS_wt);
title('Average MPS over all calls')
%% 
% Now average MPS for Trills

tiledlayout(1,3)
nexttile
IndH = logical(strcmp(CallType, 'Tr').*contains(SexDeaf, 'H')) ;
MPS_TrillH = reshape(mean(MPS(IndH,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_TrillH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title(sprintf('Normalized Trill Hearing bats (n=%d)', sum(IndH)))
nexttile
IndD = logical(strcmp(CallType, 'Tr').*contains(SexDeaf, 'D')) ;
MPS_TrillD = reshape(mean(MPS(IndD,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_TrillD./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title(sprintf('Normalized Trill Deaf bats (n=%d)', sum(IndD)))
nexttile
plot_mps2(MPS_TrillD./MPS_mean,MPS_TrillH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title('Normalized difference Trill Deaf-Hearing bats ')
%% 
% Now average MPS for Pitchy calls

tiledlayout(1,3)
nexttile
IndH = logical(strcmp(CallType, 'Pi').*contains(SexDeaf, 'H')) ;
MPS_PiH = reshape(mean(MPS(IndH,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_PiH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title(sprintf('Normalized Pitchy calls Hearing bats (n=%d)', sum(IndH)))
nexttile
IndD = logical(strcmp(CallType, 'Pi').*contains(SexDeaf, 'D')) ;
MPS_PiD = reshape(mean(MPS(IndD,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_PiD./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title(sprintf('Normalized Pitchy calls Deaf bats (n=%d)', sum(IndD)))
nexttile
plot_mps2(MPS_PiD./MPS_mean,MPS_PiH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title('Normalized difference Picthy Calls Deaf-Hearing bats ')
%% 
% Now average MPS for Barks

tiledlayout(1,3)
nexttile
IndH = logical(strcmp(CallType, 'Ba').*contains(SexDeaf, 'H'));
MPS_BaH = reshape(mean(MPS(IndH,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_BaH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title(sprintf('Normalized Barks Hearing bats (n=%d)', sum(IndH)))
nexttile
IndD = logical(strcmp(CallType, 'Ba').*contains(SexDeaf, 'D'));
MPS_BaD = reshape(mean(MPS(IndD,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_BaD./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title(sprintf('Normalized Barks Deaf bats (n=%d)', sum(IndD)))
nexttile
plot_mps2(MPS_BaD./MPS_mean,MPS_BaH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title('Normalized difference Barks Deaf-Hearing bats ')
%% 
% Now average MPS for Squeals

tiledlayout(1,3)
nexttile
IndH = logical(strcmp(CallType, 'Sq').*contains(SexDeaf, 'H'));
MPS_SqH = reshape(mean(MPS(IndH,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_SqH./MPS_mean, MPS_wf, MPS_wt,60, [-8 8]);
title(sprintf('Normalized Squeals Hearing bats (n=%d)', sum(IndH)))
nexttile
IndD = logical(strcmp(CallType, 'Sq').*contains(SexDeaf, 'D'));
MPS_SqD = reshape(mean(MPS(IndD,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_SqD./MPS_mean, MPS_wf, MPS_wt,60, [-8 8]);
title(sprintf('Normalized Squeals Deaf bats (n=%d)', sum(IndD)))
nexttile
plot_mps2(MPS_SqD./MPS_mean,MPS_SqH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title('Normalized difference Squeal Deaf-Hearing bats ')
%% 
% Now average MPS for Low Buzz

tiledlayout(1,3)
nexttile
IndH = logical(strcmp(CallType, 'Bu').*contains(SexDeaf, 'H'));
MPS_BuH = reshape(mean(MPS(IndH,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_BuH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title(sprintf('Normalized Low Buzz Hearing bats (n=%d)', sum(IndH)))
nexttile
IndD = logical(strcmp(CallType, 'Bu').*contains(SexDeaf, 'D'));
MPS_BuD = reshape(mean(MPS(IndD,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_BuD./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title(sprintf('Normalized Low Buzz Deaf bats (n=%d)', sum(IndD)))
nexttile
plot_mps2(MPS_BuD./MPS_mean,MPS_BuH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title('Normalized difference Buzz Deaf-Hearing bats ')
%% 
% Now average MPS for Panting

tiledlayout(1,3)
nexttile
IndH = logical(strcmp(CallType, 'Pa').*contains(SexDeaf, 'H'));
MPS_PaH = reshape(mean(MPS(IndH,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_PaH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title(sprintf('Normalized Panting Hearing Bats (n=%d)', sum(IndH)))
nexttile
IndD = logical(strcmp(CallType, 'Pa').*contains(SexDeaf, 'D'));
MPS_PaD = reshape(mean(MPS(IndD,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_PaD./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title(sprintf('Normalized Panting Deaf Bats (n=%d)', sum(IndD)))
nexttile
plot_mps2(MPS_PaD./MPS_mean,MPS_PaH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title('Normalized difference Panting Deaf-Hearing bats ')
%% 
% Now average MPS for Low Tuck

tiledlayout(1,3)
nexttile
IndH = logical(strcmp(CallType, 'LT').*contains(SexDeaf, 'H'));
MPS_LTH = reshape(mean(MPS(IndH,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_LTH./MPS_mean, MPS_wf, MPS_wt,60, [-8 8]);
title(sprintf('Normalized Low Tuck Hearing bats (n=%d)', sum(IndH)))
nexttile
IndD = logical(strcmp(CallType, 'LT').*contains(SexDeaf, 'D'));
MPS_LTD = reshape(mean(MPS(IndD,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_LTD./MPS_mean, MPS_wf, MPS_wt,60, [-8 8]);
title(sprintf('Normalized Low Tuck Deaf bats (n=%d)', sum(IndD)))
nexttile
plot_mps2(MPS_LTD./MPS_mean,MPS_LTH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title('Normalized difference Low Tucks Deaf-Hearing bats ')
%% 
% Now average MPS for Chuckles

tiledlayout(1,3)
nexttile
IndH = logical(strcmp(CallType, 'Ch').*contains(SexDeaf, 'H'));
MPS_ChH = reshape(mean(MPS(IndH,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_ChH./MPS_mean, MPS_wf, MPS_wt,60, [-8 8]);
title(sprintf('Normalized Chukles Hearing bats (n=%d)', sum(IndH)))
nexttile
IndD = logical(strcmp(CallType, 'Ch').*contains(SexDeaf, 'D'));
MPS_ChD = reshape(mean(MPS(IndD,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_ChD./MPS_mean, MPS_wf, MPS_wt,60, [-8 8]);
title(sprintf('Normalized Chukles Deaf bats (n=%d)', sum(IndD)))
nexttile
plot_mps2(MPS_ChD./MPS_mean,MPS_ChH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title('Normalized difference Chuckles Deaf-Hearing bats ')
%% 
% Now average MPS for Rattle

tiledlayout(1,3)
nexttile
IndH = logical(strcmp(CallType, 'Ra').*contains(SexDeaf, 'H'));
MPS_RaH = reshape(mean(MPS(IndH,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_RaH./MPS_mean, MPS_wf, MPS_wt,60, [-8 8]);
title(sprintf('Normalized Rattle Hearing bats (n=%d)', sum(IndH)))
nexttile
IndD = logical(strcmp(CallType, 'Ra').*contains(SexDeaf, 'D'));
MPS_RaD = reshape(mean(MPS(IndD,:))',length(MPS_wf), length(MPS_wt));
plot_mps(MPS_RaD./MPS_mean, MPS_wf, MPS_wt,60, [-8 8]);
title(sprintf('Normalized Rattle Deaf bats (n=%d)', sum(IndD)))
nexttile
plot_mps2(MPS_RaD./MPS_mean,MPS_RaH./MPS_mean, MPS_wf, MPS_wt,60,[-8 8]);
title('Normalized difference Rattle Deaf-Hearing bats ')
%% Project Normalized MPS in UMAP space
% First restrict the MPS spectral frequency to 3 cycles/kHz and normalized all 
% MPS

addpath(genpath('/Users/elie/Documents/CODE/GitHub/umapDistribution'))
IndWf = (MPS_wf*1000)<=3; % Max Spectral frequency taken into account in cycles/kHz
NVoc = size(MPS,1);
MPS_norm = cell(1,NVoc);
MPS_focus = cell(1,NVoc);
for vv=1:NVoc
    MPS_local = reshape(MPS(vv,:)',length(MPS_wf), length(MPS_wt));
    MPS_local = MPS_local(IndWf,:);
    MPS_norm{vv} = reshape(MPS_local./MPS_mean(IndWf,:), numel(MPS_local),1);
    MPS_focus{vv} = reshape(MPS_local, numel(MPS_local),1);
end
MPS_norm = [MPS_norm{:}]';
MPS_focus = [MPS_focus{:}]';
%% 
% You want to run these lines in the command window as the code open multiple 
% figures and that messes up with the livescripts

% default parameters
NNeigh = 30;
MinDist = 0.3;
[Reduction,UMAP]= run_umap(MPS_norm, 'metric','euclidean','n_component',3,'target_weight',0.55,'n_neighbors',NNeigh, 'min_dist',MinDist);
%%
figure()
scatter3(Reduction(~strcmp(CallType, 'Un'),1), Reduction(~strcmp(CallType, 'Un'),2),Reduction(~strcmp(CallType, 'Un'),3),20,CCT(~strcmp(CallType, 'Un'),:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of Normalized MPS (color = CallType) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
figure()
scatter3(Reduction(~strcmp(CallType, 'Un'),1), Reduction(~strcmp(CallType, 'Un'),2),Reduction(~strcmp(CallType, 'Un'),3),20,CBat(~strcmp(CallType, 'Un'),:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of Normalized MPS (color = BatID) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
% default parameters
NNeigh = 30;
MinDist = 0.3;
[Reduction,UMAP]= run_umap(MPS_focus, 'metric','euclidean','n_component',3,'target_weight',0.55,'n_neighbors',NNeigh, 'min_dist',MinDist);
%%
figure()
scatter3(Reduction(~strcmp(CallType, 'Un'),1), Reduction(~strcmp(CallType, 'Un'),2),Reduction(~strcmp(CallType, 'Un'),3),20,CCT(~strcmp(CallType, 'Un'),:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of MPS with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
% Decrease the number of neighbor to focalize more on local structure
NNeigh = 10;
MinDist = 0.3;
[Reduction,UMAP]= run_umap(MPS_norm, 'metric','euclidean','n_component',3,'target_weight',0.55,'n_neighbors',NNeigh, 'min_dist',MinDist);
%%
figure()
scatter3(Reduction(~strcmp(CallType, 'Un'),1), Reduction(~strcmp(CallType, 'Un'),2),Reduction(~strcmp(CallType, 'Un'),3),20,CCT(~strcmp(CallType, 'Un'),:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection or Normalized MPS with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
% Increase Min Dist to decrease density
NNeigh = 10;
MinDist = 0.5;
[Reduction,UMAP]= run_umap(MPS_norm, 'metric','euclidean','n_component',3,'target_weight',0.55,'n_neighbors',NNeigh, 'min_dist',MinDist);

%%
figure()
scatter3(Reduction(~strcmp(CallType, 'Un'),1), Reduction(~strcmp(CallType, 'Un'),2),Reduction(~strcmp(CallType, 'Un'),3),20,CCT(~strcmp(CallType, 'Un'),:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of normalized MPS with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
% Decrease Min Dist to increase density
NNeigh = 10;
MinDist = 0.051;
[Reduction,UMAP]= run_umap(MPS_norm, 'metric','euclidean','n_component',3,'target_weight',0.55,'n_neighbors',NNeigh, 'min_dist',MinDist);
%%
figure()
scatter3(Reduction(~strcmp(CallType, 'Un'),1), Reduction(~strcmp(CallType, 'Un'),2),Reduction(~strcmp(CallType, 'Un'),3),20,CCT(~strcmp(CallType, 'Un'),:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of Normalized MPS (Color = Call-Type) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
figure()
scatter3(Reduction(:,1), Reduction(:,2),Reduction(:,3),20,CBat,'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of Normalized MPS (Color = BatID) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
figure()
scatter3(Reduction(:,1), Reduction(:,2),Reduction(:,3),20,CSexDeaf,'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of Normalized MPS (Color = SexExp) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%% 
% *Explore the projection according to some acoustic features*
% 
% Pitch Saliency

scatter_cubehelix(Reduction, Sal_mean,2)
title(sprintf('UMAP Projection of Normalized MPS (Color = Pitch Saliency) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%% 
% Spectral entropy

scatter_cubehelix(Reduction, Spect_ent,2)
title(sprintf('UMAP Projection of Normalized MPS (Color = Spectral entropy) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%% 
% Pitch F0

scatter_cubehelix(Reduction, F0_mean,2,1,2)
title(sprintf('UMAP Projection of Normalized MPS (Color = Pitch F0) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%% 
% Time std

scatter_cubehelix(Reduction, Time_std*1000,2,1,10)
title(sprintf('UMAP Projection of Normalized MPS (Color = Time std in ms) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%% 
% Duration

scatter_cubehelix(Reduction, Duration,2,1,10)
title(sprintf('UMAP Projection of Normalized MPS (Color = duration in ms) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%% 
% RMS

scatter_cubehelix(Reduction, RMS_mean,2,2)
title(sprintf('UMAP Projection of Normalized MPS (Color = RMS) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
% Increase the number of neighbor to focalize more on global structure
NNeigh = 50;
MinDist = 0.051;
[Reduction,UMAP]= run_umap(MPS_norm, 'metric','euclidean','n_component',3,'target_weight',0.55,'n_neighbors',NNeigh, 'min_dist',MinDist);
%%
figure()
scatter3(Reduction(~strcmp(CallType, 'Un'),1), Reduction(~strcmp(CallType, 'Un'),2),Reduction(~strcmp(CallType, 'Un'),3),20,CCT(~strcmp(CallType, 'Un'),:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of Normalized MPS with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
% Try directly the MPS as an input instead of the normalized MPS
NNeigh = 10;
MinDist = 0.051;
[Reduction,UMAP]= run_umap(MPS_focus, 'metric','euclidean','n_component',3,'target_weight',0.55,'n_neighbors',NNeigh, 'min_dist',MinDist);
%%
figure()
scatter3(Reduction(~strcmp(CallType, 'Un'),1), Reduction(~strcmp(CallType, 'Un'),2),Reduction(~strcmp(CallType, 'Un'),3),20,CCT(~strcmp(CallType, 'Un'),:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of MPS with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%% Hierarchical clustering of vocalizations 
% Objective: find a partition of the acoustic space in hearing animals then 
% project the vocalizations of the deafened animals in the same space to classify 
% them accordig to the closest hearing bat partition. Within each partition explore 
% PAF between D and H bats. If there are still differences between them, then 
% quite robustly we can say that D bats have modified vocal production
% 
% First decide of an optimal representation for the MPS by using UMAP as a dimensionality 
% reduction
% 
% To inquire on the number of diemsions that might correctly capture the data, 
% let's run a PCA first

% 
[PC,Score,~, ~, VarExpl,~] = pca(MPS_norm(strcmp(Deaf, 'H'),:));
figure()
plot(cumsum(VarExpl), 'k-','LineWidth',2);
xlabel('PC#')
ylabel('% Cumulative Variance explained')
fprintf(1,'Number of PC to reach 90%% of variance = %d\n', find(cumsum(VarExpl)>90,1))
%%
fprintf(1,'******    Normalized MPS   ***********\n')
%%
[MPSReduction4Clustering,UMAP]= run_umap(MPS_norm(strcmp(Deaf, 'H'),:), 'metric','euclidean','n_component',50,'target_weight',0.55, 'n_neighbors',50, 'min_dist',0.051);
%%
% Only using MPS
Z = linkage( MPSReduction4Clustering,'ward','euclidean');
[D1, T1, ~] = dendrogram(Z,0, 'ColorThreshold','Default','Orientation','left');
set(D1,'LineWidth',2)
ClustSize = 2:9;
MeanSilhouette = nan(size(ClustSize));
MedianSilhouette = nan(size(ClustSize));
PercLowSilhouette = nan(size(ClustSize));
LowSi = 0.3;
for cc=1:length(ClustSize)
    T = cluster(Z,'maxclust',ClustSize(cc));
    [Si,h] = silhouette(MPSReduction4Clustering,T, 'Euclidean');
    MeanSilhouette(cc) = mean(Si);
    MedianSilhouette(cc) = median(Si);
    PercLowSilhouette(cc) = sum(Si<LowSi)/length(Si);
end
eva = evalclusters(MPSReduction4Clustering, 'linkage','CalinskiHarabasz','KList',1:15) %9 clusters
eva = evalclusters(MPSReduction4Clustering, 'linkage','DaviesBouldin','KList',1:15) % 2 clusters
eva = evalclusters(MPSReduction4Clustering, 'linkage','silhouette','KList',1:15) %4 clusters
%% 
% Get the UMAP projection in 3D and see how well the clustering partitioned 
% the space

NNeigh = 10;
MinDist = 0.051;
[ReductionH,UMAP]= run_umap(MPS_norm(strcmp(Deaf, 'H'),:), 'metric','euclidean','n_component',3,'target_weight',0.55,'n_neighbors',NNeigh, 'min_dist',MinDist);

%%
figure()
t=tiledlayout(1,2);
CTin9Cluster = nan(10,13);
nexttile
[D9, T9, ~] = dendrogram(Z,9, 'ColorThreshold','Default','Orientation','left');
set(D9,'LineWidth',2)
for cc=1:10
    for ct=1:length(UCT)
        CTin9Cluster(cc,ct) = sum(strcmp(CallType(T9==cc), UCT{ct}));
    end
end
nexttile
B = bar(CTin9Cluster, 'stacked');
B(8).FaceColor = ColorCode(8,:);
B(9).FaceColor = ColorCode(9,:);
B(10).FaceColor = [1 1 1];
legend(UCTFull, 'Location', 'northoutside')
%%
figure()
scatter3(ReductionH(:,1), ReductionH(:,2),ReductionH(:,3),20,ColorCode(T9,:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of MPS (Color = Acoustic clusters) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
figure()
scatter3(ReductionH(:,1), ReductionH(:,2),ReductionH(:,3),20,CCT(strcmp(Deaf, 'H'),:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of MPS (Color = Manual CT) with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
figure()
t=tiledlayout(1,2);
CTin4Cluster = nan(4,10);
nexttile
[D4, T4, ~] = dendrogram(Z,4, 'ColorThreshold','Default','Orientation','left');
set(D4,'LineWidth',2)
for cc=1:4
    for ct=1:length(UCT)
        CTin4Cluster(cc,ct) = sum(strcmp(CallType(T4==cc), UCT{ct}));
    end
end
nexttile
B= bar(CTin4Cluster, 'stacked');
B(8).FaceColor = ColorCode(8,:);
B(9).FaceColor = ColorCode(9,:);
B(10).FaceColor = [1 1 1];
legend(UCTFull, 'Location', 'northoutside')
%%
fprintf(1,'******    Raw MPS   ***********\n')
%%
[Reduction,UMAP]= run_umap(MPS_focus, 'metric','euclidean','n_component',50,'target_weight',0.55, 'n_neighbors',50, 'min_dist',0.051);
%%
% Obtain matrix of reduction MPS values for all vocal elements
MPSReduction = nan(length(VocUID), 50);
for vv=1:length(VocUID)
    MPSReduction(vv,:) = Reduction(UID==VocUID(vv),:);
end
%%
% Only using MPS
Z = linkage( MPSReduction,'ward','euclidean');
[D1, T1, ~] = dendrogram(Z,0, 'ColorThreshold','Default','Orientation','left');
set(D1,'LineWidth',2)
ClustSize = 2:9;
MeanSilhouette = nan(size(ClustSize));
MedianSilhouette = nan(size(ClustSize));
PercLowSilhouette = nan(size(ClustSize));
LowSi = 0.3;
for cc=1:length(ClustSize)
    T = cluster(Z,'maxclust',ClustSize(cc));
    [Si,h] = silhouette(MPSReduction,T, 'Euclidean');
    MeanSilhouette(cc) = mean(Si);
    MedianSilhouette(cc) = median(Si);
    PercLowSilhouette(cc) = sum(Si<LowSi)/length(Si);
end
eva = evalclusters(MPSReduction, 'linkage','CalinskiHarabasz','KList',1:30) %30 clusters
eva = evalclusters(MPSReduction, 'linkage','DaviesBouldin','KList',1:9) % 3 clusters
eva = evalclusters(MPSReduction, 'linkage','silhouette','KList',1:9) %3 clusters
%%
figure()
t=tiledlayout(1,2);
CTin9Cluster = nan(10,10);
nexttile
[D9, T9, ~] = dendrogram(Z,9, 'ColorThreshold','Default','Orientation','left');
set(D9,'LineWidth',2)
for cc=1:9
    for ct=1:length(UCT)
        CTin9Cluster(cc,ct) = sum(strcmp(CallType(T9==cc), UCT{ct}));
    end
end
nexttile
B = bar(CTin9Cluster, 'stacked');
B(8).FaceColor = ColorCode(8,:);
B(9).FaceColor = ColorCode(9,:);
B(10).FaceColor = [1 1 1];
legend(UCTFull, 'Location', 'northoutside')
%%
figure()
t=tiledlayout(1,2);
CTin4Cluster = nan(4,10);
nexttile
[D4, T4, ~] = dendrogram(Z,4, 'ColorThreshold','Default','Orientation','left');
set(D4,'LineWidth',2)
for cc=1:4
    for ct=1:length(UCT)
        CTin4Cluster(cc,ct) = sum(strcmp(CallType(T4==cc), UCT{ct}));
    end
end
nexttile
B= bar(CTin4Cluster, 'stacked');
B(8).FaceColor = ColorCode(8,:);
B(9).FaceColor = ColorCode(9,:);
B(10).FaceColor = [1 1 1];
legend(UCTFull, 'Location', 'northoutside')
%%
figure()
t=tiledlayout(1,2);
CTin3Cluster = nan(3,10);
nexttile
[D3, T3, ~] = dendrogram(Z,3, 'ColorThreshold','Default','Orientation','left');
set(D3,'LineWidth',2)
for cc=1:3
    for ct=1:length(UCT)
        CTin3Cluster(cc,ct) = sum(strcmp(CallType(T3==cc), UCT{ct}));
    end
end
nexttile
B= bar(CTin3Cluster, 'stacked');
B(8).FaceColor = ColorCode(8,:);
B(9).FaceColor = ColorCode(9,:);
B(10).FaceColor = [1 1 1];
legend(UCTFull, 'Location', 'northoutside')
%%
fprintf(1,'******    PAF    ***********\n')
DataTable = [F0_mean Sal_mean Spect_mean Spect_std Spect_ent SpectralMean_mean Time_ent Time_kurt];
DataTable1 = zscore(DataTable(~isnan(F0_mean),:));
Z = linkage( DataTable1,'ward','euclidean');
[D1, T1, ~] = dendrogram(Z,0, 'ColorThreshold','Default','Orientation','left');
set(D1,'LineWidth',2)
ClustSize = 2:9;
MeanSilhouette = nan(size(ClustSize));
MedianSilhouette = nan(size(ClustSize));
PercLowSilhouette = nan(size(ClustSize));
LowSi = 0.3;
for cc=1:length(ClustSize)
    T = cluster(Z,'maxclust',ClustSize(cc));
    [Si,h] = silhouette(DataTable1,T, 'Euclidean');
    MeanSilhouette(cc) = mean(Si);
    MedianSilhouette(cc) = median(Si);
    PercLowSilhouette(cc) = sum(Si<LowSi)/length(Si);
end
eva = evalclusters(DataTable1, 'linkage','CalinskiHarabasz','KList',1:9) %4 clusters
eva = evalclusters(DataTable1, 'linkage','DaviesBouldin','KList',1:9) % 4 clusters
eva = evalclusters(DataTable1, 'linkage','silhouette','KList',1:9) %4 clusters
%%
figure()
t=tiledlayout(1,2);
CTin3Cluster = nan(3,10);
nexttile
[D3, T3, ~] = dendrogram(Z,3, 'ColorThreshold','Default','Orientation','left');
set(D3,'LineWidth',2)
for cc=1:3
    for ct=1:length(UCT)
        CTin3Cluster(cc,ct) = sum(strcmp(CallType(T3==cc), UCT{ct}));
    end
end
nexttile
B= bar(CTin3Cluster, 'stacked');
B(8).FaceColor = ColorCode(8,:);
B(9).FaceColor = ColorCode(9,:);
B(10).FaceColor = [1 1 1];
legend(UCTFull, 'Location', 'northoutside')

%%
figure()
t=tiledlayout(1,2);
CTin4Cluster = nan(4,10);
nexttile
[D4, T4, ~] = dendrogram(Z,4, 'ColorThreshold','Default','Orientation','left');
set(D4,'LineWidth',2)
for cc=1:4
    for ct=1:length(UCT)
        CTin4Cluster(cc,ct) = sum(strcmp(CallType(T4==cc), UCT{ct}));
    end
end
nexttile
B= bar(CTin4Cluster, 'stacked');
B(8).FaceColor = ColorCode(8,:);
B(9).FaceColor = ColorCode(9,:);
B(10).FaceColor = [1 1 1];
legend(UCTFull, 'Location', 'northoutside')
%%
figure()
t=tiledlayout(1,2);
CTin5Cluster = nan(5,10);
nexttile
[D5, T5, ~] = dendrogram(Z,5, 'ColorThreshold','Default','Orientation','left');
set(D5,'LineWidth',2)
for cc=1:5
    for ct=1:length(UCT)
        CTin5Cluster(cc,ct) = sum(strcmp(CallType(T5==cc), UCT{ct}));
    end
end
nexttile
B= bar(CTin5Cluster, 'stacked');
B(8).FaceColor = ColorCode(8,:);
B(9).FaceColor = ColorCode(9,:);
B(10).FaceColor = [1 1 1];
legend(UCTFull, 'Location', 'northoutside')
%%
figure()
t=tiledlayout(1,2);
CTin9Cluster = nan(9,10);
nexttile
[D9, T9, ~] = dendrogram(Z,9, 'ColorThreshold','Default','Orientation','left');
set(D9,'LineWidth',2)
for cc=1:9
    for ct=1:length(UCT)
        CTin9Cluster(cc,ct) = sum(strcmp(CallType(T9==cc), UCT{ct}));
    end
end
nexttile
B= bar(CTin9Cluster, 'stacked');
B(8).FaceColor = ColorCode(8,:);
B(9).FaceColor = ColorCode(9,:);
B(10).FaceColor = [1 1 1];
legend(UCTFull, 'Location', 'northoutside')
%% 
% Try some UMAP projections of the PAF

NNeigh = 10;
MinDist = 0.051;
[Reduction,UMAP]= run_umap(DataTable1, 'metric','euclidean','n_component',3,'target_weight',0.55, 'n_neighbors',NNeigh, 'min_dist',MinDist);
%%
figure()
CallType_local = CallType(~isnan(F0_mean));
scatter3(Reduction(~strcmp(CallType_local, 'Un'),1), Reduction(~strcmp(CallType_local, 'Un'),2),Reduction(~strcmp(CallType_local, 'Un'),3),20,CCT(~strcmp(CallType_local, 'Un'),:,:),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
title(sprintf('UMAP Projection of PAF with Nneighbors = %d and MinDist = %.2f', NNeigh, MinDist))
%%
function [MPS4plot, Wt_local, Wf_local]=plot_mps(MPS, Wf, Wt, DBNOISE,CLim,Log, YLim, XLim, TwoDfilter)
    if nargin<4
        DBNOISE=60;
    end

    if nargin<6
        Log=1;
    end

    if nargin<7
        YLim = [0 max(Wf*10^3)];
    end
    if nargin<8
        XLim = [min(Wt) max(Wt)];
    end
    if nargin<5
        CLim = nan;
    end
    if nargin<9
        TwoDfilter=0;
    end
    
    Wf_i = logical((Wf*10^3>=YLim(1)).* (Wf*10^3<=YLim(2)));
    Wt_i = logical((Wt>=XLim(1)).* (Wt<=XLim(2)));
    MPS = MPS(Wf_i, Wt_i);
    Wf_local = Wf(Wf_i);
    Wt_local = Wt(Wt_i);
    if Log
        MPS4plot = 10*log10(MPS);
        MaxMPS = max(max(MPS4plot));
        MinMPS = MaxMPS-DBNOISE;
        MPS4plot(MPS4plot < MinMPS) = MinMPS;
    else
        MPS4plot = MPS; 
    end

    if TwoDfilter
        MPS4plot = imgaussfilt(MPS4plot,TwoDfilter);
    end
    
    %             imagesc(LogMPS, interpolation='nearest', aspect='auto', origin='lower', cmap=cmap, extent=ex)
    Im = imagesc(MPS4plot);
    axis xy
    colormap(Im.Parent,'jet');
    colorbar()
    xlabel('Temporal Frequency (Hz)')
    ylabel('Spectral Frequency (Cycles/kHz)')
    % get nice X and y tick labels
    MaxWf = max(floor(Wf_local*10^3));
    YTickLabel=0:MaxWf;
    YTick = nan(length(YTickLabel),1);
    for yy=1:length(YTick)
        YTick(yy) = find(floor(Wf_local*10^3)==YTickLabel(yy),1);
    end
    set(gca, 'YTickLabel', YTickLabel, 'YTick', YTick)
    MaxWt = max(floor(Wt_local*10^-2))*10^2;
    XTickLabel = -MaxWt:100:MaxWt;
    XTick = nan(length(XTickLabel),1);
    for xx=1:length(XTick)
        [~,XTick(xx)] = min(abs(round(Wt_local)-XTickLabel(xx)));
    end
    set(gca, 'XTickLabel',XTickLabel, 'XTick', XTick)
    if ~isnan(CLim)
        PlotMPS = gca;
        PlotMPS.CLim = CLim;
    end
end

function [MPS4plot, Wt_local, Wf_local]=plot_mps2(MPS1,MPS2, Wf, Wt, DBNOISE,CLim,Log, YLim, XLim, TwoDfilter)
    if nargin<5
        DBNOISE=60;
    end

    if nargin<7
        Log=1;
    end

    if nargin<8
        YLim = [0 max(Wf*10^3)];
    end
    if nargin<9
        XLim = [min(Wt) max(Wt)];
    end
    if nargin<6
        CLim = nan;
    end
    if nargin<10
        TwoDfilter=0;
    end
    
    Wf_i = logical((Wf*10^3>=YLim(1)).* (Wf*10^3<=YLim(2)));
    Wt_i = logical((Wt>=XLim(1)).* (Wt<=XLim(2)));
    MPS1 = MPS1(Wf_i, Wt_i);
    MPS2 = MPS2(Wf_i, Wt_i);
    Wf_local = Wf(Wf_i);
    Wt_local = Wt(Wt_i);
    if Log
        MPS4plot = 10*log10(MPS1)-10*log10(MPS2);
        MaxMPS = max(max(MPS4plot));
        MinMPS = MaxMPS-DBNOISE;
        MPS4plot(MPS4plot < MinMPS) = MinMPS;
    else
        MPS4plot = MPS1-MPS2; 
    end

    if TwoDfilter
        MPS4plot = imgaussfilt(MPS4plot,TwoDfilter);
    end
    
    %             imagesc(LogMPS, interpolation='nearest', aspect='auto', origin='lower', cmap=cmap, extent=ex)
    Im = imagesc(MPS4plot);
    axis xy
    colormap(Im.Parent,'jet');
    colorbar()
    xlabel('Temporal Frequency (Hz)')
    ylabel('Spectral Frequency (Cycles/kHz)')
    % get nice X and y tick labels
    MaxWf = max(floor(Wf_local*10^3));
    YTickLabel=0:MaxWf;
    YTick = nan(length(YTickLabel),1);
    for yy=1:length(YTick)
        YTick(yy) = find(floor(Wf_local*10^3)==YTickLabel(yy),1);
    end
    set(gca, 'YTickLabel', YTickLabel, 'YTick', YTick)
    MaxWt = max(floor(Wt_local*10^-2))*10^2;
    XTickLabel = -MaxWt:100:MaxWt;
    XTick = nan(length(XTickLabel),1);
    for xx=1:length(XTick)
        [~,XTick(xx)] = min(abs(round(Wt_local)-XTickLabel(xx)));
    end
    set(gca, 'XTickLabel',XTickLabel, 'XTick', XTick)
    if ~isnan(CLim)
        PlotMPS = gca;
        PlotMPS.CLim = CLim;
    end
end

function scatter_cubehelix(VAR, VARcol,x, MinColVal, LogToggle)
    % Getting VARZ ready for plotting.
    % Let's say you want to plot VARX VARY VARZ. If VARZ is always a non-null
    % integer then you can use it directly in cubehelix by specifying x=1
    % otherwise, the function first multiplies VARY by 10^x to make sure its
    % variability is maximized after ceiling it to an integer
    if nargin<3
        x=3;% here I'm doing *1000 as an example
    end

    if nargin<4
        MinColVal=0;
    end

    if nargin<5
        LogToggle=0;
    end
    % Make sure you don't have nan values
    VAR = VAR(~isnan(VARcol), :);
    VARcol = VARcol(~isnan(VARcol));
    if LogToggle==2
        VARcol = log2(VARcol);
    elseif LogToggle==10
        VARcol = log10(VARcol);
    end
    

    % Scale your minimum value as 0 or min Val
    Slide = 0;
    if ~MinColVal && sum(VARcol<=0)>0
        Slide = -min(VARcol)+ 10^(-x);
    elseif MinColVal
        Slide = -min(VARcol)+ 10^(-x);
    end

    
    
    VARcol_cube = ceil((VARcol+Slide)*10^x);
    %GRAD=cubehelix(max(VARZ_cube));%Note that you can choose the values of...
    ...your gradient by doing cubehelix_view and then specifying the values in...
        ...the function, for instance: GRAD=cubehelix(max(VARZ_cube),0.5, -0.8, 1.5, 1.7, [1,0]);
%     GRAD=cubehelix(max(VARcol_cube),0.5, -1.1, 1.5, 0.5, [1,0]);
    GRAD=cubehelix(max(VARcol_cube),0.5, -0.8, 1.5, 1.7, [1,0]);

    NU = length(VARcol_cube); % this is the number of points you have in your plot
    Colors = nan(NU, size(GRAD,2));
    for jj=1:NU
        Colors(jj,:) = GRAD(VARcol_cube(jj),:); % for each dot you call the line in the gradient of color that correspond to the value in Z
    end
    
    figure()
    scatter3(VAR(:,1), VAR(:,2),VAR(:,3),20,Colors,'filled')
    xlabel('Dim1')
    ylabel('Dim2')
    zlabel('Dim3')
    cc=colorbar();
    colormap(GRAD)
    
    YTL = get(cc, 'YTickLabel');
    YTL_new = max(VARcol_cube) .* str2double(YTL)/10^x - Slide;
    if LogToggle==2
        YTL_new = 2.^YTL_new;
    elseif LogToggle==10
        YTL_new = 10.^YTL_new;
    end
    YTL_new = num2str(round(YTL_new.*10^(x-1))./10^(x-1));
    set(cc, 'YTickLabel',YTL_new)%here you correct the value of the z axis that you artificially multiplied by 10^x

end
%% 
%