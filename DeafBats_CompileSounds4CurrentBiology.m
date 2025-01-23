BaseDataDir = '/Volumes/server_home-1/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
LocalDataDir = '/Users/elie/Documents/DeafBats/Data';
BoxPath = '/Users/elie/Box';
Path2Paper = fullfile(BoxPath, 'BatmanData', 'Deaf Paper');

% Loading PAF data
load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls.mat'), 'CallDate',  'MicFileNum', 'MicFileSampleID', 'MicAudioGood','OnSetOffsetSample',...
    'BatID','Time_mean','Q1Mic_mean')
MicAudioGood01 = MicAudioGood;
MicAudioGood01(isnan(MicAudioGood01)) = 0;
MicAudioGood01 = logical(MicAudioGood01);
fprintf(1,'In total, there are %d Vocal elements, %d are also of good microphone quality\n', length(MicAudioGood01), sum(MicAudioGood01))

% Get the color vector ready for call type
ColorCode = [get(groot, 'DefaultAxesColorOrder'); 0 1 1; 0.5 0.5 0.5; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0 1; 0 0 0];


% Get the vector ready for BatName, Sex and Deafness
Path2RecordingTable = fullfile(BoxPath,'JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx');
[~,~,RecTableData]=xlsread(Path2RecordingTable,2,'A1:k4','basic');
BatName = cell2mat(RecTableData(3,2:11));
BatSexDeaf = cell(size(BatName));
BatIDPaper = RecTableData(4,2:11);
CSexDeaf =  nan(length(BatID),3);
SexDeaf =  cell(length(BatID),1);
Sex =  cell(length(BatID),1);
Deaf = cell(length(BatID),1);
CBat = nan(length(BatID),3);
UCBatName = nan(length(BatName),3);
for bat = 1:length(BatName)
    ColID = 1+ bat;
    CBat(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(bat,:),sum(contains(BatID,num2str(BatName(bat)))),1);
    UCBatName(bat,:) = ColorCode(bat,:);
    if strcmp(RecTableData(1,ColID), 'K') && strcmp(RecTableData(2,ColID), 'M')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(1,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'DM'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'M'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'D'};
        BatSexDeaf{bat} = 'Kanamycin Male';
    elseif strcmp(RecTableData(1,ColID), 'S') && strcmp(RecTableData(2,ColID), 'M')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(6,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'HM'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'M'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'H'};
        BatSexDeaf{bat} = 'Saline Male';
    elseif strcmp(RecTableData(1,ColID), 'S') && strcmp(RecTableData(2,ColID), 'F')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(3,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'HF'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'F'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'H'};
        BatSexDeaf{bat} = 'Saline Female';
    elseif strcmp(RecTableData(1,ColID), 'K') && strcmp(RecTableData(2,ColID), 'F')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(2,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'DF'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'F'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'D'};
        BatSexDeaf{bat} = 'Kanamycin Female';
    end

end
USexDeaf = {'Kanamycin Male' 'Saline Male' 'Kanamycin Female' 'Saline Female'};
UCSexDeaf = ColorCode([1 6 2 3], :);
BatName = BatName([1:7 10 8:9]);
BatSexDeaf = BatSexDeaf([1:7 10 8:9]);
UCBatName = UCBatName([1:7 10 8:9],:);
UCBatNameSexDeaf = ColorCode([2 2 2 3 3 3 1 1 6 6], :);
BatIDPaper = BatIDPaper([1:7 10 8:9]);

load('WritingSounds.mat',"Sound_count","BioSoundFileNamesList","AudioFileName", "OnSetOffsetSample_Mic", "OnSetOffsetSample_AL", "OnSetOffsetSample_Good", "IndGood") % I erase that file as it was for tracking progress only
tStart = tic;
for ii=1:length(IndGood)
    if rem(ii,100)==0
        toc(tStart)
        fprintf(1, 'Vocalization %d\n',ii)
        save(fullfile(Path2Paper,'Data','WritingSounds.mat'), "ii", "Sound_count","BioSoundFileNamesList","AudioFileName", "OnSetOffsetSample_Mic", "OnSetOffsetSample_AL", "OnSetOffsetSample_Good", "IndGood")
        tStart = tic;
    end
    % Hearing Males examples (Low DFA1)
    IndEx = IndGood(ii);
    WhoFile = dir(fullfile(BaseDataDir, ['20' num2str(CallDate(IndEx))], 'audiologgers', sprintf('DeSa_%s_*_%d_%d*spec_200.pdf', num2str(CallDate(IndEx)),  MicFileNum(IndEx), MicFileSampleID(IndEx))));
    if isempty(WhoFile)
        WhoFile = dir(fullfile(BaseDataDir, ['20' num2str(CallDate(IndEx))], 'audiologgers', sprintf('*eSa_%s_*_%d_%d*_200.fig', num2str(CallDate(IndEx)),  MicFileNum(IndEx), MicFileSampleID(IndEx))));
        IndSlash = strfind(WhoFile.name, '\');
        WhoFile.name = WhoFile.name(IndSlash(end)+1:end);
    end
    Ind_ = strfind(WhoFile.name, '_');
    SetNum = WhoFile.name((Ind_(7)+1):(Ind_(8)-1));
    ExpTime = WhoFile.name((Ind_(2)+1):(Ind_(3)-1));
    DataFilename = dir(fullfile(BaseDataDir, ['20' num2str(CallDate(IndEx))], 'audiologgers', sprintf('%d_%s_VocExtractData%s_200.mat', CallDate(IndEx),ExpTime, SetNum)));
    DataFile = load(fullfile(DataFilename.folder, DataFilename.name));
    if any(cellfun(@isempty,DataFile.BioSoundFilenames(:,2)))
        IndBS = find(~cellfun(@isempty,DataFile.BioSoundFilenames(:,2)));
        CallsInd = IndBS(contains(DataFile.BioSoundFilenames(IndBS,2), sprintf('DeSa_%s_%s_voc_%d_%d', num2str(CallDate(IndEx)), ExpTime,  MicFileNum(IndEx), MicFileSampleID(IndEx))));
    else
        CallsInd = find(contains(DataFile.BioSoundFilenames(:,2), sprintf('DeSa_%s_%s_voc_%d_%d', num2str(CallDate(IndEx)), ExpTime,  MicFileNum(IndEx), MicFileSampleID(IndEx))));
    end
    CallInd = [];
    for vv=1:length(CallsInd)
        if isfield(DataFile.BioSoundCalls{CallsInd(vv),1}, 'MeanQ1t')
            if ~isnan(Q1Mic_mean(IndEx)) && any(DataFile.BioSoundCalls{CallsInd(vv),1}.MeanQ1t == Q1Mic_mean(IndEx))
                if any(DataFile.BioSoundCalls{CallsInd(vv),2}.meantime == Time_mean(IndEx)) % check how time mean is coded in Biosounds
                    CallInd = CallsInd(vv);
                    break
                else
                    keyboard
                end
            elseif isnan(Q1Mic_mean(IndEx)) && any(isnan(DataFile.BioSoundCalls{CallsInd(vv),1}.MeanQ1t))
                if any(DataFile.BioSoundCalls{CallsInd(vv),2}.meantime == Time_mean(IndEx)) % check how time mean is coded in Biosounds
                    CallInd = CallsInd(vv);
                    break
                else
                    keyboard
                end
            end
        end
    end
    if isempty(CallInd)
        keyboard
    end
    [~,Name,~] = fileparts(DataFile.BioSoundFilenames{CallInd,2});
    
    BatInd = strfind(Name,'Bat');
    % check the BatID 
    if ~strcmp(BatID{IndEx} , Name(BatInd+(3:7)))
        keyboard
    end

    % define the BatID same as in the paper
    BatIDPaper_local = BatIDPaper{BatName==str2double(BatID{IndEx})};
    
    if ii==1 || (~any(contains(BioSoundFileNamesList(1:(ii-1)), Name))) % this sound has not yet been saved
        Sound_count=Sound_count+1;
        SoundMic = DataFile.BioSoundCalls{CallInd,1}.sound;
        SoundAL = DataFile.BioSoundCalls{CallInd,2}.sound;
        SoundNameMic = fullfile(Path2Paper,'Data','VocDataBase',sprintf('%s_%s_%s_%s_%d_Mic.wav', BatIDPaper_local, Deaf{IndEx}, num2str(CallDate(IndEx)), ExpTime, Sound_count));
        if any(abs(SoundMic)>1)
            SoundMic = SoundMic ./ max(abs(SoundMic));
        end
        audiowrite(SoundNameMic,SoundMic,DataFile.BioSoundCalls{CallInd,1}.samprate)
        SoundNameAL = fullfile(Path2Paper,'Data','VocDataBase',sprintf('%s_%s_%s_%s_%d_AL.wav', BatIDPaper_local, Deaf{IndEx}, num2str(CallDate(IndEx)), ExpTime, Sound_count));
        if any(abs(SoundAL)>1)
            SoundAL = SoundAL ./ max(abs(SoundAL));
        end
        audiowrite(SoundNameAL,SoundAL,DataFile.BioSoundCalls{CallInd,2}.samprate)
        AudioFileName{ii,2}=SoundNameAL;
        AudioFileName{ii,1}=SoundNameMic;
        BioSoundFileNamesList{ii} = Name;
    else
        jj = find(contains(BioSoundFileNamesList(1:(ii-1)), Name));
        AudioFileName{ii,1}=AudioFileName{jj,1};
        AudioFileName{ii,2}=AudioFileName{jj,2};
        BioSoundFileNamesList{ii} = BioSoundFileNamesList{jj};
    end
    OnSetOffsetSample_Mic(ii,:) = OnSetOffsetSample_Good(ii,:);
    if OnSetOffsetSample_Mic(ii,1) == 0
        OnSetOffsetSample_Mic(ii,1) =1;
    end
    OnSetOffsetSample_AL(ii,:) = ceil(OnSetOffsetSample_Mic(ii,:) ./ DataFile.BioSoundCalls{CallInd,1}.samprate .* DataFile.BioSoundCalls{CallInd,2}.samprate);
end

save(fullfile(Path2Paper,'Data','WritingSounds.mat'), "ii", "Sound_count","BioSoundFileNamesList","AudioFileName", "OnSetOffsetSample_Mic", "OnSetOffsetSample_AL", "OnSetOffsetSample_Good", "IndGood")