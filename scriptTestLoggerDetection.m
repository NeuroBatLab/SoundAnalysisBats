%% Load input data
% Path2Data = '/Volumes/server_home/users/JulieE/LMC_CoEd/logger/20190623';
% Path2Audio = '/Volumes/server_home/users/JulieE/LMC_CoEd/audio/20190623';
% % Load data manually extracted
% load(fullfile(Path2Data, '190623_1401_VocExtractData_200.mat'))
% load(fullfile(Path2Data, '190623_1401_VocExtractData.mat'))
% load(fullfile(Path2Audio, '190623_1401_VocExtractTimes.mat'))


%% Load data manually extracted
Path2Data1 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audiologgers';
% Path2Data1 = '/Volumes/server_home/users/JulieE/JuvenileRecordings155/20190927/audiologgers';
load(fullfile(Path2Data1,'190927_1014_VocExtractData_200.mat'))
load(fullfile(Path2Data1,'190927_1014_VocExtractData.mat'))
% Path2Data2 = '/Volumes/server_home/users/JulieE/JuvenileRecordings155/20190927/audio';
Path2Data2 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audio';
load(fullfile(Path2Data2,'190927_1014_VocExtractTimes.mat'))

Fs_env=1000; % in Hertz, should have been saved in who_calls.m to correctly convert in time the starting and ending indices of vocalizations in IndVocStart
FS_Piezo = 50000; % could also be retrieved from Piezo_FS

%% Run the automatic detection
AllLoggers = dir(fullfile(Path2Data1, '*ogger*'));
DirFlags = [AllLoggers.isdir];
% Extract only those that are directories.
AllLoggers = AllLoggers(DirFlags);
NL = length(AllLoggers);
for ll=1:NL
    Data_directory = fullfile(AllLoggers(ll).folder,AllLoggers(ll).name, 'extracted_data');
    [SoundEvent_LoggerSamp.(sprintf('L%s',AllLoggers(ll).name(2:end))),SoundEvent_TranscTime_ms.(sprintf('L%s',AllLoggers(ll).name(2:end))),~,~] = piezo_find_calls_logger(Data_directory);
end
save(fullfile(Path2Data1, 'SoundEvent.mat'),'SoundEvent_LoggerSamp','SoundEvent_TranscTime_ms')
% load(fullfile(Path2Data1, 'SoundEvent.mat')



% Duration of the data manually analysed 100 min
% ManDur = 10*60*10^3; %in ms
ManDur = 10*10*60*10^3; %in ms

%% Get in transceiver time the onset/offset of each vocalization manually extracted for each logger
AL_ManId = fieldnames(Piezo_wave); % Names of the audioLoggers
NLoggers = length(AudioLogs);% number of audio loggers
ManCallTranscTime_ms = cell(1,NLoggers); % onset and offset time of the vocalization in ms in transceiver time
ManCallMicSamp = cell(1,NLoggers); % onset and offset time of the vocalization in microphone (raw) samples from the beginning of each individual 10min file
ManCallMicFile = cell(1,1400); % name of the raw microphone files that contain the call
NFiles = length(IndVocStart_all); % number of files that where manually cut in the first place from the microphone recording and then individually anaylzed using whocalls
voc_counter = 0;
for ff=1:NFiles
    fprintf(1,'File %d/%d\n', ff, NFiles)
    [~,File]=fileparts(VocFilename{ff});
    Ind_ = strfind(File, '_');
    FileRaw = File(1:(Ind_(end)-1));
    % onset and offset of manually detected vocalization events for all
    % loggers
    IndVocStart_local = IndVocStart_all{ff};
    IndVocStop_local = IndVocStop_all{ff};
    
    % loop through loggers and append the 2 columns matrix of onset/offset
    % of vocalization events
    for ll=1:NLoggers
        fprintf(1,'%s\n',AL_ManId{ll})
        if ~isempty(IndVocStart_local{ll})
            NewCalls = [IndVocStart_local{ll}' IndVocStop_local{ll}'] ./ Fs_env .* 10^3  + Voc_transc_time_refined(ff,1); % onset and offset time of the vocalization in ms in transceiver time
            ManCallTranscTime_ms{ll} = [ManCallTranscTime_ms{ll} ; NewCalls];
            
            
            NewCallsSample = [IndVocStart_local{ll}' IndVocStop_local{ll}'] ./ Fs_env .* FS  + Voc_samp_idx(ff,1); % onset and offset time of the vocalization in microphone (raw) samples from the beginning of each individual 10min file
            ManCallMicSamp{ll} = [ManCallMicSamp{ll} ; NewCallsSample];
            for voc=1:length(IndVocStart_local{ll})
                voc_counter = voc_counter+1;
                ManCallMicFile{voc_counter} = FileRaw;
            end
        end
    end
end
% ManCallSamples = cell2mat(ManCallSamp');
% ManCallFiles = ManCallFile(1:voc_counter);
% save('GroundTruthData.mat','ManCallFiles','ManCallSamples')
% %% Massage the input of the automatic detection
% 
% % Find the boundaries in ms in transceiver time of the section of sound that was
% % manually analysed
% Delay2FirstFile = Voc_samp_idx(1,1)./FS.*10^3; %in ms (Voc_samp_idx(1,1) is the first sample number of the first isolated section of sound in the first 10 min file analyzed here
% OnsetBoundary = Voc_transc_time(1,1) - Delay2FirstFile; % (Voc_Transc_time(1,1), is the absolute time onset in transceiver time of the first isolated section of sound in the first 10 min file analyzed here
% OffsetBoundary = OnsetBoundary + ManDur;
% 
% % Select the events in the automatically detected results that are within
% % the manually analysed window
% for ll=1:NLoggers
%     allCallTimes{ll} = cell2mat(allCallTimes{ll}').*10^(-3); % converting the cell array to a 2 column matrix and values from us to ms
%     Idx_WBound = find(sum((allCallTimes{ll}>OnsetBoundary) .* (allCallTimes{ll}<OffsetBoundary), 2)>0); % any element that start or end within the boundaries is selected
%     allCallTimes{ll} = allCallTimes{ll}(Idx_WBound,:);
% end
% 
% % %% Plot (interesting but too few vocalizations to see nicely things)
% % % order of loggers in auto detected events
% % Ord = [4 1:3];
% %
% % for ll=1:NLoggers
% %     figure(ll)
% %     for vv = 1:size(allCallTimes{Ord(ll)},1)
% %         plot(allCallTimes{Ord(ll)}(vv,:), ones(2,1)*2, 'b-', 'LineWidth',2)
% %         hold on
% %     end
% %     for vv = 1:size(ManCallTimes{ll},1)
% %         plot(ManCallTimes{ll}(vv,:), ones(2,1), 'g-', 'LineWidth',1)
% %         hold on
% %     end
% %     hold off
% % end

%% Let's loop in the dataset of manually extracted calls and see how many correct hits we get in the automatic detection
Delay = 100; %in ms error/delay between auto and man detection
AL_AutoId = fieldnames(SoundEvent_TranscTime_ms); % Names of the audioLoggers
MissedAutoDetection = [];
CorrectAutoDetection01 = zeros(size(SoundEvent_TranscTime_ms,1),1);
for ll=1:NLoggers
    for vv=1:size(ManCallTranscTime_ms{ll},1)
        fprintf(1, '*** %s vocalization %d/%d ****\n', AL_ManId{ll}, vv,size(ManCallTranscTime_ms{ll},1))
        OnOffVoc = ManCallTranscTime_ms{ll}(vv,:);
        ll_auto = contains(AL_AutoId, AL_ManId{ll});
        Idx_OnsetAuto = find(sum((SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,1)>(OnOffVoc(1)-Delay)) .* (SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,1)<(OnOffVoc(2)+Delay)), 2)>0);
        Idx_OffsetAuto = find(sum((SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,2)>(OnOffVoc(1)-Delay)) .* (SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,2)<(OnOffVoc(2)+Delay)), 2)>0);
        
        if isempty(Idx_OnsetAuto) && isempty(Idx_OffsetAuto)
            fprintf('No automatic call detected\n')
            disp(OnOffVoc)
            MissedAutoDetection = [MissedAutoDetection vv];
        else
            IdxAll = union(Idx_OnsetAuto,Idx_OffsetAuto);
            CorrectAutoDetection01(IdxAll) = ones(size(IdxAll));
            fprintf('%d events started during that vocalization\n',length(IdxAll))
            figure(1)
            clf
            plot(ManCallTranscTime_ms{ll}(vv,:), ones(2,1) * 1.1, 'g-', 'LineWidth',1)
            hold on
            for ii=1:length(IdxAll)
                plot(SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(IdxAll(ii),:), ones(2,1), 'b-', 'LineWidth',2)
                if ii==1
                    legend({'Manual', 'Automatic'})
                    legend('AutoUpdate','Off')
                end
                hold on
            end
            
            hold off
            xlabel('Time in transceiver Time')
            
            title(sprintf('%s  Voc %d/%d', AL_ManId{ll}, vv, size(ManCallTranscTime_ms{ll},1)))
            ylim([0,2])
            pause(1)
        end
    end
end




%Call I'm looking at: 4.053667319626574e10   4.053670019626574e10 of logger 10
%closest file index is 23, corresponds to 184549377 index and
%4.054173326222046e10 start
%difference of starts is -5.060065954715014e+03 ms. (off by 40 ms in loop??)
%then the index difference should be -2.530032977357507e+05
%and the call should be at 1.842963737022642e+08 which will be
%3.685927474045284e+06 on the graph
