%% Load input data
% Path2Data = '/Volumes/server_home/users/JulieE/LMC_CoEd/logger/20190623';
% Path2Audio = '/Volumes/server_home/users/JulieE/LMC_CoEd/audio/20190623';
% % Load data manually extracted
% load(fullfile(Path2Data, '190623_1401_VocExtractData_200.mat'))
% load(fullfile(Path2Data, '190623_1401_VocExtractData.mat'))
% load(fullfile(Path2Audio, '190623_1401_VocExtractTimes.mat'))
addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'))
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'))

%% Load data manually extracted
Path2Data1 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audiologgers';
Path2Results1 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audiologgers/GroundTruthResultsPipelineCheck';
% Path2Data1 = '/Volumes/server_home/users/JulieE/JuvenileRecordings155/20190927/audiologgers';
load(fullfile(Path2Results1,'190927_1014_VocExtractData_200.mat'))
load(fullfile(Path2Results1,'190927_1014_VocExtractData.mat'))
% Path2Data2 = '/Volumes/server_home/users/JulieE/JuvenileRecordings155/20190927/audio';
Path2Data2 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audio';
Path2Results2 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audio/GroundTruthResultsPipelineCheck';
load(fullfile(Path2Results2,'190927_1014_VocExtractTimes.mat'))

Fs_env=1000; % in Hertz, should have been saved in who_calls.m to correctly convert in time the starting and ending indices of vocalizations in IndVocStart
FS_Piezo = 50000; % could also be retrieved from Piezo_FS

%% Load the data saved of that script
%  load(fullfile(Path2Results1, 'SoundEvent.mat'))
%% Load the data of CoEd
Path2Data1 = '/Volumes/Julie4T/LMC_CoEd/logger/20190603';
Path2Results1 = '/Volumes/Julie4T/LMC_CoEd/logger/20190603/GroundTruthResultsRecOnly';
load(fullfile(Path2Results1,'190603_1447_VocExtractData_200.mat'))
load(fullfile(Path2Results1,'190603_1447_VocExtractData.mat'))

Path2Data2 = '/Volumes/Julie4T/LMC_CoEd/audio/20190603';
Path2Results2 = '/Volumes/Julie4T/LMC_CoEd/audio/20190603/GroundTruthResults';
load(fullfile(Path2Results2,'190603_1447_VocExtractTimes.mat'))
Fs_env=1000; % in Hertz, should have been saved in who_calls.m to correctly convert in time the starting and ending indices of vocalizations in IndVocStart
FS_Piezo = 50000; % could also be retrieved from Piezo_FS

%% Get in transceiver time the onset/offset of each vocalization manually extracted for each logger
AllLoggers = dir(fullfile(Path2Data1, '*ogger*'));
DirFlags = [AllLoggers.isdir];
% Extract only those that are directories.
AllLoggers = AllLoggers(DirFlags);
AL_ManId = fieldnames(Piezo_wave); % Names of the audioLoggers
NLoggers = length(AudioLogs);% number of audio loggers
ManCallTranscTime_ms = cell(1,NLoggers); % onset and offset time of the vocalization in ms in transceiver time
ManCallMicSamp = cell(1,NLoggers); % onset and offset time of the vocalization in microphone (raw) samples from the beginning of each individual 10min file
ManCallLogSamp = cell(1,NLoggers); % onset and offset time of the vocalization in logger samples from the beginning of each individual 10min file
SamplingFreq = cell(1, NLoggers);
ManCallMicFile = cell(1,NLoggers); % name of the raw microphone files that contain the call
NFiles = length(IndVocStart_all); % number of files that where manually cut in the first place from the microphone recording and then individually anaylzed using whocalls
voc_counter = zeros(1,NLoggers);
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
        Data_directory = fullfile(AllLoggers(ll).folder,AL_ManId{ll}, 'extracted_data');
        File = dir(fullfile(Data_directory, '*CSC0*'));
        if isempty(File)
            error('Data file not found');
        end
        Filepath = fullfile(File.folder, File.name);
        load(Filepath, 'Timestamps_of_first_samples_usec','Estimated_channelFS_Transceiver','Indices_of_first_and_last_samples')
        
        if ~isempty(IndVocStart_local{ll})
            NewCalls = [IndVocStart_local{ll}' IndVocStop_local{ll}'] ./ Fs_env .* 10^3  + Voc_transc_time_refined(ff,1); % onset and offset time of the vocalization in ms in transceiver time
            ManCallTranscTime_ms{ll} = [ManCallTranscTime_ms{ll} ; NewCalls];
            
            
            NewCallsSample = [IndVocStart_local{ll}' IndVocStop_local{ll}'] ./ Fs_env .* FS  + Voc_samp_idx(ff,1); % onset and offset time of the vocalization in microphone (raw) samples from the beginning of each individual 10min file
            ManCallMicSamp{ll} = [ManCallMicSamp{ll} ; NewCallsSample];
            for voc=1:length(IndVocStart_local{ll})
                voc_counter(ll) = voc_counter(ll)+1;
                ManCallMicFile{ll}{voc_counter(ll)} = FileRaw;
            end
            
            
            % Find the onset and offset on the logger recording
            for voc=1:length(IndVocStart_local{ll})
                VocOnset_time = NewCalls(voc,1);
                VocOffset_time = NewCalls(voc,2);
                % find the time stamp on the logger that is closest to before
                % the snippet of sound onset
                IndTSOn = find(Timestamps_of_first_samples_usec<(VocOnset_time*10^3), 1, 'Last');
                
                % find the time stamp on the logger that is closest to after
                % the snippet of sound offset
                IndTSOff = find(Timestamps_of_first_samples_usec>(VocOffset_time*10^3), 1, 'First');
                if ~isempty(IndTSOff)
                    % deduct the corresponding onset and offset samples
                    if IndTSOn<=length(Estimated_channelFS_Transceiver) && ~isnan(Estimated_channelFS_Transceiver(IndTSOn))
                        IndSampOn = round(Indices_of_first_and_last_samples(IndTSOn,1) + Estimated_channelFS_Transceiver(IndTSOn)*(10^-6)*(VocOnset_time*10^3 - Timestamps_of_first_samples_usec(IndTSOn)));
                    else
                        IndSampOn = round(Indices_of_first_and_last_samples(IndTSOn,1) + nanmean(Estimated_channelFS_Transceiver)*(10^-6)*(VocOnset_time*10^3 - Timestamps_of_first_samples_usec(IndTSOn)));
                    end
                    if IndTSOff<=length(Estimated_channelFS_Transceiver) && ~isnan(Estimated_channelFS_Transceiver(IndTSOff))
                        IndSampOff = round(Indices_of_first_and_last_samples(IndTSOff,1) - Estimated_channelFS_Transceiver(IndTSOff)*(10^-6)*(Timestamps_of_first_samples_usec(IndTSOff) - VocOffset_time*10^3));
                    else
                        IndSampOff = round(Indices_of_first_and_last_samples(IndTSOff,1) - nanmean(Estimated_channelFS_Transceiver)*(10^-6)*(Timestamps_of_first_samples_usec(IndTSOff) - VocOffset_time*10^3));
                    end
                    if IndSampOff<=IndSampOn
                        IndSampOff = round(Indices_of_first_and_last_samples(IndTSOn,1) + nanmean(Estimated_channelFS_Transceiver)*(10^-6)*(VocOffset_time*10^3 - Timestamps_of_first_samples_usec(IndTSOn)));
                    end
                else
                    % find the time stamp on the logger that is closest to before
                    % the snippet of sound offset
                    IndTSOff = find(Timestamps_of_first_samples_usec<(VocOffset_time*10^3), 1, 'Last');
                    % this vocalization is in the last recording file
                    % There is no estimation of the sample frequency for that last
                    % file. Let's estimate it as the average of the previous
                    % estimates
                    FS_local = nanmean(Estimated_channelFS_Transceiver);
                    IndSampOn = round(Indices_of_first_and_last_samples(IndTSOn,1) + FS_local*(10^-6)*(VocOnset_time*10^3 - Timestamps_of_first_samples_usec(IndTSOn)));
                    IndSampOff = round(Indices_of_first_and_last_samples(IndTSOff,1) + FS_local*(10^-6)*(VocOffset_time*10^3 - Timestamps_of_first_samples_usec(IndTSOff)));
                end
                ManCallLogSamp{ll} = [ManCallLogSamp{ll} ; IndSampOn IndSampOff];
                SamplingFreq{ll} = nanmean(Estimated_channelFS_Transceiver);
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
save(fullfile(Path2Results1, 'SoundEvent.mat'),'AL_ManId','NLoggers','ManCallTranscTime_ms','ManCallMicSamp','ManCallLogSamp', 'SamplingFreq','ManCallMicFile')








%% Run the automatic detection based on loggers with various threshold on RMS
FigOn = 0;
AllLoggers = dir(fullfile(Path2Data1, '*ogger*'));
DirFlags = [AllLoggers.isdir];
% Extract only those that are directories.
AllLoggers = AllLoggers(DirFlags);
NL = length(AllLoggers);

AllFiles = dir(fullfile(Path2Data2,'*RecOnly*.wav'));
% find the date and expstart time
Date = AllFiles(1).name(6:11);
ExpStartTime = AllFiles(1).name(13:16);
TTL_dir = dir(fullfile(Path2Data2,sprintf( '%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
TTL = load(fullfile(TTL_dir.folder, TTL_dir.name));
if strcmp(Date, '190927')
    MaxTranscTime = max(TTL.Pulse_TimeStamp_Transc(TTL.File_number==10));
elseif strcmp(Date, '190603')
    MaxTranscTime = max(TTL.Pulse_TimeStamp_Transc);
end

%RMSThresh = [1.5 2 3];
RMSThresh = 1.5;
MissedAutoDetection = cell(length(RMSThresh),1);
CorrectAutoDetection01 = cell(length(RMSThresh),1);
TotMissedAutoCall = nan(length(RMSThresh),1);
for tt=1:length(RMSThresh)
    fprintf(1,'*** RMS Thresh %d/%d *****', tt, length(RMSThresh))
    for ll=1:NL
        Data_directory = fullfile(AllLoggers(ll).folder,AllLoggers(ll).name, 'extracted_data');
        [SE_LS,SE_TT_ms,LoggerEnvelopeAll.(sprintf('L%s',AllLoggers(ll).name(2:end))),~,~] = piezo_find_calls_logger(Data_directory, RMSThresh(tt));
        % only keep detection that are within the period of analysis
        SoundEvent_TranscTime_ms.(sprintf('L%s',AllLoggers(ll).name(2:end))) = SE_TT_ms((SE_TT_ms(:,1)<MaxTranscTime),:);
        SoundEvent_LoggerSamp.(sprintf('L%s',AllLoggers(ll).name(2:end))) = SE_LS((SE_TT_ms(:,1)<MaxTranscTime),:);
    end
    save(fullfile(Path2Results1, 'SoundEvent.mat'),'SoundEvent_LoggerSamp','SoundEvent_TranscTime_ms','LoggerEnvelopeAll','-append')
    % load(fullfile(Path2Data1, 'SoundEvent.mat')
    
    
    
    
    
    %% Let's loop in the dataset of manually extracted calls and see how many correct hits we get in the automatic detection
    FS_env = 1000; %Sampling Frequency of the envelope as calculated by piezo_find_calls_logger
    Delay = 50; %in ms error/delay between auto and man detection and Delay to add before each detected call in ms
    AL_AutoId = fieldnames(SoundEvent_TranscTime_ms); % Names of the audioLoggers
    MissedAutoDetection{tt} = cell(NLoggers,1);
    TotManCall = 0;
    %CorrectAutoDetection01{tt} = struct();
    CorrectAutoDetection01{tt} = cell(NLoggers,1);
    
    for ll=1:NLoggers
        MissedAutoDetection{tt}{ll} = [];
        ll_auto = contains(AL_AutoId, AL_ManId{ll});
        %% Load the raw signal
        Data_directory = fullfile(AllLoggers(ll).folder,AL_ManId{ll}, 'extracted_data');
        File = dir(fullfile(Data_directory, '*CSC0*'));
        if isempty(File)
            error('Data file not found');
        end
        Filepath = fullfile(File.folder, File.name);
        load(Filepath, 'AD_count_int16')
        AD_count_double = double(AD_count_int16);
        clear AD_count_int16
        % Center the signal and clear the old data from memory
        Centered_piezo_signal = AD_count_double - mean(AD_count_double);
        clear AD_count_double
        %         CorrectAutoDetection01{tt}.(sprintf(AL_AutoId{ll_auto})) = zeros(size(SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto})),1),1);
        CorrectAutoDetection01{tt}{ll_auto} = zeros(size(SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto})),1),1);
        
        TotManCall = TotManCall+size(ManCallTranscTime_ms{ll},1);
        
        for vv=1:size(ManCallTranscTime_ms{ll},1)
            fprintf(1, '*** %s vocalization %d/%d ****\n', AL_ManId{ll}, vv,size(ManCallTranscTime_ms{ll},1))
            OnOffVoc = ManCallTranscTime_ms{ll}(vv,:);
            Idx_OnsetAuto = find((SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,1)>(OnOffVoc(1)-Delay)) .* (SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,1)<(OnOffVoc(2)+Delay))>0);
            Idx_OffsetAuto = find((SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,2)>(OnOffVoc(1)-Delay)) .* (SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,2)<(OnOffVoc(2)+Delay))>0);
            
            if isempty(Idx_OnsetAuto) && isempty(Idx_OffsetAuto)
                Idx_OnsetAuto = find((SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,1)<(OnOffVoc(1)-Delay)) .* (SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,2)>(OnOffVoc(1)-Delay))>0);
                Idx_OffsetAuto = find((SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,1)<(OnOffVoc(2)+Delay)) .* (SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto}))(:,2)>(OnOffVoc(2)+Delay))>0);
            end
            if isempty(Idx_OnsetAuto) && isempty(Idx_OffsetAuto)
                fprintf('No automatic call detected\n')
                disp(OnOffVoc)
                MissedAutoDetection{tt}{ll} = [MissedAutoDetection{tt}{ll} vv];
                %                 if FigOn
                figure(1)
                clf
                
                % Plot the spectrogram of the sound extract
                DBNoise = 60; % amplitude parameter for the color scale of the spectro
                FHigh = 10000; % y axis max scale for the spectrogram
                x_start = round(ManCallLogSamp{ll}(vv,1)-Delay*10^(-3)*SamplingFreq{ll});
                x_stop = round(ManCallLogSamp{ll}(vv,2)+Delay*10^(-3)*SamplingFreq{ll});
                Raw = Centered_piezo_signal(x_start:x_stop);
                Raw_ramp = cosramp(Raw-mean(Raw), SamplingFreq{ll}*10*10^-3);
                [~] = spec_only_bats(Raw_ramp,SamplingFreq{ll},DBNoise, FHigh);
                caxis('manual');
                caxis([2 70]);
                ylim([-500 10000])
                hold on
                
                yyaxis right %% There is always a problem in the plot for allignment of the envelope that I attribute to the average sample frequency estimate... Teh onset/offset detection is correctly alligned though!
                x_start_env = round((ManCallLogSamp{ll}(vv,1)/SamplingFreq{ll}-Delay*10^(-3))*FS_env);
                %             x_start_env = round(x_start/SamplingFreq{ll}*FS_env);
                x_stop_env = round((ManCallLogSamp{ll}(vv,2)/SamplingFreq{ll}+Delay*10^(-3))*FS_env);
                %             x_stop_env = round(x_stop/SamplingFreq{ll}*FS_env);
                plot(LoggerEnvelopeAll.(sprintf(AL_AutoId{ll_auto}))(x_start_env:x_stop_env), '-k','LineWidth',2)
                hold on
                %             hline((RMSfactor * Noise))
                %             hold on
                ylim([-10 300])
                ylabel('Amplitude Envelope')
                
                AllignOn = ManCallLogSamp{ll}(vv,1)-Delay*10^(-3)*SamplingFreq{ll}; % in logger samples
                x_start_man = round(Delay*10^(-3)*FS_env);
                x_stop_man = round((ManCallLogSamp{ll}(vv,2)-AllignOn)/SamplingFreq{ll}*FS_env);
                plot([x_start_man x_stop_man], ones(2,1) * 1.5, 'g-', 'LineWidth',1)
                hold off
                
                
                title(sprintf('%s  Voc %d/%d NOT DETECTED', AL_ManId{ll}, vv, size(ManCallTranscTime_ms{ll},1)))
                Player = audioplayer(Raw_ramp/std(Raw_ramp), SamplingFreq{ll});
                play(Player)
                pause()
                %                 end
                
            else
                IdxAll = union(Idx_OnsetAuto,Idx_OffsetAuto);
                %                 CorrectAutoDetection01{tt}.(sprintf(AL_AutoId{ll_auto}))(IdxAll) = ones(size(IdxAll));
                CorrectAutoDetection01{tt}{ll_auto}(IdxAll) = ones(size(IdxAll));
                fprintf('%d events started during that vocalization\n',length(IdxAll))
                if FigOn
                    figure(1)
                    clf
                    
                    % Plot the spectrogram of the sound extract
                    DBNoise = 60; % amplitude parameter for the color scale of the spectro
                    FHigh = 10000; % y axis max scale for the spectrogram
                    x_start = round(ManCallLogSamp{ll}(vv,1)-Delay*10^(-3)*SamplingFreq{ll});
                    x_stop = round(ManCallLogSamp{ll}(vv,2)+Delay*10^(-3)*SamplingFreq{ll});
                    Raw = Centered_piezo_signal(x_start:x_stop);
                    Raw_ramp = cosramp(Raw-mean(Raw), SamplingFreq{ll}*10*10^-3);
                    [~] = spec_only_bats(Raw_ramp,SamplingFreq{ll},DBNoise, FHigh);
                    caxis('manual');
                    caxis([2 70]);
                    ylim([-500 10000])
                    hold on
                    
                    yyaxis right
                    x_start_env = round((ManCallLogSamp{ll}(vv,1)/SamplingFreq{ll}-Delay*10^(-3))*FS_env);
                    %             x_start_env = round(x_start/SamplingFreq{ll}*FS_env);
                    x_stop_env = round((ManCallLogSamp{ll}(vv,2)/SamplingFreq{ll}+Delay*10^(-3))*FS_env);
                    %             x_stop_env = round(x_stop/SamplingFreq{ll}*FS_env);
                    plot(LoggerEnvelopeAll.(sprintf(AL_AutoId{ll_auto}))(x_start_env:x_stop_env), '-k','LineWidth',2)
                    hold on
                    %             hline((RMSfactor * Noise))
                    %             hold on
                    ylim([-10 300])
                    ylabel('Amplitude Envelope')
                    
                    AllignOn = ManCallLogSamp{ll}(vv,1)-Delay*10^(-3)*SamplingFreq{ll}; % in logger samples
                    x_start_man = round(Delay*10^(-3)*FS_env);
                    x_stop_man = round((ManCallLogSamp{ll}(vv,2)-AllignOn)/SamplingFreq{ll}*FS_env);
                    plot([x_start_man x_stop_man], ones(2,1) * 1.5, 'g-', 'LineWidth',1)
                    hold on
                    for ii=1:length(IdxAll)
                        OnOff = (SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(IdxAll(ii),:) - AllignOn)/SamplingFreq{ll}*FS_env;
                        plot(OnOff, ones(2,1) * (-1.5), 'b-', 'LineWidth',2)
                        if ii==1
                            legend({ 'Automatic' ,'Manual'})
                            legend('AutoUpdate','Off')
                        end
                        hold on
                    end
                    hold off
                    
                    
                    title(sprintf('%s  Voc %d/%d', AL_ManId{ll}, vv, size(ManCallTranscTime_ms{ll},1)))
                    %             Player = audioplayer(Raw_ramp/std(Raw_ramp), SamplingFreq{ll});
                    %             play(Player)
                    %             pause()
                end
            end
        end
    end
    
    
    TotMissedAutoCall(tt) = sum(cellfun('length',MissedAutoDetection{tt}));
    fprintf(1,'Missed call by the logger auto detection: %d/%d or %.1f%%\n', sum(cellfun('length',MissedAutoDetection{tt})),TotManCall,sum(cellfun('length',MissedAutoDetection{tt}))/TotManCall*100);
    fprintf(1,'Number of detected events from the logger %d, proportion of true calls %d/%d or %.1f%%\n',sum(cellfun('length',CorrectAutoDetection01{tt})),sum(cellfun(@sum,CorrectAutoDetection01{tt})),sum(cellfun('length',CorrectAutoDetection01{tt})),sum(cellfun(@sum,CorrectAutoDetection01{tt}))/sum(cellfun('length',CorrectAutoDetection01{tt}))*100);
end



%% Plot the results of various levels of RMS Thresh
TotCorrectAutoDetection01 = nan(length(RMSThresh),1);
TotAutoDetection01 = nan(length(RMSThresh),1);
for tt=1:length(RMSThresh)
    TotCorrectAutoDetection01(tt) = sum(cellfun(@sum,CorrectAutoDetection01{tt}));
    TotAutoDetection01(tt) = sum(cellfun(@length,CorrectAutoDetection01{tt}));
end
figure()
subplot(2,1,1)
plot(TotMissedAutoCall, 'LineWidth',2)
%  hold on
%  plot(TotCorrectAutoDetection01, 'LineWidth',2)
hold on
ylabel(sprintf('Number of Missed calls out of %d',TotManCall))
set(gca,'XTick', 1:length(RMSThresh))
set(gca,'XTickLabel', RMSThresh)
xlabel('RMS threshold factor on logger envelope')
hold off
ylim([0 25])
title('Logger detection')
subplot(2,1,2)
plot(TotMissedAutoCall/TotManCall*100, 'LineWidth',2)
hold on
plot(TotCorrectAutoDetection01./TotAutoDetection01.*100, 'LineWidth',2)
ylabel('Percentage')
ylim([0 10])
yyaxis right
hold on
plot(100-TotCorrectAutoDetection01./TotAutoDetection01.*100, 'LineWidth',2)
legend('Missed calls', 'Correct detection','False detection')
set(gca,'XTick', 1:length(RMSThresh))
set(gca,'XTickLabel', RMSThresh)
ylabel('Percentage')
xlabel('RMS threshold factor on logger envelope')
ylim([80 100])
hold off
%%
save(fullfile(Path2Results1, 'SoundEvent.mat'),'SoundEvent_LoggerSamp','SoundEvent_TranscTime_ms','LoggerEnvelopeAll', 'CorrectAutoDetection01','MissedAutoDetection','RMSThresh','Delay','-append')




%% Now Loop through the result of the logger detection and check what we have on the microphone
AL_AutoId = fieldnames(SoundEvent_TranscTime_ms); % Names of the audioLoggers
NL = length(AL_AutoId);
Buffer = 20; % in ms
Flow = 500;
FHigh = 10000;
[z,p,k] = butter(6,[Flow FHigh]/(192000/2),'bandpass');
sos_mic_band = zp2sos(z,p,k);
[z,p,k] = butter(6,[Flow FHigh]/(50000/2),'bandpass');
sos_logger_band = zp2sos(z,p,k);
RMS = cell(1,NL);
MaxAmp = cell(1,NL);
MeanAmp = cell(1,NL);
CC_wave = cell(1,NL);
CC_spectro = cell(1,NL);

for ll=1:NL
    fprintf(1,'Loading data of %s (%d/%d)\n',AL_AutoId{ll},ll,NL)
    % Load the raw logger signal
    Data_directory = fullfile(AllLoggers(ll).folder,AL_AutoId{ll}, 'extracted_data');
    File = dir(fullfile(Data_directory, '*CSC0*'));
    if isempty(File)
        error('Data file not found');
    end
    Filepath = fullfile(File.folder, File.name);
    load(Filepath, 'AD_count_int16')
    AD_count_double = double(AD_count_int16);
    clear AD_count_int16
    % Center the signal and clear the old data from memory
    Centered_piezo_signal = AD_count_double - mean(AD_count_double);
    clear AD_count_double
    
    % convert transceiver time to audio samp files
    Nevents = size(SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll})),1);
    [MicVoc_samp_idx,MicVoc_File]=transc_time2micsamp(Path2Data2,SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll})));
    
    % Loop through events
    RMS{ll} = nan(1,Nevents);
    MaxAmp{ll} = nan(1,Nevents);
    MeanAmp{ll} = nan(1,Nevents);
    CC_wave{ll} = nan(1,Nevents);
    CC_spectro{ll} = nan(1,Nevents);
    OldMicVoc_File = 0;
    for ee=1:Nevents
        fprintf(1,'%s event %d/%d\n', AL_AutoId{ll}, ee, Nevents)
        if OldMicVoc_File~=MicVoc_File(ee)
            RawWavDir = dir(fullfile(Path2Data2,sprintf('*mic1_%d.wav',MicVoc_File(ee))));
            [RawWav_local, FS_raw] = audioread(fullfile(RawWavDir.folder, RawWavDir.name));
            OldMicVoc_File = MicVoc_File(ee);
        end
        mic_start = round(MicVoc_samp_idx(ee,1)-Buffer*FS_raw*10^-3);
        mic_stop = round(MicVoc_samp_idx(ee,2)+Buffer*FS_raw*10^-3);
        if mic_start<0 || mic_stop<0% call occur before microphone started recording
            continue
        end
        
        if mic_stop>length(RawWav_local) % this section is cut between 2 10 min recordings
            mic_stop1 = length(RawWav_local);
            mic_stop2 = mic_stop - length(RawWav_local);
            RawWavDir2 = dir(fullfile(Path2Data2,sprintf('*mic1_%d.wav',MicVoc_File(ee)+1)));
            if ~isempty(RawWavDir2)
                [RawWav_local2, FS_raw] = audioread(fullfile(RawWavDir2.folder, RawWavDir2.name));
                Mic_Data = [RawWav_local(mic_start:mic_stop1); RawWav_local2(1:mic_stop2)];
                OldMicVoc_File = MicVoc_File(ee)+1;
                RawWav_local = RawWav_local2;
                clear RawWav_local2
            else % that was the last recording from the microphone
                if mic_start>length(RawWav_local) % call occur after microphone stopped recording
                    continue
                end
                Mic_Data = RawWav_local(mic_start:mic_stop1);
            end
        else
            Mic_Data = RawWav_local(mic_start:mic_stop);
        end
        FS_log = round(SamplingFreq{contains(AL_ManId,AL_AutoId{ll})});
        x_start = round(SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll}))(ee,1) - Buffer*FS_log*10^-3);
        x_stop = round(SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll}))(ee,2)+ Buffer*FS_log*10^-3);
        Logger_Data = Centered_piezo_signal(x_start:x_stop);
        % filter the data
        Mic_Data = Mic_Data-mean(Mic_Data);
        Logger_Data = Logger_Data-mean(Logger_Data);
        Filt_MicData = filtfilt(sos_mic_band,1,Mic_Data);
        Filt_LoggerData = filtfilt(sos_logger_band,1,Logger_Data);
        % Calculate the RMS in the microphone extract and check if it's
        % above thershold
        AmpEnv=envelope(Filt_MicData,FS_raw/FS_env,'rms');
        RMS{ll}(ee)=std(Filt_MicData(Filt_MicData~=0));
        MaxAmp{ll}(ee) = max(AmpEnv);
        MeanAmp{ll}(ee) = mean(AmpEnv);
        ResampFilt_MicData = resample(Filt_MicData,4*FHigh,FS_raw);
        ResampFilt_LogData = resample(Filt_LoggerData,4*FHigh,FS_log);
        if length(ResampFilt_MicData)~=length(ResampFilt_LogData)
            OptLength = min(length(ResampFilt_MicData),length(ResampFilt_LogData));
            ResampFilt_MicData = ResampFilt_MicData(1:OptLength);
            ResampFilt_LogData = ResampFilt_LogData(1:OptLength);
        end
        figure(10)
        clf
        title(sprintf('%s event %d/%d', AL_AutoId{ll}, ee, Nevents))
        subplot(2,1,1)
        [toM, foM, logBM,~] = spec_only_bats(ResampFilt_MicData,4*FHigh,DBNoise, FHigh,100,'Time_increment',0.005);
        subplot(2,1,2)
        [toL, foL, logBL,~] = spec_only_bats(ResampFilt_LogData,4*FHigh,DBNoise, FHigh,100,'Time_increment',0.005);
        
        
        % Do a cross correlation between the two signals and see how high
        % it is
        
        [Xcor,Lag] = xcorr(ResampFilt_MicData,ResampFilt_LogData, (2*Buffer*10^-3)*4*FHigh,'normalized'); % Running a cross correlation between the raw signal and each audio logger signal with a maximum lag equal to twice the Buffer size
        CC_wave{ll}(ee) = max(abs(Xcor));
        if max(toM)<=0.5
            Xcor2 = xcorr2(logBM, logBL); % Running a cross correlation between the raw signal and each audio logger signal
            CC_spectro{ll}(ee) = max(max(Xcor2));
        end
        
    end
end
%% Look at the features of amplitude in the microphone and correlation for correct vs incorrect detected elements in the loggers
% RMS
CorrectAutoDetection01_All = cell(size(CorrectAutoDetection01{1}'));
for ll=1:NL
    CorrectAutoDetection01_All{ll} = CorrectAutoDetection01{1}{ll}';
end
RMS = [RMS{:}];
MaxAmp = [MaxAmp{:}];
MeanAmp = [MeanAmp{:}];
CC_wave = [CC_wave{:}];
CC_spectro = [CC_spectro{:}];
CorrectAutoDetection01_All = [CorrectAutoDetection01_All{:}];

% Trim data to only look at non nans (these for which we have ground truth)
RMS = RMS(~isnan(RMS));
MaxAmp = MaxAmp(~isnan(RMS));
MeanAmp = MeanAmp(~isnan(RMS));
CC_wave = CC_wave(~isnan(RMS));
CC_spectro = CC_spectro(~isnan(RMS));
CorrectAutoDetection01_All = CorrectAutoDetection01_All(~isnan(RMS));

figure();
histogram(RMS(~CorrectAutoDetection01_All), 'Normalization', 'probability', 'BinWidth',0.0001)
hold on
histogram(RMS(logical(CorrectAutoDetection01_All)), 'Normalization', 'probability', 'BinWidth',0.0001)
hold off
legend('Noise','Vocalizations')
title('Mic RMS')

% MaxAmp could be a good criteria for sorting noise from calls with
% threshold set at MaxAmp>0.0001 (only loosing 2 calls)
figure();
histogram(MaxAmp(~CorrectAutoDetection01_All), 'Normalization', 'probability', 'BinWidth',0.001)
hold on
histogram(MaxAmp(logical(CorrectAutoDetection01_All)), 'Normalization', 'probability', 'BinWidth',0.001)
hold off
legend('Noise','Vocalizations')
title('Mic MaxAmp')

figure();
histogram(MeanAmp(~CorrectAutoDetection01_All), 'Normalization', 'probability', 'BinWidth',0.001)
hold on
histogram(MeanAmp(logical(CorrectAutoDetection01_All)), 'Normalization', 'probability', 'BinWidth',0.001)
hold off
legend('Noise','Vocalizations')
title('Mic MeanAmp')

% CC_wave could be a good measure with threshold set at
figure();
histogram(CC_wave(~CorrectAutoDetection01_All), 'Normalization', 'probability', 'BinWidth',0.001)
hold on
histogram(CC_wave(logical(CorrectAutoDetection01_All)), 'Normalization', 'probability', 'BinWidth',0.001)
hold off
legend('Noise','Vocalizations')
title('Mic CC_wave')

figure();
histogram(CC_spectro(~CorrectAutoDetection01_All), 'Normalization', 'probability')
hold on
histogram(CC_spectro(logical(CorrectAutoDetection01_All)), 'Normalization', 'probability')
hold off
legend('Noise','Vocalizations')
title('Mic CC_spectro')


%% Now loop through the detected elements and calculate biosound
WorkingDirSpectro = '/Users/elie/Documents/GroundTruthWorkDir';
mkdir(WorkingDirSpectro)
Buffer = 30;% time in ms to add before after each sound element such atht it's longer than the 23ms required for biosound to calculate fundamental and saliency parameters
F_High = 5000;
F_low = 100;
F_highSpec = 15000;
Flow = 500;
FHigh = 10000;
DBNoise = 60;
SpectroYN = 0;
[z,p,k] = butter(6,[Flow FHigh]/(192000/2),'bandpass');
sos_mic_band = zp2sos(z,p,k);
[z,p,k] = butter(6,[Flow FHigh]/(50000/2),'bandpass');
FS_env = 1000; %Sampling Frequency of the envelope as calculated by piezo_find_calls_logger
AllLoggers = dir(fullfile(Path2Data1, '*ogger*'));
DirFlags = [AllLoggers.isdir];
% Extract only those that are directories.
AllLoggers = AllLoggers(DirFlags);
NL = length(AllLoggers);

% BioSoundUniqParam = nan(21553,23);
BioSoundUniqParam = cell(1,NLoggers);
% BioSoundUniqParam = nan(24770,21);
ee_count = 0;
% BioSoundParamNames = {'stdtime' 'meantime' 'skewtime' 'entropytime'...
%         'kurtosistime' 'AmpPeriodF' 'AmpPeriodP' 'rms' 'maxAmp' 'stdspect'...
%         'meanspect' 'skewspect' 'entropyspect' 'kurtosisspect' 'q1' 'q2' 'q3'...
%         'fund' 'cvfund' 'minfund' 'maxfund' 'meansal' '01correct'};


% Turn off warnings regarding Pyton to structure conversion
% warning('off', 'MATLAB:structOnObject')

AL_AutoId = fieldnames(SoundEvent_TranscTime_ms); % Names of the audioLoggers
AL_ManId = fieldnames(Piezo_wave); % Names of the audioLoggers
% BiosoundFolder = cell(NLoggers,1);
% Data_out = fullfile(AllLoggers(ll).folder, 'BiosoundEvents');
% if ~exist(Data_out, 'dir')
%     mkdir(Data_out)
% end
parfor ll=1:NLoggers
    ll_auto = contains(AL_AutoId, AL_ManId{ll});
    fprintf(1, '*** %s %d/%d ****\n',AL_ManId{ll}, ll, NLoggers)
    % Load the raw signal
    Data_directory = fullfile(AllLoggers(ll).folder,AL_ManId{ll}, 'extracted_data');
    File = dir(fullfile(Data_directory, '*CSC0*'));
    if isempty(File)
        error('Data file not found');
    end
    Filepath = fullfile(File.folder, File.name);
    DataL=load(Filepath, 'AD_count_int16', 'Indices_of_first_and_last_samples','Estimated_channelFS_Transceiver')
    AD_count_double = double(DataL.AD_count_int16);
    DataL.AD_count_int16 = [];
    % Center the signal and clear the old data from memory
    Centered_piezo_signal = AD_count_double - mean(AD_count_double);
    AD_count_double=[];
    
    % convert transceiver time to audio samp files
    Nevents = size(SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto})),1);
    [MicVoc_samp_idx,MicVoc_File]=transc_time2micsamp(Path2Data2,SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto})));
    
    
    % Loop through sound events
    OldMicVoc_File = 0;
    TotEv = size(SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto})),1);
    BioSoundUniqParam{ll} = cell(1,TotEv);
    for ee=1:TotEv
        ee_count = ee_count+1;
        if rem(ee,100)==0
            fprintf(1, 'Event %d/%d\n', ee,TotEv)
        end
        
        % find the sampling Frequency
        OnInd1 = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,1);
        OffInd1 = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,2);
        FileIdx = find((DataL.Indices_of_first_and_last_samples(:,1)<OnInd1) .* (DataL.Indices_of_first_and_last_samples(:,2)>OffInd1));
        if isempty(FileIdx) || (length(FileIdx)~=1) || FileIdx>length(DataL.Estimated_channelFS_Transceiver)
            FS_local = round(nanmean(DataL.Estimated_channelFS_Transceiver));
        else
            FS_local = round(DataL.Estimated_channelFS_Transceiver(FileIdx));
        end
        
        % extract the sound with Buffer ms before after the sound
        OnInd_logger = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,1) - round(FS_local*Buffer*10^-3);
        OffInd_logger = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,2) + round(FS_local*Buffer*10^-3);
        if OnInd_logger<0
            OffInd_logger = OffInd_logger-OnInd_logger;
            OnInd_logger=1;
        end
        Logger_Data = Centered_piezo_signal(OnInd_logger : OffInd_logger);
        Logger_Data = Logger_Data - mean(Logger_Data);
        
        
        [BioSoundUP,~] = run_acoustic_features(Logger_Data, FS_local, 'F_High',F_High, 'F_low',F_low, 'F_highSpec',F_highSpec,'Spectro',SpectroYN);
        
        BioSoundUniqParam{ll}{ee} = [BioSoundUP'; nan(5,1)];
        
        % Add parameters regarding the microphone data
        if OldMicVoc_File~=MicVoc_File(ee)
            RawWavDir = dir(fullfile(Path2Data2,sprintf('*RecOnly_mic1_%d.wav',MicVoc_File(ee))));
            [RawWav_mic, FS_mic] = audioread(fullfile(RawWavDir.folder, RawWavDir.name));
            OldMicVoc_File = MicVoc_File(ee);
        end
        mic_start = round(MicVoc_samp_idx(ee,1)-Buffer*FS_mic*10^-3);
        mic_stop = round(MicVoc_samp_idx(ee,2)+Buffer*FS_mic*10^-3);
        if mic_start<0 || mic_stop<0% call occur before microphone started recording
            continue
        end
        
        if mic_stop>length(RawWav_mic) % this section is cut between 2 10 min recordings
            mic_stop1 = length(RawWav_mic);
            mic_stop2 = mic_stop - length(RawWav_mic);
            RawWavDir2 = dir(fullfile(Path2Data2,sprintf('*RecOnly_mic1_%d.wav',MicVoc_File(ee)+1)));
            if ~isempty(RawWavDir2)
                [RawWav_local2, FS_mic] = audioread(fullfile(RawWavDir2.folder, RawWavDir2.name));
                Mic_Data = [RawWav_mic(mic_start:mic_stop1); RawWav_local2(1:mic_stop2)];
                OldMicVoc_File = MicVoc_File(ee)+1;
                RawWav_mic = RawWav_local2;
                RawWav_local2 = [];
            else % that was the last recording from the microphone
                if mic_start>length(RawWav_mic) % call occur after microphone stopped recording
                    continue
                end
                Mic_Data = RawWav_mic(mic_start:mic_stop1);
            end
        else
            Mic_Data = RawWav_mic(mic_start:mic_stop);
        end
        % filter the data
        Mic_Data = Mic_Data-mean(Mic_Data);
        Filt_MicData = filtfilt(sos_mic_band,1,Mic_Data);
        Filt_LoggerData = filtfilt(sos_logger_band,1,Logger_Data);
        % Calculate the RMS in the microphone extract and check if it's
        % above thershold
        AmpEnv=envelope(Filt_MicData,FS_mic/FS_env,'rms');
        BioSoundUniqParam{ll}{ee}(17)=std(Filt_MicData(Filt_MicData~=0));
        BioSoundUniqParam{ll}{ee}(18) = max(AmpEnv);
        BioSoundUniqParam{ll}{ee}(19) = mean(AmpEnv);
        ResampFilt_MicData = resample(Filt_MicData,4*FHigh,FS_mic);
        ResampFilt_LogData = resample(Filt_LoggerData,4*FHigh,FS_local);
        if length(ResampFilt_MicData)~=length(ResampFilt_LogData)
            OptLength = min(length(ResampFilt_MicData),length(ResampFilt_LogData));
            ResampFilt_MicData = ResampFilt_MicData(1:OptLength);
            ResampFilt_LogData = ResampFilt_LogData(1:OptLength);
        end
        %             figure(10)
        %             clf
        %             title(sprintf('%s event %d/%d', AL_AutoId{ll}, ee, Nevents))
        %             subplot(2,1,1)
        %             [toM, foM, logBM,~] = spec_only_bats(ResampFilt_MicData,4*FHigh,DBNoise, FHigh,100,'Time_increment',0.005);
        %             subplot(2,1,2)
        %             [toL, foL, logBL,~] = spec_only_bats(ResampFilt_LogData,4*FHigh,DBNoise, FHigh,100,'Time_increment',0.005);
        %
        
        % Do a cross correlation between the two signals
        [Xcor,Lag] = xcorr(ResampFilt_MicData,ResampFilt_LogData, (2*Buffer*10^-3)*4*FHigh,'normalized'); % Running a cross correlation between the raw signal and each audio logger signal with a maximum lag equal to twice the Buffer size
        BioSoundUniqParam{ll}{ee}(20) = max(abs(Xcor));
        %kep track of ground truth
        BioSoundUniqParam{ll}{ee}(21) = CorrectAutoDetection01{1}{ll_auto}(ee);
        %         audiowrite(fullfile(Data_out, sprintf('Sound_%s_%d_%d_%d_%d.wav', (sprintf(AL_AutoId{ll_auto})), ee, OnInd, OffInd,CorrectAutoDetection01.(sprintf(AL_AutoId{ll_auto}))(ee) )),Sound, FS_local);
    end
end
BioSoundParamNames = {'MeanSaliency' 'MaxAmp' 'RMS' 'MeanTime' 'StdTime' 'KurtosisTime' 'SkewTime' 'EntropyTime' 'Q1' 'Q2' 'Q3' 'MeanSpec' 'StdSpec' 'KurtosisSpec' 'SkewSpec' 'EntropySpec' 'RMS' 'MaxAmp' 'MeanAmp' 'Corr_wave'};
% Turn back on warnings regarding Pyton to structure conversion
% warning('on', 'MATLAB:structOnObject')
%save(fullfile(Path2Data1, 'SoundEvent.mat'),'BioSoundUniqParam', 'BioSoundParamNames','AL_AutoId','AL_ManId', 'MissedAutoDetection','TotManCall','CorrectAutoDetection01','NLoggers','ManCallTranscTime_ms','ManCallMicSamp','ManCallLogSamp', 'SamplingFreq','ManCallMicFile','-append')
% Concatenate data
for ll=1:NL
    TotEv = length(BioSoundUniqParam{ll});
    for ee=1:TotEv
        if length(BioSoundUniqParam{ll}{ee})<21
            BioSoundUniqParam{ll}{ee} = [BioSoundUniqParam{ll}{ee}; nan(21-length(BioSoundUniqParam{ll}{ee}),1)];
        end
    end
    BioSoundUniqParam{ll} = [BioSoundUniqParam{ll}{:}];
end
BioSoundUniqParam = [BioSoundUniqParam{:}]';
save(fullfile(Path2Results1, 'SoundEvent2.mat'),'BioSoundUniqParam', 'BioSoundParamNames','AL_AutoId','AL_ManId', 'MissedAutoDetection','TotManCall','CorrectAutoDetection01','NLoggers','ManCallTranscTime_ms','ManCallMicSamp','ManCallLogSamp', 'SamplingFreq','ManCallMicFile','Delay','RMSThresh','SoundEvent_LoggerSamp','SoundEvent_TranscTime_ms','CorrectAutoDetection01','MissedAutoDetection','-append')

%% Now loop through the detected elements and calculate/save spectro
WorkingDirSpectro = '/Users/elie/Documents/GroundTruthWorkDir';
mkdir(WorkingDirSpectro)
Buffer = 30;% time in ms to add before after each sound element such atht it's longer than the 23ms required for biosound to calculate fundamental and saliency parameters
F_low = 100;
F_highSpec = 15000;
Flow = 500;
FHigh = 10000;
FBand = 50;
DBNOISE = 50;

[z,p,k] = butter(6,[F_low F_highSpec]/(50000/2),'bandpass');
sos_band_piezo = zp2sos(z,p,k);


FS_env = 1000; %Sampling Frequency of the envelope as calculated by piezo_find_calls_logger

Path2Data1 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audiologgers';
Path2Data2 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audio';
Path2Results1 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audiologgers/GroundTruthResultsPipelineCheck';
load(fullfile(Path2Results1, 'SoundEvent.mat'))
AllLoggers = dir(fullfile(Path2Data1, '*ogger*'));
DirFlags = [AllLoggers.isdir];
% Extract only those that are directories.
AllLoggers = AllLoggers(DirFlags);
NL = length(AllLoggers);

% BioSoundParamNames = {'stdtime' 'meantime' 'skewtime' 'entropytime'...
%         'kurtosistime' 'AmpPeriodF' 'AmpPeriodP' 'rms' 'maxAmp' 'stdspect'...
%         'meanspect' 'skewspect' 'entropyspect' 'kurtosisspect' 'q1' 'q2' 'q3'...
%         'fund' 'cvfund' 'minfund' 'maxfund' 'meansal' '01correct'};
load(fullfile(Path2Results1,'190927_1014_VocExtractData.mat'),'Piezo_wave')
% load(fullfile(Path2Results1,'190603_1447_VocExtractData.mat'),'Piezo_wave')
AL_AutoId = fieldnames(SoundEvent_TranscTime_ms); % Names of the audioLoggers
AL_ManId = fieldnames(Piezo_wave); % Names of the audioLoggers
clear Piezo_wave
SpectroFilename1 = cell(1,NLoggers);
for ll=1:NLoggers
    ll_auto = contains(AL_AutoId, AL_ManId{ll});
    fprintf(1, '*** %s %d/%d ****\n',AL_ManId{ll}, ll, NLoggers)
    % Load the raw signal
    Data_directory = fullfile(AllLoggers(ll).folder,AL_ManId{ll}, 'extracted_data');
    File = dir(fullfile(Data_directory, '*CSC0*'));
    if isempty(File)
        error('Data file not found');
    end
    Filepath = fullfile(File.folder, File.name);
    DataL=load(Filepath, 'AD_count_int16', 'Indices_of_first_and_last_samples','Estimated_channelFS_Transceiver');
    AD_count_double = double(DataL.AD_count_int16);
    DataL.AD_count_int16 = [];
    % Center the signal and clear the old data from memory
    Centered_piezo_signal = AD_count_double - mean(AD_count_double);
    AD_count_double=[];
    
    % convert transceiver time to audio samp files
    Nevents = size(SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto})),1);
    [MicVoc_samp_idx,MicVoc_File]=transc_time2micsamp(Path2Data2,SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto})));
    
    
    % Loop through sound events
    OldMicVoc_File = 0;
    TotEv = size(SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto})),1);
    SpectroFilename1{ll} = cell(1,TotEv);
    for ee=1:TotEv
        if isnan(CorrectAutoDetection01{1}{ll_auto}(ee)) % only keep data for which we have groundtruth
            continue
        end
        if rem(ee,100)==0
            fprintf(1, 'Event %d/%d\n', ee,TotEv)
        end
        
        % find the sampling Frequency
        OnInd1 = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,1);
        OffInd1 = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,2);
        FileIdx = find((DataL.Indices_of_first_and_last_samples(:,1)<OnInd1) .* (DataL.Indices_of_first_and_last_samples(:,2)>OffInd1));
        if isempty(FileIdx) || (length(FileIdx)~=1) || FileIdx>length(DataL.Estimated_channelFS_Transceiver)
            FS_local = round(nanmean(DataL.Estimated_channelFS_Transceiver));
        else
            FS_local = round(DataL.Estimated_channelFS_Transceiver(FileIdx));
        end
        
        % extract the sound with Buffer ms before after the sound
        OnInd_logger = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,1) - round(FS_local*Buffer*10^-3);
        OffInd_logger = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,2) + round(FS_local*Buffer*10^-3);
        if OnInd_logger<0
            OffInd_logger = OffInd_logger-OnInd_logger;
            OnInd_logger=1;
        end
        Logger_Data = Centered_piezo_signal(OnInd_logger : OffInd_logger);
        Logger_Data = Logger_Data - mean(Logger_Data);
        
        % Spectrogram
        Sound_filtered = filtfilt(sos_band_piezo,1,Logger_Data);

        [to, fo, logB, ~] = spec_only_bats(Sound_filtered, FS_local, DBNOISE, F_highSpec, FBand);
%         SpectroFilename{ll}{ee} = fullfile(WorkingDirSpectro,sprintf('20190603_%s_%d_%d.mat',AL_ManId{ll},ee,CorrectAutoDetection01{1}{ll_auto}(ee)));
        SpectroFilename1{ll}{ee} = fullfile(WorkingDirSpectro,sprintf('20190927_%s_%d_%d.mat',AL_ManId{ll},ee,CorrectAutoDetection01{1}{ll_auto}(ee)));
%         save(fullfile(WorkingDirSpectro,sprintf('20190927_%s_%d_%d.mat',AL_ManId{ll},ee,CorrectAutoDetection01{1}{ll_auto}(ee))), 'logB', 'to', 'fo')% CorrectAutoDetection01{1}{ll_auto}(ee) is the ground truth here
        save(SpectroFilename1{ll}{ee}, 'logB', 'to', 'fo','Sound_filtered')% CorrectAutoDetection01{1}{ll_auto}(ee) is the ground truth here
    end
    
end
SpectroFilename1 = [SpectroFilename1{:}];
save(fullfile(WorkingDirSpectro, 'UMAPSpectro.mat'),'SpectroFilename1')


Path2Results1 = '/Volumes/Julie4T/LMC_CoEd/logger/20190603/GroundTruthResultsRecOnly';
Path2Data1 = '/Volumes/Julie4T/LMC_CoEd/logger/20190603';
Path2Data2 = '/Volumes/Julie4T/LMC_CoEd/audio/20190603';
load(fullfile(Path2Results1, 'SoundEvent.mat'))
AllLoggers = dir(fullfile(Path2Data1, '*ogger*'));
DirFlags = [AllLoggers.isdir];
% Extract only those that are directories.
AllLoggers = AllLoggers(DirFlags);
NL = length(AllLoggers);

% BioSoundParamNames = {'stdtime' 'meantime' 'skewtime' 'entropytime'...
%         'kurtosistime' 'AmpPeriodF' 'AmpPeriodP' 'rms' 'maxAmp' 'stdspect'...
%         'meanspect' 'skewspect' 'entropyspect' 'kurtosisspect' 'q1' 'q2' 'q3'...
%         'fund' 'cvfund' 'minfund' 'maxfund' 'meansal' '01correct'};
% load(fullfile(Path2Results1,'190927_1014_VocExtractData.mat'),'Piezo_wave')
load(fullfile(Path2Results1,'190603_1447_VocExtractData.mat'),'Piezo_wave')
AL_AutoId = fieldnames(SoundEvent_TranscTime_ms); % Names of the audioLoggers
AL_ManId = fieldnames(Piezo_wave); % Names of the audioLoggers
clear Piezo_wave
SpectroFilename2 = cell(1,NLoggers);
for ll=1:NLoggers
    ll_auto = contains(AL_AutoId, AL_ManId{ll});
    fprintf(1, '*** %s %d/%d ****\n',AL_ManId{ll}, ll, NLoggers)
    % Load the raw signal
    Data_directory = fullfile(AllLoggers(ll).folder,AL_ManId{ll}, 'extracted_data');
    File = dir(fullfile(Data_directory, '*CSC0*'));
    if isempty(File)
        error('Data file not found');
    end
    Filepath = fullfile(File.folder, File.name);
    DataL=load(Filepath, 'AD_count_int16', 'Indices_of_first_and_last_samples','Estimated_channelFS_Transceiver');
    AD_count_double = double(DataL.AD_count_int16);
    DataL.AD_count_int16 = [];
    % Center the signal and clear the old data from memory
    Centered_piezo_signal = AD_count_double - mean(AD_count_double);
    AD_count_double=[];
    
    % convert transceiver time to audio samp files
    Nevents = size(SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto})),1);
    [MicVoc_samp_idx,MicVoc_File]=transc_time2micsamp(Path2Data2,SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto})));
    
    
    % Loop through sound events
    OldMicVoc_File = 0;
    TotEv = size(SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto})),1);
    SpectroFilename2{ll} = cell(1,TotEv);
    for ee=1:TotEv
        if isnan(CorrectAutoDetection01{1}{ll_auto}(ee)) % only keep data for which we have groundtruth
            continue
        end
        if rem(ee,100)==0
            fprintf(1, 'Event %d/%d\n', ee,TotEv)
        end
        
        % find the sampling Frequency
        OnInd1 = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,1);
        OffInd1 = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,2);
        FileIdx = find((DataL.Indices_of_first_and_last_samples(:,1)<OnInd1) .* (DataL.Indices_of_first_and_last_samples(:,2)>OffInd1));
        if isempty(FileIdx) || (length(FileIdx)~=1) || FileIdx>length(DataL.Estimated_channelFS_Transceiver)
            FS_local = round(nanmean(DataL.Estimated_channelFS_Transceiver));
        else
            FS_local = round(DataL.Estimated_channelFS_Transceiver(FileIdx));
        end
        
        % extract the sound with Buffer ms before after the sound
        OnInd_logger = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,1) - round(FS_local*Buffer*10^-3);
        OffInd_logger = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,2) + round(FS_local*Buffer*10^-3);
        if OnInd_logger<0
            OffInd_logger = OffInd_logger-OnInd_logger;
            OnInd_logger=1;
        end
        Logger_Data = Centered_piezo_signal(OnInd_logger : OffInd_logger);
        Logger_Data = Logger_Data - mean(Logger_Data);
        
        % Spectrogram
        Sound_filtered = filtfilt(sos_band_piezo,1,Logger_Data);

        [to, fo, logB, ~] = spec_only_bats(Sound_filtered, FS_local, DBNOISE, F_highSpec, FBand);
        SpectroFilename2{ll}{ee} = fullfile(WorkingDirSpectro,sprintf('20190603_%s_%d_%d.mat',AL_ManId{ll},ee,CorrectAutoDetection01{1}{ll_auto}(ee)));
%         SpectroFilename{ll}{ee} = fullfile(WorkingDirSpectro,sprintf('20190927_%s_%d_%d.mat',AL_ManId{ll},ee,CorrectAutoDetection01{1}{ll_auto}(ee)));
%         save(fullfile(WorkingDirSpectro,sprintf('20190927_%s_%d_%d.mat',AL_ManId{ll},ee,CorrectAutoDetection01{1}{ll_auto}(ee))), 'logB', 'to', 'fo')% CorrectAutoDetection01{1}{ll_auto}(ee) is the ground truth here
        save(SpectroFilename2{ll}{ee}, 'logB', 'to', 'fo','Sound_filtered')% CorrectAutoDetection01{1}{ll_auto}(ee) is the ground truth here
    end
    
end
SpectroFilename2 = [SpectroFilename2{:}];
SpectroFilename = [SpectroFilename1 SpectroFilename2];
save(fullfile(WorkingDirSpectro, 'UMAPSpectro.mat'),'SpectroFilename2','SpectroFilename', '-append')



%% Construct a matrix of dynamic time warping distance between sound events
NTotEvents = length(SpectroFilename);
DTWmat = nan(NTotEvents, NTotEvents);
for ee1=1:NTotEvents
    TempDTW = cell(1,NTotEvents-ee1);
    Spectro1 = load(fullfile(SpectroFilename{ee1}), 'logB');
    parfor ee2 = (ee1+1) : NTotEvents
        Spectro2 = load(fullfile(SpectroFilename{ee2}), 'logB');
        TempDTW{ee2} = dtw(Spectro1.logB, Spectro2.logB)
    end
    DTWmat(ee1,(ee1+1) : NTotEvents) = [TempDTW{:}];
    DTWmat((ee1+1) : NTotEvents,ee1) = [TempDTW{:}];
end
save(fullfile(WorkingDirSpectro, 'UMAPSpectro.mat'),'DTWmat','-append')

%% Draw some scatters of the parameters
DatasetVoc = BioSoundUniqParam(:,21)==1;
DatasetNoise = BioSoundUniqParam(:,21)==0;

NParam = size(BioSoundUniqParam,2);
for pp=1:(NParam-1)
    figure();
    histogram(BioSoundUniqParam(~BioSoundUniqParam(DatasetNoise,21),pp), 'Normalization', 'probability')
    hold on
    histogram(BioSoundUniqParam(logical(BioSoundUniqParam(DatasetVoc,21)),pp), 'Normalization', 'probability')
    legend('Noise','Vocalizations')
    title(BioSoundParamNames{pp})
    hold off
end

figure()
scatter(BioSoundUniqParam(:,22),BioSoundUniqParam(:,11),10, [BioSoundUniqParam(:,23) zeros(ee_count,1) ones(ee_count,1)],'filled', 'MarkerFaceAlpha',0.5)
xlabel(BioSoundParamNames{22})
ylabel(BioSoundParamNames{11})

figure()
scatter(BioSoundUniqParam(:,15),BioSoundUniqParam(:,16),10, [BioSoundUniqParam(:,23) zeros(ee_count,1) ones(ee_count,1)],'filled', 'MarkerFaceAlpha',0.5)
xlabel(BioSoundParamNames{15})
ylabel(BioSoundParamNames{16})

figure()
scatter(BioSoundUniqParam(:,15),BioSoundUniqParam(:,17),10, [BioSoundUniqParam(:,23) zeros(ee_count,1) ones(ee_count,1)],'filled', 'MarkerFaceAlpha',0.5)
xlabel(BioSoundParamNames{15})
ylabel(BioSoundParamNames{17})

figure()
scatter(BioSoundUniqParam(:,14),BioSoundUniqParam(:,12),10, [BioSoundUniqParam(:,23) zeros(ee_count,1) ones(ee_count,1)],'filled', 'MarkerFaceAlpha',0.5)
xlabel(BioSoundParamNames{14})
ylabel(BioSoundParamNames{12})

figure()
scatter(BioSoundUniqParam(:,11),BioSoundUniqParam(:,9),10, [BioSoundUniqParam(:,23) zeros(ee_count,1) ones(ee_count,1)],'filled', 'MarkerFaceAlpha',0.5)
xlabel(BioSoundParamNames{11})
ylabel(BioSoundParamNames{9})

figure()
scatter(BioSoundUniqParam(:,8),BioSoundUniqParam(:,4),10, [BioSoundUniqParam(:,23) zeros(ee_count,1) ones(ee_count,1)],'filled', 'MarkerFaceAlpha',0.5)
xlabel(BioSoundParamNames{8})
ylabel(BioSoundParamNames{4})

MinSal = 0.1; %Param 22
Q1Max=5000;% Param 15, loose 2 voc: 607, 623
Q2Max = 5000; % Param 16, loose 3: 607 623 682
KurtSpectMax = 400; % Param 14
SkewSpectMax = 8; % Param12, %loose 3 voc: 428, 538, 604
MaxMaxAmp = 5000; % Param 9
MaxMeanSpect = 5000; % Param 11
MinEntropyTime = 0.8; % Param 4 Loose 1 voc: 377
MaxRMS = 1500; % Param 8, loose 1 voc: 682

NoiseRows = logical( (BioSoundUniqParam(:,22)<MinSal) + (BioSoundUniqParam(:,15)>Q1Max) + (BioSoundUniqParam(:,16)>Q2Max) + (BioSoundUniqParam(:,14)>KurtSpectMax) + (BioSoundUniqParam(:,12)>SkewSpectMax) + (BioSoundUniqParam(:,9)>MaxMaxAmp) + (BioSoundUniqParam(:,11)>MaxMeanSpect) + (BioSoundUniqParam(:,4)< MinEntropyTime) + (BioSoundUniqParam(:,8)> MaxRMS));
fprintf(1,'%% False positive before/after restrictions on acoustic parameters: %.1f %%  and %.1f %% \n', sum(~BioSoundUniqParam(:,23))/ee_count*100, sum(~BioSoundUniqParam(~NoiseRows,23))/sum(~NoiseRows)*100)


%% Machine learning approach
% UsefulParams = [22 15 16 17 14 12 9 11 4 8 10];

DataSet190927 = ~isnan(BioSoundUniqParam(:,21));
BioSoundUniqParam_Old = BioSoundUniqParam;
BioSoundUniqParam = BioSoundUniqParam(DataSet190927,:);

UsefulParams = 1:20;
% Try a support vector machine classifier (linear) Binary SVM
SVMModel = fitcsvm(BioSoundUniqParam(:,UsefulParams),BioSoundUniqParam(:,21),'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto','Prior','Uniform');
% Cross-validate the SVM classifier. By default, the software uses 10-fold cross-validation.
CVSVMModel = crossval(SVMModel);
%Estimate the out-of-sample misclassification rate.
classLoss = kfoldLoss(CVSVMModel) % 5% error in cross-validation


% Binary Kernel classification (non-linear)
CVMdl = fitckernel(BioSoundUniqParam(:,UsefulParams),BioSoundUniqParam(:,21),'CrossVal','on','Prior','Uniform')
%CVMdl is a ClassificationPartitionedKernel model. Because fitckernel implements 10-fold cross-validation, CVMdl contains 10 ClassificationKernel models that the software trains on training-fold (in-fold) observations.
%Estimate the cross-validated classification error.
kfoldLoss(CVMdl) % 50.4% error!


% Let's try to predict data using a Binary SVM
ee_count=size(BioSoundUniqParam,1);
oosInds = unique(randi(ee_count,[round(ee_count/10) 1]));   % Out-of-sample indices
isInds = setdiff(1:ee_count, oosInds);   % In-sample indices
X_train = BioSoundUniqParam(isInds,UsefulParams);
Y_train = BioSoundUniqParam(isInds,21);
X_test = BioSoundUniqParam(oosInds,UsefulParams);
Y_test = BioSoundUniqParam(oosInds,21);
%Train an SVM classifier. Standardize the data . Conserve memory by reducing the size of the trained SVM classifier.
SVMModelSplit = fitcsvm(X_train,Y_train,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto','Prior','Uniform');
CompactSVMModelSplit = compact(SVMModelSplit);
whos('SVMModel','CompactSVMModel')

% The CompactClassificationSVM classifier (CompactSVMModel) uses less space than the ClassificationSVM classifier (SVMModel) because SVMModel stores the data.
% Estimate the optimal score-to-posterior-probability transformation function.
CompactSVMModelSplit = fitPosterior(CompactSVMModelSplit,...
    X_train,Y_train)

%The optimal score transformation function (CompactSVMModel.ScoreTransform)
%is a sigmoid  function because the classes are inseparable.
% Predict the out-of-sample labels and class posterior probabilities. Because true labels are available, compare them with the predicted labels.
[labels,PostProbs] = predict(CompactSVMModelSplit,X_test);
figure();
scatter(PostProbs(:,1), PostProbs(:,2), 40, [Y_test zeros(size(Y_test)) zeros(size(Y_test))], 'filled')
hold on
scatter(PostProbs(:,1), PostProbs(:,2), 42, [labels zeros(size(Y_test)) zeros(size(Y_test))])
hold off
xlabel(sprintf('Posterior probability class %d', CompactSVMModelSplit.ClassNames(1)))
ylabel(sprintf('Posterior probability class %d', CompactSVMModelSplit.ClassNames(2)))

fprintf(1,'Percentage of misses (vocalizations detected as noise): %.1f or %d/%d\n', sum(~labels.*Y_test)/length(labels)*100, sum(~labels.*Y_test), length(labels)) % 2.5 %
fprintf(1,'Percentage of false detection (noise detected as vocalizations): %.1f or %d/%d\n', sum(labels.*~Y_test)/length(labels)*100, sum(labels.*~Y_test), length(labels)) % 0.5%

% Now choose a less restrictive label attribution -> any sound with a
% probability of being a vocalization (class 1) above 0.1 is labeled
ProbaThresh = [0:0.001:0.04 0.05:0.05:0.5];
PercMissVoc = nan(length(ProbaThresh),1);
PercFalseDetect = nan(length(ProbaThresh),1);
for pp=1:length(ProbaThresh)
    NewLabels = PostProbs(:,2)>=ProbaThresh(pp);
    %     PercMissVoc(pp) = sum(~NewLabels.*Y_test)/length(NewLabels)*100;
    PercMissVoc(pp) = sum(~NewLabels.*Y_test)/sum(Y_test)*100;
    %     PercFalseDetect(pp) = sum(NewLabels.*~Y_test)/length(NewLabels)*100;
    PercFalseDetect(pp) = sum(NewLabels.*~Y_test)/sum(~Y_test)*100;
end
figure()
plot(ProbaThresh, PercMissVoc,'r','LineWidth',2)
hold on
plot(ProbaThresh, PercFalseDetect, 'k','LineWidth',2)
hold off
xlabel('Threshold on Posterior probability of class 1 (vocalization)')
ylabel('Percentage of error')
legend('Missed Vocalizations', 'False Detection')
% fprintf(1,'With Threshold set at %f Percentage of misses (vocalizations detected as noise): %.1f or %d/%d\n', ProbaThresh(4), PercMissVoc(4), round(length(labels)*PercMissVoc(4)/100), length(labels))
fprintf(1,'With Threshold set at %f Percentage of misses (vocalizations detected as noise): %.1f or %d/%d\n', ProbaThresh(4), PercMissVoc(4), round(sum(Y_test)*PercMissVoc(4)/100), sum(Y_test))% 1.4-1.6%
% fprintf(1,'With Threshold set at %f Percentage of false detection (noise detected as vocalizations): %.1f or %d/%d\n', ProbaThresh(4), PercFalseDetect(4), round(length(labels)*PercFalseDetect(4)/100), length(labels)) %1.4 - 1.9%
fprintf(1,'With Threshold set at %f Percentage of false detection (noise detected as vocalizations): %.1f or %d/%d\n', ProbaThresh(4), PercFalseDetect(4), round(sum(~Y_test)*PercFalseDetect(4)/100), sum(~Y_test))
fprintf(1,'With Threshold set at %f Percentage of misses (vocalizations detected as noise): %.1f or %d/%d\n', ProbaThresh(3), PercMissVoc(3), round(sum(Y_test)*PercMissVoc(3)/100), sum(Y_test))% 1.4-1.6%
fprintf(1,'With Threshold set at %f Percentage of false detection (noise detected as vocalizations): %.1f or %d/%d\n', ProbaThresh(3), PercFalseDetect(3), round(sum(~Y_test)*PercFalseDetect(3)/100), sum(~Y_test))
fprintf(1,'With Threshold set at %f Percentage of misses (vocalizations detected as noise): %.1f or %d/%d\n', ProbaThresh(2), PercMissVoc(2), round(sum(Y_test)*PercMissVoc(2)/100), sum(Y_test))% 1.4-1.6%
fprintf(1,'With Threshold set at %f Percentage of false detection (noise detected as vocalizations): %.1f or %d/%d\n', ProbaThresh(2), PercFalseDetect(2), round(sum(~Y_test)*PercFalseDetect(2)/100), sum(~Y_test))
% If threshold posterior probability set at 0.02, then false positive
% brought down to 2-10% and % of missed vocalizations brought down to 1.7%

%% Save the SVM model for use/prediction with other recordings
Path2Results1 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audiologgers/GroundTruthResultsPipelineCheck';
Data190927 = load(fullfile(Path2Results1, 'SoundEvent2.mat'));
DataSet190927 = ~isnan(Data190927.BioSoundUniqParam(:,21));
Data190927.BioSoundUniqParam = Data190927.BioSoundUniqParam(DataSet190927,:);
Path2Results1 = '/Volumes/Julie4T/LMC_CoEd/logger/20190603/GroundTruthResultsRecOnly';
Data190603 = load(fullfile(Path2Results1, 'SoundEvent2.mat'));
DataSet190603 = ~isnan(Data190603.BioSoundUniqParam(:,21));
Data190603.BioSoundUniqParam = Data190603.BioSoundUniqParam(DataSet190603,:);

BioSoundUniqParam = [Data190603.BioSoundUniqParam; Data190927.BioSoundUniqParam];
BioSoundParamNames = Data190927.BioSoundParamNames;

UsefulParams = 1:20;
% Try a support vector machine classifier (linear) Binary SVM
SVMModel = fitcsvm(BioSoundUniqParam(:,UsefulParams),BioSoundUniqParam(:,21),'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto','Prior','Uniform');
% Cross-validate the SVM classifier. By default, the software uses 10-fold cross-validation.
CVSVMModel = crossval(SVMModel);
%Estimate the out-of-sample misclassification rate.
classLoss = kfoldLoss(CVSVMModel)
CompactSVMModel = compact(SVMModel);
whos('SVMModel','CompactSVMModel')

% The CompactClassificationSVM classifier (CompactSVMModel) uses less space than the ClassificationSVM classifier (SVMModel) because SVMModel stores the data.
% Estimate the optimal score-to-posterior-probability transformation function.
CompactSVMModel = fitPosterior(CompactSVMModel,...
    BioSoundUniqParam(:,UsefulParams),BioSoundUniqParam(:,21))
save('/Users/elie/Documents/CODE/SoundAnalysisBats/SVMModelNoiseVoc_MicLog.mat', 'CompactSVMModel')



UsefulParams = 1:16;
% Try a support vector machine classifier (linear) Binary SVM
SVMModel = fitcsvm(BioSoundUniqParam(:,UsefulParams),BioSoundUniqParam(:,21),'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto','Prior','Uniform');
% Cross-validate the SVM classifier. By default, the software uses 10-fold cross-validation.
CVSVMModel = crossval(SVMModel);
%Estimate the out-of-sample misclassification rate.
classLoss = kfoldLoss(CVSVMModel)
CompactSVMModel = compact(SVMModel);
whos('SVMModel','CompactSVMModel')

% The CompactClassificationSVM classifier (CompactSVMModel) uses less space than the ClassificationSVM classifier (SVMModel) because SVMModel stores the data.
% Estimate the optimal score-to-posterior-probability transformation function.
CompactSVMModel = fitPosterior(CompactSVMModel,...
    BioSoundUniqParam(:,UsefulParams),BioSoundUniqParam(:,21))
save('/Users/elie/Documents/CODE/SoundAnalysisBats/SVMModelNoiseVoc.mat', 'CompactSVMModel')


%% Test a UMAP projection on this dataset
addpath /Users/elie/Documents/CODE/umapDistribution/umap
addpath /Users/elie/Documents/CODE/umapDistribution/util
javaaddpath('/Users/elie/Documents/CODE/umapDistribution/umap/umap.jar')
UsefulParams = [1:16 21];
fprintf(1,'\n\n ***** Distance Euclidean *****\n')
RandomSet = randperm(size(BioSoundUniqParam,1));
TrainingSet = RandomSet(1:round(0.8*size(BioSoundUniqParam,1)));
TestingSet = RandomSet((round(0.8*size(BioSoundUniqParam,1))+1) :end);
% [Reduction,UMAP,ClustID]= run_umap(BioSoundUniqParam(TrainingSet,UsefulParams),'parameter_names',BioSoundParamNames(UsefulParams(1:(end-1))),'label_column',length(UsefulParams),'save_template_file',fullfile(Path2Results1,'UMAP_templateNoiseVoc.mat'), 'metric','euclidean','label_file', fullfile(Path2Results1,'NoiseCallLabels.properties.txt'), 'target_weight',0.6);
[Reduction,UMAP]= run_umap(BioSoundUniqParam(TrainingSet,UsefulParams),'parameter_names',BioSoundParamNames(UsefulParams(1:(end-1))),'label_column',length(UsefulParams),'save_template_file',fullfile(Path2Results1,'UMAP_templateNoiseVoc.mat'), 'metric','euclidean', 'target_weight',0.6);
figure()
scatter(Reduction(:,1), Reduction(:,2),5,[BioSoundUniqParam(TrainingSet,21) zeros(size(Reduction,1),2)],'filled')
title('Euclidean distance')

[ReductionTest,UMAPTest, ClustIDTest]= run_umap(BioSoundUniqParam(TestingSet,UsefulParams(1:end-1)),'parameter_names',BioSoundParamNames(UsefulParams(1:(end-1))),'template_file',fullfile(Path2Results1,'UMAP_templateNoiseVoc.mat'), 'metric','euclidean','match_supervisors',1,'qf_tree',true,'qf_dissimilarity', true);
figure()
scatter(ReductionTest(:,1), ReductionTest(:,2),5,[BioSoundUniqParam(TestingSet,21) zeros(size(ReductionTest,1),2)],'filled')
title('Euclidean distance')



[Reduction,UMAP,ClustID]= run_umap(BioSoundUniqParam(TrainingSet,UsefulParams),'parameter_names',BioSoundParamNames(UsefulParams(1:(end-1))),'label_column',length(UsefulParams),'save_template_file',fullfile(Path2Results1,'UMAP_templateNoiseVoc.mat'), 'metric','euclidean','target_weight',0.6,'n_component',3);
figure()
scatter3(Reduction(:,1), Reduction(:,2),Reduction(:,3),5,[BioSoundUniqParam(TrainingSet,21) zeros(size(Reduction,1),2)],'filled')
title('Euclidean distance')



[Reduction,UMAP,ClustID]= run_umap(BioSoundUniqParam(:,UsefulParams),'metric','cosine');
figure()
scatter(Reduction(:,1), Reduction(:,2),5,[BioSoundUniqParam(:,21) zeros(size(Reduction,1),2)],'filled')
title('Cosyne metric')

[Reduction,UMAP,ClustID]= run_umap(BioSoundUniqParam(:,UsefulParams),'metric','correlation');
figure()
scatter(Reduction(:,1), Reduction(:,2),5,[BioSoundUniqParam(:,17) zeros(size(Reduction,1),2)],'filled')
title('Correlation metric')

%% Try with UMAP on acoustic features with same size noise/voc dataset for traning purposes
NoiseIndices = find(~BioSoundUniqParam(:,17));
SelectNoiseInd = NoiseIndices(randperm(sum(~BioSoundUniqParam(:,17)),sum(BioSoundUniqParam(:,17))));
VocIndices = find(BioSoundUniqParam(:,17));
[Reduction,UMAP,ClustID]= run_umap(BioSoundUniqParam([SelectNoiseInd VocIndices],UsefulParams));
figure()
scatter(Reduction(:,1), Reduction(:,2),5,[BioSoundUniqParam([SelectNoiseInd VocIndices],17) zeros(size(Reduction,1),2)],'filled')
title('Euclidean distance')


%% Test the calculations of spectrograms and apply UMAP with HBDSCAN as a clustering algorithm
%% Now loop through the detected elements and calculate sound duration
SoundEventDuration_ms = nan(21654,1);
ee_count = 0;
AL_AutoId = fieldnames(SoundEvent_TranscTime_ms); % Names of the audioLoggers
AL_ManId = fieldnames(Piezo_wave); % Names of the audioLoggers

for ll=1:NLoggers
    ll_auto = contains(AL_AutoId, AL_ManId{ll});
    fprintf(1, '*** %s %d/%d ****\n',AL_ManId{ll}, ll, NLoggers)
    % Load the raw signal
    Data_directory = fullfile(AllLoggers(ll).folder,AL_ManId{ll}, 'extracted_data');
    File = dir(fullfile(Data_directory, '*CSC0*'));
    if isempty(File)
        error('Data file not found');
    end
    Filepath = fullfile(File.folder, File.name);
    load(Filepath,'Indices_of_first_and_last_samples','Estimated_channelFS_Transceiver')
    
    % Loop through sound events
    TotEv = size(SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto})),1);
    for ee=1:TotEv
        ee_count = ee_count+1;
        if rem(ee,100)==0
            fprintf(1, 'Event %d/%d\n', ee,TotEv)
        end
        
        % find the sampling Frequency
        FileIdx = find((Indices_of_first_and_last_samples(:,1)<OnInd) .* (Indices_of_first_and_last_samples(:,2)>OffInd));
        if isempty(FileIdx) || (length(FileIdx)~=1) || FileIdx>length(Estimated_channelFS_Transceiver)
            FS_local = round(nanmean(Estimated_channelFS_Transceiver));
        else
            FS_local = round(Estimated_channelFS_Transceiver(FileIdx));
        end
        SoundEventDuration_ms(ee) = diff(SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,:))/FS_local*10^3;
    end
end
figure()
subplot(3,1,1)
histogram(SoundEventDuration_ms, 'BinWidth',1)
xlabel('Sound event duration ms')
subplot(3,1,2)
histogram(SoundEventDuration_ms(SoundEventDuration_ms<500), 'BinWidth',1)
xlabel('Sound event duration ms')
subplot(3,1,3)
histogram(SoundEventDuration_ms(SoundEventDuration_ms<200), 'BinWidth',1)
xlabel('Sound event duration ms')

%% Find calls using call detection on microphone
%load previous data
load(fullfile(Path2Results1, 'SoundEvent.mat'))
load(fullfile(Path2Results1, 'SoundEvent2.mat'),'TotManCall','NLoggers')
addpath(genpath('/Users/elie/Documents/CODE/neurobat-callCutting'))
FigOn = 0; % set to 1 to see the spectrogram of each vocalization with the coresponding onset/offset from the automatic and manual detection
% WD = '/Users/elie/Documents/ManipBats/LMC/20190603';
% FS = 192000;
% findcalls_session(WD,FS,'fileType','wav')
Thresh = [0.2 0.25 0.5 0.6 0.7 0.75].*10^-3; % Testing a variety of threshold on the amplitude envelope for the detection of calls in microphone data
%Thresh = 0.2*10^-3;
MissedMicAutoDetection = cell(length(Thresh),1);
MicCorrectAutoDetection01 =  cell(length(Thresh),1);
TotMissedMicCall = nan(length(Thresh),1);

for tt=1:length(Thresh)
    fprintf(1,'*** Thresh %d/%d ****',tt,length(Thresh))
    FS = 192000;
    findcalls_session(Path2Data2,FS,'fileType','wav','audio_file_filter','JuGr_190927_1014_RecOnly_mic1*','filter_raw_data',true,'thresh',Thresh(tt))
    
    % convert automatically detected Audio time to transceiver times
    [MicVoc_samp_idx,MicVoc_transcTime_ms] = mic2transc_time(Path2Data2);
    
    % save(fullfile(Path2Results1, 'SoundEvent3.mat'),'MicVoc_samp_idx','MicVoc_transcTime_ms')
    %% Let's loop in the dataset of manually extracted calls (Compare to ground truth) and see how many correct hits we get in the automatic detection
    
    Delay = 100; %in ms error/delay between auto and man detection and Delay to add before each detected call in ms
    MissedMicAutoDetection{tt} = cell(NLoggers,1);
    TotManCall = 0;
    MicCorrectAutoDetection01{tt} = zeros(size(MicVoc_transcTime_ms,1),1);
    [z,p,k] = butter(6,[100 90000]/(192000/2),'bandpass');
    sos_raw_band = zp2sos(z,p,k);
    
    for ll=1:NLoggers
        MissedMicAutoDetection{tt}{ll} = [];
        
        TotManCall = TotManCall+size(ManCallTranscTime_ms{ll},1);
        
        for vv=1:size(ManCallTranscTime_ms{ll},1)
            fprintf(1, '*** %s vocalization %d/%d ****\n', AL_ManId{ll}, vv,size(ManCallTranscTime_ms{ll},1))
            OnOffVoc = ManCallTranscTime_ms{ll}(vv,:);
            Idx_OnsetAuto = find((MicVoc_transcTime_ms(:,1)>(OnOffVoc(1)-Delay)) .* (MicVoc_transcTime_ms(:,1)<(OnOffVoc(2)+Delay))>0);
            Idx_OffsetAuto = find((MicVoc_transcTime_ms(:,2)>(OnOffVoc(1)-Delay)) .* (MicVoc_transcTime_ms(:,2)<(OnOffVoc(2)+Delay))>0);
            
            if isempty(Idx_OnsetAuto) && isempty(Idx_OffsetAuto)
                Idx_OnsetAuto = find((MicVoc_transcTime_ms(:,1)<(OnOffVoc(1)-Delay)) .* (MicVoc_transcTime_ms(:,2)>(OnOffVoc(1)-Delay))>0);
                Idx_OffsetAuto = find((MicVoc_transcTime_ms(:,1)<(OnOffVoc(2)+Delay)) .* (MicVoc_transcTime_ms(:,2)>(OnOffVoc(2)+Delay))>0);
            end
            if isempty(Idx_OnsetAuto) && isempty(Idx_OffsetAuto)
                fprintf('No automatic call detected\n')
                disp(OnOffVoc)
                MissedMicAutoDetection{tt}{ll} = [MissedMicAutoDetection{tt}{ll} vv];
                if FigOn
                    figure(1)
                    clf
                    
                    % Plot the spectrogram of the sound extract if requested
                    DBNoise = 60; % amplitude parameter for the color scale of the spectro
                    FHigh = 10000; % y axis max scale for the spectrogram
                    [Raw, FSRaw] = audioread(fullfile(Path2Data2,[ManCallMicFile{ll}{vv} '.wav']));
                    x_start = round(ManCallMicSamp{ll}(vv,1)-Delay*10^(-3)*FSRaw);
                    x_stop = round(ManCallMicSamp{ll}(vv,2)+Delay*10^(-3)*FSRaw);
                    Raw = Raw(x_start:x_stop);
                    Filt_RawVoc = filtfilt(sos_raw_band,1,Raw);
                    Raw_ramp = cosramp(Filt_RawVoc-mean(Filt_RawVoc), FSRaw*10*10^-3);
                    [~] = spec_only_bats(Raw_ramp,FSRaw,DBNoise, FHigh);
                    %             caxis('manual');
                    %             caxis([2 70]);
                    ylim([-500 10000])
                    hold on
                    
                    yyaxis right
                    ylim([-10 300])
                    AllignOn = ManCallMicSamp{ll}(vv,1)/FSRaw*10^3-Delay;
                    x_start_man = Delay;
                    x_stop_man = round(ManCallMicSamp{ll}(vv,2)/FSRaw*10^3-AllignOn);
                    plot([x_start_man x_stop_man], ones(2,1) * 1.5, 'g-', 'LineWidth',1)
                    hold off
                    
                    
                    title(sprintf('%s  Voc %d/%d NOT DETECTED', AL_ManId{ll}, vv, size(ManCallTranscTime_ms{ll},1)))
                    Player = audioplayer(Raw_ramp/std(Raw_ramp), FSRaw);
                    play(Player)
                    %             pause()
                end
                
            else
                IdxAll = union(Idx_OnsetAuto,Idx_OffsetAuto);
                MicCorrectAutoDetection01{tt}(IdxAll) = ones(size(IdxAll));
                fprintf('%d events started during that vocalization\n',length(IdxAll))
                if FigOn
                    figure(1)
                    clf
                    
                    % Plot the spectrogram of the sound extract
                    DBNoise = 60; % amplitude parameter for the color scale of the spectro
                    FHigh = 10000; % y axis max scale for the spectrogram
                    [Raw, FSRaw] = audioread(fullfile(Path2Data2,[ManCallMicFile{ll}{vv} '.wav']));
                    x_start = round(ManCallMicSamp{ll}(vv,1)-Delay*10^(-3)*FSRaw);
                    x_stop = round(ManCallMicSamp{ll}(vv,2)+Delay*10^(-3)*FSRaw);
                    Raw = Raw(x_start:x_stop);
                    Filt_RawVoc = filtfilt(sos_raw_band,1,Raw);
                    Raw_ramp = cosramp(Filt_RawVoc-mean(Filt_RawVoc), FSRaw*10*10^-3);
                    [~] = spec_only_bats(Raw_ramp,FSRaw,DBNoise, FHigh);
                    %             caxis('manual');
                    %             caxis([2 70]);
                    ylim([-500 10000])
                    hold on
                    
                    yyaxis right
                    ylim([-10 300])
                    AllignOn = ManCallMicSamp{ll}(vv,1)/FSRaw*10^3-Delay;
                    x_start_man = Delay;
                    x_stop_man = round(ManCallMicSamp{ll}(vv,2)/FSRaw*10^3-AllignOn);
                    plot([x_start_man x_stop_man], ones(2,1) * 1.5, 'g-', 'LineWidth',1)
                    hold on
                    for ii=1:length(IdxAll)
                        OnOff = (MicVoc_samp_idx(IdxAll(ii),:)/FSRaw*10^3 - AllignOn);
                        plot(OnOff, ones(2,1) * (-1.5), 'b-', 'LineWidth',2)
                        if ii==1
                            legend({ 'Manual' ,'Automatic'})
                            legend('AutoUpdate','Off')
                        end
                        hold on
                    end
                    hold off
                    
                    
                    title(sprintf('%s  Voc %d/%d', AL_ManId{ll}, vv, size(ManCallTranscTime_ms{ll},1)))
                    %             Player = audioplayer(Raw_ramp/std(Raw_ramp), FSRaw);
                    %             play(Player)
                    %             pause()
                end
            end
        end
    end
    TotMissedMicCall(tt) = sum(cellfun('length',MissedMicAutoDetection{tt}));
    fprintf(1,'Missed call by the microphone auto detection: %d/%d or %.1f%%\n', sum(cellfun('length',MissedMicAutoDetection{tt})),TotManCall,sum(cellfun('length',MissedMicAutoDetection{tt}))/TotManCall*100);
    fprintf(1,'Number of detected events from the microphone %d, proportion of true calls %d/%d or %.1f%%\n',length(MicCorrectAutoDetection01{tt}),sum(MicCorrectAutoDetection01{tt}),length(MicCorrectAutoDetection01{tt}),sum(MicCorrectAutoDetection01{tt})/length(MicCorrectAutoDetection01{tt})*100);
end
figure()
subplot(2,1,1)
plot(TotMissedMicCall, 'LineWidth',2)
hold on
%  plot(cellfun(@sum,MicCorrectAutoDetection01), 'LineWidth',2)
%  hold on
%  legend('# Missed calls')
set(gca,'XTick', 1:length(Thresh))
set(gca,'XTickLabel', Thresh)
xlabel('Threshold on mic envelope')
ylabel(sprintf('Number of missed calls out of %d', TotManCall))
hold off
ylim([0 100])
title('Microphone detection')
subplot(2,1,2)
plot(TotMissedMicCall/TotManCall*100, 'LineWidth',2)
hold on
plot(cellfun(@sum,MicCorrectAutoDetection01)./cellfun('length',MicCorrectAutoDetection01).*100, 'LineWidth',2)
hold on
ylabel('Percentage')
yyaxis right
plot(100-cellfun(@sum,MicCorrectAutoDetection01)./cellfun('length',MicCorrectAutoDetection01).*100, 'LineWidth',2)
legend('missed calls', 'correct detection', 'false detection')
set(gca,'XTick', 1:length(Thresh))
set(gca,'XTickLabel', Thresh)
xlabel('Threshold on mic envelope')
ylabel('Percentage')
ylim([80 100])
yyaxis left
ylim([0 20])
hold off

%% Internal functions
function BiosoundObj = runBiosound(Y, FS, F_high)
% Hard coded parameters for biosound
% spectrogram parameters
Spec_sample_rate = 1000; % sampling rate Hz
Freq_spacing = 50; % width of the frequency window for the FFT Hz
Min_freq = 300; % high pass filter before FFT Hz
Max_freq = 50000; % Low pass filter before FFT Hz
% temporal enveloppe parameters
Cutoff_freq = 150; % Hz
Amp_sample_rate = 1000; % Hz
if nargin<3
    % Spectrum parameters
    F_high = 50000; % frequency of Low-pass filter Hz
end
% Fundamental parameters
MaxFund = 4000;
MinFund = 300;
LowFc = 100; %100
HighFc = 18000;% 15000
MinSaliency = 0.6;
DebugFigFundest = 0;
MinFormantFreq = 2000;
MaxFormantBW = 1000; %500
WindowFormant = 0.1;
Method= 'Stack';

% create the biosound object
BiosoundObj = py.soundsig.sound.BioSound(py.numpy.array(Y),pyargs('fs',FS));
% methods(BiosoundFi, '-full') % this command plot all the methods with the available arguments

% Calculate the RMS (lhs std(varargin))
BiosoundObj.rms = BiosoundObj.sound.std();

% calculate the amplitude enveloppe
ampenv(BiosoundObj, Cutoff_freq,Amp_sample_rate);

% Calculate the periodicity of the amplitude envelope
SoundAmp = double(py.array.array('d', py.numpy.nditer(BiosoundObj.amp)));
[P,F] = pspectrum(SoundAmp,1000);
[PKS,LOCS]=findpeaks(P);
AmpPeriodF = F(LOCS(PKS == max(PKS))); % Frequency in hertz of the max peak
AmpPeriodP = max(PKS)/mean(SoundAmp.^2); % Proportion of power in the max peak of the spectrum

% calculate the spectrum (lhs spectrum(self, f_high, pyargs))
spectrum(BiosoundObj, F_high)
% calculate the spectrogram (lhs spectroCalc(self, spec_sample_rate,
% freq_spacing, min_freq, max_freq, pyargs))
try % For very short sound, the Freq_spacing is too small, doubling if error
    spectroCalc(BiosoundObj, Spec_sample_rate, Freq_spacing, Min_freq,Max_freq)
catch
    try
        spectroCalc(BiosoundObj, Spec_sample_rate, Freq_spacing.*2, Min_freq,Max_freq)
    catch
        warning('Impossible to calculate spectrogram')
        BiosoundObj.spectro = nan;
    end
end

% Calculate time varying spectralmean and spectral max
Spectro = double(BiosoundObj.spectro);
if ~isnan(Spectro)
    Fo = double(BiosoundObj.fo);
    TPoints = size(Spectro,2);
    SpectralMean = nan(1,TPoints);
    %         SpectralMax = nan(1,TPoints);
    for tt=1:TPoints
        %             SpectralMax(tt) = Fo(Spectro(:,tt)==max(Spectro(:,tt)));
        PSDSpec = Spectro(:,tt)./(sum(Spectro(:,tt)));
        SpectralMean(tt) = sum(PSDSpec' .* Fo);
    end
else
    SpectralMean = nan;
end

% calculate the fundamental and related values (lhs fundest(self, maxFund,
% minFund, lowFc, highFc, minSaliency, debugFig, pyargs)
fundest(BiosoundObj, MaxFund, MinFund,LowFc, HighFc, MinSaliency,DebugFigFundest,MinFormantFreq,MaxFormantBW,WindowFormant,Method)

% convert biosound to a strcuture
BiosoundObj = struct(BiosoundObj);
% Add some fields
BiosoundObj.AmpPeriodF = AmpPeriodF;
BiosoundObj.AmpPeriodP = AmpPeriodP;
BiosoundObj.SpectralMean = SpectralMean;
%         BiosoundObj.SpectralMax = SpectralMax;
% convert all nmpy arrays to double to be able to save as matfiles
BiosoundObj.amp = SoundAmp;
BiosoundObj.tAmp = double(BiosoundObj.tAmp);
BiosoundObj.spectro = double(BiosoundObj.spectro);
BiosoundObj.to = double(BiosoundObj.to);
BiosoundObj.fo = double(BiosoundObj.fo);
BiosoundObj.F1 = double(BiosoundObj.F1);
BiosoundObj.F2 = double(BiosoundObj.F2);
BiosoundObj.F3 = double(BiosoundObj.F3);
BiosoundObj.fpsd = double(BiosoundObj.fpsd);
BiosoundObj.psd = double(BiosoundObj.psd);
BiosoundObj.sal = double(BiosoundObj.sal);
BiosoundObj.f0 = double(BiosoundObj.f0);
BiosoundObj.f0_2 = double(BiosoundObj.f0_2);
BiosoundObj.fund = double(BiosoundObj.fund);
BiosoundObj.cvfund = double(BiosoundObj.cvfund);
BiosoundObj.fund2 = double(BiosoundObj.fund2);
BiosoundObj.minfund = double(BiosoundObj.minfund);
BiosoundObj.maxfund = double(BiosoundObj.maxfund);
BiosoundObj.sound = double(BiosoundObj.sound);
BiosoundObj.wf = double(BiosoundObj.wf);
BiosoundObj.wt = double(BiosoundObj.wt);
BiosoundObj.mps = double(BiosoundObj.mps);
end

function [MicVoc_samp_idx,MicAuto_transcTime]=mic2transc_time(RawWav_dir)
AllFiles = dir(fullfile(RawWav_dir, 'Analyzed_auto','*.mat'));
% find the date and expstart time
Date = AllFiles(1).name(6:11);
ExpStartTime = AllFiles(1).name(13:16);
TTL_dir = dir(fullfile(RawWav_dir,sprintf( '%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
TTL = load(fullfile(TTL_dir.folder, TTL_dir.name));
% loop through data
Nfiles = length(AllFiles);
MicAuto_transcTime=nan(Nfiles,2);
MicVoc_samp_idx=nan(Nfiles,2);
% Extract the transceiver time
% zscore the sample stamps
for ff=1:Nfiles
    Ind_=strfind(AllFiles(ff).name,'_');
    FileIdx=str2double(AllFiles(ff).name((Ind_(5)+1):(Ind_(6)-1)));
    Data = load(fullfile(AllFiles(ff).folder, AllFiles(ff).name));
    TTL_idx = find(unique(TTL.File_number) == FileIdx);
    MicVoc_samp_idx(ff,:) = Data.callpos;
    Voc_samp_idx_zs = (Data.callpos - TTL.Mean_std_Pulse_samp_audio(TTL_idx,1))/TTL.Mean_std_Pulse_samp_audio(TTL_idx,2);
    % calculate the transceiver times
    MicAuto_transcTime(ff,:) = TTL.Mean_std_Pulse_TimeStamp_Transc(TTL_idx,2) .* polyval(TTL.Slope_and_intercept{TTL_idx},Voc_samp_idx_zs,[], TTL.Mean_std_x{TTL_idx}) + TTL.Mean_std_Pulse_TimeStamp_Transc(TTL_idx,1);
    callpos_transcTime = MicAuto_transcTime(ff,:);
    save(fullfile(AllFiles(ff).folder, AllFiles(ff).name), 'callpos_transcTime', '-append')
end
end


function [MicVoc_samp_idx,MicVoc_File]=transc_time2micsamp(RawWav_dir,OnOffTranscTime_ms)
AllFiles = dir(fullfile(RawWav_dir, '*RecOnly*.wav'));
% find the date and expstart time
Date = AllFiles(1).name(6:11);
ExpStartTime = AllFiles(1).name(13:16);
TTL_dir = dir(fullfile(RawWav_dir,sprintf( '%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
TTL = load(fullfile(TTL_dir.folder, TTL_dir.name));

% loop through data
Nevents = size(OnOffTranscTime_ms,1);
MicVoc_samp_idx=nan(Nevents,2);
MicVoc_File = nan(Nevents,1);

% Extract the transceiver time
% zscore the sample stamps
for ff=1:Nevents
    FileNumIdx = find(TTL.Pulse_TimeStamp_Transc<OnOffTranscTime_ms(ff,1),1,'Last');
    if isempty(FileNumIdx)
        FileNumIdx = find(TTL.Pulse_TimeStamp_Transc>OnOffTranscTime_ms(ff,1),1,'First');
    end
    MicVoc_File(ff) = TTL.File_number(FileNumIdx);
    TranscTime_zs = (OnOffTranscTime_ms(ff,:) - TTL.Mean_std_Pulse_TimeStamp_Transc(MicVoc_File(ff),1))/TTL.Mean_std_Pulse_TimeStamp_Transc(MicVoc_File(ff),2);
    MicVoc_samp_idx(ff,:) =TTL.Mean_std_Pulse_samp_audio(MicVoc_File(ff),2) .* polyval(TTL.Slope_and_intercept_transc2audiosamp{MicVoc_File(ff)},TranscTime_zs,[],TTL.Mean_std_x_transc2audiosamp{MicVoc_File(ff)}) + TTL.Mean_std_Pulse_samp_audio(MicVoc_File(ff),1);
end
end

