%% Load input data
% Path2Data = '/Volumes/server_home/users/JulieE/LMC_CoEd/logger/20190623';
% Path2Audio = '/Volumes/server_home/users/JulieE/LMC_CoEd/audio/20190623';
% % Load data manually extracted
% load(fullfile(Path2Data, '190623_1401_VocExtractData_200.mat'))
% load(fullfile(Path2Data, '190623_1401_VocExtractData.mat'))
% load(fullfile(Path2Audio, '190623_1401_VocExtractTimes.mat'))
addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'))

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
% load(fullfile(Path2Data1, 'SoundEvent.mat'))

%% Run the automatic detection
AllLoggers = dir(fullfile(Path2Data1, '*ogger*'));
DirFlags = [AllLoggers.isdir];
% Extract only those that are directories.
AllLoggers = AllLoggers(DirFlags);
NL = length(AllLoggers);
for ll=1:NL
    Data_directory = fullfile(AllLoggers(ll).folder,AllLoggers(ll).name, 'extracted_data');
    [SoundEvent_LoggerSamp.(sprintf('L%s',AllLoggers(ll).name(2:end))),SoundEvent_TranscTime_ms.(sprintf('L%s',AllLoggers(ll).name(2:end))),LoggerEnvelopeAll.(sprintf('L%s',AllLoggers(ll).name(2:end))),SoundEvent_EnvSamp.(sprintf('L%s',AllLoggers(ll).name(2:end))),~] = piezo_find_calls_logger(Data_directory);
end
save(fullfile(Path2Results1, 'SoundEvent.mat'),'SoundEvent_LoggerSamp','SoundEvent_TranscTime_ms','LoggerEnvelopeAll')
% load(fullfile(Path2Data1, 'SoundEvent.mat')



% Duration of the data manually analysed 100 min
% ManDur = 10*60*10^3; %in ms
ManDur = 10*10*60*10^3; %in ms

%% Get in transceiver time the onset/offset of each vocalization manually extracted for each logger
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
save(fullfile(Path2Results1, 'SoundEvent.mat'),'AL_ManId','NLoggers','ManCallTranscTime_ms','ManCallMicSamp','ManCallLogSamp', 'SamplingFreq','ManCallMicFile','-append')
%% Let's loop in the dataset of manually extracted calls and see how many correct hits we get in the automatic detection
FS_env = 1000; %Sampling Frequency of the envelope as calculated by piezo_find_calls_logger
Delay = 100; %in ms error/delay between auto and man detection and Delay to add before each detected call in ms
AL_AutoId = fieldnames(SoundEvent_TranscTime_ms); % Names of the audioLoggers
MissedAutoDetection = cell(NLoggers,1);
TotManCall = 0;
CorrectAutoDetection01 = struct();

for ll=1:NLoggers
    MissedAutoDetection{ll} = [];
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
    CorrectAutoDetection01.(sprintf(AL_AutoId{ll_auto})) = zeros(size(SoundEvent_TranscTime_ms.(sprintf(AL_AutoId{ll_auto})),1),1);
    
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
            MissedAutoDetection{ll} = [MissedAutoDetection{ll} vv];
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
            
            
%             title(sprintf('%s  Voc %d/%d NOT DETECTED', AL_ManId{ll}, vv, size(ManCallTranscTime_ms{ll},1)))
%             Player = audioplayer(Raw_ramp/std(Raw_ramp), SamplingFreq{ll});
%             play(Player)
%             pause()
            
            
        else
            IdxAll = union(Idx_OnsetAuto,Idx_OffsetAuto);
            CorrectAutoDetection01.(sprintf(AL_AutoId{ll_auto}))(IdxAll) = ones(size(IdxAll));
            fprintf('%d events started during that vocalization\n',length(IdxAll))
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


%% Now loop through the detected elements and calculate biosound
Buffer = 30;% time in ms to add before after each sound element such atht it's longer than the 23ms required for biosound to calculate fundamental and saliency parameters
F_High = 5000;
F_low = 100;
F_highSpec = 15000;
% BioSoundUniqParam = nan(21553,23);
BioSoundUniqParam = nan(21553,17);
ee_count = 0;
% BioSoundParamNames = {'stdtime' 'meantime' 'skewtime' 'entropytime'...
%         'kurtosistime' 'AmpPeriodF' 'AmpPeriodP' 'rms' 'maxAmp' 'stdspect'...
%         'meanspect' 'skewspect' 'entropyspect' 'kurtosisspect' 'q1' 'q2' 'q3'...
%         'fund' 'cvfund' 'minfund' 'maxfund' 'meansal' '01correct'};
        
    
% Turn off warnings regarding Pyton to structure conversion
% warning('off', 'MATLAB:structOnObject')

AL_AutoId = fieldnames(SoundEvent_TranscTime_ms); % Names of the audioLoggers
AL_ManId = fieldnames(Piezo_wave); % Names of the audioLoggers
BiosoundFolder = cell(NLoggers,1);
Data_out = fullfile(AllLoggers(ll).folder, 'BiosoundEvents');
if ~exist(Data_out, 'dir')
    mkdir(Data_out)
end
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
    load(Filepath, 'AD_count_int16', 'Indices_of_first_and_last_samples','Estimated_channelFS_Transceiver')
    AD_count_double = double(AD_count_int16);
    clear AD_count_int16
    % Center the signal and clear the old data from memory
    Centered_piezo_signal = AD_count_double - mean(AD_count_double);
    clear AD_count_double
    
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
        
        % extract the sound with Buffer ms before after the sound
        OnInd = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,1) - round(FS_local*Buffer*10^-3);
        OffInd = SoundEvent_LoggerSamp.(sprintf(AL_AutoId{ll_auto}))(ee,2) + round(FS_local*Buffer*10^-3);
        Sound = Centered_piezo_signal(OnInd : OffInd);
        Sound = Sound - mean(Sound);
        
%         BioSoundCall=runBiosound(Sound, FS_local, F_High);
%         
%         % Feed biosound data into a Matrix
%             % temporal parameters (calculated on the envelope)
%             BioSoundUniqParam(ee_count,1) = BioSoundCall.stdtime;
%             BioSoundUniqParam(ee_count,2) = BioSoundCall.meantime;
%             BioSoundUniqParam(ee_count,3) = BioSoundCall.skewtime;
%             BioSoundUniqParam(ee_count,4) = BioSoundCall.entropytime;
%             BioSoundUniqParam(ee_count,5) = BioSoundCall.kurtosistime;
%             if ~isempty(BioSoundCall.AmpPeriodF)
%                 BioSoundUniqParam(ee_count,6) = BioSoundCall.AmpPeriodF;
%                 BioSoundUniqParam(ee_count,7) = BioSoundCall.AmpPeriodP;
%             end
% 
%             % Amplitude parameters calculated on the envelope
%              BioSoundUniqParam(ee_count,8) = BioSoundCall.rms;
%               BioSoundUniqParam(ee_count,9) = BioSoundCall.maxAmp;
% 
%             % Spectral parameters calculated on the spectrum
%             BioSoundUniqParam(ee_count,10) = BioSoundCall.stdspect;
%             BioSoundUniqParam(ee_count,11) = BioSoundCall.meanspect;
%             BioSoundUniqParam(ee_count,12) = BioSoundCall.skewspect;
%             BioSoundUniqParam(ee_count,13) = BioSoundCall.entropyspect;
%             BioSoundUniqParam(ee_count,14) = BioSoundCall.kurtosisspect;
%             BioSoundUniqParam(ee_count,15) = BioSoundCall.q1;
%             BioSoundUniqParam(ee_count,16) = BioSoundCall.q2;
%             BioSoundUniqParam(ee_count,17) = BioSoundCall.q3;
% 
%             % Fundamental parameters
%             if ~isempty(BioSoundCall.fund)
%                 BioSoundUniqParam(ee_count,18) = BioSoundCall.fund;
%             end
%             if ~isempty(BioSoundCall.cvfund)
%                 BioSoundUniqParam(ee_count,19) = BioSoundCall.cvfund;
%             end
%             if ~isempty(BioSoundCall.minfund)
%                 BioSoundUniqParam(ee_count,20) = BioSoundCall.minfund;
%             end
%             if ~isempty(BioSoundCall.maxfund)
%                 BioSoundUniqParam(ee_count,21) = BioSoundCall.maxfund;
%             end
%             if ~isempty(BioSoundCall.meansal)
%                 BioSoundUniqParam(ee_count,22) = BioSoundCall.meansal;
%             end
%             BioSoundUniqParam(ee_count,23) = CorrectAutoDetection01.(sprintf(AL_AutoId{ll_auto}))(ee);
        
          [BioSoundUniqParam(ee_count,1:16),AcounsticFeatureNames] = run_acoustic_features(Sound, FS_local, F_High, F_low, F_highSpec);
          BioSoundUniqParam(ee_count,17) = CorrectAutoDetection01.(sprintf(AL_AutoId{ll_auto}))(ee);
%         audiowrite(fullfile(Data_out, sprintf('Sound_%s_%d_%d_%d_%d.wav', (sprintf(AL_AutoId{ll_auto})), ee, OnInd, OffInd,CorrectAutoDetection01.(sprintf(AL_AutoId{ll_auto}))(ee) )),Sound, FS_local);
    end
end
% Turn back on warnings regarding Pyton to structure conversion
% warning('on', 'MATLAB:structOnObject')
%save(fullfile(Path2Data1, 'SoundEvent.mat'),'BioSoundUniqParam', 'BioSoundParamNames','AL_AutoId','AL_ManId', 'MissedAutoDetection','TotManCall','CorrectAutoDetection01','NLoggers','ManCallTranscTime_ms','ManCallMicSamp','ManCallLogSamp', 'SamplingFreq','ManCallMicFile','-append')
save(fullfile(Path2Results1, 'SoundEvent2.mat'),'BioSoundUniqParam', 'AcounsticFeatureNames','AL_AutoId','AL_ManId', 'MissedAutoDetection','TotManCall','CorrectAutoDetection01','NLoggers','ManCallTranscTime_ms','ManCallMicSamp','ManCallLogSamp', 'SamplingFreq','ManCallMicFile')
BioSoundParamNames = AcounsticFeatureNames;
%% Draw some scatters of the parameters
NParam = size(BioSoundUniqParam,2);
for pp=1:(NParam-1)
    figure();
    histogram(BioSoundUniqParam(~BioSoundUniqParam(:,17),pp), 'Normalization', 'probability')
    hold on
    histogram(BioSoundUniqParam(logical(BioSoundUniqParam(:,17)),pp), 'Normalization', 'probability')
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
UsefulParams = 1:16;
% Try a support vector machine classifier (linear) Binary SVM
SVMModel = fitcsvm(BioSoundUniqParam(:,UsefulParams),BioSoundUniqParam(:,17),'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto','Prior','Uniform');
% Cross-validate the SVM classifier. By default, the software uses 10-fold cross-validation.
CVSVMModel = crossval(SVMModel);
%Estimate the out-of-sample misclassification rate.
classLoss = kfoldLoss(CVSVMModel) % 9.8%-10.8% error in cross-validation


% Binary Kernel classification (non-linear)
CVMdl = fitckernel(BioSoundUniqParam(:,UsefulParams),BioSoundUniqParam(:,17),'CrossVal','on','Prior','Uniform')
%CVMdl is a ClassificationPartitionedKernel model. Because fitckernel implements 10-fold cross-validation, CVMdl contains 10 ClassificationKernel models that the software trains on training-fold (in-fold) observations.
%Estimate the cross-validated classification error.
kfoldLoss(CVMdl) % 50.4% error!


% Let's try to predict data using a Binary SVM
oosInds = unique(randi(ee_count,[round(ee_count/10) 1]));   % Out-of-sample indices
isInds = setdiff(1:ee_count, oosInds);   % In-sample indices
X_train = BioSoundUniqParam(isInds,UsefulParams);
Y_train = BioSoundUniqParam(isInds,17);
X_test = BioSoundUniqParam(oosInds,UsefulParams);
Y_test = BioSoundUniqParam(oosInds,17);
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
 % If threshold posterior probability set at 0.03, then false positive
 % brought down to 2-10% and % of missed vocalizations brought down to 1.7

 %% Save the SVM model for use/prediction with other recordings
 UsefulParams = 1:16;
% Try a support vector machine classifier (linear) Binary SVM
SVMModel = fitcsvm(BioSoundUniqParam(:,UsefulParams),BioSoundUniqParam(:,17),'Standardize',true,'KernelFunction','RBF',...
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
    BioSoundUniqParam(:,UsefulParams),BioSoundUniqParam(:,17)) 
save('/Users/elie/Documents/CODE/SoundAnalysisBats/SVMModelNoiseVoc.mat', 'CompactSVMModel')

%% Test a UMAP projection on this dataset
addpath /Users/elie/Documents/CODE/umap_1.4.1/umap
addpath /Users/elie/Documents/CODE/umap_1.4.1/util
javaaddpath('/Users/elie/Documents/CODE/umap_1.4.1/umap/umap.jar')
fprintf(1,'\n\n ***** Distance Euclidean *****\n')
[Reduction,UMAP,ClustID]= run_umap(BioSoundUniqParam(:,UsefulParams));
figure()
scatter(Reduction(:,1), Reduction(:,2),5,[BioSoundUniqParam(:,17) zeros(size(Reduction,1),2)],'filled')
title('Euclidean distance')

[Reduction,UMAP,ClustID]= run_umap(BioSoundUniqParam(:,UsefulParams),'metric','cosine');
figure()
scatter(Reduction(:,1), Reduction(:,2),5,[BioSoundUniqParam(:,17) zeros(size(Reduction,1),2)],'filled')
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
load(fullfile(Path2Results1, 'SoundEvent2.mat'),'BioSoundUniqParam', 'AcounsticFeatureNames','AL_AutoId', 'MissedAutoDetection','TotManCall','CorrectAutoDetection01','NLoggers')
addpath(genpath('/Users/elie/Documents/CODE/neurobat-callCutting'))
% WD = '/Users/elie/Documents/ManipBats/LMC/20190603';
% FS = 192000;
% findcalls_session(WD,FS,'fileType','wav')

FS = 192000;
findcalls_session(Path2Data2,FS,'fileType','wav','audio_file_filter','JuGr_190927_1014_RecOnly_mic1*','filter_raw_data',true)

% convert automatically detected Audio time to transceiver times 
[MicVoc_samp_idx,MicVoc_transcTime_ms] = mic2transc_time(Path2Data2);

% 
%% Let's loop in the dataset of manually extracted calls (Compare to ground truth) and see how many correct hits we get in the automatic detection

Delay = 100; %in ms error/delay between auto and man detection and Delay to add before each detected call in ms
MissedMicAutoDetection = cell(NLoggers,1);
TotManCall = 0;
MicCorrectAutoDetection01 = zeros(size(MicVoc_transcTime_ms,1),1);
[z,p,k] = butter(6,[100 90000]/(192000/2),'bandpass');
sos_raw_band = zp2sos(z,p,k);

for ll=1:NLoggers
    MissedMicAutoDetection{ll} = [];
    
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
            MissedMicAutoDetection{ll} = [MissedMicAutoDetection{ll} vv];
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
            hold off
            
            
            title(sprintf('%s  Voc %d/%d NOT DETECTED', AL_ManId{ll}, vv, size(ManCallTranscTime_ms{ll},1)))
            Player = audioplayer(Raw_ramp/std(Raw_ramp), FSRaw);
            play(Player)
%             pause()
            
            
        else
            IdxAll = union(Idx_OnsetAuto,Idx_OffsetAuto);
            MicCorrectAutoDetection01(IdxAll) = ones(size(IdxAll));
            fprintf('%d events started during that vocalization\n',length(IdxAll))
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
fprintf(1,'Missed call by the microphone auto detection: %d/%d or %.1f%%\n', sum(cellfun('length',MissedMicAutoDetection)),TotManCall,sum(cellfun('length',MissedMicAutoDetection))/TotManCall*100);
fprintf(1,'Number of detected events from the microphone %d, proportion of true calls %d/%d or %.1f%%\n',length(MicCorrectAutoDetection01),sum(MicCorrectAutoDetection01),length(MicCorrectAutoDetection01),sum(MicCorrectAutoDetection01)/length(MicCorrectAutoDetection01)*100);


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
    
