function [Voc_filename,Voc_samp_idx,Voc_transc_time] = voc_localize_using_piezo(Logger_dir,RawWav_dir, Date, ExpStartTime, varargin)
%% VOC_LOCALIZE_PIEZO a function to detect vocalizations based on threshold crossing in piezo recordings
%% and retrieve the position of these extracts in continuous microphone recordings

% Inputs
% Logger_dir is the folder containing the logger folders containing the
% extracted/CS0.mat files

% RawWav_dir is the folder containing the continuous recordings and the
% file *TTLPulseTimes.mat generated by align_soundmexAudio_2_logger.m
% or align_avisoft_2_logger.m in case calculating transceiver time is requested

% 'TransceiverTime' (optional input): set by default to 1 to calculate the
% onset anf offset of sound extracts in transceiver time

% Ouputs
% Voc_filename is the file list of automatic extracts 

% Voc_samp_idx is a 2 column vector that gives the onset and offset indices
% of each extract in the original recordings, same number of lines as
% Voc_filename

% Voc_transc_time is a 2 column vector that gives the expected onset and offest
% times of each extract in the piezo logger recordings in transceiver time,
% in ms
% same number of lines as Voc_filename.

% Hard coded input for finding the microphone envelope noise threshold
Dur_RMS = 0.5; % duration of the silence sample in min for the calculation of average running RMS
Fhigh_power = 20; %Hz
Fs_env = 1000; %Hz Sample frequency of the enveloppe
MicThreshNoise = 15*10^-3;
Merge_thresh = 500; % merge thershold in ms for grouping sounds in sequences if they are withinh that delay
if ~exist(fullfile(RawWav_dir, 'Detected_calls'),'dir')
    mkdir(fullfile(RawWav_dir, 'Detected_calls'))
end
    
    

%% Load data and initialize output variables
% Load the pulse times and samples
TTL_dir = dir(fullfile(RawWav_dir,sprintf( '%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
TTL = load(fullfile(TTL_dir.folder, TTL_dir.name));


%% Run the automatic detection
AllLoggers = dir(fullfile(Logger_dir, '*ogger*'));
DirFlags = [AllLoggers.isdir];
% Extract only those that are directories.
AllLoggers = AllLoggers(DirFlags);
% Loop through loggers and detect vocalizations based on thereshold
% crossings on the envelope
NL = length(AllLoggers);
Nevents = nan(NL,1);
for ll=1:NL
    Data_directory = fullfile(AllLoggers(ll).folder,AllLoggers(ll).name, 'extracted_data');
    [SoundEvent_LoggerSamp.(sprintf('L%s',AllLoggers(ll).name(2:end))),SoundEvent_TranscTime_ms.(sprintf('L%s',AllLoggers(ll).name(2:end))),LoggerEnvelopeAll.(sprintf('L%s',AllLoggers(ll).name(2:end))),SoundEvent_EnvSamp.(sprintf('L%s',AllLoggers(ll).name(2:end))),~] = piezo_find_calls_logger(Data_directory);
    Nevents = size(SoundEvent_LoggerSamp.(sprintf('L%s',AllLoggers(ll).name(2:end))),1);
end

% Sort vocalization from noise using an SVM approach on acoustic parameters
Buffer = 10;% time in ms to add before after each sound element such atht it's longer than the 23ms required for biosound to calculate fundamental and saliency parameters
F_high = 5000; % frequency low pass on the logger detected sounds for calculating temporal acoustic features and Saliency
F_low = 100;% frequency high pass n the logger detected sounds
F_highSpec = 15000;% frequency low pass on the logger detected sounds for calculating spectral parameters
TotEvents = sum(Nevents);
AcousticParams = nan(TotEvents,11);
LoggerID_unmerged = cell(TotEvents,1);
FS_logger_voc_unmerged = nan(TotEvents,1);
Voc_loggerSamp_Idx_unmerged = nan(TotEvents,2);
Voc_transc_time_unmerged = nan(TotEvents,2);
ALField_Id = fieldnames(SoundEvent_TranscTime_ms); % Names of the audioLoggers
ee_count=0;
for ll=1:NL
    fprintf(1, '*** %s %d/%d ****\n',ALField_Id{ll}, ll, NL)
    % Load the raw signal
    Data_directory = fullfile(AllLoggers(ll).folder,AllLoggers(ll).file, 'extracted_data');
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
    for ee=1:Nevents(ll)
        ee_count = ee_count+1;
        if rem(ee,100)==0
            fprintf(1, 'Event %d/%d\n', ee,Nevents(ll))
        end
        % Onset and offset of detected sound extract
        OnInd = SoundEvent_LoggerSamp.(sprintf(ALField_Id{ll}))(ee,1);
        OffInd = SoundEvent_LoggerSamp.(sprintf(ALField_Id{ll}))(ee,2);
        % find the sampling Frequency
        FileIdx = find((Indices_of_first_and_last_samples(:,1)<OnInd) .* (Indices_of_first_and_last_samples(:,2)>OffInd));
        if isempty(FileIdx) || (length(FileIdx)~=1) || FileIdx>length(Estimated_channelFS_Transceiver)
            FS_logger_voc_unmerged(ee) = round(nanmean(Estimated_channelFS_Transceiver));
        else
            FS_logger_voc_unmerged(ee) = round(Estimated_channelFS_Transceiver(FileIdx));
        end
        
        % extract the sound with Buffer ms before after the sound
        OnIndBuff = OnInd - round(FS_logger_voc_unmerged(ee)*Buffer*10^-3);
        OffIndBuff = OffInd + round(FS_logger_voc_unmerged(ee)*Buffer*10^-3);
        Sound = Centered_piezo_signal(OnIndBuff : OffIndBuff);
        Sound = Sound - mean(Sound);
        
       AcousticParams(ee_count,:) = run_acoustic_features(Sound, FS_logger_voc_unmerged(ee), F_high, F_low, F_highSpec);
       LoggerID_unmerged{ee_count} = ALField_Id{ll};
       Voc_loggerSamp_Idx_unmerged(ee_count,:) = [OnInd OffInd];
       Voc_transc_time_unmerged(ee_count,:) = SoundEvent_TranscTime_ms.(sprintf(ALField_Id{ll}))(ee,:);
    end
end
% Load the SVM compact model, predict labels according to the model and
% eliminate noise
load('SVMCompact.mat','CompactSVMModel')
[~,PostProbs] = predict(CompactSVMModel,AcousticParams);
Labels = PostProbs(:,2)>=10^-3; % Here we are being very conservative and keep events that even have a slight chance of being vocalizations
% Select presume vocalizations and order them in time
LoggerID_unmerged = LoggerID_unmerged(logical(Labels));
Voc_loggerSamp_Idx_unmerged = Voc_loggerSamp_Idx_unmerged(logical(Labels),:);
Voc_transc_time_unmerged = Voc_transc_time_unmerged(logical(Labels),:);
FS_logger_voc_unmerged = FS_logger_voc_unmerged(logical(Labels),:);

[Voc_transc_time_unmerged, OrdInd] = sort(Voc_transc_time_unmerged(:,1));
LoggerID_unmerged = LoggerID_unmerged(OrdInd);
Voc_loggerSamp_Idx_unmerged = Voc_loggerSamp_Idx_unmerged(OrdInd,:);
FS_logger_voc_unmerged = FS_logger_voc_unmerged(OrdInd,:);
TotEvents = sum(Labels);

%% Merge vocalizations into sequences if they are less than Merge_thresh appart to avoid repetition in Who Calls and add Merge_Thresh before/after
Events2Merge = [0; (Voc_transc_time_unmerged(2:end,1)-Voc_transc_time_unmerged(1:end-1,2))<= Merge_thresh];
FirstEvents2Merge = find(diff([Events2Merge; 0])==1); % onset of each sequence of events that should be merged
LastEvents2Merge = find(diff([Events2Merge; 0])==-1);% offset of each sequence of events that should be merged
Events2keep = strfind([Events2Merge' 0],[0 0]); % events that should be kept as they are
if length(FirstEvents2Merge)~=length(LastEvents2Merge)
    warning('Problem in the detection of sequences of sound events to merge')
    keyboard
end
Voc_transc_time = [Voc_transc_time_unmerged(Events2keep,:) ; [Voc_transc_time_unmerged(FirstEvents2Merge,1) Voc_transc_time_unmerged(LastEvents2Merge,2)]];
Voc_loggerSamp_Idx = [Voc_loggerSamp_Idx_unmerged(Events2keep,:) ; [Voc_loggerSamp_Idx_unmerged(FirstEvents2Merge,1) Voc_loggerSamp_Idx_unmerged(LastEvents2Merge,2)]];
LoggerID = [LoggerID_unmerged(Events2keep,:) ; [LoggerID_unmerged(FirstEvents2Merge,1) LoggerID_unmerged(LastEvents2Merge,2)]];
FS_logger_voc = [FS_logger_voc_unmerged(Events2keep,:) ; [FS_logger_voc_unmerged(FirstEvents2Merge,1) FS_logger_voc_unmerged(LastEvents2Merge,2)]];

% reorder in time
[Voc_transc_time, OrdInd] = sort(Voc_transc_time(:,1));
LoggerID = LoggerID(OrdInd);
Voc_loggerSamp_Idx = Voc_loggerSamp_Idx(OrdInd,:);
FS_logger_voc =FS_logger_voc(OrdInd,:);

% Add MergeThresh before and after each
Voc_transc_time(:,1) = Voc_transc_time(:,1) - Merge_thresh;
Voc_transc_time(:,2) = Voc_transc_time(:,2) + Merge_thresh;
Voc_loggerSamp_Idx(:,1) = Voc_loggerSamp_Idx(:,1) - Merge_thresh*10^-3.*FS_logger_voc;
Voc_loggerSamp_Idx(:,2) = Voc_loggerSamp_Idx(:,2) + Merge_thresh*10^-3.*FS_logger_voc;

%% Retrieve the Microphone file that contains the data for each detected sequence of vocalization
MeanStdAmpRawFile = nan(100,2);
MeanStdAmpRawExtract = nan(TotEvents,2);
Voc_samp_idx = nan(TotEvents,2);
Voc_filename = cell(TotEvents,1);
% Construct the frequency bandpass filter for calculating the noise rate on
% the microphone data
WavFileStruc_local = dir(fullfile(RawWav_dir, sprintf('*_%s_%s*mic*.wav',Date, ExpStartTime)));
Raw_filename = fullfile(WavFileStruc_local(1).folder, WavFileStruc_local(1).name);
Subj = WavFileStruc_local(1).name(1:4);
Info = audioinfo(Raw_filename);
FS = Info.SampleRate;
[z,p,k] = butter(6,[1000 90000]/(FS/2),'bandpass');
sos_raw_band = zp2sos(z,p,k);

for ee=1:TotEvents
    % Find the microphone file
    OnFile = find(TTL.Pulse_TimeStamp_Transc<Voc_transc_time_unmerged(ee,1),1,'Last');
    OffFile = find(TTL.Pulse_TimeStamp_Transc>Voc_transc_time_unmerged(ee,2),1,'First')-1;
    FileIdx = min(TTL.File_number(OnFile), TTL.File_number(OffFile));
    
    if isnan(MeanStdAmpRawFile(FileIdx,1)) % calculate the amplitude threshold for that file
        % load the raw file
        WavFileStruc_local = dir(fullfile(RawWav_dir, sprintf('*_%s_%s*mic*_%d.wav',Date, ExpStartTime, FileIdx)));
        Raw_filename = fullfile(WavFileStruc_local.folder, WavFileStruc_local.name);
        [Raw_10minwav, FS] = audioread(Raw_filename);
        
        % Calculate the amplitude threshold as the average amplitude on the
        % first 30 seconds of that 10 min recording file from which that file
        % come from
        % Get the average running rms in a Dur_RMS min extract in the middle of
        % the recording
        fprintf(1, 'Calculating average RMS values on a %.1f min sample of silence\n',Dur_RMS);
        SampleDur = round(Dur_RMS*60*FS);
        StartSamp = round(length(Raw_10minwav)/2);
        fprintf(1,'Calculating the amplitude threshold for file %d  ',FileIdx)
        BadSection = 1;
        while BadSection
            Filt_RawVoc = filtfilt(sos_raw_band,1,Raw_10minwav(StartSamp + (1:round(SampleDur))));
            Amp_env_Mic = running_rms(Filt_RawVoc, FS, Fhigh_power, Fs_env);
            if any(Amp_env_Mic>MicThreshNoise) % there is most likely a vocalization in this sequence look somewhere else!
                StartSamp = StartSamp + SampleDur +1;
            else
                BadSection = 0;
            end
        end
        MeanStdAmpRawFile(FileIdx,1) = mean(Amp_env_Mic);
        MeanStdAmpRawFile(FileIdx,2) = std(Amp_env_Mic);
        fprintf('-> Done\n')
    end
    MeanStdAmpRawExtract(ee,1)= MeanStdAmpRawFile(FileIdx,1);
    MeanStdAmpRawExtract(ee,2)= MeanStdAmpRawFile(FileIdx,2);
    
    % Calculate the Samples of the extract in the microphone recording
    TTL_idx = find(unique(TTL.File_number) == FileIdx);
    Voc_transc_time_zs = (Voc_transc_time - TTL.Mean_std_Pulse_TimeStamp_Transc(TTL_idx,1))/TTL.Mean_std_Pulse_TimeStamp_Transc(TTL_idx,2);
    Voc_samp_idx(ee,:) = TTL.Mean_std_Pulse_samp_audio(TTL_idx,2) .* polyval(TTL.Slope_and_intercept_transc2audiosamp{TTL_idx}, Voc_transc_time_zs,[],TTL.Mean_std_x_transc2audiosamp{TTL_idx}) + TTL.Mean_std_Pulse_samp_audio(TTL_idx,1);
    
    % Extract the wave from the microphone recording
    Raw_wave = Raw_10minwav(Voc_samp_idx(ee,1) : Voc_samp_idx(ee,2));
    
    % Save the sound as a wav file 
    Voc_filename{ee} = fullfile(RawWav_dir, 'Detected_calls',sprintf('%s_%s_%s_voc_%d_%d.wav',Subj,Date,ExpStartTime, FileIdx, Voc_samp_idx(ee,1)));
    audiowrite(Voc_filename{ee} , Raw_wave, FS)
end



%% save the calculation results
save(fullfile(RawWav_dir, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)), 'Voc_filename','Voc_samp_idx','Voc_transc_time','MeanStdAmpRawExtract','Voc_loggerSamp_Idx','LoggerID','FS_logger_voc','LoggerID_unmerged','FS_logger_voc_unmerged','Voc_loggerSamp_Idx_unmerged','Voc_transc_time_unmerged')

end
