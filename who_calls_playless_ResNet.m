function [IndVocStartRaw_merged, IndVocStopRaw_merged, IndVocStartPiezo_merge_local, IndVocStopPiezo_merge_local] = who_calls_playless_ResNet(Raw_dir, Loggers_dir, Date, ExpStartTime, MergeThresh, UseOld,CheckMicChannel, varargin)

% optional parameter: Factor_RMS_Mic, Factor by which the RMS of the
% band-pass filtered baseline signal is multiplied to obtained the
% threshold of vocalization detection on Microphone
RMS_Fig=0;
VolDenominatorLogger=5;
VolFactorMic=0.5;
pnames = {'Factor_RMS_Mic','Working_dir','Force_Save_onoffsets_mic','SaveFileType', 'Sorter'};
dflts  = {3,Loggers_dir,0,'pdf','JulieE'};
[Factor_RMS_Mic,Working_dir,Force_Save_onoffsets_mic,SaveFileType,  Sorter] = internal.stats.parseArgs(pnames,dflts,varargin{:});
% UsePrevious is a boolean to indicate if you sould only check the
% vocalizations that are new after the MergePatch

if nargin<7
    CheckMicChannel = 0;
end

if nargin<6
    UseOld = 0;
end

if nargin<5
    MergeThresh = 500; % any 2 detected call spaced by less than MergeThresh ms are merged
end
PPT = [0.96 0.9]; % Post probability threshold to consider a sequence containing a vocalizations or for sure being noise

Working_dir_read = fullfile(Working_dir, 'read');
Working_dir_write = fullfile(Working_dir, 'write');

if ~strcmp(Loggers_dir,Working_dir) && (~exist(Working_dir,'dir') || ~exist(Working_dir_read,'dir') || ~exist(Working_dir_write,'dir'))
    mkdir(Working_dir)
    mkdir(Working_dir_read)
    mkdir(Working_dir_write)
elseif strcmp(Loggers_dir,Working_dir)
    Working_dir_read = Loggers_dir;
    Working_dir_write = Loggers_dir;
end
DataFiles = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData*.mat', Date, ExpStartTime)));

ResNetDataFile = dir(fullfile(Loggers_dir, sprintf('*%s.csv', Date)));
ResNetData = readtable(fullfile(ResNetDataFile.folder, ResNetDataFile.name), 'Delimiter',',');

if isempty(DataFiles)
    warning('Vocalization data were not extracted by get_logger_data_voc.m')
else
    % select the correct files
    Gdf = zeros(length(DataFiles),1);
    for df=1:length(DataFiles)
        if length(strfind(DataFiles(df).name, '_'))==2
            Gdf(df)=1;
        end
    end
    DataFiles = DataFiles(logical(Gdf));
    % gather File indices to reorder them
    IndDataFiles = nan(length(DataFiles),1);
    for nfile = 1:length(DataFiles)
        IndData = strfind(DataFiles(nfile).name, 'Data') + length('Data');
        IndDot = strfind(DataFiles(nfile).name, '.') - 1;
        IndDataFiles(nfile) = str2double(DataFiles(nfile).name(IndData:IndDot));
    end
    [~,AscendOrd] = sort(IndDataFiles);
    InDataFilesOrdered = IndDataFiles(AscendOrd);
    DataFiles = DataFiles(AscendOrd);
    load(fullfile(Raw_dir, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)), 'MeanStdAmpRawExtract','Voc_filename')
    Nvoc_all = length(Voc_filename);
    DataFile = fullfile(DataFiles(1).folder, DataFiles(1).name);
    load(DataFile, 'VocMaxNum')
    if ~exist('VocMaxNum','var')
        VocMaxNum=1000;
    end
    if Nvoc_all>VocMaxNum
        Nvocs = [0 VocMaxNum:VocMaxNum:Nvoc_all (floor(Nvoc_all/VocMaxNum)*VocMaxNum+rem(Nvoc_all,VocMaxNum))];
    else
        Nvocs = [0 Nvoc_all];
    end
    for df=1:length(DataFiles)
        fprintf(1,'Set %d/%d\n', df, length(DataFiles))
        
        %% Identify sound elements in each vocalization extract and decide of the vocalizer
        % Output corresponds to onset and offset indices of each vocal element in
        % sound section for each vocalizing logger. A cell array of length the
        % number of sound sections identified by voc_localize or
        % voc_localize_operant. For IndVocStartRaw_merged each cell contains the onset/offset samples of each vocal element in the sound section.
        % For IndVocStartPiezo_merged each cell is a cell array of size the number of
        % audio-loggers where each cell contains the onset or offset sample of
        % individual vocal element emitted by the respective audio-loggers. Note
        % that for both these variables, vocal elements that are within MergeThresh
        % ms were merged in a single vocal element.
        %
        % Run a time window of duration 2 ms to identify who is vocalizing
        % create 2 signals for the logger, one high pass filtered above 5kHz
        % and the other one low pass filtered at 5kHz. The vocalizer will have
        % the highest energy of all vocalizers and will have more energy in the
        % lower compare to higher filtered signal
        
        Nvoc = Nvocs(df+1) - Nvocs(df);
        DataFile = fullfile(DataFiles(df).folder, DataFiles(df).name);
        
        PreviousFile = fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime, df,MergeThresh));
%         if ~isfile(PreviousFile)
%             PreviousFile = fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime,MergeThresh));
%         end
        if ~isempty(dir(PreviousFile)) && UseOld
%             load(PreviousFile, 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStart_all', 'IndVocStop_all','IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','IndNoiseStart_all','IndNoiseStop_all', 'IndNoiseStartRaw', 'IndNoiseStopRaw', 'IndNoiseStartPiezo', 'IndNoiseStopPiezo','RMSRatio_all','RMSDiff_all','vv','MicError','PiezoError','MicErrorType','PiezoErrorType');
            load(PreviousFile, 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','vv','List2WorkOnResNet','MicError','PiezoError','MicErrorType','PiezoErrorType', 'SorterName');
            if ~exist('vv','var') % There is no previous data but just data regarding piezo numbers and bats_ID
                fprintf(1, 'No Previous data...ok?')
                keyboard
                vv=1;
                clear List2WorkOnResNet
            end

        else
            vv=1;
            clear List2WorkOnResNet
        end
        
        % Check if that set of vocalizations was alread fully completed
        if Nvoc==1
            CheckDone = input('There is only one vocalization in this set, indicate if it was already processed or not as it cannot be determined programmaticaly (Redo:1 Done:0)');
            if ~CheckDone
                clear List2WorkOnResNet vv
                continue
            end
        elseif vv==Nvoc
            % All done
            clear List2WorkOnResNet vv
            continue
        end
        
        if ~strcmp(Working_dir_write,Loggers_dir) && ~isfile(fullfile(Working_dir_read,DataFiles(df).name))
            fprintf(1,'Bringing data locally from the server\n')
            [s,~]=copyfile(DataFile, Working_dir_read, 'f');
            if ~s
                fprintf(1,'File transfer did not occur correctly\n')
                keyboard
            else
                DataFile = fullfile(Working_dir_read,DataFiles(df).name);
            end
        end
        load(DataFile,'Raw_wave')
        if Nvoc ~= length(Raw_wave)
            warning('Looks like there might be an issue there!! Check variables!!')
            keyboard
        end
        minvv = 1;
%         load(DataFile,'Piezo_wave', 'AudioLogs',   'Piezo_FS',  'FS', 'DiffRMS', 'RMSLow','VocFilename','Voc_transc_time_refined');
        load(DataFile,'Piezo_wave', 'AudioLogs',   'Piezo_FS',  'FS', 'DiffRMS', 'RMSLow','VocFilename');
        Fns_AL = fieldnames(Piezo_wave);
        
        % parameters
        Consecutive_binsMic = 10; % Number of consecutive bins of the envelope difference between highpass and low pass logger signal that has to be higher than threshold to be considered as a vocalization
        Consecutive_binsPiezo = 15; % Number of consecutive bins of the envelope difference between highpass and low pass logger signal that has to be higher than threshold to be considered as a vocalization
        Factor_RMS_low = 1.5.*ones(length(AudioLogs),1); % Factor by which the RMS of the low-pass filtered baseline signal is multiplied to obtained the threshold of vocalization detection on piezos I sometimes have to manually increase that one for some loggers Logger 39 on 200131 needs ot have a threshold of 2.75...
        % Factor_AmpRatio = 1.5; % Factor by which the ratio of amplitude between low and high  pass filtered baseline signals is multiplied to obtain the threshold on calling vs hearing (when the bats call there is more energy in the lower frequency band than higher frequency band of the piezo) % used to be 3
        Factor_AmpDiff = 50; % Factor by which the ratio of amplitude between low and high  pass filtered baseline signals is multiplied to obtain the threshold on calling vs hearing (when the bats call there is more energy in the lower frequency band than higher frequency band of the piezo) % used to be 3
        DB_noise = 60; % Noise threshold for the spectrogram colormap
        FHigh_spec = 90000; % Max frequency (Hz) for the raw data spectrogram
        FHigh_spec_Logger = 10000; % Max frequency (Hz) for the raw data spectrogram
        BandPassFilter = [1000 5000 9900]; % Frequency bands chosen for digital signal processing
        Fhigh_power =50; % Frequency upper bound for calculating the envelope (time running RMS)
        Fs_env = 1000; % Sample frequency of the enveloppe
        % Buffer=100; % Maximum lag calculated for the cross correlation in ms, not used in that code for identification puposes but still calculated as of now

        if df==1 || ~exist('sos_raw_band', 'var')
            % design filters of raw ambient recording, bandpass and low pass which was
            % used for the cross correlation
            [z,p,k] = butter(6,[BandPassFilter(1) 90000]/(FS/2),'bandpass');
            sos_raw_band = zp2sos(z,p,k);
            % [z,p,k] = butter(6,BandPassFilter(1:2)/(FS/2),'bandpass');
            % sos_raw_low = zp2sos(z,p,k);
            [z,p,k] = butter(6,[100 20000]/(FS/2),'bandpass');
            sos_raw_band_listen = zp2sos(z,p,k);
        end
        
        
        if df==1 || ~exist('sos_raw_band_listen', 'var')
            % design filters of raw ambient recording, bandpass, for
            % listening 
            [z,p,k] = butter(6,[100 20000]/(FS/2),'bandpass');
            sos_raw_band_listen = zp2sos(z,p,k);
        end
        
        % Initialize variables
        
        
        if vv==1 % We need to initialize variables!
            %         Amp_env_LowPassLogVoc = cell(1,Nvoc);
            %         Amp_env_HighPassLogVoc = cell(1,Nvoc);
            %         Amp_env_Mic = cell(1,Nvoc);
            %         LowPassLogVoc = cell(1,Nvoc);
            % LowPassMicVoc = cell(1,Nvoc);
            
            IndVocStart_all = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers+microphone and for each logger the index onset of when the animal start vocalizing in the piezo recording before merge in envelope unit (FS_env)
            IndVocStop_all = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal start vocalizing in the piezo recording before merge in envelope unit (FS_env)
            %         IndNoiseStart_all = cell(1,Nvoc);
            %         IndNoiseStop_all = cell(1,Nvoc);
            IndVocStartRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc)
            % a cell array of the size the number of loggers + 1 in case only one bat without a logger
            % or +2 incase no identification possible but you want to keep onset/offset of each voc and
            % for each logger the index onset of when the animal start vocalizing in the raw recording before merge
            IndVocStartPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording before merge
            IndVocStopRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording before merge
            IndVocStopPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording before merge
            %         IndNoiseStartRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the raw recording before merge
            %         IndNoiseStartPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording before merge
            %         IndNoiseStopRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording before merge
            %         IndNoiseStopPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording before merge
            IndVocStartRaw_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the raw recording
            IndVocStopRaw_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording
            IndVocStartPiezo_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording
            IndVocStopPiezo_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording
            RMSRatio_all = cell(1,Nvoc);
            RMSDiff_all = cell(1,Nvoc);
            MicError = [0 0];% first element = number of corrections. second = number of detection (question)
            MicErrorType = [0 0];% first element false negative (detected as noise or already detected when it is a new call), second element false positive (vice versa)
            PiezoError = [0 0];% first element = number of corrections. second = number of detection (question)
            PiezoErrorType = [0 0]; % first element false negative (detected as noise or hearing when it is a new call), second element false positive (vice versa)
            SorterName = cell(1,Nvoc); % Contains for each sequence of vocalizations (Nvoc) the name of the sorter
        end
        
        
        
        %% Loop through vocalizations sequences and only keep the ones selected by resnet
        if ~exist('List2WorkOnResNet', 'var')
            List2WorkOn = vv:Nvoc;
            fprintf(1,'Only Select the sequences with a ResNet posterior probability higher than %.2f \n', PPT(1))
            List2WorkOnResNet = List2WorkOn;        
            for vvi=1:length(List2WorkOn)
                vv=List2WorkOn(vvi);
                RowsData = logical(contains(ResNetData.Var1(:), sprintf('15120%s%s_%d_Logger', Date, DataFiles(df).name(1:11), InDataFilesOrdered(df))) .*contains(ResNetData.Var1(:), sprintf('_%d.mat', vv-1)));
                pVoc = ResNetData.Var3(RowsData); % all probability of a vocalization across all 100ms windows of all loggers through that sequence
                if ~any(pVoc>PPT(1))
                    List2WorkOnResNet(List2WorkOnResNet==vv)=[];
                    IndVocStart_all{vv} = [];
                    IndVocStop_all{vv} = [];
                    %                 IndNoiseStart_all{vv} = [];
                    %                 IndNoiseStop_all{vv} = [];
                    IndVocStartRaw_merged{vv} = [];
                    IndVocStopRaw_merged{vv} = [];
                    IndVocStartPiezo_merged{vv} = [];
                    IndVocStopPiezo_merged{vv} = [];
                    SorterName{vv} = Sorter;
                    RMSRatio_all{vv} = [];
                    RMSDiff_all{vv} = [];
                    SorterName{vv} = 'ResNet';
                end
            end
            fprintf(1,'%d valid sequences found, saving data...\n', length(List2WorkOnResNet))
            if ~isempty(dir(PreviousFile))
                save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','List2WorkOnResNet','SorterName', '-append');
            else
                %                     save(fullfile(Working_dir, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','IndNoiseStart_all','IndNoiseStop_all', 'IndNoiseStartRaw', 'IndNoiseStopRaw', 'IndNoiseStartPiezo', 'IndNoiseStopPiezo','RMSRatio_all','RMSDiff_all','vv','MicError','PiezoError','MicErrorType','PiezoErrorType');
                save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','MicError','PiezoError','MicErrorType','PiezoErrorType', 'List2WorkOnResNet','SorterName');
            end
            List2WorkOn=List2WorkOnResNet;
        else
            List2WorkOn=List2WorkOnResNet(find(List2WorkOnResNet==vv):end);
        end
        if isempty(List2WorkOn)
            fprintf(1,'Nothing to change for set %d\n',df)
            clear List2WorkOnResNet vv
            continue
        end
                
                
        %% Loop through vocalizations sequences and calculate amplitude envelopes
        for vvi=1:length(List2WorkOn)
            vv=List2WorkOn(vvi);
            fprintf(1,'\n\n\n\nVoc sequence %d/%d Set %d/%d\n',vv,Nvoc, df, length(DataFiles));
            
            %% First calculate the time varying RMS of the ambient microphone
            Amp_env_LowPassLogVoc = cell(length(AudioLogs),1);
            Amp_env_HighPassLogVoc = cell(length(AudioLogs),1);
            LowPassLogVoc = cell(length(AudioLogs),1);
            % Patch for previous error in the code
            Raw_wave_nn = Raw_wave{vv - (minvv -1)};
            if length(Raw_wave_nn)<(128*FS/1000)
                warning('This extract is unexpectedly short..../n')
                keyboard
                continue
            end
              
            % retrieving file name index of the microphone
            Voc_i_start = Nvocs(df)+1;
            vv_in = vv + Voc_i_start-1;
            if ~strcmp(Voc_filename{vv_in}, VocFilename{vv})
                warning('Issues with Mic file name\n')
                keyboard
            end
            fprintf(1, 'Microphone File: %s\n', Voc_filename{vv_in})
            if Nvoc>100
                warning('Probably wrong audio file name, the code is not updated for older version of previous extraction\n')
            end
            
            if isempty(Raw_wave_nn)
                warning('Raw_wave should not be empty now!!')
                keyboard
                SaveRawWave = 1;
                [Raw_wave{vv}, FS] = audioread(VocFilename{vv});
            else
                SaveRawWave = 0;
            end
            
%             % Check that correct microphone file was saved (Trying to
%             % detect/fix bug from voc_localize_using_piezo
%             TTL_dir = dir(fullfile(Raw_dir,sprintf( '%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
%             TTL = load(fullfile(TTL_dir.folder, TTL_dir.name));
%             FileNum_u = unique(TTL.File_number);
%             OnOffTranscTime_ms = Voc_transc_time_refined(vv,:);
%             if ~any(isnan(OnOffTranscTime_ms))
%                 FileNumIdx = find(TTL.Pulse_TimeStamp_Transc<OnOffTranscTime_ms(1,1),1,'Last');
%                 if isempty(FileNumIdx)
%                     FileNumIdx = find(TTL.Pulse_TimeStamp_Transc>OnOffTranscTime_ms(1,1),1,'First');
%                 end
%                 MicVoc_File = TTL.File_number(FileNumIdx);
%                 IndFileNum = find(FileNum_u == MicVoc_File);
%                 TranscTime_zs = (OnOffTranscTime_ms - TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,1))/TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,2);
%                 MicVoc_samp_idx =round(TTL.Mean_std_Pulse_samp_audio(IndFileNum,2) .* polyval(TTL.Slope_and_intercept_transc2audiosamp{IndFileNum},TranscTime_zs,[],TTL.Mean_std_x_transc2audiosamp{IndFileNum}) + TTL.Mean_std_Pulse_samp_audio(IndFileNum,1));
%                 WavFileStruc_local = dir(fullfile(Raw_dir, sprintf('*_%s_%s*mic*_%d.wav',Date, ExpStartTime, MicVoc_File)));
%                 Raw_filename = fullfile(WavFileStruc_local.folder, WavFileStruc_local.name);
%                 [Raw_10minwav2, FS2] = audioread(Raw_filename);
%                 if MicVoc_samp_idx(1)>length(Raw_10minwav2) % This vocalization occured in the next file
%                     MicVoc_File = MicVoc_File+1;
%                     IndFileNum = find(FileNum_u == MicVoc_File);
%                     TranscTime_zs = (OnOffTranscTime_ms - TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,1))/TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,2);
%                     MicVoc_samp_idx =round(TTL.Mean_std_Pulse_samp_audio(IndFileNum,2) .* polyval(TTL.Slope_and_intercept_transc2audiosamp{IndFileNum},TranscTime_zs,[],TTL.Mean_std_x_transc2audiosamp{IndFileNum}) + TTL.Mean_std_Pulse_samp_audio(IndFileNum,1));
%                     WavFileStruc_local = dir(fullfile(Raw_dir, sprintf('*_%s_%s*mic*_%d.wav',Date, ExpStartTime, MicVoc_File)));
%                     Raw_filename = fullfile(WavFileStruc_local.folder, WavFileStruc_local.name);
%                     [Raw_10minwav2, FS2] = audioread(Raw_filename);
%                 end
%                 Raw_wave_ex = Raw_10minwav2(MicVoc_samp_idx(1) : min(MicVoc_samp_idx(2),length(Raw_10minwav2)));
%                 if length(Raw_wave_ex)<length(Raw_wave_nn)
%                     Corr(1) = corr(Raw_wave_ex,Raw_wave_nn(1:length(Raw_wave_ex)));
%                     Corr(2) = corr(Raw_wave_ex,Raw_wave_nn(end-length(Raw_wave_ex)+1:end));
%                 elseif length(Raw_wave_ex)>length(Raw_wave_nn)
%                     Corr(1) = corr(Raw_wave_nn,Raw_wave_ex(1:length(Raw_wave_nn)));
%                     Corr(2) = corr(Raw_wave_nn,Raw_wave_ex(end-length(Raw_wave_nn)+1:end));
%                 elseif length(Raw_wave_ex)==length(Raw_wave_nn)
%                     Corr = corr(Raw_wave_ex,Raw_wave_nn);
%                 end
%                 if all(Corr<0.99)
%                     warning('Error in the microphone file that was previosuly saved, fixing the issue now!\n')
%     %                 keyboard
%                     SaveRawWave = 1;
%                     Raw_wave_nn = Raw_wave_ex;
%                     Raw_wave{vv} = Raw_wave_ex;
%                     TrueVocName = fullfile(Raw_dir, 'Detected_calls',sprintf('%s_%s_%s_voc_%d_%d.wav',WavFileStruc_local.name(1:4),Date,ExpStartTime, MicVoc_File, MicVoc_samp_idx(1)));
%                     if ~strcmp(VocFilename{vv}, TrueVocName)
%                         warning('Filename was also wrong correcting %s -> %s\n',VocFilename{vv},TrueVocName)
%                         VocFilename{vv}= TrueVocName;
%                         Voc_filename{vv_in} = TrueVocName;
%                         SaveRawWaveName = 1;
%                     end
%                     audiowrite(VocFilename{vv} , Raw_wave{vv}, FS2);
%                 else
                    SaveRawWaveName = 0;
%                 end
%             end   

            % bandpass filter the ambient mic recording
            Filt_RawVoc = filtfilt(sos_raw_band,1,Raw_wave_nn);
            Amp_env_Mic = running_rms(Filt_RawVoc, FS, Fhigh_power, Fs_env);
            % Low passfilter the ambient mic recording (used for cross-correlation
            % but not done anymore
            %     LowPassMicVoc{vv} = filtfilt(sos_raw_low,1,Raw_wave{vv});
            % Plot the spectrogram of the ambient microphone
            F1=figure(1);
            clf(F1)
            ColorCode = [get(groot,'DefaultAxesColorOrder');0 0 0; 0 1 1; 1 1 0];
            subplot(length(AudioLogs)+2,1,1)
            [Raw_Spec.to, Raw_Spec.fo, Raw_Spec.logB] = spec_only_bats(Filt_RawVoc, FS, DB_noise, FHigh_spec);
            hold on
            yyaxis right
            plot((1:length(Amp_env_Mic))/Fs_env*1000, Amp_env_Mic, 'r-', 'LineWidth',2)
            ylabel(sprintf('Amp\nMic'))
            title(sprintf('Voc %d/%d Set %d/%d',vv,Nvoc, df,length(DataFiles)))
            xlabel(' ')
            set(gca, 'XTick',[],'XTickLabel',{})
        
            % Getting the posterior probability of a vocalization according
            % to resnet
            RowsData = logical(contains(ResNetData.Var1(:), sprintf('15120%s%s_%d_Logger', Date, DataFiles(df).name(1:11), InDataFilesOrdered(df))) .*contains(ResNetData.Var1(:), sprintf('_%d.mat', vv-1)));
            pVoc = ResNetData.Var3(RowsData); % all probability of a vocalization across all 100ms windows of all loggers through that sequence
            if any(pVoc>PPT(1))
                fprintf('p=%.2f ResNet predicts vocalizations within that sequence.\n', max(pVoc))
                % Listen to the microphone recording 
                pause(0.1)
                Raw_listen = filtfilt(sos_raw_band_listen,1,Raw_wave_nn);
                SampleMic = resample((Raw_listen - mean(Raw_listen))/(std(Raw_listen)/VolFactorMic),FS/4,FS);
                PlayerMic= audioplayer(SampleMic, FS/4,24); %#ok<TNMLP>
                play(PlayerMic)
                pause(length(Raw_wave_nn)/FS +1)
                % Take input from the curator if that is a vocalization
                ManCall=2;
                while (ManCall~=1) && (ManCall~=0)
                    commandwindow
                    ManCall = input('Did you hear vocalizations? yes (1) No (0) Play microphone (100)');
                    if isempty(ManCall)
                        ManCall=100;
                    end
                    if ManCall==100
                        play(PlayerMic)
                        pause(length(Raw_wave_nn)/FS +1)
                        pause(0.1)
                    end
                end
            else
                ManCall=2;
            end

            % skip to the next sequence if there is no vocalization
            if ~any(pVoc>PPT(1)) || ManCall==0
                if ManCall==0
                    fprintf(1,'Manual input enforced: Noise (0)\n');
                elseif ~any(pVoc>PPT(1))
                    fprintf('p=%.2f ResNet is not predicting any vocalization within that sequence. Skip\n', max(pVoc))
                end
                % Gather Vocalization production data:
                IndVocStart_all{vv} = [];
                IndVocStop_all{vv} = [];
                %                 IndNoiseStart_all{vv} = [];
                %                 IndNoiseStop_all{vv} = [];
                IndVocStartRaw_merged{vv} = [];
                IndVocStopRaw_merged{vv} = [];
                IndVocStartPiezo_merged{vv} = [];
                IndVocStopPiezo_merged{vv} = [];
                SorterName{vv} = Sorter;
                RMSRatio_all{vv} = [];
                RMSDiff_all{vv} = [];
                fprintf(1,'saving data...\n')
                if vvi==length(List2WorkOn) % this is the last vocalization
                    vv=Nvoc;
                end
                if ~isempty(dir(PreviousFile))
                    save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','vv','SorterName', '-append');
                else
                    %                     save(fullfile(Working_dir, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','IndNoiseStart_all','IndNoiseStop_all', 'IndNoiseStartRaw', 'IndNoiseStopRaw', 'IndNoiseStartPiezo', 'IndNoiseStopPiezo','RMSRatio_all','RMSDiff_all','vv','MicError','PiezoError','MicErrorType','PiezoErrorType');
                    save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','vv','SorterName','MicError','PiezoError','MicErrorType','PiezoErrorType');
                end
                clf(F1)
                if SaveRawWave
                    save(DataFile, 'Raw_wave','-append')
                end
                if SaveRawWaveName
                    save(DataFile, 'VocFilename','-append') %#ok<UNRCH> 
                    save(fullfile(Raw_dir, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)), 'Voc_filename','-append')
                end
                continue
            end

            %% Loop through the loggers and check the extracts length
            LengthLoggersData = nan(length(AudioLogs),1);
            for ll=1:length(AudioLogs)
                if ~isnan(Piezo_wave.(Fns_AL{ll}){vv})
                    LengthLoggersData(ll) = length(Piezo_wave.(Fns_AL{ll}){vv});
                end
            end
            %% Loop through the loggers and calculate envelopes
            Logger_Spec = cell(length(AudioLogs),1);
            ResNetPred = nan(length(AudioLogs),1);
            for ll=1:length(AudioLogs)
                if isnan(Piezo_FS.(Fns_AL{ll})(vv))
                    Piezo_FS.(Fns_AL{ll})(vv) = mean(Piezo_FS.(Fns_AL{ll}), 'omitnan');
                end
                if isempty(Piezo_wave.(Fns_AL{ll}){vv}) || any(isnan(Piezo_wave.(Fns_AL{ll}){vv}))
                    fprintf(1, 'NO DATA for Vocalization %d from %s\n', vv, Fns_AL{ll})
                else
                    % design the filters
                    [z,p,k] = butter(6,BandPassFilter(1:2)/(Piezo_FS.(Fns_AL{ll})(vv)/2),'bandpass');
                    sos_low = zp2sos(z,p,k);
                    [z,p,k] = butter(6,BandPassFilter(2:3)/(Piezo_FS.(Fns_AL{ll})(vv)/2),'bandpass');
                    sos_high = zp2sos(z,p,k);
                    % filter the loggers' signals
                    if sum(isnan(Piezo_wave.(Fns_AL{ll}){vv}))~=length(Piezo_wave.(Fns_AL{ll}){vv})
                        % replace NaN by zeros so the filter can work and
                        % check the length of extracts
                        InputPiezo = Piezo_wave.(Fns_AL{ll}){vv};
                        if sum(isnan(InputPiezo))
                            InputPiezo(isnan(InputPiezo))=0;
                            if length(InputPiezo)>min(LengthLoggersData)
                                warning('Piezo data have different durations!\n Logger %s data is truncated for analysis\n',Fns_AL{ll})
                                InputPiezo = InputPiezo(1:min(LengthLoggersData));
                            end
                            LowPassLogVoc{ll} = (filtfilt(sos_low,1,InputPiezo)); % low-pass filter the voltage trace
                            HighPassLogVoc = (filtfilt(sos_high,1,InputPiezo)); % high-pass filter the voltage trace
                        else
                            if length(InputPiezo)>min(LengthLoggersData)
                                warning('Piezo data have different durations!\n Logger %s data is truncated for analysis\n',Fns_AL{ll})
                                InputPiezo = InputPiezo(1:min(LengthLoggersData));
                            end
                            LowPassLogVoc{ll} = (filtfilt(sos_low,1,InputPiezo)); % low-pass filter the voltage trace
                            HighPassLogVoc = (filtfilt(sos_high,1,InputPiezo)); % high-pass filter the voltage trace
                        end
                        Amp_env_LowPassLogVoc{ll}=running_rms(LowPassLogVoc{ll}, Piezo_FS.(Fns_AL{ll})(vv), Fhigh_power, Fs_env);
                        Amp_env_HighPassLogVoc{ll}=running_rms(HighPassLogVoc, Piezo_FS.(Fns_AL{ll})(vv), Fhigh_power, Fs_env);
                    else
                        Amp_env_LowPassLogVoc{ll}=resample(nan(1,length(Piezo_wave.(Fns_AL{ll}){vv})), Fs_env, round(Piezo_FS.(Fns_AL{ll})(vv)));
                        Amp_env_HighPassLogVoc{ll}=resample(nan(1,length(Piezo_wave.(Fns_AL{ll}){vv})), Fs_env, round(Piezo_FS.(Fns_AL{ll})(vv)));
                    end
                    % Get ResNet posterior probability
                    RowsData = logical(contains(ResNetData.Var1(:), sprintf('15120%s%s_%d_%s_', Date, DataFiles(df).name(1:11), InDataFilesOrdered(df),Fns_AL{ll})) .*contains(ResNetData.Var1(:), sprintf('_%d.mat', vv-1)));
                    pVoclogger = ResNetData.Var3(RowsData); % all probability of a vocalization across all 100ms windows of that logger on that sequence
                    ResNetPred(ll) = max(pVoclogger);
                    if ~any(pVoclogger>PPT(2))
                        fprintf(1, 'p= %.2f No vocalization according to ResNet on %s\n', max(pVoclogger), Fns_AL{ll})
                    else
                        fprintf(1, 'p= %.2f ResNet found a vocalization on %s\n', max(pVoclogger), Fns_AL{ll})
                    end
                end
            end
            
            if ~(sum(cellfun('isempty',(Amp_env_LowPassLogVoc))) == length(AudioLogs))
                %% There is some data on the logger, extract the id of the vocalizing bat
                % Treat the case where some approximations of the calculations led to
                % slight diffreent sizes of running RMS
                if any(cellfun('isempty',(Amp_env_LowPassLogVoc))) % there is at least one missing logger
                    EmptyL = cellfun('isempty',(Amp_env_LowPassLogVoc));
                    Short = min(cellfun('length',Amp_env_LowPassLogVoc(~EmptyL)));
                    Amp_env_LowPassLogVoc{EmptyL} = nan(1,Short);
                    Amp_env_HighPassLogVoc{EmptyL} = nan(1,Short);
                else
                    Short = min(cellfun('length',Amp_env_LowPassLogVoc));
                end
                Long = max(cellfun('length',Amp_env_LowPassLogVoc));
                if (Long-Short)>1
                    error('The length of vectors of running RMS are too different than expected, please check!\n')
                elseif Long~=Short
                    Amp_env_LowPassLogVoc = cellfun(@(X) X(1:Short), Amp_env_LowPassLogVoc, 'UniformOutput',false);
                    Amp_env_HighPassLogVoc = cellfun(@(X) X(1:Short), Amp_env_HighPassLogVoc, 'UniformOutput',false);
                end
                Amp_env_LowPassLogVoc_MAT = cell2mat(Amp_env_LowPassLogVoc);
                Amp_env_HighPassLogVoc_MAT = cell2mat(Amp_env_HighPassLogVoc);
                RatioAmp = (Amp_env_LowPassLogVoc_MAT +1)./(Amp_env_HighPassLogVoc_MAT+1);
                DiffAmp = Amp_env_LowPassLogVoc_MAT-Amp_env_HighPassLogVoc_MAT;
                for ll=1:length(AudioLogs) % We are looking at actual logger!
                    Vocp_local = Amp_env_LowPassLogVoc_MAT(ll,:)>(Factor_RMS_low(ll) * RMSLow.(Fns_AL{ll})(1)); % Time points above amplitude threshold on the low-passed logger signal
                    IndVocStart_local = strfind(Vocp_local, ones(1,Consecutive_binsPiezo)); %find the first indices of every sequences of length "Consecutive_bins" higher than RMS threshold
                    if isempty(IndVocStart_local) || ResNetPred(ll)<=PPT(2) %#ok<STREMP> 
                        fprintf('\nNo vocalization detected on %s ResNet p=%.2f\n',Fns_AL{ll}, ResNetPred(ll));
                    else% Some vocalizations were detected plot the spectrogram
                        % Plot the low pass filtered signal of each logger
                        figure(1)
                        subplot(length(AudioLogs)+2,1,ll+1)
                        [Logger_Spec{ll}.to, Logger_Spec{ll}.fo, Logger_Spec{ll}.logB] = spec_only_bats(LowPassLogVoc{ll}, Piezo_FS.(Fns_AL{ll})(vv), DB_noise, FHigh_spec_Logger);
                        if ll<length(AudioLogs)
                            xlabel(' '); % supress the x label output
                            set(gca,'XTick',[],'XTickLabel',{});
                        end
                        hold on
                        yyaxis right
                        plot((1:length(Amp_env_LowPassLogVoc{ll}))/Fs_env*1000, Amp_env_LowPassLogVoc{ll}, 'b-','LineWidth', 2);
                        hold on
                        plot((1:length(Amp_env_HighPassLogVoc{ll}))/Fs_env*1000, Amp_env_HighPassLogVoc{ll}, 'r-','LineWidth',2);
                        ylabel(sprintf('Amp\n%s',Fns_AL{ll}([1:3 7:end])))
                        hold off
                        %                         if Manual
                        %                             pause(0.1)
                        %                             Player= audioplayer((Piezo_wave.(Fns_AL{ll}){vv}-mean(Piezo_wave.(Fns_AL{ll}){vv}))/std(Piezo_wave.(Fns_AL{ll}){vv}), Piezo_FS.(Fns_AL{ll})(vv)); %#ok<TNMLP>
                        %                             play(Player)
                        %                             pause(length(Raw_wave{vv})/FS +1)
                        %                         end
                    end
                end
            end

            
            %% Listen to microphone and decide if we should go anyfurther in analysis
            % Player
            ManCall=2;
            while (ManCall~=1) && (ManCall~=0)
                commandwindow
                ManCall = input('Did you hear vocalizations? yes (1) No (0) Play microphone and loggers (2) Play microphone (100) Play logger #x (x)');
                if isempty(ManCall)
                    ManCall=2;
                   continue
                end
                if ManCall==2
                    play(PlayerMic)
                    pause(length(Raw_wave_nn)/FS +1)
                    pause(0.1)
                    for ll=1:length(AudioLogs)
                        Player= audioplayer((Piezo_wave.(Fns_AL{ll}){vv}-mean(Piezo_wave.(Fns_AL{ll}){vv}))/(VolDenominatorLogger*std(Piezo_wave.(Fns_AL{ll}){vv})), Piezo_FS.(Fns_AL{ll})(vv)); %#ok<TNMLP>
                        play(Player)
                        pause(length(Raw_wave_nn)/FS +1)
                    end
                elseif ManCall==100
                    play(PlayerMic)
                    pause(length(Raw_wave_nn)/FS +1)
                    pause(0.1)
                elseif (ManCall~=1) && (ManCall~=0)
                    FieldLog = sprintf('Logger%d',ManCall);
                    ll = find(strcmp(FieldLog, Fns_AL));
                    if isempty(ll)
                        warning('Invalid Logger number, Logger %d does not exist',ManCall)
                    else
                        Player= audioplayer((Piezo_wave.(Fns_AL{ll}){vv}-mean(Piezo_wave.(Fns_AL{ll}){vv}, 'omitnan'))/(VolDenominatorLogger*std(Piezo_wave.(Fns_AL{ll}){vv}, 'omitnan')), Piezo_FS.(Fns_AL{ll})(vv)); %#ok<TNMLP> 
                        play(Player)
                        pause(length(Raw_wave_nn)/FS +1)
                    end
                end
            end
            
            % Decision
            if ManCall==0
                fprintf(1,'Manual input enforced: Noise (0)\n');
                % Gather Vocalization production data:
                IndVocStart_all{vv} = [];
                IndVocStop_all{vv} = [];
%                 IndNoiseStart_all{vv} = [];
%                 IndNoiseStop_all{vv} = [];
                IndVocStartRaw_merged{vv} = [];
                IndVocStopRaw_merged{vv} = [];
                IndVocStartPiezo_merged{vv} = [];
                IndVocStopPiezo_merged{vv} = [];
                SorterName{vv} = Sorter;
                RMSRatio_all{vv} = [];
                RMSDiff_all{vv} = [];
                fprintf(1,'saving data...\n')
                if vvi==length(List2WorkOn) % this is the last vocalization
                    vv=Nvoc;
                end
                if ~isempty(dir(PreviousFile))
                    save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','vv','SorterName', '-append');
                else
%                     save(fullfile(Working_dir, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','IndNoiseStart_all','IndNoiseStop_all', 'IndNoiseStartRaw', 'IndNoiseStopRaw', 'IndNoiseStartPiezo', 'IndNoiseStopPiezo','RMSRatio_all','RMSDiff_all','vv','MicError','PiezoError','MicErrorType','PiezoErrorType');
                    save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','vv','SorterName','MicError','PiezoError','MicErrorType','PiezoErrorType');
                end
                clf(F1)
                if SaveRawWave
                    save(DataFile, 'Raw_wave','-append')
                end
                if SaveRawWaveName
                    save(DataFile, 'VocFilename','-append') %#ok<UNRCH> 
                    save(fullfile(Raw_dir, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)), 'Voc_filename','-append')
                end
                continue
            end
                
             %% Now find onset/offset on microphone recording if there is no logger data      
            if sum(cellfun('isempty',(Amp_env_LowPassLogVoc))) == length(AudioLogs) % No logger data, just isolate onset/offset of vocalizations on the microphone
                if ~Force_Save_onoffsets_mic
                    fprintf(1,'CANNOT DETERMINE OWNERSHIP NO DETECTION OF ONSET/OFFSET\nSet Force_Save_onoffsets_mic to 1 to force saving data\n')
                else
                    fprintf(1,'CANNOT DETERMINE OWNERSHIP JUST SAVING ONSET/OFFSET OFF THE MICROPHONE\n')
                    if RMS_Fig
                        F2=figure(2); %#ok<UNRCH> 
                        clf(F2)
                        subplot(3,1,1)
                        if CheckMicChannel
                            plot(Amp_env_Mic.*10^3,'LineWidth',2, 'Color', ColorCode(ll+1,:))
                            hold on
                            plot([0 length(Amp_env_Mic)], (Factor_RMS_Mic * MeanStdAmpRawExtract(Nvocs(df)+vv,1))*ones(2,1).*10^3, 'Color',ColorCode(ll+1,:),'LineStyle','--');
                            hold on
                        end
                        legend({'Microphone' 'voc detection threshold' 'voc detection threshold'})
                    end
                    RowSize = length(AudioLogs) +2;
                
                    IndVocStartRaw{vv} = cell(RowSize,1);
                    IndVocStopRaw{vv} = cell(RowSize,1);
                    IndVocStart = cell(RowSize,1);
                    IndVocStop = cell(RowSize,1);
%                     IndNoiseStartRaw{vv} = cell(RowSize,1);
%                     IndNoiseStopRaw{vv} = cell(RowSize,1);
%                     IndNoiseStart = cell(RowSize,1);
%                     IndNoiseStop = cell(RowSize,1);
                    IndVocStartRaw_merge_local = cell(RowSize,1);
                    IndVocStopRaw_merge_local = cell(RowSize,1);
                    IndVocStartPiezo_merge_local = cell(RowSize,1);
                    IndVocStopPiezo_merge_local = cell(RowSize,1);
                
                    VocpMic = Amp_env_Mic>(Factor_RMS_Mic * MeanStdAmpRawExtract(Nvocs(df)+vv,1)); % Time points above amplitude threshold on the band-pass microphone signal
                    VocpMic = reshape(VocpMic,1,length(VocpMic));
                    IndVocStart{RowSize} = strfind(VocpMic, ones(1,Consecutive_binsMic)); %find the first indices of every sequences of length "Consecutive_bins" higher than RMS threshold
                    if isempty(IndVocStart{RowSize}) || ResNetPred(ll)<=PPT(2)
                        fprintf(1,'No vocalization detected on microphone\n');
                    else% Some vocalizations are automatically detected
                        IndVocStart_diffind = find(diff(IndVocStart{RowSize})>1);
                        IndVocStart{RowSize} = [IndVocStart{RowSize}(1) IndVocStart{RowSize}(IndVocStart_diffind +1)]; % these two lines get rid of overlapping sequences that werer detected several times
                        NV = length(IndVocStart{RowSize}); % This is the number of detected potential vocalization
                        IndVocStop{RowSize} = nan(1,NV);
                        NewCall1Noise0_temp = ones(NV,1);
                        if ManCall
                            NewCall1Noise0_man = nan(NV,1);
                        else
                            NewCall1Noise0_man = zeros(NV,1);
                        end
                        for ii=1:NV
                            IVStop = find(VocpMic(IndVocStart{RowSize}(ii):end)==0, 1, 'first');
                            if ~isempty(IVStop)
                                IndVocStop{RowSize}(ii) = IndVocStart{RowSize}(ii) + IVStop -1;
                            else
                                IndVocStop{RowSize}(ii) = length(VocpMic);
                            end
                            % Plot the localization of that sound extract on figure 3
                            % Duplicate the figure of the spectrogram for manual input purposes
                            F3 = figure(3);
                            if ii==1
                                clf(F3)
                                maxB = max(max(Raw_Spec.logB));
                                minB = maxB-DB_noise;            
                                imagesc(Raw_Spec.to*1000,Raw_Spec.fo,Raw_Spec.logB);          % to is in seconds
                                axis xy;
                                caxis('manual');
                                caxis([minB maxB]); 
                                cmap = spec_cmap();
                                colormap(cmap);
                                v_axis = axis; 
                                v_axis(3)=0; 
                                v_axis(4)=FHigh_spec;
                                axis(v_axis);                                
                                xlabel('time (ms)'), ylabel('Frequency');
                                title(sprintf('Ambient Microphone Voc %d/%d Set %d/%d',vv,Nvoc,df,length(DataFiles)))
                            end
                                %yyaxis right
                            %                                 ylabel('Logger ID')
                            %                                 set(gca, 'YTick', 1:RowSize, 'YTickLabel', [Fns_AL; 'Mic'], 'YLim', [0 (length(AudioLogs)+2)],'YDir', 'reverse')
                            %                             end
                            hold on
                            yyaxis right
                            plot([IndVocStart{RowSize}(ii)/Fs_env IndVocStop{RowSize}(ii)/Fs_env]*1000, [RowSize RowSize], 'k:', 'LineWidth',2)
                            hold off

                            % Decide if that call was already detected and ID
                            % attributed or not by looking if there is any
                            % overlap
                            %                     OverlapOnset = nan(RowSize-1,1);
                            %                     OverlapOffset = nan(RowSize-1,1);
                            %                     for ll2=1:(RowSize-1)
                            %                         OverlapOnset(ll2) = any((IndVocStart{ll2}<=IndVocStart{RowSize}(ii)) .* (IndVocStop{ll2}>=IndVocStart{RowSize}(ii)));
                            %                         OverlapOffset(ll2) = any((IndVocStop{ll2}<=IndVocStart{RowSize}(ii)) .* (IndVocStop{ll2}>=IndVocStop{RowSize}(ii)));
                            %                     end

                            %                     NewCall1OldCall0_temp(ii) = ~any([OverlapOnset; OverlapOffset]);

                            % update figure(3) with the decision
                            figure(3)
                            yyaxis right
                            hold on
                            text(IndVocStop{RowSize}(ii)/Fs_env*1000, RowSize, 'New call on Mic')
                            hold off

                            fprintf(1,'Computer guess for that sound element: New call on Mic\n');
                            if ManCall
                                TempIn = input('Indicate your choice: new call (1);    noise (0);    listen again to mic(any other number)\n');
                                if isempty(TempIn)
                                    fprintf(1, 'NO ENTRY GIVEN, playing again the sound and asking the same question\n')
                                    NewCall1Noise0_man(ii) = 2;
                                else
                                    NewCall1Noise0_man(ii) = TempIn;
                                end
                                while NewCallNoise0_man(ii)~=0 && NewCallNoise0_man(ii)~=1
                                    play(PlayerMic)
                                    pause(length(Raw_wave_nn)/FS +1)
                                    TempIn = input('Indicate your choice: new call (1);    noise (0);    listen again to mic(any other number)\n');
                                    if isempty(TempIn)
                                        fprintf(1, 'NO ENTRY GIVEN, playing again the sound and asking the same question\n')
                                        NewCall1Noise0_man(ii) = 2;
                                    else
                                        NewCall1Noise0_man(ii) = TempIn;
                                    end
                                end
                            else
                                fprintf(1,'Manual input enforced: Noise (0)\n');
                            end
                            Agree = NewCall1Noise0_man(ii) == NewCall1Noise0_temp(ii);
                            if ~Agree
                                NewCall1Noise0_temp(ii) = NewCall1Noise0_man(ii);
                                MicError = MicError + [1 1];
                                MicErrorType = MicErrorType + [NewCall1Noise0_temp(ii) ~NewCall1Noise0_temp(ii)];
                            else
                                MicError = MicError + [0 1];
                            end
                        
                            % update figure(3) with the decision
                            figure(3)
                            yyaxis right
                            hold on
                            if NewCall1Noise0_temp(ii)
                                plot([IndVocStart{RowSize}(ii)/Fs_env IndVocStop{RowSize}(ii)/Fs_env]*1000, [RowSize RowSize], 'Color',ColorCode(RowSize,:), 'LineWidth',2, 'LineStyle','-')
                                text(IndVocStop{RowSize}(ii)/Fs_env*1000, RowSize+0.2, 'New call','Color','r','FontWeight','bold')
                            else
                                plot([IndVocStart{RowSize}(ii)/Fs_env IndVocStop{RowSize}(ii)/Fs_env]*1000, [RowSize RowSize], 'k-', 'LineWidth',2)
                                text(IndVocStop{RowSize}(ii)/Fs_env*1000, RowSize+0.2, 'Noise','Color','r','FontWeight','bold')
                            end
                            hold off                        
                        end
                    
                        % Stash noise and only keep vocalizations
%                         IndNoiseStart{RowSize} = IndVocStart{RowSize}(~logical(NewCall1Noise0_temp));
%                         IndNoiseStop{RowSize} = IndVocStop{RowSize}(~logical(NewCall1Noise0_temp));
%                         IndNoiseStartRaw{vv}{RowSize} = round(IndNoiseStart{RowSize}/Fs_env*FS);
%                         IndNoiseStopRaw{vv}{RowSize} = round(IndNoiseStop{RowSize}/Fs_env*FS);
                        IndVocStart{RowSize} = IndVocStart{RowSize}(logical(NewCall1Noise0_temp));
                        IndVocStop{RowSize} = IndVocStop{RowSize}(logical(NewCall1Noise0_temp));
                        IndVocStartRaw{vv}{RowSize} = round(IndVocStart{RowSize}/Fs_env*FS);
                        IndVocStopRaw{vv}{RowSize} = round(IndVocStop{RowSize}/Fs_env*FS);
                    
                        NV = sum(NewCall1Noise0_temp);
                    
                        if NV
                            % Now merge detected vocalizations that are within MergeThresh
                            Merge01 = (IndVocStart{RowSize}(2:end)-IndVocStop{RowSize}(1:(end-1)) < (MergeThresh/1000*Fs_env));
                            NV = NV-sum(Merge01);
                            IndVocStartRaw_merge_local{RowSize} = nan(1,NV);
                            IndVocStopRaw_merge_local{RowSize} = nan(1,NV);
                            IndVocStartPiezo_merge_local{RowSize} = nan(1,NV); % There is no useful piezo data for a call detected only on Microphone, it will stay Nan
                            IndVocStopPiezo_merge_local{RowSize} = nan(1,NV);% There is no useful piezo data for a call detected only on Microphone, it will stay Nan
                            CutInd = find(~Merge01);
                            IndVocStartRaw_merge_local{RowSize}(1) = round(IndVocStart{RowSize}(1)/Fs_env*FS);
                            IndVocStopRaw_merge_local{RowSize}(end) = round(IndVocStop{RowSize}(end)/Fs_env*FS);
                            for cc=1:length(CutInd)
                                IndVocStopRaw_merge_local{RowSize}(cc) = round(IndVocStop{RowSize}( CutInd(cc))/Fs_env*FS);
                                IndVocStartRaw_merge_local{RowSize}(cc+1) = round(IndVocStart{RowSize}(CutInd(cc)+1)/Fs_env*FS);
                            end

                            % Now plot the onset/offset of each extract on the
                            % spectrograms of the microphone
                            figure(1)
                            for ii=1:NV
                                subplot(length(AudioLogs)+2,1,1)
                                hold on
                                yyaxis right
                                YLim = get(gca, 'YLim');
                                plot([IndVocStartRaw_merge_local{RowSize}(ii) IndVocStopRaw_merge_local{RowSize}(ii)]/FS*1000, [YLim(2)*4/5 YLim(2)*4/5], 'k-', 'LineWidth',2)
                                subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
                                hold on
                                plot([IndVocStartRaw_merge_local{RowSize}(ii) IndVocStopRaw_merge_local{RowSize}(ii)]/FS*1000, [RowSize RowSize], 'Color',ColorCode(RowSize,:), 'LineWidth',2, 'LineStyle','-')
                                hold off
                            end

                        end
                    end
                    [~,FileVoc]=fileparts(VocFilename{vv});
                    fprintf(1,'saving figures...\n')
                    if strcmp(SaveFileType,'pdf')
                        print(F1,fullfile(Working_dir_write,sprintf('%s_whocalls_spec_%d.pdf', FileVoc, MergeThresh)),'-dpdf','-fillpage')
                        if sum(cellfun(@length, IndVocStartRaw_merge_local)) && RMS_Fig % Only save the RMS figure if there was a vocalization
                            saveas(F2,fullfile(Working_dir_write,sprintf('%s_whocalls_RMS_%d.pdf', FileVoc, MergeThresh)),'pdf')
                        end
                    elseif strcmp(SaveFileType,'fig')
                        saveas(F1,fullfile(Working_dir_write,sprintf('%s_whocalls_spec_%d.fig', FileVoc, MergeThresh)))
                        if sum(cellfun(@length, IndVocStartRaw_merge_local)) && RMS_Fig % Only save the RMS figure if there was a vocalization
                            saveas(F2,fullfile(Working_dir_write,sprintf('%s_whocalls_RMS_%d.fig', FileVoc, MergeThresh)))
                        end
                    end
                    if ManCall
                        pause(1)
                    end
                    clf(F1)
                    if RMS_Fig
                        clf(F2) %#ok<UNRCH> 
                    end
                
                end
            else
                %% There is some data on the logger, extract the id of the vocalizing bat
                % Plot the ratio of time varying RMS, the difference in time varying
                % RMS between the high and low frequency bands and the absolute time
                % varying RMS of the low frequency band
                if RMS_Fig
                    F2=figure(2); %#ok<UNRCH> 
                    clf(F2)
                    subplot(3,1,3)
                    for ll=1:length(AudioLogs)
                        plot(RatioAmp(ll,:), 'LineWidth',2, 'Color', ColorCode(ll,:))
                        hold on
                    end
                    ylabel('Frequency bands Ratio')
                    %         hold on
                    %         for ll=1:length(AudioLogs)
                    %             plot([0 size(Amp_env_LowPassLogVoc_MAT,2)], Factor_AmpRatio * RatioRMS.(Fns_AL{ll})(1)*ones(2,1), 'Color',ColorCode(ll,:),'LineStyle','--');
                    %             hold on
                    %         end
                    %         legend({Fns_AL{:} 'calling detection threshold' 'calling detection threshold'})
                    legend(Fns_AL{:})
                    subplot(3,1,2)
                    for ll=1:length(AudioLogs)
                        plot(DiffAmp(ll,:), 'LineWidth',2, 'Color', ColorCode(ll,:))
                        hold on
                    end
                    title('Calling vs Hearing')
                    ylabel('Frequency bands Diff')
                    hold on
                    for ll=1:length(AudioLogs)
                        plot([0 size(Amp_env_LowPassLogVoc_MAT,2)], Factor_AmpDiff * DiffRMS.(Fns_AL{ll})(1)*ones(2,1), 'Color',ColorCode(ll,:),'LineStyle','--');
                        hold on
                    end
                    legend({Fns_AL{:} 'calling detection threshold' 'calling detection threshold'})
                    subplot(3,1,1)
                    for ll=1:length(AudioLogs)
                        plot(Amp_env_LowPassLogVoc_MAT(ll,:), 'LineWidth',2, 'Color', ColorCode(ll,:))
                        hold on
                    end
                    hold on
                    ylabel(sprintf('Running RMS %dHz-%dHz', BandPassFilter(1:2)))
                    title('Detection of vocalizations on each logger')
                    for ll=1:length(AudioLogs)
                        plot([0 size(Amp_env_LowPassLogVoc_MAT,2)], Factor_RMS_low(ll) * RMSLow.(Fns_AL{ll})(1)*ones(2,1), 'Color',ColorCode(ll,:),'LineStyle','--');
                        hold on
                    end

                    if CheckMicChannel
                        plot(Amp_env_Mic.*10^3,'LineWidth',2, 'Color', ColorCode(ll+1,:))
                        hold on
                        plot([0 length(Amp_env_Mic)], (Factor_RMS_Mic * MeanStdAmpRawExtract(Nvocs(df)+vv,1))*ones(2,1).*10^3, 'Color',ColorCode(ll+1,:),'LineStyle','--');
                        hold on
                    end
                    legend({Fns_AL{:} 'Microphone' 'voc detection threshold' 'voc detection threshold'})
                end
                %% Find out which calls are emitted by each individual, if checking the microphone is requested it means that a single animal did not have a collar and all microphone calls not detected on any piezo will be attributed to it
                Vocp = nan(size(DiffAmp) + [1 0]); % This is the output of the decision criterion to determine at each time point if a bat is vocalizing or not. each row=a logger, each column = a time point
                if CheckMicChannel
                    RowSize = length(AudioLogs) +1;
                else
                    RowSize = length(AudioLogs);
                end
                IndVocStartRaw{vv} = cell(RowSize,1);
                IndVocStartPiezo{vv} = cell(RowSize,1);
                IndVocStopRaw{vv} = cell(RowSize,1);
                IndVocStopPiezo{vv} = cell(RowSize,1);
                IndVocStart = cell(RowSize,1);
                IndVocStop = cell(RowSize,1);
%                 IndNoiseStartRaw{vv} = cell(RowSize,1);
%                 IndNoiseStartPiezo{vv} = cell(RowSize,1);
%                 IndNoiseStopRaw{vv} = cell(RowSize,1);
%                 IndNoiseStopPiezo{vv} = cell(RowSize,1);
%                 IndNoiseStart = cell(RowSize,1);
%                 IndNoiseStop = cell(RowSize,1);
                IndVocStartRaw_merge_local = cell(RowSize,1);
                IndVocStopRaw_merge_local = cell(RowSize,1);
                IndVocStartPiezo_merge_local = cell(RowSize,1);
                IndVocStopPiezo_merge_local = cell(RowSize,1);
                %     Xcor = cell(length(AudioLogs),1);
                RMSRatio = cell(RowSize,1);
                RMSDiff = cell(RowSize,1);
                for ll=1:RowSize
                    % This is the Microphone condition
                    if ll>length(AudioLogs) && CheckMicChannel % we are looking at the microphone data now
                        VocpMic = Amp_env_Mic>(Factor_RMS_Mic * MeanStdAmpRawExtract(Nvocs(df)+vv,1)); % Time points above amplitude threshold on the band-pass microphone signal
                        VocpMic = reshape(VocpMic,1,length(VocpMic));
                        IndVocStart{ll} = strfind(VocpMic, ones(1,Consecutive_binsMic)); %find the first indices of every sequences of length "Consecutive_bins" higher than RMS threshold
                        if isempty(IndVocStart{ll})
                            fprintf(1,'\nNo vocalization detected on microphone\n');
                        else% Some vocalizations were detected
                            IndVocStart_diffind = find(diff(IndVocStart{ll})>1);
                            IndVocStart{ll} = [IndVocStart{ll}(1) IndVocStart{ll}(IndVocStart_diffind +1)]; % these two lines get rid of overlapping sequences that werer detected several times
                            NV = length(IndVocStart{ll}); % This is the number of detected potential vocalization
                            IndVocStop{ll} = nan(1,NV);
                            NewCall1OldCall0_temp = nan(NV,1);
                            if ManCall
                                NewCall1OldCall0_man = nan(NV,1);
                            else
                                NewCall1OldCall0_man = zeros(NV,1);
                            end
                            for ii=1:NV
                                IVStop = find(VocpMic(IndVocStart{ll}(ii):end)==0, 1, 'first');
                                if ~isempty(IVStop)
                                    IndVocStop{ll}(ii) = IndVocStart{ll}(ii) + IVStop -1;
                                else
                                    IndVocStop{ll}(ii) = length(VocpMic);
                                end
                                % Plot the localization of that sound extract on figure 3
                                % Duplicate the figure of the spectrogram for manual input purposes
                                F3 = figure(3);
                                %
                                if ii==1
                                    clf(F3)
                                    maxB = max(max(Raw_Spec.logB));
                                    minB = maxB-DB_noise;
                                    imagesc(Raw_Spec.to*1000,Raw_Spec.fo,Raw_Spec.logB);          % to is in seconds
                                    axis xy;
                                    caxis('manual');
                                    caxis([minB maxB]);
                                    cmap = spec_cmap();
                                    colormap(cmap);
                                    v_axis = axis;
                                    v_axis(3)=0;
                                    v_axis(4)=FHigh_spec;
                                    axis(v_axis);
                                    xlabel('time (ms)'), ylabel('Frequency');
                                    title(sprintf('Ambient Microphone Voc %d/%d Set %d/%d',vv,Nvoc, df,length(DataFiles)))
%                                 yyaxis right
                                %                                 ylabel('Logger ID')
                                %                                 set(gca, 'YTick', 1:ll, 'YTickLabel', [Fns_AL; 'Mic'], 'YLim', [0 (length(AudioLogs)+2)],'YDir', 'reverse')
                                %                             end
                                end
                                hold on
                                yyaxis right
                                plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [ll ll], 'k:', 'LineWidth',2)
                                ylabel('Microphone')
                                set(gca, 'YTickLabel', [])
                                hold off
                            
                                % Decide if that call was already detected and ID
                                % attributed or not by looking if there is any
                                % overlap
                                OverlapOnset = nan(ll-1,1);
                                OverlapOffset = nan(ll-1,1);
                                for ll2=1:(ll-1)
                                    OverlapOnset(ll2) = any((IndVocStart{ll2}<=IndVocStart{ll}(ii)) .* (IndVocStop{ll2}>=IndVocStart{ll}(ii)));
                                    OverlapOffset(ll2) = any((IndVocStop{ll2}<=IndVocStart{ll}(ii)) .* (IndVocStop{ll2}>=IndVocStop{ll}(ii)));
                                end

                                NewCall1OldCall0_temp(ii) = ~any([OverlapOnset; OverlapOffset]);

                                % update figure(3) with the decision
                                figure(3)
                                yyaxis right
                                hold on
                                if NewCall1OldCall0_temp(ii)
                                    text(IndVocStart{ll}(ii)/Fs_env*1000, ll+0.4, 'New call')
                                else
                                    text(IndVocStart{ll}(ii)/Fs_env*1000, ll+0.4, 'Known/Noise')
                                end
                                hold off
                            
                                if NewCall1OldCall0_temp(ii)
                                    fprintf(1,'\nComputer guess for that sound element: New call on Mic\n');
                                else
                                    fprintf(1,'\nComputer guess for that sound element: known Call already attributed\n');
                                end
                                if ManCall
                                    TempIn = input('Indicate your choice: new call on Mic (1);    already known/noise (0);    listen to Mic again(any other number)\n');
                                    if isempty(TempIn)
                                        fprintf(1, 'NO ENTRY GIVEN, playing again the sound and asking the same question\n')
                                        NewCall1OldCall0_man(ii)=2;
                                    else
                                        NewCall1OldCall0_man(ii) = TempIn;
                                    end
                                    while NewCall1OldCall0_man(ii)~=0 && NewCall1OldCall0_man(ii)~=1
                                        Raw_listen = filtfilt(sos_raw_band_listen,1,Raw_wave_nn);
                                        SampleMic = resample((Raw_listen - mean(Raw_listen))/(std(Raw_listen)/VolFactorMic),FS/4,FS);
                                        PlayerMic= audioplayer(SampleMic, FS/4,24); %#ok<TNMLP>
                                        play(PlayerMic)
                                        pause(length(Raw_wave_nn)/FS +1)
                                        TempIn = input('Indicate your choice: new call on Mic (1);    already known/noise (0);    listen to mic again(any other number)\n');
                                        if isempty(TempIn)
                                            fprintf(1, 'NO ENTRY GIVEN, playing again the sound and asking the same question\n')
                                            NewCall1OldCall0_man(ii)=2;
                                        else
                                            NewCall1OldCall0_man(ii) = TempIn;
                                        end
                                    end
                                else
                                    fprintf(1,'Manual input enforced: Noise (0)\n');
                                end

                                Agree = NewCall1OldCall0_man(ii)==NewCall1OldCall0_temp(ii);
                                if ~Agree
                                    NewCall1OldCall0_temp(ii) = NewCall1OldCall0_man(ii);
                                    MicError = MicError + [1 1];
                                    MicErrorType = MicErrorType + [NewCall1OldCall0_temp(ii) ~NewCall1OldCall0_temp(ii)];
                                else
                                    MicError = MicError + [0 1];
                                end
                            
                                % update figure(3) with the decision
                                figure(3)
                                yyaxis right
                                hold on
                                if NewCall1OldCall0_temp(ii)
                                    plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
                                    text(IndVocStart{ll}(ii)/Fs_env*1000, ll+0.2, 'New call','Color','r','FontWeight','bold')
                                else
                                    plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [ll ll], 'k-', 'LineWidth',2)
                                    text(IndVocStart{ll}(ii)/Fs_env*1000, ll+0.2, 'Known/Noise','Color','r','FontWeight','bold')
                                end
                                hold off
                            
                            end
                        
                            % Stash noise and Only keep vocalizations that were
                            % not detected on piezos
%                             IndNoiseStart{ll} = IndVocStart{ll}(~logical(NewCall1OldCall0_temp));
%                             IndNoiseStop{ll} = IndVocStop{ll}(~logical(NewCall1OldCall0_temp));
%                             IndNoiseStartRaw{vv}{ll} = round(IndNoiseStart{ll}/Fs_env*FS);
%                             IndNoiseStopRaw{vv}{ll} = round(IndNoiseStop{ll}/Fs_env*FS);
                            IndVocStart{ll} = IndVocStart{ll}(logical(NewCall1OldCall0_temp));
                            IndVocStop{ll} = IndVocStop{ll}(logical(NewCall1OldCall0_temp));
                            IndVocStartRaw{vv}{ll} = round(IndVocStart{ll}/Fs_env*FS);
                            IndVocStopRaw{vv}{ll} = round(IndVocStop{ll}/Fs_env*FS);
                            NV = sum(NewCall1OldCall0_temp);
                        
                            if NV
                                % Now merge detected vocalizations that are within MergeThresh
                                Merge01 = (IndVocStart{ll}(2:end)-IndVocStop{ll}(1:(end-1)) < (MergeThresh/1000*Fs_env));
                                NV = NV-sum(Merge01);
                                IndVocStartRaw_merge_local{ll} = nan(1,NV);
                                IndVocStopRaw_merge_local{ll} = nan(1,NV);
                                IndVocStartPiezo_merge_local{ll} = nan(1,NV); % There is no useful piezo data for a call detected only on Microphone, it will stay Nan
                                IndVocStopPiezo_merge_local{ll} = nan(1,NV);% There is no useful piezo data for a call detected only on Microphone, it will stay Nan
                                CutInd = find(~Merge01);
                                IndVocStartRaw_merge_local{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*FS);
                                IndVocStopRaw_merge_local{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*FS);
                                for cc=1:length(CutInd)
                                    IndVocStopRaw_merge_local{ll}(cc) = round(IndVocStop{ll}( CutInd(cc))/Fs_env*FS);
                                    IndVocStartRaw_merge_local{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*FS);
                                end

                                % Now plot the onset/offset of each extract on the
                                % spectrograms of the microphone
                                figure(1)
                                for ii=1:NV
                                    subplot(length(AudioLogs)+2,1,1)
                                    hold on
                                    yyaxis right
                                    YLim = get(gca, 'YLim');
                                    plot([IndVocStartRaw_merge_local{ll}(ii) IndVocStopRaw_merge_local{ll}(ii)]/FS*1000, [YLim(2)*4/5 YLim(2)*4/5], 'k-', 'LineWidth',2)
                                    subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
                                    hold on
                                    plot([IndVocStartRaw_merge_local{ll}(ii) IndVocStopRaw_merge_local{ll}(ii)]/FS*1000, [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
                                    hold off
                                end

                            end
                        end
                    % Now we are looking at an actual logger    
                    elseif ll<=length(AudioLogs) % We are looking at an actual logger!
                        Vocp(ll,:) = Amp_env_LowPassLogVoc_MAT(ll,:)>(Factor_RMS_low(ll) * RMSLow.(Fns_AL{ll})(1)); % Time points above amplitude threshold on the low-passed logger signal
                        IndVocStart{ll} = strfind(Vocp(ll,:), ones(1,Consecutive_binsPiezo)); %find the first indices of every sequences of length "Consecutive_bins" higher than RMS threshold
                        if isempty(IndVocStart{ll}) || ResNetPred(ll)<=PPT(2)
                            fprintf('\nNo vocalization detected on %s\n',Fns_AL{ll});
                        else% Some vocalizations were detected
                        
                            if ManCall
                                fprintf(1,'\nSound event detected on %s.',Fns_AL{ll})
                                LoggerSound = Piezo_wave.(Fns_AL{ll}){vv}-mean(Piezo_wave.(Fns_AL{ll}){vv}, 'omitnan');
                                Player= audioplayer(LoggerSound/(VolDenominatorLogger*std(LoggerSound, 'omitnan')), Piezo_FS.(Fns_AL{ll})(vv)); %#ok<TNMLP>
                                play(Player)
                                pause(length(Raw_wave_nn)/FS +1)
                                ManCall_logger=input(sprintf('\nDid you hear any call on %s? (yes:1 ; No:0; Yes but go auto: 100 ; Yes but listen to each call: 101 ; listen again to that logger recording (any other number) )\n',Fns_AL{ll}));
                                if isempty(ManCall_logger)
                                    ManCall_logger=2;
                                end
                                while ManCall_logger~=0 && ManCall_logger~=1 && ManCall_logger~=100 && ManCall_logger~=101
                                    play(Player)
                                    pause(length(Raw_wave_nn)/FS +1)
                                    ManCall_logger = input(sprintf('\nDid you hear any call on %s? (yes:1 ; No:0  ; Yes but go auto: 100 ;Yes but listen to each call: 101 ; listen again to that logger recording (any other number) )\n',Fns_AL{ll}));
                                    if isempty(ManCall_logger)
                                        ManCall_logger=2;
                                    end
                                end
                            else
                                fprintf(1,'\nSound event detected on %s\n',Fns_AL{ll});
                                ManCall_logger = 0;
                            end
                            
                            if ManCall && ManCall_logger
                                IndVocStart_diffind = find(diff(IndVocStart{ll})>1);
                                IndVocStart{ll} = [IndVocStart{ll}(1) IndVocStart{ll}(IndVocStart_diffind +1)]; % these two lines get rid of overlapping sequences that werer detected several times
                                NV = length(IndVocStart{ll}); % This is the number of detected potential vocalization
                                IndVocStop{ll} = nan(1,NV);
                                %             Xcor{ll} = nan(NV,1);
                                RMSRatio{ll} = nan(NV,1);
                                RMSDiff{ll} = nan(NV,1);
                                Call1Hear0_temp = nan(NV,1);

                                Call1Hear0_man = nan(NV,1);
                                
                                for ii=1:NV
                                    IVStop = find(Vocp(ll,IndVocStart{ll}(ii):end)==0, 1, 'first');
                                    if ~isempty(IVStop)
                                        IndVocStop{ll}(ii) = IndVocStart{ll}(ii) + IVStop -1;
                                    else
                                        IndVocStop{ll}(ii) = length(Vocp(ll,:));
                                    end
                                    % Plot the localization of that sound extract on figure 3
                                    % Duplicate the figure of the spectrogram for manual input purposes
                                    %                             if ii==1
                                    F3 = figure(3);
                                    if ii==1
                                        clf(F3)
                                        maxB = max(max(Raw_Spec.logB));
                                        minB = maxB-DB_noise;
                                        subplot(2,1,1)
                                        imagesc(Raw_Spec.to*1000,Raw_Spec.fo,Raw_Spec.logB);          % to is in seconds
                                        axis xy;
                                        caxis('manual');
                                        caxis([minB maxB]);
                                        cmap = spec_cmap();
                                        colormap(cmap);
                                        v_axis = axis;
                                        v_axis(3)=0;
                                        v_axis(4)=FHigh_spec;
                                        axis(v_axis);
                                        xlabel('time (ms)'), ylabel('Frequency');
                                        
                                        
                                        subplot(2,1,2)
                                        imagesc(Logger_Spec{ll}.to*1000,Logger_Spec{ll}.fo,Logger_Spec{ll}.logB);          % to is in seconds
                                        maxB = max(max(Logger_Spec{ll}.logB));
                                        minB = maxB-DB_noise;
                                        axis xy;
                                        caxis('manual');
                                        caxis([minB maxB]);
                                        cmap = spec_cmap();
                                        colormap(cmap);
                                        v_axis = axis;
                                        v_axis(3)=0;
                                        v_axis(4)=FHigh_spec_Logger + 5000;
                                        axis(v_axis);
                                        xlabel('time (ms)'), ylabel('Frequency');
                                    end
%                                     yyaxis right
                                    %                                 ylabel('Logger ID')
                                    %                                 set(gca, 'YTick', 1:ll, 'YTickLabel', Fns_AL, 'YLim', [0 (length(AudioLogs)+1)],'YDir', 'reverse')
                                    %                             end
                                    figure(3);
                                    subplot(2,1,1)
                                    hold on
                                    yyaxis right
                                    plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [ll ll], 'k:', 'LineWidth',2)
                                    ylabel('Microphone')
                                    set(gca, 'YTickLabel', [])
                                    hold off
                                    
                                    subplot(2,1,2)
                                    hold on
                                    yyaxis right
                                    plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [5500 5500] , 'k:', 'LineWidth',2)
                                    hold on
                                    plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [200 200] , 'k', 'LineWidth',2)
                                    ylabel(Fns_AL{ll})
                                    ylim(v_axis(3:4))
                                    set(gca, 'YTickLabel', [])
                                    hold off
                                    
                                    %                 IndVocStartRaw{ll}(ii) = round(IndVocStart{ll}(ii)/Fs_env*FS);
                                    %                 IndVocStopRaw{ll}(ii) = round(IndVocStop{ll}(ii)/Fs_env*FS);
                                    %                 IndVocStartPiezo{ll}(ii) = round(IndVocStart{ll}(ii)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                                    %                 IndVocStopPiezo{ll}(ii) = round(IndVocStop{ll}(ii)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                                    
                                    % Cross-correlate each sound extract with the raw trace and
                                    % check if it has a value of cross correlation high enough
                                    % to be kept as a vocalization  %% NOT Working :-(
                                    %                 Resamp_Filt_Logger_wav = resample(LowPassLogVoc{vv}{ll}(IndVocStartPiezo{ll}(ii):IndVocStopPiezo{ll}(ii)), 4*BandPassFilter(2),round(Piezo_FS.(Fns_AL{ll})(vv)));
                                    %                 Resamp_Filt_Raw_wav = resample(LowPassMicVoc{vv}(IndVocStartRaw{ll}(ii):IndVocStopRaw{ll}(ii)), 4*BandPassFilter(2),FS);
                                    %                 [Xcor_local,~] = xcorr(Resamp_Filt_Raw_wav,Resamp_Filt_Logger_wav, (Buffer*2)*4*BandPassFilter(2), 'unbiased');
                                    %                 Xcor{ll}(ii) = max(Xcor_local)/(var(Resamp_Filt_Raw_wav)*var(Resamp_Filt_Logger_wav))^0.5;
                                    
                                    % Calculate the Average RMS Ratio for each sound extract
                                    % and decide about the vocalization ownership
                                    RMSRatio{ll}(ii) = mean(RatioAmp(ll,IndVocStart{ll}(ii):IndVocStop{ll}(ii)));
                                    RMSDiff{ll}(ii) = mean(DiffAmp(ll,IndVocStart{ll}(ii):IndVocStop{ll}(ii)));
                                    %                     Call1Hear0_temp(ii) = RMSRatio{ll}(ii) > Factor_AmpRatio * RatioRMS.(Fns_AL{ll})(1);
                                    Call1Hear0_temp(ii) = RMSDiff{ll}(ii) > Factor_AmpDiff * DiffRMS.(Fns_AL{ll})(1);
                                    
                                    % update figure(3) with the decision
                                    figure(3)
                                    subplot(2,1,1)
                                    yyaxis right
                                    hold on
                                    if Call1Hear0_temp(ii)
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, ll+0.4, 'Call')
                                    else
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, ll+0.5, 'Hear')
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, ll+0.4, '/Noise')
                                    end
                                    hold off
                                    
                                    subplot(2,1,2)
                                    yyaxis right
                                    hold on
                                    if Call1Hear0_temp(ii)
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, 8000, 'Call')%ll+0.4
                                    else
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, 8500, 'Hear')%ll+0.5
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, 8000, '/Noise')%ll+0.4
                                    end
                                    hold off
                                    suplabel(sprintf('Ambient Microphone and %s Voc %d/%d Set %d/%d',Fns_AL{ll},vv,Nvoc,df,length(DataFiles)), 't');
                                    
                                    
                                    if Call1Hear0_temp(ii)
                                        fprintf('Computer guess for that sound element: %s calling\n',Fns_AL{ll});
                                    else
                                        fprintf('Computer guess for that sound element: %s hearing/noise\n',Fns_AL{ll});
                                    end
                                    if ManCall &&  (ManCall_logger==1)
                                        TempIn = input('Indicate your choice: calling (1);    hearing/noise (0);    Listen again to that logger recording (any other number)\n');
                                         
                                        if isempty(TempIn)
                                            fprintf(1, 'NO ENTRY GIVEN, playing again the sound and asking the same question\n')
                                            Call1Hear0_man(ii)=2;
                                        else
                                            Call1Hear0_man(ii) = TempIn;
                                        end
                                        while Call1Hear0_man(ii)~=0 && Call1Hear0_man(ii)~=1
                                            Player= audioplayer(LoggerSound/(VolDenominatorLogger*std(LoggerSound, 'omitnan')), Piezo_FS.(Fns_AL{ll})(vv)); %#ok<TNMLP>
                                            play(Player)
                                            pause(length(Raw_wave_nn)/FS +1)
                                            TempIn = input('Indicate your choice: calling (1);    hearing/noise (0);    listen again to that logger recording (any other number)\n');
                                            if isempty(TempIn)
                                                fprintf(1, 'NO ENTRY GIVEN, playing again the sound and asking the same question\n')
                                                Call1Hear0_man(ii)=2;
                                            else
                                                Call1Hear0_man(ii) = TempIn;
                                            end
                                        end
                                    elseif ManCall &&  (ManCall_logger==101)
                                        LoggerSound_ext = LoggerSound(floor(IndVocStart{ll}(ii)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv)) : min(ceil(IndVocStop{ll}(ii)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv)), length(LoggerSound)));
                                        Player= audioplayer(LoggerSound_ext/(VolDenominatorLogger*std(LoggerSound_ext, 'omitnan')), Piezo_FS.(Fns_AL{ll})(vv)); %#ok<TNMLP>
                                        play(Player)
                                        pause(length(LoggerSound_ext)/Piezo_FS.(Fns_AL{ll})(vv) +1)
                                        TempIn = input('Indicate your choice: calling (1);    hearing/noise (0);    Listen again to that logger extract recording (any other number)\n');
                                         
                                        if isempty(TempIn)
                                            fprintf(1, 'NO ENTRY GIVEN, playing again the sound and asking the same question\n')
                                            Call1Hear0_man(ii)=2;
                                        else
                                            Call1Hear0_man(ii) = TempIn;
                                        end
                                        while Call1Hear0_man(ii)~=0 && Call1Hear0_man(ii)~=1
                                            play(Player)
                                            pause(length(LoggerSound_ext)/Piezo_FS.(Fns_AL{ll})(vv) +1)
                                            TempIn = input('Indicate your choice: calling (1);    hearing/noise (0);    listen again to that logger extract recording (any other number)\n');
                                            if isempty(TempIn)
                                                fprintf(1, 'NO ENTRY GIVEN, playing again the sound and asking the same question\n')
                                                Call1Hear0_man(ii)=2;
                                            else
                                                Call1Hear0_man(ii) = TempIn;
                                            end
                                        end
                                    elseif ManCall_logger==100
                                        fprintf(1,'Manual input enforced: Take computer guess (100)\n');
                                        Call1Hear0_man(ii)=Call1Hear0_temp(ii);
                                    else
                                        fprintf(1,'Manual input enforced: Noise (0)\n');
                                    end
                                    Agree = Call1Hear0_temp(ii)== Call1Hear0_man(ii);
                                    if ~Agree
                                        Call1Hear0_temp(ii) = Call1Hear0_man(ii);
                                        PiezoError = PiezoError + [1 1];
                                        PiezoErrorType = PiezoErrorType + [Call1Hear0_temp(ii) ~Call1Hear0_temp(ii)];
                                    else
                                        PiezoError = PiezoError + [0 1];
                                    end
                                    
                                    
                                    % update figure(3) with the decision
                                    figure(3)
                                    subplot(2,1,1)
                                    yyaxis right
                                    hold on
                                    if Call1Hear0_temp(ii)
                                        plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, ll+0.2, 'Call','Color','r','FontWeight','bold')
                                    else
                                        plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [ll ll], 'k-', 'LineWidth',2)
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, ll+0.2, 'Hear','Color','r','FontWeight','bold')
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, ll+0.1, '/Noise','Color','r','FontWeight','bold')
                                    end
                                    hold off
                                    
                                    subplot(2,1,2)
                                    yyaxis right
                                    hold on
                                    if Call1Hear0_temp(ii)
                                        plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [5500 5500], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
                                        hold on
                                        plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [200 200], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, 7500, 'Call','Color','r','FontWeight','bold')%ll+0.2
                                    else
                                        plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [200 200], 'k-', 'LineWidth',2)
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, 7500, 'Hear','Color','r','FontWeight','bold')%ll+0.2
                                        text(IndVocStart{ll}(ii)/Fs_env*1000, 7000, '/Noise','Color','r','FontWeight','bold')%ll+0.1
                                    end
                                    hold off
                                    suplabel(sprintf('Ambient Microphone and %s Voc %d/%d Set %d/%d',Fns_AL{ll},vv,Nvoc,df,length(DataFiles)), 't');
                                    
                                    
                                    %                 figure(1)
                                    %                 subplot(length(AudioLogs)+2,1,ll+1)
                                    %                 hold on
                                    %                 yyaxis right
                                    %                 YLim = get(gca, 'YLim');
                                    %                 if Call1Hear0_temp{ll}(ii)
                                    %                     plot([IndVocStart{ll}(ii) IndVocStop{ll}(ii)]/Fs_env*1000, [YLim(2)*4/5 YLim(2)*4/5], 'k-', 'LineWidth',2)
                                    %                     subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
                                    %                     hold on
                                    %                     plot([IndVocStart{ll}(ii) IndVocStop{ll}(ii)]/Fs_env*1000, [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
                                    %                     hold off
                                    %                 else
                                    %                     plot([IndVocStart{ll}(ii) IndVocStop{ll}(ii)]/Fs_env*1000, [YLim(2)*4/5 YLim(2)*4/5], 'LineStyle','-','Color',[0.8 0.8 0.8], 'LineWidth',2)
                                    %                 end
                                    %                 hold off
                                    
                                end
                                
                                % Stash noise and Only keep vocalizations that were produced by the bat
                                % according to the high RMSRatio value
                                %                             IndNoiseStart{ll} = IndVocStart{ll}(~logical(Call1Hear0_temp));
                                %                             IndNoiseStop{ll} = IndVocStop{ll}(~logical(Call1Hear0_temp));
                                %                             IndNoiseStartRaw{vv}{ll} = round(IndNoiseStart{ll}/Fs_env*FS);
                                %                             IndNoiseStartPiezo{vv}{ll} = round(IndNoiseStart{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                                %                             IndNoiseStopRaw{vv}{ll} = round(IndNoiseStop{ll}/Fs_env*FS);
                                %                             IndNoiseStopPiezo{vv}{ll} = round(IndNoiseStop{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                                IndVocStart{ll} = IndVocStart{ll}(logical(Call1Hear0_temp));
                                IndVocStop{ll} = IndVocStop{ll}(logical(Call1Hear0_temp));
                                IndVocStartRaw{vv}{ll} = round(IndVocStart{ll}/Fs_env*FS);
                                IndVocStartPiezo{vv}{ll} = round(IndVocStart{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                                IndVocStopRaw{vv}{ll} = round(IndVocStop{ll}/Fs_env*FS);
                                IndVocStopPiezo{vv}{ll} = round(IndVocStop{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                                NV = sum(Call1Hear0_temp);
                                
                                if NV
                                    % Now merge detected vocalizations that are within MergeThresh
                                    Merge01 = (IndVocStart{ll}(2:end)-IndVocStop{ll}(1:(end-1)) < (MergeThresh/1000*Fs_env));
                                    NV = NV-sum(Merge01);
                                    IndVocStartRaw_merge_local{ll} = nan(1,NV);
                                    IndVocStopRaw_merge_local{ll} = nan(1,NV);
                                    IndVocStartPiezo_merge_local{ll} = nan(1,NV);
                                    IndVocStopPiezo_merge_local{ll} = nan(1,NV);
                                    CutInd = find(~Merge01);
                                    IndVocStartRaw_merge_local{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*FS);
                                    IndVocStartPiezo_merge_local{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                                    IndVocStopRaw_merge_local{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*FS);
                                    IndVocStopPiezo_merge_local{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                                    for cc=1:length(CutInd)
                                        IndVocStopRaw_merge_local{ll}(cc) = round(IndVocStop{ll}( CutInd(cc))/Fs_env*FS);
                                        IndVocStopPiezo_merge_local{ll}(cc) = round(IndVocStop{ll}(CutInd(cc))/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                                        IndVocStartRaw_merge_local{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*FS);
                                        IndVocStartPiezo_merge_local{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                                    end
                                    
                                    % Now plot the onset/offset of each extract on the
                                    % spectrograms of the loggers
                                    figure(1)
                                    for ii=1:NV
                                        subplot(length(AudioLogs)+2,1,ll+1)
                                        hold on
                                        yyaxis right
                                        YLim = get(gca, 'YLim');
                                        plot([IndVocStartRaw_merge_local{ll}(ii) IndVocStopRaw_merge_local{ll}(ii)]/FS*1000, [YLim(2)*4/5 YLim(2)*4/5], 'k-', 'LineWidth',2)
                                        subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
                                        hold on
                                        plot([IndVocStartRaw_merge_local{ll}(ii) IndVocStopRaw_merge_local{ll}(ii)]/FS*1000, [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
                                        hold off
                                    end
                                    
                                    %             subplot(length(AudioLogs)+2,1,ll+1)
                                    %             hold on
                                    %             yyaxis right
                                    %             YLim = get(gca, 'YLim');
                                    %             for ii=1:NV
                                    %                 if Call1Hear0{ll}(ii)
                                    %                     plot(round([IndVocStartRaw{ll}(ii) IndVocStopRaw{ll}(ii)]/FS*1000), [YLim(2)*4.5/5 YLim(2)*4.5/5], 'r-', 'LineWidth',2)
                                    %                 else
                                    %                     plot(round([IndVocStartRaw{ll}(ii) IndVocStopRaw{ll}(ii)]/FS*1000), [YLim(2)*4.5/5 YLim(2)*4.5/5], '-', 'LineWidth',2, 'Color',[1 0.5 0.5])
                                    %                 end
                                    %             end
                                    %             hold off
                                end
                            else
                                
                            end
                        end
                    end
                end
                figure(1)
                subplot(length(AudioLogs)+2,1,1)
                XLIM = get(gca, 'XLim');
                subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
                set(gca, 'YTick', 1:length(AudioLogs))
                set(gca, 'YTickLabel', [Fns_AL; 'Mic'])
                set(gca, 'YLim', [0 length(AudioLogs)+1])
                set(gca, 'XLim', XLIM);
                ylabel('AL ID')
                xlabel('Time (ms)')
                if ~isempty(strfind(VocFilename{vv}, filesep))
                    [~,FileVoc]=fileparts(VocFilename{vv});
                elseif ~isempty(strfind(VocFilename{vv},'\'))
                    AllSep = strfind(VocFilename{vv},'\');
                    FileVoc = VocFilename{vv}((AllSep(end)+1):(strfind(VocFilename{vv}, '.')-1));
                elseif ~isempty(strfind(VocFilename{vv},'/'))
                    AllSep = strfind(VocFilename{vv},'/');
                    FileVoc = VocFilename{vv}((AllSep(end)+1):(strfind(VocFilename{vv}, '.')-1));
                end
            
                if strcmp(SaveFileType,'pdf')
                    if sum(cellfun(@length, IndVocStartRaw_merge_local)) % Only save the RMS and spectro figures if there was a vocalization
                        fprintf(1,'saving figures...\n')
                        try
                            print(F1,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_spec_%d.pdf',FileVoc,vv, df,MergeThresh)),'-dpdf','-fillpage')
                        catch
                            warning('Issues with saving that figure in pdf! too big?! saving as fig')
                            try 
                                saveas(F1,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_spec_%d.fig',FileVoc,vv, df,MergeThresh)))
                            catch
                                warning('cannot even save as fig...')
                            end
                        end
                        if RMS_Fig
                            saveas(F2,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_RMS_%d.pdf',FileVoc,vv, df,MergeThresh)),'pdf') %#ok<UNRCH> 
                        end
                    end
                elseif strcmp(SaveFileType,'fig')
                    if sum(cellfun(@length, IndVocStartRaw_merge_local)) % Only save the RMS and spectro figures if there was a vocalization
                        fprintf(1,'saving figures...\n')
                        saveas(F1,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_spec_%d.fig', FileVoc,vv, df,MergeThresh)))
                        if RMS_Fig
                            saveas(F2,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_RMS_%d.fig', FileVoc,vv, df,MergeThresh))) %#ok<UNRCH> 
                        end
                    end
                end
                if ManCall
                    pause(1)
                end
                clf(F1)
                if RMS_Fig
                    clf(F2) %#ok<UNRCH> 
                end
                
                % Gather Vocalization production data:
                IndVocStart_all{vv} = IndVocStart;
                IndVocStop_all{vv} = IndVocStop;
%                 IndNoiseStart_all{vv} = IndNoiseStart;
%                 IndNoiseStop_all{vv} = IndNoiseStop;
                IndVocStartRaw_merged{vv} = IndVocStartRaw_merge_local;
                IndVocStopRaw_merged{vv} = IndVocStopRaw_merge_local;
                IndVocStartPiezo_merged{vv} = IndVocStartPiezo_merge_local;
                IndVocStopPiezo_merged{vv} = IndVocStopPiezo_merge_local;
                RMSRatio_all{vv} = RMSRatio;
                RMSDiff_all{vv} = RMSDiff;
                SorterName{vv} = Sorter;
            end
            fprintf(1,'saving data...\n')
            if vvi==length(List2WorkOn) % this is the last vocalization
                vv=Nvoc;
            end
            if ~isempty(dir(PreviousFile))
%                 save(fullfile(Working_dir, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','IndNoiseStart_all','IndNoiseStop_all', 'IndNoiseStartRaw', 'IndNoiseStopRaw', 'IndNoiseStartPiezo', 'IndNoiseStopPiezo','RMSRatio_all','RMSDiff_all','vv','MicError','PiezoError','MicErrorType','PiezoErrorType','-append');
                  save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','vv','MicError','PiezoError','MicErrorType','PiezoErrorType','SorterName','-append');
            else
%                 save(fullfile(Working_dir, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','IndNoiseStart_all','IndNoiseStop_all', 'IndNoiseStartRaw', 'IndNoiseStopRaw', 'IndNoiseStartPiezo', 'IndNoiseStopPiezo','RMSRatio_all','RMSDiff_all','vv','MicError','PiezoError','MicErrorType','PiezoErrorType');
                save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','vv','SorterName','MicError','PiezoError','MicErrorType','PiezoErrorType');
            end
            if SaveRawWave
                save(DataFile, 'Raw_wave','-append')
            end
%             if SaveRawWaveName
%                 if Chunking_RawWav
%                     warning('Highly not recommended!! We are chuncking the loading, you might loose the entire Raw_wave data here')
%                     keyboard
%                 end
%                 save(DataFile, 'VocFilename','-append')
%                 save(fullfile(Raw_dir, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)), 'Voc_filename','-append')
%             end
        end
        if SaveRawWave && ~strcmp(Working_dir_read,DataFiles(df).folder)
            [s2,~]=copyfile(DataFile, fullfile(DataFiles(df).folder,DataFiles(df).name), 'f');
            if ~s2
                fprintf(1,'File transfer of %s to %s did not occur correctly\n',DataFile,fullfile(DataFiles(df).folder,DataFiles(df).name))
                keyboard
            end
            if s2  %erase local data
                [~]=rmdir(DataFile, 's');
            end
        end
        clear List2WorkOnResNet
    end
    if ~strcmp(Working_dir_write,Loggers_dir)
        fprintf(1,'Transferring data back on the server\n')
        [s1,~]=copyfile(fullfile(Working_dir_write,'*'), Loggers_dir, 'f');
        if ~s1
            fprintf(1,'File transfer did not occur correctly\n')
            keyboard
        end
        if s1  %erase local data
            [~]=rmdir(Working_dir_write, 's');
        end
        
    end
   
end
end