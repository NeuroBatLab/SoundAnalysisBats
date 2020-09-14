function [SaveRawWaveAll] = correctMicAllignment_beforeorafterWhoCalls(Raw_dir, Loggers_dir, Date, ExpStartTime, varargin)
VolFactorMic=0.5;
pnames = {'Working_dir'};
dflts  = {Loggers_dir};
[Working_dir] = internal.stats.parseArgs(pnames,dflts,varargin{:});

Manual=0; % set to 1 to listen to Mic sequences
PlotSpec = 0; % set to 1 to plot spectrograms


Working_dir_read = fullfile(Working_dir, 'read');
Working_dir_write = fullfile(Working_dir, 'write');
TTL_dir = dir(fullfile(Raw_dir,sprintf( '%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
TTL = load(fullfile(TTL_dir.folder, TTL_dir.name));
FileNum_u = unique(TTL.File_number);

if ~strcmp(Loggers_dir,Working_dir) && (~exist(Working_dir,'dir') || ~exist(Working_dir_read,'dir') || ~exist(Working_dir_write,'dir'))
    mkdir(Working_dir)
    mkdir(Working_dir_read)
    mkdir(Working_dir_write)
elseif strcmp(Loggers_dir,Working_dir)
    Working_dir_read = Loggers_dir;
    Working_dir_write = Loggers_dir;
end
DataFiles = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData*.mat', Date, ExpStartTime)));

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
    DataFiles = DataFiles(AscendOrd);
    load(fullfile(Raw_dir, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)), 'Voc_filename')
    Nvoc_all = length(Voc_filename);
    DataFile = fullfile(DataFiles(1).folder, DataFiles(1).name);
    load(DataFile, 'VocMaxNum')
    if ~exist('VocMaxNum','var')
        VocMaxNum=1000;
    end
    if (Nvoc_all>VocMaxNum) && (length(DataFiles)>1)
        Nvocs = [0 VocMaxNum:VocMaxNum:Nvoc_all (floor(Nvoc_all/VocMaxNum)*VocMaxNum+rem(Nvoc_all,VocMaxNum))];
    else
        Nvocs = [0 Nvoc_all];
    end
    SaveRawWaveAll=zeros(length(DataFiles),1);
    for df=1:length(DataFiles)
        Nvoc = Nvocs(df+1) - Nvocs(df);
        DataFile = fullfile(DataFiles(df).folder, DataFiles(df).name);
        vv=1;
        if ~strcmp(Working_dir_write,Loggers_dir) && ~isfile(fullfile(Working_dir_read,DataFiles(df).name))
            fprintf(1,'Bringing data locally from the server\n')
            [s,m,e]=copyfile(DataFile, Working_dir_read, 'f');
            if ~s
                fprintf(1,'File transfer did not occur correctly\n')
                keyboard
            else
                fprintf(1,'File transfer DONE\n')
                DataFile = fullfile(Working_dir_read,DataFiles(df).name);
            end
        elseif ~strcmp(Working_dir_write,Loggers_dir) && isfile(fullfile(Working_dir_read,DataFiles(df).name))
                fprintf(1,'File already transferred from previous code running\n')
                DataFile = fullfile(Working_dir_read,DataFiles(df).name);
        end
        load(DataFile,'Raw_wave')
        if Nvoc ~= length(Raw_wave)
            warning('Looks like there might be an issue there!! Check variables!!')
            keyboard
        end
        if Nvoc<=100
            RawWavChunking = 0;
            minvv = 1;
            maxvv = Nvoc;
            load(DataFile,'Piezo_wave', 'AudioLogs',   'Piezo_FS',  'FS','VocFilename','Voc_transc_time_refined');
        else % often problem of memory, we're going to chunck file loading
            RawWavChunking = 1;
            if ~mod(vv,100)
                minvv=floor((vv-1)/100)*100 +1;
                maxvv=ceil(vv/100)*100;
            else
                minvv = floor(vv/100)*100 +1;
                maxvv = ceil(vv/100)*100;
                if maxvv<minvv
                    maxvv = ceil((vv+1)/100)*100;
                end
            end
            Raw_wave = Raw_wave(minvv:min(maxvv, length(Raw_wave)));
            load(DataFile,'Piezo_wave', 'AudioLogs',   'Piezo_FS',  'FS', 'VocFilename','Voc_transc_time_refined');
        end
        
        Fns_AL = fieldnames(Piezo_wave);
        
        % parameters
        DB_noise = 60; % Noise threshold for the spectrogram colormap
        FHigh_spec = 90000; % Max frequency (Hz) for the raw data spectrogram
        FHigh_spec_Logger = 10000; % Max frequency (Hz) for the raw data spectrogram
        BandPassFilter = [1000 5000 9900]; % Frequency bands chosen for digital signal processing
        
        
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
        
 
        
        %% Loop through vocalizations sequences and check microphone data
        for vv=vv:Nvoc
            fprintf(1,'\n\n\n\nVoc sequence %d/%d Set %d/%d\n',vv,Nvoc, df, length(DataFiles));
            
            % Patch for previous error in the code
            if vv<=maxvv
                Raw_wave_nn = Raw_wave{vv - (minvv -1)};
            else
                RawWavChunking = 1;
                clear Raw_wave Piezo_wave
                load(DataFile,'Raw_wave')
                if ~mod(vv,100)
                    minvv=floor((vv-1)/100)*100 +1;
                    maxvv=ceil(vv/100)*100;
                else
                    minvv = floor(vv/100)*100 +1;
                    maxvv = ceil(vv/100)*100;
                    if maxvv<minvv
                        maxvv = ceil((vv+1)/100)*100;
                    end
                end
                Raw_wave = Raw_wave(minvv:min(maxvv, length(Raw_wave)));
                load(DataFile,'Piezo_wave')
                Raw_wave_nn = Raw_wave{vv - (minvv -1)};
            end
             
            
            % retrieving file name index of the microphone
            Voc_i_start = Nvocs(df)+1;
            vv_in = vv + Voc_i_start-1;
            if ~strcmp(Voc_filename{vv_in}, VocFilename{vv})
                warning('Issues with Mic file name\n')
                SaveRawWaveName = 1;
                keyboard
            else
                SaveRawWaveName = 0;
            end
            fprintf(1, 'Microphone File: %s\n', Voc_filename{vv_in})
%             if Nvoc>100
%                 warning('Probably wrong audio file name, the code is not updated for older version of previous extraction\n')
%                 keyboard
%             end
            if isempty(Raw_wave_nn)
                warning('Raw_wave should not be empty now!!')
                keyboard
                SaveRawWave = 1;
                [Raw_wave{vv - (minvv -1)}, FS] = audioread(VocFilename{vv});
                Raw_wave_nn = Raw_wave{vv - (minvv -1)};
            else
                SaveRawWave = 0;
            end
            
            % Check that correct microphone file was saved (Trying to
            % detect/fix bug from voc_localize_using_piezo
            OnOffTranscTime_ms = Voc_transc_time_refined(vv,:);
            FileNumIdx = find(TTL.Pulse_TimeStamp_Transc<OnOffTranscTime_ms(1,1),1,'Last');
            if isempty(FileNumIdx)
                FileNumIdx = find(TTL.Pulse_TimeStamp_Transc>OnOffTranscTime_ms(1,1),1,'First');
            end
            MicVoc_File = TTL.File_number(FileNumIdx);
            IndFileNum = find(FileNum_u == MicVoc_File);
            TranscTime_zs = (OnOffTranscTime_ms - TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,1))/TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,2);
            MicVoc_samp_idx =round(TTL.Mean_std_Pulse_samp_audio(IndFileNum,2) .* polyval(TTL.Slope_and_intercept_transc2audiosamp{IndFileNum},TranscTime_zs,[],TTL.Mean_std_x_transc2audiosamp{IndFileNum}) + TTL.Mean_std_Pulse_samp_audio(IndFileNum,1));
            WavFileStruc_local = dir(fullfile(Raw_dir, sprintf('*_%s_%s*mic*_%d.wav',Date, ExpStartTime, MicVoc_File)));
            Raw_filename = fullfile(WavFileStruc_local.folder, WavFileStruc_local.name);
            [Raw_10minwav2, FS2] = audioread(Raw_filename);
            if MicVoc_samp_idx(1)>length(Raw_10minwav2) % This vocalization occured in the next file
                MicVoc_File = MicVoc_File+1;
                IndFileNum = find(FileNum_u == MicVoc_File);
                TranscTime_zs = (OnOffTranscTime_ms - TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,1))/TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,2);
                MicVoc_samp_idx =round(TTL.Mean_std_Pulse_samp_audio(IndFileNum,2) .* polyval(TTL.Slope_and_intercept_transc2audiosamp{IndFileNum},TranscTime_zs,[],TTL.Mean_std_x_transc2audiosamp{IndFileNum}) + TTL.Mean_std_Pulse_samp_audio(IndFileNum,1));
                WavFileStruc_local = dir(fullfile(Raw_dir, sprintf('*_%s_%s*mic*_%d.wav',Date, ExpStartTime, MicVoc_File)));
                Raw_filename = fullfile(WavFileStruc_local.folder, WavFileStruc_local.name);
                [Raw_10minwav2, FS2] = audioread(Raw_filename);
            end
            Raw_wave_ex = Raw_10minwav2(MicVoc_samp_idx(1) : min(MicVoc_samp_idx(2),length(Raw_10minwav2)));
            if length(Raw_wave_ex)<length(Raw_wave_nn)
                Corr(1) = corr(Raw_wave_ex,Raw_wave_nn(1:length(Raw_wave_ex)));
                Corr(2) = corr(Raw_wave_ex,Raw_wave_nn(end-length(Raw_wave_ex)+1:end));
                Corr(3) = corr(Raw_wave_ex,Raw_wave{vv- (minvv -1)}(1:length(Raw_wave_ex)));
                Corr(4) = corr(Raw_wave_ex,Raw_wave{vv- (minvv -1)}(end-length(Raw_wave_ex)+1:end));
            elseif length(Raw_wave_ex)>length(Raw_wave_nn)
                Corr(1) = corr(Raw_wave_nn,Raw_wave_ex(1:length(Raw_wave_nn)));
                Corr(2) = corr(Raw_wave_nn,Raw_wave_ex(end-length(Raw_wave_nn)+1:end));
                Corr(3) = corr(Raw_wave{vv- (minvv -1)},Raw_wave_ex(1:length(Raw_wave_nn)));
                Corr(4) = corr(Raw_wave{vv- (minvv -1)},Raw_wave_ex(end-length(Raw_wave_nn)+1:end));
            elseif length(Raw_wave_ex)==length(Raw_wave_nn)
                Corr(1) = corr(Raw_wave_ex,Raw_wave_nn);
                Corr(2) = corr(Raw_wave_ex,Raw_wave{vv- (minvv -1)});
            end
            if all(Corr<0.99)
                warning('Error in the microphone file that was previosuly saved, fixing the issue now!\n')
%                 keyboard
                SaveRawWave = 1;
                Raw_wave_nn = Raw_wave_ex;
                Raw_wave{vv- (minvv -1)} = Raw_wave_ex;
                TrueVocName = fullfile(Raw_dir, 'Detected_calls',sprintf('%s_%s_%s_voc_%d_%d.wav',WavFileStruc_local.name(1:4),Date,ExpStartTime, MicVoc_File, MicVoc_samp_idx(1)));
                if ~strcmp(VocFilename{vv}, TrueVocName)
                    SaveRawWaveName = 1;
                    warning('Filename was also wrong correcting %s -> %s\n',VocFilename{vv},TrueVocName)
                    % delete old file
                    delete(VocFilename{vv})
                    VocFilename{vv}= TrueVocName;
                    Voc_filename{vv_in} = TrueVocName;
                end
                % overwrite what was wrongly saved at the time
                audiowrite(VocFilename{vv} , Raw_wave_ex, FS2);
            end

            if PlotSpec
                % bandpass filter the ambient mic recording
                Filt_RawVoc = filtfilt(sos_raw_band,1,Raw_wave_nn);

                % Plot the spectrogram of the ambient microphone
                F1=figure(1);
                clf(F1)
                subplot(length(AudioLogs)+2,1,1)
                [Raw_Spec.to, Raw_Spec.fo, Raw_Spec.logB] = spec_only_bats(Filt_RawVoc, FS, DB_noise, FHigh_spec);
                hold on
                title(sprintf('Voc %d/%d Set %d/%d',vv,Nvoc, df,length(DataFiles)))
                xlabel(' ')
                set(gca, 'XTick',[],'XTickLabel',{})
            end
        
            if Manual
                pause(0.1)
                Raw_listen = filtfilt(sos_raw_band_listen,1,Raw_wave_nn);
                SampleMic = resample((Raw_listen - mean(Raw_listen))/(std(Raw_listen)/VolFactorMic),FS/4,FS);
                PlayerMic= audioplayer(SampleMic, FS/4,24); %#ok<TNMLP>
                play(PlayerMic)
                pause(length(Raw_wave_nn)/FS +1)
            end
            
            if PlotSpec
                %% Loop through the loggers and check the extracts length
                LengthLoggersData = nan(length(AudioLogs),1);
                for ll=1:length(AudioLogs)
                    LengthLoggersData(ll) = length(Piezo_wave.(Fns_AL{ll}){vv});
                end
                %% Loop through the loggers and plot spectro
                Logger_Spec = cell(length(AudioLogs),1);
                for ll=1:length(AudioLogs)
                    if isnan(Piezo_FS.(Fns_AL{ll})(vv)) || isempty(Piezo_wave.(Fns_AL{ll}){vv})
                        fprintf(1, 'NO DATA for Vocalization %d from %s\n', vv, Fns_AL{ll})
                    else
                        % design the filters
                        [z,p,k] = butter(6,BandPassFilter(1:2)/(Piezo_FS.(Fns_AL{ll})(vv)/2),'bandpass');
                        sos_low = zp2sos(z,p,k);

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
                                LowPassLogVoc = (filtfilt(sos_low,1,InputPiezo)); % low-pass filter the voltage trace
                            else
                                if length(InputPiezo)>min(LengthLoggersData)
                                    warning('Piezo data have different durations!\n Logger %s data is truncated for analysis\n',Fns_AL{ll})
                                    InputPiezo = InputPiezo(1:min(LengthLoggersData));
                                end
                                LowPassLogVoc = (filtfilt(sos_low,1,InputPiezo)); % low-pass filter the voltage trace
                            end

                            % Plot the low pass filtered signal of each logger
                            figure(1)
                            subplot(length(AudioLogs)+2,1,ll+1)
                            [Logger_Spec{ll}.to, Logger_Spec{ll}.fo, Logger_Spec{ll}.logB] = spec_only_bats(LowPassLogVoc, Piezo_FS.(Fns_AL{ll})(vv), DB_noise, FHigh_spec_Logger);
                            if ll<length(AudioLogs)
                                xlabel(' '); % supress the x label output
                                set(gca,'XTick',[],'XTickLabel',{});
                            end

                        end
                    end
                end
            end
            
            if SaveRawWave
                if RawWavChunking
                    warning('Touchy save as we chuncked Rawwave, you might want to check\n')
%                     keyboard
                    Raw_wave_local = Raw_wave;
                    load(DataFile,'Raw_wave')
                    Raw_wave(minvv:min(maxvv, length(Raw_wave))) = Raw_wave_local;
                    save(DataFile, 'Raw_wave','-append')
                    SaveRawWaveAll(df) = 1;
                    Raw_wave = Raw_wave_local;
                else
                    save(DataFile, 'Raw_wave','-append')
                    SaveRawWaveAll(df) = 1;
                end
            end
            if SaveRawWaveName
                save(DataFile, 'VocFilename','-append')
                save(fullfile(Raw_dir, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)), 'Voc_filename','-append')
            end
        end
        
        if SaveRawWaveAll(df) && ~strcmp(Working_dir_read,DataFiles(df).folder)
            [s2,m,e]=copyfile(DataFile, fullfile(DataFiles(df).folder,DataFiles(df).name), 'f');
            if ~s2
                fprintf(1,'File transfer of %s to %s did not occur correctly\n',DataFile,fullfile(DataFiles(df).folder,DataFiles(df).name))
                keyboard
            else
                fprintf(1,'File transfer of %s to %s DONE\n',DataFile,fullfile(DataFiles(df).folder,DataFiles(df).name))
            end
            
        end
        if ~strcmp(Working_dir_read,DataFiles(df).folder)
            fprintf(1, 'Erasing local data')
            delete(DataFile);
        end
    end
end



end