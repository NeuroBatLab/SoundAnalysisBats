function audioQuality_calls_Deafs(Loggers_dir, Date, ExpStartTime)

BandPassFilter = [1000 5000];

% Set to 1 if you want to manually pause after each vocalization
ManualPause=1;


% Load data
DataFiles = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractDat*_*.mat', Date, ExpStartTime)));

% Get filters ready:
for jj = 1:length(DataFiles)
    DataFile = DataFiles(jj);
    load(fullfile(DataFile.folder, DataFile.name), 'BioSoundCalls');
    ii=1;
    while isempty(BioSoundCalls{ii,1}) && ii<=length(BioSoundCalls(:,1))
        ii=ii+1;
    end
    if ~isempty(BioSoundCalls{ii,1})
        break
    end
end

[z,p,k] = butter(6,BandPassFilter(1:2)/(BioSoundCalls{ii,1}.samprate/2),'bandpass');% a 12th order Butterworth band-pass filter; the second input argument is normalized cut-off frequency (ie. normalized to the Nyquist frequency, which is half the sampling frequency, as required by MATLAB)
sos_Raw = zp2sos(z,p,k); % obtain the second order section (biquad) filter to use as input in filtfilt

FS_ll = round(BioSoundCalls{ii,2}.samprate);
[z,p,k] = butter(6,BandPassFilter(1:2)/(FS_ll/2),'bandpass');
sos_Piezo = zp2sos(z,p,k);
clear BioSoundCalls

% Loop through the datafiles
for df=1:length(DataFiles) %1
    fprintf(1, 'Set %d/%d\n', df, length(DataFiles))
    DataFile = DataFiles(df);
    load(fullfile(DataFile.folder, DataFile.name), 'AudioGood');
    if exist('AudioGood', 'var') && (isempty(AudioGood) || ~isempty(AudioGood{end}))
        clear AudioGood
        continue
        
    elseif exist('AudioGood', 'var') && isempty(AudioGood{end})
        GoAudioGood = input(sprintf('It looks like we should start from here because the last vocalization quality is empty in AudioGood\n There is however %d empty cells\n Resume audioGood:1 skip and check the next file:0\n', sum(cellfun('isempty',AudioGood))));
        if GoAudioGood
            Rangevv = find(cellfun('isempty',AudioGood));
            load(fullfile(DataFile.folder, DataFile.name), 'BioSoundCalls','BioSoundFilenames', 'RMS', 'Duration','CorrPiezoRaw','ManualCallType','ManualAnnotationOK')
            NVoc = size(BioSoundCalls,1);
        else
            clear AudioGood
            continue
        end
    else
        load(fullfile(DataFile.folder, DataFile.name), 'BioSoundCalls','BioSoundFilenames')
        if ~exist('BioSoundCalls', 'var')
            fprintf(1, 'No Calls for that set\n')
            BioSoundCalls = [];
            AudioGood = [];
            Duration = [];
            ManualCallType = [];
            ManualAnnotationOK = [];
            save(fullfile(DataFile.folder, DataFile.name), 'Duration', 'AudioGood','BioSoundCalls','ManualCallType','ManualAnnotationOK', '-append')
            clear AudioGood BioSoundCalls ManualCallType
            continue
        end
        NVoc = size(BioSoundCalls,1);
        CorrPiezoRaw = nan(NVoc,1);
        if ManualPause
            AudioGood = cell(NVoc,1);
            ManualCallType = cell(NVoc,1);
            ManualAnnotationOK = cell(NVoc,1);
        end
        Duration = nan(NVoc,1);
        RMS = nan(NVoc,1);
        Rangevv = 1:NVoc;
    end
    for jj=1:length(Rangevv)
        vv=Rangevv(jj);
        % Open the Biosound calculation results for the piezo
        open([BioSoundFilenames{vv,2}(1:end-4) '.pdf'])
        % filter the original microphone wavfile
        if ~isfield(BioSoundCalls{vv,1}, 'sound')
            continue
        end
        Filt_Raw_wav=(filtfilt(sos_Raw,1,BioSoundCalls{vv,1}.sound)); % band-pass filter the voltage traces
        if ManualPause
            F1=figure(1);
            subplot(2,1,1)
            plot(BioSoundCalls{vv,1}.sound, '-k')
            hold on
            plot(Filt_Raw_wav, '-r')
            set(gca, 'XLim', [0 length(BioSoundCalls{vv,1}.sound)])
            title('Environmental Mic filtering and resampling')
            legend('Raw voltage trace', sprintf('BandPass %d %d Hz', BandPassFilter(1:2)))
        end
        % filter the piezo wavfile that was filtered by whatcalls between
        % 100 and 10kHz
        FS_ll = round(BioSoundCalls{vv,2}.samprate);
        Filt_Logger_wav = filtfilt(sos_Piezo,1,BioSoundCalls{vv,2}.sound); % band-pass filter the piezo sound
        if ManualPause
            figure(1)
            subplot(2,1,2)
            plot(BioSoundCalls{vv,2}.sound, '-k')
            hold on
            plot(Filt_Logger_wav, '-r')
            set(gca, 'XLim', [0 length(BioSoundCalls{vv,2}.sound)])
            title('Logger filtering and resampling')
            legend('Raw voltage trace', sprintf('BandPass %d %d Hz', BandPassFilter(1:2)))
        end
        
        
        % resample the sounds so they are at the same sample frequency of 4
        % times the low pass filter value
        Resamp_Filt_Raw_wav = resample(Filt_Raw_wav, 4*BandPassFilter(2), BioSoundCalls{vv,1}.samprate);
        Resamp_Filt_Logger_wav = resample(Filt_Logger_wav, 4*BandPassFilter(2),FS_ll);
        DiffLength = length(Resamp_Filt_Raw_wav) - length(Resamp_Filt_Logger_wav);
        if DiffLength == -1
            Resamp_Filt_Logger_wav = Resamp_Filt_Logger_wav(1:(end-1));
        elseif DiffLength ==1
            Resamp_Filt_Raw_wav = Resamp_Filt_Raw_wav(1:(end-1));
        elseif abs(DiffLength)>1
                warning('There was an isue with the extraction of microphone data, not as long as Logger data, skip vocalization')
                continue
        end
        
        if ManualPause
            figure(1)
            subplot(2,1,1)
            t2 = (0:(length(Resamp_Filt_Raw_wav)-1))*BioSoundCalls{vv,1}.samprate/(4*BandPassFilter(2));
            hold on
            plot(t2, Resamp_Filt_Raw_wav, 'g-', 'DisplayName',sprintf('BandPass + Resampled %d Hz', 4*BandPassFilter(2)))
            hold off
            subplot(2,1,2)
            t2 = (0:(length(Resamp_Filt_Logger_wav)-1))*FS_ll/(4*BandPassFilter(2));
            hold on
            plot(t2, Resamp_Filt_Logger_wav , 'g-', 'DisplayName',sprintf('BandPass + Resampled %d Hz', 4*BandPassFilter(2)))
            hold off
        end
        
        % Localize vocalizations on the microphone and logger and set to zero the
        % signal outside of putative vocalizations.
        % calculate a running RMS of the audio signal with a 1ms bin
        % resolution. a vocalization is defined as 10 consecutive time bins
        % (15ms) above the median amplitude here
        Consecutive_bins = 15;
        Fs_env = 1000; %Hz
        Vocp_logger = BioSoundCalls{vv,2}.amp>median(BioSoundCalls{vv,2}.amp);
        IndVocStart = strfind(Vocp_logger, ones(1,Consecutive_bins));
        if isempty(IndVocStart)
            keyboard
        end
        if isempty(IndVocStart)
            continue
        end
        IndVocStart_diffind = find(diff(IndVocStart)>1);
        IndVocStart = [IndVocStart(1) IndVocStart(IndVocStart_diffind +1)];
        NV = length(IndVocStart);
        IndVocStop = nan(NV,1);
        Clean_Resamp_Filt_Logger_wav = zeros(size(Resamp_Filt_Logger_wav));
        Clean_Resamp_Filt_Raw_wav = zeros(size(Resamp_Filt_Raw_wav));
        for ii=1:NV
            IVStop = find(Vocp_logger(IndVocStart(ii):end)==0, 1, 'first');
            if isempty(IVStop)
                IVStop = length(Vocp_logger(IndVocStart(ii):end));
            end
            IndVocStop(ii) = IndVocStart(ii) + IVStop;
            IndVocStart(ii) = round(IndVocStart(ii)/Fs_env*4*BandPassFilter(2));
            IndVocStop(ii) = round(IndVocStop(ii)/Fs_env*4*BandPassFilter(2));
            if IndVocStop(ii)>length(Clean_Resamp_Filt_Logger_wav) % This sound element reaches the end of the recording and downsampling messed up the exact indices values
                Clean_Resamp_Filt_Logger_wav(IndVocStart(ii):end) = Resamp_Filt_Logger_wav(IndVocStart(ii):end);
                Clean_Resamp_Filt_Raw_wav(IndVocStart(ii):end) = Resamp_Filt_Raw_wav(IndVocStart(ii):end);
            else
                Clean_Resamp_Filt_Logger_wav(IndVocStart(ii):IndVocStop(ii)) = Resamp_Filt_Logger_wav(IndVocStart(ii):IndVocStop(ii));
                Clean_Resamp_Filt_Raw_wav(IndVocStart(ii):IndVocStop(ii)) = Resamp_Filt_Raw_wav(IndVocStart(ii):IndVocStop(ii));
            end
        end
        if ManualPause
            figure(1)
            subplot(2,1,1)
            t2 = (0:(length(Resamp_Filt_Raw_wav)-1))*BioSoundCalls{vv,1}.samprate/(4*BandPassFilter(2));
            hold on
            plot(t2, Clean_Resamp_Filt_Raw_wav , 'c-', 'DisplayName','Cleaned signal for cross-correlation')
            hold off
            subplot(2,1,2)
            t2 = (0:(length(Resamp_Filt_Logger_wav)-1))*FS_ll/(4*BandPassFilter(2));
            hold on
            plot(t2, Clean_Resamp_Filt_Logger_wav , 'c-', 'DisplayName','Cleaned signal for cross-correlation')
            hold off
            
        end
        
        
        % correlate the raw data with the logger
        [XCorr, Lags]=xcorr(Clean_Resamp_Filt_Raw_wav',Clean_Resamp_Filt_Logger_wav',(10)*4*BandPassFilter(2)); % finding the optimal alignement between the 2 signals
        [~,Ind] = max(abs(XCorr));
        Lag = Lags(Ind);
        if Lag<0
            CorrPiezoRaw(vv) = abs(corr(Resamp_Filt_Raw_wav(1:(end+Lag))',Resamp_Filt_Logger_wav((1-Lag):end)')); % Running correlation with optimal alignment (mic ahead of logger signal)
        elseif Lag>0
            CorrPiezoRaw(vv) = abs(corr(Resamp_Filt_Raw_wav((1+Lag):end)',Resamp_Filt_Logger_wav(1:(end-Lag))')); % Running correlation with optimal alignment (logger ahead of mic signal)
        elseif Lag==0
            CorrPiezoRaw(vv) = abs(corr(Resamp_Filt_Raw_wav',Resamp_Filt_Logger_wav')); % % Running correlation with optimal alignment 
        end
        
        
        Duration(vv) = length(BioSoundCalls{vv,1}.amp);
        RMS(vv) = BioSoundCalls{vv,2}.rms;
        
        if ManualPause
            figure(1)
            sgtitle(sprintf('set %d/%d Voc %d/%d rho = %.2f',df,length(DataFiles), vv,NVoc,CorrPiezoRaw(vv)))
            APM=audioplayer(BioSoundCalls{vv,1}.sound./(max(abs(BioSoundCalls{vv,1}.sound))),BioSoundCalls{vv,1}.samprate);
            APP=audioplayer(BioSoundCalls{vv,2}.sound./(max(abs(BioSoundCalls{vv,2}.sound))),BioSoundCalls{vv,2}.samprate);
            
            
%             AP2=audioplayer(BioSoundCalls{vv,2}.sound./(max(abs(BioSoundCalls{vv,2}.sound))),BioSoundCalls{vv,2}.samprate);
%             play(AP2)
            % Rate the audio quality at the microphone
            commandwindow
            AudioGood{vv} = nan(size(BioSoundCalls{vv,2}.OnOffSets_elmts,1),1);
            ManualCallType{vv} = cell(size(BioSoundCalls{vv,2}.OnOffSets_elmts,1),1);
            ManualAnnotationOK{vv} = AudioGood{vv};
            for elmt = 1:length(AudioGood{vv})
                fprintf(1, 'Sound Element #%d\n', elmt)
                
                INPUT=[];
                while isempty(INPUT)
                    play(APM)
                    INPUT = input('Is it a vocalization? 1 yes, 0 No');
                    if isempty(INPUT) || ((INPUT~=1) && (INPUT~=0))
                        INPUT=[];
                    end
                end
                ManualAnnotationOK{vv}(elmt) = INPUT;
                if ~INPUT % This is not a vocalization
                    continue
                end
                
                INPUT=[];
                while isempty(INPUT)
                    play(APM)
                    INPUT = input('Good 4 Audio on Microphone? 1 yes, 0 No');
                    if isempty(INPUT) || ((INPUT~=1) && (INPUT~=0))
                        INPUT=[];
                    end
                end
                AudioGood{vv}(elmt) = INPUT;
                
                % Annotate the type of call
                INPUT=[];
                while isempty(INPUT)
                    play(APM)
                    pause(1)
                    play(APP)
                    INPUT = input('Trill (1), Bark(2), pitchy call(3), low buzz(4) panting (5) Low tuck (6) Squeal (7) Rattle (8) Chuckles (9) Unknown (0)');
                    if isempty(INPUT) || ((INPUT~=0) &&(INPUT~=1) && (INPUT~=2) && (INPUT~=3)&& (INPUT~=4) && (INPUT~=5)&& (INPUT~=6)&& (INPUT~=7)&& (INPUT~=8)&& (INPUT~=9))
                        INPUT=[];
                    end
                end
                if INPUT==1
                    ManualCallType{vv}{elmt} = 'Tr';
                elseif INPUT==2
                    ManualCallType{vv}{elmt} = 'Ba';
                elseif INPUT==3
                    ManualCallType{vv}{elmt} = 'Pi';
                elseif INPUT==4
                    ManualCallType{vv}{elmt} = 'Bu';
                elseif INPUT==5
                    ManualCallType{vv}{elmt} = 'Pa';
                elseif INPUT==6
                    ManualCallType{vv}{elmt} = 'LT';
                elseif INPUT==7
                    ManualCallType{vv}{elmt} = 'Sq';
                elseif INPUT==8
                    ManualCallType{vv}{elmt} = 'Ra';
                elseif INPUT==9
                    ManualCallType{vv}{elmt} = 'Ch';
                elseif INPUT==0
                    ManualCallType{vv}{elmt} = 'Un';
                end
                %             keyboard
                clf(F1)
                
            end
        end
    end
    if ManualPause
        save(fullfile(DataFile.folder, DataFile.name), 'CorrPiezoRaw','Duration', 'RMS', 'AudioGood','ManualCallType','ManualAnnotationOK', '-append')
    else
        save(fullfile(DataFile.folder, DataFile.name), 'CorrPiezoRaw','Duration', 'RMS', '-append')
    end
    clear BioSoundCalls CorrPiezoRaw Duration RMS AudioGood ManualCallType ManualAnnotationOK
end
end