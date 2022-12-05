function audioQuality_calls_Deafs(Loggers_dir, Date, ExpStartTime, ReDo)
BandPassFilter = [1000 5000];
VolFactorMic = 0.5;
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

[z,p,k] = butter(6,[100 20000]/(BioSoundCalls{ii,1}.samprate/2),'bandpass');
sos_raw_band_listen = zp2sos(z,p,k);
clear BioSoundCalls

% Loop through the datafiles
for df=1:length(DataFiles) 
    DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData%d_*.mat', Date, ExpStartTime,df)));
    if isempty(DataFile) % who calls was the earlier format
        DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date, ExpStartTime)));
    end
    fprintf(1,'Set %d/%d\nwith file %s\n', df, length(DataFiles), DataFile.name)
    load(fullfile(DataFile.folder, DataFile.name), 'ManualAnnotationOK', 'BioSoundFilenames');
    if ~ReDo && (~exist('BioSoundFilenames', 'var') || (exist('BioSoundFilenames', 'var') && exist('ManualAnnotationOK', 'var') && isempty(BioSoundFilenames)))
            fprintf(1, 'No Calls for that set\n')
            continue
    elseif ~ReDo && exist('BioSoundFilenames', 'var') && (exist('ManualAnnotationOK', 'var') && ~isempty(ManualAnnotationOK) && (~isempty(ManualAnnotationOK{end})))
        clear ManualAnnotationOK
        continue
        
    elseif ~ReDo && exist('ManualAnnotationOK', 'var') && ~isempty(ManualAnnotationOK) && isempty(ManualAnnotationOK{end})
        GoAudioGood = input(sprintf('It looks like we should start from here because the last vocalization quality is empty in AudioGood\n There is however %d empty cells\n Resume audioGood:1 skip and check the next file:0\n', sum(cellfun('isempty',ManualAnnotationOK))));
        if GoAudioGood
            Rangevv = find(cellfun('isempty',ManualAnnotationOK));
            load(fullfile(DataFile.folder, DataFile.name), 'BioSoundCalls','BioSoundFilenames', 'RMS', 'Duration','CorrPiezoRaw','ManualCallType','AudioGood')
            NVoc = size(BioSoundCalls,1);
        else
            clear ManualAnnotationOK
            continue
        end
    else
        load(fullfile(DataFile.folder, DataFile.name), 'BioSoundCalls')
        if ~exist('BioSoundCalls', 'var')
            fprintf(1, 'No Calls for that set\n')
            BioSoundCalls = [];
            AudioGood = [];
            Duration = [];
            ManualCallType = [];
            ManualAnnotationOK = [];
            fprintf(1, 'Saving data....')
            save(fullfile(DataFile.folder, DataFile.name), 'Duration', 'AudioGood','BioSoundCalls','ManualCallType','ManualAnnotationOK', '-append')
            fprintf(1, 'Done!\n')
            clear AudioGood BioSoundCalls ManualCallType ManualAnnotationOK
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
        fprintf(1,'set %d/%d Voc %d/%d\n',df,length(DataFiles), vv,NVoc)
        % Check that biosound could be run on the extract
        if ~isfield(BioSoundCalls{vv,1}, 'sound')
            fprintf(1,'Not Processed by BioSound for Microphone Data\n')
            system('killall Preview');
            if ManualPause && (~rem(jj,40) || (jj==length(Rangevv)))
                fprintf(1, 'Saving data....')
                save(fullfile(DataFile.folder, DataFile.name), 'CorrPiezoRaw','Duration', 'RMS', 'AudioGood','ManualCallType','ManualAnnotationOK', '-append')
                fprintf(1, 'Done!\n')
            end
            continue
        end
        % Open the Biosound calculation results for the piezo
        open([BioSoundFilenames{vv,2}(1:end-4) '.pdf'])
        % filter the original microphone wavfile
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
        end
        if abs(DiffLength)<=1 % else: 
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
                while isempty(IndVocStart) && Consecutive_bins>7
                    Consecutive_bins = Consecutive_bins-1;
                    IndVocStart = strfind(Vocp_logger, ones(1,Consecutive_bins));
                end
            end
            if ~isempty(IndVocStart)
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
            end
        else
            warning('There was an isue with the extraction of microphone data, not as long as Logger data, skip vocalization for correlation calculations')
        end
        Duration(vv) = length(BioSoundCalls{vv,1}.amp);
        RMS(vv) = BioSoundCalls{vv,2}.rms;
        
        if ManualPause
            ElmtMode = 0;
            figure(1)
            sgtitle(sprintf('set %d/%d Voc %d/%d rho = %.2f',df,length(DataFiles), vv,NVoc,CorrPiezoRaw(vv)))
            Raw_listen = filtfilt(sos_raw_band_listen,1,BioSoundCalls{vv,1}.sound);
            SampleMic = resample((Raw_listen - mean(Raw_listen))/(std(Raw_listen)/VolFactorMic),BioSoundCalls{vv,1}.samprate/4,BioSoundCalls{vv,1}.samprate);
            APM = audioplayer(SampleMic, BioSoundCalls{vv,1}.samprate/4,24);
%             APM=audioplayer(BioSoundCalls{vv,1}.sound./(max(abs(BioSoundCalls{vv,1}.sound))),BioSoundCalls{vv,1}.samprate);
            APP=audioplayer(BioSoundCalls{vv,2}.sound./(max(abs(BioSoundCalls{vv,2}.sound))),BioSoundCalls{vv,2}.samprate);
            
            
%             AP2=audioplayer(BioSoundCalls{vv,2}.sound./(max(abs(BioSoundCalls{vv,2}.sound))),BioSoundCalls{vv,2}.samprate);
%             play(AP2)
            % Rate the audio quality at the microphone
            commandwindow
            if isfield(BioSoundCalls{vv,2}, 'OnOffSets_elmts')
                AudioGood{vv} = nan(size(BioSoundCalls{vv,2}.OnOffSets_elmts,1),1);
                ManualCallType{vv} = cell(size(BioSoundCalls{vv,2}.OnOffSets_elmts,1),1);
            else
                AudioGood{vv} = nan;
                ManualCallType{vv} = cell(1,1);
            end
            ManualAnnotationOK{vv} = AudioGood{vv};
            for elmt = 1:length(AudioGood{vv})
                fprintf(1, 'Sound Element #%d\n', elmt)
                if ElmtMode
                    Raw_listenE = Raw_listen(BioSoundCalls{vv,1}.OnOffSets_elmts(elmt,1):BioSoundCalls{vv,1}.OnOffSets_elmts(elmt,2));
                    SampleMicE = resample((Raw_listenE - mean(Raw_listenE))/(std(Raw_listenE)/VolFactorMic),BioSoundCalls{vv,1}.samprate/4,BioSoundCalls{vv,1}.samprate);
                    APMe = audioplayer(SampleMicE, BioSoundCalls{vv,1}.samprate/4,24);
%             APM=audioplayer(BioSoundCalls{vv,1}.sound./(max(abs(BioSoundCalls{vv,1}.sound))),BioSoundCalls{vv,1}.samprate);
                    OOInd = BioSoundCalls{vv,2}.OnOffSets_elmts(elmt,1):BioSoundCalls{vv,2}.OnOffSets_elmts(elmt,2);
                    APPe=audioplayer(BioSoundCalls{vv,2}.sound(OOInd)./(max(abs(BioSoundCalls{vv,2}.sound(OOInd)))),BioSoundCalls{vv,2}.samprate);
                end
                    
                
                INPUT=[];
                while isempty(INPUT)
                    if ~ ElmtMode
                        play(APM)
                        pause(max(1, length(Raw_listen)/BioSoundCalls{vv,1}.samprate))
                        play(APP)
                    else
                        play(APMe)
                        pause(max(1, length(Raw_listenE)/BioSoundCalls{vv,1}.samprate))
                        play(APPe)
                    end
                    INPUT = input('Is it a vocalization? 1 yes, 0 No, -1 alignment error or piezo choice error or bad piezo recording, 100 plot the spectrogram of the element');
                    
                    if ~isempty(INPUT) && INPUT==100 && isfield(BioSoundCalls{vv,2}, 'spectro_elmts') % Plot the spectrogram of the elmt
                        ElmtMode = 1;
                        if elmt==1
                            Raw_listenE = Raw_listen(BioSoundCalls{vv,1}.OnOffSets_elmts(elmt,1):BioSoundCalls{vv,1}.OnOffSets_elmts(elmt,2));
                            SampleMicE = resample((Raw_listenE - mean(Raw_listenE))/(std(Raw_listenE)/VolFactorMic),BioSoundCalls{vv,1}.samprate/4,BioSoundCalls{vv,1}.samprate);
                            APMe = audioplayer(SampleMicE, BioSoundCalls{vv,1}.samprate/4,24);
                            %             APM=audioplayer(BioSoundCalls{vv,1}.sound./(max(abs(BioSoundCalls{vv,1}.sound))),BioSoundCalls{vv,1}.samprate);
                            OOInd = BioSoundCalls{vv,2}.OnOffSets_elmts(elmt,1):BioSoundCalls{vv,2}.OnOffSets_elmts(elmt,2);
                            APPe=audioplayer(BioSoundCalls{vv,2}.sound(OOInd)./(max(abs(BioSoundCalls{vv,2}.sound(OOInd)))),BioSoundCalls{vv,2}.samprate);
                        end
                        figure(2)
                        clf
                        DBNOISE =50;
                        f_low = 0;
                        F_high = 10000;
                        logB = BioSoundCalls{vv,2}.spectro_elmts{elmt};
                        maxB = max(max(logB));
                        minB = maxB-DBNOISE;
                        imagesc(double(BioSoundCalls{vv,2}.to_elmts{elmt})*1000,double(BioSoundCalls{vv,2}.fo_elmts{elmt}),logB);          % to is in seconds
                        axis xy;
                        caxis('manual');
                        caxis([minB maxB]);
                        cmap = spec_cmap();
                        colormap(cmap);
                        v_axis = axis;
                        v_axis(3)=f_low;
                        v_axis(4)=F_high;
                        axis(v_axis);
                        xlabel('time (ms)'), ylabel('Frequency');
                        title(sprintf('Set %d/%d Voc %d/%d Elmt %d/%d\n',df,length(DataFiles), vv,NVoc, elmt, length(AudioGood{vv})))
                        INPUT=[];
                    elseif isempty(INPUT) || ((INPUT~=1) && (INPUT~=0) && (INPUT~=-1))
                        INPUT=[];
                    end
                end
                ManualAnnotationOK{vv}(elmt) = INPUT;
                if ~INPUT || INPUT<0 % This is not a good vocalization
                    if elmt==length(AudioGood{vv})
                        system('killall Preview') 
                    end
                    if ManualPause && (~rem(jj,40) || (jj==length(Rangevv))) && (elmt==length(ManualAnnotationOK{vv}))
                        fprintf(1, 'Saving data....')
                        save(fullfile(DataFile.folder, DataFile.name), 'CorrPiezoRaw','Duration', 'RMS', 'AudioGood','ManualCallType','ManualAnnotationOK', '-append')
                        fprintf(1, 'Done!\n')
                    end
                    continue
                end
                
                INPUT=[];
                while isempty(INPUT)
                    if ~ ElmtMode
                        play(APM)
                    else
                        play(APMe)
                    end
                    INPUT = input('Good 4 Audio on Microphone? 1 yes, 0 No');
                    if isempty(INPUT) || ((INPUT~=1) && (INPUT~=0))
                        INPUT=[];
                    end
                end
                AudioGood{vv}(elmt) = INPUT;
                
                % Annotate the type of call
                INPUT=[];
                while isempty(INPUT)
                    if ~ ElmtMode
                        play(APM)
                        pause(max(1, length(Raw_listen)/BioSoundCalls{vv,1}.samprate))
                        play(APP)
                    else
                        play(APMe)
                        pause(max(1, length(Raw_listenE)/BioSoundCalls{vv,1}.samprate))
                        play(APPe)
                    end
                    INPUT = input('Trill (1), Bark(2), pitchy call(3), low buzz(4) panting (5) Low tuck (6) Squeal (7) Rattle (8) Chuckles (9) LoudPicthy (10) Unknown (0) BarkBuzz (24) PitchyCallBuzz (34) SquealBuzz (74)');
                    if isempty(INPUT) || ((INPUT~=0) &&(INPUT~=1) && (INPUT~=2) && (INPUT~=3)&& (INPUT~=4) && (INPUT~=5)&& (INPUT~=6)&& (INPUT~=7)&& (INPUT~=8)&& (INPUT~=9)&& (INPUT~=10) && (INPUT~=24) && (INPUT~=34) && (INPUT~=74))
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
                elseif INPUT==10
                    ManualCallType{vv}{elmt} = 'LPi';
                elseif INPUT==24
                    ManualCallType{vv}{elmt} = 'BB';
                elseif INPUT==34
                    ManualCallType{vv}{elmt} = 'PB';
                elseif INPUT==74
                    ManualCallType{vv}{elmt} = 'SB';
                elseif INPUT==0
                    ManualCallType{vv}{elmt} = 'Un';
                end
                %             keyboard
                clf(F1)
                
            end
        end
        system('killall Preview'); 
    
        if ManualPause && (~rem(jj,40) || (jj==length(Rangevv)))
            fprintf(1, 'Saving data....')
            save(fullfile(DataFile.folder, DataFile.name), 'CorrPiezoRaw','Duration', 'RMS', 'AudioGood','ManualCallType','ManualAnnotationOK', '-append')
            fprintf(1, 'Done!')
        end
    end
    clear BioSoundCalls CorrPiezoRaw Duration RMS AudioGood ManualCallType ManualAnnotationOK
end
end