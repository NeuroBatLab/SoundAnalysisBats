function [IndVocStartRaw, IndVocStopRaw, IndVocStartPiezo, IndVocStopPiezo] = who_calls(Loggers_dir, Date, ExpStartTime, Manual)
load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'Piezo_wave', 'Piezo_FS',  'Raw_wave','FS', 'RatioRMS', 'DiffRMS','BandPassFilter', 'AudioLogs', 'RMSHigh', 'RMSLow','VocFilename');
%% Identify sound elements in each vocalization extract
% Run a time window of duration 2 ms to identify who is vocalizing
% create 2 signals for the logger, one high pass filtered above 5kHz
% and the other one low pass filtered at 5kHz. The vocalizer will have
% the highest energy of all vocalizers and will have more energy in the
% lower compare to higher filtered signal

% parameters
Consecutive_bins = 20; % Number of consecutive bins of the envelope difference between highpass and low pass logger signal that has to be higher than threshold to be considered as a vocalization
DB_noise = 60; % Noise threshold for the spectrogram colormap
FHigh_spec = 90000; % Max frequency (Hz) for the raw data spectrogram
BandPassFilter = [1000 5000 9900]; % Frequency bands chosen for digital signal processing
Fhigh_power =20; % Frequency upper bound for calculating the envelope (time running RMS)
Fs_env = 1000; % Sample freqency of the enveloppe
% Buffer=100; % Maximum lag calculated for the cross correlation in ms, not used in that code for identification puposes but still calculated as of now
MergeThresh = 500; % any 2 detected call spaced by less than MergeThresh ms are merged

% Initialize variables
Nvoc = length(Raw_wave); %#ok<USENS>
Amp_env_LowPassLogVoc = cell(1,Nvoc);
Amp_env_HighPassLogVoc = cell(1,Nvoc);
Amp_env_Mic = cell(1,Nvoc);
LowPassLogVoc = cell(1,Nvoc);
% LowPassMicVoc = cell(1,Nvoc);
Fns_AL = fieldnames(Piezo_wave);
IndVocStart_all = cell(1,Nvoc);
IndVocStop_all = cell(1,Nvoc);
IndVocStartRaw_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the raw recording
IndVocStopRaw_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording
IndVocStartPiezo_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording
IndVocStopPiezo_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording
RMSRatio_all = cell(1,Nvoc);
RMSDiff_all = cell(1,Nvoc);


% design filters of raw ambient recording, bandpass and low pass which was
% used for the cross correlation
[z,p,k] = butter(6,[BandPassFilter(1) 90000]/(FS/2),'bandpass');
sos_raw_band = zp2sos(z,p,k);
% [z,p,k] = butter(6,BandPassFilter(1:2)/(FS/2),'bandpass');
% sos_raw_low = zp2sos(z,p,k);

%% Loop through vocalizations sequences and calculate amplitude envelopes
for vv=1:Nvoc
    %% First calculate the time varying RMS of the ambient microphone
    Amp_env_LowPassLogVoc{vv} = cell(length(AudioLogs),1);
    Amp_env_HighPassLogVoc{vv} = cell(length(AudioLogs),1);
    LowPassLogVoc{vv} = cell(length(AudioLogs),1);
    % bandpass filter the ambient mic recording
    Filt_RawVoc = filtfilt(sos_raw_band,1,Raw_wave{vv}); %#ok<IDISVAR>
    Amp_env_Mic{vv} = running_rms(Filt_RawVoc, FS, Fhigh_power, Fs_env);
    % Low passfilter the ambient mic recording (used for cross-correlation
    % but not done anymore
%     LowPassMicVoc{vv} = filtfilt(sos_raw_low,1,Raw_wave{vv});
    % Plot the spectrogram of the ambient microphone
    F1=figure(1);
    ColorCode = get(groot,'DefaultAxesColorOrder');
    subplot(length(AudioLogs)+2,1,1)
    [~] = spec_only_bats(Filt_RawVoc, FS, DB_noise, FHigh_spec);
    hold on
    yyaxis right
    plot((1:length(Amp_env_Mic{vv}))/Fs_env*1000, Amp_env_Mic{vv}, 'r-', 'LineWidth',2)
    ylabel('Amplitude')
    title(sprintf('Ambient Microphone Voc %d/%d',vv,Nvoc))
    pause(0.1)
    Player= audioplayer((Raw_wave{vv} - mean(Raw_wave{vv}))/std(Raw_wave{vv}), FS); %#ok<TNMLP>
    play(Player)
    pause(1)
    %% Loop through the loggers and calculate envelopes
    for ll=1:length(AudioLogs)
        % design the filters
        [z,p,k] = butter(6,BandPassFilter(1:2)/(Piezo_FS.(Fns_AL{ll})(vv)/2),'bandpass');
        sos_low = zp2sos(z,p,k);
        [z,p,k] = butter(6,BandPassFilter(2:3)/(Piezo_FS.(Fns_AL{ll})(vv)/2),'bandpass');
        sos_high = zp2sos(z,p,k);
        % filter the loggers' signals
        if sum(isnan(Piezo_wave.(Fns_AL{ll}){vv}))~=length(Piezo_wave.(Fns_AL{ll}){vv})
            LowPassLogVoc{vv}{ll} = (filtfilt(sos_low,1,Piezo_wave.(Fns_AL{ll}){vv})); % low-pass filter the voltage trace
            HighPassLogVoc = (filtfilt(sos_high,1,Piezo_wave.(Fns_AL{ll}){vv})); % high-pass filter the voltage trace
            Amp_env_LowPassLogVoc{vv}{ll}=running_rms(LowPassLogVoc{vv}{ll}, Piezo_FS.(Fns_AL{ll})(vv), Fhigh_power, Fs_env);
            Amp_env_HighPassLogVoc{vv}{ll}=running_rms(HighPassLogVoc, Piezo_FS.(Fns_AL{ll})(vv), Fhigh_power, Fs_env);

            % Plot the low pass filtered signal of each logger
            figure(1)
            subplot(length(AudioLogs)+2,1,ll+1)
            [~] = spec_only_bats(LowPassLogVoc{vv}{ll}, Piezo_FS.(Fns_AL{ll})(vv));
            title(sprintf('%s',Fns_AL{ll}))
            hold on
            yyaxis right
            plot((1:length(Amp_env_LowPassLogVoc{vv}{ll}))/Fs_env*1000, Amp_env_LowPassLogVoc{vv}{ll}, 'b-','LineWidth', 2);
            hold on
            plot((1:length(Amp_env_HighPassLogVoc{vv}{ll}))/Fs_env*1000, Amp_env_HighPassLogVoc{vv}{ll}, 'r-','LineWidth',2);
            ylabel('Amplitude')
            hold off
            pause(0.1)
            Player= audioplayer((Piezo_wave.(Fns_AL{ll}){vv}-mean(Piezo_wave.(Fns_AL{ll}){vv}))/std(Piezo_wave.(Fns_AL{ll}){vv}), Piezo_FS.(Fns_AL{ll})(vv)); %#ok<TNMLP>
            play(Player)
            pause(1)
        else
            Amp_env_LowPassLogVoc{vv}{ll}=resample(nan(1,length(Piezo_wave.(Fns_AL{ll}){vv})), Fs_env, round(Piezo_FS.(Fns_AL{ll})(vv)));
            Amp_env_HighPassLogVoc{vv}{ll}=resample(nan(1,length(Piezo_wave.(Fns_AL{ll}){vv})), Fs_env, round(Piezo_FS.(Fns_AL{ll})(vv)));
        end 
    end
    % Treat the case where some approximations of the calculations led to
    % slight diffreent sizes of running RMS
    Short = min(cellfun('length',Amp_env_LowPassLogVoc{vv}));
    Long = max(cellfun('length',Amp_env_LowPassLogVoc{vv}));
    if (Long-Short)>1
        error('The length of vectors of running RMS are too different than expected, please check!\n')
    elseif Long~=Short
        Amp_env_LowPassLogVoc{vv} = cellfun(@(X) X(1:Short), Amp_env_LowPassLogVoc{vv}, 'UniformOutput',false);
        Amp_env_HighPassLogVoc{vv} = cellfun(@(X) X(1:Short), Amp_env_HighPassLogVoc{vv}, 'UniformOutput',false);
    end
    Amp_env_LowPassLogVoc_MAT = cell2mat(Amp_env_LowPassLogVoc{vv});
    Amp_env_HighPassLogVoc_MAT = cell2mat(Amp_env_HighPassLogVoc{vv});
    RatioAmp = Amp_env_LowPassLogVoc_MAT./Amp_env_HighPassLogVoc_MAT;
    DiffAmp = Amp_env_LowPassLogVoc_MAT-Amp_env_HighPassLogVoc_MAT;
    
    % Plot the ratio of time varying RMS, the difference in time varying
    % RMS between the high and low frequency bands and the absolute time
    % varying RMS of the low frequency band
    F2=figure(2);
    subplot(3,1,1)
    plot(RatioAmp', 'LineWidth',2)
    title('Frequency bands Ratio')
    legend(Fns_AL)
    subplot(3,1,2)
    plot(DiffAmp', 'LineWidth',2)
    title('Frequency bands Diff')
    legend(Fns_AL)
    subplot(3,1,3)
    plot(Amp_env_LowPassLogVoc_MAT', 'LineWidth',2)
    title(sprintf('Running RMS %dHz-%dHz', BandPassFilter(1:2)))
    legend(Fns_AL)
    
    %% Find out which calls are emitted by each individual
    Vocp = nan(size(DiffAmp));
    IndVocStart = cell(length(AudioLogs),1);
    IndVocStop = cell(length(AudioLogs),1);
    IndVocStartRaw = cell(length(AudioLogs),1);
    IndVocStopRaw = cell(length(AudioLogs),1);
    IndVocStartPiezo = cell(length(AudioLogs),1);
    IndVocStopPiezo = cell(length(AudioLogs),1);
%     Xcor = cell(length(AudioLogs),1);
    RMSRatio = cell(length(AudioLogs),1);
    RMSDiff = cell(length(AudioLogs),1);
    for ll=1:length(AudioLogs)
%         Vocp(ll,:)=(DiffAmp(ll,:)>(4*DiffRMS.(Fns_AL{ll})(1))) .* (RatioAmp(ll,:)>(4*RatioRMS.(Fns_AL{ll})(1)));
        Vocp(ll,:) = Amp_env_LowPassLogVoc_MAT(ll,:)>(2*RMSLow.(Fns_AL{ll})(1));
        IndVocStart{ll} = strfind(Vocp(ll,:), ones(1,Consecutive_bins));
        if ~isempty(IndVocStart{ll}) % Some vocalizations were detected
            IndVocStart_diffind = find(diff(IndVocStart{ll})>1);
            IndVocStart{ll} = [IndVocStart{ll}(1) IndVocStart{ll}(IndVocStart_diffind +1)];
            NV = length(IndVocStart{ll});
            IndVocStop{ll} = nan(1,NV);
%             Xcor{ll} = nan(NV,1);
            RMSRatio{ll} = nan(NV,1);
            RMSDiff{ll} = nan(NV,1);
            Call1Hear0_temp = nan(NV,1);
            for ii=1:NV
                IVStop = find(Vocp(ll,IndVocStart{ll}(ii):end)==0, 1, 'first');
                if ~isempty(IVStop)
                    IndVocStop{ll}(ii) = IndVocStart{ll}(ii) + IVStop -1;
                else
                    IndVocStop{ll}(ii) = length(Vocp);
                end
                IndVocStartRaw{ll}(ii) = round(IndVocStart{ll}(ii)/Fs_env*FS);
                IndVocStopRaw{ll}(ii) = round(IndVocStop{ll}(ii)/Fs_env*FS);
                IndVocStartPiezo{ll}(ii) = round(IndVocStart{ll}(ii)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                IndVocStopPiezo{ll}(ii) = round(IndVocStop{ll}(ii)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                
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
                Call1Hear0_temp(ii) = RMSRatio{ll}(ii) > 3*RatioRMS.(Fns_AL{ll})(1);
                
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
        
            % Only keep vocalizations that were produced by the bat
            % according to the high RMSRatio value
            IndVocStart{ll} = IndVocStart{ll}(logical(Call1Hear0_temp));
            IndVocStop{ll} = IndVocStop{ll}(logical(Call1Hear0_temp));
            NV = sum(Call1Hear0_temp);
            
            if NV
                % Now merge detected vocalizations that are within MergeThresh
                Merge01 = (IndVocStart{ll}(2:end)-IndVocStop{ll}(1:(end-1)) < (MergeThresh/1000*Fs_env));
                NV = NV-sum(Merge01);
                IndVocStartRaw{ll} = nan(1,NV);
                IndVocStopRaw{ll} = nan(1,NV);
                IndVocStartPiezo{ll} = nan(1,NV);
                IndVocStopPiezo{ll} = nan(1,NV);
                CutInd = find(~Merge01);
                IndVocStartRaw{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*FS);
                IndVocStartPiezo{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                IndVocStopRaw{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*FS);
                IndVocStopPiezo{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                for cc=1:length(CutInd)
                    IndVocStopRaw{ll}(cc) = round(IndVocStop{ll}( CutInd(cc))/Fs_env*FS);
                    IndVocStopPiezo{ll}(cc) = round(IndVocStop{ll}(CutInd(cc))/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                    IndVocStartRaw{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*FS);
                    IndVocStartPiezo{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                end

                % Now plot the onset/offset of each extract on the spectrograms
                figure(1)
                for ii=1:NV 
                    subplot(length(AudioLogs)+2,1,ll+1)
                    hold on
                    yyaxis right
                    YLim = get(gca, 'YLim');
                    plot([IndVocStartRaw{ll}(ii) IndVocStopRaw{ll}(ii)]/FS*1000, [YLim(2)*4/5 YLim(2)*4/5], 'k-', 'LineWidth',2)
                    subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
                    hold on
                    plot([IndVocStartRaw{ll}(ii) IndVocStopRaw{ll}(ii)]/FS*1000, [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
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
        end
    end
    figure(1)
    subplot(length(AudioLogs)+2,1,1)
    XLIM = get(gca, 'XLim');
    subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
    set(gca, 'YTick', 1:length(AudioLogs))
    set(gca, 'YTickLabel', Fns_AL)
    set(gca, 'YLim', [0 length(AudioLogs)+1])
    set(gca, 'XLim', XLIM);
    ylabel('AL ID')
    xlabel('Time (ms)')
    [~,FileVoc]=fileparts(VocFilename{vv}); %#ok<IDISVAR,USENS>
    saveas(F1,fullfile(Loggers_dir,sprintf('%s_whocalls_spec', FileVoc)),'epsc')
    saveas(F2,fullfile(Loggers_dir,sprintf('%s_whocalls_RMS', FileVoc)),'epsc')
    if Manual
        pause()
    else
        pause(1)
    end
    clf(F1)
    clf(F2)
    
    % Gather Vocalization production data:
    IndVocStart_all{vv} = IndVocStart;
    IndVocStop_all{vv} = IndVocStop;
    IndVocStartRaw_merged{vv} = IndVocStartRaw;
    IndVocStopRaw_merged{vv} = IndVocStopRaw;
    IndVocStartPiezo_merged{vv} = IndVocStartPiezo;
    IndVocStopPiezo_merged{vv} = IndVocStopPiezo;
    RMSRatio_all{vv} = RMSRatio;
    RMSDiff_all{vv} = RMSDiff;
    
end


save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all', '-append');
end