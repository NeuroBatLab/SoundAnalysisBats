function [SoundEvent_LoggerSamp, SoundEvent_TranscTime_ms, Piezo_envelope_All, CallEnvInd_merge, CallEnvInd] = piezo_find_calls_logger(Data_directory, RMSfactor)
% Takes in the directory for the individidual logger data
% LOGGER_DIRECTORY and outputs SoundEvent_LoggerSamp, the start/stop indices of potential
% calls in the raw input data. Note that SoundEvent_LoggerSamp is a matrix where each
% row represents a sound event with the first column corresponding to onsets and
% the second column to offsets. This function should be used in combination with
% piezo_find_calls.m, which will call this function on all loggers for
% a recording session.
close all
%% Setting up variables that will be used later
draw_plots = false;

DurChunck = 10; % duration (in min) of each chunck treated in the function. This duration ensures that there is no significant time drift between the envelope calculation and the original sound
FS_env = 1000; % Sample frequency of the envelope
EnvWin = 5; % duration of the sliding window for the envelope calculation
BandPassFilter = [1000 5000]; % the frequencies we care about to identify when a call is made

% Call detection parameters
if nargin<2
    RMSfactor = 3; % how much greater the call's envelope needs to be than the noise. Subject to change implementation
end
CallLength = 0.007; % minimum length of a call in s
MergeThresh = 50e-3; % minimum length between two calls in s (50ms)

% Identify data
PathPieces = split(Data_directory, filesep);
Logger_name = PathPieces{find(contains(PathPieces,'extracted_data'))-1};
%disp(logger_name)

%it's saying the last call occurs at 3.8731, which doesn't make sense.
%Check via graph


%% Load the raw signal and raw signal frequency, then calculate ratio of
% raw signal frequency to what the envelope's frequency will be.
file = dir(fullfile(Data_directory, '*CSC0*'));
if isempty(file)
    error('Data file not found');
end
filepath = fullfile(file.folder, file.name);
load(filepath, 'AD_count_int16', 'Estimated_channelFS_Transceiver','Indices_of_first_and_last_samples','Timestamps_of_first_samples_usec')
SamplingFreq = round(nanmean(Estimated_channelFS_Transceiver));
if isnan(SamplingFreq)
    warning('There is no estimation of the sampling frequency for that logger, hardcode 50kHz value\n')
    SamplingFreq = 50000;
end
AD_count_double = double(AD_count_int16);
clear AD_count_int16

%% Center the signal and clear the old data from memory
Centered_piezo_signal = AD_count_double - mean(AD_count_double);
clear AD_count_double

%% Design the bandpass filter for the 1000-5000Hz range
[z,p,k] = butter(6,BandPassFilter / (SamplingFreq / 2),'bandpass');
Sos_low = zp2sos(z,p,k);


%% Cut the signal into chuncks of same short duration that is a multiple of FS_ratio, otherwise resampling accumulate errors
Signal_length = length(Centered_piezo_signal);
NumSampleChunck = round(DurChunck*60*SamplingFreq);
Chuncks = [1:NumSampleChunck:Signal_length Signal_length];
NChuncks = length(Chuncks)-1;
Piezo_centered_signal = cell(1,NChuncks); %

for ii = 1:NChuncks
    Piezo_centered_signal{ii} = Centered_piezo_signal(Chuncks(ii):(Chuncks(ii+1)));
end

%% Calculate the amplitude envelope independantly on every chunck
Piezo_envelope = cell(1,NChuncks); %
ErrorResample_ms = cell(1,NChuncks);
ErrorResampleLeft_ms = cell(1,NChuncks);
EnvCalcStart = tic;
parfor ii = 1:NChuncks %%% parallelize
    fprintf(1,'Start Envelope chunck %d/%d\n', ii, NChuncks)
    LocalStart = tic;
    Filtered_piezo_sample = filtfilt(Sos_low, 1, Piezo_centered_signal{ii});
    Filtered_sample_envelope = envelope(Filtered_piezo_sample, round(EnvWin*10^(-3) * SamplingFreq), 'rms');
    Piezo_envelope{ii} = resample(Filtered_sample_envelope, FS_env, SamplingFreq);
    % check the error by resample and randomly eliminate the corresponding
    % number of points in the envelope to suppress cumulative error of
    % allignement between envellope calculation and raw data 
    ErrorResample_ms{ii} = (length(Piezo_envelope{ii})/FS_env - length(Piezo_centered_signal{ii})/SamplingFreq)*10^3;
    if round(ErrorResample_ms{ii})>=1
        for jj=1:round(ErrorResample_ms{ii})
            Piezo_envelope{ii}(randi(length(Piezo_envelope{ii}),1))=[];
        end
    end
    ErrorResampleLeft_ms{ii} = ErrorResample_ms{ii}-round(ErrorResample_ms{ii});
    fprintf(1,'Done Envelope chunck %d/%d in %.1fs\n', ii, NChuncks, toc(LocalStart))
end
EnvCalcStop = toc(EnvCalcStart);
fprintf(1,'Done Calculating Envelope in %.1fs\n', EnvCalcStop)
ErrorResample_ms = cell2mat(ErrorResample_ms);
ErrorResampleLeft_ms = cell2mat(ErrorResampleLeft_ms);
% adjust envelope again!
CS_envError_ms = cumsum(ErrorResampleLeft_ms);
while any(CS_envError_ms>=1)
    ii = find(CS_envError_ms>=1, 1,'First');
    Piezo_envelope{ii}(randi(length(Piezo_envelope{ii}),1))=[];
    CS_envError_ms = CS_envError_ms-1;
end
% clear Centered_piezo_signal

%% Find the noise threshold for the given logger taking into account all the
% recording
fprintf(1, 'Calculating Noise threshold\n')
LocalStart = tic;
EnvLengths = cellfun('length',Piezo_envelope);
Piezo_envelope_All = reshape([Piezo_envelope{:}],1,sum(EnvLengths))';
NumSamples4Noise = 50; % this is the number of samples of the distribution of noise values
Noise = getAverageNoise(Piezo_envelope_All, FS_env, NumSamples4Noise);
fprintf(1,'Done with calculating noise: %s in %.1fs\n', num2str(Noise),toc(LocalStart))

if draw_plots
    %Display the noise threshold
    figure(2)
    subplot(3,1,1)
    plot((1:length(Piezo_envelope_All)), Piezo_envelope_All, 'Color',[0,0.5,0.9])
    hold on
    line([1, length(Piezo_envelope_All)], [Noise * RMSfactor, Noise * RMSfactor], 'Color','red','LineStyle','--', 'LineWidth',2)
    title(Logger_name)
    xlabel('Time (ms or sample)')
    legend({'Envelope', 'Sound detection threshold'})
    hold on
end

%% Detect sound events: Loop through chuncks and Create a logical vector:
% 1-> every time the data point is above the noise
% threshold and 0 -> every time it isn't
% StartTimes = cell(NChuncks,1);
% StopTimes = cell(NChuncks,1);
% EnvLengths_ = [0 EnvLengths];
% for ii = 1:NChuncks %%% parallelize
%     fprintf(1,'Find sound events in chunck %d/%d\n', ii, NChuncks)
%     CallIndicator = Piezo_envelope{ii} > (RMSfactor * Noise);
%     
%     % Start/Stop times is any 1ms bin there is a change from 1 to 0 or 0 to 1
%     StartTimes{ii} = find(diff(CallIndicator) > 0)+1 + sum(EnvLengths_(1:ii));
%     StopTimes{ii} = find(diff(CallIndicator) < 0) + sum(EnvLengths_(1:ii));
%         
% end
% 
% StartTimes = [StartTimes{:}];
% StopTimes = [StopTimes{:}];

StartTimes = (find(diff(Piezo_envelope_All>(RMSfactor * Noise))>0)+1)';
StopTimes = (find(diff(Piezo_envelope_All>(RMSfactor * Noise))<0))';
% Check edge case 1
if StopTimes(1) < StartTimes(1)
    StartTimes = [1, StartTimes];
end

% Check edge case 2
if StartTimes(end) > StopTimes(end)
    StopTimes = [StopTimes, length(Piezo_envelope_All)];
end


if length(StartTimes) ~= length(StopTimes)
    warning('Start and Stop times are not the same size');
    keyboard
end
if any(StartTimes > StopTimes)
    warning('Start time should not be greater than Stop time');
    keyboard
end



% Find the potential call times, the requirements of which are that the duration of the sound event above threshold
% is greater than callLength and also merge any calls that are less
% than mergethresh apart
LongSoundEvents = (StopTimes-StartTimes)>= (FS_env .* CallLength);
CallEnvInd = [StartTimes(LongSoundEvents)' StopTimes(LongSoundEvents)'];


if draw_plots
    %visually inspect that the previous step is correct
    figure(2)
    subplot(3,1,1)
    hold on
    legend('AutoUpdate','Off')
    plot(CallEnvInd(:,1), -250*ones(size(CallEnvInd(:,1))), 'g*')
    plot(CallEnvInd(:,2), -500*ones(size(CallEnvInd(:,2))),'m*')
    hold off
    subplot(3,1,2)
    histogram(StopTimes-StartTimes,'BinWidth',1)
    hold on
    histogram(CallEnvInd(:,2)-CallEnvInd(:,1),'BinWidth',1)
    legend('Before threshold on duration','After threshold on duration')
    hold off
    xlabel('detected sound duration (ms or samples)')
    ylabel('# sounds')
    
    subplot(3,1,3)
    histogram(StopTimes-StartTimes,'BinWidth',1)
    hold on
    histogram(CallEnvInd(:,2)-CallEnvInd(:,1),'BinWidth',1)
    legend('Before threshold on duration','After threshold on duration')
    hold off
    xlabel('detected sound duration (ms or samples)')
    ylabel('# sounds')
    xlim([-10 20])
    
end

%% Check that the detection is working properly
TotalNumSoundEvent = size(CallEnvInd,1);
if TotalNumSoundEvent == 0
    fprintf(1,'No events detected for %s\n', Logger_name);
    SoundEvent_LoggerSamp = [];
    SoundEvent_TranscTime_ms = [];
    CallEnvInd_merge = [];
else
    if draw_plots
        Step=1;
        Delay = 100; % delay to add before after each call in ms
        DBNoise = 60; % amplitude parameter for the color scale of the spectro
        FHigh = 10000; % y axis max scale for the spectrogram
        %visually inspect that the previous step is correct
        for ii=1:Step:TotalNumSoundEvent
                Fig3=figure(3);
                clf(Fig3)
                Call = CallEnvInd(ii,:);
                x_start = Call(1)-Delay;
                x_stop = Call(2)+Delay;
                Raw = Centered_piezo_signal(round(x_start/FS_env*SamplingFreq):round(x_stop/FS_env*SamplingFreq));
                Raw_ramp = cosramp(Raw-mean(Raw), SamplingFreq*10*10^-3);
                [~] = spec_only_bats(Raw_ramp,SamplingFreq,DBNoise, FHigh);
                caxis('manual');
                caxis([2 70]);
                ylim([-500 10000])
                hold on

                yyaxis right
                plot(Piezo_envelope_All(x_start:x_stop), '-k','LineWidth',2)
                hold on
                HL = hline((RMSfactor * Noise));
                HL.LineWidth=2;
                hold on
                %             line([call(1)-x_start, call(2)-x_start], max(piezo_envelope(x_start:x_stop))*ones(2,1), 'Color','g','LineStyle', '-', 'LineWidth',4)
                line([Call(1)-x_start, Call(2)-x_start], -5*ones(2,1), 'Color','g','LineStyle', '-', 'LineWidth',4)
                ylim([-10 300])
                ylabel('Amplitude Envelope')
                hold off
                title(sprintf('detection %d/%d',ii,TotalNumSoundEvent))
                Player = audioplayer(Raw_ramp/std(Raw_ramp), SamplingFreq);
                play(Player)
                pause(length(Raw_ramp)/SamplingFreq+1)
        end

    end

    %% Merge sound events that are separated by less than mergethresh
    %find which successive events we can merge
    Events2Merge = [0; (CallEnvInd(2:end,1)-CallEnvInd(1:end-1,2))<= (MergeThresh * FS_env)];
    FirstEvents2Merge = find(diff([Events2Merge; 0])==1); % onset of each sequence of events that should be merged
    LastEvents2Merge = find(diff([Events2Merge; 0])==-1);% offset of each sequence of events that should be merged
    Events2keep = strfind([Events2Merge' 0],[0 0]); % events that should be kept as they are
    if length(FirstEvents2Merge)~=length(LastEvents2Merge)
        warning('Problem in the detection of sequences of sound events to merge')
        keyboard
    end
    CallEnvInd_merge = [CallEnvInd(Events2keep,:) ; [CallEnvInd(FirstEvents2Merge,1) CallEnvInd(LastEvents2Merge,2)]];
    % reorder in increasing samp value
    [~,IndOrd]=sort(CallEnvInd_merge(:,1));

    CallEnvInd_merge = CallEnvInd_merge(IndOrd,:);
    TotalNumSoundEvent_merge = size(CallEnvInd_merge,1);

    % check that the merge work properly
    if draw_plots
        figure(5);
        subplot(2,1,1)
        histogram((CallEnvInd(2:end,1)-CallEnvInd(1:end-1,2))/FS_env*10^3, 'BinWidth',1)
        hold on
        histogram((CallEnvInd_merge(2:end,1)-CallEnvInd_merge(1:end-1,2))/FS_env*10^3, 'BinWidth',1)
        legend('Before merge', 'After merge')
        xlabel('Intercall intervals (ms)')
        subplot(2,1,2)
        histogram((CallEnvInd(2:end,1)-CallEnvInd(1:end-1,2))/FS_env*10^3, 'BinWidth',1)
        hold on
        histogram((CallEnvInd_merge(2:end,1)-CallEnvInd_merge(1:end-1,2))/FS_env*10^3, 'BinWidth',1)
        legend('Before merge', 'After merge')
        xlabel('Intercall intervals (ms)')
        xlim([0 MergeThresh*2*10^3])


        Delay = 100; % delay to add before after each call in ms
        DBNoise = 60; % amplitude parameter for the color scale of the spectro
        FHigh = 10000; % y axis max scale for the spectrogram
        Step = 1;
        %visually inspect that the previous step is correct
        for ii=1:Step:size(CallEnvInd_merge,1)
                Fig4=figure(4);
                clf(Fig4)
                Call = CallEnvInd_merge(ii,:);
                x_start = Call(1)-Delay;
                x_stop = Call(2)+Delay;
                Raw = Centered_piezo_signal(round(x_start/FS_env*SamplingFreq):round(x_stop/FS_env*SamplingFreq));
                Raw_ramp = cosramp(Raw-mean(Raw), SamplingFreq*10*10^-3);
                [~] = spec_only_bats(Raw_ramp,SamplingFreq,DBNoise, FHigh);
                caxis('manual');
                caxis([2 70]);
                ylim([-500 10000])
                hold on
                yyaxis right
                plot(Piezo_envelope_All(x_start:x_stop), '-k','LineWidth',2)
                hold on
                HL = hline((RMSfactor * Noise));
                HL.LineWidth=2;
                hold on
                %             line([call(1)-x_start, call(2)-x_start], max(piezo_envelope(x_start:x_stop))*ones(2,1), 'Color','g','LineStyle', '-', 'LineWidth',4)
                line([Call(1)-x_start, Call(2)-x_start], -5*ones(2,1), 'Color','g','LineStyle', '-', 'LineWidth',4)
                ylim([-10 300])
                ylabel('Amplitude Envelope')
                hold off
                title(sprintf('detection %d/%d after merge',ii,TotalNumSoundEvent_merge))
                Player = audioplayer(Raw_ramp/std(Raw_ramp), SamplingFreq);
                play(Player)
                pause(length(Raw_ramp)/SamplingFreq+1)
        end

    end

    %% Calculate the best estimate of the onset/offset of each sound event...
    % in the original data centered_piezo_signal or AD_countint16

    SoundEvent_LoggerSamp = round(CallEnvInd_merge/FS_env*SamplingFreq);

    %% Convert in transceiver time
    NEvents = size(SoundEvent_LoggerSamp,1);
    SoundEvent_TranscTime_ms = nan(size(SoundEvent_LoggerSamp));
    for ee=1:NEvents
        % find the File onset stamp on the logger that is closest to before
        % the snippet of sound onset
        IndTSOn = find(Indices_of_first_and_last_samples(:,1)<SoundEvent_LoggerSamp(ee,1), 1, 'Last');

        % find the File onset stamp on the logger that is closest to after
        % the snippet of sound offset
        IndTSOff = find(Indices_of_first_and_last_samples(:,2)>SoundEvent_LoggerSamp(ee,2), 1, 'First');

        if ~isempty(IndTSOff)
            if IndTSOn<=length(Estimated_channelFS_Transceiver) && ~isnan(Estimated_channelFS_Transceiver(IndTSOn))
                SoundEvent_TranscTime_ms(ee,1) = Timestamps_of_first_samples_usec(IndTSOn)*10^-3 + (SoundEvent_LoggerSamp(ee,1)-Indices_of_first_and_last_samples(IndTSOn,1)) / Estimated_channelFS_Transceiver(IndTSOn)*10^3;
            elseif ~isnan(nanmean(Estimated_channelFS_Transceiver))
                SoundEvent_TranscTime_ms(ee,1) = Timestamps_of_first_samples_usec(IndTSOn)*10^-3 + (SoundEvent_LoggerSamp(ee,1)-Indices_of_first_and_last_samples(IndTSOn,1)) / nanmean(Estimated_channelFS_Transceiver)*10^3;
            else
                SoundEvent_TranscTime_ms(ee,1) = Timestamps_of_first_samples_usec(IndTSOn)*10^-3 + (SoundEvent_LoggerSamp(ee,1)-Indices_of_first_and_last_samples(IndTSOn,1)) / SamplingFreq*10^3;
            end

            if IndTSOff<=length(Estimated_channelFS_Transceiver) && ~isnan(Estimated_channelFS_Transceiver(IndTSOff))
                SoundEvent_TranscTime_ms(ee,2) = Timestamps_of_first_samples_usec(IndTSOff)*10^-3 - (SoundEvent_LoggerSamp(ee,2)-Indices_of_first_and_last_samples(IndTSOff,1)) / Estimated_channelFS_Transceiver(IndTSOff)*10^3;
            elseif ~isnan(nanmean(Estimated_channelFS_Transceiver))
                SoundEvent_TranscTime_ms(ee,2) = Timestamps_of_first_samples_usec(IndTSOff)*10^-3 - (SoundEvent_LoggerSamp(ee,2)-Indices_of_first_and_last_samples(IndTSOff,1)) / nanmean(Estimated_channelFS_Transceiver)*10^3;
            else
                SoundEvent_TranscTime_ms(ee,2) = Timestamps_of_first_samples_usec(IndTSOff)*10^-3 - (SoundEvent_LoggerSamp(ee,2)-Indices_of_first_and_last_samples(IndTSOff,1)) / SamplingFreq*10^3;
            end

            if diff(SoundEvent_TranscTime_ms(ee,:))<=0
                if ~isnan(nanmean(Estimated_channelFS_Transceiver))
                    SoundEvent_TranscTime_ms(ee,2) = SoundEvent_TranscTime_ms(ee,1) + diff(SoundEvent_LoggerSamp(ee,:)) / nanmean(Estimated_channelFS_Transceiver)*10^3;
                else
                    SoundEvent_TranscTime_ms(ee,2) = SoundEvent_TranscTime_ms(ee,1) + diff(SoundEvent_LoggerSamp(ee,:)) / SamplingFreq*10^3;
                end
            end
        else
            % find the time stamp on the logger that is closest to before
            % the snippet of sound offset
            IndTSOff = find(Indices_of_first_and_last_samples(:,1)<SoundEvent_LoggerSamp(ee,2), 1, 'Last');
            % this vocalization is in the last recording file
            % There is no estimation of the sample frequency for that last
            % file. Let's estimate it as the average of the previous
            % estimates
            FS_local = nanmean(Estimated_channelFS_Transceiver);
            if isnan(FS_local)
                FS_local = SamplingFreq;
            end
            SoundEvent_TranscTime_ms(ee,1) = Timestamps_of_first_samples_usec(IndTSOn)*10^-3 + (SoundEvent_LoggerSamp(ee,1)-Indices_of_first_and_last_samples(IndTSOn,1)) / FS_local*10^3;
            SoundEvent_TranscTime_ms(ee,2) = Timestamps_of_first_samples_usec(IndTSOff)*10^-3 + (SoundEvent_LoggerSamp(ee,2)-Indices_of_first_and_last_samples(IndTSOff,1)) / FS_local*10^3;
        end
    end
end
%% SAVE
%
% %     save('CallTimes.mat', 'callTimes', 'piezo_envelope', 'noise', 'samplingFreq', 'AD_count_double', '-v7.3')
%
%
%    intervalTimes = zeros(1, length(callTimes) - 1);
%    for ii = 1:length(callTimes) - 1
%        nextCall = callTimes{ii + 1};
%        currCall = callTimes{ii};
%        intervalTimes(ii) = nextCall(1) - currCall(2);
%    end
%
%    %looks like 50-100ms is the sweet spot
%    histogram(intervalTimes, [0:10:1000])
%    ylabel("counts")
%    xlabel("interval (ms)")

%     load('CallTimes.mat', 'callTimes', 'piezo_envelope', 'noise', 'samplingFreq')
%     callTimes = broadenCallLength(callTimes, piezo_envelope, noise);

%     %Notes: If Formant0 exists seems to be a decent discriminator...
%     %nevermind. Call 90, 91 is a mix of a call and lots of noise so could be
%     %useful to extract information from. Can possibly view it as a mixture.
%     %In general, calls have redundant structure at three ranges, noise is
%     %usually more. Calls 1-38 exhibit this excellent noise structure, 39-59 is
%     %different. 60-65 has the fluttering sound that shows as continuous streak
%     %and Julie said was an artifact. Call 68 has beautiful structure, the
%     %noise ones in 70-75 have no solid structure and weak curves. 76-80
%     %have some more interesting mixed structures and weaker but still there
%     %vocalization curves
%
%     % Plot the biosound features of the "calls"
%     sound_directory = '/Users/ryanmoughan/Research/Data/20190623/audiologgers/logger12/extracted_data/extracted_calls/';
%     files = dir(sound_directory);
%     for ii = 90:length(files)
%        file_name = strcat('cut_call_', int2str(ii), '.wav');
%        if contains(file_name, 'wav')
%            file_path = strcat(sound_directory, file_name);
%            [Y] = audioread(file_path);
%            Y = Y - mean(Y);
%            [z,p,k] = butter(6,[1000, 5000] / (samplingFreq / 2),'bandpass');
%            sos_low = zp2sos(z,p,k);
%            filtered_Y = filtfilt(sos_low, 1, Y);
%            BiosoundObj = runBiosound(filtered_Y, samplingFreq, 10000);
%            plotBiosound(BiosoundObj, 10000, 1)
%            title(file_name)
%            pause()
%            close all
%        end
%     end

%     %plotting different features to look at
%     vocalizations = [68, 69, 76:91];
%     false_positives = (1:66);
%     feature1 = blah;
%     feature2 = bah;
%     plot(feature1(vocalizations), feature2(vocalizations), 'filled', 'd')
%     plot(feature1(false_positives), feature2(false_positives), 'red')

%     cut_calls(data_directory, callTimes)

end

% function BiosoundObj = runBiosound(Y, FS, F_high)
% % Hard coded parameters for biosound
% % spectrogram parameters
% Spec_sample_rate = 1000; % sampling rate Hz
% Freq_spacing = 50; % width of the frequency window for the FFT Hz
% Min_freq = 300; % high pass filter before FFT Hz
% Max_freq = 50000; % Low pass filter before FFT Hz
% % temporal enveloppe parameters
% Cutoff_freq = 150; % Hz
% Amp_sample_rate = 1000; % Hz
% if nargin<3
%     % Spectrum parameters
%     F_high = 50000; % frequency of Low-pass filter Hz
% end
% % Fundamental parameters
% MaxFund = 4000;
% MinFund = 300;
% LowFc = 100; %100
% HighFc = 18000;% 15000
% MinSaliency = 0.6;
% DebugFigFundest = 0;
% MinFormantFreq = 2000;
% MaxFormantBW = 1000; %500
% WindowFormant = 0.1;
% Method= 'Stack';
% 
% % create the biosound object
% BiosoundObj = py.soundsig.sound.BioSound(py.numpy.array(Y),pyargs('fs',FS));
% % methods(BiosoundFi, '-full') % this command plot all the methods with the available arguments
% 
% % Calculate the RMS (lhs std(varargin))
% BiosoundObj.rms = BiosoundObj.sound.std();
% 
% % calculate the amplitude enveloppe
% ampenv(BiosoundObj, Cutoff_freq,Amp_sample_rate);
% 
% % Calculate the periodicity of the amplitude envelope
% SoundAmp = double(py.array.array('d', py.numpy.nditer(BiosoundObj.amp)));
% [P,F] = pspectrum(SoundAmp,1000);
% [PKS,LOCS]=findpeaks(P);
% AmpPeriodF = F(LOCS(PKS == max(PKS))); % Frequency in hertz of the max peak
% AmpPeriodP = max(PKS)/mean(SoundAmp.^2); % Proportion of power in the max peak of the spectrum
% 
% % calculate the spectrum (lhs spectrum(self, f_high, pyargs))
% spectrum(BiosoundObj, F_high)
% % calculate the spectrogram (lhs spectroCalc(self, spec_sample_rate,
% % freq_spacing, min_freq, max_freq, pyargs))
% try % For very short sound, the Freq_spacing is too small, doubling if error
%     spectroCalc(BiosoundObj, Spec_sample_rate, Freq_spacing, Min_freq,Max_freq)
% catch
%     spectroCalc(BiosoundObj, Spec_sample_rate, Freq_spacing.*2, Min_freq,Max_freq)
% end
% 
% % Calculate time varying spectralmean and spectral max
% Spectro = double(BiosoundObj.spectro);
% Fo = double(BiosoundObj.fo);
% TPoints = size(Spectro,2);
% SpectralMean = nan(1,TPoints);
% %         SpectralMax = nan(1,TPoints);
% for tt=1:TPoints
%     %             SpectralMax(tt) = Fo(Spectro(:,tt)==max(Spectro(:,tt)));
%     PSDSpec = Spectro(:,tt)./(sum(Spectro(:,tt)));
%     SpectralMean(tt) = sum(PSDSpec' .* Fo);
% end
% 
% % calculate the fundamental and related values (lhs fundest(self, maxFund,
% % minFund, lowFc, highFc, minSaliency, debugFig, pyargs)
% fundest(BiosoundObj, MaxFund, MinFund,LowFc, HighFc, MinSaliency,DebugFigFundest,MinFormantFreq,MaxFormantBW,WindowFormant,Method)
% 
% % convert biosound to a strcuture
% BiosoundObj = struct(BiosoundObj);
% % Add some fields
% BiosoundObj.AmpPeriodF = AmpPeriodF;
% BiosoundObj.AmpPeriodP = AmpPeriodP;
% BiosoundObj.SpectralMean = SpectralMean;
% %         BiosoundObj.SpectralMax = SpectralMax;
% % convert all nmpy arrays to double to be able to save as matfiles
% BiosoundObj.amp = SoundAmp;
% BiosoundObj.tAmp = double(BiosoundObj.tAmp);
% BiosoundObj.spectro = double(BiosoundObj.spectro);
% BiosoundObj.to = double(BiosoundObj.to);
% BiosoundObj.fo = double(BiosoundObj.fo);
% BiosoundObj.F1 = double(BiosoundObj.F1);
% BiosoundObj.F2 = double(BiosoundObj.F2);
% BiosoundObj.F3 = double(BiosoundObj.F3);
% BiosoundObj.fpsd = double(BiosoundObj.fpsd);
% BiosoundObj.psd = double(BiosoundObj.psd);
% BiosoundObj.sal = double(BiosoundObj.sal);
% BiosoundObj.f0 = double(BiosoundObj.f0);
% BiosoundObj.f0_2 = double(BiosoundObj.f0_2);
% BiosoundObj.sound = double(BiosoundObj.sound);
% BiosoundObj.wf = double(BiosoundObj.wf);
% BiosoundObj.wt = double(BiosoundObj.wt);
% BiosoundObj.mps = double(BiosoundObj.mps);
% end
% 
% 
% function plotBiosound(BiosoundObj, F_high, FormantPlot)
% if nargin<3
%     FormantPlot=1;
% end
% % Plot the results of biosound calculations
% subplot(2,1,1)
% ColorCode = get(groot,'DefaultAxesColorOrder');
% DBNOISE =12;
% f_low = 0;
% logB = - 20*log10(abs(double(BiosoundObj.spectro)));
% maxB = max(max(logB));
% minB = maxB-DBNOISE;
% 
% imagesc(double(BiosoundObj.to)*1000,double(BiosoundObj.fo),logB);          % to is in seconds
% axis xy;
% caxis('manual');
% caxis([minB maxB]);
% cmap = spec_cmap();
% colormap(cmap);
% %         colorbar()
% 
% v_axis = axis;
% v_axis(3)=f_low;
% v_axis(4)=F_high;
% axis(v_axis);
% xlabel('time (ms)'), ylabel('Frequency');
% 
% % Plot the fundamental and formants if they were calculated
% %     if double(BiosoundFi.sal)>MinSaliency
% Legend = {'F0' 'Formant1' 'Formant2' 'Formant3'};
% IndLegend = [];
% if ~isempty(double(BiosoundObj.f0))
%     hold on
%     plot(double(BiosoundObj.to)*1000,double(BiosoundObj.f0),'r-','LineWidth',2)
%     IndLegend = [1 IndLegend];
% end
% if FormantPlot
%     hold on
%     plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F1),'Color',ColorCode(4,:),'LineWidth',2)
%     hold on
%     plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F2),'Color',ColorCode(2,:),'LineWidth',2)
%     hold on
%     if any(~isnan(double(BiosoundObj.F3)))
%         plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F3),'Color',ColorCode(3,:),'LineWidth',2)
%         IndLegend = [IndLegend 2:4];
%     else
%         IndLegend = [IndLegend 2:3];
%     end
% end
% legend(Legend(IndLegend))
% hold off
% subplot(2,1,2)
% yyaxis left
% plot((1:length(double(BiosoundObj.sound)))/BiosoundObj.samprate*1000,double(BiosoundObj.sound), 'k-','LineWidth',2)
% hold on
% YLIM = get(gca,'YLim');
% YLIM = max(abs(YLIM)).*[-1 1];
% set(gca, 'YLim', YLIM)
% SoundAmp = double(py.array.array('d', py.numpy.nditer(BiosoundObj.amp)));
% yyaxis right
% plot(double(BiosoundObj.tAmp)*1000,double(SoundAmp), 'r-', 'LineWidth',2)
% YLIM = get(gca,'YLim');
% YLIM = max(abs(YLIM)).*[-1 1];
% set(gca, 'YLim', YLIM)
% set(gca, 'XLim', v_axis(1:2))
% xlabel('Time (ms)')
% title(sprintf('AmpPeriodicity = %.3f AmpPF = %.1f Hz',BiosoundObj.AmpPeriodP, BiosoundObj.AmpPeriodF))
% hold off
% end
% 
% function [callTimes] = broadenCallLength(callTimes, piezo_envelope, noise)
% % Recovers the beginning and end of a call in CALLTIMES by looking
% % below the noise threshold NOISE of PIEZO_ENVELOPE for where the call
% % likely began and ended
% 
% for ii = 1:length(callTimes)
%     call = callTimes{ii};
%     callStart = call(1);
%     callEnd = call(2);
%     while piezo_envelope(callStart) >= noise * 1.5
%         callStart = callStart - 1;
%     end
%     while piezo_envelope(callEnd) >= noise * 1.5
%         callEnd = callEnd + 1;
%     end
%     newCall = [callStart, callEnd];
%     callTimes{ii} = newCall;
%     
% end
% 
% end

function [Average_noise] = getAverageNoise(Piezo_envelope, FS_env, Num_samples)
% Takes in PIEZO_ENVELOPE, FS_ENV, and NUM_SAMPLES and computes the
% average noise over the PIEZO_ENVELOPE. This is done by calling
% getNoiseLevel NUM_SAMPLES of time.
rng('shuffle') % reinitialize the random number generator with a seed based on the current time
Average_noise = 0;
for ns = 1:Num_samples
    Average_noise = Average_noise + getNoiseLevel(Piezo_envelope, FS_env);
end

Average_noise = Average_noise / Num_samples;
end

function [Average_noise] = getNoiseLevel(Piezo_envelope, FS_env)
% Takes in PIEZO_ENVELOPE and FS_ENV and computes the mean noise over
% an interval by first randomly selecting an interval using
% getRandomIndex, then by checking if all values in that interval are
% less than 50, and lastly returning the average of all those values.

SampleLength = FS_env * 1; % takes 1 second
Envelope_sample = 100;
while any(Envelope_sample > 50)
    % find a random index in PIEZO_ENVELOPE that will have at least SAMPLELENGTH length
    Rand_indx = round(rand(1) * (length(Piezo_envelope)-SampleLength));
    % Extract the corresponding envelope sample
    Envelope_sample = Piezo_envelope(Rand_indx: Rand_indx + SampleLength);
end
Average_noise = mean(Envelope_sample);
end


