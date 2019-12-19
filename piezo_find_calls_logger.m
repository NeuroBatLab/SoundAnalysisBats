function [SoundEvent_LoggerSamp] = piezo_find_calls_logger(data_directory)
    % Takes in the directory for the individidual logger data
    % LOGGER_DIRECTORY and outputs SoundEvent_LoggerSamp, the start/stop indices of potential
    % calls in the raw input data. Note that SoundEvent_LoggerSamp is a matrix where each
    % row represents a sound event with the first column corresponding to onsets and
    % the second column to offsets. This function should be used in combination with
    % piezo_find_calls.m, which will call this function on all loggers for 
    % a recording session. 
    
    %% Setting up variables that will be used later
    draw_plots = true;
    DurChunck = 10; % duration of each chunck in min. This duration ensures that there is no significant time drift between the envelope calculation and the original sound
    Fhigh_power =75; % Frequency upper bound for calculating the envelope (time running RMS)
    RMSfactor = 2; % how much greater the call's envelope needs to be than the noise. Subject to change implementation
    callLength = 0.007; % minimum length of a call in s
    mergeThresh = 5e-3; % minimum length between two calls in s
    FS_env = 1000; % Sample frequency of the envelope
    BandPassFilter = [1000 5000]; % the frequencies we care about to identify when a call is made
    PathPieces = split(data_directory, filesep);
    logger_name = PathPieces{find(contains(PathPieces,'extracted_data'))-1};
    disp(logger_name)
    
    %it's saying the last call occurs at 3.8731, which doesn't make sense.
    %Check via graph
    
   
    %% Load the raw signal and raw signal frequency, then calculate ratio of
    % raw signal frequency to what the envelope's frequency will be.
    file = dir(fullfile(data_directory, '*CSC0*'));
    if isempty(file)
        error('Data file not found');
    end  
    filepath = fullfile(file.folder, file.name);
    load(filepath, 'AD_count_int16', 'Estimated_channelFS_Transceiver')
    SamplingFreq = nanmean(Estimated_channelFS_Transceiver);
    FS_ratio = round(SamplingFreq / FS_env); 
    AD_count_double = double(AD_count_int16);
    clear AD_count_int16
    
    %% Center the signal and clear the old data from memory
    centered_piezo_signal = AD_count_double - mean(AD_count_double);
    clear AD_count_double
    
%     centered_piezo_signal = centered_piezo_signal(1.842963737022642e+08 - 10000: 1.842977237022640e+08 + 10000);
%     centered_piezo_signal = centered_piezo_signal(1: 1.842977237022640e+08);
%     if draw_plots
%         %Debug plotting
%         figure(1)
%         plot((1:length(centered_piezo_signal)), centered_piezo_signal)
%     end

    %% Design the bandpass filter for the 1000-5000Hz range
    [z,p,k] = butter(6,BandPassFilter / (SamplingFreq / 2),'bandpass');
    sos_low = zp2sos(z,p,k);    

    
    %% Cut the signal into chuncks of same short duration
    Signal_length = length(centered_piezo_signal);
    NumSampleChunck = round(DurChunck*60*SamplingFreq);
    Chuncks = [1:NumSampleChunck:Signal_length Signal_length];
    NChuncks = length(Chuncks)-1;
    piezo_centered_signal = cell(1,NChuncks); % 
    
    for ii = 1:NChuncks
        piezo_centered_signal{ii} = centered_piezo_signal(Chuncks(ii):(Chuncks(ii+1)));
    end
    
    %% Calculate the amplitude envelope independantly on every chunck
    piezo_envelope = cell(1,NChuncks); %
    piezo_envelope2 = cell(1,NChuncks); %
    EnvCalcStart = tic;
    parfor ii = 1:NChuncks %%% parallelize
        fprintf(1,'Start Envelope chunck %d/%d\n', ii, NChuncks)
        LocalStart = tic;
        filtered_piezo_sample = filtfilt(sos_low, 1, piezo_centered_signal{ii});
        filtered_sample_envelope = envelope(filtered_piezo_sample, round(1e-3 * SamplingFreq), 'rms');
        piezo_envelope2{ii} = resample(filtered_sample_envelope, 1, FS_ratio);
        piezo_envelope{ii} = running_rms(filtered_piezo_sample,SamplingFreq, Fhigh_power,FS_env);
        fprintf(1,'Done Envelope chunck %d/%d in %.1fs\n', ii, NChuncks, toc(LocalStart))
    end
    EnvCalcStop = toc(EnvCalcStart);
    fprintf(1,'Done Calculating Envelope in %.1fs\n', EnvCalcStop)
    clear centered_piezo_signal
    
    %% Find the noise for the given logger taking into account all the
    % recording
    EnvLength = sum(cellfun('length',piezo_envelope));
    piezo_envelope_All = reshape([piezo_envelope{:}],1,EnvLength)';
    noise = getAverageNoise(piezo_envelope_All, FS_env, 50);
    disp(['Done with calculating noise: ', num2str(noise)])
    
    if draw_plots
        %Display the noise threshold
        figure(2)
        plot((1:length(piezo_envelope_All)), piezo_envelope_All, 'Color',[0,0.5,0.9])
        hold on
        line([1, length(piezo_envelope_All)], [noise * RMSfactor, noise * RMSfactor], 'Color','red','LineStyle','--', 'LineWidth',2)
        title(logger_name)
        hold on
    end
 
    %% Detect sound events: Loop through chuncks and Create a logical vector:
    % 1-> every time the data point is above the noise 
    % threshold and 0 -> every time it isn't
    callTimes = cell(NChuncks,1);
    for ii = 1:NChuncks %%% parallelize
        fprintf(1,'Find sound events in chunck %d/%d\n', ii, NChuncks)
        callIndicator = piezo_envelope{ii} > (RMSfactor * noise);
        
        % Start/Stop times is any 1ms bin there is a change from 1 to 0 or 0 to 1
        startTimes = find(diff(callIndicator) > 0);
        stopTimes = find(diff(callIndicator) < 0);
        
        % Check edge case 1
        if stopTimes(1) < startTimes(1)
            startTimes = [1, startTimes];
        end
        
        % Check edge case 2
        if startTimes(end) > stopTimes(end)
            stopTimes = [stopTimes, length(piezo_envelope{ii})];
        end
        
%         if draw_plots
%             %visually inspect that the previous step is correct
%             figure(2)
%             x_values = [startTimes stopTimes];
%             y_values = repmat((100),1,size(x_values,1));
%             scatter(x_values, y_values)
%             hold on
%         end
        
        if length(startTimes) ~= length(stopTimes)
            error('Start and Stop times are not the same size');
        end
        if any(startTimes > stopTimes)
            error('Start time should not be greater than Stop time');
        end
        
        % Find the potential call times, the requirements of which are that the duration of the sound event above threshold
        % is greater than callLength and also merge any calls that are less
        % than mergethresh apart
        LongSoundEvents = (stopTimes-startTimes)>= FS_env .* callLength;
        callTimes{ii} = [startTimes(LongSoundEvents)' stopTimes(LongSoundEvents)'];
    end
    
    
    %% Check that the detection is working properly within each chunck
    TotalNumSoundEvent = sum(cellfun(@numel,callTimes))/2;
    if draw_plots
        Delay = 100; % delay to add before after each call in ms
        DBNoise = 60; % amplitude parameter for the color scale of the spectro
        FHigh = 10000; % y axis max scale for the spectrogram
        %visually inspect that the previous step is correct
        for ct = 1:length(callTimes)
            for ii=1:50:size(callTimes{ct},1)
                Fig3=figure(3);
                clf(Fig3)
                call = callTimes{ct}(ii,:);
                x_start = call(1)-Delay;
                x_stop = call(2)+Delay;
                Raw = piezo_centered_signal{ct}(round(x_start/FS_env*SamplingFreq):round(x_stop/FS_env*SamplingFreq));
                [~] = spec_only_bats(Raw,SamplingFreq,DBNoise, FHigh);
                caxis('manual');
                caxis([2 70]);
                ylim([-500 10000])
                hold on
                yyaxis right
                plot(piezo_envelope{ct}(x_start:x_stop), '-k','LineWidth',2)
                hold on
                plot(piezo_envelope2{ct}(x_start:x_stop), ':r','LineWidth',2)
                hold on
                %             line([call(1)-x_start, call(2)-x_start], max(piezo_envelope(x_start:x_stop))*ones(2,1), 'Color','g','LineStyle', '-', 'LineWidth',4)
                line([call(1)-x_start, call(2)-x_start], -5*ones(2,1), 'Color','g','LineStyle', '-', 'LineWidth',4)
                ylim([-10 300])
                hold off
                title(sprintf('detection %d/%d',sum(cellfun(@numel,callTimes(1:(ct-1))))/2+ii,TotalNumSoundEvent))
                pause(1)
            end
        end
        
    end
    
    %% Calculate the best estimate of the onset/offset of each sound event...
    % in the original data centered_piezo_signal or AD_countint16
    Length_raw_chuncks = [0 cellfun('length',piezo_centered_signal)];
    SoundEvent_LoggerSamp_local = cell(length(callTimes),1);
    for ct = 1:length(callTimes)
        SoundEvent_LoggerSamp_local{ct} = (round(callTimes{ct}/FS_env*SamplingFreq) + sum(Length_raw_chuncks(1:ct)))';
    end
    SoundEvent_LoggerSamp_local = reshape([SoundEvent_LoggerSamp_local{:}],2,TotalNumSoundEvent)';
    
     %% Merge sound events that are separated by less than mergethresh
     %find which successive events we can merge
     Events2Merge = [0; (SoundEvent_LoggerSamp_local(2:end,1)-SoundEvent_LoggerSamp_local(1:end-1,2))<= (mergeThresh * SamplingFreq)];
     FirstEvents2Merge = find([diff(Events2Merge); 0]==1); % onset of each sequence of events that should be merged
     LastEvents2Merge = find([diff(Events2Merge); 0]==-1);% offset of each sequence of events that should be merged
     Events2keep = strfind([Events2Merge' 0],[0 0]); % events that should be kept as they are
     if length(FirstEvents2Merge)~=length(LastEvents2Merge)
         error('Problem in the detection of sequences of sound events to merge')
     end
     SoundEvent_LoggerSamp = [SoundEvent_LoggerSamp_local(Events2keep,:) ; [SoundEvent_LoggerSamp_local(FirstEvents2Merge,1) SoundEvent_LoggerSamp_local(LastEvents2Merge,2)]];
     % reorder in increasing samp value
     [~,IndOrd]=sort(SoundEvent_LoggerSamp(:,1));
     SoundEvent_LoggerSamp = SoundEvent_LoggerSamp(IndOrd,:);
     
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
        spectroCalc(BiosoundObj, Spec_sample_rate, Freq_spacing.*2, Min_freq,Max_freq)
    end

    % Calculate time varying spectralmean and spectral max
    Spectro = double(BiosoundObj.spectro);
    Fo = double(BiosoundObj.fo);
    TPoints = size(Spectro,2);
    SpectralMean = nan(1,TPoints);
    %         SpectralMax = nan(1,TPoints);
    for tt=1:TPoints
        %             SpectralMax(tt) = Fo(Spectro(:,tt)==max(Spectro(:,tt)));
        PSDSpec = Spectro(:,tt)./(sum(Spectro(:,tt)));
        SpectralMean(tt) = sum(PSDSpec' .* Fo);
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
    BiosoundObj.sound = double(BiosoundObj.sound);
    BiosoundObj.wf = double(BiosoundObj.wf);
    BiosoundObj.wt = double(BiosoundObj.wt);
    BiosoundObj.mps = double(BiosoundObj.mps);
end


function plotBiosound(BiosoundObj, F_high, FormantPlot)
    if nargin<3
        FormantPlot=1;
    end
    % Plot the results of biosound calculations
    subplot(2,1,1)
    ColorCode = get(groot,'DefaultAxesColorOrder');
    DBNOISE =12;
    f_low = 0;
    logB = - 20*log10(abs(double(BiosoundObj.spectro)));
    maxB = max(max(logB));
    minB = maxB-DBNOISE;

    imagesc(double(BiosoundObj.to)*1000,double(BiosoundObj.fo),logB);          % to is in seconds
    axis xy;
    caxis('manual');
    caxis([minB maxB]);
    cmap = spec_cmap();
    colormap(cmap);
    %         colorbar()

    v_axis = axis;
    v_axis(3)=f_low;
    v_axis(4)=F_high;
    axis(v_axis);
    xlabel('time (ms)'), ylabel('Frequency');

    % Plot the fundamental and formants if they were calculated
    %     if double(BiosoundFi.sal)>MinSaliency
    Legend = {'F0' 'Formant1' 'Formant2' 'Formant3'};
    IndLegend = [];
    if ~isempty(double(BiosoundObj.f0))
        hold on
        plot(double(BiosoundObj.to)*1000,double(BiosoundObj.f0),'r-','LineWidth',2)
        IndLegend = [1 IndLegend];
    end
    if FormantPlot
        hold on
        plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F1),'Color',ColorCode(4,:),'LineWidth',2)
        hold on
        plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F2),'Color',ColorCode(2,:),'LineWidth',2)
        hold on
        if any(~isnan(double(BiosoundObj.F3)))
            plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F3),'Color',ColorCode(3,:),'LineWidth',2)
            IndLegend = [IndLegend 2:4];
        else
            IndLegend = [IndLegend 2:3];
        end
    end
    legend(Legend(IndLegend))
    hold off
    subplot(2,1,2)
    yyaxis left
    plot((1:length(double(BiosoundObj.sound)))/BiosoundObj.samprate*1000,double(BiosoundObj.sound), 'k-','LineWidth',2)
    hold on
    YLIM = get(gca,'YLim');
    YLIM = max(abs(YLIM)).*[-1 1];
    set(gca, 'YLim', YLIM)
    SoundAmp = double(py.array.array('d', py.numpy.nditer(BiosoundObj.amp)));
    yyaxis right
    plot(double(BiosoundObj.tAmp)*1000,double(SoundAmp), 'r-', 'LineWidth',2)
    YLIM = get(gca,'YLim');
    YLIM = max(abs(YLIM)).*[-1 1];
    set(gca, 'YLim', YLIM)
    set(gca, 'XLim', v_axis(1:2))
    xlabel('Time (ms)')
    title(sprintf('AmpPeriodicity = %.3f AmpPF = %.1f Hz',BiosoundObj.AmpPeriodP, BiosoundObj.AmpPeriodF))
    hold off
end

function [callTimes] = broadenCallLength(callTimes, piezo_envelope, noise)
    % Recovers the beginning and end of a call in CALLTIMES by looking 
    % below the noise threshold NOISE of PIEZO_ENVELOPE for where the call 
    % likely began and ended

    for ii = 1:length(callTimes)
        call = callTimes{ii};
        callStart = call(1);
        callEnd = call(2);
        while piezo_envelope(callStart) >= noise * 1.5
            callStart = callStart - 1;
        end
        while piezo_envelope(callEnd) >= noise * 1.5
            callEnd = callEnd + 1;
        end
        newCall = [callStart, callEnd];
        callTimes{ii} = newCall;
        
    end

end

function [average_noise] = getAverageNoise(piezo_envelope, FS_env, num_samples)
    % Takes in PIEZO_ENVELOPE, FS_ENV, and NUM_SAMPLES and computes the
    % average noise over the PIEZO_ENVELOPE. This is done by calling
    % getNoiseLevel NUM_SAMPLES of time. 
    
    average_noise = 0;
    for i = 1:num_samples
        average_noise = average_noise + getNoiseLevel(piezo_envelope, FS_env);
    end
    
    average_noise = average_noise / num_samples;

end

function [average_noise] = getNoiseLevel(piezo_envelope, FS_env)
    % Takes in PIEZO_ENVELOPE and FS_ENV and computes the mean noise over
    % an interval by first randomly selecting an interval using
    % getRandomIndex, then by checking if all values in that interval are
    % less than 50, and lastly returning the average of all those values.
    
    sampleLength = FS_env * 1;
    rand_indx = uint64(getRandomIndex(piezo_envelope, sampleLength));    
    repeat = true;
    while repeat
        total = 0;
        repeat = false;
        envelope_sample = piezo_envelope(rand_indx: rand_indx + sampleLength);
        if any(envelope_sample > 50)
            rand_indx = uint64(getRandomIndex(piezo_envelope, sampleLength));
            repeat = true;
        else
            total = sum(envelope_sample);
        end
    end
    
    average_noise = total / sampleLength;
end

function [rand_indx] = getRandomIndex(piezo_envelope, sampleLength)
    % Takes in PIEZO_ENVELOPE and SAMPLELENGTH and returns a random index
    % in PIEZO_ENVELOPE that will have at least SAMPLELENGTH length
    
    rand_indx = rand * length(piezo_envelope);
    while (rand_indx + sampleLength) > length(piezo_envelope)
        rand_indx = rand * length(piezo_envelope);
    end
end

