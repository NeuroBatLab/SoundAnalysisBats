function [callTimes] = piezo_find_calls_logger(data_directory)
    % Takes in the directory for the individidual logger data
    % LOGGER_DIRECTORY and outputs CALLTIMES, the start/stop times of potential
    % calls. Note that CallTIMES is currently a cell of vectors where each
    % vector represents a call and the first element of the vector is the 
    % start of the call and the second element of the vector is the end of
    % the call. This function should be used in combination with
    % piezo_find_calls.m, which will call this function on all loggers for 
    % a recording session. 
    
    % Setting up variables that will be used later
    draw_plots = true;
    RMSfactor = 2; % how much greater the call's envelope needs to be than the noise. Subject to change implementation
    callLength = 0.007; % minimum length of a call in s
    mergeThresh = 5e-3; % minimum length between two calls in s
    FS_env = 1000; % Sample frequency of the envelope
    BandPassFilter = [1000 5000]; % the frequencies we care about to identiy when a call is made
    logger_name = data_directory(end - 23 : end - 16); % I'm assuming all names of the form loggerXX
    disp(logger_name)
    
    %it's saying the last call occurs at 3.8731, which doesn't make sense.
    %Check via graph
    
   
    % Loop through the files in the directory and find the data file. 
    % Load the raw signal and raw signal frequency, then calculate ratio of
    % raw signal frequency to what the envelope's frequency will be.
    files = dir(data_directory);
    found = 0;
    for ii = 1:length(files)
       file = files(ii);
       if contains(file.name, "CSC0")
           filepath = strcat(data_directory, file.name);
           load(filepath, 'AD_count_int16', 'Estimated_channelFS_Transceiver')
           samplingFreq = nanmean(Estimated_channelFS_Transceiver);
           FS_ratio = round(samplingFreq / FS_env); 
           found = 1;
       end
    end
    AD_count_double = double(AD_count_int16);
    clear AD_count_int16
    
    if found == 0
        ME = MException('Data file not found');
        throw(ME)
    end    

    % Center the signal and clear the old data from memory
    centered_piezo_signal = AD_count_double - mean(AD_count_double);
    
%     AD_count_double = AD_count_double(1.842963737022642e+08 - 10000: 1.842977237022640e+08 + 10000);
%     centered_piezo_signal = centered_piezo_signal(1: 1.842977237022640e+08);
    %Debug plotting
    figure(1)
    plot((1:length(AD_count_double)), AD_count_double)

    %Design the bandpass filter for the 1000-5000Hz range
    [z,p,k] = butter(6,BandPassFilter / (samplingFreq / 2),'bandpass');
    sos_low = zp2sos(z,p,k);    

    % Filter the signal and compute the envelope (using the rms). Once
    % again, remove old data from memory
    signal_length = length(centered_piezo_signal);
    piezo_envelope = zeros(1, floor(signal_length / (FS_ratio - 1)));
    envelope_length = 0;
    for ii = 1:4
        tic;
        sampleStart = 1;
        sampleEnd = floor(signal_length / 4);
        if ii == 4
            sampleEnd = length(centered_piezo_signal);
        end
        sample = centered_piezo_signal(sampleStart:sampleEnd);
        centered_piezo_signal(sampleStart: sampleEnd) = [];
        filtered_piezo_sample = filtfilt(sos_low, 1, sample);
        filtered_sample_envelope = envelope(filtered_piezo_sample, 1e-3 * 50000, 'rms');
        filtered_sample_envelope = resample(filtered_sample_envelope, 1, FS_ratio);
        piezo_envelope((envelope_length + 1) : (envelope_length + length(filtered_sample_envelope))) = filtered_sample_envelope;
        envelope_length = envelope_length + length(filtered_sample_envelope);
        toc;
    end
    piezo_envelope = piezo_envelope(1:envelope_length);
    clear centered_piezo_signal
    
    %Find the noise for the given logger
    noise = getAverageNoise(piezo_envelope, FS_env, 50);
    disp(['Done with calculating noise: ', num2str(noise)])
    
    if draw_plots
        %Display the noise threshold
        figure(2)
        plot((1:length(piezo_envelope)), piezo_envelope, 'Color',[0,0.5,0.9])
        line([1, length(piezo_envelope)], [noise * RMSfactor, noise * RMSfactor], 'Color','red','LineStyle','--')
        title(logger_name)
        hold on
    end
 
    % Create a vector of 1s every time the data point is above the noise 
    % threshold and 0s every time it isn't
    callIndicator = piezo_envelope > (RMSfactor * noise);
    
    % Start/Stop times is any time there is a change from 1 to 0 or 0 to 1
    startTimes = find(diff(callIndicator) > 0);
    stopTimes = find(diff(callIndicator) < 0);
                
    % Check edge case 1
    if stopTimes(1) < startTimes(1)
        startTimes = [1, startTimes];
    end
        
    % Check edge case 2
    if startTimes(end) > stopTimes(end)
        stopTimes = [stopTimes, length(piezo_envelope)];
    end
        
    if draw_plots
        %visually inspect that the previous step is correct
        figure(2)
        x_values = [startTimes, stopTimes];
        y_values = repmat((100),1,length(x_values));
        scatter(x_values, y_values) 
        hold on
    end

    if length(startTimes) ~= length(stopTimes)
        ME = MException('Start and Stop times are not the same size');
        throw(ME) 
    end
        
    % Find the call times, the requirements of which are that the call time 
    % is greater than callLength and also merge any calls that are less 
    % than mergethresh apart
    callTimes = cell(1, length(startTimes));
    index = 1;
    for ii = 1:length(startTimes)
        start = startTimes(ii);
        stop = stopTimes(ii);

        if start > stop
            ME = MException('Start time should not be greater than Stop time');
            throw(ME)
        end

        if stop - start >= FS_env * callLength
            callTimes{index} = [start, stop];
            index = index + 1;
%             if index ~= 1
%                 %check to see if we can merge with last call 
%                 %just this start - previous stop < mergethresh * FS_env
%                 previous_call = callTimes{index - 1};
%                 if start - previous_call(2) <= mergeThresh * FS_env
%                     callTimes{end} = [previous_call(1), stop];
%                 else
%                     callTimes{index} = [start, stop];
%                     index = index + 1;
%                 end
%             else
%                 callTimes{index} = [start, stop];
%                 index = index + 1;
%             end
        end     
    end
    callTimes = callTimes(1 : index - 1);

    if draw_plots
        %visually inspect that the previous step is correct
        for ii = 1:length(callTimes)
            figure(2)
            call = callTimes{ii};
            x_start = call(1);
            x_stop = call(2);
            line([x_start, x_stop], [1000,1000], 'Color','g','Marker', '+')
            hold on
        end
        hold off
    end
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

