function [callTimes] = piezo_find_calls_logger(data_directory)
%     % Takes in the directory for the individidual logger data
%     % LOGGER_DIRECTORY and outputs CALLTIMES, the start/stop times of potential
%     % calls. Note that CallTIMES is currently a cell of vectors where each
%     % vector represents a call and the first element of the vector is the 
%     % start of the call and the second element of the vector is the end of
%     % the call. This function should be used in combination with
%     % piezo_find_calls.m, which will call this function on all loggers for 
%     % a recording session. 
%     
%     % Setting up variables that will be used later
%     draw_plots = true;
%     RMSfactor = 2; % how much greater the call's envelope needs to be than the noise. Subject to change implementation
%     callLength = 0.007; % minimum length of a call in s
%     mergeThresh = 5e-3; % minimum length between two calls in s
%     FS_env = 1000; % Sample frequency of the envelope
%     BandPassFilter = [1000 5000]; % the frequencies we care about to identiy when a call is made
%     logger_name = data_directory(end - 23 : end - 16); % I'm assuming all names of the form loggerXX
%     disp(logger_name)
%     
%    
%     % Loop through the files in the directory and find the data file. 
%     % Load the raw signal and raw signal frequency, then calculate ratio of
%     % raw signal frequency to what the envelope's frequency will be.
%     files = dir(data_directory);
%     found = 0;
%     for ii = 1:length(files)
%        file = files(ii);
%        if contains(file.name, "CSC0")
%            filepath = strcat(data_directory, file.name);
%            load(filepath, 'AD_count_int16', 'Estimated_channelFS_Transceiver')
%            samplingFreq = nanmean(Estimated_channelFS_Transceiver);
%            FS_ratio = round(samplingFreq / FS_env); 
%            found = 1;
%        end
%     end
%     AD_count_double = double(AD_count_int16);
%     clear AD_count_int16
%     
%     if found == 0
%         ME = MException('Data file not found');
%         throw(ME)
%     end    
%     
% %     %Ryan debug info
% %     %6,254,719 : 6,255,619
% %     %really does not sound like anything but noise
% %     data = AD_count_int16(6229719 : 6378800);
% %     audiowrite("DebugCall.wav", data, round(samplingFreq))
% %     figure(100)
% %     [~]=spec_only_bats(data, 50000, 60, 10000);
% %     hold off
% 
% 
%     % Center the signal and clear the old data from memory
%     centered_piezo_signal = AD_count_double - mean(AD_count_double);
% 
%     %Design the bandpass filter for the 1000-5000Hz range
%     [z,p,k] = butter(6,BandPassFilter / (samplingFreq / 2),'bandpass');
%     sos_low = zp2sos(z,p,k);    
% 
%     % Filter the signal and compute the envelope (using the rms). Once
%     % again, remove old data from memory
%     signal_length = length(centered_piezo_signal);
%     piezo_envelope = zeros(1, floor(signal_length / (FS_ratio - 1)));
%     envelope_length = 0;
%     for ii = 1:4
%         tic;
%         sampleStart = 1;
%         sampleEnd = floor(signal_length / 4);
%         if ii == 4
%             sampleEnd = length(centered_piezo_signal);
%         end
%         sample = centered_piezo_signal(sampleStart:sampleEnd);
%         centered_piezo_signal(sampleStart: sampleEnd) = [];
%         filtered_piezo_sample = filtfilt(sos_low, 1, sample);
%         filtered_sample_envelope = envelope(filtered_piezo_sample, 1e-3 * 50000, 'rms');
%         filtered_sample_envelope = resample(filtered_sample_envelope, 1, FS_ratio);
%         piezo_envelope(envelope_length + 1 : envelope_length + length(filtered_sample_envelope)) = filtered_sample_envelope;
%         envelope_length = envelope_length + length(filtered_sample_envelope);
%         toc;
%     end
%     piezo_envelope = piezo_envelope(1:envelope_length);
%     clear centered_piezo_signal
%     
%     %Find the noise for the given logger
%     noise = getAverageNoise(piezo_envelope, FS_env, 50);
%     disp(['Done with calculating noise: ', num2str(noise)])
%     
%     if draw_plots
%         %Display the noise threshold
%         plot((1:length(piezo_envelope)), piezo_envelope, 'Color',[0,0.5,0.9])
%         line([1, length(piezo_envelope)], [noise * RMSfactor, noise * RMSfactor], 'Color','red','LineStyle','--')
%         title(logger_name)
%         hold on
%     end
%  
%     % Create a vector of 1s every time the data point is above the noise 
%     % threshold and 0s every time it isn't
%     callIndicator = piezo_envelope > (RMSfactor * noise);
%     
%     % Start/Stop times is any time there is a change from 1 to 0 or 0 to 1
%     startTimes = find(diff(callIndicator) > 0);
%     stopTimes = find(diff(callIndicator) < 0);
%                 
%     % Check edge case 1
%     if stopTimes(1) < startTimes(1)
%         startTimes = [1, startTimes];
%     end
%         
%     % Check edge case 2
%     if startTimes(end) > stopTimes(end)
%         stopTimes = [stopTimes, length(piezo_envelope)];
%     end
%         
%     if draw_plots
%         %visually inspect that the previous step is correct
%         x_values = [startTimes, stopTimes];
%         y_values = repmat((100),1,length(x_values));
%         scatter(x_values, y_values) 
%         hold on
%     end
% 
%     if length(startTimes) ~= length(stopTimes)
%         ME = MException('Start and Stop times are not the same size');
%         throw(ME) 
%     end
%         
%     % Find the call times, the requirements of which are that the call time 
%     % is greater than callLength and also merge any calls that are less 
%     % than mergethresh apart
%     callTimes = cell(1, length(startTimes));
%     index = 1;
%     for ii = 1:length(startTimes)
%         start = startTimes(ii);
%         stop = stopTimes(ii);
% 
%         if start > stop
%             ME = MException('Start time should not be greater than Stop time');
%             throw(ME)
%         end
% 
%         if stop - start >= FS_env * callLength
%             callTimes{index} = [start, stop];
%             index = index + 1;
% %             if index ~= 1
% %                 %check to see if we can merge with last call 
% %                 %just this start - previous stop < mergethresh * FS_env
% %                 previous_call = callTimes{index - 1};
% %                 if start - previous_call(2) <= mergeThresh * FS_env
% %                     callTimes{end} = [previous_call(1), stop];
% %                 else
% %                     callTimes{index} = [start, stop];
% %                     index = index + 1;
% %                 end
% %             else
% %                 callTimes{index} = [start, stop];
% %                 index = index + 1;
% %             end
%         end     
%     end
%     callTimes = callTimes(1 : index - 1);
% 
%     if draw_plots
%         %visually inspect that the previous step is correct
%         for ii = 1:length(callTimes)
%             call = callTimes{ii};
%             x_start = call(1);
%             x_stop = call(2);
%             line([x_start, x_stop], [1000,1000], 'Color','g','Marker', '+')
%             hold on
%         end
%         hold off
%     end
%     
%     save('CallTimes.mat', 'callTimes', 'piezo_envelope', 'noise', 'samplingFreq', 'AD_count_double', '-v7.3')
    
% 
%    intervalTimes = zeros(1, length(callTimes) - 1);
%    for ii = 1:length(callTimes) - 1
%        nextCall = callTimes{ii + 1};
%        currCall = callTimes{ii};
%        intervalTimes(ii) = nextCall(1) - currCall(2);
%    end
%  
% 
%     %ones I don't know about: 3, 4, 5, 9, 26, 29, 31-36,
%     labels = [6,7,11:21,23];
%    
%    noise_indx = [1,2,3,4,5,9,11,90,134];
%    noise_intervals = intervalTimes(noise_indx);
%    calls_indx = (150:170);
%    calls_intervals = intervalTimes(calls_indx);
%    histogram(noise_intervals)
%    ylabel("counts")
%    xlabel("interval (ms)")
%    hold on
%    histogram(calls_intervals)
   
   %looks like 50-100ms is the sweet spot
%    histogram(intervalTimes, [0:10:1000])
%    ylabel("counts")
%    xlabel("interval (ms)")

    load('CallTimes.mat', 'callTimes', 'piezo_envelope', 'noise', 'samplingFreq', 'AD_count_double')
    callTimes = broadenCallLength(callTimes, piezo_envelope, noise);
    
    % Call I was looking at earlier is at index 76 for reference
    % Ask Julie about call 61
    for ii = 1:length(callTimes)
        call = callTimes{ii};
        % Calculate the dynamic pitch saliency of the sound event
        [sal, t] = salEstimator(AD_count_double(call(1)*50:call(2)*50), samplingFreq);
        title(strcat(int2str(ii), '/', int2str(length(callTimes))))
        % Calculate the spectrum of the sound event
        figure(8)
        [p,f] = pspectrum(AD_count_double(call(1)*50:call(2)*50), samplingFreq, 'Leakage', 0.85);
        disp(size(p))
        plot(f, pow2db(p), 'LineWidth',2,'k-')
        xlabel('Frequency (Hz)')
        ylabel('Power (dB)')
        % Find quartile power
        cum_power = cumsum(p);
        tot_power = sum(p);
        quartile_freq = zeros(1,3);
        quartile_values = [0.25, 0.5, 0.75];
        nfreqs = length(cum_power);
        iq = 1;
        for ifreq=1:nfreqs
            if (cum_power(ifreq) > quartile_values(iq)*tot_power)
                quartile_freq(iq) = ifreq;
                iq = iq+1;
                if (iq > 3)
                    break;
                end
            end
        end
        hold on
        YLimFig8 = get(gca, 'YLim');
        for iq=1:length(quartile_freq)
            plot(quartile_freq(iq)*ones(2,1), YLimFig8, 'g--')
            hold on\
        end
        MeanSpect = (p/sum(p)).*f/length(p);
        plot(MeanSpect*ones(2,1), YLimFig8, 'r--')
        pause()
        close all
    end
%     disp([fund, sal])
    
    
    
    
%     cut_calls(data_directory, callTimes)
    
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

