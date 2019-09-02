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
    RMSfactor = 10; % how much greater the call's envelope needs to be than the noise. Subject to change implementation
    callLength = 0.015; % minimum length of a call in s
    mergeThresh = 5e-3; % minimum length between two calls in s
    FS_env = 1000; % Sample frequency of the envelope
    BandPassFilter = [1000 5000]; % the frequencies we care about to identiy when a call is made
    logger_name = data_directory(end - 23 : end - 16); % I'm assuming all names of the form loggerXX
    disp(logger_name)
    
   
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
           samplingFreq = mean(Estimated_channelFS_Transceiver);
           FS_ratio = round(samplingFreq / FS_env); 
           found = 1;
       end
    end
    
    if found == 0
        ME = MException('Data file not found');
        throw(ME)
    end    

    % Center the signal and clear the old data from memory
    centered_piezo_signal = AD_count_int16 - mean(AD_count_int16);
    centered_piezo_signal = double(centered_piezo_signal);
    clear AD_count_int16

    %Design the bandpass filter for the 1000-5000Hz range
    [z,p,k] = butter(6,BandPassFilter / (samplingFreq / 2),'bandpass');
    sos_low = zp2sos(z,p,k);    

    % Filter the signal and compute the envelope (using the rms). Once
    % again, remove old data from memory
    signal_length = length(centered_piezo_signal);
    piezo_envelope = zeros(1, floor(signal_length / (FS_ratio - 1)));
    envelope_length = 0;
    for i = 1:4
        tic;
        sampleStart = 1;
        sampleEnd = floor(signal_length / 4);
        if i == 4
            sampleEnd = length(centered_piezo_signal);
        end
        sample = centered_piezo_signal(sampleStart:sampleEnd);
        centered_piezo_signal(sampleStart: sampleEnd) = [];
        filtered_piezo_sample = filtfilt(sos_low, 1, sample);
        filtered_sample_envelope = envelope(filtered_piezo_sample, 1e-3 * 50000, 'rms');
%         disp(length(filtered_sample_envelope))
        %can likely get rid of, envelope window length takes care of this
        filtered_sample_envelope = resample(filtered_sample_envelope, 1, 50);
%         disp(length(filtered_sample_envelope))
        piezo_envelope(envelope_length + 1 : envelope_length + length(filtered_sample_envelope)) = filtered_sample_envelope;
        envelope_length = envelope_length + length(filtered_sample_envelope);
        toc;
    end
    piezo_envelope = piezo_envelope(1:envelope_length);
    clear centered_piezo_signal
    
    %Find the noise for the given logger
    noise = getAverageNoise(piezo_envelope, FS_env, 50);
    disp(['Done with calculating noise: ', num2str(noise)])
    
    %Display the nosie threshold
    plot((1:length(piezo_envelope)), piezo_envelope, 'Color',[0,0.5,0.9])
    line([1, length(piezo_envelope)], [noise * RMSfactor, noise * RMSfactor], 'Color','red','LineStyle','--')
    title(logger_name)
    hold on
 
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
        
    %visually inspect that the previous step is correct
    x_values = [startTimes, stopTimes];
    y_values = repmat((100),1,length(x_values));
    scatter(x_values, y_values) 
    hold on
        

    if length(startTimes) ~= length(stopTimes)
        ME = MException('Start and Stop times are not the same size');
        throw(ME) 
    end
        
    % Find the call times, the requirements of which are that the call time 
    % is greater than callLength and also merge any calls that are less 
    % than mergethresh apart
    callTimes = cell(1, length(startTimes));
    index = 1;
    for i = 1:length(startTimes)
        start = startTimes(i);
        stop = stopTimes(i);

        if start > stop
            ME = MException('Start time should not be greater than Stop time');
            throw(ME)
        end

        if stop - start >= FS_env * callLength
            if index ~= 1
                %check to see if we can merge with last call 
                %just this start - previous stop < mergethresh * FS_env
                previous_call = callTimes{index - 1};
                if start - previous_call(2) <= mergeThresh * FS_env
                    callTimes{end} = [previous_call(1), stop];
                else
                    callTimes{index} = [start, stop];
                    index = index + 1;
                end
            else
                callTimes{index} = [start, stop];
                index = index + 1;
            end
        end     
    end
    callTimes = callTimes(1 : index - 1);

    %visually inspect that the previous step is correct
    for i = 1:length(callTimes)
        call = callTimes{i};
        x_start = call(1);
        x_stop = call(2);
        line([x_start, x_stop], [1000,1000], 'Color','g','Marker', '+')
        hold on
    end
    hold off
    
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

