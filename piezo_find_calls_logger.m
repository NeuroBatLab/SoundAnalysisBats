function [] = piezo_find_calls_logger(directory)
    % We will see what this ends up doing...
    
    RMSfactor = 10;
    callLength = 0.015; % in s
    mergeThresh = 5e-3; % in s
    Fhigh_power =50; % Frequency upper bound for calculating the envelope (time running RMS)
    FS_env = 1000; % Sample frequency of the envelope
    BandPassFilter = [1000 5000]; %the frequencies we care about to identiy when a call is made


    % AD_count_int16 is the raw_data variable in the extracted CSC.mat file
    % Esimated_channelFS_Transceiver is the estimated sample frequencies
    load(strcat(directory, '/extracted_data/71306_20180907_CSC0.mat'), 'AD_count_int16', 'Estimated_channelFS_Transceiver')
    samplingFreq = mean(Estimated_channelFS_Transceiver);
    recordingLength = length(AD_count_int16) / samplingFreq / 60;  % in minutes
    sampleSize = 2; % in minutes



    % Center the signal and filter it
    centered_piezo_signal = AD_count_int16 - mean(AD_count_int16);
    centered_piezo_signal = double(centered_piezo_signal);

    %Design the bandpass filter for the 1000-5000Hz range
    [z,p,k] = butter(6,BandPassFilter / (samplingFreq / 2),'bandpass');
    sos_low = zp2sos(z,p,k);    

    
    size = length(centered_piezo_signal);
    piezo_envelope = zeros(1, floor(size / 49)); %note this is hardcoded
    envelope_length = 0;
    for i = 1:4
        tic;
        sampleStart = 1;
        sampleEnd = floor(size / 4);
        if i == 4
            sampleEnd = length(centered_piezo_signal);
        end
        sample = centered_piezo_signal(sampleStart:sampleEnd);
        centered_piezo_signal(sampleStart: sampleEnd) = [];
        filtered_piezo_sample = filtfilt(sos_low, 1, sample);
        filtered_sample_envelope = envelope(filtered_piezo_sample, 1e-3 * 50000, 'rms');
        filtered_sample_envelope = resample(filtered_sample_envelope, 1, 50);
        piezo_envelope(envelope_length + 1 : envelope_length + length(filtered_sample_envelope)) = filtered_sample_envelope;
        envelope_length = envelope_length + length(filtered_sample_envelope);
        toc;
    end
    piezo_envelope = piezo_envelope(1:envelope_length);
    clear centered_piezo_signal
    
    
    %try graphing noise as a mixture model
    noise = getAverageNoise(piezo_envelope, FS_env, 50);
    disp(['Done with calculating noise: ', num2str(noise)])
    %still varies more than I like (6-9 ish), could increase sample size
    
%     tic;            
    plot((1:length(piezo_envelope)), piezo_envelope, 'Color',[0,0.5,0.9])
    line([1, length(piezo_envelope)], [noise * RMSfactor, noise * RMSfactor], 'Color','red','LineStyle','--')
    hold on
 
    %List of a bunch of 1s and 0s indicating calls
    callIndicator = piezo_envelope > (RMSfactor * noise);
        
    %Start/Stop times is any time there is a change
    startTimes = find(diff(callIndicator) > 0);
    stopTimes = find(diff(callIndicator) < 0);
                
    %check edge cases
    if stopTimes(1) < startTimes(1)
        startTimes = [1, startTimes];
    end
        
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
        
    %Find the call times, the requirements of which are that the call 
    %is greater than callLength and also merge any calls that are less 
    %than mergethresh apart
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
%     toc;

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
    average_noise = 0;
    for i = 1:num_samples
        average_noise = average_noise + getNoiseLevel(piezo_envelope, FS_env);
    end
    
    average_noise = average_noise / num_samples;

end

function [average_noise] = getNoiseLevel(piezo_envelope, FS_env)
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
    rand_indx = rand * length(piezo_envelope);
    while (rand_indx + sampleLength) > length(piezo_envelope)
        rand_indx = rand * length(piezo_envelope);
    end
end

