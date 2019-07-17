function [] = piezo_find_calls(directory)
    % We will see what this ends up doing...
    RMSfactor = 10;
    callLength = 0.015; % in s
    mergethresh = 5e-3; % in s
    Fhigh_power =50; % Frequency upper bound for calculating the envelope (time running RMS)
    FS_env = 1000; % Sample frequency of the envelope
    BandPassFilter = [1000 5000]; %the frequencies we care about to identiy when a call is made


    % AD_count_int16 is the raw_data variable in the extracted CSC.mat file
    % Esimated_channelFS_Transceiver is the estimated sample frequencies
    load(strcat(directory, '/logger12/extracted_data/71306_20180907_CSC0.mat'), 'AD_count_int16', 'Estimated_channelFS_Transceiver')
    samplingFreq = mean(Estimated_channelFS_Transceiver);
    recordingLength = length(AD_count_int16) / samplingFreq / 60;  % in minutes
    sampleSize = 2; % in minutes



    % Center the signal and filter it
    centered_piezo_signal = AD_count_int16 - mean(AD_count_int16);
    centered_piezo_signal = double(centered_piezo_signal);

    %Design the bandpass filter for the 1000-5000Hz range
    [z,p,k] = butter(6,BandPassFilter / (samplingFreq / 2),'bandpass');
    sos_low = zp2sos(z,p,k);

    %get the baseline noise for the RMS of the piezo data
    noise = getNoiseLevel(centered_piezo_signal, samplingFreq, sos_low, Fhigh_power, FS_env);
    disp(['Done with calculating noise: ', num2str(noise)])
    %still varies more than I like (6-9 ish), could increase sample size

    %loop after designing filter but before applying it
    %seems to be about 15 seconds per 2 minute sample to filter and calculate
    %the envelope. 

    for j = 0:sampleSize: 2%recordingLength
        sampleStart = uint64(j * 60 * samplingFreq + 1);
        sampleEnd = uint64((j + sampleSize) * 60 * samplingFreq);
        if sampleEnd > length(centered_piezo_signal)
           sampleEnd = length(centered_piezo_signal);
        end
        sample = centered_piezo_signal(sampleStart : sampleEnd);
        
        filtered_piezo_sample = filtfilt(sos_low, 1, sample);
        piezo_envelope = running_rms(filtered_piezo_sample, samplingFreq, Fhigh_power, FS_env);
        

        x_start = j * (length(piezo_envelope) / sampleSize) + 1;
        x_end = (j + sampleSize) * (length(piezo_envelope) / sampleSize);
        x_values = (x_start : x_end);
        plot(x_values, piezo_envelope, 'Color',[0,0.5,0.9])
        line([x_start, x_end], [noise * RMSfactor, noise * RMSfactor], 'Color','red','LineStyle','--')
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
        
%         %visually inspect that the previous step is correct
%         x_values = [startTimes + j * (length(piezo_envelope) / sampleSize), stopTimes + j * (length(piezo_envelope) / sampleSize)];
%         y_values = repmat((100),1,length(x_values));
%         scatter(x_values, y_values) 
%         hold on
        

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
                    if start - previous_call(2) <= mergethresh * FS_env
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
            x_start = call(1) + j * (length(piezo_envelope) / sampleSize);
            x_stop = call(2) + j * (length(piezo_envelope) / sampleSize);
            line([x_start, x_stop], [1000,1000], 'Color','g','LineWidth', 5.0)
            hold on
        end

    end
    hold off
 
end

function [average_noise] = getNoiseLevel(centered_piezo_signal, samplingFreq, filter, Fhigh_power, FS_env)
    rand_indx = uint64(getRandomIndex(centered_piezo_signal, samplingFreq));    
    repeat = true;
    while repeat
        total = 0;
        repeat = false;
        sample = centered_piezo_signal(rand_indx: rand_indx + samplingFreq);
        filtered_piezo_sample = filtfilt(filter, 1, sample);
        piezo_envelope = running_rms(filtered_piezo_sample, samplingFreq, Fhigh_power, FS_env);
        if any(piezo_envelope >50)
            rand_indx = uint64(getRandomIndex(centered_piezo_signal, samplingFreq));
            repeat = true;
        else
            total = sum(piezo_envelope);
        end
    end
    
    average_noise = total / length(piezo_envelope);
end

function [rand_indx] = getRandomIndex(centered_piezo_signal, samplingFreq)
    rand_indx = rand * length(centered_piezo_signal);
    while (rand_indx + samplingFreq) > length(centered_piezo_signal)
        rand_indx = rand * length(centered_piezo_signal);
    end
end

