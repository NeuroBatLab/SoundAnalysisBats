function [] = piezo_find_calls_logger_optimization(Data_directory)
    % Takes in the directory for the individidual logger data
    % LOGGER_DIRECTORY and outputs SoundEvent_LoggerSamp, the start/stop indices of potential
    % calls in the raw input data. Note that SoundEvent_LoggerSamp is a matrix where each
    % row represents a sound event with the first column corresponding to onsets and
    % the second column to offsets. This function should be used in combination with
    % piezo_find_calls.m, which will call this function on all loggers for 
    % a recording session. 
    
    %% Setting up variables that will be used later
    DurChuncksList = 1:1:10; % Duration of chuncks we want to test in min. 10 min ensures that there is no significant time drift between the envelope calculation and the original sound
    Fhigh_power =75; % Frequency upper bound for calculating the envelope (time running RMS)
    FS_env = 1000; % Sample frequency of the envelope
    BandPassFilter = [1000 5000]; % the frequencies we care about to identify when a call is made
    PathPieces = split(Data_directory, filesep);
    logger_name = PathPieces{find(contains(PathPieces,'extracted_data'))-1};
    disp(logger_name)
    
    %it's saying the last call occurs at 3.8731, which doesn't make sense.
    %Check via graph
    
   
    %% Load the raw signal and raw signal frequency, then calculate ratio of
    % raw signal frequency to what the envelope's frequency will be.
    file = dir(fullfile(Data_directory, '*CSC0*'));
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
    
    %% Design the bandpass filter for the 1000-5000Hz range
    [z,p,k] = butter(6,BandPassFilter / (SamplingFreq / 2),'bandpass');
    sos_low = zp2sos(z,p,k);    

    %% Loop through the different chunck sizes calculating speed of calculation for running rms
    Env1CalcStop=nan(length(DurChuncksList),1);
    NChuncks1 = nan(length(DurChuncksList),1);
    for cc=1:length(DurChuncksList)
        DurChunck = DurChuncksList(cc);
        
        %% Cut the signal into chuncks of same short duration
        Signal_length = length(centered_piezo_signal);
        NumSampleChunck = round(DurChunck*60*SamplingFreq);
        Chuncks = [1:NumSampleChunck:Signal_length Signal_length];
        NChuncks1(cc) = length(Chuncks)-1;
        piezo_centered_signal = cell(1,NChuncks1(cc)); % 

        for ii = 1:NChuncks1(cc)
            piezo_centered_signal{ii} = centered_piezo_signal(Chuncks(ii):(Chuncks(ii+1)));
        end
        
        %% Calculate the amplitude envelope independantly on every chunck
        piezo_envelope = cell(1,NChuncks1(cc)); %
        NChuncks_local = NChuncks1(cc);
        EnvCalcStart = tic;
        parfor ii = 1:NChuncks_local %%% parallelize
            fprintf(1,'Start Envelope chunck %d/%d\n', ii, NChuncks_local)
            LocalStart = tic;
            filtered_piezo_sample = filtfilt(sos_low, 1, piezo_centered_signal{ii});
            piezo_envelope{ii} = running_rms(filtered_piezo_sample,SamplingFreq, Fhigh_power,FS_env);
            fprintf(1,'Done Envelope chunck %d/%d in %.1fs\n', ii, NChuncks_local, toc(LocalStart))
        end
        Env1CalcStop(cc) = toc(EnvCalcStart);
        fprintf(1,'Done Calculating Envelope with running_rms in %.1fs\n', Env1CalcStop)
    end
    
    %% Loop through the different chunck sizes calculating speed of calculation for envelope
    Env2CalcStop=nan(length(DurChuncksList),1);
    NChuncks2 = nan(length(DurChuncksList),1);
    for cc=1:length(DurChuncksList)
        DurChunck = DurChuncksList(cc);
        
        %% Cut the signal into chuncks of same short duration
        Signal_length = length(centered_piezo_signal);
        NumSampleChunck = round(DurChunck*60*SamplingFreq);
        Chuncks = [1:NumSampleChunck:Signal_length Signal_length];
        NChuncks2(cc) = length(Chuncks)-1;
        piezo_centered_signal = cell(1,NChuncks2(cc)); % 

        for ii = 1:NChuncks2(cc)
            piezo_centered_signal{ii} = centered_piezo_signal(Chuncks(ii):(Chuncks(ii+1)));
        end
        
        %% Calculate the amplitude envelope independantly on every chunck
        piezo_envelope2 = cell(1,NChuncks2(cc)); %
        EnvCalcStart = tic;
        NChuncks_local = NChuncks2(cc);
        parfor ii = 1:NChuncks_local %%% parallelize
            fprintf(1,'Start Envelope chunck %d/%d\n', ii, NChuncks_local)
            LocalStart = tic;
            filtered_piezo_sample = filtfilt(sos_low, 1, piezo_centered_signal{ii});
            filtered_sample_envelope = envelope(filtered_piezo_sample, round(1e-3 * SamplingFreq), 'rms');
            piezo_envelope2{ii} = resample(filtered_sample_envelope, 1, FS_ratio);
            fprintf(1,'Done Envelope chunck %d/%d in %.1fs\n', ii, NChuncks_local, toc(LocalStart))
        end
        Env2CalcStop(cc) = toc(EnvCalcStart);
        fprintf(1,'Done Calculating Envelope with enveloppe in %.1fs\n', Env2CalcStop)
    end
    clear centered_piezo_signal
    
    figure()
    subplot(2,1,1)
    plot(DurChuncksList,Env1CalcStop*SamplingFreq/Signal_length,'k','LineWidth',2)
    hold on
    plot(DurChuncksList,Env2CalcStop*SamplingFreq/Signal_length,':r','LineWidth',2)
    legend({'running rms' 'envelope'})
    xlabel('Duration of chuncks (min)')
    ylabel('Speed (s of calculation/s of signal)')
    title('Envelope calculation speed optimization')
    hold off
    subplot(2,1,2)
    plot(NChuncks1,Env1CalcStop*SamplingFreq/Signal_length,'k','LineWidth',2)
    hold on
    plot(NChuncks2,Env2CalcStop*SamplingFreq/Signal_length,':r','LineWidth',2)
    legend({'running rms' 'envelope'})
    xlabel('Number of loops')
    ylabel('Speed (s of calculation/s of signal)')
    title(sprintf('Mean Running rms = %dms Mean Envelope = %dms',round(nanmean(Env1CalcStop*SamplingFreq/Signal_length)*10^3),round(nanmean(Env2CalcStop*SamplingFreq/Signal_length)*10^3)))
    hold off
end
