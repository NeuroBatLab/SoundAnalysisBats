function [logger_power_bands,logger_xcorr,success] = get_AL_stats(cut_call_data,tsData,audio2nlg,manual_noise_flag)
addpath('C:\Users\phyllo\Documents\GitHub\SoundAnalysisBats\')
if nargin < 4
    manual_noise_flag = false;
end    

al_fs = 50e3;
wav_fs = 250e3;
ds_factor = wav_fs/al_fs;
envelope_window_s = 10e-3;
nLogger = length(tsData);

if manual_noise_flag
    cut_call_data = cut_call_data(~[cut_call_data.noise]);
end
callPos = vertcat(cut_call_data.corrected_callpos);
nCalls = size(callPos,1);

logger_power_nfft = 256;

logger_power_freqs = linspace(1,10e3,logger_power_nfft);
f_bounds = [1 5;5 10]*1e3;
n_freq_band = 2;
f_idx = cell(1,n_freq_band);
for band_k = 1:n_freq_band
    [~, f_idx{band_k}] = inRange(logger_power_freqs,f_bounds(band_k,:));
end
[b,a]=butter(3,f_bounds(1,:)/(50e3/2),'bandpass');
success = true(1,nCalls);
logger_power_bands = nan(n_freq_band,nLogger,nCalls);
logger_xcorr = nan(nLogger,nCalls);
out_of_bounds = false;
for k = 1:nCalls
    
        call_ts = callPos(k,:);
        call_ts_nlg = (call_ts + audio2nlg.first_nlg_pulse_time)*1e3;
        
        for logger_k = 1:nLogger
            call_ts_nlg_sample = get_voltage_samples_for_Nlg_timestamps(call_ts_nlg(:),...
                tsData(logger_k).Indices_of_first_and_last_samples(:,1)',...
                tsData(logger_k).Timestamps_of_first_samples_usec,...
                tsData(logger_k).Samples_per_channel_per_file,...
                1e6/nanmean(tsData(logger_k).Estimated_channelFS_Transceiver));
            if isempty(call_ts_nlg_sample)
                success(k) = false;
                break
            end
            piezo_bout_file_callpos = reshape(call_ts_nlg_sample,[],2);
            
            piezo_sample_idx = piezo_bout_file_callpos(1):piezo_bout_file_callpos(2);
            
            if any(piezo_sample_idx<0 | piezo_sample_idx>length(tsData(logger_k).AD_count_int16))
                out_of_bounds = true;
                break
            end
            call_logger_chunk = double(tsData(logger_k).AD_count_int16(piezo_sample_idx));
            
            for band_k = 1:n_freq_band
                logger_power_bands(band_k,logger_k,k) = bandpower(call_logger_chunk,al_fs,f_bounds(band_k,:));
            end
            
            call_wav_env = zscore(downsample(envelope(cut_call_data(k).cut,envelope_window_s*wav_fs,'rms'),ds_factor));
            call_al_env = envelope(filtfilt(b,a,call_logger_chunk),envelope_window_s*al_fs,'rms');
            call_wav_env_L = length(call_wav_env);
            call_al_env_L = length(call_al_env);
            N = min(call_al_env_L,call_wav_env_L);

            r = xcorr(call_al_env(1:N),call_wav_env(1:N),1,'unbiased');
            logger_xcorr(logger_k,k) = r(2);
            
        end
        if out_of_bounds
            success(k) = false;
            out_of_bounds = false; 
        end
end