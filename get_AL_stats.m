function [logger_power_bands,logger_xcorr,success] = get_AL_stats(cut_call_data,tsData,audio2nlg,manual_noise_flag)
addpath('C:\Users\phyllo\Documents\GitHub\SoundAnalysisBats\')
if nargin < 4
    manual_noise_flag = false;
end    

al_fs = round(1/(1e-6*unique([tsData.Sampling_period_usec_Logger])));
wav_fs = unique([cut_call_data.fs]);
[resample_N,resample_D] = rat(wav_fs/al_fs);
envelope_window_s = 10e-3;
nLogger = length(tsData);

if manual_noise_flag
    cut_call_data = cut_call_data(~[cut_call_data.noise]);
end
callPos = vertcat(cut_call_data.corrected_callpos);
nCalls = size(callPos,1);

logger_power_nfft = 256;

logger_power_freqs = linspace(1,10e3,logger_power_nfft);
f_bounds = [1 3;5 10]*1e3;
n_freq_band = 2;
f_idx = cell(1,n_freq_band);
for band_k = 1:n_freq_band
    [~, f_idx{band_k}] = inRange(logger_power_freqs,f_bounds(band_k,:));
end
[b,a]=butter(3,f_bounds(1,:)/(al_fs/2),'bandpass');
success = true(1,nCalls);
logger_power_bands = nan(n_freq_band,nLogger,nCalls);
logger_xcorr = nan(nLogger,nCalls);
for k = 1:nCalls
    
    call_ts = callPos(k,:);
    call_wav_env = zscore(resample(envelope(cut_call_data(k).cut,envelope_window_s*wav_fs,'rms'),resample_D,resample_N));
    call_wav_env_L = length(call_wav_env);
    
    for logger_k = 1:nLogger
        [call_logger_chunk,success(k)] = get_AL_chunk(tsData(logger_k),call_ts,audio2nlg);
        
        if ~success(k)
            continue
        end
        
        for band_k = 1:n_freq_band
            logger_power_bands(band_k,logger_k,k) = bandpower(call_logger_chunk,al_fs,f_bounds(band_k,:));
        end
        
        call_al_env = zscore(envelope(filtfilt(b,a,call_logger_chunk),envelope_window_s*al_fs,'rms'));
        call_al_env_L = length(call_al_env);
        N = min(call_al_env_L,call_wav_env_L);
        
        r = xcorr(call_al_env(1:N),call_wav_env(1:N),1,'unbiased');
        logger_xcorr(logger_k,k) = r(2);
        
    end
end