function [call_bat_IDs, logger_power_bands] = get_bat_call_ID(cut_call_data,tsData,audio2nlg)
al_fs = 50e3;
nLogger = length(tsData);
bc = bout_call_data;
bc.min_calls_in_bout = 1;
bc.min_bout_length = 0;

% cut_call_data = cut_call_data(~[cut_call_data.noise]);
callPos = vertcat(cut_call_data.corrected_callpos);
file_callPos = vertcat(cut_call_data.callpos);

nCalls = size(file_callPos,1);

logger_power_nfft = 256;

logger_power_freqs = linspace(1,10e3,logger_power_nfft);
f_bounds = [1 5;5 10]*1e3;
n_freq_band = 2;
f_idx = cell(1,n_freq_band);
for band_k = 1:n_freq_band
    [~, f_idx{band_k}] = inRange(logger_power_freqs,f_bounds(band_k,:));
end
call_bat_IDs = cell(1,nCalls);
logger_power_bands = zeros(n_freq_band,nLogger,nCalls);
out_of_bounds = false;
for k = 1:nCalls
    
        call_ts = callPos(k,:);
        call_ts_nlg = (call_ts + audio2nlg.first_nlg_pulse_time)*1e3;
        
        for logger_k = 1:nLogger
            call_ts_nlg_sample = get_voltage_samples_for_Nlg_timestamps(call_ts_nlg(:),...
                tsData(logger_k).Indices_of_first_and_last_samples(:,1)',...
                tsData(logger_k).Timestamps_of_first_samples_usec,...
                tsData(logger_k).Sampling_period_usec_Logger);
            if isempty(call_ts_nlg_sample)
                call_bat_IDs{k} = NaN;
                break
            end
            piezo_bout_file_callpos = reshape(call_ts_nlg_sample,[],2);
            
            piezo_sample_idx = piezo_bout_file_callpos(1):piezo_bout_file_callpos(2);
            
            if any(piezo_sample_idx<0 | piezo_sample_idx>length(tsData(logger_k).AD_count_int16))
                out_of_bounds = true;
                break
            end
            
            for band_k = 1:n_freq_band
                logger_power_bands(band_k,logger_k,k) = bandpower(double(tsData(logger_k).AD_count_int16(piezo_sample_idx)),al_fs,f_bounds(band_k,:));
            end
        end
        if ~out_of_bounds
            [~,logger_identity_idx] = max(squeeze(logger_power_bands(1,:,k)./logger_power_bands(2,:,k)));
            call_bat_IDs{k} = tsData(logger_identity_idx).Bat_id;
        else
            call_bat_IDs{k} = NaN;
            out_of_bounds = false; 
        end
end