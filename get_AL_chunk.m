function [call_logger_chunk,success] = get_AL_chunk(tsData,call_ts,audio2nlg)
success = true;
call_logger_chunk = [];
call_ts_nlg = (call_ts + audio2nlg.first_nlg_pulse_time)*1e3;
call_ts_nlg_sample = get_voltage_samples_for_Nlg_timestamps(call_ts_nlg(:),...
    tsData.Indices_of_first_and_last_samples(:,1)',...
    tsData.Timestamps_of_first_samples_usec,...
    tsData.Samples_per_channel_per_file,...
    1e6/nanmean(tsData.Estimated_channelFS_Transceiver));
if isempty(call_ts_nlg_sample)
    success = false;
end
piezo_bout_file_callpos = reshape(call_ts_nlg_sample,[],2);

piezo_sample_idx = piezo_bout_file_callpos(1):piezo_bout_file_callpos(2);

if any(piezo_sample_idx<0 | piezo_sample_idx>length(tsData.AD_count_int16) | isnan(piezo_sample_idx))
    success = false;
    return
end
call_logger_chunk = double(tsData.AD_count_int16(piezo_sample_idx));
end