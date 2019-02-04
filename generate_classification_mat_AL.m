function [logger_power_bands,logger_xcorr,classMat,batNums,usableIdx] = generate_classification_mat_AL(cut_call_data,tsData,audio2nlg)

wav_fs = 250e3;
[logger_power_bands,logger_xcorr,success] = get_AL_stats(cut_call_data,tsData,audio2nlg,false);
nLogger = size(logger_power_bands,2);
batNums = cell(1,nLogger);
max_n_logger = 6;

if nLogger < max_n_logger
    nan_data_append = nan(size(logger_power_bands,1),max_n_logger - nLogger,size(logger_power_bands,3));
    logger_power_bands = cat(2,logger_power_bands,nan_data_append);
end


classMat = cell(1,nLogger);
usableIdx = success';
callRMS = cellfun(@rms,{cut_call_data.cut})';
callLength = cellfun(@length,{cut_call_data.cut})'/wav_fs;
for k = 1:nLogger
    batNums{k} = tsData(k).Bat_id;
    lowPower = squeeze(log10(logger_power_bands(1,k,:)));
    highPower = squeeze(log10(logger_power_bands(2,k,:)));
    powerRatio = lowPower./highPower;
    powerDifference = lowPower - highPower;
    crossCorr = logger_xcorr(k,:)';
    classMat{k} = table(lowPower,highPower,powerRatio,powerDifference,crossCorr);
    usableIdx = usableIdx & ~any(isnan(table2array(classMat{k})),2) & ~any(isinf(table2array(classMat{k})),2);
end


callRMS = callRMS(usableIdx);
callLength = callLength(usableIdx);

for k = 1:nLogger
    classMat{k} = classMat{k}(usableIdx,:);
    classMat{k} = addvars(classMat{k},callLength,callRMS,'NewVariableNames',{'callLength','callRMS'});
end

end
