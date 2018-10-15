function [logger_power_bands,classMat,usableIdx] = generate_classification_mat_AL(expDate)

[cut_call_data,tsData,audio2nlg] = get_AL_data(expDate);

[logger_power_bands,success] = get_AL_power_bands(cut_call_data,tsData,audio2nlg,false);
%%
if length(logger_nums) < max_n_logger
    nan_data_append = nan(size(logger_power_bands,1),max_n_logger - length(logger_nums),size(logger_power_bands,3));
    logger_power_bands = cat(2,logger_power_bands,nan_data_append);
end
%%
lowPower = squeeze(log10(max(logger_power_bands(1,:,:),[],2)));
highPower = squeeze(log10(max(logger_power_bands(2,:,:),[],2)));
powerRatio = lowPower./highPower;
powerDifference = lowPower - highPower;
classMat = table(lowPower,highPower,powerRatio,powerDifference);
usableIdx = success' & ~any(isnan(table2array(classMat)),2) & ~any(isinf(table2array(classMat)),2);
classMat = classMat{usableIdx,:};