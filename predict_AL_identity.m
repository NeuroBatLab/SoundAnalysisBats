function logger_ID = predict_AL_identity(AL_info)

s = load('C:\Users\phyllo\Documents\Maimon\acoustic_recording\analysis_results\call_noise_prediction_mdl_ALs.mat');
call_noise_mdl = s.call_noise_prediction_mdl_ALs;

noise_pred = nan(size(AL_info.classMat{1},1),length(AL_info.classMat));
for k = 1:length(AL_info.classMat)
    noise_pred(:,k) = call_noise_mdl.predictFcn(AL_info.classMat{k});
end

[ID_rows,ID_col] = find(~noise_pred);
logger_ID = cell(1,size(noise_pred,1));
for k = 1:size(noise_pred,1)
    logger_ID{k} = ID_col(ID_rows==k);
end


end