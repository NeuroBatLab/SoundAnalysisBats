function [cut_call_data,tsData,audio2nlg] = get_AL_data(expDate)

addpath('C:\Users\phyllo\Documents\GitHub\SoundAnalysisBats\')
dataDir = 'D:\acoustic_recording\';
T = readtable([dataDir 'recording_logs.csv']);
T = T(T.Date==expDate,:);
logger_base_dir = 'Z:\users\Maimon\acoustic_recording\video_and_AL\';
audio_base_dir = 'Z:\users\Maimon\acoustic_recording\audio\';
date_str_format = 'mmddyyyy';

dateStr = datestr(expDate,date_str_format);

audio_dir = [audio_base_dir dateStr filesep 'audio\ch1\'];
audio2nlg = load([audio_dir 'audio2nlg_fit.mat']);
cut_call_data = load([audio_dir 'cut_call_data.mat']);
cut_call_data = cut_call_data.cut_call_data;

logger_nums = T{1,4:2:14};
logger_nums = logger_nums(~isnan(logger_nums));
logger_nums = setdiff(logger_nums,str2double(T.malfunction_loggers{1}));
loggerDir = [logger_base_dir dateStr '\audiologgers\'];
for logger_k = 1:length(logger_nums)
    data_fname = dir([loggerDir 'logger' num2str(logger_nums(logger_k)) filesep 'extracted_data' filesep '*CSC0.mat']);
    tsData(logger_k) = load(fullfile(data_fname.folder,data_fname.name));
end

end