function [cut_call_data,tsData,audio2nlg,AL_info,audio_dir,logger_dir,logger_nums] = get_AL_data(expDate,expType,varargin)


pnames = {'sessionType', 'boxNum'};
dflts  = {'communication', NaN};
[sessionType,boxNum] = internal.stats.parseArgs(pnames,dflts,varargin{:});

addpath('C:\Users\phyllo\Documents\GitHub\SoundAnalysisBats\')

switch expType
    case 'acousticRecording'
        dataDir = 'C:\Users\phyllo\Documents\Maimon\acoustic_recording\';
        logger_base_dir = 'Z:\users\Maimon\acoustic_recording\video_and_AL\';
        audio_base_dir = 'Z:\users\Maimon\acoustic_recording\audio\';
        audio_dir_str = 'audio\ch1\';
        T = readtable([dataDir 'recording_logs.csv']);
        T = T(T.Date==expDate,:);
    case 'adultRecording'
        dataDir = 'Y:\users\maimon\adult_recording\';
        logger_base_dir = dataDir;
        audio_base_dir = dataDir;
        audio_dir_str = 'audio\communication\ch1\';
        T = readtable([dataDir 'documents\recording_logs.csv']);
        dateIdx = T.Date==expDate & strcmp(T.Session,sessionType);
        T = T(dateIdx,:);
    case 'operant_adult_recording'
        dataDir = 'Y:\users\maimon\adult_operant_recording\';
        logger_base_dir = dataDir;
        audio_base_dir = dataDir;
        T = readtable([dataDir 'documents\recording_logs.csv']);
        switch sessionType
            case 'communication'
                audio_dir_str = 'audio\communication\ch1';
                dateIdx = T.Date==expDate & strcmp(T.Session,sessionType);
            case 'operant'
                audio_dir_str = ['operant\box' boxNum];
                dateIdx = T.Date==expDate & strcmp(T.Session,sessionType) & T.Box == str2double(boxNum);
        end
        
        T = T(dateIdx,:);
        
end

date_str_format = 'mmddyyyy';

dateStr = datestr(expDate,date_str_format);

audio_dir = fullfile(audio_base_dir,dateStr,audio_dir_str);
audio2nlg = load(fullfile(audio_dir,'audio2nlg_fit.mat'));
cut_call_data = load(fullfile(audio_dir,'cut_call_data.mat'));
cut_call_data = cut_call_data.cut_call_data;

if isfile(fullfile(audio_dir,'AL_class_info.mat'))
    AL_info = load(fullfile(audio_dir,'AL_class_info.mat'));
else
    AL_info = [];
end
AL_idx = contains(T.Properties.VariableNames,'AL');
logger_nums = T{1,AL_idx};

logger_nums = logger_nums(~isnan(logger_nums));
% logger_nums = setdiff(logger_nums,cast(T.malfunction_loggers{:},'double'));

logger_dir = [logger_base_dir dateStr '\audiologgers\'];
for logger_k = 1:length(logger_nums)
    data_fname = dir([logger_dir 'logger' num2str(logger_nums(logger_k)) filesep 'extracted_data' filesep '*CSC0.mat']);
    tsData(logger_k) = load(fullfile(data_fname.folder,data_fname.name));
end

end