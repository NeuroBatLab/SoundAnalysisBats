function [cut_call_data,tsData,audio2nlg,AL_info,audio_dir,logger_dir,logger_nums] = get_AL_data(expDate,expType,varargin)


pnames = {'sessionType', 'boxNum'};
dflts  = {'communication', NaN};
[sessionType,boxNum] = internal.stats.parseArgs(pnames,dflts,varargin{:});

% addpath('C:\Users\phyllo\Documents\GitHub\SoundAnalysisBats\')
[~,net_use_str]= dos('net use');
net_use_str = strsplit(net_use_str,'\n');

net_use_str = net_use_str{contains(net_use_str,'yartsev_server3')};
remote_drive_letter = regexp(net_use_str,'[A-Z]{1}\:','match');
remote_drive_letter = remote_drive_letter{1};
remote_base_dir = fullfile(remote_drive_letter,'users\maimon');

switch expType
    case 'acousticRecording'
        dataDir = 'C:\Users\phyllo\Documents\Maimon\acoustic_recording\';
        logger_base_dir = fullfile(remote_base_dir,'\acoustic_recording\video_and_AL\');
        audio_base_dir = fullfile(remote_base_dir,'acoustic_recording\audio\');
        audio_dir_str = 'audio\ch1\';
        T = readtable([dataDir 'recording_logs.csv']);
        T = T(T.Date==expDate,:);
    case 'adultRecording'
        dataDir = fullfile(remote_base_dir,'adult_recording');
        logger_base_dir = dataDir;
        audio_base_dir = dataDir;
        audio_dir_str = 'audio\communication\ch1\';
        T = readtable([dataDir 'documents\recording_logs.csv']);
        dateIdx = T.Date==expDate & strcmp(T.Session,sessionType);
        T = T(dateIdx,:);
    case 'adult_social'
        dataDir = fullfile(remote_base_dir,'adult_social_recording');
        logger_base_dir = dataDir;
        audio_base_dir = dataDir;
        audio_dir_str = ['audio\' sessionType '\ch1\'];
        T = get_rec_logs;
        dateIdx = T.Date==expDate & strcmp(T.Session,sessionType);
        T = T(dateIdx,:);
    case 'operant_adult_recording'
        dataDir = fullfile(remote_base_dir,'adult_operant_recording');
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

logger_dir = fullfile(logger_base_dir,dateStr,'audiologgers');
for logger_k = 1:length(logger_nums)
    data_fname = dir(fullfile(logger_dir, ['logger' num2str(logger_nums(logger_k))],'extracted_data','*CSC0.mat'));
    tsData(logger_k) = load(fullfile(data_fname.folder,data_fname.name));
end

end