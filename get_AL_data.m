function [cut_call_data,tsData,audio2nlg,AL_info,audio_dir,logger_dir,logger_nums] = get_AL_data(expDate,expType,varargin)


pnames = {'sessionType', 'boxNum'};
dflts  = {'communication', NaN};
[sessionType,boxNum] = internal.stats.parseArgs(pnames,dflts,varargin{:});

remote_drive_letter = get_server_letter;
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
bat_num_idx = contains(T.Properties.VariableNames,'Bat_');
logger_nums = T{1,AL_idx};
batNums = T{1,bat_num_idx};

logger_nums = logger_nums(~isnan(logger_nums));
nLogger = length(logger_nums);
% logger_nums = setdiff(logger_nums,cast(T.malfunction_loggers{:},'double'));

logger_dir = fullfile(logger_base_dir,dateStr,'audiologgers');
tsData = cell(1,nLogger);
missing_loggers = false(1,nLogger);
for logger_k = 1:nLogger
    data_fname = dir(fullfile(logger_dir, ['logger' num2str(logger_nums(logger_k))],'extracted_data','*CSC0.mat'));
    if ~isempty(data_fname)
        tsData{logger_k} = matfile(fullfile(data_fname.folder,data_fname.name),'Writable',false);
    else
        missing_loggers(logger_k) = true;
        fprintf('Could not find data for logger #%d\n',logger_nums(logger_k));
    end
end

if any(missing_loggers)
   fieldNames = fieldnames(tsData{find(~missing_loggers,1)});
   empty_ts_data_struct = cell(length(fieldNames),1);
   empty_ts_data_struct = cell2struct(empty_ts_data_struct,fieldNames);
   for logger_k = find(missing_loggers)
       tsData{logger_k} = empty_ts_data_struct;
       tsData{logger_k}.logger_serial_number = num2str(logger_nums(logger_k));
       tsData{logger_k}.Bat_id = num2str(batNums(logger_k));
       tsData{logger_k}.Indices_of_first_and_last_samples = zeros(0,2);
       tsData{logger_k}.logger_type = 'Audi';
   end
end

end