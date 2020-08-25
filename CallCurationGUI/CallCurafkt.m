function CallCurafkt(action)
global vv Nvoc df redo DataFiles ManCall FhGUI;
global submith noCallh redoh starth oldvv olddf;
global redoEditVoch redoEditSeth checkboxh stopclick;


switch action
    
    case 'Start'
        set(starth,'enable','off')
        initFiles
        loadnextfile(vv)
        
    case 'NoCall'
        disableEvals
        drawnow;
        newmessage('Manual input enforced: Noise (=NoCall)');
        ManCall=0;
        savingData
        if redo
            vv=oldvv;
            df=olddf;
            redo=0;
        end
        vv=vv+1;
        if vv<=Nvoc
            loadnextfile(vv)
        else
            eraseData
            if df<=length(DataFiles)
                grabNewDatafile
                loadnextfile(vv)
            else
                newmessage('Annotation done!');
                initFiles
                loadnextfile(vv)
            end
        end
        
    case 'Redo'
        disableEvals
        drawnow;
        oldvv=vv-1;
        olddf=df;
        %grab which file to reevaluate
        vv=str2num(get(redoEditVoch,'string'));
        df=str2num(get(redoEditSeth,'string'));
        if vv>Nvoc || df>length(DataFiles)
            newmessage('Incorrect Voc# or Set#');
            set([redoh redoEditVoch redoEditSeth],'enable','on')
            df=olddf; vv=oldvv+1;
            set(redoEditVoch,'String',num2str(vv))
            set(redoEditSeth,'String',num2str(df))
        elseif vv>oldvv+1 && df>=olddf
            newmessage('Do not evaluate into the future!')
            set([redoh redoEditVoch redoEditSeth],'enable','on')
            df=olddf; vv=oldvv+1;
            set(redoEditVoch,'String',num2str(vv))
            set(redoEditSeth,'String',num2str(df))
        elseif df>olddf
            newmessage('Do not evaluate into the future!')
            set([redoh redoEditVoch redoEditSeth],'enable','on')
            df=olddf; vv=oldvv+1;
            set(redoEditVoch,'String',num2str(vv))
            set(redoEditSeth,'String',num2str(df))
        else
            newmessage('Redoing evaluation');
            newmessage(['Grabbing Voc#' num2str(vv) ' and Set#' num2str(df) '...']);
            if df~=olddf
                grabNewDatafile
            end
            loadnextfile(vv)
            redo=1;
        end
        
    case 'Checkbox'
        stopclick=0;
        set(checkboxh,'String','X')
        set(checkboxh,'BackgroundColor',[88 117 88]./255)
        set([submith noCallh redoh],'enable','on');
        
        
    case 'Submit'
        disableEvals
        drawnow;
        set(submith,'String','Submitted')
        set(submith,'BackgroundColor',[204 88 88]./255)
        %save submit
        ManCall=1;
        evaluationDone(vv)
        savingData
        if redo
            vv=oldvv;
            df=olddf;
            redo=0;
        end
        %load next file
        vv=vv+1;
        if vv<=Nvoc
            loadnextfile(vv)
        else
            eraseData
            if df<=length(DataFiles)
                grabNewDatafile
                loadnextfile(vv)
            else
                newmessage('Annotation done!');
                initFiles
                loadnextfile(vv)
            end
        end
        set(submith,'String','Submit')
        
    case 'EvalLog1'
        evaluatingCalls(vv,1)
    case 'EvalLog2'
        evaluatingCalls(vv,2)
    case 'EvalLog3'
        evaluatingCalls(vv,3)
    case 'EvalLog4'
        evaluatingCalls(vv,4)
    case 'EvalLog5'
        evaluatingCalls(vv,5)
    case 'EvalLog6'
        evaluatingCalls(vv,6)
    case 'EvalLog7'
        evaluatingCalls(vv,7)
    case 'EvalLog8'
        evaluatingCalls(vv,8)
    case 'EvalLog9'
        evaluatingCalls(vv,9)
    case 'EvalLog10'
        evaluatingCalls(vv,10)
        
    case 'Quit'
        close(FhGUI);        close all;
        clear all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initFiles
global DoneListDetect BaseDataDir NExpe DoneListWho ParamFile;
global Date ExpStartTime Logger_dir AudioDataPath ee WorkingDir;

newmessage('Grabbing new session...');
checkSession=1;
while checkSession && ee<=NExpe
    ee=ee+1;
    BatsID = DoneListDetect{1}{ee};
    Date = DoneListDetect{2}{ee};
    Time = DoneListDetect{3}{ee};
    ParamFile = dir(fullfile(BaseDataDir,['20' Date],'audio',sprintf('%s_%s_%s*RecOnly_param.txt', BatsID, Date, Time)));
    fprintf(1, '\n\n\n Date: %s, experiment %d/%d\n%s\n', Date,ee,NExpe,ParamFile.name)
    % Check that the file was not already set aside or done
    if ~isempty(DoneListWho)
        Done = sum(contains(DoneListWho{1},BatsID) .* contains(DoneListWho{2},Date) .* contains(DoneListWho{3},Time));
    else
        Done=0;
    end
    if Done
        fprintf(1, '   -> Data already processed\n')
        checkSession=1;
    else
        checkSession=0;
    end
end
Filepath = fullfile(ParamFile.folder, ParamFile.name);
% NCalls = result_reconly_DbatsWho(Filepath, WorkingDir);

% In case the identification of bats was already done but you want to re-do it again
ForceWhoID = 1;

% Get the recording date
[AudioDataPath, DataFile ,~]=fileparts(Filepath);
Date = DataFile(6:11);
ExpStartTime = DataFile(13:16);
Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'audiologgers');
fprintf(1,'*** Check the clock drift correction of the logger ***\n')
LoggersDir = dir(fullfile(Logger_dir, 'logger*'));
Check = zeros(length(LoggersDir)+1,1);
for ll=1:length(LoggersDir)
    FigCD = open(fullfile(LoggersDir(ll).folder, LoggersDir(ll).name,'extracted_data','CD_correction0.fig'));
    failsafe=1;
    while failsafe
        Check(ll) = 1;%input('Is everything ok? (yes ->1, No -> 0): ');
        failsafe=0;
        if isempty(Check(ll))
            failsafe=1;
            disp('Entry is empty, please repeat!')
        end
    end
    fprintf('\n')
    close(FigCD)
end
fprintf(1,'*** Check the allignement of the TTL pulses ***\n')
AllignmentPath = fullfile(AudioDataPath,sprintf('%s_%s_CD_correction_audio_piezo.fig', Date, ExpStartTime));
FigAP = open(AllignmentPath);
failsafe=1;
while failsafe
    Check(length(LoggersDir)+1) = 1;%input('Is everything ok? (yes ->1, No -> 0): ');
    failsafe=0;
    if isempty(Check(ll))
        failsafe=1;
        disp('Entry is empty, please repeat!')
    end
end
fprintf('\n')
close(FigAP)
if any(~Check)
    %NCalls = nan;
    fprintf(1,'\n****** Error in allignement reported ******\n')
else
    
    %% Identify who is calling
    fprintf('\n*** Identify who is calling ***\n')
    WhoCall_dir = dir(fullfile(Logger_dir, sprintf('*%s_%s*whocalls*', Date, ExpStartTime)));
    if isempty(WhoCall_dir) || ForceWhoID
        who_calls_playless_init('Working_dir',WorkingDir)
    else
        fprintf('\n*** ALREADY DONE: Identify who is calling ***\n')
        initFiles
    end
    %NCalls = sum(cellfun('length',IndVocStartRawMerged));
    %fprintf(FidWho, '%s\t%s\t%s\t%d\n',ParamFile.name(1:4),ParamFile.name(6:11),ParamFile.name(13:16),NCalls);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function who_calls_playless_init(varargin)
global DataFiles Working_dir Logger_dir Date ExpStartTime;
global MeanStdAmpRawExtract Voc_filename AudioDataPath DataFile Nvoc Nvocs;
global Factor_RMS_Mic Force_Save_onoffsets_mic SaveFileType pnames;
global MergeThresh df vv;
global Working_dir_read Working_dir_write PreviousFile;
global IndVocStartRaw_merged IndVocStopRaw_merged IndVocStartPiezo_merged;
global IndVocStopPiezo_merged IndVocStartRaw IndVocStopRaw IndVocStartPiezo;
global IndVocStopPiezo IndVocStart_all IndVocStop_all RMSRatio_all RMSDiff_all;
global MicError PiezoError MicErrorType PiezoErrorType Raw_wave;
global sos_raw_band_listen FS minvv maxvv UseOld;
global sos_raw_band BandPassFilter Piezo_wave AudioLogs Piezo_FS DiffRMS RMSLow VocFilename;
global Fns_AL Consecutive_binsMic Consecutive_binsPiezo Factor_RMS_low Factor_AmpDiff;
global DB_noise FHigh_spec FHigh_spec_Logger Fhigh_power Fs_env;

dflts  = {3,Logger_dir,0,'pdf'};
[Factor_RMS_Mic,Working_dir,Force_Save_onoffsets_mic,SaveFileType] = internal.stats.parseArgs(pnames,dflts,varargin{:});

Working_dir_read = fullfile(Working_dir, 'read');
Working_dir_write = fullfile(Working_dir, 'write');

if ~strcmp(Logger_dir,Working_dir) && (~exist(Working_dir,'dir') || ~exist(Working_dir_read,'dir') || ~exist(Working_dir_write,'dir'))
    mkdir(Working_dir)
    mkdir(Working_dir_read)
    mkdir(Working_dir_write)
elseif strcmp(Logger_dir,Working_dir)
    Working_dir_read = Logger_dir;
    Working_dir_write = Logger_dir;
end

DataFiles = dir(fullfile(Logger_dir, sprintf('%s_%s_VocExtractData*.mat', Date, ExpStartTime)));
if isempty(DataFiles)
    warning('Vocalization data were not extracted by get_logger_data_voc.m')
    initFiles
    loadnextfile(vv)
else
    % select the correct files
    Gdf = zeros(length(DataFiles),1);
    for df=1:length(DataFiles)
        if length(strfind(DataFiles(df).name, '_'))==2
            Gdf(df)=1;
        end
    end
    DataFiles = DataFiles(logical(Gdf));
    % gather File indices to reorder them
    IndDataFiles = nan(length(DataFiles),1);
    for nfile = 1:length(DataFiles)
        IndData = strfind(DataFiles(nfile).name, 'Data') + length('Data');
        IndDot = strfind(DataFiles(nfile).name, '.') - 1;
        IndDataFiles(nfile) = str2double(DataFiles(nfile).name(IndData:IndDot));
    end
    [~,AscendOrd] = sort(IndDataFiles);
    DataFiles = DataFiles(AscendOrd);
    load(fullfile(AudioDataPath, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)),...
        'MeanStdAmpRawExtract','Voc_filename')
    Nvoc_all = length(Voc_filename);
    DataFile = fullfile(DataFiles(1).folder, DataFiles(1).name);
    load(DataFile, 'VocMaxNum')
    if ~exist('VocMaxNum','var')
        VocMaxNum=1000;
    end
    if Nvoc_all>VocMaxNum
        Nvocs = [0 VocMaxNum:VocMaxNum:Nvoc_all (floor(Nvoc_all/VocMaxNum)*VocMaxNum+rem(Nvoc_all,VocMaxNum))];
    else
        Nvocs = [0 Nvoc_all];
    end
    check=1;
    df=6;
    while check && df<=length(DataFiles)
        
        %% Identify sound elements in each vocalization extract and decide of the vocalizer
        % Output corresponds to onset and offset indices of each vocal element in
        % sound section for each vocalizing logger. A cell array of length the
        % number of sound sections identified by voc_localize or
        % voc_localize_operant. For IndVocStartRaw_merged each cell contains the onset/offset samples of each vocal element in the sound section.
        % For IndVocStartPiezo_merged each cell is a cell array of size the number of
        % audio-loggers where each cell contains the onset or offset sample of
        % individual vocal element emitted by the respective audio-loggers. Note
        % that for both these variables, vocal elements that are within MergeThresh
        % ms were merged in a single vocal element.
        %
        % Run a time window of duration 2 ms to identify who is vocalizing
        % create 2 signals for the logger, one high pass filtered above 5kHz
        % and the other one low pass filtered at 5kHz. The vocalizer will have
        % the highest energy of all vocalizers and will have more energy in the
        % lower compare to higher filtered signal
        
        Nvoc = Nvocs(df+1) - Nvocs(df);
        DataFile = fullfile(DataFiles(df).folder, DataFiles(df).name);
        
        PreviousFile = fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat',...
            Date, ExpStartTime, df,MergeThresh));
        if ~isfile(PreviousFile)
            PreviousFile = fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData_%d.mat',...
                Date, ExpStartTime,MergeThresh));
        end
        if ~isempty(dir(PreviousFile)) && UseOld
            load(PreviousFile, 'IndVocStartRaw_merged', 'IndVocStopRaw_merged',...
                'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', ...
                'IndVocStartRaw','IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo',...
                'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all',...
                'vv','MicError','PiezoError','MicErrorType','PiezoErrorType');
            if ~exist('vv','var') % There is no previous data but just data regarding piezo numbers and bats_ID
                vv=1;
            end
            
        else
            vv=1;
        end
        
        % Check if that set of vocalizations was alread fully completed
        if vv==Nvoc
            % All done
            df=df+1;
            vv=1;
        else
            check=0;
            if ~strcmp(Working_dir_write,Logger_dir) && ~isfile(fullfile(Working_dir_read,DataFiles(df).name))
                fprintf(1,'Bringing data locally from the server\n')
                [s,m,e]=copyfile(DataFile, Working_dir_read, 'f');
                if ~s
                    fprintf(1,'File transfer did not occur correctly\n')
                    keyboard
                else
                    DataFile = fullfile(Working_dir_read,DataFiles(df).name);
                end
            end
            load(DataFile,'Raw_wave')
            if Nvoc ~= length(Raw_wave)
                warning('Looks like there might be an issue there!! Check variables!!')
                keyboard
            end
            if Nvoc<=100
                minvv = 1;
                maxvv = Nvoc;
                load(DataFile,'Piezo_wave', 'AudioLogs',   'Piezo_FS',  'FS', 'DiffRMS', 'RMSLow','VocFilename');
            else % often problem of memory, we're going to chunck file loading
                if ~mod(vv,100)
                    minvv=floor((vv-1)/100)*100 +1;
                    maxvv=ceil(vv/100)*100;
                else
                    minvv = floor(vv/100)*100 +1;
                    maxvv = ceil(vv/100)*100;
                    if maxvv<minvv
                        maxvv = ceil((vv+1)/100)*100;
                    end
                end
                Raw_wave = Raw_wave(minvv:min(maxvv, length(Raw_wave)));
                load(DataFile,'Piezo_wave', 'AudioLogs',   'Piezo_FS',  'FS', 'DiffRMS', 'RMSLow','VocFilename');
            end
        end
    end
    Fns_AL = fieldnames(Piezo_wave);
    
    % parameters
    Consecutive_binsMic = 10; % Number of consecutive bins of the envelope difference between highpass and low pass logger signal that has to be higher than threshold to be considered as a vocalization
    Consecutive_binsPiezo = 15; % Number of consecutive bins of the envelope difference between highpass and low pass logger signal that has to be higher than threshold to be considered as a vocalization
    Factor_RMS_low = 1.5.*ones(length(AudioLogs),1); % Factor by which the RMS of the low-pass filtered baseline signal is multiplied to obtained the threshold of vocalization detection on piezos
    Factor_AmpDiff = 50; % Factor by which the ratio of amplitude between low and high  pass filtered baseline signals is multiplied to obtain the threshold on calling vs hearing (when the bats call there is more energy in the lower frequency band than higher frequency band of the piezo) % used to be 3
    DB_noise = 60; % Noise threshold for the spectrogram colormap
    FHigh_spec = 90000; % Max frequency (Hz) for the raw data spectrogram
    FHigh_spec_Logger = 10000; % Max frequency (Hz) for the raw data spectrogram
    BandPassFilter = [1000 5000 9900]; % Frequency bands chosen for digital signal processing
    Fhigh_power =50; % Frequency upper bound for calculating the envelope (time running RMS)
    Fs_env = 1000; % Sample frequency of the enveloppe
    
    % if df==1 || ~exist('sos_raw_band', 'var')
    % design filters of raw ambient recording, bandpass and low pass which was
    % used for the cross correlation
    [z,p,k] = butter(6,[BandPassFilter(1) 90000]/(FS/2),'bandpass');
    sos_raw_band = zp2sos(z,p,k);
    % [z,p,k] = butter(6,BandPassFilter(1:2)/(FS/2),'bandpass');
    % sos_raw_low = zp2sos(z,p,k);
    [z,p,k] = butter(6,[100 20000]/(FS/2),'bandpass');
    sos_raw_band_listen = zp2sos(z,p,k);
    %   end
    
    
    % if df==1 || ~exist('sos_raw_band_listen', 'var')
    % design filters of raw ambient recording, bandpass, for
    % listening
    [z,p,k] = butter(6,[100 20000]/(FS/2),'bandpass');
    sos_raw_band_listen = zp2sos(z,p,k);
    % end
    
    % Initialize variables
    
    if vv==1 % We need to initialize variables!
        IndVocStart_all = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers+microphone and for each logger the index onset of when the animal start vocalizing in the piezo recording before merge in envelope unit (FS_env)
        IndVocStop_all = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal start vocalizing in the piezo recording before merge in envelope unit (FS_env)
        IndVocStartRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc)
        % a cell array of the size the number of loggers + 1 in case only one bat without a logger
        % or +2 incase no identification possible but you want to keep onset/offset of each voc and
        % for each logger the index onset of when the animal start vocalizing in the raw recording before merge
        IndVocStartPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording before merge
        IndVocStopRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording before merge
        IndVocStopPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording before merge
        IndVocStartRaw_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the raw recording
        IndVocStopRaw_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording
        IndVocStartPiezo_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording
        IndVocStopPiezo_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording
        RMSRatio_all = cell(1,Nvoc);
        RMSDiff_all = cell(1,Nvoc);
        MicError = [0 0];% first element = number of corrections. second = number of detection (question)
        MicErrorType = [0 0];% first element false negative (detected as noise or already detected when it is a new call), second element false positive (vice versa)
        PiezoError = [0 0];% first element = number of corrections. second = number of detection (question)
        PiezoErrorType = [0 0]; % first element false negative (detected as noise or hearing when it is a new call), second element false positive (vice versa)
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savingData
global PreviousFile Working_dir_write Date ExpStartTime df MergeThresh vv;
global IndVocStartRaw_merged IndVocStopRaw_merged IndVocStartPiezo_merged;
global IndVocStopPiezo_merged IndVocStartRaw IndVocStopRaw IndVocStartPiezo;
global IndVocStopPiezo IndVocStart_all IndVocStop_all RMSRatio_all RMSDiff_all;
global MicError PiezoError MicErrorType PiezoErrorType SaveRawWave Raw_wave;

while 0
    newmessage('Saving data...')
    if ~isempty(dir(PreviousFile))
        save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)),...
            'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', ...
            'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo',...
            'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all',...
            'vv','MicError','PiezoError','MicErrorType','PiezoErrorType','-append');
    else
        save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)),...
            'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged',...
            'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo',...
            'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all',...
            'vv','MicError','PiezoError','MicErrorType','PiezoErrorType');
    end
    while 0
        if SaveRawWave
            warning('Make sure you want to change that variable!! NOT recommended here as we are chuncking the loading!!!')
            keyboard
            save(DataFile, 'Raw_wave','-append')
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grabNewDatafile
global df Nvocs Nvoc DataFile DataFiles PreviousFile Working_dir_write;
global Date ExpStartTime MergeThresh UseOld vv Logger_dir Working_dir_read;
global Raw_wave  minvv  maxvv;
global Piezo_wave AudioLogs Piezo_FS FS DiffRMS RMSLow VocFilename;
global IndVocStartRaw_merged IndVocStopRaw_merged IndVocStartPiezo_merged;
global IndVocStopPiezo_merged IndVocStartRaw IndVocStopRaw IndVocStartPiezo;
global IndVocStopPiezo IndVocStart_all IndVocStop_all RMSRatio_all RMSDiff_all;
global MicError PiezoError MicErrorType PiezoErrorType;

newmessage('Grabbing new set...');
df=df+1;
Nvoc = (df+1) - Nvocs(df);%not sure why this line?
DataFile = fullfile(DataFiles(df).folder, DataFiles(df).name);

PreviousFile = fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat',...
    Date, ExpStartTime, df,MergeThresh));
if ~isfile(PreviousFile)
    PreviousFile = fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData_%d.mat',...
        Date, ExpStartTime,MergeThresh));
end
if ~isempty(dir(PreviousFile)) && UseOld
    load(PreviousFile, 'IndVocStartRaw_merged', 'IndVocStopRaw_merged',...
        'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', ...
        'IndVocStartRaw','IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo',...
        'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all',...
        'vv','MicError','PiezoError','MicErrorType','PiezoErrorType');
    % There is no previous data but just data regarding piezo numbers and bats_ID
    if ~exist('vv','var')
        vv=1;
    end
    
else
    vv=1;
end

if ~strcmp(Working_dir_write,Logger_dir) && ~isfile(fullfile(Working_dir_read,DataFiles(df).name))
    fprintf(1,'Bringing data locally from the server\n')
    [s,m,e]=copyfile(DataFile, Working_dir_read, 'f');
    if ~s
        fprintf(1,'File transfer did not occur correctly\n')
        keyboard
    else
        DataFile = fullfile(Working_dir_read,DataFiles(df).name);
    end
end
load(DataFile,'Raw_wave')
if Nvoc ~= length(Raw_wave)
    warning('Looks like there might be an issue there!! Check variables!!')
    keyboard
end
if Nvoc<=100
    minvv = 1;
    maxvv = Nvoc;
    load(DataFile,'Piezo_wave', 'AudioLogs',   'Piezo_FS',  'FS', 'DiffRMS', 'RMSLow','VocFilename');
else % often problem of memory, we're going to chunck file loading
    if ~mod(vv,100)
        minvv=floor((vv-1)/100)*100 +1;
        maxvv=ceil(vv/100)*100;
    else
        minvv = floor(vv/100)*100 +1;
        maxvv = ceil(vv/100)*100;
        if maxvv<minvv
            maxvv = ceil((vv+1)/100)*100;
        end
    end
    Raw_wave = Raw_wave(minvv:min(maxvv, length(Raw_wave)));
    load(DataFile,'Piezo_wave', 'AudioLogs',   'Piezo_FS',  'FS', 'DiffRMS', 'RMSLow','VocFilename');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadnextfile(vv)
global Nvoc df DataFiles AudioLogs Amp_env_LowPassLogVoc filenameh;
global Date ee NExpe ParamFile redoEditVoch redoEditSeth;

newmessage('Grabbing new vocalization...');
fprintf(1, '\n\n\n Date: %s, experiment %d/%d\n%s\n', Date,ee,NExpe,ParamFile.name)
checkforerror(vv)
flindx=strfind(ParamFile.name,'_');
set(filenameh,'String',[ParamFile.name(1:flindx(3)-1) ':      Voc sequence ' num2str(vv) '/' num2str(Nvoc)...
    ' Set ' num2str(df) '/' num2str(length(DataFiles))])
set(redoEditVoch,'String',num2str(vv))
set(redoEditSeth,'String',num2str(df))
grabAmbientMic(vv)
grabLoggers(vv)
% No logger data, just isolate onset/offset of vocalizations on the microphone
if sum(cellfun('isempty',(Amp_env_LowPassLogVoc))) == length(AudioLogs)
    newmessage('CANNOT DETERMINE OWNERSHIP NO DETECTION OF ONSET/OFFSET')
    %ForceSaveOnOffSetMic(vv)
else
    DataOnLogger_prep
    prepfindCaller(vv)
    enableEvals
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eraseData
global DataFile SaveRawWave DataFiles df;

while 0
    if SaveRawWave
        [s2,m,e]=copyfile(DataFile, fullfile(DataFiles(df).folder,DataFiles(df).name), 'f');
        if ~s2
            fprintf(1,'File transfer did not occur correctly\n')
            keyboard
        end
        if s2  %erase local data
            [sdel,mdel,edel]=rmdir(DataFile, 's');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkforerror(vv)
global minvv maxvv Raw_wave DataFile Piezo_wave Raw_wave_nn;
global SaveRawWave VocFilename FS;

% Patch for previous error in the code
if vv<=maxvv
    Raw_wave_nn = Raw_wave{vv - (minvv -1)};
else
    clear Raw_wave Piezo_wave
    load(DataFile,'Raw_wave')
    if ~mod(vv,100)
        minvv=floor((vv-1)/100)*100 +1;
        maxvv=ceil(vv/100)*100;
    else
        minvv = floor(vv/100)*100 +1;
        maxvv = ceil(vv/100)*100;
        if maxvv<minvv
            maxvv = ceil((vv+1)/100)*100;
        end
    end
    Raw_wave = Raw_wave(minvv:min(maxvv, length(Raw_wave)));
    load(DataFile,'Piezo_wave')
    Raw_wave_nn = Raw_wave{vv - (minvv -1)};
end
while 0
    if isempty(Raw_wave_nn)
        warning('That should not be empty now!!')
        keyboard
        SaveRawWave = 1;
        [Raw_wave{vv}, FS] = audioread(VocFilename{vv});
    else
        SaveRawWave = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grabAmbientMic(vv)
global DataFiles Raw_wave_nn Nvoc df sos_raw_band Amp_env_Mic;
global Filt_RawVoc FS Fhigh_power Fs_env DB_noise FHigh_spec ColorCode;
global Raw_Spec sliderLefth;
global plotmich;

%% First calculate the time varying RMS of the ambient microphone
% bandpass filter the ambient mic recording
Filt_RawVoc = filtfilt(sos_raw_band,1,Raw_wave_nn);
Amp_env_Mic = running_rms(Filt_RawVoc, FS, Fhigh_power, Fs_env);
% Plot the spectrogram of the ambient microphone
ColorCode = [get(groot,'DefaultAxesColorOrder');1 1 1; 0 1 1; 1 1 0];
[Raw_Spec.to, Raw_Spec.fo, Raw_Spec.logB] = ...
    spec_only_bats_gui(Filt_RawVoc, FS, DB_noise, FHigh_spec);
maxB = max(max(Raw_Spec.logB));
minB = maxB-DB_noise;
axes(plotmich);
cla(plotmich ,'reset')
hold on;
imagesc(plotmich,Raw_Spec.to*1000,Raw_Spec.fo,Raw_Spec.logB);          % to is in seconds
hold on;
axis xy;
caxis('manual');
caxis([minB maxB]);
cmap = spec_cmap();
colormap(cmap);

v_axis = axis;
v_axis(3)=0;
v_axis(4)=FHigh_spec+5e3;
axis(v_axis);
set(gca,'xlim',[0 Raw_Spec.to(end)*1000],'ylim',[0 FHigh_spec+5e3],...
    'ytick',[0 FHigh_spec/2 FHigh_spec]);

xlabel('time (ms)'), ylabel('Frequency');
yyaxis left;
plotmich.YColor='w';

hold on
yyaxis right
plot(plotmich,(1:length(Amp_env_Mic))/Fs_env*1000, Amp_env_Mic, 'r-', 'LineWidth',2)
%ylabel(sprintf('Amp\nMic'))
title(sprintf('Voc %d/%d Set %d/%d',vv,Nvoc, df,length(DataFiles)))
xlabel(' ')
set(gca, 'XTick',[],'XTickLabel',{})
set(gca,'ytick','')
drawnow;
hold off;
set(sliderLefth,'SliderStep', [1/(length(Amp_env_Mic)-1), 10/(length(Amp_env_Mic)-1)], ...
    'Min', 1, 'Max', length(Amp_env_Mic), 'Value', 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grabLoggers(vv)
global AudioLogs Piezo_wave Piezo_FS Fns_AL BandPassFilter;
global Fhigh_power Fs_env Amp_env_LowPassLogVoc;
global Amp_env_HighPassLogVoc LowPassLogVoc;

Amp_env_LowPassLogVoc = cell(length(AudioLogs),1);
Amp_env_HighPassLogVoc = cell(length(AudioLogs),1);
LowPassLogVoc = cell(length(AudioLogs),1);
%% Loop through the loggers and check the extracts length
LengthLoggersData = nan(length(AudioLogs),1);
for ll=1:length(AudioLogs)
    LengthLoggersData(ll) = length(Piezo_wave.(Fns_AL{ll}){vv});
end
%% Loop through the loggers and calculate envelopes
for ll=1:length(AudioLogs)
    if isnan(Piezo_FS.(Fns_AL{ll})(vv)) || isempty(Piezo_wave.(Fns_AL{ll}){vv})
        newmessage( 'NO DATA for Vocalization %d from %s\n', vv, Fns_AL{ll})
    else
        % design the filters
        [z,p,k] = butter(6,BandPassFilter(1:2)/(Piezo_FS.(Fns_AL{ll})(vv)/2),'bandpass');
        sos_low = zp2sos(z,p,k);
        [z,p,k] = butter(6,BandPassFilter(2:3)/(Piezo_FS.(Fns_AL{ll})(vv)/2),'bandpass');
        sos_high = zp2sos(z,p,k);
        % filter the loggers' signals
        if sum(isnan(Piezo_wave.(Fns_AL{ll}){vv}))~=length(Piezo_wave.(Fns_AL{ll}){vv})
            % replace NaN by zeros so the filter can work and
            % check the length of extracts
            InputPiezo = Piezo_wave.(Fns_AL{ll}){vv};
            if sum(isnan(InputPiezo))
                InputPiezo(isnan(InputPiezo))=0;
                if length(InputPiezo)>min(LengthLoggersData)
                    warning('Piezo data have different durations!\n Logger %s data is truncated for analysis\n',Fns_AL{ll})
                    InputPiezo = InputPiezo(1:min(LengthLoggersData));
                end
                % low-pass filter the voltage trace
                LowPassLogVoc{ll} = (filtfilt(sos_low,1,InputPiezo));
                % high-pass filter the voltage trace
                HighPassLogVoc = (filtfilt(sos_high,1,InputPiezo));
            else
                if length(InputPiezo)>min(LengthLoggersData)
                    warning('Piezo data have different durations!\n Logger %s data is truncated for analysis\n',Fns_AL{ll})
                    InputPiezo = InputPiezo(1:min(LengthLoggersData));
                end
                % low-pass filter the voltage trace
                LowPassLogVoc{ll} = (filtfilt(sos_low,1,InputPiezo));
                % high-pass filter the voltage trace
                HighPassLogVoc = (filtfilt(sos_high,1,InputPiezo));
            end
            Amp_env_LowPassLogVoc{ll}=running_rms(LowPassLogVoc{ll}, ...
                Piezo_FS.(Fns_AL{ll})(vv), Fhigh_power, Fs_env);
            Amp_env_HighPassLogVoc{ll}=running_rms(HighPassLogVoc, ...
                Piezo_FS.(Fns_AL{ll})(vv), Fhigh_power, Fs_env);
            
            % Plot the low pass filtered signal of each logger
            if ll<=10
                plotLogger(vv,ll,1)
            else
                newmessage('Logger # higher than 10, not plotted!')
            end
        else
            Amp_env_LowPassLogVoc{ll}=resample(nan(1,length(Piezo_wave.(Fns_AL{ll}){vv})),...
                Fs_env, round(Piezo_FS.(Fns_AL{ll})(vv)));
            Amp_env_HighPassLogVoc{ll}=resample(nan(1,length(Piezo_wave.(Fns_AL{ll}){vv})),...
                Fs_env, round(Piezo_FS.(Fns_AL{ll})(vv)));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prepfindCaller(vv)
global DiffAmp AudioLogs CheckMicChannel RowSize;
global IndVocStartRaw IndVocStartPiezo IndVocStopRaw IndVocStopPiezo ;
global IndVocStart IndVocStop IndVocStartRaw_merge_local IndVocStopRaw_merge_local;
global IndVocStartPiezo_merge_local IndVocStopPiezo_merge_local;
global RMSRatio RMSDiff;
global Vocp;

%% Find out which calls are emitted by each individual,
%if checking the microphone is requested it means that a single animal did not
%have a collar and all microphone calls not detected on any piezo will be attributed to it
% This is the output of the decision criterion to determine at each time point
%if a bat is vocalizing or not. each row=a logger, each column = a time point
Vocp = nan(size(DiffAmp) + [1 0]);
if CheckMicChannel
    RowSize = length(AudioLogs) +1;
else
    RowSize = length(AudioLogs);
end
IndVocStartRaw{vv} = cell(RowSize,1);
IndVocStartPiezo{vv} = cell(RowSize,1);
IndVocStopRaw{vv} = cell(RowSize,1);
IndVocStopPiezo{vv} = cell(RowSize,1);
IndVocStart = cell(RowSize,1);
IndVocStop = cell(RowSize,1);
IndVocStartRaw_merge_local = cell(RowSize,1);
IndVocStopRaw_merge_local = cell(RowSize,1);
IndVocStartPiezo_merge_local = cell(RowSize,1);
IndVocStopPiezo_merge_local = cell(RowSize,1);
RMSRatio = cell(RowSize,1);
RMSDiff = cell(RowSize,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataOnLogger_prep
global Amp_env_LowPassLogVoc Amp_env_HighPassLogVoc;
global AudioLogs ColorCode Fns_AL Factor_AmpDiff DiffRMS BandPassFilter;
global Factor_RMS_low RMSLow;
global DiffAmp RatioAmp;
global Amp_env_LowPassLogVoc_MAT Amp_env_HighPassLogVoc_MAT;
global F2 FhGUI;

%% There is some data on the logger, extract the id of the vocalizing bat
% Treat the case where some approximations of the calculations led to
% slight diffreent sizes of running RMS
Short = min(cellfun('length',Amp_env_LowPassLogVoc));
Long = max(cellfun('length',Amp_env_LowPassLogVoc));
if (Long-Short)>1
    error('The length of vectors of running RMS are too different than expected, please check!\n')
elseif Long~=Short
    Amp_env_LowPassLogVoc = cellfun(@(X) X(1:Short), Amp_env_LowPassLogVoc, 'UniformOutput',false);
    Amp_env_HighPassLogVoc = cellfun(@(X) X(1:Short), Amp_env_HighPassLogVoc, 'UniformOutput',false);
end
Amp_env_LowPassLogVoc_MAT = cell2mat(Amp_env_LowPassLogVoc);
Amp_env_HighPassLogVoc_MAT = cell2mat(Amp_env_HighPassLogVoc);
RatioAmp = (Amp_env_LowPassLogVoc_MAT +1)./(Amp_env_HighPassLogVoc_MAT+1);
DiffAmp = Amp_env_LowPassLogVoc_MAT-Amp_env_HighPassLogVoc_MAT;
while 0
    F2=figure(2);
    % Plot the ratio of time varying RMS, the difference in time varying
    % RMS between the high and low frequency bands and the absolute time
    % varying RMS of the low frequency band
    subplot(3,1,3)
    for ll=1:length(AudioLogs)
        plot(RatioAmp(ll,:), 'LineWidth',2, 'Color', ColorCode(ll,:))
        hold on
    end
    ylabel('Frequency bands Ratio')
    legend(Fns_AL{:})
    subplot(3,1,2)
    for ll=1:length(AudioLogs)
        plot(DiffAmp(ll,:), 'LineWidth',2, 'Color', ColorCode(ll,:))
        hold on
    end
    title('Calling vs Hearing')
    ylabel('Frequency bands Diff')
    hold on
    for ll=1:length(AudioLogs)
        plot([0 size(Amp_env_LowPassLogVoc_MAT,2)], Factor_AmpDiff * DiffRMS.(Fns_AL{ll})(1)*ones(2,1),...
            'Color',ColorCode(ll,:),'LineStyle','--');
        hold on
    end
    legend({Fns_AL{:} 'calling detection threshold' 'calling detection threshold'})
    subplot(3,1,1)
    for ll=1:length(AudioLogs)
        plot(Amp_env_LowPassLogVoc_MAT(ll,:), 'LineWidth',2, 'Color', ColorCode(ll,:))
        hold on
    end
    hold on
    ylabel(sprintf('Running RMS %dHz-%dHz', BandPassFilter(1:2)))
    title('Detection of vocalizations on each logger')
    for ll=1:length(AudioLogs)
        plot([0 size(Amp_env_LowPassLogVoc_MAT,2)], Factor_RMS_low(ll) * RMSLow.(Fns_AL{ll})(1)*ones(2,1),...
            'Color',ColorCode(ll,:),'LineStyle','--');
        hold on
    end
    legend({Fns_AL{:} 'Microphone' 'voc detection threshold' 'voc detection threshold'})
    winontopch=WinOnTop( FhGUI );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableEvals
global Amp_env_LowPassLogVoc_MAT AudioLogs Factor_RMS_low RMSLow Fns_AL;
global Consecutive_binsPiezo Vocp IndVocStart evalb playb;
global noCallh redoh redoEditVoch redoEditSeth sliderLefth sliderRighth;
global playMich submith plotmicevalh plotlogevalh evalbon;

for ll=1:length(AudioLogs)
    if ll<=10
        % Time points above amplitude threshold on the low-passed logger signal
        Vocp(ll,:) = Amp_env_LowPassLogVoc_MAT(ll,:)>(Factor_RMS_low(ll) * RMSLow.(Fns_AL{ll})(1));
        %find the first indices of every sequences of length "Consecutive_bins" higher than RMS threshold
        IndVocStart{ll} = strfind(Vocp(ll,:), ones(1,Consecutive_binsPiezo));
        if isempty(IndVocStart{ll})
            fprintf('\nNo vocalization detected on %s\n',Fns_AL{ll});
            set(evalb{ll},'enable','off')
            evalbon(ll)=0;
        else
            set(evalb{ll},'enable','on')
            evalbon(ll)=1;
        end
        set(evalb{ll},'String',['Eval' Fns_AL{ll}([1:3 7:end])])
        set(playb{ll},'enable','on','String',['Play' Fns_AL{ll}([1:3 7:end])])
    end
end
set([submith noCallh redoh redoEditVoch...
    redoEditSeth sliderLefth sliderRighth playMich],'enable','on')
set(submith,'BackgroundColor',[19 159 255]./255)
axes(plotmicevalh);
cla(plotmicevalh ,'reset')
axes(plotlogevalh);
cla(plotlogevalh ,'reset')
newmessage('Ready for evaluation')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disableEvals
global submith evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h evalLog6h;
global evalLog7h evalLog8h evalLog9h checkboxh;
global noCallh redoh redoEditVoch redoEditSeth sliderLefth sliderRighth;
global playMich playLog1h playLog2h playLog3h playLog4h playLog5h playLog6h;
global playLog7h playLog8h  playLog9h playLog10h playMicEvalh playLogEvalh evalLog10h;

set([evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h submith evalLog6h...
    evalLog7h evalLog8h evalLog9h evalLog10h noCallh redoh redoEditVoch...
    redoEditSeth playMich playLog1h playLog2h...
    playLog3h playLog4h playLog5h playLog6h playLog7h playLog8h playMicEvalh...
    playLog9h playLog10h playLogEvalh sliderLefth sliderRighth checkboxh],'enable','off')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallOnLogger(vv,ll)
global DiffAmp Fns_AL IndVocStart IndVocStop FHigh_spec_Logger FHigh_spec;
global RatioAmp RMSRatio RMSDiff Amp_env_Mic sliderRighth logdone Call1Listen0_temp;
global Vocp Factor_AmpDiff DiffRMS Fs_env Call1Hear0_temp plotlogevalh plotmicevalh;

logdone=0;
if ~isempty(IndVocStart{ll})
    IndVocStart_diffind = find(diff(IndVocStart{ll})>1);
    % these two lines get rid of overlapping sequences that werer detected several times
    IndVocStart{ll} = [IndVocStart{ll}(1) IndVocStart{ll}(IndVocStart_diffind +1)];
    % This is the number of detected potential vocalization
    NV = length(IndVocStart{ll});
    IndVocStop{ll} = nan(1,NV);
    RMSRatio{ll} = nan(NV,1);
    RMSDiff{ll} = nan(NV,1);
    Call1Hear0_temp = zeros(NV,1);
    Call1Listen0_temp= zeros(NV,1);
    plotMic(vv,3)
    plotLogger(vv,ll,3)
    for ii=1:NV
        IVStop = find(Vocp(ll,IndVocStart{ll}(ii):end)==0, 1, 'first');
        if ~isempty(IVStop)
            IndVocStop{ll}(ii) = IndVocStart{ll}(ii) + IVStop -1;
        else
            IndVocStop{ll}(ii) = length(Vocp(ll,:));
        end
        
        % Calculate the Average RMS Ratio for each sound extract
        % and decide about the vocalization ownership
        RMSRatio{ll}(ii) = mean(RatioAmp(ll,IndVocStart{ll}(ii):IndVocStop{ll}(ii)));
        RMSDiff{ll}(ii) = mean(DiffAmp(ll,IndVocStart{ll}(ii):IndVocStop{ll}(ii)));
        Call1Hear0_temp(ii) = RMSDiff{ll}(ii) > Factor_AmpDiff * DiffRMS.(Fns_AL{ll})(1);
        
        % update figure(3) with the decision
        
        hold on
        yyaxis left
        line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
            [FHigh_spec_Logger FHigh_spec_Logger]+3e3,'linewidth',20,'color',[0 0 0])
        yyaxis left
        line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
            [FHigh_spec_Logger FHigh_spec_Logger]-3e3,'linewidth',20,'color',[0 0 0])
        
        if Call1Hear0_temp(ii)
            %computer guess is calling
            line(plotmicevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                [FHigh_spec FHigh_spec]-3e3,'linewidth',20,'color',[.6 0 .7])
        else
            %computer guess is hearing/noise
            line(plotmicevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                [FHigh_spec FHigh_spec]-3e3,'linewidth',20,'color',[0 .7 .8])
        end
        hold off
        set(sliderRighth,'SliderStep', [1/(length(Amp_env_Mic)-1), 10/(length(Amp_env_Mic)-1)], ...
            'Min', 1, 'Max', length(Amp_env_Mic), 'Value', 1)
        
    end
else
    newmessage('Evaluation done for this logger')
    logdone=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function evaluatingCalls(vv,ll)
global IndVocStart Fs_env IndVocStop FHigh_spec_Logger checkboxh playll;
global submith noCallh redoh stopclick logdone Call1Listen0_temp;
global Call1Hear0_temp PiezoError PiezoErrorType Fns_AL;
global playMicEvalh playLogEvalh plotmicevalh plotlogevalh;
global evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h evalLog6h;
global evalLog7h evalLog8h evalLog9h evalLog10h evalbon;

playll=ll;
CallOnLogger(vv,ll)
if logdone==0
    NV = length(IndVocStart{ll});
    Call1Hear0_man = zeros(NV,1);
    Call1Listen0_man = zeros(NV,1);
    set(playLogEvalh,'enable','on','String',['Play' Fns_AL{ll}([1:3 7:end])])
    set([playMicEvalh playLogEvalh],'enable','on')
    set([submith noCallh redoh],'enable','off');
    set(checkboxh,'String','V')
    set(checkboxh,'BackgroundColor',[0 255 0]./255)
    set(checkboxh,'enable','on')
    set([evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h submith evalLog6h...
        evalLog7h evalLog8h evalLog9h evalLog10h],'enable','off')
    hold on;
    stopclick=1;
    axes(plotlogevalh);
    axorigl = axis;
    axes(plotmicevalh);
    axorigm = axis;
    while  stopclick%|| clickxv<arrowfieldl || clickyv>= arrowfield
        axes(plotlogevalh);
        [clickxv,clickyv,b]=ginput(1);%_ax(axcl,1);
        if b==2
            ax= axis;width=ax(2)-ax(1);
            axis([clickxv-width/2 clickxv+width/2 ax(3) ax(4)]);
            zoom xon;
            zoom(2)
            ax=axis;
            axes(plotmicevalh)
            axis([ax(1:2) axorigm(3) axorigm(4)])
            axes(plotlogevalh)
            axis([ax(1:2) axorigl(3) axorigl(4)])
        elseif b==3
            axes(plotlogevalh);
            axis([axorigl]);
            axes(plotmicevalh);
            axis([axorigm]);
        else
            check=1;
            ii=0;
            while check && ii<NV
                ii=ii+1;
                if clickxv>=(IndVocStart{ll}(ii)/Fs_env)*1e3  && clickxv<=(IndVocStop{ll}(ii)/Fs_env)*1e3
                    if clickyv<FHigh_spec_Logger
                        if Call1Hear0_man(ii)==1
                            Call1Hear0_man(ii)=0;
                            line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                                [FHigh_spec_Logger FHigh_spec_Logger]-3e3,'linewidth',20,'color',[0 0 0])
                        else
                            Call1Hear0_man(ii)=1;
                            line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                                [FHigh_spec_Logger FHigh_spec_Logger]-3e3,'linewidth',20,'color',[1 0 0])
                            line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                                [FHigh_spec_Logger FHigh_spec_Logger]+3e3,'linewidth',20,'color',[0 0 0])
                        end
                    elseif clickyv>FHigh_spec_Logger+2e3 && clickyv<FHigh_spec_Logger+8e3
                        if Call1Listen0_man(ii)==1
                            Call1Listen0_man(ii)=0;
                            line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                                [FHigh_spec_Logger FHigh_spec_Logger]+3e3,'linewidth',20,'color',[0 0 0])
                        else
                            Call1Listen0_man(ii)=1;
                            line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                                [FHigh_spec_Logger FHigh_spec_Logger]+3e3,'linewidth',20,'color',[0 0 1])
                            line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                                [FHigh_spec_Logger FHigh_spec_Logger]-3e3,'linewidth',20,'color',[0 0 0])
                        end
                    else
                        newmessage('Not valid click field')
                        newmessage('No choice made')
                    end
                    check=0;
                end
            end
        end
        pause(.1)
        if clickxv<0 && strcmp(get(checkboxh,'string'),'X')
            drawnow;
            stopclick=0;
        end
    end
    for ii=1:NV
        Agree = Call1Hear0_temp(ii)== Call1Hear0_man(ii);
        if ~Agree
            Call1Hear0_temp(ii) = Call1Hear0_man(ii);
            PiezoError = PiezoError + [1 1];
            PiezoErrorType = PiezoErrorType + [Call1Hear0_temp(ii) ~Call1Hear0_temp(ii)];
        else
            PiezoError = PiezoError + [0 1];
        end
        Agree = Call1Hear0_temp(ii)== Call1Listen0_man(ii);
        if ~Agree
            Call1Listen0_temp(ii) = Call1Listen0_man(ii);
        end
    end
    updateLoggerEval(vv,ll)
    evalbon(ll)=0;
    enableEvals_again
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableEvals_again
global noCallh redoh redoEditVoch redoEditSeth sliderLefth sliderRighth;
global playMich submith plotmicevalh plotlogevalh evalb evalbon;

for ll=1:length(evalbon)
    if ll<=10
        if evalbon(ll)
            set(evalb{ll},'enable','on')
        end
    end
end
set([submith noCallh redoh redoEditVoch...
    redoEditSeth sliderLefth sliderRighth playMich],'enable','on')
set(submith,'BackgroundColor',[19 159 255]./255)
axes(plotmicevalh);
cla(plotmicevalh ,'reset')
axes(plotlogevalh);
cla(plotlogevalh ,'reset')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function evaluationDone(vv)
global AudioLogs;
global IndVocStart IndVocStop IndVocStartRaw_merge_local IndVocStopRaw_merge_local;
global IndVocStartPiezo_merge_local IndVocStopPiezo_merge_local Fns_AL;
global RMSRatio RMSDiff Working_dir_write df MergeThresh FileVoc;
global ManCall SaveFileType VocFilename plotmich plotevalh Figcopy F2;
global IndVocStartRaw_merged IndVocStopRaw_merged IndVocStartPiezo_merged;
global IndVocStopPiezo_merged IndVocStart_all IndVocStop_all RMSRatio_all RMSDiff_all;

axes(plotmich)
XLIM = get(gca, 'XLim');
axes(plotevalh)
set(gca, 'YTick', 1:length(AudioLogs))
set(gca, 'YTickLabel', [Fns_AL; 'Mic'])
set(gca, 'YLim', [0 length(AudioLogs)+1])
set(gca, 'XLim', XLIM);
ylabel('AL ID')
xlabel('Time (ms)')
[~,FileVoc]=fileparts(VocFilename{vv});


while 0
    % Only save the RMS and spectro figures if there was a vocalization
    if ManCall
        copyFigure
        if strcmp(SaveFileType,'pdf')
            fprintf(1,'saving figures...\n')
            print(Figcopy,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_spec_%d.pdf',...
                FileVoc,vv, df,MergeThresh)),'-dpdf','-fillpage')
            saveas(F2,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_RMS_%d.pdf',...
                FileVoc,vv, df,MergeThresh)),'pdf')
        elseif strcmp(SaveFileType,'fig')
            fprintf(1,'saving figures...\n')
            saveas(Figcopy,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_spec_%d.fig',...
                FileVoc,vv, df,MergeThresh)))
            saveas(F2,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_RMS_%d.fig',...
                FileVoc,vv, df,MergeThresh)))
        end
    end
end

%clf(Figcopy)
%clf(F2)
% Gather Vocalization production data:
IndVocStart_all{vv} = IndVocStart;
IndVocStop_all{vv} = IndVocStop;
IndVocStartRaw_merged{vv} = IndVocStartRaw_merge_local;
IndVocStopRaw_merged{vv} = IndVocStopRaw_merge_local;
IndVocStartPiezo_merged{vv} = IndVocStartPiezo_merge_local;
IndVocStopPiezo_merged{vv} = IndVocStopPiezo_merge_local;
RMSRatio_all{vv} = RMSRatio;
RMSDiff_all{vv} = RMSDiff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function copyFigure
global Figcopy AudioLogs plotb;

Figcopy=figure(10);
for ll=1:length(AudioLogs)+2
    axcopy=subplot(10,1,ll);
    copyobj(plotb{ll}.Children,axcopy)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateLoggerEval(vv,ll)
global Call1Hear0_temp;
global IndVocStartRaw IndVocStartPiezo IndVocStopRaw IndVocStopPiezo ;
global IndVocStart IndVocStop IndVocStartRaw_merge_local IndVocStopRaw_merge_local;
global IndVocStartPiezo_merge_local IndVocStopPiezo_merge_local Fns_AL;
global ColorCode FS Piezo_FS Logger_Spec;
global Fs_env MergeThresh;
global plotb;

%!!!!!!!!!!!!!!!!!Still add evaluation of listener in Call1Listen0_temp!!!!!!!!!!!!!!!!!!!

% Stash noise and Only keep vocalizations that were produced by the bat
% according to the high RMSRatio value
IndVocStart{ll} = IndVocStart{ll}(logical(Call1Hear0_temp));
IndVocStop{ll} = IndVocStop{ll}(logical(Call1Hear0_temp));
IndVocStartRaw{vv}{ll} = round(IndVocStart{ll}/Fs_env*FS);
IndVocStartPiezo{vv}{ll} = round(IndVocStart{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
IndVocStopRaw{vv}{ll} = round(IndVocStop{ll}/Fs_env*FS);
IndVocStopPiezo{vv}{ll} = round(IndVocStop{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
NV = sum(Call1Hear0_temp);

if NV
    % Now merge detected vocalizations that are within MergeThresh
    Merge01 = (IndVocStart{ll}(2:end)-IndVocStop{ll}(1:(end-1)) < (MergeThresh/1000*Fs_env));
    NV = NV-sum(Merge01);
    IndVocStartRaw_merge_local{ll} = nan(1,NV);
    IndVocStopRaw_merge_local{ll} = nan(1,NV);
    IndVocStartPiezo_merge_local{ll} = nan(1,NV);
    IndVocStopPiezo_merge_local{ll} = nan(1,NV);
    CutInd = find(~Merge01);
    IndVocStartRaw_merge_local{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*FS);
    IndVocStartPiezo_merge_local{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
    IndVocStopRaw_merge_local{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*FS);
    IndVocStopPiezo_merge_local{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
    for cc=1:length(CutInd)
        IndVocStopRaw_merge_local{ll}(cc) = round(IndVocStop{ll}( CutInd(cc))/Fs_env*FS);
        IndVocStopPiezo_merge_local{ll}(cc) = round(IndVocStop{ll}(CutInd(cc))/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
        IndVocStartRaw_merge_local{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*FS);
        IndVocStartPiezo_merge_local{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
    end
    
    % Now plot the onset/offset of each extract on the
    % spectrograms of the loggers
    for ii=1:NV
        axes(plotb{ll+1})
        hold on
        yyaxis right
        YLim = get(gca, 'YLim');
        plot(plotb{ll+1},[IndVocStartRaw_merge_local{ll}(ii) IndVocStopRaw_merge_local{ll}(ii)]/FS*1000,...
            [YLim(2)*4/5 YLim(2)*4/5], 'k-', 'LineWidth',2)
        hold on
        plot(plotb{end},[IndVocStartRaw_merge_local{ll}(ii) IndVocStopRaw_merge_local{ll}(ii)]/FS*1000,...
            [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
        axes(plotb{end})
        set(gca,'xlim',[0 Logger_Spec{1}.to(end)*1000])
        drawnow;
        hold off
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMic(vv,FigN)
global Raw_Spec  Nvoc df DataFiles FHigh_spec;
global DB_noise plotmich plotmicevalh;

if FigN==1
    axpl=plotmich;
else
    axpl=plotmicevalh;
end
maxB = max(max(Raw_Spec.logB));
minB = maxB-DB_noise;
axes(axpl);
cla(axpl ,'reset')
imagesc(axpl,Raw_Spec.to*1000,Raw_Spec.fo,Raw_Spec.logB);% to is in seconds
hold on;
axis xy;
caxis('manual');
caxis([minB maxB]);
cmap = spec_cmap();
colormap(cmap);
v_axis = axis;
v_axis(3)=0;
v_axis(4)=FHigh_spec+5e3;
axis(v_axis);
set(gca,'xlim',[0 Raw_Spec.to(end)*1000],'ylim',[0 FHigh_spec+5e3],...
    'ytick',[0 FHigh_spec/2 FHigh_spec])
axpl.YColor='w';
xlabel('time (ms)'), ylabel('Frequency');
if FigN==1
    title(sprintf('Ambient Microphone Voc %d/%d Set %d/%d',vv,Nvoc,df,length(DataFiles)))
end

drawnow;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function plotLogger(vv,ll,FigN)
global LowPassLogVoc Piezo_FS Fns_AL DB_noise FHigh_spec_Logger;
global Amp_env_LowPassLogVoc Amp_env_HighPassLogVoc Fs_env AudioLogs;
global plotb plotlogevalh plotevalh Logger_Spec;

if FigN==1
    axes(plotevalh);
    cla(plotevalh ,'reset')
    axpl=plotb{ll+1};
    [Logger_Spec{ll}.to, Logger_Spec{ll}.fo, Logger_Spec{ll}.logB] = ...
        spec_only_bats_gui(LowPassLogVoc{ll}, Piezo_FS.(Fns_AL{ll})(vv),...
        DB_noise, FHigh_spec_Logger);
    val=0;
else
    axpl=plotlogevalh;
    val=5e3;
end
maxB = max(max(Logger_Spec{ll}.logB));
minB = maxB-DB_noise;
axes(axpl);
cla(axpl ,'reset')
hold on;

imagesc(axpl,Logger_Spec{ll}.to*1000,Logger_Spec{ll}.fo,Logger_Spec{ll}.logB);
hold on;% to is in seconds
axis xy;
caxis('manual');
caxis([minB maxB]);
cmap = spec_cmap();
colormap(cmap);

v_axis = axis;
v_axis(3)=0;
v_axis(4)=FHigh_spec_Logger+val;
axis(v_axis);
set(gca,'xlim',[0 Logger_Spec{ll}.to(end)*1000],'ylim',[0 FHigh_spec_Logger+val],...
    'ytick',[0 FHigh_spec_Logger/2 FHigh_spec_Logger])
xlabel('time (ms)'), ylabel('Frequency');
yyaxis left;
axpl.YColor='w';
yyaxis right;
set(gca,'ytick','')
if FigN==1
    hold on;
    if ll<length(AudioLogs)
        xlabel(' '); % supress the x label output
        set(gca,'XTick',[],'XTickLabel',{});
    end
    yyaxis right
    plot(axpl,(1:length(Amp_env_LowPassLogVoc{ll}))/Fs_env*1000,...
        Amp_env_LowPassLogVoc{ll}, 'b-','LineWidth', 2);
    hold on
    plot(axpl,(1:length(Amp_env_HighPassLogVoc{ll}))/Fs_env*1000,...
        Amp_env_HighPassLogVoc{ll}, 'r-','LineWidth',2);
    %ylabel(sprintf('Amp\n%s',Fns_AL{ll}([1:3 7:end])))
    hold off
else
    axpl.XColor='w';
end

hold off;
zoom on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newmessage(mstring)
global message mh Use_AppDesigner;

message={message{2};message{3};message{4};message{5};message{6};...
    mstring};
if Use_AppDesigner
    set(mh,'value',message);
else
    set(mh,'string',message);
end