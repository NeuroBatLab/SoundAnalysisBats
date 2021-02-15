function CallCurafkt(action)
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



global Nvocs df vv redo DataFiles string_handle string_handle2;
global submith noCallh redoh starth oldvv olddf;
global redoEditVoch redoEditSeth checkboxh stopclick;
global RowSize


switch action
    
    case 'Start'
        set(starth,'enable','off')
        fprintf(1,'\n*** Identify who is calling ***\n')
        newmessage('STARTING...')
        % Initializing output, finding where we left at (finding df and vv)
        % and grabing the correct set
        who_calls_playless_init
        % Loading the next vocalization identified by
        % Who_calls_playless_init
        loadnextfile
        
    case 'MaybeCall'
        loadnextfile2(vv)
        
    case 'NoCall'
        disableEvals
        drawnow;
        newmessage('Manual input enforced: Noise (=NoCall)');
        % saving data to file
        savingData(vv)
        if redo
            redo=0;
            vv=oldvv;
            if df~=olddf
                df=olddf;
                grabNewDatafile(df);
            else
                df=olddf;
            end
        end
        
        % Switching to next sound extract
        vv=vv+1;
        Nvoc = Nvocs(df+1) - Nvocs(df);
        if vv<=Nvoc % sound extract is in the same set
            % loading the next sound extract
            loadnextfile
        else % sound extract is in the next set
            % if rawwave was updated in the original data, transfer and erase the data we imported 
            transferEraseSetData
            
            if df<=length(DataFiles)
                % switch to the next set
                df=df+1;
                % grab the new dataset
                [vv,~] = grabNewDatafile(df);
                % loading the sound extract
                loadnextfile
            else % we did all sets of that recording session
                newmessage('Annotation done!');
                transferresults
                newmessage('Transfer done!\nRestart the GUI to start the next session\n')
            end
        end
        
    case 'Redo'
        % disable all evaluation butons
        disableEvals
        drawnow;
        
        
        %grab which file to reevaluate
        vv_requested=str2double(get(redoEditVoch,string_handle2));
        df_requested=str2double(get(redoEditSeth,string_handle2));
        Nvoc = Nvocs(df+1) - Nvocs(df);
        if vv_requested>Nvoc || df_requested>length(DataFiles)
            newmessage('Incorrect Voc# or Set#');
            set([redoh redoEditVoch redoEditSeth],'enable','on')
            set(redoEditVoch,string_handle2,num2str(vv))
            set(redoEditSeth,string_handle2,num2str(df))
        elseif (vv_requested>(vv+1) && df_requested>=df) || (df_requested>df)
            newmessage('Do not evaluate into the future!')
            set([redoh redoEditVoch redoEditSeth],'enable','on')
            set(redoEditVoch,string_handle2,num2str(vv))
            set(redoEditSeth,string_handle2,num2str(df))
        else
            redo=1;
            newmessage('Redoing evaluation');
            newmessage(['Grabbing Voc#' num2str(vv_requested) ' and Set#' num2str(df_requested) '...']);
            fprintf(1,'Redoing evaluation');
            fprintf(1,['Grabbing Voc#' num2str(vv_requested) ' and Set#' num2str(df_requested) '...']);
            % keep in memory the previous sound extract reference to make sure
            % we get back to where we are after redoing the requested extract
            if vv_requested~=1
                if vv_requested~=vv
                    oldvv=vv-1;
                    olddf=df;
                elseif vv_requested==vv
                    oldvv=vv;
                    olddf=df;
                end
            else % we switch to the next set!
                oldvv=Nvocs(df) - Nvocs(df-1);
                olddf=df-1;
            end
            if df_requested~=df
                df = df_requested;
                [~]=grabNewDatafile(df_requested);
                vv=vv_requested;
            else
                vv=vv_requested;
            end
            loadnextfile
        end
        
    case 'Checkbox' %(Selection Done)
        stopclick=0;
        set(checkboxh,string_handle,'X')
        set(checkboxh,'BackgroundColor',[88 117 88]./255)
        set([submith noCallh redoh],'enable','on');
        
        
    case 'Submit'
        disableEvals
        drawnow;
        set(submith,string_handle,'Submitted')
        set(submith,'BackgroundColor',[204 88 88]./255)
        % save figures and gather manual curation results only if
        % vocalizations were detected
        evaluationDone(vv)
        % saving data to file
        savingData(vv)
        if redo
            redo=0;
            vv=oldvv;
            if df~=olddf
                df=olddf;
                grabNewDatafile(df);
            else
                df=olddf;
            end
        end
        
        %load next sound extract
        vv=vv+1;
        Nvoc = Nvocs(df+1) - Nvocs(df);
        if vv<=Nvoc % next sound extract is in the same set
            loadnextfile
        else % next sound extract is in the next set
            % if rawwave was updated in the original data, transfer and erase the data we imported 
            transferEraseSetData
            
            if df<=length(DataFiles)
                % switch to the next set
                df = df +1;
                % Grabng this new set
                [vv,~]=grabNewDatafile(df);
                % Loading the next file
                loadnextfile
            else % we did all sets of that recording session
                newmessage('Annotation done!');
                transferresults
                newmessage('Transfer done!\nRestart the GUI to access the new session\n');
            end
        end
        set(submith,string_handle,'Submit')
        
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
    case 'EvalMic'
        evaluatingCalls(vv,RowSize)
        
    case 'Quit'
        close all;
        clearvars;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function who_calls_playless_init
global WorkingDir Filepath BatID Date ExpStartTime DataFiles AudioDataPath Logger_dir;
global MeanStdAmpRawExtract Voc_filename Nvocs;
global Working_dir_read Working_dir_write;
global sos_raw_band_listen FS;
global sos_raw_band BandPassFilter;
global Consecutive_binsMic Consecutive_binsPiezo Factor_RMS_low Factor_RMS_Mic Factor_AmpDiff;
global DB_noise FHigh_spec FHigh_spec_Logger Fhigh_power Fs_env;
global df vv



% parameters of the detection
Consecutive_binsMic = 10; % Number of consecutive bins of the envelope difference between highpass and low pass logger signal that has to be higher than threshold to be considered as a vocalization
Consecutive_binsPiezo = 15; % Number of consecutive bins of the envelope difference between highpass and low pass logger signal that has to be higher than threshold to be considered as a vocalization
Factor_RMS_low = 1.5; % Factor by which the RMS of the low-pass filtered baseline signal is multiplied to obtained the threshold of vocalization detection on piezos
Factor_RMS_Mic = 3; % Factor by which the RMS of the band-pass filtered baseline signal is multiplied to obtained the threshold of vocalization detection on Microphone
Factor_AmpDiff = 50; % Factor by which the ratio of amplitude between low and high  pass filtered baseline signals is multiplied to obtain the threshold on calling vs hearing (when the bats call there is more energy in the lower frequency band than higher frequency band of the piezo) % used to be 3
DB_noise = 60; % Noise threshold for the spectrogram colormap
FHigh_spec = 90000; % Max frequency (Hz) for the raw data spectrogram
FHigh_spec_Logger = 10000; % Max frequency (Hz) for the raw data spectrogram
BandPassFilter = [1000 5000 9900]; % Frequency bands chosen for digital signal processing
Fhigh_power =50; % Frequency upper bound for calculating the envelope (time running RMS)
Fs_env = 1000; % Sample frequency of the enveloppe




%% Defining Path to Data and specific parameters of that session
Working_dir_read = fullfile(WorkingDir, 'read');
Working_dir_write = fullfile(WorkingDir, 'write');

[AudioDataPath,ParamFile,~]=fileparts(Filepath);
BatID = ParamFile(1:4);
Date = ParamFile(6:11);
ExpStartTime = ParamFile(13:16);

if ~strcmp(Logger_dir,WorkingDir) && (~exist(WorkingDir,'dir') || ~exist(Working_dir_read,'dir') || ~exist(Working_dir_write,'dir'))
    mkdir(WorkingDir)
    mkdir(Working_dir_read)
    mkdir(Working_dir_write)
elseif strcmp(Logger_dir,WorkingDir)
    Working_dir_read = Logger_dir;
    Working_dir_write = Logger_dir;
end

%% Grabbing the Data coresponding to the particular date and find the correct file at which the curation was stopped
DataFiles = dir(fullfile(Logger_dir, sprintf('%s_%s_VocExtractData*.mat', Date, ExpStartTime)));
if isempty(DataFiles)
    error('Vocalization data were not extracted by get_logger_data_voc.m\nData cannot be found\n')
else
    % select the correct files
    Gdf = zeros(length(DataFiles),1);
    for dfi=1:length(DataFiles)
        if length(strfind(DataFiles(dfi).name, '_'))==2
            Gdf(dfi)=1;
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
   
    % Load the list of sound events to manually cure, and the threshold on
    % the Microphone
    load(fullfile(AudioDataPath, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)),...
        'MeanStdAmpRawExtract','Voc_filename')
    
    % Find the expected number of vocalizations Nvoc in each file
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
    
    % Loop through sets (files) to find where we left at and load the
    % previous data if requested (UseOld=1) or start from scratch (UseOld=0)
    for df=1:length(DataFiles)
        [vv,Success] = grabNewDatafile(df);
        if Success
            break
        end
    end
    if (df==length(DataFiles)) && ~Success
        newmessage('All Sets Done!');
        transferresults
    else
        %% % design filters of raw ambient recording
        % bandpass filter for detection 
        % and plotting (now that FS was loaded by grabNewDatafile
        [z,p,k] = butter(6,[BandPassFilter(1) 90000]/(FS/2),'bandpass');
        sos_raw_band = zp2sos(z,p,k);

        % design filters of raw ambient recording, bandpass, for
        % listening to microphone recordings
        [z,p,k] = butter(6,[100 20000]/(FS/2),'bandpass');
        sos_raw_band_listen = zp2sos(z,p,k);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savingData(vv)
% Save the data to files
global PreviousFile Working_dir_write Date ExpStartTime df MergeThresh AudioDataPath;
global IndVocStartRaw_merged IndVocStopRaw_merged IndVocStartPiezo_merged;
global IndVocStopPiezo_merged IndVocStartRaw IndVocStopRaw IndVocStartPiezo;
global IndVocStopPiezo IndVocStart_all IndVocStop_all RMSRatio_all RMSDiff_all;
global IndHearStartRaw IndHearStopRaw IndHearStartPiezo;
global IndHearStopPiezo IndHearStart_all IndHearStop_all;
global MicError PiezoError MicErrorType PiezoErrorType SaveRawWave Raw_wave;
global Chunking_RawWav SaveRawWaveName VocFilename Voc_filename redo

newmessage(sprintf('Saving data Voc %d...', vv))
fprintf('Saving data Voc %d...', vv)
if ~isempty(dir(PreviousFile))
    save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)),...
        'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', ...
        'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo',...
        'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all',...
        'IndHearStartRaw', 'IndHearStopRaw', 'IndHearStartPiezo', 'IndHearStopPiezo',...
        'IndHearStart_all', 'IndHearStop_all',...
        'MicError','PiezoError','MicErrorType','PiezoErrorType','-append');
else
    save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)),...
        'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged',...
        'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo',...
        'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all',...
        'IndHearStartRaw', 'IndHearStopRaw', 'IndHearStartPiezo', 'IndHearStopPiezo',...
        'IndHearStart_all', 'IndHearStop_all',...
        'MicError','PiezoError','MicErrorType','PiezoErrorType');
end
if ~redo
    save(fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat', Date, ExpStartTime,df, MergeThresh)),...
        'vv','-append');
end
if SaveRawWave
    if Chunking_RawWav
        warning('Make sure you want to change that variable!! NOT recommended here as we are chuncking the loading!!!')
        keyboard
    end
    save(DataFile, 'Raw_wave','-append')
end
if SaveRawWaveName
    if Chunking_RawWav
        warning('Highly not recommended!! We are chuncking the loading, you might loose the entire Raw_wave data here')
        keyboard
    end
    save(DataFile, 'VocFilename','-append')
    save(fullfile(AudioDataPath, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)), 'Voc_filename','-append')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vv,Success] = grabNewDatafile(df_local)
% input
global Nvocs DataFiles PreviousFile Working_dir_write redo;
global Date ExpStartTime MergeThresh UseOld Logger_dir Working_dir_read;
% output variables
global Nvoc DataFile 
global Raw_wave  minvv  maxvv;
global Piezo_wave AudioLogs Piezo_FS FS DiffRMS RMSLow VocFilename Fns_AL;
global IndVocStartRaw_merged IndVocStopRaw_merged IndVocStartPiezo_merged;
global IndVocStopPiezo_merged IndVocStartRaw IndVocStopRaw IndVocStartPiezo;
global IndVocStopPiezo IndVocStart_all IndVocStop_all RMSRatio_all RMSDiff_all;
global IndHearStart_all IndHearStop_all IndHearStartRaw IndHearStartPiezo IndHearStopRaw IndHearStopPiezo;
global MicError PiezoError MicErrorType PiezoErrorType;
global Chunking_RawWav;

newmessage(sprintf('Grabbing new set #%d...', df_local));
fprintf('Grabbing new set #%d...', df_local);
Nvoc = Nvocs(df_local+1) - Nvocs(df_local);
DataFile = fullfile(DataFiles(df_local).folder, DataFiles(df_local).name);

PreviousFile = fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData%d_%d.mat',...
    Date, ExpStartTime, df_local,MergeThresh));
if ~isfile(PreviousFile)
    % Compatibility with old version of the code
    PreviousFile = fullfile(Working_dir_write, sprintf('%s_%s_VocExtractData_%d.mat',...
        Date, ExpStartTime,MergeThresh));
end

if ~isempty(dir(PreviousFile)) && UseOld
    load(PreviousFile, 'IndVocStartRaw_merged', 'IndVocStopRaw_merged',...
        'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', ...
        'IndVocStartRaw','IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo',...
        'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all',...
        'vv','MicError','PiezoError','MicErrorType','PiezoErrorType');
    
    if ~exist('vv','var')
        % There is no previous data but just data regarding piezo numbers and bats_ID
        vv=1;
    end
else
    vv=1;
end

if vv~=Nvoc || redo % This is the file that we need to complete
    if ~strcmp(Working_dir_read,Logger_dir) && ~isfile(fullfile(Working_dir_read,DataFiles(df_local).name))
        fprintf(1,'Bringing data locally from the server\n')
        [s,m,e]=copyfile(DataFile, Working_dir_read, 'f'); %#ok<ASGLU>
        if ~s
            fprintf(1,'File transfer did not occur correctly\n')
            keyboard
        else
            DataFile = fullfile(Working_dir_read,DataFiles(df_local).name);
        end
    end
    load(DataFile,'Raw_wave')
    if Nvoc ~= length(Raw_wave)
        warning('Looks like there might be an issue there!! Check variables!!')
        keyboard
    end
    if Nvoc<=100
        Chunking_RawWav = 0;
        minvv = 1;
        maxvv = Nvoc;
        load(DataFile,'Piezo_wave', 'AudioLogs',   'Piezo_FS',  'FS', 'DiffRMS', 'RMSLow','VocFilename');
    else % often problem of memory, we're going to chunck file loading
        Chunking_RawWav = 1;
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
    
    Fns_AL = fieldnames(Piezo_wave);
    % Initialize variables
    if vv==1 && ~redo % We need to initialize variables!
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
        
        IndHearStart_all = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers+microphone and for each logger the index onset of when the animal start HEARING in the piezo recording before merge in envelope unit (FS_env)
        IndHearStop_all = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal start HEARING in the piezo recording before merge in envelope unit (FS_env)
        IndHearStartRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc)
        % a cell array of the size the number of loggers + 1 in case only one bat without a logger
        % or +2 incase no identification possible but you want to keep onset/offset of each voc and
        % for each logger the index onset of when the animal start HEARING in the raw recording before merge
        IndHearStartPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start HEARING in the piezo recording before merge
        IndHearStopRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop HEARING the raw recording before merge
        IndHearStopPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop HEARING in the piezo recording before merge
        
        
        RMSRatio_all = cell(1,Nvoc);
        RMSDiff_all = cell(1,Nvoc);
        MicError = [0 0];% first element = number of corrections. second = number of detection (question)
        MicErrorType = [0 0];% first element false negative (detected as noise or already detected when it is a new call), second element false positive (vice versa)
        PiezoError = [0 0];% first element = number of corrections. second = number of detection (question)
        PiezoErrorType = [0 0]; % first element false negative (detected as noise or hearing when it is a new call), second element false positive (vice versa)
    end
    Success = 1;
    if redo
        newmessage('Set sucessfully loaded for doing again manual annotation!');
    else
        newmessage(sprintf('Set %d sucessfully loaded!', df_local));
    end
elseif vv==Nvoc && ~redo
    newmessage(sprintf('Set %d already done!', df_local));
    Success = 0;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadnextfile
global Nvoc df vv DataFiles filenameh;
global Filepath redoEditVoch redoEditSeth string_handle string_handle2;
global noCallh maybeCallh plotb AudioLogs playMich redoh;
% Loading sound extract vv calculating microphone
% spectrogram + playing microphone

newmessage(sprintf('Grabbing next vocalization #%d...',vv));

% Find the correct wavfile
checkforerror(vv)

% Edit the Gui Info
[~,ParamFile,~]=fileparts(Filepath);
flindx=strfind(ParamFile,'_');
set(filenameh,string_handle,[ParamFile(1:flindx(3)-1) ':      Voc sequence ' num2str(vv) '/' num2str(Nvoc)...
    ' Set ' num2str(df) '/' num2str(length(DataFiles))])
set(redoEditVoch,string_handle2,num2str(vv))
set(redoEditSeth,string_handle2,num2str(df))

% Load, Check and plot the microphone file on the left pannel
Success = grabAmbientMic(vv);
if Success
    for ll=1:(length(AudioLogs))
        axes(plotb{ll+1});
        cla(plotb{ll+1} ,'reset')
    end
    axes(plotb{end});
    cla(plotb{end} ,'reset')
    set([noCallh maybeCallh playMich redoh redoEditVoch...
        redoEditSeth],'enable','on')
else
    vv=vv+1;
    loadnextfile
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadnextfile2(vv)
global AudioLogs Amp_env_LowPassLogVoc noCallh maybeCallh;
% Getting ready for evaluation (calculating
% spectrograms of loggers, plotting them,finding calls...)
set([noCallh maybeCallh],'enable','off')
newmessage(sprintf('Finishing loading #%d...',vv));

% Load piezo data, calculate envelope of high pass and low pass fltered
% signals and plot the piezo signals on the left pannel
grabLoggers(vv)


if sum(cellfun('isempty',(Amp_env_LowPassLogVoc))) == length(AudioLogs)
    % No logger data, just isolate onset/offset of vocalizations on the microphone
    newmessage('CANNOT DETERMINE OWNERSHIP NO DETECTION OF ONSET/OFFSET')
    %ForceSaveOnOffSetMic(vv)
    %!!!! This option is not availabe as of now!!!!
else
    % Get ready variables
    % Detect vocalizations on Loggers and enable buttons for manual evaluation
    % of loggers that have sound events
    prepfindCaller(vv)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function transferEraseSetData
global DataFile SaveRawWave DataFiles df Working_dir_read;

if SaveRawWave && ~strcmp(Working_dir_read,DataFiles(df).folder)
        [s2,m,e]=copyfile(DataFile, fullfile(DataFiles(df).folder,DataFiles(df).name), 'f'); %#ok<ASGLU>
        if ~s2
            warning('File transfer of %s to %s did not occur correctly\n',DataFile,fullfile(DataFiles(df).folder,DataFiles(df).name))
            keyboard
        end
        if s2  %erase local data
%             [sdel,mdel,edel]=rmdir(DataFile, 's');
            [~]=rmdir(DataFile, 's');
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function transferresults
global Working_dir_write Logger_dir FidWho IndVocStartRaw_merged BatID Date ExpStartTime
%%% TRANSFER DATA BACK ON SERVER AND KEEP TRACK THAT THIS EXPERIMENT DAY IS DONE%%%%
if ~strcmp(Working_dir_write,Logger_dir)
    newmessage('Transferring data back on the server\n')
    fprintf(1,'Transferring data back on the server\n')
    [s1,m,e]=copyfile(fullfile(Working_dir_write,'*'), Logger_dir, 'f');
    if ~s1
        newmessage('!!!! File transfer did not occur correctly!!!!\n')
        newmessage('!!!! Contact JULIE Asap via Slack, Stop here!!!!\n')
        error('!!!! File transfer did not occur correctly!!!!\n!!!! Contact JULIE Asap via Slack, Stop here!!!!\n')
    end
    if s1  %erase local data
        [sdel,mdel,edel]=rmdir(Working_dir_write, 's');
    end
    
end
NCalls = sum(cellfun('length',IndVocStartRaw_merged));
fprintf(FidWho, '%s\t%s\t%s\t%d\n',BatID,Date,ExpStartTime,NCalls);
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkforerror(vv)
global minvv maxvv Raw_wave DataFile Piezo_wave Raw_wave_nn Chunking_RawWav;
global SaveRawWave VocFilename FS;

% Patch for previous error in the code
if vv<=maxvv
    Raw_wave_nn = Raw_wave{vv - (minvv -1)};
else
    clear Raw_wave Piezo_wave
    load(DataFile,'Raw_wave')
    Chunking_RawWav = 1;
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

if isempty(Raw_wave_nn)
    warning('That should not be empty now!!')
    keyboard
    SaveRawWave = 1;
    [Raw_wave{vv}, FS] = audioread(VocFilename{vv});
else
    SaveRawWave = 0;
end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Success = grabAmbientMic(vv)
% Load, Check and plot the microphone file for sound event vv
% global DataFile AudioDataPath Date ExpStartTime
global Raw_wave_nn Nvoc Nvocs df sos_raw_band Amp_env_Mic;
global Voc_filename VocFilename
global Filt_RawVoc FS Fhigh_power Fs_env DB_noise FHigh_spec ColorCode;
global Raw_Spec;
% global SaveRawWaveName SaveRawWave

%% Check that correct microphone file was saved (Trying to
% detect/fix bug from voc_localize_using_piezo

% retrieving file name index of the microphone
Voc_i_start = Nvocs(df)+1;
vv_in = vv + Voc_i_start-1;
if ~strcmp(Voc_filename{vv_in}, VocFilename{vv})
    warning('Issues with Mic file name\n')
    keyboard
end
fprintf(1, 'Microphone File: %s\n', Voc_filename{vv_in})
if Nvoc>100
    warning('Probably wrong audio file name, the code is not updated for older version of previous extraction\n')
end

% Section that was checking the correct extraction of microphone data
% fprintf(1, 'Checking correct mic file was selected\n')
% TimerMicCheck = tic;
% load(DataFile,'Voc_transc_time_refined');
% TTL_dir = dir(fullfile(AudioDataPath,sprintf( '%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
% TTL = load(fullfile(TTL_dir.folder, TTL_dir.name));
% FileNum_u = unique(TTL.File_number);
% OnOffTranscTime_ms = Voc_transc_time_refined(vv,:);
% FileNumIdx = find(TTL.Pulse_TimeStamp_Transc<OnOffTranscTime_ms(1,1),1,'Last');
% MicVoc_File = TTL.File_number(FileNumIdx);
% IndFileNum = find(FileNum_u == MicVoc_File);
% TranscTime_zs = (OnOffTranscTime_ms - TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,1))/TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,2);
% MicVoc_samp_idx =round(TTL.Mean_std_Pulse_samp_audio(IndFileNum,2) .* polyval(TTL.Slope_and_intercept_transc2audiosamp{IndFileNum},TranscTime_zs,[],TTL.Mean_std_x_transc2audiosamp{IndFileNum}) + TTL.Mean_std_Pulse_samp_audio(IndFileNum,1));
% WavFileStruc_local = dir(fullfile(AudioDataPath, sprintf('*_%s_%s*mic*_%d.wav',Date, ExpStartTime, MicVoc_File)));
% Raw_filename = fullfile(WavFileStruc_local.folder, WavFileStruc_local.name);
% [Raw_10minwav2, FS2] = audioread(Raw_filename);
% if MicVoc_samp_idx(1)>length(Raw_10minwav2) % This vocalization occured in the next file
%     MicVoc_File = MicVoc_File+1;
%     IndFileNum = find(FileNum_u == MicVoc_File);
%     TranscTime_zs = (OnOffTranscTime_ms - TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,1))/TTL.Mean_std_Pulse_TimeStamp_Transc(IndFileNum,2);
%     MicVoc_samp_idx =round(TTL.Mean_std_Pulse_samp_audio(IndFileNum,2) .* polyval(TTL.Slope_and_intercept_transc2audiosamp{IndFileNum},TranscTime_zs,[],TTL.Mean_std_x_transc2audiosamp{IndFileNum}) + TTL.Mean_std_Pulse_samp_audio(IndFileNum,1));
%     WavFileStruc_local = dir(fullfile(AudioDataPath, sprintf('*_%s_%s*mic*_%d.wav',Date, ExpStartTime, MicVoc_File)));
%     Raw_filename = fullfile(WavFileStruc_local.folder, WavFileStruc_local.name);
%     [Raw_10minwav2, FS2] = audioread(Raw_filename);
% end
% Raw_wave_ex = Raw_10minwav2(MicVoc_samp_idx(1) : min(MicVoc_samp_idx(2),length(Raw_10minwav2)));
% if length(Raw_wave_ex)<length(Raw_wave_nn)
%     Corr(1) = corr(Raw_wave_ex,Raw_wave_nn(1:length(Raw_wave_ex)));
%     Corr(2) = corr(Raw_wave_ex,Raw_wave_nn(end-length(Raw_wave_ex)+1:end));
% elseif length(Raw_wave_ex)>length(Raw_wave_nn)
%     Corr(1) = corr(Raw_wave_nn,Raw_wave_ex(1:length(Raw_wave_nn)));
%     Corr(2) = corr(Raw_wave_nn,Raw_wave_ex(end-length(Raw_wave_nn)+1:end));
% elseif length(Raw_wave_ex)==length(Raw_wave_nn)
%     Corr = corr(Raw_wave_ex,Raw_wave_nn);
% end
% if all(Corr<0.99)
%     warning('Error in the microphone file that was previosuly saved, fixing the issue now!\n')
%     %                 keyboard
%     SaveRawWave = 1;
%     Raw_wave_nn = Raw_wave_ex;
%     Raw_wave{vv} = Raw_wave_ex;
%     TrueVocName = fullfile(AudioDataPath, 'Detected_calls',sprintf('%s_%s_%s_voc_%d_%d.wav',WavFileStruc_local.name(1:4),Date,ExpStartTime, MicVoc_File, MicVoc_samp_idx(1)));
%     if ~strcmp(VocFilename{vv}, TrueVocName)
%         warning('Filename was also wrong correcting %s -> %s\n',VocFilename{vv},TrueVocName)
%         VocFilename{vv}= TrueVocName;
%         Voc_filename{vv_in} = TrueVocName;
%         SaveRawWaveName = 1;
%     end
%     audiowrite(VocFilename{vv} , Raw_wave{vv}, FS2);
% end
% fprintf(1, '-> DONE in %.1f s\n', toc(TimerMicCheck))
% clear TimerMicCheck Raw_10minwav2 Raw_wave_ex

%% First calculate the time varying RMS of the ambient microphone
% bandpass filter the ambient mic recording
if length(Raw_wave_nn)/FS>=0.1
    Filt_RawVoc = filtfilt(sos_raw_band,1,Raw_wave_nn);
    Amp_env_Mic = running_rms(Filt_RawVoc, FS, Fhigh_power, Fs_env);
    
    %% Plot the spectrogram of the ambient microphone on the left pannel
    ColorCode = [get(groot,'DefaultAxesColorOrder');1 1 1; 0 1 1; 1 1 0];
    [Raw_Spec.to, Raw_Spec.fo, Raw_Spec.logB] = ...
        spec_only_bats_gui(Filt_RawVoc, FS, DB_noise, FHigh_spec);
    plotMic(vv,1)
    PlayCallfkt('PlayMic')
    Success=1;
else
    warning('That sound is too short to be processed correctly.\n')
    Success=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grabLoggers(vv)
% Load piezo data, calculate envelope of high pass and low pass fltered
% signals and plot the piezo signals
global AudioLogs Piezo_wave Piezo_FS Fns_AL BandPassFilter;
global Fhigh_power Fs_env Amp_env_LowPassLogVoc;
global Amp_env_HighPassLogVoc LowPassLogVoc Logger_Spec;

Amp_env_LowPassLogVoc = cell(length(AudioLogs),1);
Amp_env_HighPassLogVoc = cell(length(AudioLogs),1);
LowPassLogVoc = cell(length(AudioLogs),1);

%% Loop through the loggers and check the extracts length
LengthLoggersData = nan(length(AudioLogs),1);
for ll=1:length(AudioLogs)
    LengthLoggersData(ll) = length(Piezo_wave.(Fns_AL{ll}){vv});
end
%% Loop through the loggers and calculate envelopes
Logger_Spec = cell(length(AudioLogs),1);
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
                newmessage(sprintf('Logger %d higher than 10, not plotted!',ll))
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
global Amp_env_LowPassLogVoc Amp_env_HighPassLogVoc Amp_env_Mic;
global AudioLogs ColorCode Fns_AL Factor_AmpDiff DiffRMS BandPassFilter;
global Factor_RMS_low RMSLow;
global DiffAmp RatioAmp;
global Amp_env_LowPassLogVoc_MAT Amp_env_HighPassLogVoc_MAT;
global RMSFig FhGUI;

global CheckMicChannel RowSize Nvocs df Factor_RMS_Mic MeanStdAmpRawExtract;
global IndVocStartRaw IndVocStartPiezo IndVocStopRaw IndVocStopPiezo ;
global IndVocStart IndVocStop IndVocStartRaw_merge_local IndVocStopRaw_merge_local;
global IndVocStartPiezo_merge_local IndVocStopPiezo_merge_local;
global IndHearStartRaw IndHearStartPiezo IndHearStopRaw IndHearStopPiezo ;
global IndHearStart IndHearStop
global RMSRatio RMSDiff PlotRMSFig;
global Vocp;

global Consecutive_binsPiezo Consecutive_binsMic  evalLog evalMich playb string_handle;
global noCallh redoh redoEditVoch redoEditSeth sliderLefth sliderRighth;
global playMich submith plotmicevalh plotlogevalh evalbon evalmicbon ;

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

if CheckMicChannel
    try
        Short = min([cellfun('length',Amp_env_LowPassLogVoc)' length(Amp_env_Mic)]);
    catch
        Short = min([cellfun('length',Amp_env_LowPassLogVoc) length(Amp_env_Mic)]);
    end
    try
        Long = max([cellfun('length',Amp_env_LowPassLogVoc)' length(Amp_env_Mic)]);
    catch
        Long = max([cellfun('length',Amp_env_LowPassLogVoc) length(Amp_env_Mic)]);
    end
    if (Long-Short)>1
        error('The length of vectors of running RMS are too different than expected, please check!\n')
    elseif Long~=Short
        Amp_env_LowPassLogVoc = cellfun(@(X) X(1:Short), Amp_env_LowPassLogVoc, 'UniformOutput',false);
        Amp_env_HighPassLogVoc = cellfun(@(X) X(1:Short), Amp_env_HighPassLogVoc, 'UniformOutput',false);
        Amp_env_Mic = Amp_env_Mic(1:Short);
    end
    Amp_env_LowPassLogVoc_MAT = cell2mat(Amp_env_LowPassLogVoc);
    Amp_env_HighPassLogVoc_MAT = cell2mat(Amp_env_HighPassLogVoc);
    RatioAmp = (Amp_env_LowPassLogVoc_MAT +1)./(Amp_env_HighPassLogVoc_MAT+1);
    DiffAmp = Amp_env_LowPassLogVoc_MAT-Amp_env_HighPassLogVoc_MAT;
end

%% plot the RMS values if requested
if PlotRMSFig
    RMSFig=figure(2);
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
    legend([Fns_AL(:)' {'calling detection threshold' 'calling detection threshold'}])
    subplot(3,1,1)
    for ll=1:length(AudioLogs)
        plot(Amp_env_LowPassLogVoc_MAT(ll,:), 'LineWidth',2, 'Color', ColorCode(ll,:))
        hold on
    end
    hold on
    ylabel(sprintf('Running RMS %dHz-%dHz', BandPassFilter(1:2)))
    title('Detection of vocalizations on each logger')
    for ll=1:length(AudioLogs)
        plot([0 size(Amp_env_LowPassLogVoc_MAT,2)], Factor_RMS_low * RMSLow.(Fns_AL{ll})(1)*ones(2,1),...
            'Color',ColorCode(ll,:),'LineStyle','--');
        hold on
    end
    legend([Fns_AL(:)' {'Microphone' 'voc detection threshold' 'voc detection threshold'}])
    WinOnTop( FhGUI );
end

%% Initialize variables to find out which calls are emitted by each individual,
%if checking the microphone is requested it means that a single animal did not
%have a collar and all microphone calls not detected on any piezo will be attributed to it
% This is the output of the decision criterion to determine at each time point
%if a bat is vocalizing or not. each row=a logger, each column = a time point
if CheckMicChannel
    RowSize = length(AudioLogs) +1;
    Vocp = nan(length(AudioLogs) +1,size(DiffAmp,2));
else
    RowSize = length(AudioLogs);
    Vocp = nan(length(AudioLogs),size(DiffAmp,2));
end
evalbon = nan(length(AudioLogs),1); % This vector keeps track of which
% logger would still need to be evaluated (potential calls)
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
IndHearStartRaw{vv} = cell(RowSize,1);
IndHearStartPiezo{vv} = cell(RowSize,1);
IndHearStopRaw{vv} = cell(RowSize,1);
IndHearStopPiezo{vv} = cell(RowSize,1);
IndHearStart = cell(RowSize,1);
IndHearStop = cell(RowSize,1);

RMSRatio = cell(RowSize,1);
RMSDiff = cell(RowSize,1);


%% Detect vocalizations on Loggers and enable buttons for manual evaluation
for ll=1:length(AudioLogs)
    if ll<=10
        % Time points above amplitude threshold on the low-passed logger signal
        Vocp(ll,:) = Amp_env_LowPassLogVoc_MAT(ll,:)>(Factor_RMS_low * RMSLow.(Fns_AL{ll})(1));
        %find the first indices of every sequences of length "Consecutive_bins" higher than RMS threshold
        IndVocStart{ll} = strfind(Vocp(ll,:), ones(1,Consecutive_binsPiezo));
        if isempty(IndVocStart{ll})
            fprintf(1, '\nNo vocalization detected on %s\n',Fns_AL{ll});
            set(evalLog{ll},'enable','off')
            evalbon(ll)=0;
        else
            set(evalLog{ll},'enable','on')
            evalbon(ll)=1;
        end
        set(evalLog{ll},string_handle,['Eval' Fns_AL{ll}([1:3 7:end])])
        set(playb{ll},'enable','on',string_handle,['Play' Fns_AL{ll}([1:3 7:end])])
    end
end

if CheckMicChannel
    VocpMic = Amp_env_Mic>(Factor_RMS_Mic * MeanStdAmpRawExtract(Nvocs(df)+vv,1)); % Time points above amplitude threshold on the band-pass microphone signal
    Vocp(length(AudioLogs)+1,:) = reshape(VocpMic,1,length(VocpMic));
    IndVocStart{length(AudioLogs)+1} = strfind(Vocp(length(AudioLogs)+1,:), ones(1,Consecutive_binsMic)); %find the first indices of every sequences of length "Consecutive_bins" higher than RMS threshold
    if isempty(IndVocStart{length(AudioLogs)+1})
        fprintf(1,'\nNo vocalization detected on microphone\n');
        set(evalMich,'enable','off')
        evalmicbon=0;
    else% Some vocalizations were detected
        set(evalMich,'enable','on')
        evalmicbon=1;
    end
    set(evalMich,string_handle,'Eval Mic')
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
global submith evalLog evalMich;
global  checkboxh;
global noCallh maybeCallh redoh redoEditVoch redoEditSeth sliderLefth sliderRighth;
global playMich playLog1h playLog2h playLog3h playLog4h playLog5h playLog6h;
global playLog7h playLog8h  playLog9h playLog10h playMicEvalh playLogEvalh;

set([evalLog{1} evalLog{2} evalLog{3} evalLog{4} evalLog{5} submith evalLog{6}...
    evalLog{7} evalLog{8} evalLog{9} evalLog{10} evalMich noCallh maybeCallh redoh redoEditVoch...
    redoEditSeth playMich playLog1h playLog2h...
    playLog3h playLog4h playLog5h playLog6h playLog7h playLog8h playMicEvalh...
    playLog9h playLog10h playLogEvalh sliderLefth sliderRighth checkboxh],'enable','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallOnLogger(vv,ll)
% finding the offset of each sound extract in that logger or on Mic (ll>length(AudioLogs)); plot the
% spectrograms of the microphone and logger ll (ll<=length(AudioLogs)) in the evaluation area;
% plot the sound extracts in the evaluation area 

global DiffAmp Fns_AL IndVocStart IndVocStop FHigh_spec_Logger  FHigh_spec Hline AudioLogs;
global RatioAmp RMSRatio RMSDiff Amp_env_Mic sliderRighth logdone;
global Vocp Factor_AmpDiff DiffRMS Fs_env NoiseCall01_ComputerPredict plotlogevalh ;
% global FHigh_spec plotmicevalh ;

logdone=0;
if ~isempty(IndVocStart{ll}) % Some vocalizations were detected for that logger or Mic track
    IndVocStart_diffind = find(diff(IndVocStart{ll})>1);
    % these two lines get rid of overlapping sequences that werer detected several times
    IndVocStart{ll} = [IndVocStart{ll}(1) IndVocStart{ll}(IndVocStart_diffind +1)];
    % This is the number of detected potential vocalization
    NV = length(IndVocStart{ll});
    IndVocStop{ll} = nan(1,NV);
    NoiseCall01_ComputerPredict = zeros(NV,1);
    if ll<=length(AudioLogs) % this is a logger
        RMSRatio{ll} = nan(NV,1);
        RMSDiff{ll} = nan(NV,1);
    end
    
    % Plot the spectrogram of the microphone in the upper evaluation area
    plotMic(vv,3)
    % Plot the spectrogram of the logger or Microphone in the lower evaluation area
    if ll<=length(AudioLogs)
        plotLogger(vv,ll,3)
    else
        plotMic(vv,2)
    end
    Hline = cell(NV,1); % This contains the handle to the lines to change...
    % their colors in evaluatingCalls
    for ii=1:NV
        IVStop = find(Vocp(ll,IndVocStart{ll}(ii):end)==0, 1, 'first');
        if ~isempty(IVStop)
            IndVocStop{ll}(ii) = IndVocStart{ll}(ii) + IVStop -1;
        else
            IndVocStop{ll}(ii) = length(Vocp(ll,:));
        end
        
        % Calculate the Average RMS Ratio for each sound extract
        % and decide about the vocalization ownership
        if ll<=length(AudioLogs)
            RMSRatio{ll}(ii) = mean(RatioAmp(ll,IndVocStart{ll}(ii):IndVocStop{ll}(ii)));
            RMSDiff{ll}(ii) = mean(DiffAmp(ll,IndVocStart{ll}(ii):IndVocStop{ll}(ii)));
            NoiseCall01_ComputerPredict(ii) = RMSDiff{ll}(ii) > Factor_AmpDiff * DiffRMS.(Fns_AL{ll})(1);
        end
        
        % update the logger evaluation pannel (right) with the position of
        % the sound extracts in black
        hold on
        yyaxis left
        if ll<=length(AudioLogs)
            % Computer prediction of calling behavior
            if NoiseCall01_ComputerPredict(ii)
                %computer guess is calling
                Hline{ii} = line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                [FHigh_spec_Logger FHigh_spec_Logger]-2e3,'linewidth',20,'color',[1 0 0]);
                line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                [FHigh_spec_Logger FHigh_spec_Logger]+1e3,'linewidth',20,'color',[.6 0 .7]);
            else
                %computer guess is not calling
                Hline{ii} = line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                [FHigh_spec_Logger FHigh_spec_Logger]-2e3,'linewidth',20,'color',[0 0 0]);
                line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
                [FHigh_spec_Logger FHigh_spec_Logger]+1e3,'linewidth',20,'color',[0 0 0]);
            end
        else
            Hline{ii} = line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
            [FHigh_spec FHigh_spec]-15e3,'linewidth',20,'color',[0 0 0]);
            line(plotlogevalh,[IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000,...
            [FHigh_spec FHigh_spec]-5e3,'linewidth',20,'color',[0 0 0]);
        end
        hold off
        set(sliderRighth,'SliderStep', [1/(length(Amp_env_Mic)-1), 10/(length(Amp_env_Mic)-1)], ...
            'Min', 1, 'Max', length(Amp_env_Mic), 'Value', 1)
        
    end
elseif isempty(IndVocStart{ll}) && ll<=length(AudioLogs)
    newmessage('Evaluation done for this logger: no vocalizations')
    logdone=1;
elseif isempty(IndVocStart{ll}) && ll==length(AudioLogs)+1
    newmessage('Evaluation done for the microphone: no vocalizations')
    logdone=1;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function evaluatingCalls(vv,ll)
global IndVocStart Fs_env IndVocStop Hline FHigh_spec checkboxh playll;
global submith noCallh redoh stopclick logdone AudioLogs;
global  Fns_AL NoiseCall01_ComputerPredict;
global playMicEvalh playLogEvalh plotmicevalh plotlogevalh;
global evalLog evalMich evalbon evalmicbon string_handle;
% global FHigh_spec_Logger NoiseCall01_ComputerPredict PiezoError PiezoErrorType

playll=ll;

% finding the offset of each sound extract in that logger; plot the
% spectrograms of the microphone and logger ll in the evaluation area;
% plot the sound extracts in the evaluation area 
CallOnLogger(vv,ll)

% Now enabling Play buttons of microphone and logger in the evaluation pannel
% disable evaluation button on the left pannel and gather clicking
% data
if logdone==0
    NV = length(IndVocStart{ll});
    NoiseCallListen012_man = NoiseCall01_ComputerPredict; % Initialize all sound elements as noise
    if ll<=length(AudioLogs)
        set(playLogEvalh,'enable','on',string_handle,['Play' Fns_AL{ll}([1:3 7:end])])
        set(playLogEvalh,'enable','on')
    end
    set(playMicEvalh,'enable','on')
    set([submith noCallh redoh],'enable','off');
    set(checkboxh,string_handle,'V')
    set(checkboxh,'BackgroundColor',[0 255 0]./255)
    set(checkboxh,'enable','on')
    set([evalLog{1} evalLog{2} evalLog{3} evalLog{4} evalLog{5} submith evalLog{6}...
        evalLog{7} evalLog{8} evalLog{9} evalLog{10} evalMich],'enable','off')
    hold on;
    stopclick=1;
    axes(plotlogevalh);
    axorigl = axis;
    axes(plotmicevalh);
    axorigm = axis;
    while  stopclick%|| clickxv<arrowfieldl || clickyv>= arrowfield
        axes(plotlogevalh); %#ok<*LAXES>
        [clickxv,clickyv,b]=ginput(1);%_ax(axcl,1);
        if b==2 % Click on middle button of the mouse or Shift + click zoom in on the xaxis
            ax= axis;
            width=ax(2)-ax(1);
            axis([clickxv-width/2 clickxv+width/2 ax(3) ax(4)]);
            zoom xon;
            zoom(2)
            ax=axis;
            axes(plotmicevalh)
            axis([ax(1:2) axorigm(3) axorigm(4)])
            axes(plotlogevalh)
            axis([ax(1:2) axorigl(3) axorigl(4)])
        elseif b==3 % Right click recover the original axis
            axes(plotlogevalh);
            axis(axorigl);
            axes(plotmicevalh);
            axis(axorigm);
        else % left click do a selection
            % find the indices of the sound element that was selected
            ii = find((clickxv>=(IndVocStart{ll}/Fs_env)*1e3) .* (clickxv<=(IndVocStop{ll}/Fs_env)*1e3));
            if isempty(ii) %This was an incorrect click
                continue
            elseif clickyv>0 && clickyv<FHigh_spec + 3e3
                NoiseCallListen012_man(ii) = NoiseCallListen012_man(ii)+1;
                if NoiseCallListen012_man(ii)>2 % Make sure we are iterating through the 3 possibilities (noise ->0, vocalize->1, hear->2)
                    NoiseCallListen012_man(ii)=0;
                end
                if NoiseCallListen012_man(ii)==0 % Noise
                    Hline{ii}.Color = [0 0 0];
                elseif NoiseCallListen012_man(ii)==1 % vocalize
                    Hline{ii}.Color = [1 0 0];
                elseif NoiseCallListen012_man(ii)==2 %hear
                    Hline{ii}.Color = [0 0 1];
                end
            end
        end
        pause(.1)
        if clickxv<0 && strcmp(get(checkboxh,string_handle),'X')
            drawnow;
            stopclick=0; % this sounds useless to me as the checkboxh callback function already set stopclick to 0  
        end
    end
    
    
    
%     for ii=1:NV
%         Agree = (NoiseCall01_ComputerPredict(ii)) && (NoiseCallListen012_man(ii)==1);
%         if ~Agree
%             NoiseCall01_ComputerPredict(ii) = NoiseCallListen012_man(ii);
%             PiezoError = PiezoError + [1 1];
%             PiezoErrorType = PiezoErrorType + [NoiseCallListen012_man(ii) ~NoiseCallListen012_man(ii)];
%         else
%             PiezoError = PiezoError + [0 1];
%         end
%     end

    % Save production and hearing data for that logger, merge produced
    % vocalizations, update bottom plot on left panel with produced
    % vocalizations.
    updateLoggerEval(vv,ll,NoiseCallListen012_man)
    
    % Indicate evaluation as done
    if ll<= length(AudioLogs)
        evalbon(ll)=0;
    else
        evalmicbon = 0;
    end
        
    
    % Enable evaluation buttons of each logger that still need to be evaluated
    % according to variable evalbon and reset evaluation panel plots
    enableEvals_again
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableEvals_again

% Enable evaluation buttons of each logger that still need to be evaluated
% according to variable evalbon and reset evaluation panel plots
global noCallh redoh redoEditVoch redoEditSeth sliderLefth sliderRighth;
global playMich submith plotmicevalh plotlogevalh evalLog evalMich evalbon evalmicbon;

for ll=1:length(evalbon)
    if ll<=10
        if evalbon(ll)
            set(evalLog{ll},'enable','on')
        end
    end
end
if evalmicbon
    set(evalMich, 'enable', 'on')
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
global SaveFileType VocFilename plotmich plotevalh plotb RMSFig PlotRMSFig;
global IndVocStartRaw_merged IndVocStopRaw_merged IndVocStartPiezo_merged;
global IndVocStopPiezo_merged IndVocStart_all IndVocStop_all RMSRatio_all RMSDiff_all;
global IndHearStart IndHearStop IndHearStart_all IndHearStop_all

% Update the figure on the left panel, in particular scale the bottom
% summary plot
axes(plotmich)
XLIM = get(gca, 'XLim');
axes(plotevalh)
set(gca, 'YTick', 1:(length(AudioLogs)+1))
set(gca, 'YTickLabel', [Fns_AL; 'Mic'])
set(gca, 'YLim', [0 length(AudioLogs)+1])
set(gca, 'XLim', XLIM);
ylabel('AL ID')
xlabel('Time (ms)')
set(gca, 'YColor', 'w');
set(gca, 'XColor', 'w');
[~,FileVoc]=fileparts(VocFilename{vv});

% check if evaluation produced some call detection, if not no need to save
% figures and data
if sum(cellfun(@isempty,IndVocStartRaw_merge_local)) == length(IndVocStartRaw_merge_local)
    % No vocalization heard or produced
else
    % Save the RMS and spectro figures
    Figcopy = copyFigure(plotb,vv); % copy left pannel in a figure
    if strcmp(SaveFileType,'pdf')
        fprintf(1,'saving figures...\n')
        print(Figcopy,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_spec_%d.pdf',...
            FileVoc,vv, df,MergeThresh)),'-dpdf','-fillpage')
        if PlotRMSFig
            saveas(RMSFig,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_RMS_%d.pdf',...
            FileVoc,vv, df,MergeThresh)),'pdf')
            clf(RMSFig)
        end
    elseif strcmp(SaveFileType,'fig')
        fprintf(1,'saving figures...\n')
        saveas(Figcopy,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_spec_%d.fig',...
            FileVoc,vv, df,MergeThresh)))
        if PlotRMSFig
            saveas(RMSFig,fullfile(Working_dir_write,sprintf('%s_%d_%d_whocalls_RMS_%d.fig',...
            FileVoc,vv, df,MergeThresh)))
            clf(RMSFig)
        end
    end

    clf(Figcopy)

    % Gather Vocalization production data:
    IndVocStart_all{vv} = IndVocStart;
    IndVocStop_all{vv} = IndVocStop;
    IndVocStartRaw_merged{vv} = IndVocStartRaw_merge_local;
    IndVocStopRaw_merged{vv} = IndVocStopRaw_merge_local;
    IndVocStartPiezo_merged{vv} = IndVocStartPiezo_merge_local;
    IndVocStopPiezo_merged{vv} = IndVocStopPiezo_merge_local;
    RMSRatio_all{vv} = RMSRatio;
    RMSDiff_all{vv} = RMSDiff;
    % Gather vocalization hearing data
    IndHearStart_all{vv} = IndHearStart;
    IndHearStop_all{vv} = IndHearStop;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function Figcopy = copyFigure(plotb, vv)
% Create a copy of the left pannel to save it as a pdf
global AudioLogs Fns_AL;

Figcopy=figure(10);
clf(Figcopy)
Axcopy = subplot(length(AudioLogs)+2,1,1);
plotMic(vv,0,Axcopy);
for ll=1:length(AudioLogs)
    Axcopy = subplot(length(AudioLogs)+2,1,ll+1);
    plotLogger(vv,ll,0,Axcopy);
    yyaxis left
    Axcopy.YColor = 'k';
    yyaxis right
    Axcopy.YColor = 'k';
end
AxcopyLast=subplot(length(AudioLogs)+2,1,length(AudioLogs)+2);
copyobj(plotb{end}.Children,AxcopyLast)
AxcopyLast.XLim = Axcopy.XLim;
AxcopyLast.YTick = 1:(length(AudioLogs)+1);
AxcopyLast.YTickLabel = [Fns_AL; 'Mic'];
AxcopyLast.YLim = [0 length(AudioLogs)+1];
AxcopyLast.YColor = 'k';
AxcopyLast.XColor = 'k';
ylabel('AL ID')
xlabel('Time (ms)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateLoggerEval(vv,ll, NoiseCallListen012_man)
% Save production and hearing data for that logger, merge produced
% vocalizations, update bottom plot on left panel with produced
% vocalizations.
global IndVocStartRaw IndVocStartPiezo IndVocStopRaw IndVocStopPiezo ;
global IndVocStart IndVocStop IndVocStartRaw_merge_local IndVocStopRaw_merge_local;
global IndVocStartPiezo_merge_local IndVocStopPiezo_merge_local;
global IndHearStartRaw IndHearStartPiezo IndHearStopRaw IndHearStopPiezo ;
global IndHearStart IndHearStop
global Fns_AL ColorCode FS Piezo_FS Logger_Spec;
global Fs_env MergeThresh;
global plotb AudioLogs;

% Keep vocalizations that were heared by the bat
IndHearStart{ll} = IndVocStart{ll}(logical(NoiseCallListen012_man==2));
IndHearStop{ll} = IndVocStop{ll}(logical(NoiseCallListen012_man==2));
IndHearStartRaw{vv}{ll} = round(IndHearStart{ll}/Fs_env*FS);
IndHearStopRaw{vv}{ll} = round(IndHearStop{ll}/Fs_env*FS);
if ll<=length(AudioLogs) % This is a vocalization detected on the logger
    IndHearStartPiezo{vv}{ll} = round(IndHearStart{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
    IndHearStopPiezo{vv}{ll} = round(IndHearStop{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
end

% Stash noise and Only keep vocalizations that were produced by the bat
% according to the high RMSRatio value
IndVocStart{ll} = IndVocStart{ll}(logical(NoiseCallListen012_man==1));
IndVocStop{ll} = IndVocStop{ll}(logical(NoiseCallListen012_man==1));
IndVocStartRaw{vv}{ll} = round(IndVocStart{ll}/Fs_env*FS);
IndVocStopRaw{vv}{ll} = round(IndVocStop{ll}/Fs_env*FS);
if ll<=length(AudioLogs) % This is a vocalization detected on the logger
    IndVocStartPiezo{vv}{ll} = round(IndVocStart{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
    IndVocStopPiezo{vv}{ll} = round(IndVocStop{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
end
NV = sum(NoiseCallListen012_man==1);

if NV
    % Now merge detected vocalizations that are within MergeThresh
    Merge01 = (IndVocStart{ll}(2:end)-IndVocStop{ll}(1:(end-1)) < (MergeThresh/1000*Fs_env));
    NV = NV-sum(Merge01);
    IndVocStartRaw_merge_local{ll} = nan(1,NV);
    IndVocStopRaw_merge_local{ll} = nan(1,NV);
    if ll<=length(AudioLogs) % This is a vocalization detected on the logger
        IndVocStartPiezo_merge_local{ll} = nan(1,NV);
        IndVocStopPiezo_merge_local{ll} = nan(1,NV);
    end
    CutInd = find(~Merge01);
    IndVocStartRaw_merge_local{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*FS);
    IndVocStopRaw_merge_local{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*FS);
    if ll<=length(AudioLogs) % This is a vocalization detected on the logger
        IndVocStartPiezo_merge_local{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
        IndVocStopPiezo_merge_local{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
    end
    for cc=1:length(CutInd)
        IndVocStopRaw_merge_local{ll}(cc) = round(IndVocStop{ll}( CutInd(cc))/Fs_env*FS);
        IndVocStartRaw_merge_local{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*FS);
        if ll<=length(AudioLogs) % This is a vocalization detected on the logger
            IndVocStopPiezo_merge_local{ll}(cc) = round(IndVocStop{ll}(CutInd(cc))/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
            IndVocStartPiezo_merge_local{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
        end
    end
    
    % Now plot the onset/offset of each extract on the
    % spectrograms of the loggers
    if ll<=length(AudioLogs)
        Plotll = ll+1;
    elseif ll==length(AudioLogs)+1 %detection on microphone of bat not wearing logger
        Plotll = 1;
    end
    for ii=1:NV
        axes(plotb{Plotll})
        hold on
        yyaxis right
        YLim = get(gca, 'YLim');
        plot(plotb{Plotll},[IndVocStartRaw_merge_local{ll}(ii) IndVocStopRaw_merge_local{ll}(ii)]/FS*1000,...
            [YLim(2)*4/5 YLim(2)*4/5], 'k-', 'LineWidth',2)
        hold off
        axes(plotb{end})
        hold on
        plot(plotb{end},[IndVocStartRaw_merge_local{ll}(ii) IndVocStopRaw_merge_local{ll}(ii)]/FS*1000,...
            [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
        set(gca,'xlim',[0 Logger_Spec{1}.to(end)*1000])
        drawnow;
        hold off
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMic(vv,FigN, axpl)
global Raw_Spec  Nvoc df DataFiles FHigh_spec;
global DB_noise plotmich plotmicevalh plotlogevalh;
global Amp_env_Mic Fs_env sliderLefth

if FigN==1
    axpl=plotmich;
elseif FigN==3
    axpl=plotmicevalh;
elseif FigN ==2
    axpl=plotlogevalh;
elseif FigN == 0
    if nargin<3
        error('Provide the axis handle to wear the figure should be plotted')
    end
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
if FigN==0
    set(gca,'xlim',[0 Raw_Spec.to(end)*1000],'ylim',[0 FHigh_spec],...
    'ytick',[0 FHigh_spec/2 FHigh_spec])
    xlabel(' ')
    set(gca, 'XTick',[],'XTickLabel',{})
else
    set(gca,'xlim',[0 Raw_Spec.to(end)*1000],'ylim',[0 FHigh_spec+5e3],...
    'ytick',[0 FHigh_spec/2 FHigh_spec])
    xlabel('time (ms)');
end
ylabel('Frequency');

if FigN==1 || FigN==0
    hold on
    yyaxis right
    plot(axpl,(1:length(Amp_env_Mic))/Fs_env*1000, Amp_env_Mic, 'r-', 'LineWidth',2)
    %ylabel(sprintf('Amp\nMic'))
    title(sprintf('Voc %d/%d Set %d/%d',vv,Nvoc, df,length(DataFiles)))
    xlabel(' ')
    set(gca, 'XTick',[],'XTickLabel',{})
    set(gca,'ytick','')
    drawnow;
    hold off;
    if FigN==1
        yyaxis left
        axpl.YColor='w';
        set(sliderLefth,'SliderStep', [1/(length(Amp_env_Mic)-1), 10/(length(Amp_env_Mic)-1)], ...
        'Min', 1, 'Max', length(Amp_env_Mic), 'Value', 1)
    end
else
    yyaxis left
    axpl.YColor='w';
    drawnow;
    hold off;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function plotLogger(vv,ll,FigN, axpl)
global LowPassLogVoc Piezo_FS Fns_AL DB_noise FHigh_spec_Logger;
global Amp_env_LowPassLogVoc Amp_env_HighPassLogVoc Fs_env AudioLogs;
global plotb plotlogevalh plotevalh Logger_Spec;

val=0;
if FigN==1
    axes(plotevalh);
    cla(plotevalh ,'reset')
    axpl=plotb{ll+1};
    [Logger_Spec{ll}.to, Logger_Spec{ll}.fo, Logger_Spec{ll}.logB] = ...
        spec_only_bats_gui(LowPassLogVoc{ll}, Piezo_FS.(Fns_AL{ll})(vv),...
        DB_noise, FHigh_spec_Logger);
elseif FigN==3
    axpl=plotlogevalh;
    val=3e3;
elseif nargin<3
    error('Provide the axis handle to wear the figure should be plotted')
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
if FigN==3
    axpl.XColor='w';
else
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
    if FigN==0
        ylabel(sprintf('Amp\n%s',Fns_AL{ll}([1:3 7:end])))
    end
    hold off
end

hold off;
if FigN~=0
    zoom on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newmessage(mstring)
global message mh string_handle2;

message={message{2};message{3};message{4};message{5};message{6};...
    mstring};

set(mh,string_handle2,message);
