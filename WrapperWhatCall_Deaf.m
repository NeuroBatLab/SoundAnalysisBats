%% This script run through files that have been manually cured to apply what_calls.m
%% These are specific to the dataset and computer
BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
BaseCodeDir = 'C:\Users\tobias\Documents\GitHub\';
% WorkingDir = 'C:\Users\tobias\Documents\DeafWhoWorkDir\';
WhoLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWho.txt');
WhatLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWhat.txt');

%% Paths to code
% you should have pulled from github the last versions of
% LMC,LoggerDataProcessing and SoundAnalysisBats
addpath(genpath(fullfile(BaseCodeDir,'LMC')))
addpath(genpath(fullfile(BaseCodeDir, 'LoggerDataProcessing')))
addpath(genpath(fullfile(BaseCodeDir,'SoundAnalysisBats')))

%% Parameters of WhatCalls



%% Get the name of the next experiment that needs to be manually curated
ee=0;
% these are the list of files to run (replace with WhoLog)
if ~exist(WhoLog, 'file')
    error('Cannot find the list of file to run in: %s \n',WhoLog);
else
    FidWho = fopen(WhoLog, 'r');
    Header = textscan(FidWho,'%s\t%s\t%s\t%s\t%s\n',1);
    DoneListDetect = textscan(FidWho,'%s\t%s\t%s\t%.1f\t%d');
    fclose(FidWho);
end

% This is the list of files that has been processed (to replace with the new
% WhatLog)
if ~exist(WhoLog, 'file')
    FidWho = fopen(WhoLog, 'a');
    fprintf(FidWho, 'Subject\tDate\tTime\n');
    DoneListWho = [];
else
    FidWho = fopen(WhoLog, 'r');
    Header = textscan(FidWho,'%s\t%s\t%s\n',1);
    DoneListWho = textscan(FidWho,'%s\t%s\t%s');
    fclose(FidWho);
    FidWho = fopen(WhoLog, 'a');
end


%% Run what_calls on the selected data
% 1- Do a loop through all the files that need to run
% 2- Indisde the loop Check that the files listed in WhoLog are not already listed in
% WhatLog (if they are listed in WhatLog it means they've already been
% processed by WhatCalls and we might not want to run them and overwrite
% the data)
% 3- Once you've identified an experimental date that has not been run
% through whatcalls, apply whatcalls
% what_calls(Logger_dir, Date, ExpStartTime)

if ~exist(WhoLog, 'file')
    error('Cannot find the list of file to run in: %s \n',ExpLog);
else
    FidExp = fopen(WhoLog, 'r');
end
    
while checkFiles
    if ~isempty(WhatLog)
        Done = sum(contains(DoneListWho{1},BatsID) .* contains(DoneListWho{2},Date) .* contains(DoneListWho{3},ExpStartTime)); 
    else 
        Done = 0;
    end 

    if Done
        checkFiles=1;
    elseif ~Done && ~isempty(WhatLog)
        if WhatLog
            ParamFile = dir(fullfile(BaseDataDir,['20' Date],'audio',sprintf('%s_%s_%s*RecOnly_param.txt', BatsID, Date, ExpStartTime)));
            Logger_dir = fullfile(ParamFile.folder(1:(strfind(ParamFile.folder, 'audio')-1)), 'audiologgers');
            checkFiles=0;
        else
            checkFiles=1;
         end
    elseif ~Done && isempty(WhatLog)
        fprintf(1, '   -> Starting new session\n')
    end
end

fprintf(1, 'This is where I would call what_calls with Logger_dir being:\n %s\n Date: %s\n ExpStartTime:%s\n', Logger_diir, Date, ExpStartTime)
pause()

%% OLD CODE TO DELETE
%Grabbing a new Session
fprintf(1,'Grabbing the session...');
NExpe = length(DoneListDetect{1});
checkSession=1;
AlliOk=[];
while checkSession && ee<=NExpe
    ee=ee+1;
    BatsID = DoneListDetect{1}{ee};
    Date = DoneListDetect{2}{ee};
    ExpStartTime = DoneListDetect{3}{ee};
    ParamFile = dir(fullfile(BaseDataDir,['20' Date],'audio',sprintf('%s_%s_%s*RecOnly_param.txt', BatsID, Date, ExpStartTime)));
    fprintf(1, '\n\n\n Date: %s, experiment %d/%d\n%s\n', Date,ee,NExpe,ParamFile.name)
    
    % Check that the file was not already set aside or done
    if ~isempty(DoneListWho)
        Done = sum(contains(DoneListWho{1},BatsID) .* contains(DoneListWho{2},Date) .* contains(DoneListWho{3},ExpStartTime));
    else
        Done=0;
    end
    
    if ~isempty(ListAlliOk)
        AlliOkInd = find(contains(ListAlliOk{1},BatsID) .* contains(ListAlliOk{2},Date) .* contains(ListAlliOk{3},ExpStartTime));
        AlliOk = ListAlliOk{4}(AlliOkInd);
    else
        AlliOk=[];
    end
    
    if Done
        fprintf(1, '   -> Data already processed\n')
        checkSession=1;
    elseif ~Done && ~isempty(AlliOk)
        if AlliOk
            fprintf(1, '   -> Starting from where we left on this session\n')
            Logger_dir = fullfile(ParamFile.folder(1:(strfind(ParamFile.folder, 'audio')-1)), 'audiologgers');
            checkSession=0;
        else
            fprintf(1, '   -> Session flagged as not alligned correctly\n')
            checkSession=1;
        end
    elseif ~Done && isempty(AlliOk)
        fprintf(1, '   -> Starting new session\n')
        % Check that the clocks drifts were correctly corrected
        fprintf(1,'*** Check the clock drift correction of the logger ***\n')
        Logger_dir = fullfile(ParamFile.folder(1:(strfind(ParamFile.folder, 'audio')-1)), 'audiologgers');
        LoggersDir = dir(fullfile(Logger_dir, 'logger*'));
        Check = zeros(length(LoggersDir)+1,1);
        for ll=1:length(LoggersDir)
            FigCD = open(fullfile(LoggersDir(ll).folder, LoggersDir(ll).name,'extracted_data','CD_correction0.fig'));
            failsafe=1;
            while failsafe
                Check(ll) = input('Is everything ok? (yes ->1, No -> 0): ');
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
        AllignmentPath = fullfile(ParamFile.folder,sprintf('%s_%s_CD_correction_audio_piezo.fig', Date, ExpStartTime));
        FigAP = open(AllignmentPath);
        failsafe=1;
        while failsafe
            Check(length(LoggersDir)+1) = input('Is everything ok? (yes ->1, No -> 0): ');
            failsafe=0;
            if isempty(Check(ll))
                failsafe=1;
                disp('Entry is empty, please repeat!')
            end
        end
        fprintf('\n')
        close(FigAP)
        if any(~Check)
            AlliOk=0;
            fprintf(FidAlli, '%s\t%s\t%s\t%d\n',BatsID,Date,ExpStartTime,AlliOk);
            fprintf(1,'\n****** Error in allignement reported ******\n')
            checkSession=1;
        else
            AlliOk=1;
            fprintf(FidAlli, '%s\t%s\t%s\t%d\n',BatsID,Date,ExpStartTime,AlliOk);
            fprintf(1,'\n****** Allignement reported as good! ******\n')
            checkSession=0;
        end
    end
    
end
% This is the name to the experiment that needs to be analyzed
Filepath = fullfile(ParamFile.folder, ParamFile.name);
fprintf(1, '\n\n\n Date: %s, experiment %d/%d\n%s\n', Date,ee,NExpe,ParamFile.name)


