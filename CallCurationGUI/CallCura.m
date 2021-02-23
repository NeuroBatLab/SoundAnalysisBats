%% Script calling the Gui CallCuraGui (code in CallCurafkt) for manual
% curation of vocalizations extracted from microphone and loggers. It is
% the faithful Gui version of Who_calls_playless

% This code gets ready the Gui input

close all;
clear all;

global WorkingDir Filepath Logger_dir FidWho
global VolDenominatorLogger VolFactorMic Manual MergeThresh;
global Force_Save_onoffsets_mic SaveFileType;
global CheckMicChannel UseOld  PlotRMSFig FhGUI;
global evalMich evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h submith string_handle2;
global evalLog6h evalLog7h evalLog8h starth evalLog plotb string_handle;
global filenameh noCallh maybeCallh redoh redoEditVoch redoEditSeth sliderLefth sliderRighth;
global playMich playLog1h playLog2h playLog3h playLog4h playLog5h playLog6h;
global playLog7h playLog8h playMicEvalh playLogEvalh redo checkboxh playLog9h;
global playLog10h evalLog9h evalLog10h playb message mh traineeNameh;
global plotmich plotmicevalh plotlog1h plotlog2h plotlog3h plotlog4h plotlog5h;
global plotlog6h plotlog7h plotlog8h plotlog9h plotlog10h plotevalh plotlogevalh;



%% These are specific to the dataset and computer
BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
BaseCodeDir = 'C:\Users\tobias\Documents\GitHub\';
WorkingDir = 'C:\Users\tobias\Documents\DeafWhoWorkDir\';
ExpLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSal.txt');
WhoLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWho.txt');
AlliLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalAllignement.txt');

%% Paths to code
% you should have pulled from github the last versions of
% LMC,LoggerDataProcessing and SoundAnalysisBats
addpath(genpath(fullfile(BaseCodeDir,'LMC')))
addpath(genpath(fullfile(BaseCodeDir, 'LoggerDataProcessing')))
addpath(genpath(fullfile(BaseCodeDir,'SoundAnalysisBats')))

%% Parameters of WhoCalls

VolDenominatorLogger=5;
VolFactorMic=0.5;

Force_Save_onoffsets_mic = 0; % Currently useless
SaveFileType = 'pdf';

MergeThresh=200;
Manual=1;
UseOld=1;% Set to 1 if you want to append to existing manually curated data
CheckMicChannel=0; % If set to 1 authorize evaluation of Microphone channel (case of an individual without a collar)
PlotRMSFig = 0; % Set to 1 if you want to plot and save the figure showing the RMS

%% Get the name of the next experiment that needs to be manually curated
ee=0;
if ~exist(ExpLog, 'file')
    error('Cannot find the list of file to run in: %s \n',ExpLog);
else
    FidExp = fopen(ExpLog, 'r');
    Header = textscan(FidExp,'%s\t%s\t%s\t%s\t%s\n',1);
    DoneListDetect = textscan(FidExp,'%s\t%s\t%s\t%.1f\t%d');
    fclose(FidExp);
end

if ~exist(WhoLog, 'file')
    FidWho = fopen(WhoLog, 'a');
    fprintf(FidWho, 'Subject\tDate\tTime\tNCalls\n');
    DoneListWho = [];
else
    FidWho = fopen(WhoLog, 'r');
    Header = textscan(FidWho,'%s\t%s\t%s\t%s\n',1);
    DoneListWho = textscan(FidWho,'%s\t%s\t%s\t%d');
    fclose(FidWho);
    FidWho = fopen(WhoLog, 'a');
end

if ~exist(AlliLog, 'file')
    FidAlli = fopen(AlliLog, 'a');
    fprintf(FidAlli, 'Subject\tDate\tTime\tAllignement\n');
    ListAlliOk = [];
else
    FidAlli = fopen(AlliLog, 'r');
    Header = textscan(FidAlli,'%s\t%s\t%s\t%s\n',1);
    ListAlliOk = textscan(FidAlli,'%s\t%s\t%s\t%d');
    fclose(FidAlli);
    FidAlli = fopen(AlliLog, 'a');
end

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

%% Starting the GUI
% Initializing variables
redo=0;
message={'';'';'';'';'';''};
% Use_AppDesigner=0;
% if Use_AppDesigner
%     string_handle='Text';
%     string_handle2='Value';
%     FhGUI=CallCuraGui_App;
%     starth=FhGUI.start;
%     noCallh=FhGUI.noCall;
%     filenameh=FhGUI.filename;
%     mh=FhGUI.message;
%     plotmich=FhGUI.plotMic;
%     plotmicevalh=FhGUI.plotMicEval;
%     plotlog1h=FhGUI.plotLog1;
%     plotlog2h=FhGUI.plotLog2;
%     plotlog3h=FhGUI.plotLog3;
%     plotlog4h=FhGUI.plotLog4;
%     plotlog5h=FhGUI.plotLog5;
%     plotlog6h=FhGUI.plotLog6;
%     plotlog7h=FhGUI.plotLog7;
%     plotlog8h=FhGUI.plotLog8;
%     plotlog9h=FhGUI.plotLog9;
%     plotlog10h=FhGUI.plotLog10;
%     plotevalh=FhGUI.plotEval;
%     plotlogevalh=FhGUI.plotLogEval;
%     sliderLefth=FhGUI.sliderLeft;
%     sliderRighth=FhGUI.sliderRight;
%     redoh=FhGUI.redo;
%     redoEditVoch=FhGUI.redoEditVoc;
%     redoEditSeth=FhGUI.redoEditSet;
%     submith=FhGUI.submit;
%     checkboxh=FhGUI.checkbox;
%     evalMich = FhGUI.evalMic;
%     evalLog1h=FhGUI.evalLog1;
%     evalLog2h=FhGUI.evalLog2;
%     evalLog3h=FhGUI.evalLog3;
%     evalLog4h=FhGUI.evalLog4;
%     evalLog5h=FhGUI.evalLog5;
%     evalLog6h=FhGUI.evalLog6;
%     evalLog7h=FhGUI.evalLog7;
%     evalLog8h=FhGUI.evalLog8;
%     evalLog9h=FhGUI.evalLog9;
%     evalLog10h=FhGUI.evalLog10;
%     playMich=FhGUI.playMic;
%     playMicEvalh=FhGUI.playMicEval;
%     playLogEvalh=FhGUI.playLogEval;
%     playLog1h=FhGUI.playLog1;
%     playLog2h=FhGUI.playLog2;
%     playLog3h=FhGUI.playLog3;
%     playLog4h=FhGUI.playLog4;
%     playLog5h=FhGUI.playLog5;
%     playLog6h=FhGUI.playLog6;
%     playLog7h=FhGUI.playLog7;
%     playLog8h=FhGUI.playLog8;
%     playLog9h=FhGUI.playLog9;
%     playLog10h=FhGUI.playLog10;
% else
    string_handle='String';
    string_handle2='String';
    FhGUI=CallCuraGui;
    starth=findobj(FhGUI,'tag','start');
    noCallh=findobj(FhGUI,'tag','noCall');
    maybeCallh=findobj(FhGUI,'tag','maybeCall');
    filenameh=findobj(FhGUI,'tag','filename');
    mh=findobj(FhGUI,'tag','message');
    plotmich=findobj(FhGUI,'tag','plotMic');
    plotmicevalh=findobj(FhGUI,'tag','plotMicEval');
    plotlog1h=findobj(FhGUI,'tag','plotLog1');
    plotlog2h=findobj(FhGUI,'tag','plotLog2');
    plotlog3h=findobj(FhGUI,'tag','plotLog3');
    plotlog4h=findobj(FhGUI,'tag','plotLog4');
    plotlog5h=findobj(FhGUI,'tag','plotLog5');
    plotlog6h=findobj(FhGUI,'tag','plotLog6');
    plotlog7h=findobj(FhGUI,'tag','plotLog7');
    plotlog8h=findobj(FhGUI,'tag','plotLog8');
    plotlog9h=findobj(FhGUI,'tag','plotLog9');
    plotlog10h=findobj(FhGUI,'tag','plotLog10');
    plotevalh=findobj(FhGUI,'tag','plotEval');
    plotlogevalh=findobj(FhGUI,'tag','plotLogEval');
    sliderLefth=findobj(FhGUI,'tag','sliderLeft');
    sliderRighth=findobj(FhGUI,'tag','sliderRight');
    redoh=findobj(FhGUI,'tag','redo');
    redoEditVoch=findobj(FhGUI,'tag','redoEditVoc');
    redoEditSeth=findobj(FhGUI,'tag','redoEditSet');
    traineeNameh=findobj(FhGUI,'tag','trainee_name');
    submith=findobj(FhGUI,'tag','submit');
    checkboxh=findobj(FhGUI,'tag','checkbox');
    evalLog = cell(10,1);
    evalLog{1}=findobj(FhGUI,'tag','evalLog1');
    evalLog{2}=findobj(FhGUI,'tag','evalLog2');
    evalLog{3}=findobj(FhGUI,'tag','evalLog3');
    evalLog{4}=findobj(FhGUI,'tag','evalLog4');
    evalLog{5}=findobj(FhGUI,'tag','evalLog5');
    evalLog{6}=findobj(FhGUI,'tag','evalLog6');
    evalLog{7}=findobj(FhGUI,'tag','evalLog7');
    evalLog{8}=findobj(FhGUI,'tag','evalLog8');
    evalLog{9}=findobj(FhGUI,'tag','evalLog9');
    evalLog{10}=findobj(FhGUI,'tag','evalLog10');
    evalMich =findobj(FhGUI,'tag','evalMic');
    playMich=findobj(FhGUI,'tag','playMic');
    playMicEvalh=findobj(FhGUI,'tag','playMicEval');
    playLogEvalh=findobj(FhGUI,'tag','playLogEval');
    playLog1h=findobj(FhGUI,'tag','playLog1');
    playLog2h=findobj(FhGUI,'tag','playLog2');
    playLog3h=findobj(FhGUI,'tag','playLog3');
    playLog4h=findobj(FhGUI,'tag','playLog4');
    playLog5h=findobj(FhGUI,'tag','playLog5');
    playLog6h=findobj(FhGUI,'tag','playLog6');
    playLog7h=findobj(FhGUI,'tag','playLog7');
    playLog8h=findobj(FhGUI,'tag','playLog8');
    playLog9h=findobj(FhGUI,'tag','playLog9');
    playLog10h=findobj(FhGUI,'tag','playLog10');
% end

plotb{1}=plotmich;
plotb{2}=plotlog1h;
plotb{3}=plotlog2h;
plotb{4}=plotlog3h;
plotb{5}=plotlog4h;
plotb{6}=plotlog5h;
plotb{7}=plotlog6h;
plotb{8}=plotlog7h;
plotb{9}=plotlog8h;
plotb{10}=plotlog9h;
plotb{11}=plotlog10h;
plotb{12}=plotevalh;

playb{1}=playLog1h;
playb{2}=playLog2h;
playb{3}=playLog3h;
playb{4}=playLog4h;
playb{5}=playLog5h;
playb{6}=playLog6h;
playb{7}=playLog7h;
playb{8}=playLog8h;
playb{9}=playLog9h;
playb{10}=playLog10h;

set([evalMich evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h submith evalLog6h...
    evalLog7h evalLog8h evalLog9h evalLog10h noCallh maybeCallh redoh redoEditVoch...
    redoEditSeth sliderLefth sliderRighth playMich playLog1h playLog2h...
    playLog3h playLog4h playLog5h playLog6h playLog7h playLog8h playLog9h playLog10h...
    playMicEvalh playLogEvalh checkboxh],'enable','off')
