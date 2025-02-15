%% Script calling the Gui CallCuraGui (code in CallCurafkt) for manual
% curation of vocalizations extracted from microphone and loggers. It is
% the faithful Gui version of Who_calls_playless

% This code gets ready the Gui input

close all;
clear all;

global WorkingDir Filepath df Logger_dir
global VolDenominatorLogger VolFactorMic Manual MergeThresh;
global Force_Save_onoffsets_mic SaveFileType;
global CheckMicChannel UseOld  PlotRMSFig FhGUI;
global evalMich evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h submith string_handle2;
global evalLog6h evalLog7h evalLog8h starth evalLog plotb string_handle;
global filenameh noCallh maybeCallh redoh redoEditVoch redoEditSeth sliderLefth sliderRighth;
global traineeNameh
global playMich playLog1h playLog2h playLog3h playLog4h playLog5h playLog6h;
global playLog7h playLog8h playMicEvalh playLogEvalh redo checkboxh playLog9h;
global playLog10h evalLog9h evalLog10h playb message mh;
global plotmich plotmicevalh plotlog1h plotlog2h plotlog3h plotlog4h plotlog5h;
global plotlog6h plotlog7h plotlog8h plotlog9h plotlog10h plotevalh plotlogevalh;



%% These are specific to the dataset and computer
BaseDataDir = 'Z:\users\JulieE\LMC\LMC_CoEd\audio';
% BaseDataDir = '/Volumes/JulieE4T/LMC_CoEd/audio';
BaseCodeDir = 'C:\Users\tobias\Documents\GitHub\';
% BaseCodeDir = '/Users/elie/Documents/CODE/';
WorkingDir = 'C:\Users\tobias\Documents\TrainWhoWorkDir\';
% WorkingDir = '/Users/elie/Documents/TrainWhoWorkDir/';

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
%Grabbing a new Session
fprintf(1,'Grabbing the session...');

% This is the name to the experiment that needs to be analyzed
BatsID = 'CoEd';
Date = '190612';
ExpStartTime = '1438';
ParamFile = dir(fullfile(BaseDataDir,['20' Date],sprintf('%s_%s_%s*RecOnly_param.txt', BatsID, Date, ExpStartTime)));
Filepath = fullfile(ParamFile.folder, ParamFile.name);
fprintf(1, '\n\n\n Date: %s, experiment:\n%s\n', Date,ParamFile.name)
Logger_dir = fullfile(ParamFile.folder(1:(strfind(ParamFile.folder, 'audio')-1)), 'logger',['20' Date]);
% this is the set # we want to focus on
df = 4;
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
    FhGUI=CallCura_TrainGui;
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
