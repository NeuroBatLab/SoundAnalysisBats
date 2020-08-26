clear all;
close all;

global BaseDataDir WorkingDir VolDenominatorLogger VolFactorMic Manual MergeThresh;
global CheckMicChannel pnames NExpe UseOld FhGUI ee DoneListDetect DoneListWho;
global evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h submith string_handle2;
global evalLog6h evalLog7h evalLog8h starth evalb plotb string_handle;
global filenameh noCallh redoh redoEditVoch redoEditSeth sliderLefth sliderRighth;
global playMich playLog1h playLog2h playLog3h playLog4h playLog5h playLog6h;
global playLog7h playLog8h playMicEvalh playLogEvalh redo checkboxh playLog9h;
global playLog10h evalLog9h evalLog10h playb message mh;
global plotmich plotmicevalh plotlog1h plotlog2h plotlog3h plotlog4h plotlog5h;
global plotlog6h plotlog7h plotlog8h plotlog9h plotlog10h plotevalh plotlogevalh;




BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';
BaseCodeDir = 'C:\Users\BatLab\Documents\GitHub\';
WorkingDir = 'C:\Users\BatLab\Documents\DeafWhoWorkDir\';

addpath(genpath(fullfile(BaseCodeDir,'LMC')))
addpath(genpath(fullfile(BaseCodeDir, 'LoggerDataProcessing')))
addpath(genpath(fullfile(BaseCodeDir,'SoundAnalysisBats')))

DatesDir = dir(fullfile(BaseDataDir,'20*'));
NDates = length(DatesDir);
ExpLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSal.txt');
WhoLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWho.txt');

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

%general variables
NExpe = length(DoneListDetect{1});
% optional parameter: Factor_RMS_Mic, Factor by which the RMS of the
% band-pass filtered baseline signal is multiplied to obtained the
% threshold of vocalization detection on Microphone
VolDenominatorLogger=5;
VolFactorMic=0.5;
pnames = {'Factor_RMS_Mic','Working_dir','Force_Save_onoffsets_mic','SaveFileType'};
MergeThresh=200;
Manual=1;
UseOld=1;
CheckMicChannel=0;
redo=0;

message={'';'';'';'';'';''};
Use_AppDesigner=1;
if Use_AppDesigner
    string_handle='Text';
    string_handle2='Value';
    FhGUI=CallCuraGui_App;
    starth=FhGUI.start;
    noCallh=FhGUI.noCall;
    filenameh=FhGUI.filename;
    mh=FhGUI.message;
    plotmich=FhGUI.plotMic;
    plotmicevalh=FhGUI.plotMicEval;
    plotlog1h=FhGUI.plotLog1;
    plotlog2h=FhGUI.plotLog2;
    plotlog3h=FhGUI.plotLog3;
    plotlog4h=FhGUI.plotLog4;
    plotlog5h=FhGUI.plotLog5;
    plotlog6h=FhGUI.plotLog6;
    plotlog7h=FhGUI.plotLog7;
    plotlog8h=FhGUI.plotLog8;
    plotlog9h=FhGUI.plotLog9;
    plotlog10h=FhGUI.plotLog10;
    plotevalh=FhGUI.plotEval;
    plotlogevalh=FhGUI.plotLogEval;
    sliderLefth=FhGUI.sliderLeft;
    sliderRighth=FhGUI.sliderRight;
    redoh=FhGUI.redo;
    redoEditVoch=FhGUI.redoEditVoc;
    redoEditSeth=FhGUI.redoEditSet;
    submith=FhGUI.submit;
    checkboxh=FhGUI.checkbox;
    evalLog1h=FhGUI.evalLog1;
    evalLog2h=FhGUI.evalLog2;
    evalLog3h=FhGUI.evalLog3;
    evalLog4h=FhGUI.evalLog4;
    evalLog5h=FhGUI.evalLog5;
    evalLog6h=FhGUI.evalLog6;
    evalLog7h=FhGUI.evalLog7;
    evalLog8h=FhGUI.evalLog8;
    evalLog9h=FhGUI.evalLog9;
    evalLog10h=FhGUI.evalLog10;
    playMich=FhGUI.playMic;
    playMicEvalh=FhGUI.playMicEval;
    playLogEvalh=FhGUI.playLogEval;
    playLog1h=FhGUI.playLog1;
    playLog2h=FhGUI.playLog2;
    playLog3h=FhGUI.playLog3;
    playLog4h=FhGUI.playLog4;
    playLog5h=FhGUI.playLog5;
    playLog6h=FhGUI.playLog6;
    playLog7h=FhGUI.playLog7;
    playLog8h=FhGUI.playLog8;
    playLog9h=FhGUI.playLog9;
    playLog10h=FhGUI.playLog10;
else
    string_handle='String';
     string_handle2='String';
    FhGUI=CallCuraGui;
    starth=findobj(FhGUI,'tag','start');
    noCallh=findobj(FhGUI,'tag','noCall');
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
    submith=findobj(FhGUI,'tag','submit');
    checkboxh=findobj(FhGUI,'tag','checkbox');
    evalLog1h=findobj(FhGUI,'tag','evalLog1');
    evalLog2h=findobj(FhGUI,'tag','evalLog2');
    evalLog3h=findobj(FhGUI,'tag','evalLog3');
    evalLog4h=findobj(FhGUI,'tag','evalLog4');
    evalLog5h=findobj(FhGUI,'tag','evalLog5');
    evalLog6h=findobj(FhGUI,'tag','evalLog6');
    evalLog7h=findobj(FhGUI,'tag','evalLog7');
    evalLog8h=findobj(FhGUI,'tag','evalLog8');
    evalLog9h=findobj(FhGUI,'tag','evalLog9');
    evalLog10h=findobj(FhGUI,'tag','evalLog10');
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
end

evalb{1}=evalLog1h;
evalb{2}=evalLog2h;
evalb{3}=evalLog3h;
evalb{4}=evalLog4h;
evalb{5}=evalLog5h;
evalb{6}=evalLog6h;
evalb{7}=evalLog7h;
evalb{8}=evalLog8h;
evalb{9}=evalLog9h;
evalb{10}=evalLog10h;


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

set([evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h submith evalLog6h...
    evalLog7h evalLog8h evalLog9h evalLog10h noCallh redoh redoEditVoch...
    redoEditSeth sliderLefth sliderRighth playMich playLog1h playLog2h...
    playLog3h playLog4h playLog5h playLog6h playLog7h playLog8h playLog9h playLog10h...
    playMicEvalh playLogEvalh checkboxh],'enable','off')
