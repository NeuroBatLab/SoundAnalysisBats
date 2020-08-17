clear all;
close all;

global VolDenominatorLogger VolFactorMic Manual MergeThresh;
global CheckMicChannel pnames NExpe UseOld FhGUI;
global evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h submith;
global evalLog6h evalLog7h evalLog8h leftPloth rightPloth starth evalb;
global filenameh noCallh redoh redoEditVoch redoEditSeth sliderLefth sliderRighth;
global playMich playLog1h playLog2h playLog3h playLog4h playLog5h playLog6h;
global playLog7h playLog8h playMicEvalh playLogEvalh redo;


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

%GUI variables
FhGUI=CallCuraGui;
%set (FhGUI, 'position', [6.8000   36.4615  153.4000   43.7692]);
starth=findobj(FhGUI,'tag','start');
noCallh=findobj(FhGUI,'tag','noCall');
filenameh=findobj(FhGUI,'tag','filename');
leftPloth=findobj(FhGUI,'tag','leftPlot');
rightPloth=findobj(FhGUI,'tag','rightPlot');
sliderLefth=findobj(FhGUI,'tag','sliderLeft');
sliderRighth=findobj(FhGUI,'tag','sliderRight');
redoh=findobj(FhGUI,'tag','redo');
redoEditVoch=findobj(FhGUI,'tag','redoEditVoc');
redoEditSeth=findobj(FhGUI,'tag','redoEditSet');
submith=findobj(FhGUI,'tag','submit');

evalLog1h=findobj(FhGUI,'tag','evalLog1');
evalLog2h=findobj(FhGUI,'tag','evalLog2');
evalLog3h=findobj(FhGUI,'tag','evalLog3');
evalLog4h=findobj(FhGUI,'tag','evalLog4');
evalLog5h=findobj(FhGUI,'tag','evalLog5');
evalLog6h=findobj(FhGUI,'tag','evalLog6');
evalLog7h=findobj(FhGUI,'tag','evalLog7');
evalLog8h=findobj(FhGUI,'tag','evalLog8');
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

evalb{1}=evalLog1h;
evalb{2}=evalLog2h;
evalb{3}=evalLog3h;
evalb{4}=evalLog4h;
evalb{5}=evalLog5h;
evalb{6}=evalLog6h;
evalb{7}=evalLog7h;
evalb{8}=evalLog8h;
set([evalLog1h evalLog2h evalLog3h evalLog4h evalLog5h submith evalLog6h...
    evalLog7h evalLog8h noCallh redoh redoEditVoch...
    redoEditSeth sliderLefth sliderRighth playMich playLog1h playLog2h...
    playLog3h playLog4h playLog5h playLog6h playLog7h playLog8h playMicEvalh...
    playLogEvalh],'enable','off')
