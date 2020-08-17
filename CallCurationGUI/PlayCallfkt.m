function PlayCallfkt(action)
global playll vv;

switch action
    
    case 'PlayMic'
        Player=prepPlayMic;
        play(Player)
        updateSliderLeft
        
    case 'PlayMicEval'
        Player=prepPlayMic;
        play(Player)
        updateSliderRight
        
    case 'PlayLog1'
        [Player]=getplaylogger(1,vv);
        play(Player)
        updateSliderLeft
        
    case 'PlayLog2'
        [Player]=getplaylogger(2,vv);
        play(Player)
        updateSliderLeft
        
    case 'PlayLog3'
        [Player]=getplaylogger(3,vv);
        play(Player)
        updateSliderLeft
        
    case 'PlayLog4'
        [Player]=getplaylogger(4,vv);
        play(Player)
        updateSliderLeft
        
    case 'PlayLog5'
        [Player]=getplaylogger(5,vv);
        play(Player)
        updateSliderLeft
        
    case 'PlayLog6'
        [Player]=getplaylogger(6,vv);
        play(Player)
        updateSliderLeft
        
    case 'PlayLog7'
        [Player]=getplaylogger(7,vv);
        play(Player)
        updateSliderLeft
        
    case 'PlayLog8'
        [Player]=getplaylogger(8,vv);
        play(Player)
        updateSliderLeft
        
    case PlayLogEval
        [Player]=getplaylogger(playll,vv);
        play(Player)
        updateSliderRight
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Player]=getplaylogger(LogN,vv)
global Piezo_wave Fns_AL VolDenominatorLogger Piezo_FS;

Player= audioplayer((Piezo_wave.(Fns_AL{LogN}){vv}-...
    nanmean(Piezo_wave.(Fns_AL{LogN}){vv}))/(VolDenominatorLogger*nanstd(Piezo_wave.(Fns_AL{LogN}){vv})), ...
    Piezo_FS.(Fns_AL{LogN})(vv));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Player=prepPlayMic
global Raw_wave_nn sos_raw_band_listen VolFactorMic FS Raw_listen;

Raw_listen = filtfilt(sos_raw_band_listen,1,Raw_wave_nn);
SampleMic = resample((Raw_listen - mean(Raw_listen))/(std(Raw_listen)/VolFactorMic),FS/4,FS);
Player= audioplayer(SampleMic, FS/4,24);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateSliderLeft
global sliderLefth Raw_wave_nn FS;

startpos=get(sliderLefth,'Value');
tic;
while toc<=length(Raw_wave_nn)/FS-(startpos/FS)
    set(sliderLefth,'Value', round(toc*FS)+startpos)
end
set(sliderLefth,'Value', 1)
pause(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateSliderRight
global sliderRighth Raw_wave_nn FS;

startpos=get(sliderRighth,'Value');
tic;
while toc<=length(Raw_wave_nn)/FS-(startpos/FS)
    set(sliderRighth,'Value', round(toc*FS)+startpos)
end
set(sliderRighth,'Value', 1)
pause(1);

