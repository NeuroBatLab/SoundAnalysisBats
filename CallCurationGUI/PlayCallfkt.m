function PlayCallfkt(action)
global playll vv sliderLefth sliderRighth;

switch action
    
    case 'PlayMic'
        Player=prepPlayMic;
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayMicEval'
        Player=prepPlayMic;
        startpos=get(sliderRighth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderRight(startpos)
        
    case 'PlayLog1'
        [Player]=getplaylogger(1,vv);
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayLog2'
        [Player]=getplaylogger(2,vv);
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayLog3'
        [Player]=getplaylogger(3,vv);
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayLog4'
        [Player]=getplaylogger(4,vv);
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayLog5'
        [Player]=getplaylogger(5,vv);
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayLog6'
        [Player]=getplaylogger(6,vv);
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayLog7'
        [Player]=getplaylogger(7,vv);
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayLog8'
        [Player]=getplaylogger(8,vv);
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayLog9'
        [Player]=getplaylogger(9,vv);
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayLog10'
        [Player]=getplaylogger(10,vv);
        startpos=get(sliderLefth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderLeft(startpos)
        
    case 'PlayLogEval'
        [Player]=getplaylogger(playll,vv);
        startpos=get(sliderRighth,'Value');
        play(Player,round((startpos/1e3)*Player.SampleRate))
        updateSliderRight(startpos)
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
function updateSliderLeft(startpos)
global sliderLefth Raw_wave_nn FS;

tic;
while toc*1e3<=(length(Raw_wave_nn)/FS)*1e3-(startpos)
    set(sliderLefth,'Value', round(toc*1e3)+startpos)
    drawnow;
    pause(.02);
end
set(sliderLefth,'Value', 1)
pause(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateSliderRight(startpos)
global sliderRighth Raw_wave_nn FS;

tic;
while toc*1e3<=(length(Raw_wave_nn)/FS)*1e3-(startpos)
    set(sliderRighth,'Value', round(toc*1e3)+startpos)
    drawnow;
    pause(.02);
end
set(sliderRighth,'Value', 1)
pause(1);

