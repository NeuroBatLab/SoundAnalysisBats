% This script was meant to check again call-Type manually annotated and
% correct some misclassifications
% we are using the sounds, spectros and call-type collected by DeafBats_PythonPlot_dataGather.mlx
% and comparing with the call type collected by DeafBats_CatCalls1.mlx
% (just to be on the safe side)
% Note that only MicAudioGood calls have been collected by
% DeafBats_PythonPlot_dataGather.mlx so we are only checking again these
% high quality vocalizations
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
LocalDataDir = '/Users/elie/Documents/DeafBats/Data';
load(fullfile(LocalDataDir, 'MicData4_DeafBats_PythonPlot.mat'), 'BatID', 'ALID', 'CallType','CallDate','MicFileNum','MicFileSampleID',...
    'OnSetOffsetSample', 'Sound_Mic', ...
    'Spectro_mic', 'Spectro_mic_to','Spectro_mic_fo', 'SampleRate_Mic');
load(fullfile(LocalDataDir, 'PiezoData4_DeafBats_PythonPlot.mat'),'Spectro', 'Spectro_to', 'Spectro_fo',...
    'Sound','SampleRate');
D = load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls.mat'), 'BatID', 'ALID', 'CallType','CallDate','MicFileNum','MicFileSampleID',...
    'OnSetOffsetSample','MicAudioGood');

VolFactorMic = 0.5;

MicAudioGood01 = logical(D.MicAudioGood);
D.BatID = D.BatID(MicAudioGood01);
D.ALID = D.ALID(MicAudioGood01);
D.CallType = D.CallType(MicAudioGood01);
D.CallDate = D.CallDate(MicAudioGood01);
D.MicFileNum = D.MicFileNum(MicAudioGood01);
D.MicFileSampleID = D.MicFileSampleID(MicAudioGood01);
D.OnSetOffsetSample = D.OnSetOffsetSample(find(MicAudioGood01),:);
DBNOISE = 60;
FHigh = 50000;
FHigh_piezo = 10000;
SampRate = 192000;
SampRate_piezo = 50000;

NVoc = sum(MicAudioGood01);
% Check that vocalizations are in the same order
if sum(str2double(D.BatID) == BatID)~=NVoc || sum(strcmp(D.ALID, ALID))~=NVoc...
        || sum(strcmp(D.CallType, CallType))~=NVoc || sum(D.CallDate == CallDate)~=NVoc...
        || sum(D.MicFileNum == MicFileNum)~=NVoc || sum(D.MicFileSampleID == MicFileSampleID)~=NVoc...
        || sum(D.OnSetOffsetSample(:,1) == OnSetOffsetSample(:,1))~=NVoc...
        || sum(D.OnSetOffsetSample(:,2) == OnSetOffsetSample(:,2))~=NVoc
    error('Issue of allignment')
end
CallType_Old = CallType; % keep for record the old cold type

% Prepare file output
if ~exist(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_CallType2ndCheckMicAudioGood.mat'), 'file')
    save(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_CallType2ndCheckMicAudioGood.mat'),'BatID', 'ALID', 'CallType','CallDate','MicFileNum','MicFileSampleID',...
    'OnSetOffsetSample')
    vvstart=1;
    CallType = cell(size(CallType));
else
    load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_CallType2ndCheckMicAudioGood.mat'), 'vv', 'CallType')
    if ~exist('vv', 'var')
        vv=1;
    end
    vvstart = vv;
end

FIG=figure(1);
for vv=vvstart:NVoc
    clf(FIG)
   fprintf(1, '---------------------------------------\n');
    fprintf(1, 'Call %d/%d\n', vv, NVoc);
   % Plot the spectrograms of the extract
   subplot(1,2,1)
   logB = Spectro_mic{vv};
   maxB = max(max(logB));
   minB = maxB-DBNOISE;
   IM=imagesc(Spectro_mic_to{vv}*1000,Spectro_mic_fo{vv},logB);          % to is in seconds
   axis xy;
   caxis('manual');
   caxis([minB maxB]);
   cmap = spec_cmap();
   colormap(IM.Parent, cmap);
   v_axis = axis;
   v_axis(3)=0;
   v_axis(4)=FHigh;
   axis(v_axis);
   xlabel('Time ms')
   ylabel('Frequency')
   title(sprintf('Spectro Mic Call %d/%d CallType = %s\n', vv, NVoc, CallType_Old{vv}))
   
   subplot(1,2,2)
   logB = Spectro{vv};
   maxB = max(max(logB));
   minB = maxB-DBNOISE;
   IM=imagesc(Spectro_to{vv}*1000,Spectro_fo{vv},logB);          % to is in seconds
   axis xy;
   caxis('manual');
   caxis([minB maxB]);
   cmap = spec_cmap();
   colormap(IM.Parent, cmap);
   v_axis = axis;
   v_axis(3)=0;
   v_axis(4)=FHigh_piezo;
   axis(v_axis);
   title(sprintf('Spectro Piezo Call %d/%d CallType = %s\n', vv, NVoc, CallType_Old{vv}))
   xlabel('Time ms')
   ylabel('Frequency')
   suplabel(sprintf('Call %d/%d\n', vv, NVoc), 't');
   
   % prepare playback of the extract
   APM = audioplayer(Sound_Mic{vv}, SampRate/4,24);
            %             APM=audioplayer(BioSoundCalls{vv,1}.sound./(max(abs(BioSoundCalls{vv,1}.sound))),BioSoundCalls{vv,1}.samprate);
   APP=audioplayer(Sound{vv}./(max(abs(Sound{vv}))),SampRate_piezo);
   
   % Correct the call type if needed
   INPUT=[];
   while isempty(INPUT)
       play(APM)
       pause(max(1, length(Sound_Mic{vv})/SampRate/4))
       play(APP)
       commandwindow
       INPUT = input(sprintf('*** Type of call was %s ***\n  Correct (11)\n Bad Mic Quality (99)\n  Trill (1), Bark(2), pitchy call(3), low buzz(4) panting (5)\n    Low tuck (6) Squeal (7) Rattle (8) Chuckles (9) LoudPicthy (10)\n    Cricket (12) Unknown (0)\n BarkBuzz (24) PitchyCallBuzz (34) SquealBuzz (74)', CallType_Old{vv}));
       if isempty(INPUT) || ((INPUT~=0) &&(INPUT~=1) && (INPUT~=2) && (INPUT~=3)&& (INPUT~=4) && (INPUT~=5)&& (INPUT~=6)&& (INPUT~=7)&& (INPUT~=8)&& (INPUT~=9)&& (INPUT~=10) && (INPUT~=11) && (INPUT~=12) && (INPUT~=24) && (INPUT~=34) && (INPUT~=74) && (INPUT~=99))
           INPUT=[];
       end
   end
   if INPUT==1
       CallType{vv} = 'Tr';
   elseif INPUT==2
       CallType{vv} = 'Ba';
   elseif INPUT==3
       CallType{vv} = 'Pi';
   elseif INPUT==4
       CallType{vv} = 'Bu';
   elseif INPUT==5
       CallType{vv} = 'Pa';
   elseif INPUT==6
       CallType{vv} = 'LT';
   elseif INPUT==7
       CallType{vv} = 'Sq';
   elseif INPUT==8
       CallType{vv} = 'Ra';
   elseif INPUT==9
       CallType{vv} = 'Ch';
   elseif INPUT==10
       CallType{vv} = 'LPi';
   elseif INPUT==24
       CallType{vv} = 'BB';
   elseif INPUT==34
       CallType{vv} = 'PB';
   elseif INPUT==74
       CallType{vv} = 'SB';
   elseif INPUT==0
       CallType{vv} = 'Un';
   elseif INPUT==11
       CallType{vv} = CallType_Old{vv};
   elseif INPUT==12
       CallType{vv} = 'Cr';
   elseif INPUT==99
       CallType{vv} = 'OO';
   end
   
  save(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_CallType2ndCheckMicAudioGood.mat'),  'CallType','vv', '-append')
end