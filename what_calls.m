function what_calls(Loggers_dir, Date, ExpStartTime)

% Hard coded parameters for biosound
% spectrogram parameters
Spec_sample_rate = 1000; % sampling rate Hz
Freq_spacing = 50; % width of the frequency window for the FFT Hz
Min_freq = 300; % high pass filter before FFT Hz
Max_freq = 50000; % Low pass filter before FFT Hz
% temporal enveloppe parameters (currenty unused, I don't find a way of
% passing them to the function under Matlab
Cutoff_freq = 75; % Hz
Amp_sample_rate = 1000; % Hz
% Spectrum parameters
F_high = 50000; % frequency of Low-pass filter Hz
% Fundamental parameters
MaxFund = 4000;
MinFund = 300;
LowFc = 500;
HighFc = 10000;
MinSaliency = 0.5;
DebugFigFundest = 0;
MinFormantFreq = 500;
MaxFormantBW = 500;
Method= 'Stack';



DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date, ExpStartTime)));
load(fullfile(DataFile.folder, DataFile.name), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'BatID','LoggerName');
load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'FS','Piezo_wave','Raw_wave', 'Piezo_FS','VocFilename');
% Number of call sequences with identified vocalizations
VocInd = find(~cellfun('isempty',IndVocStartRaw_merged));
NV = length(VocInd);
Fns_AL = fieldnames(Piezo_wave);
% create the output directoty
Path2Wav = fullfile(Loggers_dir, 'VocExtracts');
mkdir(Path2Wav);

% Count the number of vocalization cuts for preallocation of space
VocCall = 0;
for vv=1:NV
    for ll=1:length(IndVocStartRaw_merged{VocInd(vv)})
        VocCall = VocCall + length(IndVocStartRaw_merged{VocInd(vv)}{ll});
    end
end

BioSoundFilenames = cell(VocCall,1);

%% Loop through calls and save them as wav files
Ncall = nan(NV,1);

NVocFile = 0;
for vv=1:NV
    [~,FileVoc]=fileparts(VocFilename{VocInd(vv)});
    for ll=1:length(IndVocStartRaw_merged{VocInd(vv)})
        % Logger number
        AL_local = Fns_AL{ll};
        ALNum = AL_local(7:end);
        % ID of the bat
        ALIndex = contains(LoggerName, 'AL') .* contains(LoggerName, ALNum);
        BatID_local =BatID{find(ALIndex)};
        Ncall(vv) = length(IndVocStartRaw_merged{VocInd(vv)}{ll});
        if Ncall(vv)
            for nn=1:Ncall(vv)
                % Extract the sound of the microphone that
                % correspond to the data
                IndOn = IndVocStartRaw_merged{VocInd(vv)}{ll}(nn);
                IndOff = IndVocStopRaw_merged{VocInd(vv)}{ll}(nn);
                WL = Raw_wave{VocInd(vv)}(IndOn:IndOff);
                NVocFile = NVocFile +1;
                BioSoundFilenames{NVocFile} = fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Raw.wav',FileVoc, BatID_local,ALNum,nn));
                audiowrite(BioSoundFilenames{NVocFile},WL,FS);
                
                
                % Extract the sound of the audio-logger that
                % correspond to the data
                IndOn = IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn);
                IndOff = IndVocStopPiezo_merged{VocInd(vv)}{ll}(nn);
                WL = Piezo_wave.(Fns_AL{ll}){VocInd(vv)}(IndOn:min(IndOff, length(Piezo_wave.(Fns_AL{ll}){VocInd(vv)})));
                WL = WL - mean(WL); % center the piezo data around 0
                if any(abs(WL)>=1)
                    WL = WL./max(abs(WL)); % scale between 0 and 1 if exceeding 1
                end
                NVocFile = NVocFile +1;
                BioSoundFilenames{NVocFile} =fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Piezo.wav',FileVoc,BatID_local,ALNum,nn));
                audiowrite(BioSoundFilenames{NVocFile},WL,round(Piezo_FS.(Fns_AL{ll})(VocInd(vv))));
            end
        end
    end
end

%% Run the biosound on all the extracts
py.importlib.import_module('soundsig')

for vv=1:NVocFile
    [Y] =audioread(BioSoundFilenames{vv});
    BiosoundFi = py.soundsig.sound.BioSound(py.numpy.array(Y),pyargs('fs',50000));
    % methods(BiosoundFi, '-full') % this command plot all the methods with the available arguments
    
    % Calculate the RMS (lhs std(varargin))
    BiosoundFi.rms = BiosoundFi.sound.std();
    
    % calculate the amplitude enveloppe
    ampenv(BiosoundFi, Cutoff_freq,Amp_sample_rate);
    
    % calculate the spectrum (lhs spectrum(self, f_high, pyargs))
    spectrum(BiosoundFi, F_high)
    
    % calculate the spectrogram (lhs spectroCalc(self, spec_sample_rate,
    % freq_spacing, min_freq, max_freq, pyargs))
    spectroCalc(BiosoundFi, Spec_sample_rate, Freq_spacing, Min_freq,Max_freq)
    
    % calculate the fundamental and related values (lhs fundest(self, maxFund,
    % minFund, lowFc, highFc, minSaliency, debugFig, pyargs)
    fundest(BiosoundFi, MaxFund, MinFund,LowFc, HighFc, MinSaliency,DebugFigFundest,MinFormantFreq,MaxFormantBW,Method)
    
    % export/save the values!
    
end



end
