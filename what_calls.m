function what_calls(Loggers_dir, Date, ExpStartTime)
% Hard coded parameters for the calculation of the spectrum in biosound
F_high_Raw = 50000;
F_high_Piezo = 10000;

% Set to 1 if you want to manually pause after each vocalization and listen
% to them
ManualPause=1;

% Import biosound library
py.importlib.import_module('soundsig')

% Load data
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

BioSoundFilenames = cell(VocCall,2);
BioSoundCalls = cell(VocCall,2);

%% Loop through calls, save them as wav files and run biosound
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
                NVocFile = NVocFile +1;
                fprintf(1,'%d/%d Vocalization\n',NVocFile,VocCall)
                % Extract the sound of the microphone that
                % correspond to the data
                IndOn = IndVocStartRaw_merged{VocInd(vv)}{ll}(nn);
                IndOff = IndVocStopRaw_merged{VocInd(vv)}{ll}(nn);
                WL = Raw_wave{VocInd(vv)}(IndOn:IndOff);
                BioSoundFilenames{NVocFile,1} = fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Raw.wav',FileVoc, BatID_local,ALNum,nn));
                audiowrite(BioSoundFilenames{NVocFile,1},WL,FS);
                BioSoundCalls{NVocFile,1} = runBiosound(WL, FS, F_high_Raw);
                % Plot figures of biosound results for Microphone data
                Fig1=figure(1);
                clf
                title(sprintf('%d/%d Vocalization',NVocFile,VocCall))
                plotBiosound(BioSoundCalls{NVocFile,1}, F_high_Raw)
                % Play the sound
                if ManualPause
                    AP=audioplayer(WL,FS);
                    play(AP)
                end
                print(Fig1,fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Raw.pdf', FileVoc, BatID_local,ALNum,nn)),'-dpdf','-fillpage')
                
                
                % Extract the sound of the audio-logger that
                % correspond to the data
                IndOn = IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn);
                IndOff = IndVocStopPiezo_merged{VocInd(vv)}{ll}(nn);
                WL = Piezo_wave.(Fns_AL{ll}){VocInd(vv)}(IndOn:min(IndOff, length(Piezo_wave.(Fns_AL{ll}){VocInd(vv)})));
                WL = WL - mean(WL); % center the piezo data around 0
                if any(abs(WL)>=1)
                    WL = WL./max(abs(WL)); % scale between 0 and 1 if exceeding 1
                end
                BioSoundFilenames{NVocFile,2} =fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Piezo.wav',FileVoc,BatID_local,ALNum,nn));
                FSpiezo = round(Piezo_FS.(Fns_AL{ll})(VocInd(vv)));
                audiowrite(BioSoundFilenames{NVocFile,2},WL,FSpiezo);
                BioSoundCalls{NVocFile,2} = runBiosound(WL, FSpiezo, F_high_Piezo);
                % Plot figures of biosound results for piezo data
                Fig2=figure(2);
                clf
                title(sprintf('%d/%d Vocalization',NVocFile,VocCall))
                plotBiosound(BioSoundCalls{NVocFile,2}, F_high_Piezo,0)
                % Play the sound
                if ManualPause
                    AP=audioplayer(WL,FSpiezo);
                    play(AP)
                end
                print(Fig2,fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Piezo.pdf', FileVoc, BatID_local,ALNum,nn)),'-dpdf','-fillpage')
                
                % Plot figures of dynamic jointly evaluated by piezo and
                % microphone data
                Fig3 = figure(3);
                clf
                title(sprintf('%d/%d Vocalization',NVocFile,VocCall))
                plotCallDynamic(BioSoundCalls{NVocFile,1}, BioSoundCalls{NVocFile,2})
                print(Fig3,fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Dyn.pdf', FileVoc, BatID_local,ALNum,nn)),'-dpdf','-fillpage')
                if ManualPause
                    pause()
                end
            end
        end
    end
end
% save the values!
save(fullfile(DataFile.folder, DataFile.name), 'BioSoundCalls','-append');

%% Internal functions


    function BiosoundObj = runBiosound(Y, FS, F_high)
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
        if nargin<3
            % Spectrum parameters
            F_high = 50000; % frequency of Low-pass filter Hz
        end
        % Fundamental parameters
        MaxFund = 4000;
        MinFund = 300;
        LowFc = 500;
        HighFc = 10000;
        MinSaliency = 0.6;
        DebugFigFundest = 0;
        MinFormantFreq = 2000;
        MaxFormantBW = 2000;
        Method= 'Stack';
        
        % create the biosound object
        BiosoundObj = py.soundsig.sound.BioSound(py.numpy.array(Y),pyargs('fs',FS));
        % methods(BiosoundFi, '-full') % this command plot all the methods with the available arguments
        
        % Calculate the RMS (lhs std(varargin))
        BiosoundObj.rms = BiosoundObj.sound.std();
        
        % calculate the amplitude enveloppe
        ampenv(BiosoundObj, Cutoff_freq,Amp_sample_rate);
        
        % calculate the spectrum (lhs spectrum(self, f_high, pyargs))
        spectrum(BiosoundObj, F_high)
        % calculate the spectrogram (lhs spectroCalc(self, spec_sample_rate,
        % freq_spacing, min_freq, max_freq, pyargs))
        try % For ver short sound, the Freq_spacing is too small, doubling if error
            spectroCalc(BiosoundObj, Spec_sample_rate, Freq_spacing, Min_freq,Max_freq)
        catch
            spectroCalc(BiosoundObj, Spec_sample_rate, Freq_spacing.*2, Min_freq,Max_freq)
        end
        
        % calculate the fundamental and related values (lhs fundest(self, maxFund,
        % minFund, lowFc, highFc, minSaliency, debugFig, pyargs)
        fundest(BiosoundObj, MaxFund, MinFund,LowFc, HighFc, MinSaliency,DebugFigFundest,MinFormantFreq,MaxFormantBW,Method)
    end

    function plotBiosound(BiosoundObj, F_high, FormantPlot)
        if nargin<3
            FormantPlot=1;
        end
        % Plot the results of biosound calculations
        subplot(2,1,1)
        ColorCode = get(groot,'DefaultAxesColorOrder');
        DBNOISE =12;
        f_low = 0;
        logB = - 20*log10(abs(double(BiosoundObj.spectro)));
        maxB = max(max(logB));
        minB = maxB-DBNOISE;
        
        imagesc(double(BiosoundObj.to)*1000,double(BiosoundObj.fo),logB);          % to is in seconds
        axis xy;
        caxis('manual');
        caxis([minB maxB]);
        cmap = spec_cmap();
        colormap(cmap);
%         colorbar()
        
        v_axis = axis;
        v_axis(3)=f_low;
        v_axis(4)=F_high;
        axis(v_axis);
        xlabel('time (ms)'), ylabel('Frequency');
        
        % Plot the fundamental and formants if they were calculated
        %     if double(BiosoundFi.sal)>MinSaliency
        Legend = {'F0' 'Formant1' 'Formant2' 'Formant3'};
        IndLegend = [];
        if ~isempty(double(BiosoundObj.f0))
            hold on
            plot(double(BiosoundObj.to)*1000,double(BiosoundObj.f0),'r-','LineWidth',2)
            IndLegend = [1 IndLegend];
        end
        if FormantPlot
            hold on
            plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F1),'Color',ColorCode(4,:),'LineWidth',2)
            hold on
            plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F2),'Color',ColorCode(2,:),'LineWidth',2)
            hold on
            if any(~isnan(double(BiosoundObj.F3)))
                plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F3),'Color',ColorCode(3,:),'LineWidth',2)
                IndLegend = [IndLegend 2:4];
            else
                IndLegend = [IndLegend 2:3];
            end
        end
        legend(Legend(IndLegend))
        hold off
        subplot(2,1,2)
        yyaxis left
        plot((1:length(double(BiosoundObj.sound)))/BiosoundObj.samprate*1000,double(BiosoundObj.sound), 'k-','LineWidth',2)
        hold on
        SoundAmp = double(py.array.array('d', py.numpy.nditer(BiosoundObj.amp)));
        yyaxis right
        plot(double(BiosoundObj.tAmp)*1000,double(SoundAmp), 'r-', 'LineWidth',2)
        YLIM = get(gca,'YLim');
        YLIM(1) = -YLIM(2);
        set(gca, 'YLim', YLIM)
        xlabel('Time (ms)')
        hold off
    end


    function plotCallDynamic(BiosoundRaw, BiosoundPiezo)
        Span = 9;% Span is an unevennumber. smooth has a default span of 5 points = 5ms However end points are unchanged...
        HalfSpan = (Span-1)/2;
        % Plot the pitch saliency vs amplitude on microphone
        subplot(3,1,1)
        Saliency = mysmooth(double(BiosoundRaw.sal), Span);
        SoundAmp = double(py.array.array('d', py.numpy.nditer(BiosoundRaw.amp)));
        TimeSound = double(BiosoundRaw.to)*1000;
        TimeSound = TimeSound./max(TimeSound);
        cmap = colormap('jet');
        ncolors = length(cmap);
        nx = length(Saliency);
        
        for ii=HalfSpan:nx-HalfSpan
            segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
            plot([Saliency(ii), Saliency(ii+1)], [SoundAmp(ii), SoundAmp(ii+1)], "Color",segcolor, "LineWidth",2);
            hold on;
        end
        set(gca,'XLim',[0 1]);
        xlabel(sprintf('Pitch Saliency %.1f', double(BiosoundRaw.meansal)))
        ylabel('Amplitude')
        
        % Plot the difference of formants (Mic data) vs sound amplitude (Mic
        % Data)
        subplot(3,1,2)
        SoundSpeed = 350;
        F1 = double(BiosoundRaw.F1);
        F2 = double(BiosoundRaw.F2);
        FormantDisp = mysmooth(SoundSpeed./(2*(F2 - F1))*1000, Span);
        nx = length(FormantDisp);
        
        
        for ii=HalfSpan:nx-HalfSpan
            segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
            plot([FormantDisp(ii), FormantDisp(ii+1)], [SoundAmp(ii), SoundAmp(ii+1)], "Color",segcolor, "LineWidth",2);
            hold on;
        end
        set(gca,'XLim',[10 130])
        xlabel('1/Formant disp (vocal tract length (mm))')
        ylabel('Amplitude')
        
        % Plot the ratio of formants (Mic data) vs fundamental (Piezo
        % Data)
        subplot(3,1,3)
        SoundFund = mysmooth(double(BiosoundPiezo.f0), Span);
        if ~isempty(SoundFund)
            for ii=HalfSpan:nx-HalfSpan
                segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
                plot([SoundFund(ii), SoundFund(ii+1)], [SoundAmp(ii), SoundAmp(ii+1)], "Color",segcolor, "LineWidth",2);
                hold on;
            end
            ylabel('Amplitude')
            xlabel(sprintf('Fundamental (Hz), %.1f Hz', double(BiosoundPiezo.fund)))
            set(gca,'XLim',[0 3000])
        end
        
        
    end

    function outyy = mysmooth(yy,Span)
        if nargin<2
            Span = 5;
        end
        outyy=nan(size(yy));
        for ii=1:length(yy)
            if ii==1 || ii==length(yy)
                outyy(ii) = yy(ii);
            elseif ii<((Span-1)/2)
                HalfSpan = ii-1;
                outyy(ii) = nanmean(yy(1:(ii+HalfSpan)));
            elseif (length(yy)-ii) < ((Span-1)/2)
                HalfSpan = length(yy)-ii;
                outyy(ii) = nanmean(yy((ii-HalfSpan):end));
            else
                outyy(ii) = nanmean(yy((ii-HalfSpan):(ii+HalfSpan)));
            end
        end
    end



end


