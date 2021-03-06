function run_biosound_wav(InputDir,F_high_Raw,PlotFormant, SaveBiosoundperFile, SaveFig)

if nargin<5
    SaveFig=0;
end
if nargin<4
    SaveBiosoundperFile = 0;
end

if nargin<3
    PlotFormant=1;
end

if nargin<2
    % Hard coded parameters for the calculation of the spectrum in biosound
    F_high_Raw = 50000;
end

% Set to 1 if you want to manually pause after each vocalization and listen
% to them
ManualPause=0;

% Import biosound library
py.importlib.import_module('soundsig')

% Load data
% WaveFiles1 = dir(fullfile(InputDir, 'Kanamycin*.wav'));
% WaveFiles2 = dir(fullfile(InputDir, 'Saline-Deaf*.wav'));
% WaveFiles = [WaveFiles1; WaveFiles2];
WaveFiles = dir(fullfile(InputDir, '*.wav'));
if isempty(WaveFiles)
    warning('**** No vocalization in %s ****\n',InputDir)
else
    
    % Find out the sampling frequency from one file
    Info = audioinfo(fullfile(WaveFiles(1).folder, WaveFiles(1).name));
    FS = Info.SampleRate;
    
    % Filter for the Mic signal
    [z,p,k] = butter(3,100/(FS/2),'high');
    sos_high_raw = zp2sos(z,p,k);
    
    % create the output directoty
    OutputDir = fullfile(InputDir, 'BiosoundOutput');
    mkdir(OutputDir);
    
    % Number of vocalization cuts for preallocation of space
    NV = length(WaveFiles);
    BioSoundFilenames = cell(NV,1);
    BatID = cell(NV,1);
    Condition = cell(NV,1);
    if ~SaveBiosoundperFile
        BioSoundUniqParam = nan(NV,22);
    end
    BioSoundParamNames = {'stdtime' 'meantime' 'skewtime' 'entropytime'...
        'kurtosistime' 'AmpPeriodF' 'AmpPeriodP' 'rms' 'maxAmp' 'stdspect'...
        'meanspect' 'skewspect' 'entropyspect' 'kurtosisspect' 'q1' 'q2' 'q3'...
        'fund' 'cvfund' 'minfund' 'maxfund' 'meansal'};
        
    % Turn off warnings regarding Pyton to structure conversion
    warning('off', 'MATLAB:structOnObject')
    
    
    %% Loop through calls, save them as wav files and run biosound
    for vv=1:NV
        fprintf(1,'%d/%d Vocalization\n',vv,NV)
        % load the file
        [WL, FS] = audioread(fullfile(WaveFiles(vv).folder, WaveFiles(vv).name));
        % filter and center the data
        FiltWL = filtfilt(sos_high_raw,1,WL);
        FiltWL = FiltWL-mean(FiltWL);
        %parse out info
        BioSoundFilenames{vv} = fullfile(WaveFiles(vv).folder, WaveFiles(vv).name);
        %Ind_ = strfind(WaveFiles(vv).name,'_');
        %BatID{vv} = WaveFiles(vv).name((Ind_(1)-5):(Ind_(1)-1));
        %Condition{vv} = WaveFiles(vv).name(1:(Ind_(1)-9));
        BatID{vv} = WaveFiles(vv).name((end-8):(end-4));
        Condition{vv} = WaveFiles(vv).name(1:5);
        
        % run biosound
        BioSoundCall = runBiosound(FiltWL, FS, F_high_Raw);
        
        if ~SaveBiosoundperFile
            % Feed biosound data into a Matrix
            % temporal parameters (calculated on the envelope)
            BioSoundUniqParam(vv,1) = BioSoundCall.stdtime;
            BioSoundUniqParam(vv,2) = BioSoundCall.meantime;
            BioSoundUniqParam(vv,3) = BioSoundCall.skewtime;
            BioSoundUniqParam(vv,4) = BioSoundCall.entropytime;
            BioSoundUniqParam(vv,5) = BioSoundCall.kurtosistime;
            if ~isempty(BioSoundCall.AmpPeriodF)
                BioSoundUniqParam(vv,6) = BioSoundCall.AmpPeriodF;
                BioSoundUniqParam(vv,7) = BioSoundCall.AmpPeriodP;
            end

            % Amplitude parameters calculated on the envelope
             BioSoundUniqParam(vv,8) = BioSoundCall.rms;
              BioSoundUniqParam(vv,9) = BioSoundCall.maxAmp;

            % Spectral parameters calculated on the spectrum
            BioSoundUniqParam(vv,10) = BioSoundCall.stdspect;
            BioSoundUniqParam(vv,11) = BioSoundCall.meanspect;
            BioSoundUniqParam(vv,12) = BioSoundCall.skewspect;
            BioSoundUniqParam(vv,13) = BioSoundCall.entropyspect;
            BioSoundUniqParam(vv,14) = BioSoundCall.kurtosisspect;
            BioSoundUniqParam(vv,15) = BioSoundCall.q1;
            BioSoundUniqParam(vv,16) = BioSoundCall.q2;
            BioSoundUniqParam(vv,17) = BioSoundCall.q3;

            % Fundamental parameters
            if ~isempty(BioSoundCall.fund)
                BioSoundUniqParam(vv,18) = BioSoundCall.fund;
            end
            if ~isempty(BioSoundCall.cvfund)
                BioSoundUniqParam(vv,19) = BioSoundCall.cvfund;
            end
            if ~isempty(BioSoundCall.minfund)
                BioSoundUniqParam(vv,20) = BioSoundCall.minfund;
            end
            if ~isempty(BioSoundCall.maxfund)
                BioSoundUniqParam(vv,21) = BioSoundCall.maxfund;
            end
            if ~isempty(BioSoundCall.meansal)
                BioSoundUniqParam(vv,22) = BioSoundCall.meansal;
            end
        end
        
        % Plot figures of biosound results for Microphone data
        Fig1=figure(1);
        clf
        title(sprintf('%d/%d Vocalization',vv,NV))
        plotBiosound(BioSoundCall, F_high_Raw,PlotFormant)
        % Play the sound
        if ManualPause
            AP=audioplayer(FiltWL./(max(abs(FiltWL))),FS);
            play(AP)
        end
        if SaveFig
            print(Fig1,fullfile(OutputDir,sprintf('%s_biosound.pdf', WaveFiles(vv).name)),'-dpdf','-fillpage')
        end
        if SaveBiosoundperFile
            save(fullfile(OutputDir,sprintf('%s_biosound.mat', WaveFiles(vv).name)),'BioSoundCall')
        end
    end
    % Turn back on warnings regarding Pyton to structure conversion
    warning('on', 'MATLAB:structOnObject')

    if ~SaveBiosoundperFile
        % save all the values as a matrix!
        save(fullfile(OutputDir, 'BioSoundMatrix.mat'),'BioSoundFilenames','BioSoundUniqParam','BioSoundParamNames', 'BatID','Condition');
    end
end
end



%% Internal functions

    function BiosoundObj = runBiosound(Y, FS, F_high)
        % Hard coded parameters for biosound
        % Spectrum parameters
        if nargin<3
            F_high = 50000; % frequency of Low-pass filter Hz
        end
        Quartile_values = [0.25, 0.5, 0.75];
        DBNOISE = 60; % Threshold on amplitude of spectrogram to extract time varying spectral mean and quartiles
        
        % spectrogram parameters
        Spec_sample_rate = 1000; % sampling rate Hz
        Freq_spacing = 50; % width of the frequency window for the FFT Hz
        Min_freq = 300; % high pass filter before FFT Hz
        Max_freq = F_high; % Low pass filter before FFT Hz
        % temporal enveloppe parameters
        Cutoff_freq = 150; % Hz
        Amp_sample_rate = 1000; % Hz
        % Fundamental parameters
        MaxFund = 4000;
        MinFund = 300;
        LowFc = 100; %100
        HighFc = 18000;% 15000
        MinSaliency = 0.6;
        DebugFigFundest = 0;
        MinFormantFreq = 2000;
        MaxFormantBW = 1000; %500
        WindowFormant = 0.1;
        Method= 'Stack';
        
        % create the biosound object
        BiosoundObj = py.soundsig.sound.BioSound(py.numpy.array(Y),pyargs('fs',FS));
        % methods(BiosoundFi, '-full') % this command plot all the methods with the available arguments
        
        % Calculate the RMS (lhs std(varargin))
        BiosoundObj.rms = BiosoundObj.sound.std();
        
        % calculate the amplitude enveloppe
        ampenv(BiosoundObj, Cutoff_freq,Amp_sample_rate);
        
        % Calculate the periodicity of the amplitude envelope
        SoundAmp = double(py.array.array('d', py.numpy.nditer(BiosoundObj.amp)));
        [P,F] = pspectrum(SoundAmp,1000);
        [PKS,LOCS]=findpeaks(P);
        AmpPeriodF = F(LOCS(PKS == max(PKS))); % Frequency in hertz of the max peak
        AmpPeriodP = max(PKS)/mean(SoundAmp.^2); % Proportion of power in the max peak of the spectrum
        
        % calculate the spectrum (lhs spectrum(self, f_high, pyargs))
        spectrum(BiosoundObj, F_high)
        % calculate the spectrogram (lhs spectroCalc(self, spec_sample_rate,
        % freq_spacing, min_freq, max_freq, pyargs))
        try % For very short sound, the Freq_spacing is too small, doubling if error
            spectroCalc(BiosoundObj, Spec_sample_rate, Freq_spacing, Min_freq,Max_freq)
        catch
            spectroCalc(BiosoundObj, Spec_sample_rate, Freq_spacing.*2, Min_freq,Max_freq)
        end
        
        % Calculate time varying spectralmean and spectral quartiles
        Spectro = double(BiosoundObj.spectro);
        maxB = max(max(Spectro));
        minB = maxB-DBNOISE; 
        Spectro(Spectro<minB) =minB;
        Spectro = Spectro-minB;
        Fo = double(BiosoundObj.fo);
        TPoints = size(Spectro,2);
        SpectralMean = nan(1,TPoints);
        Quartile_freq = nan(3,TPoints);
        %         SpectralMax = nan(1,TPoints);
        for tt=1:TPoints
            %             SpectralMax(tt) = Fo(Spectro(:,tt)==max(Spectro(:,tt)));
            PSDSpec = Spectro(:,tt)./(sum(Spectro(:,tt)));
            SpectralMean(tt) = sum(PSDSpec' .* Fo);
            % Find quartile power
            Cum_power = cumsum(Spectro(:,tt));
            Tot_power = sum(Spectro(:,tt));
            Nfreqs = length(Cum_power);
            iq = 1;
            for ifreq=1:Nfreqs
                if (Cum_power(ifreq) > Quartile_values(iq)*Tot_power)
                    Quartile_freq(iq,tt) = Fo(ifreq);
                    iq = iq+1;
                    if (iq > 3)
                        break;
                    end
                end
            end
        end
        
        % calculate the fundamental and related values (lhs fundest(self, maxFund,
        % minFund, lowFc, highFc, minSaliency, debugFig, pyargs)
        fundest(BiosoundObj, MaxFund, MinFund,LowFc, HighFc, MinSaliency,DebugFigFundest,MinFormantFreq,MaxFormantBW,WindowFormant,Method)
        
        % convert biosound to a strcuture
        BiosoundObj = struct(BiosoundObj);
        % Add some fields
        BiosoundObj.AmpPeriodF = AmpPeriodF;
        BiosoundObj.AmpPeriodP = AmpPeriodP;
        BiosoundObj.SpectralMean = SpectralMean;
        %         BiosoundObj.SpectralMax = SpectralMax;
        BiosoundObj.Quartile_freq = Quartile_freq;
        % convert all nmpy arrays to double to be able to save as matfiles
        BiosoundObj.amp = SoundAmp;
        BiosoundObj.tAmp = double(BiosoundObj.tAmp);
        BiosoundObj.spectro = double(BiosoundObj.spectro);
        BiosoundObj.to = double(BiosoundObj.to);
        BiosoundObj.fo = double(BiosoundObj.fo);
        BiosoundObj.F1 = double(BiosoundObj.F1);
        BiosoundObj.F2 = double(BiosoundObj.F2);
        BiosoundObj.F3 = double(BiosoundObj.F3);
        BiosoundObj.fpsd = double(BiosoundObj.fpsd);
        BiosoundObj.psd = double(BiosoundObj.psd);
        BiosoundObj.sal = double(BiosoundObj.sal);
        BiosoundObj.f0 = double(BiosoundObj.f0);
        BiosoundObj.f0_2 = double(BiosoundObj.f0_2);
        BiosoundObj.fund = double(BiosoundObj.fund);
        BiosoundObj.cvfund = double(BiosoundObj.cvfund);
        BiosoundObj.fund2 = double(BiosoundObj.fund2);
        BiosoundObj.minfund = double(BiosoundObj.minfund);
        BiosoundObj.maxfund = double(BiosoundObj.maxfund);
        BiosoundObj.sound = double(BiosoundObj.sound);
        BiosoundObj.wf = double(BiosoundObj.wf);
        BiosoundObj.wt = double(BiosoundObj.wt);
        BiosoundObj.mps = double(BiosoundObj.mps);
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
        YLIM = get(gca,'YLim');
        YLIM = max(abs(YLIM)).*[-1 1];
        set(gca, 'YLim', YLIM)
        SoundAmp = double(py.array.array('d', py.numpy.nditer(BiosoundObj.amp)));
        yyaxis right
        plot(double(BiosoundObj.tAmp)*1000,double(SoundAmp), 'r-', 'LineWidth',2)
        YLIM = get(gca,'YLim');
        YLIM = max(abs(YLIM)).*[-1 1];
        set(gca, 'YLim', YLIM)
        set(gca, 'XLim', v_axis(1:2))
        xlabel('Time (ms)')
        title(sprintf('AmpPeriodicity = %.3f AmpPF = %.1f Hz',BiosoundObj.AmpPeriodP, BiosoundObj.AmpPeriodF))
        hold off
    end


    function plotCallDynamic(BiosoundRaw, BiosoundPiezo)
        Span = 9;% Span is an unevennumber. smooth has a default span of 5 points = 5ms However end points are unchanged...
        HalfSpan = (Span-1)/2;
        % Plot the pitch saliency vs amplitude on microphone
        subplot(4,1,1)
        Saliency = mysmooth(double(BiosoundRaw.sal), Span);
        TimeSound = double(BiosoundRaw.to)*1000;
        TimeSound = TimeSound./max(TimeSound);
        cmap = colormap('jet');
        ncolors = length(cmap);
        nx = length(Saliency);
        
        for ii=HalfSpan:nx-HalfSpan
            segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
            plot([Saliency(ii), Saliency(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
            hold on;
        end
        set(gca,'XLim',[0 1]);
        xlabel(sprintf('Pitch Saliency %.1f', double(BiosoundRaw.meansal)))
        ylabel('Amplitude')
        
        % Plot the difference of formants (Mic data) vs sound amplitude (Mic
        % Data)
        subplot(4,1,2)
        SoundSpeed = 350;
        F1 = double(BiosoundRaw.F1);
        F2 = double(BiosoundRaw.F2);
        FormantDisp = mysmooth(SoundSpeed./(2*(F2 - F1))*1000, Span);
        nx = length(FormantDisp);
        
        
        for ii=HalfSpan:nx-HalfSpan
            segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
            plot([FormantDisp(ii), FormantDisp(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
            hold on;
        end
        set(gca,'XLim',[10 150])
        xlabel('1/Formant disp (vocal tract length (mm))')
        ylabel('Amplitude')
        
        % Plot the amplitude (Mic data) vs fundamental (Piezo
        % Data)
        subplot(4,1,3)
        SoundFund = mysmooth(double(BiosoundPiezo.f0), Span);
        if ~isempty(SoundFund)
            for ii=HalfSpan:nx-HalfSpan
                segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
                plot([SoundFund(ii), SoundFund(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
                hold on;
            end
            ylabel('Amplitude')
            xlabel(sprintf('Fundamental (Hz), %.1f Hz', double(BiosoundPiezo.fund)))
            set(gca,'XLim',[200 3000])
        end
        
        % Plot the Amplitude (Mic data) vs SpectralMean (Mic
        % Data)
        subplot(4,1,4)
        SoundSpecMean = mysmooth(double(BiosoundRaw.SpectralMean), Span);
        if ~isempty(SoundSpecMean)
            for ii=HalfSpan:nx-HalfSpan
                segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
                plot([SoundSpecMean(ii), SoundSpecMean(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
                hold on;
            end
            ylabel('Amplitude')
            xlabel(sprintf('Spectral Mean (Hz), %.1f Hz', nanmean(double(BiosoundRaw.SpectralMean))))
            set(gca,'XLim',[25000 30000])
        end
        
        %         % Plot the Amplitude (Mic data) vs Spectral Max (Piezo
        %         % Data)
        %         subplot(5,1,5)
        %         SoundSpecMax = mysmooth(double(BiosoundPiezo.SpectralMax), Span);
        %         if ~isempty(SoundSpecMax)
        %             for ii=HalfSpan:nx-HalfSpan
        %                 segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
        %                 plot([SoundSpecMax(ii), SoundSpecMax(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
        %                 hold on;
        %             end
        %             ylabel('Amplitude')
        %             xlabel(sprintf('Spectral Max (Hz), %.1f Hz', nanmean(double(BiosoundPiezo.SpectralMax))))
        %             set(gca,'XLim',[0 10000])
        %         end
        
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





