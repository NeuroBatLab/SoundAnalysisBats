function what_calls(Loggers_dir, Date, ExpStartTime)
% Hard coded parameters for the calculation of the spectrum in biosound
F_high_Raw = 50000;
F_high_Piezo = 10000;

% Set to 1 if you want to manually pause after each vocalization and listen
% to them
ManualPause=0;

% Import biosound library
py.importlib.import_module('soundsig')

% Load data
Data1 = fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime));
if ~isfile(Data1)
    warning('No vocalization data extracted by who_calls.m or get_logger_data_voc.m')
else
    DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date, ExpStartTime)));
    if length(DataFile)>1
        warning('2 potential data files where found')
        DataFile
        warning('we work with the first one %s',DataFile(1).name)
    end
    % bringing the file back on the local computer (we're going to write
    % pretty often to it)
    WorkDir = ['~' filesep 'WorkingDirectoryWhat'];
    fprintf(1,'Transferring data from the server %s\n on the local computer %s\n', DataFile(1).folder, WorkDir);
    mkdir(WorkDir)
    [s,m,e]=copyfile(fullfile(DataFile(1).folder, DataFile(1).name), WorkDir, 'f');
    if ~s
        m %#ok<NOPRT>
        e %#ok<NOPRT>
        error('File transfer did not occur correctly for %s\n', fullfile(DataFile(1).folder, DataFile(1).name));
    end
    load(fullfile(WorkDir, DataFile(1).name), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'BatID','LoggerName');
    load(Data1, 'FS','Piezo_wave','Raw_wave', 'Piezo_FS','VocFilename');
    
    
    
    % Number of call sequences with identified vocalizations
    VocInd = find(~cellfun('isempty',IndVocStartRaw_merged));
    NV = length(VocInd);
    Fns_AL = fieldnames(Piezo_wave);
    
    % Filter for the Mic signal
    [z,p,k] = butter(3,100/(FS/2),'high');
    sos_high_raw = zp2sos(z,p,k);
    
    % Filter for the Piezo signal
    PFS = round(Piezo_FS.(Fns_AL{1})(1));
    [z,p,k] = butter(3,100/(PFS/2),'high');
    sos_high_piezo = zp2sos(z,p,k);
    
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
    
    try
        load(fullfile(WorkDir, DataFile(1).name), 'BioSoundCalls','BioSoundFilenames','NVocFile','vv');
        if exist('BioSoundCalls','var')
            PrevData = input('Do you want to use previous data?');
        else
            fprintf(1, 'No previous data, starting from scratch');
            PrevData = 0;
        end
    catch
        fprintf(1, 'No previous data, starting from scratch');
        PrevData = 0;
    end
    if ~PrevData
        BioSoundFilenames = cell(VocCall,2);
        BioSoundCalls = cell(VocCall,2);
        Firstcall = 1;
        NVocFile = 0;
    else
        Firstcall=vv;
        NVocFile = NVocFile-1;
    end
    
    %% Loop through calls, save them as wav files and run biosound
    Ncall = nan(NV,1);
    
    % Turn off warning notifications for python 2 struct conversion
    warning('off', 'MATLAB:structOnObject')
    
    
    for vv=Firstcall:NV
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
                    IndOff = min(length(Raw_wave{VocInd(vv)}),IndVocStopRaw_merged{VocInd(vv)}{ll}(nn)); % we take the min here as sometimes the rounding procedures gets numbers outisde of wave length
                    if IndOn>=IndOff
                        keyboard
                    end
                    WL = Raw_wave{VocInd(vv)}(IndOn:IndOff);
                    FiltWL = filtfilt(sos_high_raw,1,WL);
                    FiltWL = FiltWL-mean(FiltWL);
                    BioSoundFilenames{NVocFile,1} = fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Raw.wav',FileVoc, BatID_local,ALNum,nn));
                    audiowrite(BioSoundFilenames{NVocFile,1},FiltWL,FS);
                    BioSoundCalls{NVocFile,1} = runBiosound(FiltWL, FS, F_high_Raw);
                    % Plot figures of biosound results for Microphone data
                    Fig1=figure(1);
                    clf
                    title(sprintf('%d/%d Vocalization',NVocFile,VocCall))
                    plotBiosound(BioSoundCalls{NVocFile,1}, F_high_Raw)
                    % Play the sound
                    if ManualPause
                        AP=audioplayer(FiltWL./(max(abs(FiltWL))),FS);
                        play(AP)
                    end
                    print(Fig1,fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Raw.pdf', FileVoc, BatID_local,ALNum,nn)),'-dpdf','-fillpage')
                    
                    
                    % Extract the sound of the audio-logger that
                    % correspond to the data
                    IndOn = IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn);
                    IndOff = IndVocStopPiezo_merged{VocInd(vv)}{ll}(nn);
                    if IndOn>=IndOff
                        keyboard
                    end
                    WL = Piezo_wave.(Fns_AL{ll}){VocInd(vv)}(IndOn:min(IndOff, length(Piezo_wave.(Fns_AL{ll}){VocInd(vv)})));
                    WL = WL - mean(WL); % center the piezo data around 0
                    if any(abs(WL)>=1)
                        WL = WL./max(abs(WL)); % scale between 0 and 1 if exceeding 1
                    end
                    FiltWL = filtfilt(sos_high_piezo,1,WL);
                    BioSoundFilenames{NVocFile,2} =fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Piezo.wav',FileVoc,BatID_local,ALNum,nn));
                    FSpiezo = round(Piezo_FS.(Fns_AL{ll})(VocInd(vv)));
                    audiowrite(BioSoundFilenames{NVocFile,2},FiltWL,FSpiezo);
                    BioSoundCalls{NVocFile,2} = runBiosound(FiltWL, FSpiezo, F_high_Piezo);
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
%                     print(Fig3,fullfile(Path2Wav,sprintf('%s_Bat%d_AL%s_Elmt%d_Dyn.pdf', FileVoc, BatID_local,ALNum,nn)),'-dpdf','-fillpage')
                    
                    
                    % Guess for the call category
                    try double(BioSoundCalls{NVocFile,1}.AmpPeriodP)
                        if (BioSoundCalls{NVocFile,1}.AmpPeriodF<40.5) && (BioSoundCalls{NVocFile,1}.AmpPeriodF>34) && (BioSoundCalls{NVocFile,1}.AmpPeriodP>0.075)
                            Guess ='Tr';
                        else
                            Guess ='Ba';
                        end
                        
                    catch
                        Guess = 'Ba';
                    end
                    %                 if ManualPause
                    %                     Resp = input(sprintf('Is this a Trill (t) or a Bark (b)? Computer guess: %s. Leave empty if you agree',Guess),'s');
                    Resp = [];
                    if isempty(Resp)
                        BioSoundCalls{NVocFile,1}.type = Guess;
                        BioSoundCalls{NVocFile,2}.type = Guess;
                    elseif strcmp(Resp, 't')
                        BioSoundCalls{NVocFile,1}.type = 'Tr';
                        BioSoundCalls{NVocFile,2}.type = 'Tr';
                    elseif strcmp(Resp, 'b')
                        BioSoundCalls{NVocFile,1}.type = 'Ba';
                        BioSoundCalls{NVocFile,2}.type = 'Ba';
                    end
                    %                 end
                end
            end
        end
        
        % save the values!
        if length(DataFile)>1
            save(fullfile(WorkDir, DataFile(1).name), 'BioSoundCalls','BioSoundFilenames','NVocFile','-append');
        else
            save(fullfile(WorkDir, DataFile.name), 'BioSoundCalls','BioSoundFilenames','NVocFile','-append');
        end
    end
    % Turn off warning notifications for python 2 struct conversion
    warning('on', 'MATLAB:structOnObject')
    % Transfer data bacn on the server
    fprintf(1,'Transferring data from the local computer %s\n back on the server %s\n', WorkDir, DataFile(1).folder);
    [s,m,e]=copyfile(fullfile(WorkDir, DataFile(1).name), DataFile(1).folder, 'f');
    if ~s
        TicTransfer = tic;
        while toc(TicTransfer)<30*60
            [s,m,e]=copyfile(fullfile(WorkDir, DataFile(1).name), DataFile(1).folder, 'f');
            if s
                return
            end
        end
        if ~s
            s %#ok<NOPRT>
            m %#ok<NOPRT>
            e %#ok<NOPRT>
            error('File transfer did not occur correctly for %s\n Although we tried for 30min\n', DataFile(1).folder);
        else
            fprintf('Data transfered back on server in:\n%s\n',  DataFile(1).folder);
        end
    else
        fprintf('Data transfered back on server in:\n%s\n',  DataFile(1).folder);
    end
    if s  %erase local data
        [sdel,mdel,edel]=rmdir(WorkDir, 's');
        if ~sdel
            TicErase = tic;
            while toc(TicErase)<30*60
                [sdel,mdel,edel]=rmdir(WorkDir, 's');
                if sdel
                    return
                end
            end
        end
        if ~sdel
            sdel %#ok<NOPRT>
            mdel %#ok<NOPRT>
            edel %#ok<NOPRT>
            error('File erase did not occur correctly for %s\n Although we tried for 30min\n', WorkDir);
        end
    end
end

%% Internal functions


    function BiosoundObj = runBiosound(Y, FS, F_high)
        % Hard coded parameters for biosound
        % spectrogram parameters
        Spec_sample_rate = 1000; % sampling rate Hz
        Freq_spacing = 50; % width of the frequency window for the FFT Hz
        Min_freq = 300; % high pass filter before FFT Hz
        Max_freq = 50000; % Low pass filter before FFT Hz
        % temporal enveloppe parameters
        Cutoff_freq = 150; % Hz
        Amp_sample_rate = 1000; % Hz
        if nargin<3
            % Spectrum parameters
            F_high = 50000; % frequency of Low-pass filter Hz
        end
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
        
        % Calculate time varying spectralmean and spectral max
        Spectro = double(BiosoundObj.spectro);
        Fo = double(BiosoundObj.fo);
        TPoints = size(Spectro,2);
        SpectralMean = nan(1,TPoints);
        %         SpectralMax = nan(1,TPoints);
        for tt=1:TPoints
            %             SpectralMax(tt) = Fo(Spectro(:,tt)==max(Spectro(:,tt)));
            PSDSpec = Spectro(:,tt)./(sum(Spectro(:,tt)));
            SpectralMean(tt) = sum(PSDSpec' .* Fo);
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
        BiosoundObj.sound = double(BiosoundObj.sound);
        BiosoundObj.wf = double(BiosoundObj.wf);
        BiosoundObj.wt = double(BiosoundObj.wt);
        BiosoundObj.mps = double(BiosoundObj.mps);
        BiosoundObj = rmfield(BiosoundObj,'emitter');
        BiosoundObj.hashid = double(BiosoundObj.hashid);
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



end


