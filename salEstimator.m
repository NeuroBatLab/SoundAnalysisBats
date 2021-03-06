function [sal,t] = salEstimator(SoundIn_filtered, FS, minFund, maxFund,RMSThresh)
DebugFig=0;
% Estimates the pitch
% saliency (sal) 
% soundIn is the sound pressure waveform.
% FS is the sampling rate

% Some user parameters (should be part of the function at some time)
if nargin<5
    RMSThresh=0.1;      % Minimum relative threshold on Max RMS to calculate pitch saliency
end
if nargin<3
    minFund = 600;       % Minimum fundamental frequency expected
end
if nargin<5
    maxFund = 4000;      % Maximum fundamental frequency expected
end

if DebugFig
    DBNOISE = 60;         % dB in Noise for the log compression in the spectrogram calculation - values below will be set to zero.
    f_high=10000;         % Upper frequency bound to get average amplitude in spectrogram
    fband = 100;            % Size of each frequency band in spectrogram, also determine time resolution
end

if DebugFig
    % Calculate and plot the spectrogram    
    figure(154)
    clf
    [ ~] = spec_only_bats(SoundIn_filtered, FS, DBNOISE, f_high, fband);
end


% Initializations and useful variables
Tt = 1; % period in ms, getting a value per ms
if size(SoundIn_filtered,2)==1
    SoundIn_filtered = SoundIn_filtered';
end
soundLen = length(SoundIn_filtered);
t = (0:Tt:(floor(soundLen/FS*10^3)-1))*10^-3; % getting a value per ms, note that t is in seconds
nt=length(t);
soundRMS = zeros(1,nt);
fund = zeros(1,nt);
fund2 = zeros(1,nt);
sal = zeros(1,nt);

% Calculate the size of the window for the auto-correlation
alpha = 5; % Number of sd in the Gaussian window
winLen = fix((2*alpha/minFund)*FS);    % Length of Gaussian window based on minFund
if (mod(winLen,2) == 0)  % Make a symmetric window
    winLen = winLen+1;
end
w = gausswin(winLen, alpha);
maxlags = 2*ceil((FS/minFund));


%% First calculate the rms in each window
for it = 1:nt
    tval = t(it)+Tt/2*10^-3;   % Center of window in time
    tind = fix(tval*FS);  % Center of window in ind
    tstart = tind - (winLen-1)/2;
    tend = tind + (winLen-1)/2;
    
    if tstart < 1
        winstart = 2 - tstart;
        tstart = 1;
    else
        winstart = 1;
    end
    
    if tend > soundLen
        windend = winLen - (tend-soundLen);
        tend = soundLen;
    else
        windend = winLen;
    end
    
    soundWin = SoundIn_filtered(tstart:tend).*w(winstart:windend)';
    soundRMS(it) = std(soundWin);
    
end

soundRMSMax = max(soundRMS);

%% Calculate the auto-correlation in windowed segments and the pitch saliency
soundlen = 0;
for it = 1:nt
    fund(it) = NaN;
    sal(it) = NaN;
    fund2(it) = NaN;
    if (soundRMS(it) < soundRMSMax*RMSThresh)
        continue;
    end
    soundlen = soundlen + 1;
    tval = t(it)+Tt/2*10^-3;   % Center of window in time
    tind = fix(tval*FS);  % Center of window in ind
    tstart = tind - (winLen-1)/2;
    tend = tind + (winLen-1)/2;
    
    if tstart < 1
        winstart = 2 - tstart;
        tstart = 1;
    else
        winstart = 1;
    end
    
    if tend > soundLen
        windend = winLen - (tend-soundLen);
        tend = soundLen;
    else
        windend = winLen;
    end
    
    soundWin = SoundIn_filtered(tstart:tend).*w(winstart:windend)';
    
    [autoCorr, lags] = xcorr(soundWin, maxlags, 'unbiased');
    Timelags = (lags/FS)*10^3;
    ind0 = find(lags == 0);
    
    % find peaks
    [pksCorr, indPeaksCorr] = findpeaks(autoCorr, 'MINPEAKPROMINENCE', max(autoCorr)./10);
    
    if DebugFig
        figure(155)
        clf
        plot(Timelags, autoCorr, 'LineWidth',2)
        hold on; plot(Timelags(indPeaksCorr),pksCorr,'r*')
    end
        
    % Eliminate center peak and all peaks too close to middle
    pksCorr(abs(indPeaksCorr-ind0) < FS/maxFund ) = [];
    indPeaksCorr(abs(indPeaksCorr-ind0) < FS/maxFund ) = [];
    if DebugFig
        hold on; plot(Timelags(indPeaksCorr),pksCorr,'g*')
    end
    
    % Find max peak
    if isempty(pksCorr)
        pitchSaliency = 0.1; 
    else
        indIndMax = find(pksCorr == max(pksCorr), 1, 'first');
        indMax = indPeaksCorr(indIndMax);   
        pitchSaliency = autoCorr(indMax)./autoCorr(ind0);
        if DebugFig
            hold on; plot(Timelags(indMax),max(pksCorr),'mo', 'MarkerSize',16)
            vline(10^3/maxFund)
            hold on
            vline(-10^3/maxFund)
            hold off
            title(sprintf('Autocorr for time point %d/%d Sal = %.2f',it,nt,pitchSaliency))
            xlabel('Time lag (ms)')
            ylabel('Autocorrelation')
            pause(1)
        end
    end

    sal(it) = pitchSaliency;
    
 
end

if DebugFig   
    figure(154)
    % plot the pitch saliency on top
    hold on
    yyaxis right
    plot(t*10^3,sal, 'r--', 'LineWidth',2)
    ylabel('Pitch Saliency')
    ylim([0 1])
    hold off
    pause(1)
end


end



