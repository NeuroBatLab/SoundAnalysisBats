function [sal,t] = salEstimator(SoundIn, FS)
% Estimates the pitch
% saliency (sal) 
% soundIn is the sound pressure waveform.
% FS is the sampling rate

% Some user parameters (should be part of the function at some time)
maxFund = 10000;      % Maximum fundamental frequency
minFund = 100;       % Minimum fundamental frequency
LowFc = 100;         % Low frequency cut-off for band-passing the signal prior to auto-correlation.
HighFc = 24000;       % High frequency cut-off
DBNOISE = 60;         % dB in Noise for the log compression in the spectrogram calculation - values below will be set to zero.
f_high=20000;         % Upper frequency bound to get average amplitude in spectrogram
fband = 100;            % Size of each frequency band in spectrogram, also determine time resolution

% Band-pass filtering signal prior to auto-correlation
[z,p,k] = butter(6,[LowFc HighFc]/(FS/2),'bandpass');
sos_band = zp2sos(z,p,k);
soundIn_filtered = filtfilt(sos_band,1,SoundIn);

% Calculate and plot the spectrogram    
figure(9)
[ ~] = spec_only_bats(soundIn_filtered, FS, DBNOISE, f_high, fband);


% Initializations and useful variables
soundLen = length(SoundIn);
t = (0:1:floor(length(SoundIn)/FS*10^3))*10^-3; % getting a value per ms, note that t is in seconds
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
    tval = t(it);   % Center of window in time
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
    
    soundWin = SoundIn(tstart:tend).*w(winstart:windend)';
    soundRMS(it) = std(soundWin);
    
end

soundRMSMax = max(soundRMS);

%% Calculate the auto-correlation in windowed segments and the pitch saliency
soundlen = 0;
for it = 1:nt
    fund(it) = NaN;
    sal(it) = NaN;
    fund2(it) = NaN;
    if (soundRMS(it) < soundRMSMax*0.1)
        continue;
    end
    soundlen = soundlen + 1;
    tval = t(it);   % Center of window in time
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
    
    soundWin = SoundIn(tstart:tend).*w(winstart:windend)';
    
    [autoCorr, lags] = xcorr(soundWin, maxlags, 'unbiased');
    ind0 = find(lags == 0);
    
    % find peaks
    [pksCorr, indPeaksCorr] = findpeaks(autoCorr, 'MINPEAKHEIGHT', max(autoCorr)./10);
    
    % Eliminate center peak and all peaks too close to middle
    
    pksCorr(abs(indPeaksCorr-ind0) < FS/maxFund ) = [];
    indPeaksCorr(abs(indPeaksCorr-ind0) < FS/maxFund ) = [];
    
    % Find max peak
    if isempty(pksCorr)
        pitchSaliency = 0.1; 
    else
        indIndMax = find(pksCorr == max(pksCorr), 1, 'first');
        indMax = indPeaksCorr(indIndMax);   
        pitchSaliency = autoCorr(indMax)./autoCorr(ind0);
    end

    sal(it) = pitchSaliency;
    
 
end

figure(9)
hold on
yyaxis left
plot(t*10^3, fund, 'k-', 'LineWidth',2)
hold on
yyaxis right
plot(t*10^3,sal, 'r--', 'LineWidth',2)
ylabel('Pitch Saliency')
ylim([0 1])
hold off
pause()


end



