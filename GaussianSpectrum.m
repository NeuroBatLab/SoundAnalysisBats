function [s, to, fo, pg, tError, fError] = GaussianSpectrum(input, increment, winLength, samprate)
%
% Gaussian spectrum
% 	s = GuassianSpectrum(input, increment, winLength, samprate)
% 	Compute the guassian spectrogram of an input signal with a guassian
% 	window that is given by winLength. The standard deviation of that
% 	guassian is 1/6th of winLength.
%	Each time frame is [winLength]-long and
%	starts [increment] samples after previous frame's start.
%	Only zero and the positive frequencies are returned.
%   to and fo are the time and frequency for each bin in s and Hz
%   pg is a rumming rms.
%   Etf is the mean square error between marginals and the input^2 and
%   frequency ^2

%%%%%%%%%%%%%%%%%%%%%%%
% Massage the input
%%%%%%%%%%%%%%%%%%%%%%%
debugfig = 0;

% Enforce even winLength to have a symmetric window
if rem(winLength, 2) == 1
    winLength = winLength +1;
end

% Make input it into a row vector if it isn't
if size(input, 1) > 1,
	input = input';
end;

% Padd the input with zeros
pinput = zeros(1,length(input)+winLength);
pinput(winLength/2+1:winLength/2+length(input)) = input;
inputLength = length(pinput);

% The number of time points in the spectrogram
frameCount = floor((inputLength-winLength)/increment)+1;

% The window of the fft
fftLen = winLength;

% Temporal power
x2 = pinput.*pinput;
ti = (-winLength/2:inputLength-1-winLength/2)./samprate;

% Fourier power
fval = fft(pinput)/sqrt(inputLength);   % The matlab fft is not symmetric
f2 = real(fval.*conj(fval));   
if rem(inputLength, 2)   % winLength is odd
    select = 1:(inputLength+1)/2;
else
    select = 1:inputLength/2+1;
end
fi = (select-1)'*samprate/inputLength;


%%%%%%%%%%%%%%%%%%%%%%%%
% Guassian window 
%%%%%%%%%%%%%%%%%%%%%%%%
nstd = 6;                   % Number of standard deviations in one window.
wx2 = ((1:winLength)-((winLength+1)/2)).^2;
wvar = (winLength/nstd)^2;
ws = exp(-0.5*(wx2./wvar));  
sumws = sum(ws);
ws = ws.*(increment/sumws);      % nomalize so that area is one * increment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize output "s" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rem(fftLen, 2)   
    % winLength is odd
    s = zeros((fftLen+1)/2+1, frameCount);
else
    % winLength is even 
    s = zeros(fftLen/2+1, frameCount);
end

pg = zeros(1, frameCount);
for i=1:frameCount
    start = (i-1)*increment + 1;
    last = start + winLength - 1;
    f = zeros(fftLen, 1);
    f(1:winLength) = ws.*pinput(start:last);
    pg(i) = std(f(1:winLength));

    specslice = fft(f);  % Here we don't need normalization because we did it with the window
    if rem(fftLen, 2)   % winLength is odd
        s(:,i) = specslice(1:((fftLen+1)/2+1));
    else
        s(:,i) = specslice(1:(fftLen/2+1));
    end
    %s(:,i) = specslice(1:(fftLen/2+1));
end

% Assign frequency_label
if rem(fftLen, 2)   % winLength is odd
    select = 1:(fftLen+1)/2;
else
    select = 1:fftLen/2+1;
end
fo = (select-1)'*samprate/fftLen;

% assign time_label
to = ((1:size(s,2))-1)'.*(increment/samprate);    % time starts at zero because we padded the signal with 1/2 window of zeros on both sides

% Calculate the marginals
absS = abs(s).^2;
tMarginal = timeseries(sum(absS,1)'.*2, to);      % Here both are multiplied by 2 because they both come from the "upper" part of the spectrogram
fMarginal = timeseries(sum(absS,2).*2, fo);
% Let's smooth tactual to take into account min res of our spectrogram.
% x2pow = sum(x2);
% if (rem(increment, 2) )
%     x2 = smooth(x2, increment);
% else
%     x2 = smooth(x2, increment-1);
% end
% x2 = x2.*(x2pow/sum(x2));
x2 = x2';
tActual = timeseries(x2, ti);
fActual = timeseries(f2(1:length(fi))'.*2, fi);   % Here we multiply by 2 because we are only taking on half of the frequency spectrum

% Match the sampling rate of all these vectors but preserve power.
tMarginalPower = sum(tMarginal.data);
tActualPower = sum(tActual.data);
fMarginalPower = sum(fMarginal.data);
fActualPower = sum(fActual.data);

[tMarginal, tActual] = synchronize(tMarginal, tActual, 'Uniform', 'Interval', ti(2)-ti(1));
[fMarginal, fActual] = synchronize(fMarginal, fActual, 'Uniform', 'Interval', fi(2)-fi(1));

% The overall power should be matched here but it is not...
tMarginal.data = tMarginal.data.*(tMarginalPower/sum(tMarginal.data));
tActual.data = tActual.data.*(tActualPower/sum(tActual.data));
fMarginal.data = fMarginal.data.*(fMarginalPower/sum(fMarginal.data));
fActual.data = fActual.data.*(fActualPower/sum(fActual.data));

% So we match it here...
tMarginal.data = tMarginal.data.*(tActualPower/sum(tMarginal.data));
fMarginal.data = fMarginal.data.*(fActualPower/sum(fMarginal.data));

tError = sum(abs(tActual.data-tMarginal.data))./(length(input)/samprate);   % The error is normalized by duration of signal - if we want to compare different stims...
fError = sum(abs(fActual.data-fMarginal.data))./(length(input)/samprate);

% Plot marginals on new figure for debugging purpose.
if debugfig
    % Checking Parseval's theorem
    fprintf(1, 'Parseval check t^2 = %f f^2 =%f\n', sum(tActual.data), sum(fActual.data));
    fprintf(1, 'Temporal marginal = %f Spectral marginal = %f\n', sum(tMarginal.data), sum(fMarginal.data)); 

    cfh = gcf;  % Get a handle to current figure
    figure();
    subplot(211);   
    plot(tActual, 'k');
    hold on;
    plot(tMarginal, 'r', 'LineWidth', 2);
    xlabel('Time (s)');
    title(sprintf('Temporal Marginal Error %f', tError));
    hold off;
    subplot(212);  
    plot(fActual, 'k');
    hold on;
    plot(fMarginal, 'r', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    title(sprintf('Spectral Marginal Error %f', fError));
    hold off;
    figure(cfh);  % Restore current figure
end



return

