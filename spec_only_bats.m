function [to, fo, logB, pg, tError, fError] = spec_only_bats(sound_in, fband, samprate, dBScale, f_high)
% The Theunissen Lab Spectrogram with Guassian window.  Plots the
% spectrogram and oscillogram of the sound

% Parameters for the Spectrogram
nstd = 6;
twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
winLength = fix(winLength/2)*2;            % Enforce even window length
increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at approximately 1000 Hz
f_low=0;                                 % Lower frequency bounds to get average amplitude in spectrogram
%f_high=15000;                               % Upper frequency bound to get average amplitude in spectrogram
DBNOISE = dBScale;                          % dB in Noise for the log compression - values below will be set to zero.

% Calculate and plot the spectrogram    
[s, to, fo, pg, tError, fError] = GaussianSpectrum(sound_in, increment, winLength, samprate); 
logB = 20*log10(abs(s));
maxB = max(max(logB));
minB = maxB-DBNOISE;            

imagesc(to*1000,fo,logB);          % to is in seconds
axis xy;
caxis('manual');
caxis([minB maxB]); 
cmap = spec_cmap();
colormap(cmap);

v_axis = axis; 
v_axis(3)=f_low; 
v_axis(4)=f_high;
axis(v_axis);                                

xlabel('time (ms)'), ylabel('Frequency');

end