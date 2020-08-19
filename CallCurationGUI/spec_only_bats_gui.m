function [to, fo, logB, pg, tError, fError] = spec_only_bats_gui(sound_in, samprate, DBNOISE, f_high, fband, varargin)
% The Theunissen Lab Spectrogram with Guassian window.  Plots the
% spectrogram and oscillogram of the sound

if nargin<3
    DBNOISE = 60;                          % dB in Noise for the log compression - values below will be set to zero.
end

if nargin<4
    f_high=60000;                               % Upper frequency bound to get average amplitude in spectrogram
end

if nargin<5
    fband = 100;                            % Size of each frequency band
end

pnames = {'Time_increment'};
dflts  = {0.0001};
[Time_increment] = internal.stats.parseArgs(pnames,dflts,varargin{:});

% Parameters for the Spectrogram
nstd = 6;
twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
winLength = fix(winLength/2)*2;            % Enforce even window length
increment = fix(Time_increment*samprate);           % Sampling rate of spectrogram in number of points - set at approximately 1000 Hz
f_low=0;                                 % Lower frequency bounds to get average amplitude in spectrogram

% Calculate and plot the spectrogram    
[s, to, fo, pg, tError, fError] = GaussianSpectrum(sound_in, increment, winLength, samprate); 
logB = 20*log10(abs(s));
maxB = max(max(logB));
minB = maxB-DBNOISE;            

imagesc(axh,to*1000,fo,logB);          % to is in seconds
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