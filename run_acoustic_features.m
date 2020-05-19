function [AcousticFeatureValues,AcounsticFeatureNames, Spectro] = run_acoustic_features(Sound, FS,  varargin)
%% RUN_ACOUSTIC_FEATURES calculate 16 bioacoustic parameters on the sound input Sound
% the pitch saliency
% Q1, Q2, Q3, spectral mean
pnames = {'F_high', 'F_low', 'F_highSpec', 'Spectro'};
dflts  = {5000, 100, 15000, 0};
[F_high, F_low, F_highSpec,Spectro] = internal.stats.parseArgs(pnames,dflts,varargin{:});

%    F_high = 5000 is good value for logger data to focus on sound emitted by the bat
%    F_low = 100;  remove noise
%    F_highSpec = 15000; % Good vaue to catch noise artefact (high power in high frequencies)


% Bandpass the input signal in the range of freq of interest prior applying calculations
[z,p,k] = butter(6,[F_low F_high]/(FS/2),'bandpass');
sos_band = zp2sos(z,p,k);
Sound_filtered = filtfilt(sos_band,1,Sound);

% Calculate pitch saliency
[Sal,~] = salEstimator(Sound_filtered, FS);
MeanSal = nanmean(Sal);

% Temporal enveloppe features
Cutoff_freq = 150; % Hz
Amp_sample_rate = 1000; % Hz
[MaxAmp, RMS, MeanTime, StdTime, KurtosisTime, SkewTime, EntropyTime] = temp_env_paramEstimator(Sound_filtered, FS, Cutoff_freq,Amp_sample_rate);

% Spectral enveloppe features
[z,p,k] = butter(6,[F_low F_highSpec]/(FS/2),'bandpass');
sos_band = zp2sos(z,p,k);
Sound_filtered = filtfilt(sos_band,1,Sound);
[QuartileSpec, MeanSpec, StdSpec, KurtosisSpec, SkewSpec, EntropySpec] = spec_env_paramEstimator(Sound_filtered, FS, F_high);

% Spectrogram
FBand = 50;
DBNOISE = 50;
if Spectro
    [to, fo, logB, ~] = spec_only_bats(Sound_filtered, FS, DBNOISE, F_highSpec, FBand);
    Spectro = struct();
    Spectro.to = to;
    Spectro.fo = fo;
    Spectro.logB = logB;
end
% Return values in a vector
AcousticFeatureValues = [MeanSal, MaxAmp, RMS, MeanTime, StdTime, KurtosisTime, SkewTime, EntropyTime, QuartileSpec(1), QuartileSpec(2), QuartileSpec(3), MeanSpec, StdSpec, KurtosisSpec, SkewSpec, EntropySpec];
AcounsticFeatureNames = {'MeanSaliency' 'MaxAmp' 'RMS' 'MeanTime' 'StdTime' 'KurtosisTime' 'SkewTime' 'EntropyTime' 'Q1' 'Q2' 'Q3' 'MeanSpec' 'StdSpec' 'KurtosisSpec' 'SkewSpec' 'EntropySpec'};
end

