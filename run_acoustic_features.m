function [AcousticFeatureValues,AcounsticFeatureNames] = run_acoustic_features(Sound, FS, F_high, F_low, F_highSpec)
%% RUN_ACOUSTIC_FEATURES calculate 16 bioacoustic parameters on the sound input Sound
% the pitch saliency
% Q1, Q2, Q3, spectral mean


if nargin<3
    F_high = 5000; % good value for logger data to focus on sound emitted by the bat
end
if nargin<4
    F_low = 100; % remove noise
end
if nargin<5
    F_highSpec = 15000; % Good vaue to catch noise artefact (high power in high frequencies)
end

% Bandpass the input signal in the range of freq of interest prior applying calculations
[z,p,k] = butter(6,[F_low F_high]/(FS/2),'bandpass');
sos_band = zp2sos(z,p,k);
Sound_filtered = filtfilt(sos_band,1,Sound);

% Calculate pitch saliency
[Sal,Sal_t] = salEstimator(Sound_filtered, FS, F_low, F_high);

% Temporal enveloppe features
Cutoff_freq = 150; % Hz
Amp_sample_rate = 1000; % Hz
[MaxAmp, RMS, MeanTime, StdTime, KurtosisTime, SkewTime, EntropyTime] = temp_env_paramEstimator(Sound_filtered, FS, Cutoff_freq,Amp_sample_rate);

% Spectral enveloppe features
[z,p,k] = butter(6,[F_low F_highSpec]/(FS/2),'bandpass');
sos_band = zp2sos(z,p,k);
Sound_filtered = filtfilt(sos_band,1,Sound);
[QuartileSpec, MeanSpec, StdSpec, KurtosisSpec, SkewSpec, EntropySpec] = spec_env_paramEstimator(Sound_filtered, FS, F_high);

        
% Return values in a vector
AcousticFeatureValues = [Sal, Sal_t, MaxAmp, RMS, MeanTime, StdTime, KurtosisTime, SkewTime, EntropyTime, QuartileSpec(1), QuartileSpec(2), QuartileSpec(3), MeanSpec, StdSpec, KurtosisSpec, SkewSpec, EntropySpec];
AcounsticFeatureNames = {'Saliency' 'SaliencyT' 'MaxAmp' 'RMS' 'MeanTime' 'StdTime' 'KurtosisTime' 'SkewTime' 'EntropyTime' 'Q1' 'Q2' 'Q3' 'MeanSpec' 'StdSpec' 'KurtosisSpec' 'SkewSpec' 'EntropySpec'};
end

