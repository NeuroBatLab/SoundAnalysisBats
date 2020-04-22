function [maxAmp, RMS, MeanTime, StdTime, KurtosisTime, SkewTime, EntropyTime] = temp_env_paramEstimator(Filtered_SoundIn, FS,Cutoff_freq, FS_env)

if nargin<4
    FS_env=1000;% sampling frequency of the envelope in Hz
end
% Calculating the RMS
RMS = std(Filtered_SoundIn(Filtered_SoundIn~=0));

% Calculate the temporal envelope

[Amp_env, ~]=running_rms(Filtered_SoundIn, FS, Cutoff_freq, FS_env);

% Max Amplitude
maxAmp = max(Amp_env);

if length(Amp_env)~=fix(length(Filtered_SoundIn)*FS_env/FS)
    if abs(length(Amp_env)-fix(length(Filtered_SoundIn)*FS_env/FS))>1
        warning('Issue in envelope calculation, the envelope is different by more than one sample from the expected length\n')
        keyboard
    end
end

% Calculate the momentum of the temporal envelope
tAmp_env = 0:round(length(Filtered_SoundIn)*FS_env/FS);
Amp_env = Amp_env./sum(Amp_env);
tAmp_env = tAmp_env(1:length(Amp_env));
MeanTime = sum(tAmp_env.*Amp_env);
StdTime = sqrt(sum(Amp_env.*((tAmp_env-MeanTime).^2)));
SkewTime = sum(Amp_env.*(tAmp_env-MeanTime).^3);
SkewTime = SkewTime./(StdTime.^3);
KurtosisTime = sum(Amp_env.*(tAmp_env-MeanTime).^4);
KurtosisTime = KurtosisTime./(StdTime.^4);
indpos = find(Amp_env>0);
EntropyTime = -sum(Amp_env(indpos).*log2(Amp_env(indpos)))/log2(length(indpos));
end
