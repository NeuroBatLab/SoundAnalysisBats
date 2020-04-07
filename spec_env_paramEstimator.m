function [Quartile_freq, MeanSpect, StdSpect, KurtosisSpect, SkewSpect, EntropySpect] = spec_env_paramEstimator(Sound_filtered, FS, F_high)

Noverlap = 256;
NFFT = 512;% the spectral resolution is then F_high/NFFT
% Calculating the spectrum
[Psd,Freqs] = pwelch(Sound_filtered,[],Noverlap,NFFT,FS);

% Find quartile power
    Cum_power = cumsum(Psd);
    Tot_power = sum(Psd);
    Quartile_freq = nan(1,3);
    Quartile_values = [0.25, 0.5, 0.75];
    Nfreqs = length(Cum_power);
    iq = 1;
    for ifreq=1:Nfreqs
        if (Cum_power(ifreq) > Quartile_values(iq)*Tot_power)
            Quartile_freq(iq) = Freqs(ifreq);
            iq = iq+1;
            if (iq > 3)
                break;
            end
        end
    end
    
    % Calculate the momentum of the spectral envelope
% Find skewness, kurtosis and entropy for power spectrum below
    % f_high
    Ind_fmax = Nfreqs;
    for ifreq=1:Nfreqs
        if (Freqs(ifreq) > F_high )
            Ind_fmax = ifreq;
            break;
        end
    end
    
    % Description of spectral enveloppe
    Spectdata = Psd(1:Ind_fmax);
    Freqdata = Freqs(1:Ind_fmax);
    Spectdata = Spectdata./sum(Spectdata);
    MeanSpect = sum(Freqdata.*Spectdata);
    StdSpect = sqrt(sum(Spectdata.*((Freqdata-MeanSpect).^2)));
    SkewSpect = sum(Spectdata.*(Freqdata-MeanSpect).^3);
    SkewSpect = SkewSpect./(StdSpect.^3);
    KurtosisSpect = sum(Spectdata.*(Freqdata-MeanSpect).^4);
    KurtosisSpect = KurtosisSpect./(StdSpect.^4);
    EntropySpect = -sum(Spectdata.*log2(Spectdata))/log2(Ind_fmax);
end
    