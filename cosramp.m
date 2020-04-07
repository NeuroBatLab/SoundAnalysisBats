function Signal_out = cosramp(Signal_in, Ramp_nsamples)

% Apply half-raised-cosine ramp to a signal begining and end

Ramp_nsamples = round(Ramp_nsamples);
Cos_mask = hanning(2*Ramp_nsamples-1);
Mask = ones(size(Signal_in));
Mask(1)=0;
Mask(2:Ramp_nsamples) = Cos_mask(1:Ramp_nsamples-1);
Mask(end) = 0;
Mask(end-Ramp_nsamples+1:end-1) = Cos_mask(Ramp_nsamples+1:end);

Signal_out = Signal_in .* Mask;
end