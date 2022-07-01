%This function masks a 3D wave spectrum, setting values to NaN's outside
%the wavenumber limits and frequency limits specified by the inputs.



function Spectrum_out = mask_wave_spectrum(Spectrum,wavenumberLimits,frequencyLimits)


Spectrum_out = Spectrum;
P = Spectrum.power_Spectrum;

K = sqrt(Spectrum.Kx_3D.^2 + Spectrum.Ky_3D.^2);

kmask = and(K>=min(abs(wavenumberLimits)),K<=max(abs(wavenumberLimits)));
omega_mask = and(Spectrum.W_3D>=min(abs(frequencyLimits)),Spectrum.W_3D<=max(abs(frequencyLimits)));
maskTot = kmask.*omega_mask;

P(~maskTot) = NaN;
Spectrum_out.power_Spectrum = P;


end