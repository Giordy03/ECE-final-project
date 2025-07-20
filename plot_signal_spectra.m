function plot_signal_spectra(f, spectra, title_text)
% PLOT_SIGNAL_SPECTRA Plots the power spectral density (PSD) of a signal.
%
% This function visualizes the spectral characteristics of a given signal
% using two different representations: 
%   - Normalized PSD in dB (logarithmic scale).
%   - Magnitude squared spectrum (linear scale).
%
% INPUTS:
%   f          - Frequency axis (Hz)
%   spectra    - Fourier Transform of the signal, representing its spectrum.
%   title_text - Title for the plot, describing the signal.
%
% OUTPUT:
%   A figure with two subplots:
%   - Left: Normalized PSD in dB.
%   - Right: Magnitude squared spectrum.


    % Compute normalization factor for dB conversion
    normalization = 10 * log10(max(abs(spectra).^2));

    % Create figure and set overall title
    figure
    sgtitle(title_text) 

    % First subplot: Normalized PSD in dB scale
    subplot(1, 2, 1)
    plot(f, 10 * log10(abs(spectra).^2) - normalization, "r")
    title("Normalized Power Spectral Density (PSD)")
    ylabel("|X(f)|^2 (dB)") 
    xlabel("Frequency (Hz)")
    grid on;

    % Second subplot: Magnitude squared spectrum (linear scale)
    subplot(1, 2, 2)
    plot(f, abs(spectra).^2, "r", "LineWidth", 1.2)
    title("Magnitude Squared Spectrum")
    ylabel("|X(f)|^2") 
    xlabel("Frequency (Hz)")
    grid on;
   
end
