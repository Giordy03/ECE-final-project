function plot_SIR(sir, freqs, maxSIR, optimalFreq)
% PLOT_SIR Plots the Signal-to-Interference Ratio (SIR) vs. cutoff frequency
% 
% Inputs:
%   sir        - A vector containing SIR values (in dB) for each cutoff frequency
%   freqs      - A vector containing the corresponding cutoff frequencies (in Hz)
%   maxSIR     - The maximum SIR value (in dB) achieved
%   optimalFreq - The cutoff frequency (in Hz) at which the maximum SIR occurs
% 
% Output:
%   A plot of SIR vs. cutoff frequency with the optimal point highlighted

% Create a new figure
figure
hold on
grid on
% Plot SIR values against cutoff frequencies
plot(freqs, sir, '-o')
xlabel('Cutoff Frequency (Hz)')
ylabel('SIR (dB)')
title('SIR vs. Cutoff Frequency')

% Highlight the optimal point (cutoff frequency with maximum SIR)
plot(optimalFreq, maxSIR, 'ro', 'MarkerSize', 10, 'LineWidth', 2)
text(optimalFreq, maxSIR, sprintf('Max SIR\n%.2f Hz, %.2f dB', optimalFreq, maxSIR), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center')

hold off
end