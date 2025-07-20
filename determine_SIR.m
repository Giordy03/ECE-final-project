function [SIR] = determine_SIR(x_out, x_ref)
% DETERMINE_SIR Calculates the Signal-to-Interference Ratio (SIR) in decibels (dB)
%
% Inputs:
%   x_out  - The output signal (e.g., filtered or processed signal)
%   x_ref  - The reference signal (e.g., original or desired signal)
%
% Output:
%   SIR    - The Signal-to-Interference Ratio (SIR) in decibels (dB)

% Calculate the power of the output signal
signalEnergy = sum(abs(x_out).^2); 

% Calculate the power of the interference (difference between x_out and x_ref)
interferenceEnergy = sum(abs(x_out - x_ref).^2);  

% Calculate the Signal-to-Interference Ratio (SIR) in decibels (dB)
SIR = 10 * log10(signalEnergy / interferenceEnergy);
% SIR = abs(SIR);
end