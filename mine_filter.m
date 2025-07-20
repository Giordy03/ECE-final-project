function [x_out, f_pole, maxSIR, sir, A_best1, B_best1] = mine_filter(x_in, x_ref, pole_freq_range, Fs, B, r)
    % mine_filter: Optimizes a second-order IIR filter by selecting the best pole frequency
    % to maximize the Signal-to-Interference Ratio (SIR).
    %
    % Inputs:
    %   x_in            - Input noisy signal
    %   x_ref           - Reference clean signal
    %   pole_freq_range - Range of pole frequencies to test
    %   Fs              - Sampling frequency
    %   B               - Numerator coefficients of the filter
    %   r               - Pole radius
    %
    % Outputs:
    %   x_out   - Best filtered output signal
    %   f_pole  - Best pole frequency that maximizes SIR
    %   maxSIR  - Maximum achieved SIR value
    %   sir     - Array of SIR values for each tested frequency
    %   A_best1 - Denominator coefficients of the optimal filter
    %   B_best1 - Numerator coefficients of the optimal filter

    maxSIR = 0;
    f_pole = 0;
    sir = zeros(size(pole_freq_range));
    i = 0;
    for poleFreq = pole_freq_range
        i = i+1;
        % Design poles
        p1 = r * exp(1j*2*pi*poleFreq/Fs);         
        p2 = conj(p1);                  
        A = poly([p1, p2]);             
        
        DC_gain = sum(B) / sum(A);      % Gain at low freq = 0 dB 
        B = B / DC_gain;     % Scale numerator
        xFiltered1 = filter(B, A, x_in);
        
        delay = finddelay(x_ref, xFiltered1);    
        xFiltered1 = circshift(xFiltered1, -delay);
        sir(i) = determine_SIR(2*xFiltered1, x_ref);
        
        if sir(i) > maxSIR
            maxSIR = sir(i);
            f_pole = poleFreq;
            x_out = 2*xFiltered1;
            A_best1 = A;
            B_best1 = B;
        end
    end

