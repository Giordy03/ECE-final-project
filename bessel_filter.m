function [x_out, cut_off, maxSIR, sir] = bessel_filter(x_in, x_ref, cut_off_range, filterOrder, Fs)
    cut_off = 0;
    maxSIR = 0;
    sir = zeros(1, length(cut_off_range));
    i = 0;
    x_out = -inf;

    for cutoffFrequency = cut_off_range
        i = i + 1;      % update iteration index
        
        % Design Bessel filter (the function besself output the numerator B and
        % denumerator A for an analogue filter
        [B, A] = besself(filterOrder, 2 * pi * cutoffFrequency);  
        B = double(B);
        A = double(A);
        [Bz, Az] = bilinear(B, A, Fs);  % Convert to digital filter
    
        % Apply filter to signal
        xFiltered = filter(Bz, Az, x_in);
        
        % consider delay caused by filtering
        delay = finddelay(x_ref, xFiltered);
        xFiltered = circshift(xFiltered, -delay);
    
        % Calculate SIR
        sir(i) = determine_SIR(2*xFiltered, x_ref);
    
        % Check for maximum SIR
        if sir(i) > maxSIR
            maxSIR = sir(i);
            cut_off = cutoffFrequency;
            x_out = 2*xFiltered;
            A_bessel_best = Az; 
            B_bessel_best = Bz;
        end
    end
    % A = A_bessel_best;
    % B = B_bessel_best;
end