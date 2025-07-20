function song_modulated = zero_padding(song_modulated, Nch)
    % zero_padding: Adjusts the length of the modulated song signal to match Nch.
    %
    % Inputs:
    %   song_modulated - Input signal to be adjusted (column vector)
    %   Nch            - Desired length of the output signal
    %
    % Output:
    %   song_modulated - Padded or truncated signal with length Nch
    if length(song_modulated) > Nch
        song_modulated = song_modulated(1:Nch);
    elseif length(song_modulated) < Nch
        song_modulated = [song_modulated; zeros((Nch) - length(song_modulated), 1)];
    end
end
