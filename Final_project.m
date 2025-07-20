clc, clear, close all;
%% Task 0
SONG = 1;
    % SONG == 1: IMAGINE-John Lennon
    % SONG == 2: MAMMA MIA-Abba 

% % Load reference signal (withount any noise and/or distortion)
switch SONG
    case 1
        [data, Fs] = audioread("John Lennon - Imagine.mp3");
        title_song = "Imagine";
    case 2
        [data, Fs] = audioread("ABBA - Mamma Mia.mp3");
        title_song = "Mamma mia";
end

% data is a stereo signal; make it a mono audio (consider only one of the
% two column
song_ref = data(:, 1);    
Nsamples = length(song_ref);

% definition of axis
Ts = 1 / Fs;
t = 0:Ts:(Nsamples-1)*Ts;
t = t';         % we need column vector
DeltaF = Fs / Nsamples;
f = -Fs/2:DeltaF:Fs/2-DeltaF;
f = f';

% compute and plot Fourier transform of the reference signal (emplying fft)
Song_ref = fftshift(fft(song_ref))*Ts;
plot_signal_spectra(f, Song_ref, "Reference signal: "+title_song)

% Load signal with noise
switch SONG
    case 1
        song_noise= load("SONG1_disturbance.mat");
    case 2
        song_noise = load("SONG2_disturbance.mat");
end
% loaded data is a structure, take the relevant data
song_noise = song_noise.xtot;       
Song_noise = fftshift(fft(song_noise))*Ts;
plot_signal_spectra(f, Song_noise, "Noisy signal: "+title_song)

SIR_initial = determine_SIR(Song_noise, Song_ref);
fprintf("The SIR for the original song is %.2f \n", SIR_initial)

% compare reference and noisy signal in the time domain
figure
sgtitle("Time domain signal: "+title_song)
subplot(2, 1, 1)
plot(t, song_ref, "r")
xlabel("time (s)")
ylabel("amplitude")
title("Reference signal")
subplot(2, 1, 2)
plot(t, song_noise, "r")
xlabel("time (s)")
ylabel("amplitude")
title("Noisy signal")

% to find the frequency of the disturbance use the command 
% [maxValue, index] = max(abs(Song_noise)), then the frequency is f(index)
[noise, ind_noise] = max(abs(Song_noise));
f_noise = abs(f(ind_noise));
% the absolute value is needed to ensure that the f_noise is positive (the
% negative frequency has the same value, but for calculation purposes we
% need the positive one)
% Example disturbance frequency in Hz
fprintf("The frequency of the sinusoids is: %.2f kHz \n", f_noise)
%% TASK 1, BESSEL filter
filterOrder = 4;  % 4th-order filter
cutoffFrequencies_bessel = 3000:10:4000;  % Range of cutoff frequencies

% Initialize variables to store results
optimalCutoff_bessel = 0;
maxSIR_bessel = 0;
sir_bessel = zeros(1, length(cutoffFrequencies_bessel));
i = 0;

% Loop over cutoff frequencies
for cutoffFrequency = cutoffFrequencies_bessel
    i = i + 1;      % update iteration index

    % Design Bessel filter (the function besself output the numerator B and
    % denumerator A for an analogue filter
    [B, A] = besself(filterOrder, 2 * pi * cutoffFrequency);  
    B = double(B);
    A = double(A);
    [Bz, Az] = bilinear(B, A, Fs);  % Convert to digital filter

    % Apply filter to signal
    xFiltered = filter(Bz, Az, song_noise);
    
    % consider delay caused by filtering
    delay = finddelay(song_ref, xFiltered);
    xFiltered = circshift(xFiltered, -delay);

    % Calculate SIR
    sir_bessel(i) = determine_SIR(xFiltered, song_ref);

    % Check for maximum SIR
    if sir_bessel(i) > maxSIR_bessel
        maxSIR_bessel = sir_bessel(i);
        optimalCutoff_bessel = cutoffFrequency;
        song_filtered_bessel = xFiltered;
        A_bessel_best = Az; 
        B_bessel_best = Bz;
    end
end

% Visualize Frequency Response% Compute the frequency response
[H, F] = freqz(B_bessel_best, A_bessel_best, f, Fs);

% Plot the magnitude response
figure
subplot(2, 1, 1);
hold on
grid on
plot(F, 20*log10(abs(H)), "LineWidth",  1.2) % Convert magnitude to dB
xline(optimalCutoff_bessel, "r", "LineWidth", 1.2)
xline(-optimalCutoff_bessel, "r", "LineWidth", 1.2)
title("Magnitude Response")
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")

% Plot the phase response
subplot(2, 1, 2);

hold on
grid on
plot(F, unwrap(angle(H)), "LineWidth", 1.2) % Unwrap phase for better visualization
xline(optimalCutoff_bessel, "r", "LineWidth", 1.2)
xline(-optimalCutoff_bessel, "r", "LineWidth", 1.2)
title("Phase Response")
xlabel("Frequency (Hz)")
ylabel("Phase (radians)")

% find filter transfer function (Z domain)
z = tf("z", Ts);
H_bessel = zpk(tf(B_bessel_best, A_bessel_best, Ts))

figure
pzplot(H_bessel)

% Display optimal cutoff frequency and maximum SIR
fprintf("Optimal Cutoff Frequency: %.2f Hz\n", optimalCutoff_bessel);
fprintf("Maximum SIR: %.2f dB\n", maxSIR_bessel);

% Display the filtered signal in the frequency domain
Song_filtered_bessel = fftshift(fft(song_filtered_bessel))*Ts;
plot_signal_spectra(f, Song_filtered_bessel, "Signal filtered using " + ...
    "4^{th} order Bessel filter: "+title_song)

plot_SIR(sir_bessel, cutoffFrequencies_bessel, maxSIR_bessel, ...
    optimalCutoff_bessel)

%% Task1, BUTTERWORTH (optional)
cutoffFrequencies_butter = 2450:10:3000;
optimalCutoff_butter = 0;
maxSIR_butter = 0;
sir_butter = zeros(1, length(cutoffFrequencies_butter));
i = 0;

for cutoffFrequency = cutoffFrequencies_butter
    i = i + 1;

    % Design Butterworth filter (by default LPF, digital)
    % cut off frequency must be normalized (between 0, 1)
    [B_butter, A_butter] = butter(filterOrder, cutoffFrequency/(Fs/2));  
    
    % Apply filter to signal
    xFiltered_butter = filter(B_butter, A_butter, song_noise);

    delay = finddelay(song_ref, xFiltered_butter);
    xFiltered_butter = circshift(xFiltered_butter, -delay);

    sir_butter(i) = determine_SIR(xFiltered_butter, song_ref);

    if sir_butter(i) > maxSIR_butter
        maxSIR_butter = sir_butter(i);
        optimalCutoff_butter = cutoffFrequency;
        song_filtered_butter = xFiltered_butter;
        A_butter_best = A_butter; 
        B_butter_best = B_butter;
    end
end

figure
freqz(B_butter_best, A_butter_best, f, Fs);

H_butter = zpk(tf(B_butter_best, A_butter_best, Ts))

fprintf("Optimal Cutoff Frequency: %.2f Hz\n", optimalCutoff_butter);
fprintf("Maximum SIR: %.2f dB\n", maxSIR_butter);

Song_filtered_butter = fftshift(fft(song_filtered_butter))*Ts;
plot_signal_spectra(f, Song_filtered_butter, "Signal filtered using " + ...
    "4^{th} order Butterworth filter: "+title_song)

plot_SIR(sir_butter, cutoffFrequencies_butter, maxSIR_butter, ...
    optimalCutoff_butter)

%% Task1: IIR custom filter 
r = 0.99; % Pole radius

% defines the zero to eliminate at best the disturbance noise
z1 = exp(1j*2*pi*f_noise/Fs);
z2 = conj(z1);              % zeroes are complex conjugate
B = poly([z1, z2]); % Zeros
i = 0;
maxSIR_mine = 0 ;
poleFrequencies = 5500:10:6000;
sir_mine = zeros(1, length(poleFrequencies));

for poleFrequengy = poleFrequencies
    i = i+1;
    % calculate the poles
    p1 = r * exp(1j*2*pi*poleFrequengy/Fs);
    p2 = conj(p1);          % poles are complex conjugate

    A = poly([p1, p2]);     % determines the A polynomial
    xFiltered = filter(B, A, song_noise);
    
    delay = finddelay(song_ref, xFiltered);    
    xFiltered = circshift(xFiltered, -delay);
    sir_mine(i) = determine_SIR(xFiltered, song_ref);
    
    if sir_mine(i) > maxSIR_mine
        maxSIR_mine = sir_mine(i);
        optimalPole = poleFrequengy;
        song_filtered_mine = xFiltered;
        best_freq = poleFrequengy;
        A_best_mine = A;
    end    
end

figure
freqz(B, A_best_mine, f, Fs);

H_mine = zpk(tf(B, A_best_mine, Ts))
figure
pzplot(H_mine)
axis("equal")

fprintf("Best frequency to place the poles of the filter: %.2f Hz \n", optimalPole)
fprintf("Maximum SIR obtained by my own filter: %.2f dB \n", maxSIR_mine)

Song_filtered_mine = fftshift(fft(song_filtered_mine))*Ts;
plot_signal_spectra(f, Song_filtered_mine, "Singal filtered using " + ...
    "filter design by placing zeroes and poles: "+title_song)

plot_SIR(sir_mine, poleFrequencies, maxSIR_mine, optimalPole)

%% Task2: share the channel
fprintf("Execution has been paused, please press any key to continue \n")
pause
clc, clear all, close all
[data1, Fs] = audioread("John Lennon - Imagine.mp3");
[data2, Fs] = audioread("ABBA - Mamma Mia.mp3");

% data is a stereo signal; make it a mono audio (consider only one of the
% two column
song1 = data1(:, 1);
song2 = data2(:, 1);
% N1 = length(song1);
% N2 = length(song2);

channel = load("channel.mat");
channel = channel.channel; % Combined channel signal
Ts = 1/Fs;
% Plot the channel spectrum
N = length(channel);
f = (-N/2:N/2-1)*(Fs/N); % Frequency axis
f = f';
t = 0:Ts:(N-1)*Ts;
t = t';
Channel = abs(fftshift(fft(channel)));

plot_channel(f, Channel, "Original")

song1 = zero_padding(song1, N);
song2 = zero_padding(song2, N);
% Amplitude modulation
f_c1 = 5000; % Carrier frequency for song1 (Hz)
f_c2 = 19000; % Carrier frequency for song2 (Hz)
carrier = cos(2*pi*f_c1*t);
carrier2 = cos(2*pi*f_c2*t);
modulated_song1 = song1 .* carrier;
modulated_song2 = song2 .* carrier2;

modulated_song1 = zero_padding(modulated_song1, N);
modulated_song2 = zero_padding(modulated_song2, N);

channel_final = channel + modulated_song1 + modulated_song2;
Channel_final = fftshift(fft(channel_final));
plot_channel(f, Channel_final, "Channel with both the signals")
y1 = channel_final.*carrier;
y2 = channel_final.*carrier2;
cutoffFrequencies_bessel = 2000:10:2500;  % Range of cutoff frequencies

[song_filtered_bessel1, optimalCutoff_bessel1, sir_bessel1, maxSIR_bessel1] = bessel_filter(y1, song1, cutoffFrequencies_bessel, filterOrder, Fs);

plot_SIR(sir_bessel1, cutoffFrequencies_bessel, maxSIR_bessel1, optimalCutoff_bessel1)
Song_filtered_bessel1 = fftshift(fft(song_filtered_bessel1));
plot_signal_spectra(f, Song_filtered_bessel1, "")

cutoffFrequencies_bessel2 = 5000:10:7000;
[song_filtered_bessel2, optimalCutoff_bessel2, sir_bessel2, maxSIR_bessel2] = bessel_filter(y2, song2, cutoffFrequencies_bessel2, filterOrder, Fs);
plot_SIR(sir_bessel2, cutoffFrequencies_bessel2, maxSIR_bessel2, optimalCutoff_bessel2)

Song_filtered_bessel2 = fftshift(fft(song_filtered_bessel2));
plot_signal_spectra(f, Song_filtered_bessel2, "")

fprintf("Song1 has been modulates using a carrier with frequency " + ...
    "f_c = %dHz, when demodulating and filtering with a Bessel filter, " + ...
    "the obtained signal has a SIR of %.2fdB\n", f_c1, maxSIR_bessel1);

fprintf("Songs has been modulates using a carrier with frequency " + ...
    "f_c = %dHz, when demodulating and filtering with a Bessel filter, " + ...
    "the obtained signal has a SIR of %.2fdnB\n", f_c2, maxSIR_bessel2);

%% Task2: Bessel filter, one song at a time
fprintf("Execution has been paused, please press any key to continue \n")
pause
clc, close all
fc = 16000;
carrier = cos(2*pi*fc*t);
channel_song1 = channel + song1.*carrier;
cut_off_range = 5000:10:6000;
channel_song1_demod = channel_song1 .* carrier;
[song1_alone, cut_off_1alone, maxSIR_1alone, sir_1alone] = bessel_filter(channel_song1_demod, song1, cut_off_range, filterOrder, Fs);

Song1_alone = fftshift(fft(song1_alone));
plot_signal_spectra(f, Song1_alone, "Song1")
plot_SIR(sir_1alone,cut_off_range, maxSIR_1alone, cut_off_1alone)

%% Task2: IIR custom filter,both song together
fprintf("Execution has been paused, please press any key to continue \n")
pause
clc, close all;

f_c1 = 11000;
f_c2 = 19000;
carrier = cos(2*pi*f_c1*t);
carrier2 = cos(2*pi*f_c2*t);
modulated_song1 = song1 .* carrier;
modulated_song2 = song2 .* carrier2;
channel_final = channel + modulated_song1 + modulated_song2;
Channel_final = fftshift(fft(channel_final));
% plot_channel(f, Channel_final, "Channel with both the signals")
y1 = channel_final.*carrier;
y2 = channel_final.*carrier2;
Y1 =fftshift(fft(y1));
Y2 =fftshift(fft(y2));
% plot_channel(f, Y1, "Song1")
% plot_channel(f, Y2, "Song2")

z1 = -1;
z2 = conj(z1);
B = poly([z1, z2]);

poleFrequencies = 600:10:1600;
r = 0.9;
[song1_filt, optimalPole, maxSIR_mine, sir, A_best1, B_best1] = mine_filter(y1, song1, poleFrequencies, Fs, B, r);

plot_SIR(sir, poleFrequencies, maxSIR_mine, optimalPole)

figure
freqz(B_best1, A_best1, f, Fs)

Song1_filt = fftshift(fft(song1_filt));
plot_signal_spectra(f, Song1_filt, "Song1, f_c=11kHz")

H_mine1 = zpk(tf(B_best1, A_best1, Ts));
figure
pzplot(H_mine1)
axis("equal")

poleFrequencies = 2000:10:3000;
r = 0.75;
[song2_filt, optimalPole, maxSIR_mine, sir, A_best2, B_best2] = mine_filter(y2, song2, poleFrequencies, Fs, B, r);

plot_SIR(sir, poleFrequencies, maxSIR_mine, optimalPole)

figure
freqz(B_best2, A_best2, f, Fs)

Song2_filt = fftshift(fft(song2_filt));
plot_signal_spectra(f, Song2_filt, "Song2, f_c=19kHz")

H_mine2 = zpk(tf(B_best2, A_best2, Ts));
figure
pzplot(H_mine2)
axis("equal")

%% Task2: IIR custom filter, one song at a time
fprintf("Execution has been paused, please press any key to continue \n")
pause
clc, close all
cut_off_range = 7500:10:8500;

channel_song2 = channel + song2.*carrier;
channel_song2_demod = channel_song2 .* carrier;
[song2_alone, cut_off_2alone, maxSIR_2alone, sir_2alone] = bessel_filter(channel_song2_demod, song2, cut_off_range, filterOrder, Fs);

Song2_alone = fftshift(fft(song2_alone));
plot_signal_spectra(f, Song2_alone, "Song2")
plot_SIR(sir_2alone, cut_off_range, maxSIR_2alone, cut_off_2alone)


song1 = zero_padding(song1, N);
song2 = zero_padding(song2, N);

fc = 15000;
carrier = cos(2*pi*fc*t);
channel_song1 = channel + song1.*carrier;
channel_song1_demod = channel_song1 .* carrier;

channel_song2 = channel + song2.*carrier;
channel_song2_demod = channel_song2 .* carrier;
Channel_song1_demod = fftshift(fft(channel_song1_demod));
Channel_song2_demod = fftshift(fft(channel_song2_demod));
figure
sgtitle("Demodulated songs")
subplot(2, 1, 1)
plot(f, abs(Channel_song1_demod))
title("Song1"); xlabel("f (Hz)"); ylabel("magnitude")
subplot(2, 1, 2)
plot(f, abs(Channel_song2_demod))
title("Song2"); xlabel("f (Hz)"); ylabel("magnitude")

f_zero = 7150;
z1 = -exp(1j*2*pi*f_zero/Fs);
z2 = conj(z1); 
B = poly([z1, z2]);            

% Initialize SIR tracking
poleFrequencies = 3300:10:4600; 
r = 0.75; % Test pole radii
[song1_filt, optimalPole1, maxSIR_mine1, sir1, A_best1, B_best1] = mine_filter(channel_song1_demod, song1, poleFrequencies, Fs, B, r);
[song_filt, optimalPole2, maxSIR_mine2, sir2, A_best2, B_best2] = mine_filter(channel_song2_demod, song1, poleFrequencies, Fs, B, r);

plot_SIR(sir1, poleFrequencies, maxSIR_mine1, optimalPole1)
plot_SIR(sir2, poleFrequencies, maxSIR_mine2, optimalPole2)

figure
freqz(B_best1, A_best1, f, Fs)
figure
freqz(B_best2, A_best2, f, Fs)

Song1_filt = fftshift(fft(song1_filt));
Song2_filt = fftshift(fft(song2_filt));
plot_signal_spectra(f, Song1_filt, "Song1, f_c=15kHz")
plot_signal_spectra(f, Song2_filt, "Song2, f_c=15kHz")


H_mine1 = zpk(tf(B_best1, A_best1, Ts));
figure
pzplot(H_mine1)
axis("equal")

H_mine2 = zpk(tf(B_best2, A_best2, Ts));
figure
pzplot(H_mine2)
axis("equal")
