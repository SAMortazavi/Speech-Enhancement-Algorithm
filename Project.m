close all
clear 
clc
%% part 1
desiredSampleRate = 44100;
[data, originalSampleRate] = audioread("car_clean_lom.wav");

% If the original sample rate is different, resample the data
if originalSampleRate ~= desiredSampleRate
    data = resample(data, desiredSampleRate, originalSampleRate);
end
len=length(data);
rounded_len=1000*floor(len/1000);
zero = zeros([rounded_len 1]);
zero(1:rounded_len) = data(1:rounded_len);
data = zero;
% Play the audio
% sound(data, desiredSampleRate);
%Noise
white_noise=0.1*randn(rounded_len,1);
noisy_data=data+white_noise;
% sound(noisy_data, desiredSampleRate);

%% Part 2
window_length = 1000;
overlap = window_length / 2;
% ST FFT for noisy data
[S_noisy, F_noisy, T_noisy] = spectrogram(noisy_data, hamming(window_length), overlap, [], desiredSampleRate);
amplitude_noisy_data = abs(S_noisy);
phase_noisy_data = angle(S_noisy);
% ST FFT for noise
[S_noise, F_noise, T_noise] = spectrogram(white_noise, hamming(window_length), overlap, [], desiredSampleRate);
amplitude_noise = abs(S_noise);
phase_noise = angle(S_noise);
% finding Data differences and setting the negative samples to 0
data_differences=amplitude_noisy_data-amplitude_noise;
data_differences(data_differences<0)=0;
% finding ISTFFT for the data and use its real  part for the rest of code
Data_reconstructed = inverseSTFT(sqrt(data_differences).*exp(1i*phase_noisy_data), hamming(window_length), overlap, window_length);
real_data=real(Data_reconstructed);
% sound(real_data, desiredSampleRate);

%% Part 3
wavelet_type = 'db4';  
level = 1;             % Set the decomposition level

% Perform wavelet decomposition
[c, l] = wavedec(real_data, level, wavelet_type);

% Extract detail coefficients
detail_coefficients = detcoef(c, l, level);
lamba=mean(detail_coefficients);

len_co=length(detail_coefficients);
i=1;
% thresholding
while i<len_co
    if abs(detail_coefficients(i))< lamba
        detail_coefficients(i)=0;
    elseif abs(detail_coefficients(i))> lamba
        detail_coefficients(i)=sign(detail_coefficients(i))*abs(abs(detail_coefficients(i)) - lamba);
    end
i=i+1;
end


reconstructed_signal = wrcoef('a', c, l, wavelet_type, level) + wrcoef('d', c, l, wavelet_type, level, detail_coefficients);
sound(reconstructed_signal, desiredSampleRate);



%% plots
figure;
subplot(3,1,1)
plot(data)
subtitle('Original Data')
subplot(3,1,2)
plot(noisy_data)
subtitle('Noisy Data')
subplot(3,1,3)
plot(reconstructed_signal)
subtitle('Reconstructed Data')




%% functions
% because of limitations of ISTFFT of matlab I got some errors so I
% fond another function from the mathwork.com and used it here.

function x = inverseSTFT(S, window, over_size, fft_length)
    % Perform inverse Short-Time Fourier Transform
    % Input:
    %   S: Spectrogram matrix
    %   window: Analysis window
    %   over_size: overlap size
    %   fft_length: Length of the FFT
    % Output:
    %   x: Reconstructed signal

    % Calculate the number of time steps
    num_frames = size(S, 2);

    % Initialize the reconstructed signal
    x = zeros((num_frames - 1) * over_size + fft_length, 1);

    % Perform inverse STFT
    for col = 1:num_frames
        index = (col - 1) * over_size + 1;
        x(index:index + fft_length - 1) = x(index:index + fft_length - 1) + ifft(S(:, col), fft_length) .* window;
    end

    % Normalize the result
    x = x / sum(window.^2);
end


