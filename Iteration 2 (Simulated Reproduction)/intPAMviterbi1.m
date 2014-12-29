%% Andrew Apollonsky
% This script simulates a PAM signal coupled with an interference signal
% passing through an AWGN channel. An interference-ignorant detector
% (nearest point on constellation chosen) is used to demodulate the
% received signal. As far as I know, this script is fully functional and
% provides the correct SER vs SNR vs SIR plot for a given scenario. It
% could be easily modified for different SIRs, SNRs, constellation sizes or
% modulation types. 

% close all;
clear all;
clc;

%% Configuration
m1 = 4;   % PAM Constellation Size of main signal
m2 = 4;
len = 50000;    % Signal length
slen1 = log2(m1); % Bits per symbol
slen2 = log2(m2);

EbNo1 = (-3-10*log10(slen1)):2:(30-3-10*log10(slen1));    % Signal Eb/No.
% Weird numbers for good-looking plot bounds.
sir = [0 3 6 12 60]'; 
snr1 = EbNo1 + 3 + 10*log10(slen1); % Signal SNR Calculated

%% Lessen later workload
y2 = zeros(len, length(sir));
w = zeros(len, length(sir), length(EbNo1));
z = zeros(len, length(sir), length(EbNo1));
bervec = zeros(length(sir), length(EbNo1), 3);

%% Prepare Communications things
% PAMModulator / PAMDemodulator
% RectangularQAMModulator / RectangularQAMDemodulator

h1 = comm.PAMModulator(m1); % Modulator
h2 = comm.PAMModulator(m2);
h3 = comm.PAMDemodulator(m1); %  Demodulator
enc1 = comm.ConvolutionalEncoder;
enc2 = comm.ViterbiDecoder('InputFormat', 'Hard');
hchan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)'); ...
    % AWGN Channel with noise set by SNR
hErrorCalc = comm.ErrorRate; % Thing that will calculate error

%% Signal Generation
x1 = randi([0 m1-1], len, 1); % Generate primary signal
x2 = randi([0 m2-1], len, 1); % Generate interference

y1 = step(h1, step(enc1, reshape(de2bi(x1), 1, [])')); % Modulate primary signal
hchan1.SignalPower = (real(y1)' * real(y1))/length(real(y1)); ...
    % Find Signal Power of Primary Signal

for k = 1:length(sir) % Modulate interference; amplitude depends on SIR
    y2(:,k) = step(h2, x2) / (10^(sir(k)/20)); ...
end

%% Demodulation
for k = 1:length(snr1)
    hchan1.SNR = snr1(k); % Set SNR of primary channel
    for n = 1:length(sir)
        reset(hErrorCalc);   % Reset the thing that will calculate error
        w(:, n, k) = step(hchan1, real(y1)); % Use for PAM
%         w(:, n, k) = step(hchan1, y1); % Received main signal. Use for QAM
        z(:, n, k) = step(h3, complex(w(:, n, k) + y2(:,n)));% Use for PAM
%         z(:, n, k) = step(h3, w(:, n, k) + y2(:,n)); %Use for QAM
        z2(:, n, k) = bi2de(step(enc2, reshape(z(:, n, k), 2, [])'));
        bervec(n, k, :) = step(hErrorCalc, x1, z(:, n, k)); % Calculate SER
    end
end

%% Plot
figure;
hold all;
for k = 1:length(sir)
    plot(snr1, bervec(k, :, 1), 'LineWidth', 2);
end
xlabel('SNR');
ylabel('SER');
grid on;
set(gca,'yscale','log');
legend('SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB', 'SIR = 12dB', 'SIR = 60dB');
title('SNR vs SER for SIR Variations with Interference-Ignorant Detector');
        
        