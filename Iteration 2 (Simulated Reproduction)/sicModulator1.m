%% Andrew Apollonsky
% This script simulates a PAM signal coupled with an interference signal
% passing through an AWGN channel. The SIC detector is used to
% demodulate the signal.
close all;
clear all;
clc;

%% Configuration
m1 = 4;   % PAM Constellation Size of main signal
m2 = 4;
len = 10000;    % Signal length
slen1 = log2(m1); % Bits per symbol
slen2 = log2(m2);

EbNo1 = (-3-10*log10(slen1)):2:(30-3-10*log10(slen1));    % Signal Eb/No
sir = [0 3 6 12 60]'; 
snr1 = EbNo1 + 3 + 10*log10(slen1); % Signal SNR Calculated

%% Lessen later workload
g2 = zeros(1, length(sir));
y2 = zeros(len);
w = zeros(len, length(sir), length(EbNo1));
z = zeros(len, length(sir), length(EbNo1));
z2 = zeros(len, length(sir), length(EbNo1));
z3 = zeros(len, length(sir), length(EbNo1));
bervec = zeros(length(sir), length(EbNo1), 3);

%% Prepare Communications things
% PAMModulator / PAMDemodulator
% RectangularQAMModulator / RectangularQAMDemodulator

h1 = comm.PAMModulator(m1); % PAM Modulator
h2 = comm.PAMModulator(m2);
h3 = comm.PAMDemodulator(m1); % PAM Demodulator
h4 = comm.PAMDemodulator(m2); % PAM Demodulator for interference
hchan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)'); ...
    % AWGN Channel with noise set by SNR
hErrorCalc = comm.ErrorRate; % Thing that will calculate error

%% Signal Generation
x1 = randi([0 m1-1], len, 1); % Generate primary signal
x2 = randi([0 m2-1], len, 1); % Generate interference

y1 = step(h1, x1); % Modulate primary signal
hchan1.SignalPower = (real(y1)' * real(y1))/length(real(y1)); ...
    % Find Signal Power of Primary Signal
g1 = 1;   % Gain of primary signal channel
for k = 1:length(sir) % Modulate interference; amplitude depends on SIR
    y2 = step(h2, x2); % Sent signal, gain not included.
    g2(k) = 10^(-sir(k)/20); % Gain of interference
end

%% Demodulation SIC2 (one level further than normal SIC, better for high SIR).
% Might be pretty much identical to II...
% for k = 1:length(snr1)
%     hchan1.SNR = snr1(k); % Set SNR of primary channel
%     for n = 1:length(sir)
%         reset(hErrorCalc);   % Reset the thing that will calculate error
%         w(:, n, k) = step(hchan1, real(y1)) + y2*g2(n); % Pass signal through channel
%         z(:, n, k) = step(h3, complex(w(:, n, k)));% Demodulate
%         z2(:, n, k) = step(h4, complex(w(:, n, k) - step(h1, z(:, n, k)))); %Guess int
%         z3(:, n, k) = step(h3, w(:, n, k) - step(h2, z2(:, n, k)));
%         bervec(n, k, :) = step(hErrorCalc, x1, z3(:, n, k)); % Calculate SER
%     end
% end

%% Demodulation 
for k = 1:length(snr1)
    hchan1.SNR = snr1(k); % Set SNR of primary channel
    for n = 1:length(sir)
        reset(hErrorCalc);   % Reset the thing that will calculate error
        w(:, n, k) = step(hchan1, real(y1)) + y2*g2(n); % Pass signal through channel
        z(:, n, k) = step(h4, complex(w(:, n, k)));% Call this interference
        z2(:, n, k) = step(h3, complex(w(:, n, k) - step(h2, z(:, n, k)))); %Signal
        bervec(n, k, :) = step(hErrorCalc, x1, z2(:, n, k)); % Calculate SER
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
title('SNR vs SER for SIR Variations using Successive Interference Cancellation');
        
        