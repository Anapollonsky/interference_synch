%% Andrew Apollonsky
% This script simulates a PAM signal coupled with an interference signal
% passing through an AWGN channel. The joint-MD detector is used to
% demodulate the signal. Note that as of now, while accurate, 
close all;
clear all;
clc;

%% Configuration
m1 = 2;   % PAM Constellation Size of main signal
m2 = 2;
len = 15000;    % Signal length
slen1 = log2(m1); % Bits per symbol
slen2 = log2(m2);

EbNo1 = (-3-10*log10(slen1)):2:(30-3-10*log10(slen1));    % Signal Eb/No
sir = [-5 -3 -1 1 3 5]'; 
snr1 = EbNo1 + 3 + 10*log10(slen1); % Signal SNR Calculated

%% Lessen later workload
g2 = zeros(1, length(sir));
y2 = zeros(len);
w = zeros(len, length(sir), length(EbNo1));
z = zeros(len, length(sir), length(EbNo1));
bervec = zeros(length(sir), length(EbNo1), 3);

%% Prepare Communications things
% PAMModulator / PAMDemodulator
% RectangularQAMModulator / RectangularQAMDemodulator

h1 = comm.PAMModulator(m1); % Modulator
h2 = comm.PAMModulator(m2);
h3 = comm.PAMDemodulator(m1); % Demodulator
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

%% Preparation for Demodulation
for k = 1:length(snr1)
    hchan1.SNR = snr1(k); % Set SNR of primary channel
    for n = 1:length(sir)
        w(:, n, k) = step(hchan1, real(y1)) + y2*g2(n); % Received signal
    end
end

%% Demodulation
minimat = zeros(len, m1, m2);
mini = zeros(len, m1, m2, length(sir));
minintime = zeros(1, len);

%This subsection generates the h1x1 + h2x2 term. There needs to be a value
%for every combination of x1 and x2, as well as one for every SIR due to
%the changing h2 (gain, here seen as g2).
for n = 1:length(sir)
    for q = 0:(m1-1) % y1
        for a = 0:(m2-1) % y2
            mini(:, q+1, a+1, n) = g1*step(h1, q) + g2(n)*step(h2, a);
        end
    end
end

% This subsection generates minimat, a matrix that stores matrix
% representing the entire y - h1x1 - h2x2 term for every received signal.
% For every term in time, SIR and SNR, the combination of y1 and y2 values
% that maximizes the chances of the given signal being received is chosen
% by finding the values that minimize the term. The value of y1 which fits
% this is the one that is estimated to be correct. SER is then calculated. 
for k = 1:length(snr1)
    for n = 1:length(sir)
        reset(hErrorCalc);   % Reset the thing that will calculate error
        minimat = abs(repmat(w(:, n, k), [1 m1 m2]) - mini(:, :, :, n));
        minintime = min(min(minimat, [], 2), [], 3);   % Pre-calculates 
        % of minimat 
        for p = 1:len
            [r, c, v] = find(shiftdim(minimat(p, :, :), 1) == minintime(p));
            z(p, n, k) = r(1)-1;
        end
        bervec(n, k, :) = step(hErrorCalc, x1, z(:, n, k)); % Calculate SER
    end
end

%% Old, inefficient demodulation
% Easier to visualize than the mess above.

% for k = 1:length(snr1)
%     for n = 1:length(sir)
%         reset(hErrorCalc);   % Reset the thing that will calculate error
%         for p = 1:len % For every moment in time
%             for q = 0:(m1-1) % For every possible value of y1
%                 for a = 0:(m2-1) % For every possible value of y2
%                     minimat(q+1, a+1) = abs(...  % Calculate the 'term'
%                         w(p, n, k)...
%                         - g1*step(h1, q) ...
%                         - g2(n)*step(h2, a)...
%                         );
%                 end
%             end
%         [r, c, v] = find(minimat == min(min(minimat))); % Find minimum
%         z(p, n, k) = r(1)-1;    % Assign responsible y1 value as estimate
%         end
%         bervec(n, k, :) = step(hErrorCalc, x1, z(:, n, k)); % Calculate SER
%     end
% end

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
legend('SIR = -5dB', 'SIR = -3dB', 'SIR = -1dB', 'SIR = 1dB', 'SIR = 3dB', 'SIR = 5dB');
title('SNR vs SER for SIR Variations using joint-MD');
        
        