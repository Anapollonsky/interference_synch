%% Andrew Apollonsky
close all;
clc;
clear all;

%% Simulation Parameters


ofdmlen = 48;
datalen = 2000;
fftpt = datalen;
sirs = [-5 -3 -1 1 3 5];
snrs = [0:2:14];
mu = 16;
m1 = 4;   
m2 = 4;

data1 = randi([0 m1-1], ofdmlen * datalen, 1);
mod1 = comm.RectangularQAMModulator(m1);
demod1 = comm.RectangularQAMDemodulator(m1);
errcalc1 = comm.ErrorRate;
chan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');

data2 = randi([0 m1-1], ofdmlen * datalen, 1);
mod2 = comm.RectangularQAMModulator(m2);
demod2 = comm.RectangularQAMDemodulator(m2);
errcalc2 = comm.ErrorRate;
chan2 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
    
g1 = 1;

minimat = zeros(m1, m2, datalen * ofdmlen);


for siriter = 1:6  % Iterate through SIR's
    g2 = 10^(-sirs(siriter)/20);   
    for p = 1:8  % Iterate through SNRs

        %% Signal 1
        datamod1 = step(mod1, data1);
        datamodpar1 = reshape(datamod1, ofdmlen, []);
        datamodparfreq1 = ifft(datamodpar1, fftpt, 2);   
        datamodparfreqcyc1 = [datamodparfreq1(:, 1:mu) datamodparfreq1(:, :)];
        datamodfreqcyc1 = reshape(datamodparfreqcyc1, [], 1);
        chan1.SignalPower = (real(datamodfreqcyc1)' * real(datamodfreqcyc1))/length(datamodfreqcyc1);  
        chan1.SNR = snrs(p);
        rec1 = step(chan1, datamodfreqcyc1);

        %% Interference
        datamod2 = step(mod1, data1);
        datamodpar2 = reshape(datamod2, ofdmlen, []);
        datamodparfreq2 = ifft(datamodpar2, fftpt, 2);   
        datamodparfreqcyc2 = [datamodparfreq2(:, 1:mu) datamodparfreq2(:, :)];
        datamodfreqcyc2 = reshape(datamodparfreqcyc2, [], 1);
        chan2.SignalPower = (real(datamodfreqcyc2)' * real(datamodfreqcyc2))/length(datamodfreqcyc2);  
        chan2.SNR = snrs(p)-sirs(siriter);
        rec2 = step(chan2, datamodfreqcyc2);

        %% Combined
        rectot = rec1*g1 + rec2*g2;

        recpar = reshape(rectot, ofdmlen, []);
        datapar = recpar(:, mu+1:end);

        dataftot = fft(datapar, fftpt, 2);

        datafser = reshape(dataftot, [], 1);

        for t = 1:length(datafser)
            for q = 0:(m1-1) % For every possible value of y1
                for a = 0:(m2-1) % For every possible value of y2
                    minimat(q+1, a+1, :) = abs(...  % Calculate the 'term'
                        datafser(t)...
                        - g1*step(mod1, q) ...
                        - g2*step(mod2, a)...
                        );
                end
            end
            [r, c, v] = find(minimat(:, :, t) == min(min(minimat(:, :, t)))); % Find minimum
    %         z(p, n, k) = r(1)-1;    % Assign responsible y1 value as estimate
            datarec(t) = r(1)-1;
        end


    %     datarec = step(demod1, datafser);
        error = step(errcalc1, data1, datarec.');
        ber(siriter, p) = error(1);
        reset(errcalc1);
    end
end
hold all;
for k = 1:length(sirs)
    plot(snrs, ber(k, :), 'LineWidth', 2);
end
legend('SIR = -5', 'SIR = -3', 'SIR = -1', 'SIR = 1', 'SIR = 3', 'SIR = 5');
set(gca,'yscale','log');
ylabel('SER');
xlabel('SNR');
grid on;
title('OFDM JMD SIR & SNR vs SER');
