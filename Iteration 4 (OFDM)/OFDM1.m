%% Andrew Apollonsky
close all;
clc;
clear all;

%% Simulation Parameters


ofdmlen = 48;
datalen = 6912;
fftpt = datalen;
mu = 16;
m1 = 2;   
  

data1 = randi([0 m1-1], ofdmlen * datalen, 1);
mod1 = comm.RectangularQAMModulator(m1);
demod1 = comm.RectangularQAMDemodulator(m1);
errcalc1 = comm.ErrorRate;
chan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');    
snr = -10:20;
%% Signal
for p = 1:31
    
    datamod1 = step(mod1, data1);
    datamodpar1 = reshape(datamod1, ofdmlen, []);
    
    datamodparfreq1 = ifft(datamodpar1, fftpt, 2);
    
    datamodparfreqcyc1 = [datamodparfreq1(:, 1:mu) datamodparfreq1(:, :)];
    datamodfreqcyc1 = reshape(datamodparfreqcyc1, [], 1);
    chan1.SignalPower = (real(datamodfreqcyc1)' * real(datamodfreqcyc1))/length(datamodfreqcyc1);  
    chan1.SNR = snr(p);
    rec1 = step(chan1, datamodfreqcyc1);
    recpar1 = reshape(rec1, ofdmlen, []);
    datapar1 = recpar1(:, mu+1:end);
    
    dataf1 = fft(datapar1, fftpt, 2);
    
    datafser1 = reshape(dataf1, [], 1);
    datarec1 = step(demod1, datafser1);
    error = step(errcalc1, data1, datarec1);
    ber(p) = error(1);
    reset(errcalc1);
end

plot(snr, ber);
set(gca,'yscale','log');
ylabel('SER');
xlabel('SNR');
grid on;
title('16-QAM OFDM SNR vs SER');
