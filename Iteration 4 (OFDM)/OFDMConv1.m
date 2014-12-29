%% Andrew Apollonsky
close all;
clc;
clear all;

%% Simulation Parameters


ofdmlen = 48;
datalen = 128;
fftpt = datalen;
mu = 16;
m1 = 16;   
rateadj = 2;
bitadj = log2(m1);
mod1 = comm.RectangularQAMModulator(m1, 'BitInput', true);
demod1 = comm.RectangularQAMDemodulator(m1, 'BitOutput', true);
conenc = comm.ConvolutionalEncoder;
chan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
condec = comm.ViterbiDecoder('InputFormat', 'Hard');
delay = condec.TracebackDepth*log2(condec.TrellisStructure.numInputSymbols);
errcalc1 = comm.ErrorRate('ReceiveDelay', delay);

%% Signal
for p = 1:30
    data1 = randi([0 1], datalen * ofdmlen, 1);
    encdata1 = step(conenc, data1);
    modata1 = step(mod1, encdata1);
    modatapar1 = reshape(modata1, ofdmlen, []);
    
    modatapart1 = ifft(modatapar1, fftpt*rateadj/bitadj, 2);
    
    modatatcyc1 = [modatapart1(:, 1:mu) modatapart1(:, :)];
    modatat1 = reshape(modatatcyc1, [], 1);
    chan1.SignalPower = (real(modatat1)'*real(modatat1))/length(modatat1);
    chan1.SNR = p;
    rec1 = step(chan1, modatat1);
    recpar1 = reshape(rec1, ofdmlen, []);
    recdecycpar1 = recpar1(:, mu+1:end);
    
    dataf1 = fft(recdecycpar1, fftpt*rateadj/bitadj, 2);
    
    datafser1 = reshape(dataf1, [], 1);
    demodata1 = step(demod1, datafser1);
    recbits1 = step(condec, demodata1);
    err = step(errcalc1, data1, recbits1);
    ber(p) = err(1);
end


plot(1:30, ber);
set(gca,'yscale','log');
ylabel('SER');
xlabel('SNR');
grid on;
title('4-PAM OFDM SNR vs SER');


% mod1 = comm.PAMModulator(m1, 'BitInput', true);
% demod1 = comm.PAMDemodulator(m1, 'BitOutput', true);
% conenc = comm.ConvolutionalEncoder;
% condec = comm.ViterbiDecoder('InputFormat', 'Hard');
% delay = condec.TracebackDepth*log2(condec.TrellisStructure.numInputSymbols);
% hError = comm.ErrorRate('ComputationDelay', 3, 'ReceiveDelay', delay);
% x = randi([0 1], 30, 1);
% y = step(conenc, x);
% p = step(mod1, y);
% q = step(demod1, p);
% z = step(condec, q);
% err = step(hError, x, z);
% err(1);

