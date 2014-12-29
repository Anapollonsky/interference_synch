%% Andrew Apollonsky
clc;
clear all;
close all;
hold all;
%% Configuration
% Basic
datalen = 1000;
ebnos = 0:2:15;
iters = 50;
m1 = 2;

% Channel
numTx = 2;
numRx = 2;
rS = 1e6;
maxDopp = 30;
% recCorrMat = [1 0; 0 1];

%% Object Creation
% Modulation
bpskmod = comm.BPSKModulator;
bpskdemod = comm.BPSKDemodulator;

% Convolutional Encoding
conenc12 = comm.ConvolutionalEncoder;
condec12 = comm.ViterbiDecoder('InputFormat', 'Hard');
delay12 = condec12.TracebackDepth*log2(condec12.TrellisStructure.numInputSymbols);

% MIMO
ostbcenc = comm.OSTBCEncoder('NumTransmitAntennas', numTx);
ostbcdec = comm.OSTBCCombiner('NumTransmitAntennas', numTx, 'NumReceiveAntennas', numRx);

% Convolutional Interleaver/Deinterleaver 
convint = comm.ConvolutionalInterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [1 1 1]');
convdeint = comm.ConvolutionalDeinterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [1 1 1]');
intdelay = 12;

% Error
errcalc = comm.ErrorRate;

% Automatic Gain Control
AGC1 = comm.AGC;

% Channel
chan1 = comm.MIMOChannel(...
    'SampleRate', rS, ...
    'MaximumDopplerShift', maxDopp, ...
    'NumTransmitAntennas', numTx, ...
    'NumReceiveAntennas',  numRx, ...
    'PathGainsOutputPort', true);
%     'ReceiveCorrelationMatrix', recCorrMat, ...
    
chan2 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');

%% Beginning Simulations
wB = waitbar(0,'Initializing waitbar...');

conenc = conenc12;
condec = condec12;
delay = delay12;

%% No Convolutional Encoding
for n = 1:length(ebnos)
    for p = 1:iters
        % Generation
        data1 = randi([0 1], datalen * min(numTx, numRx), 1);
        dataenc1 = data1;
        datamod1 = step(bpskmod, dataenc1);
        dataoenc1 = step(ostbcenc, datamod1);

        % Channel Implementation
        ebno = ebnos(n);
        chan2.SNR = ebno + 10*log10(log2(m1));
        chan2.SignalPower = 1;
        g1 = 1;
        [dataRay1, pathG1] = step(chan1, dataoenc1);
        rxSig1 = step(chan2, dataRay1);

        %Decoding
        datadec1 = step(ostbcdec, rxSig1, squeeze(pathG1));
        dademod1 = step(bpskdemod, datadec1);
        ber = step(errcalc, data1, dademod1);
        ser(n, p) = ber(1);
        reset(errcalc);
    end
    prog = n/length(ebnos)/3;
    waitbar(prog, wB, strcat(int2str(round(prog*100)), '%'));
end

%% Plot
ser2 = mean(ser, 2);
plot(ebnos, ser2);
set(gca ,'yscale','log');
ylabel('SER');
xlabel('EbNo');
grid on;
title('BPSK SNR vs SER');

%% Convolutional Encoding
release(chan1);
release(chan2);
for n = 1:length(ebnos)
    for p = 1:iters
        % Generation
        data1 = randi([0 1], datalen * min(numTx, numRx), 1);
        dataenc1 = convenc(data1, conenc);
%         dataenc1 = data1;
        datamod1 = step(bpskmod, dataenc1);
        dataoenc1 = step(ostbcenc, datamod1);

        % Channel Implementation
        ebno = ebnos(n);
        chan2.SNR = ebno + 10*log10(log2(m1));
        chan2.SignalPower = 1;
        g1 = 1;
        [dataRay1, pathG1] = step(chan1, dataoenc1);
        rxSig1 = step(chan2, dataRay1);

        %Decoding
        datadec1 = step(ostbcdec, rxSig1, squeeze(pathG1));
        dademod1 = step(bpskdemod, datadec1);
        datadec1 = convdec(dademod1, condec, delay);
        ber = step(errcalc, data1(1:end-delay), datadec1(1:end-delay));
        ser(n, p) = ber(1);
        reset(errcalc);
    end
    prog = .33 + n/length(ebnos)/3;
    waitbar(prog, wB, strcat(int2str(round(prog*100)), '%'));
end

%% Plot
ser2 = mean(ser, 2);
plot(ebnos, ser2);

%% Convolutional Encoding + Interleaving
release(chan1);
release(chan2);
for n = 1:length(ebnos)
    for p = 1:iters
        % Generation
        data1 = randi([0 1], datalen * min(numTx, numRx), 1);
        dataenc1 = convenc(data1, conenc);
        dataint1 = inter(dataenc1, convint, intdelay);
        datamod1 = step(bpskmod, dataint1);
        dataoenc1 = step(ostbcenc, datamod1);

        % Channel Implementation
        ebno = ebnos(n);
        chan2.SNR = ebno + 10*log10(log2(m1));
        chan2.SignalPower = 1;
        g1 = 1;
        [dataRay1, pathG1] = step(chan1, dataoenc1);
        rxSig1 = step(chan2, dataRay1);

        %Decoding
        datadec1 = step(ostbcdec, rxSig1, squeeze(pathG1));
        dademod1 = step(bpskdemod, datadec1);
        datadeint1 = deinter(dademod1, convdeint);
        datadec1 = convdec(datadeint1, condec, delay);
        ber = step(errcalc, data1(1:end-delay-intdelay), datadec1(1:end-delay-intdelay));
        ser(n, p) = ber(1);
        reset(errcalc);
    end
    prog = .67 + n/length(ebnos)/3;
    waitbar(prog, wB, strcat(int2str(round(prog*100)), '%'));
end
close(wB)
%% Plot
ser2 = mean(ser, 2);
plot(ebnos, ser2);
legend('No Encoding', 'Convolutional Encoding', 'Encoding + Interleaving');
