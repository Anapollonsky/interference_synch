%% Andrew Apollonsky
clc;
clear all;
close all;

%% Configuration
% Basic
datalen = 10000;
ebnos = 0:2:26;
iters = 150;
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

% Convolutional Interleaver/Deinterleaver 
convint = comm.ConvolutionalInterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [1 1 1]');
convdeint = comm.ConvolutionalDeinterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [1 1 1]');
intdelay = 12;

% Error
errcalc = comm.ErrorRate;

% MIMO
mimoenc = MIMOEnc;
mimodec = comm.SphereDecoder('Constellation', constellation(bpskmod), ...
    'BitTable', [0; 1], 'DecisionType', 'Hard');

% Channel
chan1 = comm.MIMOChannel(...
    'SampleRate', rS, ...
    'MaximumDopplerShift', maxDopp, ...
    'NumTransmitAntennas', numTx, ...
    'NumReceiveAntennas',  numRx, ...
    'PathGainsOutputPort', true);
%     'ReceiveCorrelationMatrix', recCorrMat, ...
    
chan2 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
wB = waitbar(0,'Initializing waitbar...');

%% Simulation
for n = 1:length(ebnos)
    for p = 1:iters
        % Generation
        data1 = randi([0 1], datalen * min(numTx, numRx), 1);
        dataenc1 = data1;
        datamod1 = step(bpskmod, dataenc1);
        dataoenc1 = step(mimoenc, datamod1, numTx);

        % Channel Implementation
        ebno = ebnos(n);
        chan2.SNR = ebno + 10*log10(log2(m1));
        chan2.SignalPower = 1;
        g1 = 1;
        [dataRay1, pathG1] = step(chan1, squeeze(dataoenc1));
        rxSig1 = step(chan2, dataRay1);
        
        %Decoding
        datadec1 = step(mimodec, rxSig1, squeeze(pathG1));
        for k = 1:numRx
            datadec2(k:numRx:2*length(datadec1), :) = datadec1(:, k);
        end
%         dademod1 = step(bpskdemod, datadec2);
        ber = step(errcalc, data1, double(datadec2));
        ser(n, p) = ber(1);
        reset(errcalc);
    end
    prog = n/length(ebnos);
    waitbar(prog, wB, strcat(int2str(round(prog*100)), '%'));
end
close(wB);
hold all;
ser2 = mean(ser, 2);
plot(ebnos, ser2);
set(gca ,'yscale','log');
ylabel('SER');
xlabel('EbNo');
grid on;
title('BPSK SNR vs SER');
