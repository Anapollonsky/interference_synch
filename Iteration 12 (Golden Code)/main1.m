%% Andrew Apollonsky
% if exist('wB', 'var')
%     delete(wB);
% %     close(wB);
% end
clc;
clear all;
close all;
% hold all;

%% Configuration
% Basic
datalen = 30;
ebnos = 0:2:20;
iters = 1500;
m1 = 2;

% Channel
numTx = 2;
numRx = 2;
rS = 1e6;
maxDopp = 30;
chanresetiteration = 10;
%% Object Creation
% Modulation
global bpskmod bpskdemod;
bpskmod = comm.BPSKModulator;
bpskdemod = comm.BPSKDemodulator;

% Error
errcalc = comm.ErrorRate;

% MIMO
goldenc = GoldEnc;
dec1 = ML;
decs = {dec1};

% Channel
chan1 = comm.MIMOChannel(...
    'SampleRate', rS, ...
    'MaximumDopplerShift', maxDopp, ...
    'NumTransmitAntennas', numTx, ...
    'NumReceiveAntennas',  numRx, ...
    'PathGainsOutputPort', true);
%     'ReceiveCorrelationMatrix', recCorrMat, ...
    
chan2 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
wB = waitbar(0,'0%');

%% Simulation
hold all;
set(gca ,'yscale','log');
ylabel('SER');
xlabel('EbNo');
grid on;
title('FDFR BPSK SNR vs SER');

tests = [1];
numdon = 0;
for reciter = tests
    for n = 1:length(ebnos)
        ebno = ebnos(n);
        for p = 1:iters
            % Generation
            data1 = randi([0 1], datalen * min(numTx, numRx), 1);
            datamod1 = step(bpskmod, data1);
            dataenc1 = step(goldenc, datamod1);
            
            % Channel Implementation
            if mod(p ,chanresetiteration) == 0
                reset(chan1);
            end
            [dataRay1, pathG1] = step(chan1, dataenc1);          
            chan2.SNR = ebno + 10*log10(log2(m1));
            siz = size(dataRay1);
            chan2.SignalPower = sum(sum(abs(dataRay1.^2)))/(siz(1)*siz(2));
            rxSig1 = step(chan2, dataRay1);

            %Decoding
            datadec1 = step(dec1, rxSig1, squeeze(pathG1));   
            ber = step(errcalc, data1, double(datadec1));
            ser(n, p) = ber(1);
            reset(errcalc);
        end
        prog = numdon / length(tests) + n/(length(ebnos) * length(tests));
        waitbar(prog, wB, strcat(int2str(round(prog*100)), '%'));
    end
    numdon = numdon + 1;
    ser2 = mean(ser, 2);
    plot(ebnos, ser2);
end
delete(wB);
legend('ML');

