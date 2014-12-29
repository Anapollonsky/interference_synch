%% Andrew Apollonsky
clc;
clear all;
close all force;
% hold all;
%% Configuration
% Basic
datalen = 40;
ebnos = 0:4:20;
iters = 6;
M = 2;
% Channel
numTx = 2;
numRx = 2;
rS = 1e6;
maxDopp = 30;
chanresetiteration = 1;
if numRx == numTx
    recCorrMat = eye(numRx);
    tranCorrMat = eye(numRx);
else
    tranCorrMat = [1 0; 0 1; 1 0; 0 1];
    recCorrMat = [1 0; 0 1; 1 0; 0 1].';
end
%% Object Creation
% Modulation
bpskmod = comm.BPSKModulator;
bpskdemod = comm.BPSKDemodulator;

qam4mod = comm.RectangularQAMModulator(4);
qam4demod = comm.RectangularQAMDemodulator(4);
qam4modb = comm.RectangularQAMModulator(4, 'BitInput', true);
qam4demodb = comm.RectangularQAMDemodulator(4, 'BitOutput', true);

modu = bpskmod;
if modu == qam4mod
    modb = qam4modb;
    demod = qam4demod;
elseif modu == bpskmod
    modb = bpskdemod;
    demod = bpskdemod;
end

m = length(constellation(modu));

% Error
errcalc = comm.ErrorRate;

% MIMO
mimoenc = MIMOEnc;
mimodec1 = comm.SphereDecoder('Constellation', constellation(bpskmod), ...
    'BitTable', [0; 1], 'DecisionType', 'Hard');
mimodec2 = JMDMIMODec;
mimodec3 = ZFMIMODec;
mimodec4 = MMSEMIMODec;
mimodec5 = MMSESICMIMODec;
mimodec6 = LLLMMSEMIMODec;
mimodec7 = QRMMLDMIMODec;

decs = {mimodec1, mimodec2, mimodec3, mimodec4, mimodec5, mimodec6, mimodec7};
needemod = [0 0 1 1 0 1 1];


% Channel
chan1 = comm.MIMOChannel(...
    'SampleRate', rS, ...
    'MaximumDopplerShift', maxDopp, ...
    'NumTransmitAntennas', numTx, ...
    'NumReceiveAntennas',  numRx, ...
    'PathGainsOutputPort', true, ...
    'TransmitCorrelationMatrix', tranCorrMat, ...
    'ReceiveCorrelationMatrix', recCorrMat);
    
chan2 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
wB = waitbar(0,'0%');

%% Simulation
hold all;
set(gca ,'yscale','log');
ylabel('SER');
xlabel('EbNo');
grid on;
title('2x2 MIMO BPSK SNR vs SER');

tests = [6];
numdon = 0;
for reciter = tests
    for n = 1:length(ebnos)
        ebno = ebnos(n);
        for p = 1:iters
            % Generation
            data1 = randi([0 1], datalen * min(numTx, numRx), 1);
%             dataenc1 = convenc(data1, conenc12);
            datamod1 = step(modu, data1);
            dataoenc1 = step(mimoenc, datamod1, numTx);

            % Channel Implementation
            if mod(p ,chanresetiteration) == 0
                reset(chan1);
            end
            [dataRay1, pathG1] = step(chan1, dataoenc1);          
            chan2.SNR = ebno + 10*log10(log2(m));
            siz = size(dataRay1);
            chan2.SignalPower = sum(sum(abs(dataRay1.^2)))/(siz(1)*siz(2));
            rxSig1 = step(chan2, dataRay1);

            %Decoding
            if reciter == 1 
                datadec1 = step(decs{reciter}, rxSig1, squeeze(pathG1));
            elseif reciter == 2
                datadec1 = step(decs{reciter}, rxSig1, ...
                    squeeze(pathG1), modu, modb, numTx, numRx);
            elseif reciter == 3
                datadec1 = step(decs{reciter}, rxSig1, squeeze(pathG1), numTx);
            elseif reciter == 4 
                datadec1 = step(decs{reciter}, rxSig1, squeeze(pathG1), ...
                    var(rxSig1(:, 1) - dataRay1(:, 1)), numTx);
            elseif reciter == 5
                datadec1 = step(decs{reciter}, rxSig1, squeeze(pathG1), ...
                    var(rxSig1(:, 1) - dataRay1(:, 1)), bpskmod, bpskdemod,...
                    numTx);
            elseif reciter == 6
                datadec1 = step(decs{reciter}, rxSig1, squeeze(pathG1), ...
                    var(rxSig1(:, 1) - dataRay1(:, 1)));
            elseif reciter == 7
                datadec1 = step(decs{reciter}, rxSig1, squeeze(pathG1), ...
                    numTx, modu, M);
            end
            
            
            for k = 1:numRx
                if needemod(reciter) == 1
                    datadec3(k:numRx:numRx*length(datadec1), :)...
                        = datadec1(:, k);
                else
                    datadec2(k:numRx:numRx*length(datadec1), :) = datadec1(:, k);
                end
            end
            
            if needemod(reciter) == 1
                datadec4 = step(bpskdemod, datadec3);
            else
                datadec4 = datadec2;
            end
           
            ber = step(errcalc, data1, double(datadec4));
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
close(wB);

% legend('Sphere Decoder', 'JMD', 'ZF', 'MMSE', 'MMSE-SIC', 'LLL-MMSE');

legend('Sphere Decoder', 'JMD', 'ZF', 'MMSE', 'MMSE-SIC', strcat('QRM-MLD (M = ', int2str(M), ')'));

