%% Andrew Apollonsky
% This script will produce a plot that compares performance of various MIMO
% decoding schemes in a 2x2, QPSK environment. Other configurations
% partially supported; to change to BPSK, replace 'modu = qam4dmod' with
% 'modu = bpskmod'. To change to higher order modulations, change 'numTx'
% and 'numRx'. Not all implementations support other antenna
% configurations. Requires the associated decoder files.
clc;
clear all;
close all force;
%% Configuration
% Basic
datalen = 80;
ebnos = 0:2:20;
numiter = 100;

% Channel
numTx = 2;
numRx = 2;
rS = 1;
maxDopp = 1e-1;
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

% modu = qam4mod;
modu = qam4mod;
if modu == qam4mod
    demod = qam4demod;
    modb = qam4modb;
elseif modu == bpskmod
    demod = bpskdemod;
    modb = bpskmod;
end

m = length(constellation(modu));

% MIMO
mimoenc = MIMOEnc;
mimodec1 = ZFMIMODec;
mimodec2 = MMSEMIMODec;
mimodec3 = MMSESICMIMODec;
mimodec4 = JMDMIMODec;
mimodec5 = comm.SphereDecoder('Constellation', constellation(modu), ...
    'BitTable', [0 0; 1 0; 0 1; 1 1], 'DecisionType', 'Hard');
% [0 1].'

decs = {mimodec1, mimodec2, mimodec3, mimodec4, mimodec5};
needemod = [1 1 0 0 2 1 1];


% Channel
chan1 = comm.MIMOChannel(...
    'SampleRate', rS, ...
    'MaximumDopplerShift', maxDopp, ...
    'PathGainsOutputPort', true, ...
    'TransmitCorrelationMatrix', tranCorrMat, ...
    'ReceiveCorrelationMatrix', recCorrMat);
chanresetiteration = 5;

chan2 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
wB = waitbar(0,'0%');

%% Simulation
hold all;
set(gca ,'yscale','log');
ylabel('SER');
xlabel('E_b/N_o');
grid on;
title('2x2 MIMO QPSK E_b/N_o vs SER');

tests = [length(needemod)];
prog = 0;
datadec2 = zeros(datalen*numTx, 1);
for reciter = 1:tests
    for n = 1:length(ebnos)
        ebno = ebnos(n);
        for p = 1:numiter
            % Generation
            data1 = randi([0 m-1], datalen * min(numTx, numRx), 1);
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
                datadec1 = step(decs{reciter}, rxSig1, squeeze(pathG1), numTx);
            elseif reciter == 2   
                datadec1 = step(decs{reciter}, rxSig1, squeeze(pathG1), ...
                    var(rxSig1(:, 1) - dataRay1(:, 1)), numTx);
            elseif reciter == 3
                datadec1 = step(decs{reciter}, rxSig1, squeeze(pathG1), ...
                    var(rxSig1(:, 1) - dataRay1(:, 1)), modu, demod,...
                    numTx);
            elseif reciter == 4 
                datadec1 = step(decs{reciter}, rxSig1, ...
                    squeeze(pathG1), modu, modb, numTx, numRx);
            elseif reciter == 5
                datadec1 = step(decs{reciter}, rxSig1, squeeze(pathG1));
            elseif reciter == 6
                datadec1 = qrm_mld(rxSig1, squeeze(pathG1), modu, 2);
            elseif reciter == 7
                datadec1 = qrm_mld(rxSig1, squeeze(pathG1), modu, 4);
            end
            
            %% Serialize
            if needemod(reciter) == 1
                datadec3 = reshape(datadec1.', [], 1);
                datadec4 = step(demod, datadec3);
            elseif needemod(reciter) == 0
                datadec4 = reshape(datadec1.', [], 1);
            elseif needemod(reciter) == 2
                datadec2 = reshape(datadec1, log2(m), []).';
                datadec3 = bi2de(datadec2, log2(m));
                datadec4 = zeros(length(datadec1), 1);
                datadec4(1:numRx:length(datadec2)) = datadec3(1:length(datadec2)/2);
                datadec4(2:numRx:length(datadec2)) = datadec3(length(datadec2)/2+1:length(datadec2));
            end
            ber(n, p) = mean(reshape(de2bi(data1, log2(m)), [], 1) ~= reshape(de2bi((datadec4+0), log2(m)), [], 1));
        end
        prog = prog + 1/(length(ebnos) * tests);
        waitbar(prog, wB, strcat(int2str(round(prog*100)), '%'));
    end
    plot(ebnos, mean(ber, 2), 'LineWidth', 2);
end
close(wB);


% legend('ZF', 'MMSE', 'MMSE-SIC', 'JMD', 'Sphere Decoder');
legend('ZF', 'MMSE', 'MMSE-SIC', 'Minimum-Distance', 'Sphere Decoder', 'QRM-MLD, M = 2', 'QRM-MLD, M = 4');

