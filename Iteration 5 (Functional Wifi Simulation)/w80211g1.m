%% Andrew Apollonsky
close all;
clc;
clear all;

%% Simulation Parameters
ofdmlen = 48;
datalen = 100000; %In Bits
datatrans = 0;
maxbin = 6000;
% framedatalen = (2312); %In octets, 2312 max
m1 = 4;
datalengthinbits = 14; %Number of bits to allocate for frame data length)

data = randi([0 1], datalen, 1);

bpskmod = comm.BPSKModulator(0);
bpskdemod = comm.BPSKDemodulator(0);
qpskmod = comm.RectangularQAMModulator(4);
qpskdemod = comm.RectangularQAMDemodulator(4);
qam16kmod = comm.RectangularQAMModulator(16);
qam16demod = comm.RectangularQAMDemodulator(16);
qam64mod = comm.RectangularQAMModulator(64);
qam64demod = comm.RectangularQAMDemodulator(64);

conenc12 = comm.ConvolutionalEncoder;
condec12 = comm.ViterbiDecoder('InputFormat', 'Hard');
delay = condec12.TracebackDepth*log2(condec12.TrellisStructure.numInputSymbols);
errcalc12 = comm.ErrorRate('ReceiveDelay', delay);

convint = comm.ConvolutionalInterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [-1 -2 -3]');
convdeint = comm.ConvolutionalDeinterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [-1 -2 -3]');

% preamble1 = [trainingsymbols 
trainingsymbolser = randi([0 1], ofdmlen * 12, 1);
trainingsymbolpar = reshape(trainingsymbolser, 12, []);

% signalinfo = [0 1 

numfr1 = ceil(datalen / maxbin);

%% Fragment data into serialized form for frames
for p = 1:numfr1
    if p ~= numfr1
        dataser(:, 1, p) = data(1+maxbin*(p-1):maxbin*(p));
        framedatalen(p) = maxbin;
    else
        dataser(1:datalen - maxbin*(p-1), 1, p) = data(1+maxbin*(p-1):end);
        framedatalen(p) = length(data(1+maxbin*(p-1):end));
    end
end
framedatalenbin = dec2bin(framedatalen, datalengthinbits);
serfrmdata = [repmat([0 1], [numfr1 1])...   % Modulation Type (QPSK)
     repmat([0 0], [numfr1 1])...          % Code Rate (1/2)
     0+framedatalenbin-48 ...
     zeros(numfr1, 24-datalengthinbits-4)];

for p = 1:numfr1
    serfrmdataenc(p, :) = step(conenc12, serfrmdata(p, :).');
    datapar(:, :, p) = reshape(dataser(:, :, p), [], ofdmlen);
end

beg = [repmat(trainingsymbolpar, [1 1 numfr1]); ...
    permute(serfrmdataenc, [3 2 1]) ];   % Concatenate training symbol and...
                                         % signal data in time-subcar-frame
                                         % domain.
                                         
fullmsg = [beg; datapar];
 
% preamble = repmat(trainingsymbols, [numfr1 1])...  % Training Symbols
%     + framedatalenbin...                           % Length of data
%     + repmat([0 1], [numfr1 1])...                 % Modulation Type (QPSK)
%     + repmat([0 1], [numfr1 1]);                   % Code Rate (1/2)






















