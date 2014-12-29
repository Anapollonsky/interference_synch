%% Andrew Apollonsky
clc;
clear all;
close all;
%% Random variables
mu=16; % Cyclic Prefix Length
%% Object Creation
% Modulators/Demodulators
global bpskmod bpskdemod qam4mod qam4demod qam16mod qam16demod qam64mod qam64demod;
bpskmod = comm.BPSKModulator(0);
bpskdemod = comm.BPSKDemodulator(0);
qam4mod = comm.RectangularQAMModulator(4,'BitInput',true);
qam4demod = comm.RectangularQAMDemodulator(4,'BitOutput',true);
qam16mod = comm.RectangularQAMModulator(16,'BitInput',true);
qam16demod = comm.RectangularQAMDemodulator(16,'BitOutput',true);
qam64mod = comm.RectangularQAMModulator(64,'BitInput',true);
qam64demod = comm.RectangularQAMDemodulator(64,'BitOutput',true);

global qam4decmod qam4decdemod qam16decmod qam16decdemod qam64decmod qam64decdemod;
qam4decmod = comm.RectangularQAMModulator(4);
qam4decdemod = comm.RectangularQAMDemodulator(4);
qam16decmod = comm.RectangularQAMModulator(16);
qam16decdemod = comm.RectangularQAMDemodulator(16);
qam64decmod = comm.RectangularQAMModulator(64);
qam64decdemod = comm.RectangularQAMDemodulator(64);

% Convolutional Encoders/Decoders
global conenc12 condec12 conenc23 condec23 conenc34 condec34 errcalc...
    delay12 delay23 delay34;
conenc12 = comm.ConvolutionalEncoder;
condec12 = comm.ViterbiDecoder('InputFormat', 'Hard');
delay12 = condec12.TracebackDepth*log2(condec12.TrellisStructure.numInputSymbols);
errcalc = comm.ErrorRate;

conenc23 = comm.ConvolutionalEncoder;
condec23 = comm.ViterbiDecoder('InputFormat', 'Hard');
conenc23.PuncturePatternSource = 'Property';
conenc23.PuncturePattern = [1;1;1;0];
condec23.PuncturePatternSource = 'Property';
delay23 = condec23.TracebackDepth*log2(condec23.TrellisStructure.numInputSymbols);
condec23.PuncturePattern = conenc23.PuncturePattern;

conenc34 = comm.ConvolutionalEncoder;
condec34 = comm.ViterbiDecoder('InputFormat', 'Hard');
conenc34.PuncturePatternSource = 'Property';
conenc34.PuncturePattern = [1;1;0;1;1;0];
condec34.PuncturePatternSource = 'Property';
delay34 = condec34.TracebackDepth*log2(condec34.TrellisStructure.numInputSymbols);
condec34.PuncturePattern = conenc34.PuncturePattern;

% Convolutional Interleaver/Deinterleaver 
global convint convdeint intdelay mimoenc mimodec;
convint = comm.ConvolutionalInterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [1 1 1]');
convdeint = comm.ConvolutionalDeinterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [1 1 1]');
intdelay = 12;

% Custom Things
ofdmenc = OFDMPre;
ofdmdec = OFDMPost;

% MIMO Processing
mimoenc = comm.OSTBCEncoder('NumTransmitAntennas', 2);
mimodec = comm.OSTBCCombiner('NumTransmitAntennas',2,'NumReceiveAntennas',2);
%% Magic Looping
snrs = 15:.5:25;

for n = 1:length(snrs)
    %% Frame Creation
    data1 = randi([0 1], 1728*8*16, 1);
    training = randi([0 1], 12, 48);
    frame1 = wififrame(data1, 8, training);

    %% OFDM Pre-Processing
    datat1 = step(ofdmenc, frame1.SignalParallel);
    %% Channel Processing
    snr = snrs(n);
    chan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
    g1 = 1;
    datat1ser = reshape(datat1, [], 1);
    chan1.SignalPower = (real(datat1ser*g1)'*real(datat1ser*g1))/length(datat1ser);
    chan1.SNR = snr;
    %% OFDM Post-Processing
    recfreq = step(ofdmdec, step(chan1, g1*datat1));
%     recfreq = recfreq * exp(j*3);
    %% Equalization
    trainingmod = reshape(wifimod(reshape(training, [], 1), 1), [], 48);
    rectrainingmod = recfreq(1:12, :);
    warning('off','MATLAB:lscov:RankDefDesignMat');
    t = lscov(rectrainingmod, trainingmod);
    for k = 1:48
      t(k, :) = lscov(rectrainingmod(:, k), trainingmod(:, k));
    end
    w = mean(t, 1);
    for k = 1:48
      recfreq(13:end, k) = recfreq(13:end, k)*w(k); %* (1/abs(w(k)));
    end
    
%       t = trainingmod/rectrainingmod;
%     recfreq(13:end, :) = recfreq(13:end, :) * t;

%     t = mean(trainingmod ./ rectrainingmod);

    for k = 1:48
        t(k, :) = lsqlin(rectrainingmod(:, k), trainingmod(:, k), 0, 1);
    end
    w = mean(t, 1);
    for k = 1:48
      recfreq(13:end, k) = recfreq(13:end, k)*w(k); %* (1/abs(w(k)));
    end
%     for k = 1:48
%         t(k) =  rectrainingmod(:, k) \ trainingmod(:, k) ;
%     end
%     recfreq(13:end, :) = recfreq(13:end, :) .* ...
%         repmat(t, length(recfreq(13:end, :)'), 1);
    %% Received Signal Deconstruction
    demoheader = wifidemod(reshape(recfreq(13:15, :), [], 1), 1);
    binheader = wifidec(wifideint(demoheader), 1);
    [TransType, ModOrder, Rate, BinLength] = GetHeader(binheader);
    
%     dataserdemod1 = wifidec(wifideint(wifidemod(reshape(...
%         recfreq(16:end, :), [], 1), TransType)), ...
%         TransType);
    
    dataserdemod1 = wifidec(wifideint(wifidemod(reshape(...
        recfreq(16:end, :), [], 1), frame1.TransmissionRate)), ...
        frame1.TransmissionRate);

%     if TransType == frame1.TransmissionRate
%         ber = step(errcalc, frame1.Data(1:end-34), dataserdemod1(1:end-34));
%         ser(n) = ber(1);
%     else
%         ser(n) = 1;
%     end
    ber = step(errcalc, frame1.Data(1:end-34), dataserdemod1(1:end-34));
    ser(n) = ber(1);
    reset(errcalc);
end
plot(snrs, ser);

set(gca,'yscale','log');
ylabel('SER');
xlabel('SNR');
grid on;
title('QPSK 1/2 802.11g SNR vs SER');

hold off;