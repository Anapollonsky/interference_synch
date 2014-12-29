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
global conenc12 condec12 errcalc12 conenc23 condec23...
    errcalc23 conenc34 condec34 errcalc34;
conenc12 = comm.ConvolutionalEncoder;
condec12 = comm.ViterbiDecoder('InputFormat', 'Hard');
delay12 = condec12.TracebackDepth*log2(condec12.TrellisStructure.numInputSymbols);
errcalc12 = comm.ErrorRate('ReceiveDelay', delay12, 'ComputationDelay',3);

conenc23 = comm.ConvolutionalEncoder;
condec23 = comm.ViterbiDecoder('InputFormat', 'Hard');
conenc23.PuncturePatternSource = 'Property';
conenc23.PuncturePattern = [1;1;1;0];
condec23.PuncturePatternSource = 'Property';
condec23.PuncturePattern = conenc23.PuncturePattern;
delay23 = condec23.TracebackDepth*log2(condec23.TrellisStructure.numInputSymbols);
errcalc23 = comm.ErrorRate('ReceiveDelay', delay23, 'ComputationDelay',3);

conenc34 = comm.ConvolutionalEncoder;
condec34 = comm.ViterbiDecoder('InputFormat', 'Hard');
conenc34.PuncturePatternSource = 'Property';
conenc34.PuncturePattern = [1;1;0;1;1;0];
condec34.PuncturePatternSource = 'Property';
condec34.PuncturePattern = conenc34.PuncturePattern;
delay34 = condec34.TracebackDepth*log2(condec34.TrellisStructure.numInputSymbols);
errcalc34 = comm.ErrorRate('ReceiveDelay', delay34, 'ComputationDelay',3);

% Convolutional Interleaver/Deinterleaver 
global convint convdeint intdelay mimoenc mimodec;
convint = comm.ConvolutionalInterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [-1 -2 -3]');
convdeint = comm.ConvolutionalDeinterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [-1 -2 -3]');
intdelay = -12;

% Custom Things
ofdmenc = OFDMPre;
ofdmdec = OFDMPost;

% MIMO Processing
mimoenc = comm.OSTBCEncoder('NumTransmitAntennas', 2);
mimodec = comm.OSTBCCombiner('NumTransmitAntennas',2,'NumReceiveAntennas',2);
%% Magic Looping
hold all;
snrs = -3:3:9;
sirs = [50];
for p = 1:length(sirs)
    for n = 1:length(snrs)
    %% Frame Creation
    data1 = randi([0 1], 1728*8, 1);
    data2 = randi([0 1], 1728*8, 1);
    training = randi([0 1], 12, 48);
    frame1 = wififrame(data1, 3, training);
    frame2 = wififrame(data2, 3, training);
    
    %% OFDM Pre-Processing
    datat1 = step(ofdmenc, frame1.SignalParallel);
    datat2 = step(ofdmenc, frame2.SignalParallel);
    %% Channel Processing
    snr = snrs(n);
    sir = sirs(p);
    chan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
    g1 = 1;
    g2 = 10^(-sir/20);
    datat1ser = reshape(datat1, [], 1);
    chan1.SignalPower = (real(datat1ser*g1)'*real(datat1ser*g1))/length(datat1ser);
    chan1.SNR = snr;

    
    %% Equalization
    
    %% OFDM Post-Processing
    recfreq = step(ofdmdec, step(chan1, g1*datat1) + g2*datat2);
    %% Received Signal Deconstruction
    [dataserdemod1 dataserdemod2] = ...
        wifidemodjmd(reshape(recfreq(14:end, :), [], 1), ...
        frame1.TransmissionRate, frame2.TransmissionRate, g1, g2);
    
    ber = step(errcalc12, frame1.Data, dataserdemod1);
    ser(n) = ber(1);
    reset(errcalc12);
    end
    plot(snrs, ser);
end
set(gca,'yscale','log');
ylabel('SER');
xlabel('SNR');
grid on;
title('QPSK 1/2 802.11g SNR vs SER');
% legend('SIR = -5', 'SIR = -3', 'SIR = -1', 'SIR = 1', 'SIR = 3', 'SIR = 5');
hold off;