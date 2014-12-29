%% Andrew Apollonsky
clc;
clear all;
close all;
% hold all;
%% Random variables
cdrt = 1;

if cdrt < 3
    m1 = 2;
elseif cdrt < 5
    m1 = 4;
elseif cdrt < 7
    m1 = 16;
else
    m1 = 64;
end

if cdrt == 1 || cdrt == 3 || cdrt == 5
    cr1 = 1/2;
elseif cdrt == 7
    cr1 = 2/3;
else
    cr1 = 3/4;
end

snrmod1 = log2(m1) * cr1;

mu=16; % Cyclic Prefix Length
ebnos = [0:1:10];
numtx = 2;
numrx = 2;
frmlen = 64;
numfr = 300;
numofdm = 48;

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
for n = 1:length(ebnos)
    %% Preallocation
    ebno = ebnos(n);
    chan2 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
    g1 = 1;
    for p = 1:numfr
        %% Data Creation
        data1 = randi([0 1], frmlen*numofdm, 1);
        datamod1 = ...
            reshape(wifimod(wifint(wifienc(data1, cdrt)), cdrt), [], numofdm);
        chan2.SNR = ebno + 10*log10(snrmod1) +...
            10*log10(numofdm / length(datamod1'));
%         chan2.SNR = ebno + 10*log10(snrmod1);
        datat1 = step(ofdmenc, datamod1);
        
        %% Channel Processing
        datat1ser = reshape(datat1, [], 1);
%         datamod1ser = reshape(datamod1, [], 1);
        chan2.SignalPower = (real(datat1ser*g1)'*real(datat1ser*g1))/length(datat1ser);
%         chan2.SignalPower = (real(datamod1ser*g1)'*real(datamod1ser*g1))/length(datamod1ser);
        %% OFDM Post-Processing
        recfreq = step(ofdmdec, step(chan2, g1*datat1));
        
%         recfreq = step(chan2, g1*datamod1);
%         dataserdemod1 = wifidec(wifideint(wifidemod(reshape(recfreq...
%             , [], 1), cdrt)), cdrt);
        %% Received Signal Deconstruction
        dataserdemod1 = wifidec(wifideint(wifidemod(reshape(recfreq...
            , [], 1), cdrt)), cdrt);
        
        ber = step(errcalc, data1(34:end-34), dataserdemod1(34:end-34));
        ser(n, p) = ber(1);
        reset(errcalc);
    end
end
ser2 = mean(ser, 2);
plot(ebnos, ser2);

set(gca,'yscale','log');
ylabel('SER');
xlabel('EbNo');
grid on;
title('BPSK 1/2 802.11g SNR vs SER');

hold off;