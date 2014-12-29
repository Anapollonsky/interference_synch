%% Andrew Apollonsky
% Simulates JMD with two signals of (the same length/time interfering in an
% AWGN channel. Signals simulate wifi-, with 48 subcarrier OFDM, maximum
% data frame length of 2304 octets, and a header consisting of 12/subcar
% training symbols (short/long distinction not simulated) and 24-bits
% including length of data and type of modulation used. Receiver does not
% decode the header, however; it is assumed that it knows everything about
% the transmitted data beforehand. Implements functional convolutional
% encoding, interleaving.

close all;
clc;
clear all;

%% Simulation Parameters
subcardatanum = 48;
framaxoct = 2304; %Maximum frame data length in octets
framelen1 = 2304;  %2304
framelen2 = 2304;
modtype1 = 1;
modtype2 = 1;
snrthing = [-4:1:8];
sirs = [-5 -2 2 5];

if modtype1 < 3
    m1 = 2;
elseif modtype1 < 5
    m1 = 4;
elseif modtype1 < 7
    m1 = 16;
else
    m1 = 64;
end

%% Objects 
% Modulators/Demodulators
bpskmod = comm.BPSKModulator(0);
bpskdemod = comm.BPSKDemodulator(0);
qam4mod = comm.RectangularQAMModulator(4,'BitInput',true);
qam4demod = comm.RectangularQAMDemodulator(4,'BitOutput',true);
qam16mod = comm.RectangularQAMModulator(16,'BitInput',true);
qam16demod = comm.RectangularQAMDemodulator(16,'BitOutput',true);
qam64mod = comm.RectangularQAMModulator(64,'BitInput',true);
qam64demod = comm.RectangularQAMDemodulator(64,'BitOutput',true);

bpskdecmod = bpskmod;
bpskdecdemod = bpskdemod;
qam4decmod = comm.RectangularQAMModulator(4);
qam4decdemod = comm.RectangularQAMDemodulator(4);
qam16decmod = comm.RectangularQAMModulator(16);
qam16decdemod = comm.RectangularQAMDemodulator(16);
qam64decmod = comm.RectangularQAMModulator(64);
qam64decdemod = comm.RectangularQAMDemodulator(64);

% Convolutional Encoders/Decoders
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
convint = comm.ConvolutionalInterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [-1 -2 -3]');
convdeint = comm.ConvolutionalDeinterleaver('NumRegisters', 3, ...
    'RegisterLengthStep', 2, 'InitialConditions', [-1 -2 -3]');
intdelay = -12;

% MIMO Processing
mimoenc = comm.OSTBCEncoder('NumTransmitAntennas', 2);
mimodec = comm.OSTBCCombiner('NumTransmitAntennas',2,'NumReceiveAntennas',2);
%% MIMO Channel Specification (Realistic 802.11n)

S = RandStream('swb2712', 'Seed', 12345); % Set a local random number stream
hModem = modem.pskmod(m1);   % 2-PSK modulator object

Rsym = 10e3;                % Input symbol rate
Rbit = Rsym * log2(m1);      % Input bit rate
Nos  = 4;                   % Oversampling factor
Rs   = Rbit * Nos;          % Input sample rate


tau = [0 10 20 30 40 50 60 70 80] * 1e-9; % Path delays, in seconds

% Average path gains of cluster 1, in dB
pdb1 = [0 -5.4 -10.8 -16.2 -21.7 -inf -inf -inf -inf];
% Average path gains of cluster 2, in dB
pdb2 = [-inf -inf -3.2 -6.3 -9.4 -12.5 -15.6 -18.7 -21.8];
% Total average path gains for both clusters, in dB
pdb = 10*log10(10.^(pdb1/10)+10.^(pdb2/10));

fd = 3;             % Maximum Doppler shift for all paths (identical)
ds = doppler.bell;  % Bell doppler spectrum, with default parameters

Nt = 2;             % Number of transmit antennas
Nr = 2;             % Number of receive antennas

% Element spacing at the transmit and receive antennas (normalized by the
% wavelength)
TxSpacing = 0.5;
RxSpacing = 0.5;

% Spatial parameters on transmitter side:
%   Angular spreads - Cluster 1
AS_Tx_C1 = [14.4 14.4 14.4 14.4 14.4 -inf -inf -inf -inf];
%   Angular spreads - Cluster 2
AS_Tx_C2 = [-inf -inf 25.4 25.4 25.4 25.4 25.4 25.4 25.4];
%   Mean angles of departure - Cluster 1
AoD_C1 = [225.1 225.1 225.1 225.1 225.1 -inf -inf -inf -inf];
%   Mean angles of departure - Cluster 2
AoD_C2 = [-inf -inf 106.5 106.5 106.5 106.5 106.5 106.5 106.5];

% Spatial parameters on receiver side:
%   Angular spreads - Cluster 1
AS_Rx_C1 = [14.4 14.4 14.4 14.4 14.4 -inf -inf -inf -inf];
%   Angular spreads - Cluster 2
AS_Rx_C2 = [-inf -inf 25.2 25.2 25.2 25.2 25.2 25.2 25.2];
%   Mean angles of arrival - Cluster 1
AoA_C1 = [4.3 4.3 4.3 4.3 4.3 -inf -inf -inf -inf];
%   Mean angles of arrival - Cluster 2
AoA_C2 = [-inf -inf 118.4 118.4 118.4 118.4 118.4 118.4 118.4];

% Calculation of transmit and receive correlation arrays
[TxCorrelationMatrix, RxCorrelationMatrix] = ...
    calculateCorrMatrix(Nt, Nr, pdb1, pdb2, TxSpacing, RxSpacing, ...
    AS_Tx_C1, AS_Tx_C2, AoD_C1, AoD_C2, ...
    AS_Rx_C1, AS_Rx_C2, AoA_C1, AoA_C2);

chan2 = comm.MIMOChannel( ...
        'SampleRate',                Rs, ...
        'PathDelays',                tau, ...
        'AveragePathGains',          pdb, ...
        'MaximumDopplerShift',       fd, ...
        'DopplerSpectrum',           ds, ...
        'NumTransmitAntennas',       Nt, ...
        'NumReceiveAntennas',        Nr, ...
        'TransmitCorrelationMatrix', TxCorrelationMatrix, ...
        'ReceiveCorrelationMatrix',  RxCorrelationMatrix, ...
        'RandomStream',              'mt19937ar with seed', ...
        'Seed',                      99, ...
        'PathGainsOutputPort',       true);
%% MIMO Channel Specification (Simplistic)
N=2;
M=2;
chan2 = comm.MIMOChannel('MaximumDopplerShift',       0, ...
                         'NumTransmitAntennas',       N, ...
                         'NumReceiveAntennas',        M, ...
                         'TransmitCorrelationMatrix', eye(N), ...
                         'ReceiveCorrelationMatrix',  eye(M), ...
                         'PathGainsOutputPort',       true);
%%
hold all;
for p = 1:4
    for n = 1:13
        sir = sirs(p);
%         sir = 10;
        snr1 = snrthing(n);
%         snr1 = 40;
        %% Signal Generation 
        %Signal 1

        % Training Symbols
        trainingsymbolser1 = randi([0 1], 12*subcardatanum, 1);
        trainingsymbolpar1 = step(bpskmod, trainingsymbolser1);
        trainingsymbolparmod1 = reshape(trainingsymbolpar1, 12, []);
        % Preamble - Currently includes 3-bit modulation type & 16-bit octet length
        pream1 = [dec2bin(modtype1-1, 3)-48 dec2bin(framelen1,16)-48 ...
            randi([0, 1], 1, 5)] + 0;
        preamod1 = wifimod(pream1.',1,bpskmod,qam4mod,qam16mod,qam64mod,...
            conenc12, conenc23, conenc34, convint, intdelay);
        % Entire header. Includes Training Symbols and Preamble
        header1 = [trainingsymbolparmod1; preamod1.'];
        % Data
        data1 = randi([0 1], framelen1 * 8, 1);
        datamod1 = wifimod(data1,modtype1,bpskmod,qam4mod,qam16mod,qam64mod,...
            conenc12, conenc23, conenc34, convint, intdelay);
        datamodpar1 = reshape(datamod1, [], subcardatanum);
        % Entire signal. Header + Data
        sigf1 = [header1; datamodpar1];
        sigf1ser = reshape(sigf1, [], 1);

        %Signal 2

         % Training Symbols
        trainingsymbolser2 = randi([0 1], 12*subcardatanum, 1);
        trainingsymbolpar2 = step(bpskmod, trainingsymbolser1);
        trainingsymbolparmod2 = reshape(trainingsymbolpar2, 12, []);
        % Preamble
        pream2 = [dec2bin(modtype1-1, 3)-48 dec2bin(framelen1,16)-48 ...
            randi([0, 1], 1, 5)] + 0;
        preamod2 = wifimod(pream2.',1,bpskmod,qam4mod,qam16mod,qam64mod, ...
            conenc12, conenc23, conenc34, convint, intdelay);
        % Entire header. Includes Training Symbols and Preamble
        header2 = [trainingsymbolparmod1; preamod2.'];
        % Data
        data2 = randi([0 1], framelen2 * 8, 1);
        datamod2 = wifimod(data2,modtype2,bpskmod,qam4mod,qam16mod,qam64mod,...
            conenc12, conenc23, conenc34, convint, intdelay);
        datamodpar2 = reshape(datamod2, [], subcardatanum);
        % Entire signal. Header + Data
        sigf2 = [header2; datamodpar2];
        sigf2ser = reshape(sigf2, [], 1);

        %% OFDM Pre-Processing

        % Signal 1
        datalen1 = length(sigf1');
        fftpt1 = datalen1;
        mu1 = 16; % Cyclic Prefix Length
        sigt1 = ifft(sigf1, fftpt1, 1);
        sigtcyc1 = [sigt1(1:mu1, :); sigt1(:,:)];
        sigtcycser1 = reshape(sigtcyc1, [], 1);

        % Signal 2
        datalen2 = length(sigf2');
        fftpt2 = datalen2;
        mu2 = 16; % Cyclic Prefix Length
        sigt2 = ifft(sigf2, fftpt2, 1);
        sigtcyc2 = [sigt2(1:mu2, :); sigt2(:,:)];
        sigtcycser2 = reshape(sigtcyc2, [], 1);

        %% Channel Things
        chan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
        g1 = 1;
        g2 = 10^(-sir/20);
        chan1.SignalPower = (real(sigtcycser1*g1)'*real(sigtcycser1*g1))/ ...
            length(sigtcycser1);
        chan1.SNR = snr1;
        rec1 = step(chan1, sigtcycser1);


        %% OFDM Post-Processing (Single Signal)
        recparcyc1 = reshape(rec1, [], subcardatanum);
        recpar1 = recparcyc1(mu1+1:end, :);
        sigf1rec = fft(recpar1, fftpt1, 1);
        
        recparcyc2 = reshape(sigtcycser2*g2, [], subcardatanum);
        recpar1 = recparcyc1(mu1+1:end, :);
        sigf1rec = fft(recpar1, fftpt1, 1);

        %% Received Signal Deconstruction (Single Signal)
        trainingserdemod1 = step(bpskdemod, reshape(sigf1rec(1:12, :), [], 1));
        preambleserdemod1 = wifidemod(reshape(sigf1rec(13, :), [], 1),1,bpskdemod,...
            qam4demod,qam16demod,qam64demod,condec16,condec78);
        dataserdemod1 = wifidemod(reshape(sigf1rec(14:end, :), [], 1),1,bpskdemod,...
            qam4demod,qam16demod,qam64demod,condec16,condec78);

        %% OFDM Post-Processing (Two Signals)
%         recparcyc = reshape(g1*rec1 + g2*sigtcycser2, [], subcardatanum);
%         recpar = recparcyc(mu1+1:end, :);
%         sigfrec = fft(recpar, fftpt1, 1);

        %% Received Signal Deconstruction (Two Signals)
% 
%         [dataserdemod1 dataserdemod2] = ...
%             wifidemodjmd(reshape(sigfrec(14:end, :), [], 1)...
%             , modtype1, modtype2, g1, g2, bpskmod, qam4mod, qam16mod, qam64mod,...
%             bpskdecmod, qam4decmod, qam16decmod, qam64decmod,...
%             bpskdemod, qam4demod, qam16demod, qam64demod,...
%             condec12, condec23, condec34, convdeint);
        
        %% Plot Stuff
        ber = step(errcalc12, data1, dataserdemod1);
        ser(n) = ber(1);
        reset(errcalc12);
    end
        plot(snrthing, ser);
end

set(gca,'yscale','log');
ylabel('SER');
xlabel('SNR');
grid on;
title('QPSK 1/2 802.11G SNR vs SER');
legend('SIR = -5', 'SIR = -2', 'SIR = 2', 'SIR = 5');
hold off;
