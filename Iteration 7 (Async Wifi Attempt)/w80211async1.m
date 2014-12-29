%% Andrew Apollonsky
clc;
clear all;
close all;
%% Configurables
mu=16; % Cyclic Prefix Length (16 used no matter what for most things)
snrs = -3:3:9;
sirs = [50];
delay = 50;
numsigs = 1;
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
intdelay = -12;

% Custom Things
ofdmenc = OFDMPre;
ofdmdec = OFDMPost;

% MIMO Processing
mimoenc = comm.OSTBCEncoder('NumTransmitAntennas', 2);
mimodec = comm.OSTBCCombiner('NumTransmitAntennas',2,'NumReceiveAntennas',2);
%% Magic Looping
% hold all;

%% Frame Creation
data1 = randi([0 1], 1728*8, 1);
data2 = randi([0 1], 1728*8, 1);
training = randi([0 1], 12, 48);
frame1 = wififrame(data1, 1, training);
frame2 = wififrame(data2, 1, training);

%% OFDM Pre-Processing
datat1 = step(ofdmenc, [frame1.SignalParallel;...
    zeros(length(frame2.SignalParallel') + delay -...
    length(frame1.SignalParallel'), 48)]);
datat2 = step(ofdmenc, [zeros(delay, 48); frame2.SignalParallel]);

%% Channel Processing
% snr = snrs(n);
% sir = sirs(p);
snr = 90; sir = 90;
chan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
g1 = 1;
g2 = 10^(-sir/20);
datat1ser = reshape(datat1, [], 1);
datat2ser = reshape(datat2, [], 1);
chan1.SignalPower = (real(datat1ser*g1)'*real(datat1ser*g1))/length(datat1ser);
chan1.SNR = snr;

%% Signal Combination
phaseoffset1 = 1;
phaseoffset2 = 1;
% rectime = zeros(max(length(datat1'), length(datat2') + delay), frame1.OFDMLen);
% rectime(1:length(datat1'), :) = datat1 * g1 * phaseoffset1;
% rectime(delay+1:end, :) = rectime(delay+1:end, :) + g2 * datat2 * phaseoffset2;
rectime = step(chan1, datat1*g1) + datat2*g2;

%% OFDM Post-Processing
recfreqnorm = step(ofdmdec, rectime);      % OFDM Post-Processing

%% Equalization
% trainingmod = reshape(wifimod(reshape(training, [], 1), 1), [], 48);
% rectrainingmod = recfreqnorm(1:12, :);
% warning('off','MATLAB:lscov:RankDefDesignMat');
% eqmat1 = lscov(rectrainingmod, trainingmod);
% 
% for k = 1:48
%     eqmat1(k, :) = lscov(rectrainingmod(:, k), trainingmod(:, k));
% end
% w = mean(eqmat1, 1);
% for k = 1:48
%     recfreqnorm(16:end, k) = recfreqnorm(16:end, k)*w(k);
% end

controldemod = wifidec(wifideint(wifidemod(...
    reshape(recfreqnorm(13:15, :), [], 1), 1)), 1);
[tt1, m1, r1, len1] = GetHeader(controldemod);
% tt1 = 1;
    
convrecdemod = wifidec(wifideint(wifidemod(...
    reshape(recfreqnorm(16:end-50, :), [], 1), tt1)), tt1);

realdataenc = reshape(wifimod(wifideint(wifienc(frame1.Data, tt1)), tt1), [], 48);

comrecfreqnorm = recfreqnorm(16:end-50, :);
err = mean(abs(realdataenc(1:end, :) - comrecfreqnorm(1:end, :)), 2);


% plot(1:length(convrecdemod)-34, err(1:100));
plot(1:length(err), err);

%% Received Signal Deconstruction
% [dataserdemod1 dataserdemod2] = ...
%     wifidemodjmd(reshape(recfreq(14:end, :), [], 1), ...
%     frame1.TransmissionRate, frame2.TransmissionRate, g1, g2);

% dataserdemod1 = convrecdemod;
% dataserdemod1 = reshape(frame1.DataMod, [], 1);
ber = step(errcalc, frame1.Data(1:end-34), convrecdemod(1:end-34))
% ser(n) = ber(1);
% reset(errcalc12);
% 
% plot(snrs, ser);
% 
% set(gca,'yscale','log');
% ylabel('SER');
% xlabel('SNR');
% grid on;
% title('QPSK 1/2 802.11g SNR vs SER');
% % legend('SIR = -5', 'SIR = -3', 'SIR = -1', 'SIR = 1', 'SIR = 3', 'SIR = 5');
% hold off;


% data1 = randi([0 1], 50, 48);
% data2 = randi([0 1], 50, 48);
% data3 = step(ofdmenc, [data1; zeros(50, 48)]) ...
%     + step(ofdmenc, [zeros(50, 48); data2]);
% data4 = step(ofdmdec, data3);
% datacomp1 = reshape(data1, [], 1);
% datacomp41 = reshape(data4(1:50, :), [], 1);
% sum(abs(datacomp1-datacomp41))
