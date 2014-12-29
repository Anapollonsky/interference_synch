%% Andrew Apollonsky
close all force;
clc;
clear all;

%% Simulation Parameters
len = 12;
numofdm = 1;

m1 = 4;
m2 = 4;
mu = 5; % Cyclic Prefix

ebno = 50;
sir = 0;
rS = 2e7;
maxDopp = 50;
rotFac = 1;
%% Object Initialization
bpskmod = comm.BPSKModulator(0);
bpskdemod = comm.BPSKDemodulator(0);
qam4mod = comm.RectangularQAMModulator(4);
qam4demod = comm.RectangularQAMDemodulator(4);
qam16mod = comm.RectangularQAMModulator(16);
qam16demod = comm.RectangularQAMDemodulator(16);
qam64mod = comm.RectangularQAMModulator(64);
qam64demod = comm.RectangularQAMDemodulator(64);
qam256mod = comm.RectangularQAMModulator(256);
qam256demod = comm.RectangularQAMDemodulator(256);

if m1 == 2
    mod1 = bpskmod;
    demod1 = bpskdemod;
elseif m1 == 4
    mod1 = qam4mod;
    demod1 = qam4demod;
elseif m1 == 16
    mod1 = qam16mod;
    demod1 = qam16demod;
elseif m1 == 64
    mod1 = qam64mod;
    demod1 = qam64demod;
elseif m1 == 256
    mod1 = qam256mod;
    demod1 = qam256demod;
end

if m2 == 2
    mod2 = bpskmod;
    demod2 = bpskdemod;
elseif m2 == 4
    mod2 = qam4mod;
    demod2 = qam4demod;
elseif m2 == 16
    mod2 = qam16mod;
    demod2 = qam16demod;
elseif m2 == 64
    mod2 = qam64mod;
    demod2 = qam64demod;
elseif m2 == 256
    mod2 = qam256mod;
    demod2 = qam256demod;
end

errcalc = comm.ErrorRate;

% OFDM
ofdmenc = OFDMPre;
ofdmdec = OFDMPost;

% Channels
chan1 = comm.MIMOChannel(...
    'SampleRate', rS, ...
    'MaximumDopplerShift', maxDopp, ...
    'NumTransmitAntennas', 1, ...
    'NumReceiveAntennas',  1, ...
    'PathGainsOutputPort', true, ...
    'TransmitCorrelationMatrix', 1, ...
    'ReceiveCorrelationMatrix', 1);

chan2 = comm.MIMOChannel(...
    'SampleRate', rS, ...
    'MaximumDopplerShift', maxDopp, ...
    'NumTransmitAntennas', 1, ...
    'NumReceiveAntennas',  1, ...
    'PathGainsOutputPort', true, ...
    'TransmitCorrelationMatrix', 1, ...
    'ReceiveCorrelationMatrix', 1);
    
chan3 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');

chan3.SNR = ebno + 10*log10(log2(m1));
sir2 = sir - 10*log10(log2(m1)) + 10*log10(log2(m2));
g1 = 1;
g2 = 10^(-sir2/20);

%% Data Creation
x1 = randi([0 m1-1], len * numofdm, 1);
x2 = randi([0 m2-1], len * numofdm, 1);
y1 = step(mod1, x1);
y2 = step(mod2, x2);
y1par = reshape(y1, [], numofdm);
y2par = reshape(y2, [], numofdm);
y1ofdm = step(ofdmenc, y1par, mu);
y2ofdm = step(ofdmenc, y2par, mu);

%% Channel Processing
chan3.SignalPower = (real(y1*g1)'*real(y1*g1))/length(y1);
[rx1, pG1] = OFDMChan(y1ofdm, chan1);
[rx2, pG2] = OFDMChan(y2ofdm, chan2);

for t = 1:length(rx2)
    rx2(t) = rx2(t) * exp(t*1i*rotFac);
    pG2(t) = pG2(t) * exp(t * 1i * rotFac);
end

for k = 1:numofdm
    rx1(:, k) = step(chan3, rx1(:, k));
end

pG1 = pG1(mu+1:end, :);
pG2 = pG2(mu+1:end, :);
rxSigOFDM = g1*rx1 + g2*rx2;
rxSig = step(ofdmdec, rxSigOFDM, mu);
%% Decoding
% rxSigSer = reshape(rxSig, [], 1);
% pG1ser = reshape(pG1, [], 1);

[out1, wts2, iters] = gen_al1(rxSig, reshape(y2, len, []) , mod1, pG1);

if isreal(out1)
    out1demod = step(demod1, complex(out1));
else
    out1demod = step(demod1, reshape(out1, [], 1)); 
end

out1demodbin = reshape(de2bi(out1demod, log2(m1)), [], 1);
x1bin = reshape(de2bi(x1, log2(m1)), [], 1);

sert1 = step(errcalc, x1bin, out1demodbin);
reset(errcalc);
iters
BER = [sert1(1)]
% rotfac
%% Plot

subplot(211), hold all;
plot(1:len, abs(wts2(:, 1))), 
plot(1:len, abs(pG2(:, 1) .* g2)), 
title('Filter 2 Weight Magnitude');
legend('Estimated', 'Real'); 
xlabel('Time'); ylabel('Magnitude')
hold off;

subplot(212), hold all,
plot(1:len, rad2deg(angle(wts2(:, 1)))),
plot(1:len, rad2deg(angle(pG2(:, 1)))), 
axis([0 len -180 180]);
title('Filter 2 Weight Angle'),
legend('Estimated', 'Real');
xlabel('Time'); ylabel('Angle (^\circ)');
set(gca,'YTick', -180:90:180);
hold off;
% 
% subplot(313), plot(1:len, abs(errsmooth)), title('Error');
% xlabel('Time'); ylabel('Error');
% set(gca ,'yscale','log');
