%% Andrew Apollonsky
close all force;
clc;
clear all;

%% Simulation Parameters
len = 10;
m1 = 4;
m2 = 4;

ebno = 8;
sir = 2;
rS = 2e7;
maxDopp = .8;
rotFac = 1/12000;
%% Object Initialization
bpskmod = comm.BPSKModulator(0);
bpskdemod = comm.BPSKDemodulator(0);
qam4mod = comm.RectangularQAMModulator(4);
qam4demod = comm.RectangularQAMDemodulator(4);
qam16mod = comm.RectangularQAMModulator(16);
qam16demod = comm.RectangularQAMDemodulator(16);
qam64mod = comm.RectangularQAMModulator(64);
qam64demod = comm.RectangularQAMDemodulator(64);

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
elseif m1 == 64
    mod2 = qam64mod;
    demod2 = qam64demod;
end

errcalc = comm.ErrorRate;

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
g1 = 1;
g2 = 10^(-sir/20);

%% Data Creation
x1 = randi([0 m1-1], len, 1);
x2 = randi([0 m2-1], len, 1);
y1 = step(mod1, x1);
y2 = step(mod2, x2);
chan3.SignalPower = (real(y1*g1)'*real(y1*g1))/length(y1);

[dataRay1, pathG1] = step(chan1, y1);
[rx2, pathG2] = step(chan2, y2);
% pathG1 = 1;
% pathG2 = ones(len, 1);
% dataRay1 = y1;
% rx2 = y2;
% for t = 1:len
%     rx2(t) = rx2(t) * exp(t*1i*rotFac);
%     pathG2(t) = pathG2(t) * exp(t * 1i * rotFac);
% end

rx1 = step(chan3, dataRay1);
rxSig = g1*rx1 + g2*rx2;

% [out1, out2, err1, err2] = JMDDemod2(rxSig, y1, y2, mod1, mod2, m1, m2);
% [out1, out2, err1, err2, wts1, wts2] = SoftDecisionML1(rxSig, y1, y2, mod1, mod2, m1, m2);
% [out1, out2, err1, err2, wts1, wts2] = SoftDecisionZF(rxSig, y1, y2, demod1, demod2);
% [out1, out2, err1, err2, wts1, wts2] = SoftDecisionML2(rxSig, y1, y2, mod1, mod2, m1, m2);
[out1, out2, err1, wts1, wts2] = HeavyML1(rxSig, y1, y2, mod1, mod2, m1, m2);
out1 = step(demod1, out1); out2 = step(demod2, out2); 
% errsmooth = smooth(abs(err1) + abs(err2), 31);
errsmooth = err1;
% dec1 = step(demod1, out);
sert1 = step(errcalc, out1(end/2:end), x1(end/2:end));
reset(errcalc);
sert2 = step(errcalc, out2(end/2:end), x2(end/2:end));
reset(errcalc);
ser = [sert1(1) sert2(1)]
%% Plot


subplot(231), hold all, 
plot(1:len, abs(wts1(1:len))), 
plot(1:len, abs(pathG1 .* g1)), 
title('Filter 1 Weight Magnitude');
% legend('Adaptive Weight', 'Real Weight');
% set(gca ,'yscale','log'), 
hold off;

subplot(234), hold all,
plot(1:len, angle(wts1(1:len))),
plot(1:len, angle(pathG1)), 
title('Filter 1 Weight Angle')
% legend('Adaptive Weight', 'Real Weight');
hold off;

subplot(232), hold all;
plot(1:len, abs(wts2(1:len))), 
plot(1:len, abs(pathG2 .* g2)), 
title('Filter 2 Weight Magnitude');
% legend('Adaptive Weight', 'Real Weight');
% set(gca ,'yscale','log'), 
hold off;

subplot(235), hold all,
plot(1:len, angle(wts2(1:len))),
plot(1:len, angle(pathG2)), 
title('Filter 2 Weight Angle'),
% legend('Adaptive Weight', 'Real Weight');
hold off;

subplot(236), plot(1:len, abs(errsmooth)), title('Error');
set(gca ,'yscale','log'),  