%% Andrew Apollonsky
close all force;
clc;
clear all;

%% Simulation Parameters
len = 300;
m1 = 4;
m2 = 4;

ebno = 40;
sir = 0;
rS = 2e7;
maxDopp = 500;
rotFac = 1/100;
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
x1 = randi([0 m1-1], len, 1);
x2 = randi([0 m2-1], len, 1);
y1 = step(mod1, x1);
y2 = step(mod2, x2);
chan3.SignalPower = (real(y1*g1)'*real(y1*g1))/length(y1);

[dataRay1, pathG1] = step(chan1, y1);
[rx2, pathG2] = step(chan2, y2);
for t = 1:len
    rx2(t) = rx2(t) * exp(t*1i*rotFac);
    pathG2(t) = pathG2(t) * exp(t * 1i * rotFac);
end

rx1 = step(chan3, dataRay1);
rxSig = g1*rx1 + g2*rx2;

% [out1, out2, wts1, wts2] = ItJDec1(rxSig, pathG1(1)*g1, pathG2(1)*g2, ...
%     mod1, mod2, demod1, demod2);

[out1, out2, wts1, wts2] = JMDec1(rxSig, pathG1(1)*g1, pathG2(1)*g2, ...
    mod1, mod2);

% if isreal(out1)
%     out1demod = step(demod1, complex(out1));
% else   
%     out1demod = step(demod1, out1); 
% end
% 
% if isreal(out2)
%     out2demod = step(demod2, complex(out2));
% else   
%     out2demod = step(demod2, out2); 
% end

out1demod = out1;
out2demod = out2;
errsmooth1 = smooth(abs(x1 - out1demod), len/5);
errsmooth2 = smooth(abs(x2 - out2demod), len/5);

sert1 = step(errcalc, out1demod(end/2:end), x1(end/2:end));
reset(errcalc);
sert2 = step(errcalc, out2demod(end/2:end), x2(end/2:end));
reset(errcalc);
ser = [sert1(1) sert2(1)]
%% Plot

subplot(321), hold all;
plot(1:len, abs(wts1(1:len))), 
plot(1:len, abs(pathG1 .* g1)), 
title('Filter 1 Weight Magnitude'); 
xlabel('Time'); ylabel('Magnitude')
hold off;

subplot(323), hold all,
plot(1:len, rad2deg(angle(wts1(1:len)))),
plot(1:len, rad2deg(angle(pathG1))), 
axis([0 len -180 180]);
title('Filter 1 Weight Angle'),
xlabel('Time'); ylabel('Angle (^\circ)');
set(gca,'YTick', -180:90:180);
hold off;

subplot(325), plot(1:len, abs(errsmooth1)), title('Error 1');
xlabel('Time'); ylabel('Error');
set(gca ,'yscale','log');

subplot(322), hold all;
plot(1:len, abs(wts2(1:len))), 
plot(1:len, abs(pathG2 .* g2)), 
title('Filter 2 Weight Magnitude'); 
xlabel('Time'); ylabel('Magnitude')
hold off;

subplot(324), hold all,
plot(1:len, rad2deg(angle(wts2(1:len)))),
plot(1:len, rad2deg(angle(pathG2))), 
axis([0 len -180 180]);
title('Filter 2 Weight Angle'),
xlabel('Time'); ylabel('Angle (^\circ)');
set(gca,'YTick', -180:90:180);
hold off;

subplot(326), plot(1:len, abs(errsmooth2)), title('Error 2');
xlabel('Time'); ylabel('Error');
set(gca ,'yscale','log');
