%% Andrew Apollonsky

close all;
clc;
clear all;

%% Simulation Parameters
len = 500;

%% Signal 1 Parameters/Generation
m1 = 4;
len1 = 200;
slen1 = 10*log10(m1);

ebno1 = [15];
snr1 = ebno1 + 3 + 10*log10(slen1);

mod1 = comm.PSKModulator(m1);
mod1.PhaseOffset = 0;
demod1 = comm.PSKDemodulator(m1);
demod1.PhaseOffset = 0;
errcalc1 = comm.ErrorRate;

% Signal Generation & Modulation
x1 = randi([0 m1-1], len1, 1);
y1 = step(mod1, x1);


%% Channel Parameters/Processing
chan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
chan1.SignalPower = (real(y1)'*real(y1))/length(real(y1));

dist1 = [.8 .2 - .14i .1 + .08i];
% dist1 = [1];

chan1.SNR = snr1;
w1 = step(chan1, filter(dist1, 1, y1));

% Phase Rotation
for t = 1:len1
    w1(t) = w1(t) * exp(t*i/50);
end
%% LMS Equalization 
% Pilot
step1 = .2;
taplen1 = 6;
e1 = zeros(1, len1);
tapin1 = cell(1, len1);
tapest1 = cell(1, len1);
w1eq = zeros(len1, 1);
% tapest1 = deal(zeros(plen1, 1));
% tapest1{1} = [0];

mvavgerlen = 5;
mvavger = 99;
isconv = 0;

for t = 1:len1
    if t > taplen1
        tapin1(t) = {fliplr(w1(t-taplen1 + 1:t))}; % If time > WN, use WN
    else
        tapin1(t) = {fliplr(w1(1:t))}; % If time < weight numbers, use smaller #
            tapest1{t} = [tapest1{t}; 0]; % If time < WN, increase size of next WN
    end
    e1(t) = y1(t) -  tapest1{t}' * tapin1{t};
    tapest1{t+1} = tapest1{t} + step1*tapin1{t} * e1(t)';

    %Approx avg for convergence in polynomial for QPSK and ebno < 30;
    errorthresh = .3;
    if t > mvavgerlen
        mvavger = mean(abs(e1(t-mvavgerlen:t)));
    end
    if mvavger < errorthresh && isconv == 0
        isconv = t;
    end
    

    w1eq(t) = sum(tapest1{t}' * w1(t));
end

%% Demodulation
z1 = step(demod1, w1eq(50:end));
ber = step(errcalc1, z1, x1(50:end))

% plot([1:plen1], ep1);
subplot(221),plot([1:len1],y1),ylabel('Desired Signal'),
subplot(222),plot([1:len1],w1),ylabel('Input Signal+Noise'),axis([0 len1 -1 1]),
subplot(223),plot([1:len1],abs(e1)),ylabel('Error'),axis([0 len1 -1 1]),
subplot(224),plot([1:len1],w1eq),ylabel('Adaptive Desired output'),
axis([0 len1 -1 1]);

% h = scatterplot(yp1, 1, plen1, 'bx');
% scatterplot(yd, 1, plen1, 'k*', h);


