%% Andrew Apollonsky

close all;
clc;
clear all;

%% Simulation Parameters
len = 500;
numsigs = 2;
sigbd = zeros(numsigs, 2);

%% Signal 1 Parameters/Generation
m1 = 2;
delay1 = 0;
len1 = 300;
slen1 = 10*log10(m1);

ebno1 = [15];
snr1 = ebno1 + 3 + 10*log10(slen1);

mod1 = comm.PSKModulator(m1);
mod1.PhaseOffset = pi/2;
demod1 = comm.PSKDemodulator(m1);
demod1.PhaseOffset = 0;
errcalc1 = comm.ErrorRate;

g1 = 1; % Gainnonz = find(y1);

% Signal Generation
x1 = [randi([0 m1-1], len1, 1)];
y1 = [zeros(delay1, 1); step(mod1, x1); zeros(len-len1-delay1, 1)];

nonz1 = find(y1);
sigbd(1, :) = [nonz1(1) nonz1(end)]; % Simulates rec knowing sig start/end t

% Channel Effects
chan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
chan1.SignalPower = (real(y1*g1)'*real(y1*g1))/length(real(y1*g1));
dist1 = [.8 .2 - .14i .1 + .08i]; % Signal 1 distortion
chan1.SNR = snr1;
w1 = step(chan1, filter(dist1, 1, y1));

% Phase Rotation
for t = 1:len1
    w1(t) = w1(t) * exp(t*i/70);
end
%% Signal 2 Parameters/Generation
s12r = 2; % Like SIR, ratio of power of signal 1 to power of signal 2
g2 = 10^(-s12r/20);

m2 = 2;
delay2 = 300;
len2 = 200;
slen2 = 10*log10(m2);

ebno2 = [10];
snr2 = ebno2 + 3 + 10*log10(slen2);

mod2 = comm.PSKModulator(m2);
mod2.PhaseOffset = 0;
demod2 = comm.PSKDemodulator(m2);
demod2.PhaseOffset = 0;
errcalc2 = comm.ErrorRate;

% Signal Generation
x2 = [randi([0 m2-1], len2, 1)];
y2 = [zeros(delay2, 1); step(mod2, x2); zeros(len -len2 - delay2, 1)];

nonz2 = find(y2);
sigbd(2, :) = [nonz2(1) nonz2(end)];

% Channel Effects
chan2 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
chan2.SignalPower = (real(y2*g2)'*real(y2*g2))/length(real(y2*g2));
dist2 = [.9 + .1i .1+.1i .05 - .09i]; % Signal 2 distortion
chan2.SNR = snr2;
w2 = step(chan2, filter(dist2, 1, y2));

% Phase Rotation
for t = 1:len1
    w2(t) = w2(t) * exp(t*i/90);
end

%% Receiver
w = w1 + w2;

%  LMS Equalization
stepg = .2;

taplen = 6;
err = zeros(1, len);

tapin = cell(1, len);
tapest = cell(1, len);

tapin1 = cell(1, len);
tapin2 = cell(1, len);

weq = zeros(len, 1);
mvavgerlen = 5;
mvavger = 100;
isconv = 0;
numsigs = 1;
minimat = zeros(m1, m2);
z1 = zeros(len, 1);
z2 = zeros(len, 1);
for t = 1:len
    if numsigs == 1  % Beginning: One signal, no convergence
        if t > taplen
            tapin(t) = {fliplr(w(t-taplen + 1:t))}; % If time > WN, use WN
        else
            tapin(t) = {fliplr(w(1:t))}; % If time < weight numbers, use smaller #
            tapest{t} = [tapest{t}; 0]; % If time < WN, increase size of next WN
        end
        err(t) = y1(t) -  tapest{t}' * tapin{t};
        tapest{t+1} = tapest{t} + stepg*tapin{t} * err(t)';
        errthresh = .3;
        if t > mvavgerlen
            mvavger = mean(abs(err(t-mvavgerlen:t)));
        end
        if mvavger < errthresh && isconv == 0
            isconv = t;
        end
        
        weq(t) = sum(tapest{t}' * w(t)); % Equalized signal
        if isreal(weq(t))
            z1(t) = step(demod1, complex(weq(t)));
        else
            z1(t) = step(demod1, weq(t));
        end
        
        z2(t) = -1;
    end
    

        
    if mvavger > errthresh && isconv > 1
        numsigs = 2;
        isconv = 0;
    end
    if numsigs == 2 && isconv == 0
        tapin1(t) = {fliplr(step(mod1, z1(t-taplen:t-1).'))};
        temp = length(find(z2 + 1));
        if temp > taplen
            tapin2(t) = {fliplr(step(mod2, z2(t-taplen:t-1)))};
        else
            tapin2(t) = {fliplr(step(mod2, [z2(t-temp:t-1); zeros(taplen-temp, 1)]))};
        end
        err(t) = y1(t) + y2(t) - tapest{t}'*tapin1{t} - tapest{t}'*tapin2{t};
        tapest{t+1} = tapest{t} + stepg*tapin1{t} * err(t)' + stepg*tapin2{t} * err(t)';
              
        weq(t) = sum(tapest{t}' * w(t));
        
        for q = 0:(m1-1) % For every possible y1
            for a = 0:(m2-1) %For every possible y2
                minimat(q+1, a+1) = ...
                    abs(weq(t) - g1*step(mod1, q) - g2*step(mod2, a));
            end
        end
        [r, c, v] = find(minimat == min(min(minimat))); % Find minimum
        
        z1(t) = r(1)-1;
        z2(t) = c(1)-1;
    end
end
        
        
        
% %% Demodulation
% z1 = step(demod1, w1eq(50:end));
% ber = step(errcalc1, z1, x1(50:end))
% 
% plot([1:plen1], ep1);
subplot(221),plot([1:len],y1),ylabel('Desired Signal'),
subplot(222),plot([1:len],w),ylabel('Input Signal+Noise'),axis([0 len1 -1 1]),
subplot(223),plot([1:len],abs(err)),ylabel('Error'),axis([0 len1 -1 1]),
subplot(224),plot([1:len],weq),ylabel('Adaptive Desired output'),
axis([0 len1 -1 1]);
% 
% % h = scatterplot(yp1, 1, plen1, 'bx');
% % scatterplot(yd, 1, plen1, 'k*', h);


