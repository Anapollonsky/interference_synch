%% Andrew Apollonsky

close all;
clc;
clear all;

%% Simulation Parameters
len = 2000;

%% Signal 1 Parameters/Generation
m1 = 2;
delay1 = 0;
len1 = 2000;
slen1 = 10*log10(m1);

ebno1 = [15];
snr1 = ebno1 + 3 + 10*log10(slen1);

mod1 = comm.PSKModulator(m1);
mod1.PhaseOffset = pi/4;
demod1 = comm.PSKDemodulator(m1);
demod1.PhaseOffset = pi/4;
errcalc1 = comm.ErrorRate;

g1 = 1; % Gainnonz = find(y1);

% Signal Generation
x1 = randi([0 m1-1], len1, 1);
y1 = [zeros(delay1, 1); step(mod1, x1); zeros(len-len1-delay1, 1)];

% Channel Effects
chan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
chan1.SignalPower = (real(y1*g1)'*real(y1*g1))/length(y1);
% dist1 = [.8 .2 - .14i .1 + .08i]; % Signal 1 distortion
dist1 = 1;
chan1.SNR = snr1;
w1 = step(chan1, filter(dist1, 1, y1));

% Phase Rotation
for t = 1:len1
    w1(t) = w1(t) * exp(t*1i/500);
end
%% Signal 2 Parameters/Generation
s12r = 0; % Like SIR, ratio of power of signal 1 to power of signal 2
g2 = 10^(-s12r/20);

m2 = 2;
delay2 = 300;
len2 = 1700;
slen2 = 10*log10(m2);

ebno2 = [15];
snr2 = ebno2 + 3 + 10*log10(slen2);

mod2 = comm.PSKModulator(m2);
mod2.PhaseOffset = -pi/4;
demod2 = comm.PSKDemodulator(m2);
demod2.PhaseOffset = -pi/4;
errcalc2 = comm.ErrorRate;

% Signal Generation
x2 = randi([0 m2-1], len2, 1);
y2 = [zeros(delay2, 1); step(mod2, x2); zeros(len -len2 - delay2, 1)];

% Channel Effects
chan2 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
chan2.SignalPower = (real(y2*g2)'*real(y2*g2))/length(y2);
% dist2 = [.9 + .1i .1+.1i .05 - .09i]; % Signal 2 distortion
dist2 = 1;
chan2.SNR = snr2;
w2 = step(chan2, filter(dist2, 1, y2));

% Phase Rotation
for t = 1:len
    w2(t) = w2(t) * exp(t*1i/300);
end

%% Receiver
w = w1 + w2;

%  LMS Equalization
lstep = .1;             % Step of LMS filter. 
taplen = 12;             % Length of LMS filter.
err = zeros(1, len);        % Error vector for LMS.
err1 = zeros(1, len);        % Error vector for LMS.
err2 = zeros(1, len);        % Error vector for LMS.
tapin1 = zeros(taplen, len);   % Undemodulated signal received for LMS
tapin2 = zeros(taplen, len);
tapest = zeros(taplen, len);   % Tap weight vector
w1eq = zeros(len, 1);   % Undemodulated estimate for signal1
errmvvarlen = 20;       % Declares length of errmvvar
errmvvar = ones(len, 1);   % Keeps track of error variance over time.
numsigs = ones(len, 1);     % Number of signals detected
convthreshup = .1;   % Variance Threshold to declare that de-converged
convthreshdown = .02;
isconv = zeros(len, 1); % Keeps track of convergence;
minimat = zeros(m1, m2, len);
z1 = zeros(len, 1);
z2 = zeros(len, 1);
weq = zeros(len, 1);

for t = 1:len
    if t > 1
        numsigs(t) = numsigs(t-1);
        isconv(t) = isconv(t-1);
    end
    % Keep track of whether converged, detected # of signals.
    if isconv(t) == 0 && t > 5 && errmvvar(t-1) < convthreshdown
        isconv(t) = 1;
    end
    
    if isconv(t) == 1 && errmvvar(t-1) > convthreshup
        isconv(t) = 0;
        numsigs(t) = numsigs(t-1) + 1;
    end
    
    % If one signal, simple LMS detection.
    if numsigs(t) == 1
        if t >= taplen
            tapin1(:, t) = flipud(w(t-taplen+1:t));
        else
            tapin1(:, t) = [flipud(w(1:t)); zeros(taplen-t, 1)];
        end
        err(t) = y1(t) - tapest(:, t)' * tapin1(:, t);
        tapest(:, t+1) = tapest(:, t) + lstep * tapin1(:, t) * err(t)';
        w1eq(t) = sum(tapest(:, t)' * w(t));
        if isreal(w1eq(t))
            z1(t) = step(demod1, complex(w1eq(t)));
        else
            z1(t) = step(demod1, w1eq(t));
        end
    end
    
    % If two signals, LMS + MD
    if numsigs(t) > 1
        switch12temp = find(numsigs == 2);
        switch12 = switch12temp(1);
        
%         if t < (switch12 + taplen)
            tapin1(:, t) = flipud(step(mod1, z1(t-taplen:t-1)));
%             tapin1(:, t) = flipud(w(t-taplen+1:t));
            tapin2(:, t) = flipud(step(mod2, z2(t-taplen:t-1)));
%         else
%             for q = 0:(m1-1) % y1
%                 tapin1(:, t) = tapin1(:, t) + ...
%                     flipud(squeeze(sum(minimat(q+1, :, t-taplen:t-1), 2) * step(mod1, q)));
%                 angle(:, t) = 
%             end
%             
%             for a = 0:(m2-1) %y2
%                   tapin2(:, t) = tapin2(:, t) + ...
%                     flipud(squeeze(sum(minimat(:, a+1, t-taplen:t-1), 1) * step(mod2, a)));  
%             end
%         end
%         err(t) = 1*(y1(t-1) + y2(t-1)) - tapest(:, t)'*tapin1(:, t) - tapest(:, t)'*tapin2(:, t);
%         tapest(:, t+1) = tapest(:, t) + lstep * err(t)' * (tapin1(:, t)+tapin2(:, t));
        err1(t) = y1(t-1) - tapest(:, t)'*tapin1(:, t);
%         err2(t) = y2(t-1) - tapest(:, t)'*tapin2(:, t);
%         tapest(:, t+1) = tapest(:, t) + lstep * err1(t)' * tapin1(:, t) + lstep * err2(t)' * tapin2(:, t);
        err(t) = err1(t);% + err2(t);
        tapest(:, t+1) = tapest(:, t) + lstep * err(t)' * (tapin1(:, t));
        weq(t) = sum(tapest(:, t)' * w(t));
        for q = 0:(m1-1) % For every possible y1
            for a = 0:(m2-1) %For every possible y2
                minimat(q+1, a+1, t) = ...
                    abs(weq(t) - g1*step(mod1, q) - g2*step(mod2, a));
            end
        end
        [r, c, v] = find(minimat(:, :, t) == min(min(minimat(:, :, t)))); % Find minimum
        minimat(:, :, t) = 1./(minimat(:, :, t).^2);
        minimat(:, :, t) = minimat(:, :, t)/sum(sum(minimat(:, :, t)));
        
        z1(t) = r(1)-1;
        z2(t) = c(1)-1;
    end
    
    if t <= errmvvarlen
        errmvvar(t) = var(abs(err(1:t)));
    else
        errmvvar(t) = var(abs(err(t-errmvvarlen:t)));
    end

end

w1eq = step(mod1, z1);
w2eq = step(mod2, z2);
z1 = z1(delay1 + 1:delay1 + len1);
z2 = z2(delay2 + 1:delay2 + len2);

%% Error Calculation 
ber1 = step(errcalc1, z1(500:1000), x1(500:1000))
ber2 = step(errcalc2, z2(100:500), x2(100:500))

%% Plot
subplot(331),plot(1:len,y1),title('Desired Signal 1');
subplot(332),plot(1:len,w1),title('Distorted Signal 1'),axis([0 len1 -1 1]),
subplot(333),plot(1:len,w1eq),title('Detected Signal 1'),axis([0 len1 -1 1]),
subplot(334),plot(1:len,y2),title('Desired Signal 2');
subplot(335),plot(1:len,w2),title('Distorted Signal 2'),axis([0 len1 -1 1]),
subplot(336),plot(1:len,w2eq),title('Detected Signal 2'),axis([0 len1 -1 1]),
subplot(337),plot(1:len,abs(err)),title('LMS Error'),%axis([0 len1 0 5]),
subplot(338),plot(1:len, errmvvar),title('Local Error Variance');



