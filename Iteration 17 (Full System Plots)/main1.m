%% Andrew Apollonsky
% This script generates plots for random 
close all force;
clc;
clear all;

%% Simulation Parameters
len1 = 12; % Training length
len2 = 60; % Frame length (discounting training)
m1 = 4;
m2 = 4;

ebnos = [0:4:20];
sirs = [10];
numiter = 350;

ts = 4e-6; %x/s
freq = 2.4; %(GHz)

maxDopp = 1;
rotfacreal = .05;
% 
% temp = (rotfacguess - rotfacreal)*180/pi;
% freqerr = (temp / ts) / 360 / 1e3; 
% freqerr2 = mean(abs(freqerr), 3);

realfreqoffset = rotfacreal /(2*pi) / ts / 1e3; %KHz

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

g1 = 1;
bersgen = nan(length(sirs), length(ebnos), numiter);
berstand = nan(length(sirs), length(ebnos), numiter);
pgsguess = nan(length(sirs), length(ebnos), numiter);
rotfacguess = nan(length(sirs), length(ebnos), numiter);
pgsreal = nan(length(sirs), length(ebnos), numiter);
numiters = nan(length(sirs), length(ebnos), numiter);
delayerr = nan(length(sirs), length(ebnos), numiter);
out1ber = nan(length(sirs), length(ebnos), numiter);
out2ber = nan(length(sirs), length(ebnos), numiter);
out3ber = nan(length(sirs), length(ebnos), numiter);
sigpresentest = nan(length(sirs), length(ebnos), numiter);
sigpresentreal = nan(length(sirs), length(ebnos), numiter);


wB = waitbar(0,'0%');
prog = 0;
%% Calculations
[evm1, evm2] = errcalc(sirs, ebnos + 10*log10(log2(m1)), 1e5, mod1, demod1, bpskmod);
tic;
for siriter = 1:length(sirs)
    sir2 = sirs(siriter) - 10*log10(log2(m1)) + 10*log10(log2(m2));
    g2 = 10^(-sir2/20);
    for ebnoiter = 1:length(ebnos)
        SNR = ebnos(ebnoiter) + 10*log10(log2(m1));
        for iter = 1:numiter
            sig2delay = randi([1 len2-len1 - 1]);
            totlen = sig2delay + len1 + len2;
            
            x1 = [randi([0 m1-1], len2, 1); zeros(sig2delay + len1, 1)];
            x2 = [zeros(sig2delay, 1); randi([0 1], len1, 1);...
                randi([0 m2-1], len2, 1)]; 
            
            y1 = [step(mod1, x1(1:len2)); zeros(sig2delay + len1, 1)];
            y2 = [zeros(sig2delay, 1); ...
                step(bpskmod, x2(sig2delay+1:sig2delay+len1)); ...
                step(mod2, x2(sig2delay+len1+1:end))];
            
            pathG1 = [ones(len2, 1) * exp(rand(1) * 2*pi*1i) * g1; zeros(len1 + sig2delay, 1)];
            pathG2 = [zeros(sig2delay, 1); ones(len1 + len2, 1) * exp(rand(1) * 2*pi*1i) * g2];
            
            for t = sig2delay+1:totlen
                pathG2(t) = pathG2(t) * exp((t-sig2delay-1) * 1i * rotfacreal);
            end
            
            % Apply path gains to data
            rx1 = pathG1 .* y1;
            rx2 = pathG2 .* y2;

            % Add AWGN
            rx1 = awgn(rx1, SNR, 'measured');

            % Sum signals
            if rand(1) > .5
                rxSig = rx1 + rx2;
                sigpresreal = 1;
            else
                rxSig = rx1;
                sigpresreal = 0;
            end
            
            % Find if signal present, and if so, delay
            [delayest, sigpres] = delayfind(rxSig, mod1,...
                demod1, pathG1, len2, evm1(siriter, ebnoiter), evm2(siriter, ebnoiter));
            
            if sigpres == 1 && sigpresreal == 1
                if sig2delay == delayest
                    % Synchronize
                    [out1gen, wts2, iters, pgbase, rotfac] = ...
                        gen_al4(rxSig(delayest + 1:delayest + len1), ...
                        y2(delayest+1:delayest + len1), mod1, ...
                        pathG1(delayest+1:delayest + len1));

                    % Combine data to estimate future signal 2 path gain for
                    % non-training data
                    pg2 = pgbase * exp(1i * rotfac * [0:len1 + len2-1]).';

                    % Joint Detection
                    [out1, out2] = JMD(rxSig(delayest+len1+1:len2),...
                        mod1, mod2, pathG1(delayest+len1+1:len2), pg2(len1+1:len2-delayest));

                    % Finish Detecting 2
                    out3demod = step(demod2, rxSig(len2+1:end)./ pg2(len2-delayest+1:end));

                    % Demodulate synchronization things
                    out1demodgen = step(demod1, out1gen); 
                    out1demodstand = step(demod1, rxSig(delayest+1:delayest+len1) ./ pathG1(delayest+1:delayest+len1));

                    % Convert to binary
                    out1demodgenbin = reshape(de2bi(out1demodgen, log2(m1)), [], 1);
                    out1demodstandbin = reshape(de2bi(out1demodstand, log2(m1)), [], 1);
                    x1binsynch = reshape(de2bi(x1(sig2delay+1:sig2delay+len1), log2(m1)), [], 1);

                    out1bin = reshape(de2bi(out1, log2(m1)), [], 1);
                    out2bin = reshape(de2bi(out2, log2(m2)), [], 1);
                    out3bin = reshape(de2bi(out3demod, log2(m2)), [], 1);

                    jmdx1b = reshape(de2bi(x1(sig2delay + len1 + 1:len2), log2(m1)), [], 1);
                    jmdx2b = reshape(de2bi(x2(sig2delay + len1 + 1:len2), log2(m2)), [], 1);
                    normx2endb = reshape(de2bi(x2(len2+1:end), log2(m2)), [], 1);

                    % Find BER
                    synchgenber = sum(out1demodgenbin ~= x1binsynch)/length(x1binsynch);
                    synchstandber = sum(out1demodstandbin ~= x1binsynch)/length(x1binsynch);
                    out1ber1 = sum(out1bin ~= jmdx1b)/length(jmdx1b);
                    out2ber1 = sum(out2bin ~= jmdx2b)/length(jmdx2b);
                    out3ber1 = sum(out3bin ~= normx2endb)/length(out3bin);
                end
            end
                        
            %% Assign to Global
     
            sigpresentreal(siriter, ebnoiter, iter) = sigpresreal;
            sigpresentest(siriter, ebnoiter, iter) = sigpres;
            
            if (sigpres == 1 && sigpresreal == 1)
                delayerr(siriter, ebnoiter, iter) = (sig2delay ~= delayest);
                if delayerr(siriter, ebnoiter, iter) == 0           
                    bersgen(siriter, ebnoiter, iter) = synchgenber;
                    berstand(siriter, ebnoiter, iter) = synchstandber;
                    pgsguess(siriter, ebnoiter, iter) = pgbase;
                    rotfacguess(siriter, ebnoiter, iter) = rotfac;
                    pgsreal(siriter, ebnoiter, iter) = pathG2(sig2delay+1);
                    numiters(siriter, ebnoiter, iter) = iters;
                    out1ber(siriter, ebnoiter, iter) = out1ber1;
                    out2ber(siriter, ebnoiter, iter) = out2ber1;
                    out3ber(siriter, ebnoiter, iter) = out3ber1;
                end
            end
            
         prog = prog + 1/(length(sirs)*length(ebnos)*numiter);
         waitbar(prog, wB, strcat(int2str(round(prog*100)), '%'));
        end

    end
end
toc;
%% Plot
temp = sum(~isnan(bersgen), 3);

%Plain Demodulation BER
subplot(341), hold all;
berstand2 = berstand;
berstand2(isnan(berstand2)) = 0;
berstand3 = sum(berstand2, 3) ./ temp;
plot(ebnos, berstand3)
title('BER with Interference-Ignorant Receiver');
xlabel('E_b/N_o'); ylabel('ber')
set(gca ,'yscale','log');
grid on;
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB');
hold off;

% False Alarm of signal
subplot(342), hold all;
plot(ebnos, mean((sigpresentest-sigpresentreal) > 0, 3))
title('False Alarm Detection of Interferer');
xlabel('E_b/N_o'); ylabel('Error Rate')
set(gca ,'yscale','log');
grid on;
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB');
hold off;

% Non-Detection Error
subplot(343), hold all;
plot(ebnos, mean((sigpresentest-sigpresentreal) < 0, 3))
title('Non-Detection of Interferer');
xlabel('E_b/N_o'); ylabel('Error Rate')
set(gca ,'yscale','log');
grid on;
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB');
hold off;

% Delay error
subplot(344), hold all;
temp2 = sum(~isnan(delayerr), 3);
delayerr2 = delayerr;
delayerr2(isnan(delayerr)) = 0;
delayerr3 = sum(delayerr2, 3) ./ temp2;
plot(ebnos, delayerr3)
title('Delay Estimation Error Rate');
xlabel('E_b/N_o'); ylabel('Delay Error Rate')
set(gca ,'yscale','log');
grid on;
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB');
hold off;


% Joint BER
subplot(345), hold all;
bersgen2 = bersgen;
bersgen2(isnan(bersgen2)) = 0;
bersgen3 = sum(bersgen2, 3) ./ temp;
plot(ebnos, bersgen3);
title('BER of Primary Signal during Synchronization');
xlabel('E_b/N_o'); ylabel('ber')
set(gca ,'yscale','log');
grid on;
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB');
hold off;


%out1 BER
subplot(346), hold all;
out1ber2 = out1ber;
out1ber2(isnan(out1ber2)) = 0;
out1ber3 = sum(out1ber2, 3) ./ temp;
plot(ebnos, out1ber3)
title('BER of Primary Signal during Joint Detection');
xlabel('E_b/N_o'); ylabel('ber')
set(gca ,'yscale','log');
grid on;
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB');
hold off;

%out2 BER
subplot(347), hold all;
out2ber2 = out2ber;
out2ber2(isnan(out2ber2)) = 0;
out2ber3 = sum(out2ber2, 3) ./ temp;
plot(ebnos, out2ber3)
title('BER of Interferer during Joint Detection');
xlabel('E_b/N_o'); ylabel('ber')
set(gca ,'yscale','log');
grid on;
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB');
hold off;

%out3 BER
subplot(348), hold all;
out3ber2 = out3ber;
out3ber2(isnan(out3ber2)) = 0;
out3ber3 = sum(out3ber2, 3) ./ temp;
plot(ebnos, out3ber3)
title('BER with Plain Demodulation of Interferer');
xlabel('E_b/N_o'); ylabel('ber')
set(gca ,'yscale','log');
grid on;
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB');
hold off;

% Channel Gain MSE
temp2 = abs((pgsguess - pgsreal).^2);
temp2(isnan(temp2)) = 0;
pgerr = sum(temp2, 3) ./ temp;
subplot(349), hold all, grid on;
plot(ebnos, pgerr)
title('E_b/N_o Channel Gain Mean-Square Error');
xlabel('E_b/N_o'); ylabel('Channel Gain Mean-Square Error')
set(gca ,'yscale','log');
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB');
hold off;

% Frequency Estimation Error
rotfacerr = abs(rotfacguess - rotfacreal)*180/pi;
freqerr = (rotfacerr / ts) / 360 / 1e3; 
freqerr(isnan(freqerr)) = 0;
freqerr2 = sum(freqerr, 3) ./ temp;
subplot(3,4,10), hold all, grid on;
plot(ebnos, freqerr2);
title('Frequency Estimation Error');
xlabel('E_b/N_o'); ylabel('Frequency Error (KHz)')
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB');
hold off;

subplot(3,4,11), hold all;
title('Hard-Decoding BER without Interference');
berfad = berawgn(ebnos, 'qam', 4);
plot(ebnos, berfad);
xlabel('E_b/N_o'); ylabel('ber')
set(gca ,'yscale','log');
grid on;
hold off;

close(wB);
