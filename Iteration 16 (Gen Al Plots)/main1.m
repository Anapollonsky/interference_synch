%% Andrew Apollonsky
% This script generates plots for random 
close all force;
clc;
clear all;
tic;
%% Simulation Parameters
len = 12;
m1 = 4;
m2 = 2;

ebnos = [-5:4:15];
sirs = [0];
numiter = 10;

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
bersgen = zeros(length(sirs), length(ebnos), numiter);
berstand = zeros(length(sirs), length(ebnos), numiter);
pgsguess = zeros(length(sirs), length(ebnos), numiter);
rotfacguess = zeros(length(sirs), length(ebnos), numiter);
pgsreal = zeros(length(sirs), length(ebnos), numiter);
numiters = zeros(length(sirs), length(ebnos), numiter);

wB = waitbar(0,'0%');
prog = 0;
%% Calculations
for siriter = 1:length(sirs)
    sir2 = sirs(siriter) - 10*log10(log2(m1)) + 10*log10(log2(m2));
    g2 = 10^(-sir2/20);
    for ebnoiter = 1:length(ebnos)
        SNR = ebnos(ebnoiter) + 10*log10(log2(m1));
        for iter = 1:numiter
            [out1demodgenbin, out1demodstandbin, x1bin, pgbase, rotfac, pathG2, iters] ...
             = processdat(mod1, mod2, demod1, m1, m2, g1, g2, SNR, len, rotfacreal);
            %% Assign to Global
            temp = sum(out1demodgenbin ~= x1bin)/length(x1bin);
            temp2 = sum(out1demodstandbin ~= x1bin)/length(x1bin);
            
            bersgen(siriter, ebnoiter, iter) = temp;
            berstand(siriter, ebnoiter, iter) = temp2;
            pgsguess(siriter, ebnoiter, iter) = pgbase;
            rotfacguess(siriter, ebnoiter, iter) = rotfac;
            pgsreal(siriter, ebnoiter, iter) = pathG2(1);
            numiters(siriter, ebnoiter, iter) = iters;
        end
         prog = prog + 1/(length(sirs)*length(ebnos));
         waitbar(prog, wB, strcat(int2str(round(prog*100)), '%'));
    end
end
%% Plot

% Joint BER
subplot(331), hold all;
plot(ebnos, mean(bersgen, 3))
% plot(ebnos, mean(bersgen, 3))
title('E_b/N_o vs BER with joint synchronization/detection');
xlabel('E_b/N_o'); ylabel('ber')
set(gca ,'yscale','log');
grid on;
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB', 'SIR = InfdB');
hold off;

%Plain Demodulation BER
subplot(332), hold all;
plot(ebnos, mean(berstand, 3))
title('E_b/N_o vs BER with plain demodulation');
xlabel('E_b/N_o'); ylabel('ber')
set(gca ,'yscale','log');
grid on;
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB', 'SIR = InfdB');
hold off;

% Channel Gain MSE
pgerr = mean(abs((pgsguess - pgsreal).^2), 3);
subplot(334), hold all, grid on;
plot(ebnos, pgerr)
title('E_b/N_o Channel Gain Mean-Square Error');
xlabel('E_b/N_o'); ylabel('Channel Gain Mean-Square Error')
set(gca ,'yscale','log');
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB', 'SIR = InfdB');
hold off;

% SIR Estimate Error
SirEstDb = 20*log10(abs(pgsguess)/g1);
SirErrDb = abs(repmat(sirs.', [1 length(ebnos)]) - flipud(mean(SirEstDb, 3)));
subplot(335), hold all, grid on;
plot(ebnos, SirErrDb)
title('E_b/N_o vs Mean SIR Estimate Error');
xlabel('E_b/N_o'); ylabel('SIR Estimate Error (dB)')
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB', 'SIR = InfdB');
hold off;

% Initial Channel Phase Change Error
InitAngleErr = mean(min((2 * pi) - abs(angle(pgsguess) - angle(pgsreal)), ...
    abs(angle(pgsguess) - angle(pgsreal))), 3) * 180/pi;
subplot(336), hold all; grid on;
plot(ebnos, InitAngleErr)
title('E_b/N_o vs Mean Initial Angle Error');
xlabel('E_b/N_o'); ylabel('Initial Angle Error (^o)')
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB', 'SIR = InfdB');
hold off;

% Rotation Angle Mean Error
rotfacerr = mean((abs(rotfacguess - rotfacreal)*180/pi), 3);
subplot(337), hold all; grid on;
plot(ebnos, rotfacerr);
title('E_b/N_o Rotation Angle Error');
xlabel('E_b/N_o'); ylabel('Rotation Angle Error (^o)')
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB', 'SIR = InfdB');
hold off;

% Frequency Estimation Error
temp = (rotfacguess - rotfacreal)*180/pi;
freqerr = (temp / ts) / 360 / 1e3; 
freqerr2 = mean(abs(freqerr), 3);
subplot(338), hold all, grid on;
plot(ebnos, freqerr2);
title('Frequency Estimation Error');
xlabel('E_b/N_o'); ylabel('Frequency Error (KHz)')
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB', 'SIR = InfdB');
hold off;

% Genetic Algorithm Iterations
subplot(339), hold all; grid on;
plot(ebnos, mean(numiters, 3));
title('Iterations for Convergence');
xlabel('E_b/N_o'); ylabel('Number of Iterations')
legend('SIR = -6dB', 'SIR = -3dB', 'SIR = 0dB', 'SIR = 3dB', 'SIR = 6dB', 'SIR = InfdB');
hold off;

close(wB);
toc