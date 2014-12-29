%% Andrew Apollonsky
close all;
clc;
clear all;
m1 = 2;
m2 = 4;
tlen = 200;
len = 5000;
poffset1 = 0;


for iter = 1
    slen1 = log2(m1);
    slen2 = log2(m2);

    EbNo1 = [10];
    sir = 5;
    snr1 = EbNo1(iter) + 3 + 10*log10(slen1);

    mod1 = comm.RectangularQAMModulator(m1);
    mod2 = comm.RectangularQAMModulator(m2);
    demod1 = comm.RectangularQAMDemodulator(m1);
    demod2 = comm.RectangularQAMDemodulator(m2);

    hchan1 = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
    hErrorCalc1 = comm.ErrorRate;
    hErrorCalc2 = comm.ErrorRate;

    x1 = randi([0 m1-1], len, 1); % Actual data
    x2 = randi([0 m2-1], len, 1);
    y1 = step(mod1, x1); % Modulated data
    y2 = step(mod2, x2);
    g1 = 1;               % Gain
    g2 = 10^(-sir/20);

    hchan1.SignalPower = (real(y1)'*real(y1))/length(real(y1));
    hchan1.SNR = snr1;


    %% Channel Distortion
    dist1 = [.986; .845; .237; .123+.31i];
    dist2 = [.97; .5; .2];

    w1 = step(hchan1, filter(dist1, 1, y1));
%     w1 = filter(dist1, 1, y1);
    w2 = filter(dist2, 1, y2);

    %% Equalize
%     eq1 = lineareq(8, lms(.01));
%     eq1.SigConst = step(mod1, (0:m1-1)')';
%     [symbolest, yd] = equalize(eq1, w1, y1(1:tlen));
    eq1 = adaptfilt.lsl(32, .995, 1);
    [symbolest, yd] = filter(eq1, w1(1:tlen), y1(1:tlen));

%     h = scatterplot(w1, 1, tlen, 'bx');
%     hold on;
%     scatterplot(symbolest, 1, tlen, 'g.', h);
%     scatterplot(eq1.SigConst, 1, 0, 'k*', h);

    %% Normal Demodulation
    demnoeq = step(demod1, w1);
    demeq = step(demod1, symbolest);
    sernoeq = step(hErrorCalc1, x1(tlen+1:end), demnoeq(tlen+1:end));
    reset(hErrorCalc1);
    sereq = step(hErrorCalc1, x1(tlen+1:end), demeq(tlen+1:end));
    sernoeq1(iter) = sernoeq(1)
    sereq1(iter) = sereq(1)
end

% axis([-1.5 1.5 -1.5 1.5]);
% hold all;
% plot(EbNo1, sernoeq1);
% plot(EbNo1, sereq1);
% legend('No Equalization', 'Equalization');
% hold off;
% xlabel('E_b/N_o');
% ylabel('SER');