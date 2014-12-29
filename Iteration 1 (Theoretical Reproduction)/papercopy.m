%% Andrew Apollonsky

clear all;
close all;
clc;

%% Configuration
% Constellation sizes. Some approximations don't utilize them. Q-functions
% assume m1 = m2 = 2.
m1 = 2;
m2 = 2;

% dB values for signal, noise and interference. Used instead of explicit
% SNR/SIR due to SINR requirement for interference-ignorant plots.
sigdb = 0;  
noidb = [0:.5:-40];
intdb = [5 3 1 -1 -3 -5];

minser = 10^-8; %Smallest SER to plot. Affects the y-axis range.
%% Other
% Automated things; finding snr in both db and not, etc.
sig = 10.^(sigdb/10);
noi = 10.^(noidb/10);
int = 10.^(intdb/10);

snrdb = sigdb - noidb;
sirdb = sigdb - intdb;
snr = 10.^(snrdb/10);
sir = 10.^(sirdb/10);


%% Interference Ignorant
arrii = zeros(length(sir), length(snr));
for ind = 1:length(sir)
    sinr = sig./(int(ind) + noi);
%     arrii(ind, :) = intigserr(m1, sinr); %Uses approximation.
    
    arrii(ind, :) = intserrqfunc(snr, snr./sir(ind), sir); %Uses Q-function
end

%% Successive Interference Cancellation
arrsic = zeros(length(sir), length(snr));
for ind = 1:length(sir)
      arrsic(ind, :) = sicserrqfunc(snr, snr./sir(ind), sir); % Uses Q-function.
end

%% Joint MD
arrjmd = zeros(length(sir), length(snr));
for ind = 1:length(sir)
    
%      arrjmd(ind, :) = jointmd1serrep(m1, m2, snr);   %Assumes
%      repetitive coding.

%     arrjmd (ind, :) = jointmd2serr(m1, m2, snr, snr./sir(ind)); %
%     Assumes estimation of both signals, uses approximation.

      arrjmd(ind, :) =  jointmd1serrqfunc(snr, snr./sir(ind), sir); % Uses Q-function.
end

%% Plot
figure;
hold all;
for ind = 1:length(sir)
    plot(snrdb, arrii(ind, :));
%     plot(snrdb, arrsic(ind, :));
%     plot(snrdb, arrjmd(ind, :));
end
hold off;
legend('-5', '-3', '-1', '1', '3', '5');
set(gca,'yscale','log');
title('SER for various detectors');
xlabel('SNR');
ylabel('SER');
axis([min(snrdb) max(snrdb) minser 1])
