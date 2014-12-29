SIRS = [-10:.2:10];
SNRS = [0:.2:30];

[evm1, evm2] = errcalc(SIRS, SNRS, 1e3, mod1, demod1, bpskmod);


close all;
subplot(1, 3, 1);
h1 = pcolor(SNRS, SIRS, evm1);
ylabel('SIR (dB)');
xlabel('SNR (dB)');
title('EVM_1');
colorbar;
% colormap(gray)
set(h1, 'EdgeColor', 'none');

subplot(1, 3, 2);
h2 = pcolor(SNRS, SIRS, evm2);
ylabel('SIR (dB)');
xlabel('SNR (dB)');
title('EVM_2');
colorbar;
% colormap(gray)
set(h2, 'EdgeColor', 'none');

subplot(1, 3, 3);
h3 = pcolor(SNRS, SIRS, evm2 - evm1);
ylabel('SIR (dB)');
xlabel('SNR (dB)');
title('EVM_2 - EVM_1');
colorbar;
% colormap(gray)
set(h3, 'EdgeColor', 'none');