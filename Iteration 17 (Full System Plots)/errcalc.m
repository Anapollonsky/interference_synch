function [evm1, evm2] = errcalc(sirs, snrs, iters, mod1, demod1, bpskmod)
    m1 = length(constellation(mod1));
    m2 = 2;
    g1 = 1;
    
    x1 = randi([0 m1-1], iters, 1);
    x2 = randi([0 1], iters, 1);
    y1 = step(mod1, x1);
    y2 = step(bpskmod, x2);
    
    evm1 = zeros(length(sirs), length(snrs));
    evm2 = zeros(length(sirs), length(snrs));
    for k = 1:length(sirs)
        g2 = 10^(-sirs(k)/20);
        for m = 1:length(snrs)
            pathG1 = exp(rand(iters, 1) * 2*pi*1i) * g1;
            z1 = y1 .* pathG1;
            z2 = y2 .* exp(rand(iters, 1) * 2*pi*1i) * g2;
            rx1 = awgn(z1, snrs(m), 'measured');
            rx2 = awgn(z1, snrs(m), 'measured') + z2;
            evm1(k, m) = mean(abs(rx1./pathG1 - step(mod1, step(demod1, rx1./pathG1))));
            evm2(k, m) = mean(abs(rx2./pathG1 - step(mod1, step(demod1, rx2./pathG1))));
        end
    end
    
end