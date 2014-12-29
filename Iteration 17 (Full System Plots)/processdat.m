
% Generate data
x1 = randi([0 m1-1], len, 1);
x2 = randi([0 m2-1], len, 1);

% Modulate
y1 = step(mod1, x1);
y2 = step(mod2, x2);

% Generate path gains
pathG1 = ones(len, 1) * exp(rand(1) * 2*pi*1i) * g1;
pathG2 = ones(len, 1) * exp(rand(1) * 2*pi*1i) * g2;

% Apply rotation to 2nd path gain
for t = 2:len
    pathG2(t) = pathG2(t) * exp((t-1) * 1i * rotfacreal);
end

% Apply path gains to data
rx1 = pathG1 .* y1;
rx2 = pathG2 .* y2;

% Add AWGN
rx1 = awgn(rx1, SNR);

% Sum signals
rxSig = rx1 + rx2;

% Apply genetic algorithm
[out1gen, wts2, iters, pgbase, rotfac] = gen_al4(rxSig, y2, mod1, pathG1);

% Demodulate
out1demodgen = step(demod1, out1gen); 
out1demodstand = step(demod1, rxSig ./ pathG1);

% Convert to binary
out1demodgenbin = reshape(de2bi(out1demodgen, log2(m1)), [], 1);
out1demodstandbin = reshape(de2bi(out1demodstand, log2(m1)), [], 1);
x1bin = reshape(de2bi(x1, log2(m1)), [], 1);

end

