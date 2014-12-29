function [delay, sig] = delayfind(rx, mod1, demod1, pathG1, len2, evm1, evm2)
    k = rx(1:len2)./pathG1(1:len2);
    k2 = step(mod1, step(demod1, k));
    evmvec = abs(k2 - k);
    msqerr = zeros(len2+1, 1);
    for delayiter = 1:len2
        msqerr(delayiter) = mean(...
            [abs((evmvec(1:delayiter) - evm1).^2); ...
            abs((evmvec(delayiter + 1:min(len2, delayiter+12)) - evm2).^2)]);
    end
    msqerr(end) = mean(abs((evmvec(1:len2)-evm1).^2));
    if msqerr(end) == min(msqerr(1:len2));
        sig = 0;
        delay = nan(1);
    else
        [~, delay] = min(msqerr(1:len2));
        sig = 1;
    end

end