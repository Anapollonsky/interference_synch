function [out1, PG2] = LMS1(data, d2, mod1, pg)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
    m1 = length(constellation(mod1));
    mu1 = .1;
    mu2 = .1;
    numiter = 5;
    pgest = 0;
    dataest = zeros(length(data), 1);
    err = zeros(length(data), numiter);
    pgest = zeros(length(data), numiter);
    
    for t = 1:length(data)
        for k = 1:numiter
            err(t, k) = data(t) - pgest(t, k)'*d2(t) - dataest(t)*pg(t);
            if k == numiter
                pgest(t + 1, 1) = pgest(t, k) + mu1 * err(t, k)';
            else
                pgest(t, k + 1) = pgest(t, k) + mu1 * err(t, k)';
            end
            dataest(t) = dataest(t) + mu2 * err(t, k)';
        end
    end
    out1 = dataest;
    PG2 = pgest(:, numiter);
end
        