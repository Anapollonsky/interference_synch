function [out1, PG2] = KMeans2(data, d2, mod1, demod1, pg)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
    m1 = length(constellation(mod1));
    fgt = .98;
    numiter = 20;
    pgs = ones(numiter, 1);
    assigns = zeros(length(data), numiter);
    assignsmod = zeros(length(data), numiter);
     
    for k = 1:numiter
        for t = 1:length(data)
            assigns(:, k) = (data - pgs(k) * d2) ./ pg;
            assignsmod(:, k) = step(mod1, step(demod1, assigns(:, k)));
            pgs(k+1) = mean((data - assignsmod(:, k) .* pg) ./ d2);
        end
    end
    out1 = assignsmod(:, numiter);
    PG2 = repmat(pgs(end), [length(data) 1]);
end

