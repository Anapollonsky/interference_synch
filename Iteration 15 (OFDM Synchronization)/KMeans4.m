function [out1, PG2, rotfac] = KMeans4(data, d2, mod1, demod1, pg)
% SoftDecisionML1 Performs joint detection, including adaptive equalization.
    numiter1 = 100;
    pgs = ones(length(data), numiter1);
    assigns = zeros(length(data), numiter1);
    assignsmod = zeros(length(data), numiter1);
    rotfac = 1;
    rotthing = ones(length(data) - 1, 1);
    rotangle = ones(length(data) - 1, 1);
    rotvec = ones(length(data), 1);
    
    
    for k = 1:numiter1
        %% Guess demodulation
        assigns(:, k) = (data - pgs(:, k) .* d2) ./ pg;
        assignsmod(:, k) = step(mod1, step(demod1, assigns(:, k)));
        %% Guess Rotation Factor
        for t = 2:length(data)
            rotthing(t-1) = ((data(t) - assignsmod(t, k) .* pg(t)) ./ d2(t)) / ...
                            ((data(t-1) - assignsmod(t-1, k) .* pg(t-1)) ./ d2(t-1));
%             rotthing(t-1) = rotthing(t-1)/abs(rotthing(t-1));
            rotangle(t-1) = angle(rotthing(t-1));
        end
% rotangle
        rotfac = exp(1i*mean(rotangle));
%         rotfac = mean(rotthing);
        %% Guess Beginning
        for t = 1:length(data)
            rotvec(t) = rotfac^(t-1);
        end


        pginit = mean(((data - assignsmod(:, k) .* pg) ./ d2) ./ rotvec);
        pgs(:, k+1) = pginit * rotvec;
    end
    out1 = assignsmod(:, end);
    PG2 = pgs(:, end);
end

