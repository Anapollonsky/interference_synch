function [out1, PG2, rotfac] = KMeans5(data, d2, mod1, demod1, pg)
% SoftDecisionML1 Performs joint detection, including adaptive equalization.
    numiter1 = 300;
    numiter2 = 30;
    pgs = ones(length(data), numiter1);
    assigns = zeros(length(data), numiter1);
    assignsmod = zeros(length(data), numiter1);
    rotfac = 1;
    rotthing = ones(length(data) - 1, 1);
    rotvec = ones(length(data), 1);
   
    for k = 1:numiter1
        %% Guess demodulation
        assigns(:, k) = (data - pgs(:, k) .* d2) ./ pg;
        assignsmod(:, k) = step(mod1, step(demod1, assigns(:, k)));
        %% Guess Rotation Factor
        for t = 2:length(data)
            rotthing(t-1) = ((data(t) - assignsmod(t, k) .* pg(t)) ./ d2(t)) / ...
                            ((data(t-1) - assignsmod(t-1, k) .* pg(t-1)) ./ d2(t-1));
        end
%         rotthing
        rotfac = mean(rotthing);
        rotfac = rotfac / abs(rotfac);
        %% Guess Beginning
        for t = 1:length(data)
            rotvec(t) = rotfac^(t-1);
        end
        
        pginit = mean(((data - assignsmod(:, k) .* pg) ./ d2) ./ rotvec);
        pgs(:, k+1) = pginit * rotvec;
    end
    
    for k = 1:length(data)
        pginit = mean
        
    
    out1 = assignsmod(:, end);
    PG2 = pgs(:, end);
end

