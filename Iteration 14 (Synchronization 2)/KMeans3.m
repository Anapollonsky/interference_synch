function [out1, PG2, rotfac] = KMeans3(data, d2, mod1, demod1, pg)
% SoftDecisionML1 Performs joint detection, including adaptive equalization.
    numiter1 = 10;
    numiter2 = 10;
    pgs = ones(length(data), numiter1);
    assigns = zeros(length(data), numiter1);
    assignsmod = zeros(length(data), numiter1);
    rotfac = 0;
    rotthing = zeros(length(data) - 1, 1);
    rotvec = zeros(1, length(data));
    
    
    for k = 1:numiter1
        for t = 1:length(data)
            assigns(:, k) = (data - pgs(:, k) .* d2) ./ pg;
            assignsmod(:, k) = step(mod1, step(demod1, assigns(:, k)));
        end
        for p = 1:numiter2
            
            for t = 1:length(data)
                rotvec(t) = rotfac ^ (t-1);
            end
            
%             pgs(:, k+1) = rotvec .* mean((data - assignsmod(:, k) .* pg) ./ d2);
            
            pgs(:, k+1) = mean( ((data - assignsmod(:, k) .* pg) ./ d2) ./ rotvec.')...
                * rotvec.';
            for t = 2:length(data)
%                 size(((data(t) - assignsmod(t, k) .* pg) ./ d2(t)))
                
                rotthing(t-1) = ((data(t) - assignsmod(t, k) .* pg(t)) ./ d2(t)) / ...
                                ((data(t-1) - assignsmod(t-1, k) .* pg(t)) ./ d2(t-1));
                rotfac = mean(rotthing);
            end
            
        end
    end
    out1 = assignsmod(:, numiter1);
    PG2 = repmat(pgs(end), [length(data) 1]);
end

