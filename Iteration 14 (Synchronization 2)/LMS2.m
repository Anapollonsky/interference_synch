function [out1, PG2, err] = LMS2(data, d2, pg)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
    mu1 = .002;
    numiter = 30;
    err = zeros(length(data), 1);
    pgest = ones(length(data), numiter);
    symguess = zeros(length(data), 1);
    
    
    guesstep = .8;
    
    for t = 1:length(data)
        for k = 1:numiter
            %% Estimate PG
            err(t) = data(t) - symguess(t) * pg(t) - pgest(t, k)'*d2(t);
            if k == numiter
                pgest(t + 1, 1) = pgest(t, k) + mu1 * d2(t) * err(t)';
            else
                pgest(t, k + 1) = pgest(t, k) + mu1 * d2(t) * err(t)';
            end
            
            %% Estimate Symbol
            if k == numiter
                newguess = (data(t) - pgest(t + 1, 1)'*d2(t))/pg(t);
            else
                newguess = (data(t) - pgest(t, k + 1)'*d2(t))/pg(t);
            end
            symguess(t) = symguess(t) + (newguess - symguess(t)) * guesstep;
        end
    end
    
    out1 = symguess;
    PG2 = pgest(:, numiter);
end
        