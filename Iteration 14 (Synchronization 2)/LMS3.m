function [out1, PG2, err] = LMS3(data, d2, pg)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
    mu1 = .02;
    numiter = 30;

    
    totiter = length(data) * numiter;
    err = zeros(totiter, 1);
    pgest = ones(totiter, 1);
    symguess = zeros(length(data), 1);
    guesstep = .8;
    
    for t = 1:length(data)
        for k = 1:numiter
            %% Estimate PG
            iterid = (t-1)*numiter + k;
            err(iterid) = data(t) - symguess(t) * pg(t) - pgest(iterid)'*d2(t);
            pgest(iterid + 1) = pgest(iterid) + mu1 * d2(t) * err(iterid)';
            newguess = (data(t) - pgest(iterid + 1)'*d2(t))/pg(t);
            symguess(t) = symguess(t) + (newguess - symguess(t)) * guesstep;        
        end
    end
    
    out1 = symguess;
    PG2 = pgest(numiter:numiter:end);
    err = err(numiter:numiter:end);
end
        