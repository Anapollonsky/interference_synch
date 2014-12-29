function [out1, PG2, err] = LMS5(data, d2, pg)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
%     mu1 = [.001 .001 .002 .002 .003 .003 .004 .005 .007 .01];
%     mu1 = [.005 .005 .007 .01 .013 .016 .02 .025 .03 .05];
    numiter = 1;
%     mu1 = [.005 * ones(1, 10) .01 * ones(1, 10) .03 * ones(1, 10) .1 * ones(1, 10)];
mu1 = .2;
%     mu1 = ones(1, numiter) * .1;


    
    totiter = length(data) * numiter;
    err = zeros(totiter, 1);
    pgest = ones(totiter, 1);
    symguess = zeros(length(data), 1);
    guesstep = 1;
    
    for t = 1:length(data)
        for k = 1:numiter
            %% Estimate PG
            iterid = (t-1)*numiter + k;   
            err(iterid) = data(t) - symguess(t) * pg(t) - pgest(iterid)'*d2(t);
            pgest(iterid + 1) = pgest(iterid) + mu1(numiter) * d2(t) * err(iterid)';
            newguess = (data(t) - pgest(iterid + 1)'*d2(t))/pg(t);
            symguess(t) = symguess(t) + (newguess - symguess(t)) * guesstep;  
        end
    end
    
    out1 = symguess;
    PG2 = pgest(numiter:numiter:end)'.';
    err = err(numiter:numiter:end);
end
        