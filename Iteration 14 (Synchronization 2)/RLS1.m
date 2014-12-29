function [out1, PG2, err] = RLS1(data, d2, pg)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
    delta = 10;
    lambda = .98;
    numiter = 10;  
    
    w0 = 1;
    P0 = eye(1)/delta;
    
    totiter = length(data) * numiter;
    P = zeros(1, 1, totiter);
    pgest = ones(totiter, 1);
    pi1 = zeros(totiter, 1);
    K = zeros(totiter, 1);
    err = zeros(totiter, 1);
    symguess = zeros(length(data), 1);
    
    guesstep = .8;
    
    for t = 1:length(data)
        for k = 1:numiter
            %% Estimate PG
            iterid = (t-1)*numiter + k;
            if iterid == 1
                pi1(iterid) = P0 * d2(t);
                K(iterid) = pi1(iterid) / (lambda + d2(t)'*pi1(iterid));
                err(iterid) = data(t)  - symguess(t) * pg(t) - w0'*d2(t);
                pgest(iterid) = w0 + K(iterid)*err(iterid)';
                P(:, :, iterid) = P0 / lambda - K(iterid) * d2(t)' * P0 / lambda;
            else
                pi1(iterid) = P(:, :, iterid-1) * d2(t);
                K(iterid) = pi1(iterid) / (lambda + d2(t)'*pi1(iterid));
                err(iterid) = data(t)  - symguess(t) * pg(t) - pgest(iterid-1)'*d2(t);
                pgest(iterid) = pgest(iterid-1) + K(iterid)*err(iterid)';
                P(:, :, iterid) = P(:, :, iterid-1) / lambda - ...
                K(iterid) * d2(t)' * P(:, :, iterid-1) / lambda;
            end
                
            %% Estimate Symbol
            newguess = (data(t) - pgest(iterid)'*d2(t))/pg(t);
            symguess(t) = symguess(t) + (newguess - symguess(t)) * guesstep;
        end
    end
    out1 = symguess;
    PG2 = pgest(numiter:numiter:end);
    err = err(numiter:numiter:end);
end
        