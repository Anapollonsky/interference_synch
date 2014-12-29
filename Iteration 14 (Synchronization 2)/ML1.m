function [out1, PG2, err] = ML1(data, d2, mod1, demod1, pg)
% SoftDecisionML1 Performs joint detection, including adaptive equalization.
    mu1 = .001;
    m1 = length(constellation(mod1));
    pgs = ones(length(data), 1);
    assigns = zeros(length(data), length(data));
    err = zeros(length(data), 1);
%     rotfac = 1;
%     rotthing = ones(length(data) - 1, 1);
%     rotangle = ones(length(data) - 1, 1);
%     rotvec = ones(length(data), 1);
    
    pgs(1) = mean(data ./ d2);
%     pgs(1) = 1;

    for t = 1:length(data)
        for t2 = 1:t
            minimat  = pg(t2) * step(mod1, (0:m1-1).');
            maximat = abs(data(t2) - pgs(t2).*d2(t2) - minimat); 
            
            [~, minind] = min(maximat);
            assigns(t2, t) = minind-1;
        end
        err(t) = mean((data(1:t) - step(mod1, assigns(1:t, t)) * pg(t) - pgs(t)'*d2(1:t)).^2);
        pgs(t+1) = pgs(t) + mu1 * d2(t) * err(t)'; 
    end
    out1 = assigns(:, end);
    PG2 = pgs;
       
end

