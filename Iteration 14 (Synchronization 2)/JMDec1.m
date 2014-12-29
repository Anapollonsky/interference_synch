function [out1, out2, pg1, pg2] = JMDec1(data, pg11, pg12, ...
    mod1, mod2)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
    m1 = length(constellation(mod1));
    m2 = length(constellation(mod2));
    mu = .3;
    out1 = zeros(length(data), 1);
    out2 = zeros(length(data), 1);
%     out1mod = zeros(length(data), 1);
%     out2mod = zeros(length(data), 1);
    pg1 = zeros(length(data), 1);
    pg1(1) = pg11;
    pg2 = zeros(length(data), 1); 
    pg2(1) = pg12;
    minimat = zeros(m1, m2, length(data));
    maximat = zeros(m1, m2, length(data));
    
    amat = zeros(m1, length(data));
    bmat = zeros(m2, length(data));
    
%     numiter = 5;
    
    tempiter = 1;
    for t = 1:length(data)
        sirl = pg1(tempiter) / pg2(tempiter);
        mu1 = mu * (sirl / (sirl + 1));
        mu2 = mu * (1/(sirl + 1));
        for a = 0:m1-1
            aa = step(mod1, a);
            for b = 0:m2-1
                bb = step(mod2, b);
                minimat(a+1, b+1, t) = pg1(tempiter) * aa + pg2(tempiter) * bb;
            end
        end
        
        maximat(:, :, t) = abs(repmat(data(t), [m1 m2]) - minimat(:, :, t)).^2;
        amat(:, t) = sum(maximat(:, :, t), 2);
        bmat(:, t) = sum(maximat(:, :, t), 1).';
        
%         [~, minind] = min(reshape(maximat(:, :, t), [], 1));
%         [a, b] = ind2sub(size(maximat(:, :, 1)), minind);
%         out1(t) = a-1;
%         out2(t) = b-1;
%         out1mod = step(mod1, out1(t));
%         out2mod = step(mod2, out2(t));
        
        [~, a] = min(amat(:, t));
        [~, b] = min(bmat(:, t));
        out1(t) = a-1;
        out2(t) = b-1;
        out1mod = step(mod1, out1(t));
        out2mod = step(mod2, out2(t));



        newguess = (data(t) - pg2(tempiter)*out2mod)/out1mod;
        pg1(t) = pg1(tempiter) + mu1 * (newguess - pg1(tempiter));
        newguess = (data(t) - pg1(tempiter)*out1mod)/out2mod;
        pg2(t) = pg2(tempiter) + mu2 * (newguess - pg2(tempiter));
        
%         pg1(t) = mean(data - ass
%         pg1(:) = mean(
        
        
        tempiter = t;
    end
end

