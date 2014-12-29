function [out1, out2, error1, error2] = JMDDemod2(data, d1, d2, mod1, mod2, m1, m2)
%JMDDemod1 Performs joint detection, including adaptive equalization.
    taplen = 1;
    mu = .0002;

    minimat = zeros(m1, m2);
    out1 = zeros(length(data), 1);
    out2 = zeros(length(data), 1);
    maximat = zeros(m1, m2, length(data));
    tapest1 = zeros(taplen, length(data));
    tapest1(1, 1) = 1;
    tapest2 = zeros(taplen, length(data));
    tapest2(1, 1) = .4;
    
    error1 = zeros(length(data), 1);
    error2 = zeros(length(data), 1);
    tapin1 = zeros(taplen, 1);
    tapin2 = zeros(taplen, 1);

    for t = 1:length(data)
        tapin1 = [data(t); tapin1(1:taplen-1)];
        tapin2 = [data(t); tapin2(1:taplen-1)];
 
        for a = 0:m1-1
            aa = step(mod1, a);
            for b = 0:m2-1
                bb = step(mod2, b);
                minimat(a+1, b+1) = tapest1(t)* aa + tapest2(t)* bb;
            end
        end

        maximat(:, :, t) = abs(repmat(data(t), [m1 m2]) - minimat);
        
        [~, minind] = min(reshape(maximat(:, :, t), [], 1));
        [a, b] = ind2sub(size(maximat(:, :, 1)), minind);
        out1(t) = a-1;
        out2(t) = b-1;
        
        comp1 = step(mod1, a-1);
        comp2 = step(mod2, b-1);
        error1(t) = d1(t) - tapest1(:, t)'*comp1;
        error2(t) = d2(t) - tapest2(:, t)'*comp2;
        tapest1(:, t+1) = tapest1(:, t) + mu*comp1*error1(t)';
        tapest2(:, t+1) = tapest2(:, t) + mu*comp2*error2(t)';
    end
end

