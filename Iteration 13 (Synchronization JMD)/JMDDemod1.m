function [out1, out2, error] = JMDDemod1(data, d1, d2, mod1, mod2, m1, m2, g1, g2)
%JMDDemod1 Performs joint detection, including adaptive equalization.
    taplen = 1;
    mu = .002;

    minimat = zeros(m1, m2);
    out1 = zeros(length(data), 1);
    out2 = zeros(length(data), 1);
    data2 = zeros(length(data), 1);
    maximat = zeros(m1, m2, length(data));
    tapest = zeros(taplen, length(data));
    tapest(1, 1) = 1;
    error = zeros(length(data), 1);
    tapin1 = zeros(taplen, 1);
%     tapin2 = zeros(taplen, 1);
    
    for a = 0:m1-1
        aa = step(mod1, a);
        for b = 0:m2-1
            bb = step(mod2, b);
            minimat(a+1, b+1) = g1* aa + g2* bb;
        end
    end
 
    for t = 1:length(data)
        tapin1 = [data(t); tapin1(1:taplen-1)];
        error(t) = d1(t) - tapest(:, t)'*tapin1;
        tapest(:, t+1) = tapest(:, t) + mu*(tapin1)*error(t)';          
        
%         movdat = [zeros(taplen-t+1, 1); data(t-taplen+1:t)];
%         data2(t) = tapest(:, t)'*flipud(movdat);
        if t <= taplen
            data2(t) = tapest(:, t)'*flipud([zeros(taplen-t, 1); data(1:t)]);
        else
            data2(t) = tapest(:, t)'*flipud(data(t-taplen+1:t));
        end
        
        
        maximat(:, :, t) = abs(repmat(data2(t), [m1 m2]) - minimat);
        
        [~, minind] = min(reshape(maximat(:, :, t), [], 1));
        [a, b] = ind2sub(size(maximat(:, :, 1)), minind);
        out1(t) = a-1;
        out2(t) = b-1;
    end
end

