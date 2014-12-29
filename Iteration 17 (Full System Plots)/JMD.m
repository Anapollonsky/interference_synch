function [out1, out2] = JMD(rx, mod1, mod2, pathG1, pg2)
    m1 = length(constellation(mod1));
    m2 = length(constellation(mod2));
    out1 = zeros(length(rx), 1);
    out2 = zeros(length(rx), 1);
    mini = zeros(length(out1), m1, m2);
    
    for q = 0:(m1-1) % y1
        for p = 0:(m2-1) % y2
            mini(:, q+1, p+1) = pathG1.*step(mod1, q) + pg2.*step(mod2, p);
        end
    end

    minimat = abs(repmat(rx, [1 m1 m2]) - mini);

    for p = 1:length(rx)
        [~, minind] = min(reshape(minimat(p, :, :), [], 1));
        [~, a, b] = ind2sub(size(minimat(p, :, :)), minind);
        out1(p) = a-1;
        out2(p) = b-1;
    end 
end