function [out1, PG2, err] = LMS6(data, d2, mod1, demod1, pg)
% SoftDecisionML1 Performs joint detection, including adaptive equalization.
    mu1 = .02;
    pgs = ones(length(data), 1);
    assigns = zeros(length(data), length(data));
    assignsmod = zeros(length(data), length(data));
    err = zeros(length(data), 1);
%     rotfac = 1;
%     rotthing = ones(length(data) - 1, 1);
%     rotangle = ones(length(data) - 1, 1);
%     rotvec = ones(length(data), 1);
    
    pgs(1) = mean(data ./ d2);
%     pgs(1) = 1;
    for t = 1:length(data)
       assigns(1:t, t) = (data(1:t) - pgs(t) * d2(1:t)) ./ pg(1:t);
       assignsmod(1:t, t) = step(mod1, step(demod1, assigns(1:t, t)));
%        assignsmod = assigns;
       err(t) = mean(data(1:t) - assignsmod(1:t, t) .* pg(1:t) - pgs(t) * d2(1:t));
       pgs(t + 1) = pgs(t) + mu1 * d2(t) * err(t)';
    end
    out1 = assignsmod(:, end);
    PG2 = pgs;
end

