function [out1, out2, error1, error2, tapest1, tapest2] = ...
    SoftDecisionZF(data, d1, d2, demod1, demod2)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
%     taplen = 1;
%     mu = .04;

    tapest1 = zeros(length(data) + 1);
    tapest1(:) = 1;
    tapest2 = zeros(length(data) + 1);
    tapest2(:) = .7;
    
    out1 = zeros(length(data), 1);
    out2 = zeros(length(data), 1);
    error1 = zeros(length(data), 1);
    error2 = zeros(length(data), 1);

    for t = 1:length(data)
        k = [tapest1(t); tapest2(t)];
        W = (k'*k)\k;
        out = W' * data(t);
        out1(t) = step(demod1, out(1));
        out2(t) = step(demod2, out(2));
        
        error1(t) = d1(t) - out(1);
        error2(t) = d2(t) - out(2);
%         tapest1(t+1) = tapest1(t) + mu*error1(t);
%         tapest2(t+1) = tapest2(t) + mu*error2(t);
    end
end


