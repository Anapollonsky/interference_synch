function [perr] = jointmd2serr(m1,m2,snr,inr)
% jointmdserr calculates the symbol error for the joint MD estimator for
% both signals.
aub1 = ((m1 - 1)^2 .* m2)./(6*snr) + 2*(1 - 1/(sqrt(m2))) * ((m2 - 1)./(3*inr));
aub2 = ((m2 - 1)^2 .* m1)./(6*inr) + 2*(1 - 1/(sqrt(m1))) * ((m1 - 1)./(3*snr));
perr = min(aub1, aub2);
end

