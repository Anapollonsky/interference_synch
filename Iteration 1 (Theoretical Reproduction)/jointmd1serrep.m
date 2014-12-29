function [perr] = jointmd1serrep(m1,m2,snr)
% jointmdserr calculates the symbol error for the joint MD estimator for
% one signal with repetition coding.
perr = 4 * (m1 - 1)^3 * m2^2 ./ (9 * snr.^2)
end

