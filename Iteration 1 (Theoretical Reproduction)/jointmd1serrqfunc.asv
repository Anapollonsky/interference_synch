function [perr] = jointmd1serrqfunc(snr, inr, sir)
% jointmdserr calculates the symbol error for the joint MD estimator for
% a single signal with the Q function.
    if sir > 1
        perr = 1/2 * qfunc(sqrt(snr) - sqrt(inr));
    elseif sir > 1/4 
        perr = 1/2 * qfunc(sqrt(inr) - sqrt(snr));
    else
        perr = qfunc(sqrt(snr));
    end
end

