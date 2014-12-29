function [perr] = intserrqfunc(snr, inr, sir)
% jointmdserr calculates the symbol error for the interference ignorant estimator for
% a single signal with the Q function.
    if sir > 1
        perr = 1/2 * qfunc(sqrt(snr) - sqrt(inr));
    else
        perr = 1/2 * (1 - qfunc(sqrt(inr) - sqrt(snr)));
    end
end

